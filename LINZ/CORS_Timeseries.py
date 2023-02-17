from __future__ import annotations

import sys
import os
import os.path
import re
import sqlite3
import json
import logging

import numpy as np
import datetime as dt
import pandas as pd
from scipy import stats

from LINZ.Geodetic.Ellipsoid import GRS80

from typing import Any, Callable, NewType


def robustStandardError(obs: np.ndarray, percentile: float = 95.0) -> np.ndarray:
    """
    Estimate the standard error for components of a time series.

    Standard errors are estimated using the differences between consecutive elements
    of the time series.  The 95 percentile of the absoluted differences is used as a
    robust estimator for the 95% cumulative probability (two tailed), so is scaled by
    1.96 to get standard error.
    """

    ppf = stats.norm.ppf((1.0 + percentile / 100.0) / 2.0)
    errors = [0] * 3
    for axis in range(3):
        diffs = np.abs(obs[1:, axis] - obs[:-1, axis])
        # Note sqrt(2) accounts for the fact that these are differences
        # between two observations
        se = np.percentile(diffs[~np.isnan(diffs)], percentile) / (ppf * np.sqrt(2))
        errors[axis] = se
    return np.array(errors)


def findOutliers(
    enu_data: pd.DataFrame,
    ndays: int = 10,
    tolerance: int | float = 5.0,
    percentile: int | float = 95.0,
    goodrows: bool = False,
) -> pd.Index:
    """
    Detect outliers in the time series and return either the indexes of outliers
    (goodrows=False) or of the non-outliers (goodrows=True).  Assumes that enu_data
    is a dataframe with columns e, n, u and a datetime index.

    Based on the difference between the observated value and a local median of
    values within ndays of the value being tested.  Outliers are points with an
    E,N, or U value more than tolerance times the robust standard error from the
    median.  The robust standard error is calculated with a specific percentage.

    Note: calculating the medians is slow!
    """
    testrange = pd.DateOffset(days=ndays)
    idx = enu_data.index
    obs = np.vstack((enu_data.e, enu_data.n, enu_data.u)).T
    rse = robustStandardError(obs, percentile)
    medians = np.array(
        [
            np.median(
                obs[np.logical_and(idx >= d - testrange, idx <= d + testrange)], axis=0
            )
            for d in idx
        ]
    )
    if goodrows:
        rows = np.where(np.all(np.abs(obs - medians) <= rse * tolerance, axis=1))
    else:
        rows = np.where(np.any(np.abs(obs - medians) > rse * tolerance, axis=1))
    return idx[rows]


def normalizeSolutionType(solutiontype: str | list[str]) -> str:
    if isinstance(solutiontype, list):
        solutiontype = "+".join(solutiontype)
    return solutiontype


TransformFunc = Callable[[np.ndarray], np.ndarray]


class Timeseries(object):
    _loaded: bool
    _isfunction: bool
    _code: str
    _solutiontype: str
    _xyz0: np.ndarray | None
    _xyzenu: np.ndarray | None
    _transform: TransformFunc | None
    _before: dt.datetime | None
    _after: dt.datetime | None
    _normalize: bool
    _sourcedata: pd.DataFrame | None

    def __init__(
        self,
        code: str,
        solutiontype: str = "default",
        data: pd.DataFrame | None = None,
        dates: np.ndarray = None,
        xyz: np.ndarray = None,
        xyz0: np.ndarray = None,
        xyzenu: np.ndarray = None,
        transform: TransformFunc | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
    ):
        """
        Create a time series. Parameters are:

        code            the station code being loaded
        solutiontype    a string identifying the solution type of the time series
        data            a pandas data frame with columns indexed on date and with
                        columns x, y, z
        dates, xyz      Alternative format for loading, dates is an array of dates,
                        xyz is a numpy array with shape (n,3)
        xyz0            Reference xyz coordinate for calculated enu components
        xyzenu          Reference xyz coordinate for calculated enu components
        transform       A transformation function applied to the XYZ coordinates
        after           The earliest date of interest
        before          The latest date of interest
        normalize       Normalize the observation dates
        """

        self._loaded = False
        self._isfunction = False
        self._code = code
        self._solutiontype = normalizeSolutionType(solutiontype)
        self._xyz0 = None if xyz0 is None else np.array(xyz0)
        self._xyzenu = None if xyzenu is None else np.array(xyzenu)
        self._transform = transform
        self._before = None
        self._after = None
        self._normalize = normalize
        self.setDateRange(after=after, before=before)
        if xyz is not None and dates is not None:
            data = pd.DataFrame(xyz, columns=("x", "y", "z"))
            data.set_index(pd.to_datetime(dates), inplace=True)
        self._sourcedata = data

    def _loadData(self) -> pd.DataFrame:
        """
        Abstract class, overridden in base classes
        """
        assert self._sourcedata is not None
        return self._sourcedata.copy()

    def _load(self) -> None:

        if self._loaded:
            return

        # Call provider specific data load
        # Should return a dataframe with x,y,z indexed on the date/time of the each point

        self._loaded = True
        data = self._loadData()
        if data is None or data.shape[0] == 0:
            message = [
                f"No {self._solutiontype} timeseries data found for {self._code}"
            ]
            if self._after is not None:
                message.append(f"after {self._after:%Y-%m-%d}")
            if self._before is not None:
                if self._after:
                    message.append("and")
                message.append(f"before {self._before:%Y-%m-%d}")
            raise RuntimeError(" ".join(message))

        if self._after is not None:
            data = data[data.index > self._after]

        if self._before is not None:
            data = data[data.index < self._before]

        if self._normalize:
            data.set_index(data.index.normalize(), inplace=True)

        data.sort_index(inplace=True)
        xyz = np.vstack((data.x, data.y, data.z)).T
        if self._transform:
            xyz = np.array([self._transform(x) for x in xyz])
            data["x"] = xyz[:, 0]
            data["y"] = xyz[:, 1]
            data["z"] = xyz[:, 2]

        xyz0 = self._xyz0
        xyz0 = np.array(xyz0) if xyz0 is not None else xyz[0]
        xyzenu = self._xyzenu
        xyzenu = np.array(xyzenu) if xyzenu is not None else xyz0
        lon, lat = GRS80.geodetic(xyzenu)[:2]
        enu_axes = GRS80.enu_axes(lon, lat)

        diff = xyz - xyz0
        enu = (xyz - xyz0).dot(enu_axes.T)
        data["e"] = enu[:, 0]
        data["n"] = enu[:, 1]
        data["u"] = enu[:, 2]

        self._xyz0 = xyz0
        self._data = data

    def setXyzTransform(
        self,
        xyz0: np.ndarray | None = None,
        xyzenu: np.ndarray | None = None,
        transform: TransformFunc | None = None,
    ) -> None:
        xyz0 = None if xyz0 is None else np.array(xyz0)
        xyzenu = None if xyzenu is None else np.array(xyzenu)
        if xyz0 is not None and (self._xyz0 is None or (self._xyz0 != xyz0).any()):
            self._xyz0 = xyz0
            self._loaded = False
        if xyzenu is not None and (
            self._xyzenu is None or (self._xyzenu != xyzenu).any()
        ):
            self._xyzenu = xyzenu
            self._loaded = False
        if transform != self._transform:
            self._transform = transform
            self._loaded = False

    def setDateRange(
        self,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        reset: bool = False,
    ) -> None:
        if reset:
            self._after = None
            self._before = None
        if after is not None and after != self._after:
            if isinstance(after, str):
                after = dt.datetime.strptime(after, "%Y-%m-%d")
            self._after = after
            self._loaded = False
        if before is not None and before != self._before:
            if isinstance(before, str):
                before = dt.datetime.strptime(before, "%Y-%m-%d")
            self._before = before
            self._loaded = False

    def xyz0(self) -> np.ndarray:
        if self._xyz0 is None:
            self._load()
        assert self._xyz0 is not None
        return self._xyz0

    def xyzenu(self) -> np.ndarray:
        if self._xyzenu is None:
            self._load()
        assert self._xyzenu is not None
        return self._xyzenu

    def code(self) -> str:
        return self._code

    def solutiontype(self) -> str:
        return self._solutiontype

    def setName(self, code: str) -> None:
        self._code = code

    def trend(self, columns: list[str] | str = "e n u".split()) -> pd.DataFrame:
        import matplotlib.dates as mdates

        self._load()
        columns = [columns] if isinstance(columns, str) else list(columns)
        series = self._data[columns]
        trends = []
        days = mdates.date2num(self._data.index.to_pydatetime())
        pfit = np.polyfit(days, series, 1)
        trends = [np.poly1d(pfit[:, i])(days) for i, c in enumerate(columns)]
        return pd.DataFrame(
            np.vstack(trends).T, index=self._data.index, columns=columns
        )

    def getObs(
        self,
        enu: bool = True,
        detrend: bool = False,
        normalize: bool = False,
        index: pd.Index = None,
    ) -> tuple[list[dt.datetime], np.ndarray]:
        """
        Returns the time series as date and enu/xyz arrays

        Parameters:
            enu         Select the ENU data rather than XYZ (default True)
            detrend     Remove trend from observations
            index       Pandas data frame index to extract subset of data
            normalize   Normalize dates in time series

        Returns two values:
            dates    An array of dates for the time series
            enu      An numpy array of east,north,up values relative to the
                     model reference coordinate, or xyz coordinates if enu is False
        """
        self._load()
        cols = ("e", "n", "u") if enu else ("x", "y", "z")
        data = self.getData(enu=enu, detrend=detrend, normalize=normalize, index=index)
        return data.index.to_pydatetime(), np.vstack([data[x] for x in cols]).T

    def getData(
        self,
        enu: bool = True,
        detrend: bool = False,
        index: bool = None,
        normalize: bool = False,
    ) -> pd.DataFrame:
        """
        Returns the time series as a pandas DataFrame.

        Parameters:
            enu         Select the ENU data rather than XYZ (default True)
            detrend     Remove trend from observations
            index       Pandas data frame index to extract subset of data
            normalize   Normalize dates in time series

        Returns:
            data        A pandas dataframe with columns e,n,u or x,y,z and a date index.
        """
        self._load()
        columns = ["e", "n", "u"] if enu else ["x", "y", "z"]
        result = self._data[columns]
        if detrend:
            trend = self.trend(columns)
            result -= trend
        if index is not None:
            result = result.loc[index]
        if normalize:
            result.set_index(result.index.normalize(), inplace=True)
        return result

    def dates(self) -> pd.Index:
        """
        Returns the dates of the time series
        """
        self._load()
        return self._data.index

    PlotDefDict = dict[str, Any]

    def plot(
        self,
        detrend: bool = False,
        enu: bool = True,
        independent: bool = False,
        samescale: bool = False,
        mmunits: bool = False,
        symbol: str | None = None,
        baseplot: PlotDefDict = None,
        figtitle: bool = True,
        **kwds: dict[str, Any],
    ) -> PlotDefDict:
        """
        Plot the time series onto 3 separate graphs

           detrend=True to remove the trend from the plots
           samescale=True to force the X,Y,Z axes to share the same scale
           independent=True to calculate trends independently from the baseplot
           mmunits=True to use millimetres for units rather than metres
           symbol='' to plot using a specified symbol
           figtitle to specify the figure title (None or False for no title)
           baseplot to plot againt a specified baseplot

           Additional keywords are passed to the pyplot.subplots function
           call or pyplot.plot call (depending on the keyword, the most common
           plot keywords are passed to plot).
        """

        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates

        default_symbols = ["b+", "r+", "g+", "m+"]
        default_lines = ["b-", "r-", "g-", "m-"]

        plotkwds = {}
        for kw in kwds:
            if (
                kw in "aa antialiased alpha c color lw ls mec mew mfc mfcalt ms".split()
                or kw.startswith("marker")
                or kw.startswith("line")
            ):
                plotkwds[kw] = kwds[kw]
        for kw in plotkwds:
            del kwds[kw]

        havebase = True
        if baseplot is None:
            havebase = False
            baseplot = {}
            title = "{0} {1} timeseries".format(self._code, self._solutiontype)
            if detrend:
                title = title + " detrended"

            settings = {
                "sharex": True,
                "sharey": samescale,
                "figsize": (8, 6),
                "dpi": 100,
            }
            settings.update(kwds)
            fig, plots = plt.subplots(3, 1, **settings)
            if isinstance(figtitle, str):
                title = figtitle
            if figtitle is not None and figtitle is not False:
                fig.suptitle(title)
            baseplot["figure"] = fig
            baseplot["plots"] = plots
            baseplot["trends"] = [None, None, None]
            baseplot["mmunits"] = mmunits
            baseplot["enu"] = enu
            baseplot["symbols"] = []
        else:
            enu = baseplot["enu"]
            mmunits = baseplot["mmunits"]
            self.setXyzTransform(xyz0=baseplot["xyz0"], xyzenu=baseplot["xyzenu"])

        if enu:
            axis_labels = ("East", "North", "Up")
            columns = ("e", "n", "u")
        else:
            axis_labels = ("X", "Y", "Z")
            columns = ("x", "y", "z")

        self._load()
        data = self._data
        if not havebase:
            baseplot["xyz0"] = self._xyz0
            baseplot["xyzenu"] = self._xyzenu
        plots = baseplot["plots"]

        if symbol is None:
            defaults = default_lines if self._isfunction else default_symbols
            colours = [s[:1] for s in baseplot["symbols"]]
            possible = None
            for s in defaults:
                if s[:1] not in colours:
                    symbol = s
                    baseplot["symbols"].append(s)
                    break
                if s not in baseplot["symbols"]:
                    possible = s
            if symbol is None:
                symbol = possible
            if symbol is None:
                symbol = defaults[0]

        for i, axis in enumerate(columns):
            series = data[axis]
            ylabel = axis_labels[i]
            trendp = None
            days = None
            if detrend and (independent or not havebase):
                days = mdates.date2num(data.index.to_pydatetime())
                trendp = np.poly1d(np.polyfit(days, series, 1))
                baseplot["trends"][i] = trendp
            else:
                trendp = baseplot["trends"][i]
            if trendp is not None:
                if days is None:
                    days = mdates.date2num(data.index.to_pydatetime())
                trend = trendp(days)
                series = series - trend
            if mmunits:
                series = series * 1000
                ylabel = ylabel + " mm"
            plots[i].plot(
                data.index,
                series,
                symbol,
                label=self._code + " " + self._solutiontype,
                picker=5,
                **plotkwds,
            )
            if not havebase:
                plots[i].set_ylabel(ylabel)
                plots[i].tick_params(labelsize=8)
                plots[i].format_xdata = mdates.DateFormatter("%d-%m-%Y")
        return baseplot

    def subtract(
        self,
        other: Timeseries,
        newcode: str | None = None,
        newtype: str | None = None,
        normalize: bool = False,
    ) -> Timeseries:
        """
        Returns a time series subtracting other from the current series.  Other
        can be another time series, or a function that takes a date as input
        and returns an XYZ coordinate

        By default requires that both series have the same code.  Will
        return an error if not.  Over-ride by including a newcode parameter
        """

        d1 = self.getData(enu=False, normalize=normalize)
        if isinstance(other, Timeseries):
            if newcode is None:
                if self._code != other._code:
                    newcode = other._code + "-" + self._code
                else:
                    newcode = self._code
            if newtype is None:
                if self._solutiontype != other._solutiontype:
                    newtype = other._solutiontype + "-" + self._solutiontype
                else:
                    newtype = self._solutiontype
            d2 = other.getData(enu=False, normalize=normalize)
            join = d1.join(d2, rsuffix="2", how="inner")
            data = pd.DataFrame(
                data={
                    "x": join.x - join.x2,
                    "y": join.y - join.y2,
                    "z": join.z - join.z2,
                }
            )
        elif callable(other):
            newcode = newcode or self._code
            newtype = newtype or self._solutiontype
            index = d1.index
            diffdata = []
            for i, d in zip(index, index.to_pydatetime()):
                xyz1 = d1.loc[i]
                xyz2 = other(d)
                diffdata.append([xyz1.x - xyz2[0], xyz1.y - xyz2[1], xyz1.z - xyz2[2]])
            data = pd.DataFrame(
                data=diffdata, index=index.copy(), columns=("x", "y", "z")
            )
        else:
            raise RuntimeError("Invalid value to subtract in Timeseries.subtract")
        xyz0 = np.array([0, 0, 0])
        xyzenu = self._xyzenu if self._xyzenu is not None else self._xyz0
        return Timeseries(newcode, newtype, data=data, xyz0=xyz0, xyzenu=xyzenu)

    def resample(
        self, rule: str | int, normalize: bool = True, how: str | None = None
    ) -> Timeseries:
        """
        Resamples the time series using pd.DataFrame.resample.

        rule can be one of the pandas.tseries rules such as M (month), Y (year), or an integer
        number, meaning a number of days

        how by default is mean.  Good alternatives are 'median'
        """

        loffset = None
        if type(rule) == int:
            ndays = rule
            rule = pd.tseries.offsets.Day(ndays)
            loffset = pd.tseries.offsets.Day(int(ndays / 2))
        d1 = self.getData(enu=False, normalize=normalize).resample(
            rule, how=how, loffset=loffset
        )
        return Timeseries(
            self._code,
            self._solutiontype,
            data=d1,
            xyz0=self._xyz0,
            xyzenu=self._xyzenu,
        )

    def robustStandardError(self, percentile: float | int = 95.0) -> np.ndarray:
        """
        Returns the "robust standard error" of the time series E,N,U components
        as a numpy array.
        """
        return robustStandardError(self.getObs()[1], percentile)

    def findOutliers(
        self,
        ndays: int = 10,
        tolerance: float | int = 5.0,
        percentile: float | int = 95.0,
        goodrows: bool = False,
    ) -> pd.Index:
        """
        Return index of outliers or good rows.
        """
        return findOutliers(
            self.getData(),
            ndays=ndays,
            tolerance=tolerance,
            percentile=percentile,
            goodrows=goodrows,
        )

    def withoutOutliers(
        self,
        ndays: int = 10,
        tolerance: float | int = 5.0,
        percentile: float | int = 95.0,
    ) -> Timeseries:
        """
        Return a new version of the time series without the outliers
        """
        index = self.findOutliers(
            ndays=ndays, tolerance=tolerance, percentile=percentile, goodrows=True
        )
        return self.filtered(index=index)

    def filtered(
        self,
        index: pd.Index | None = None,
        before: dt.datetime | None = None,
        after: dt.datetime | None = None,
    ) -> Timeseries:
        """
        Return a subset of the data filtered to a set of index values
        """
        data = self.getData(enu=False)
        filter = data.index.copy()
        if index is not None:
            filter = filter.intersection(index)
        if before is not None:
            filter = filter[filter <= before]
        if after is not None:
            filter = filter[filter >= after]
        data = data.loc[filter]
        return self.usingData(data)

    def offsetBy(
        self, offset: pd.DataFrame | np.ndarray | list[float] | float
    ) -> Timeseries:
        """
        Returns a new Timeseries created by adding an offset to the
        time series.

        The offset can be an [x,y,z] vector, or a time indexed
        data frame of x, y, z values.  In the latter case common times
        are included in the final data.
        """
        data = self.getData(enu=False)
        diff = data + offset
        diff = diff[np.isfinite(diff.x)]
        return self.usingData(diff)

    def usingData(self, data: np.ndarray) -> Timeseries:
        """
        Returns a copy of the current time series but using data
        specified.  Keeps code, solutiontype, xyz0, xyzenu.
        Discards transform.
        """
        result = Timeseries(
            self._code,
            data=data,
            solutiontype=self._solutiontype,
            xyz0=self._xyz0,
            xyzenu=self._xyzenu,
            transform=None,
        )
        result._isfunction = self._isfunction
        return result

    class Comparison(object):
        ts1: Timeseries
        ts2: Timeseries
        diff: Timeseries

        def __init__(
            self,
            ts1: Timeseries,
            ts2: Timeseries,
            newcode: str | None = None,
            newtype: str | None = None,
            normalize: bool = False,
        ):
            """
            Represents a comparison of two time series, with elements ts1, ts2, diff
            Diff is in the sense ts1-ts2
            """
            self.ts1 = ts1
            self.ts2 = ts2
            self.diff = ts1.subtract(
                ts2, newcode=newcode, newtype=newtype, normalize=normalize
            )

        def plot(self, detrend: bool = True) -> None:
            """
            Plot the two time series overlaid
            """
            p1 = self.ts1.plot(detrend=detrend)
            p1 = self.ts2.plot(baseplot=p1)

        def plotDiff(self) -> None:
            """
            Plot the difference
            """
            self.diff.plot()

        def stats(self) -> pd.DataFrame:
            """
            Return the stats of the difference
            """
            return self.diff.getData().describe()

        def printDetrendedStats(self) -> None:
            """
            Print the stats for each time series.
            """
            ts1 = self.ts1
            ts2 = self.ts2
            rsets1 = ts1.robustStandardError()
            rsets2 = ts2.robustStandardError()
            dts1 = ts1.getData() - ts1.trend()
            dts2 = ts2.getData() - ts2.trend()
            print(self.ts1.solutiontype(), "solution detrended stats")
            print(dts1.describe())
            print("Robust SE: ", rsets1)
            print(self.ts2.solutiontype(), "solution detrended stats")
            print(dts2.describe())
            print("Robust SE: ", rsets2)


class SqliteTimeseries(Timeseries):
    _dbfile: str

    _sql = """
        select epoch, code as code, X as x, Y as y, Z as z 
        from mark_coordinate m
        where code=? and solution_type=?
        {when}
        order by epoch
        """

    _sqlMultiType = """
        select
             m.epoch, m.code as code, m.solution_type, m.X as x, m.Y as y, m.Z as z 
        from
             mark_coordinate m,
             (select
                 epoch,
                 min({case}) as version
              from 
                 mark_coordinate
              where code=? and
                    {case} < ?
              group by
                 epoch
              ) as v
        where 
            m.code = ? and
            m.epoch=v.epoch and
            {case} = v.version
            {when}
        order by m.epoch
        """

    _sqlList = """
        select distinct code, solution_type
        from mark_coordinate
        order by code, solution_type
        """

    @staticmethod
    def _openDb(dbfile: str) -> sqlite3.Connection:
        if not os.path.exists(dbfile):
            raise RuntimeError(str(dbfile) + " does not exist")
        db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES)
        return db

    @staticmethod
    def seriesList(
        dbfile,
        solutiontype: str | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
    ) -> list[SqliteTimeseries]:
        db = SqliteTimeseries._openDb(dbfile)
        seriescodes = pd.read_sql(SqliteTimeseries._sqlList, db)
        db.close()
        series = []
        found = []
        stypes = solutiontype.split("+") if solutiontype else None
        for i in seriescodes.index:
            code, solntype = (seriescodes.code[i], seriescodes.solution_type[i])
            if stypes is None or solntype in stypes:
                solntype = solutiontype or solntype
                if (code, solntype) in found:
                    continue
                series.append(
                    SqliteTimeseries(
                        dbfile,
                        code,
                        solntype,
                        after=after,
                        before=before,
                        normalize=normalize,
                    )
                )
                found.append((code, solntype))
        return series

    def __init__(
        self,
        dbfile: str,
        code: str,
        solutiontype: str = "default",
        xyz0: np.ndarray | None = None,
        transform: TransformFunc | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
    ):
        Timeseries.__init__(
            self,
            code,
            solutiontype=solutiontype,
            xyz0=xyz0,
            transform=transform,
            after=after,
            before=before,
            normalize=normalize,
        )
        self._dbfile = dbfile

    def _loadData(self) -> None:
        solntype = self.solutiontype()
        code = self.code()
        logging.debug(
            f"CORS_Timeseries: Loading SQLite data from {self._dbfile} for {code}:{solntype}"
        )
        db = SqliteTimeseries._openDb(self._dbfile)
        params: list = []
        if "+" not in solntype:
            sql = SqliteTimeseries._sql
            params = [code, solntype]
        else:
            types = solntype.split("+")
            casesql = "CASE solution_type"
            for i, t in enumerate(types):
                casesql = casesql + " WHEN ? THEN " + str(i)
            casesql = casesql + " ELSE " + str(len(types)) + " END"
            sql = SqliteTimeseries._sqlMultiType
            sql = sql.replace("{case}", casesql)
            params = []
            params.extend(types)
            params.append(code)
            params.extend(types)
            params.append(len(types))
            params.append(code)
            params.extend(types)

        when = ""
        if self._before is not None:
            when = " and m.epoch < ?"
            params.append(self._before.strftime("%Y-%m-%d"))
        if self._after is not None:
            when = when + " and m.epoch > ?"
            params.append(self._after.strftime("%Y-%m-%d"))
        sql = sql.replace("{when}", when)
        logging.debug(f"SQL: {sql}")
        logging.debug(f"Params: {params}")
        data = pd.read_sql(sql, db, params=params, index_col="epoch")
        db.close()
        # , parse_dates=['epoch'])
        data.set_index(pd.to_datetime(data.index), inplace=True)
        return data


class PgTimeseries(Timeseries):
    """
    Postgres time series.  Assumes data are held in a table coords with columns
    solution_id (identifies calc date SINEX file from which solutions loaded)
    solution_type (string defining solution type)
    code (mark id)
    epoch (date/time of calculated coordinate)
    x,y,z (ordinates)

    Assumes connection details are defined in environment variables.
    """

    GetCoordinatesSql = """
        with p as (
            select 
               %s as solution_type,
               %s as code,
               %s::timestamp without time zone as min_epoch,
               %s::timestamp without time zone as max_epoch
        )
        select DISTINCT ON (date(epoch))
           to_char(epoch,'YYYY-MM-DD') as epoch, 
           X as x, Y as y, Z as z 
        from coords c, p
        where 
           c.code=p.code and 
           c.solution_type=p.solution_type and
           (p.min_epoch IS NULL OR epoch >= p.min_epoch) and
           (p.max_epoch IS NULL OR epoch <= p.max_epoch)
        order by date(epoch),solution_id DESC
        """

    GetMultiSolutionCoordinatesSql = """
            with p as (
            select 
               %s as solution_type,
               %s as code,
               %s::timestamp without time zone as min_epoch,
               %s::timestamp without time zone as max_epoch
        ),
        stypes as (
            select
                s as solution_type,
                array_position(string_to_array((select solution_type from p),'+'),s) as priority
            from
                unnest(string_to_array((select solution_type from p),'+')) s
        ),
        solns as (
            select 
                epoch, 
                c.solution_type, 
                solution_id,
                X as x, Y as y, Z as z, 
                stypes.priority
            from p, coords c join stypes on c.solution_type=stypes.solution_type
            where 
            c.code=p.code and 
           (p.min_epoch IS NULL OR epoch >= p.min_epoch) and
           (p.max_epoch IS NULL OR epoch <= p.max_epoch)
        )
        select 
            distinct on (date(epoch))
            to_char(epoch,'YYYY-MM-DD') as epoch,
            solution_type,
            x, y, z
        from
            solns
        where 
            priority is not NULL
        order by
            date(epoch),priority,solution_id DESC
        """

    GetCodeSolutionTypeSql = """
        select distinct code, solution_type
        from coords
        order by code, solution_type
        """

    @staticmethod
    def _openDb(source: str | None) -> pd._SqlConnection:
        import psycopg2

        source = source or ""
        srcprm: dict = {}
        if source.startswith("pg:"):
            source = source[3:]
        if source == "":
            dbname = os.environ.get("PGDATABASE")
            if dbname is None:
                raise RuntimeError(
                    "Postgres database not defined by environment variables"
                )
            source = "dbname=" + dbname
        try:
            if source.startswith("{"):
                srcprm = json.loads(source)
            else:
                if "=" not in source:
                    source = "dbname=" + source
                srcprm = {"dsn": source}
            db = psycopg2.connect(**srcprm)
        except:
            if "password" in srcprm:
                srcprm["password"] = "..."
            raise RuntimeError(f"Cannot connect to database {srcprm}")
        return db

    @staticmethod
    def seriesList(
        source: str | None = None,
        solutiontype: str | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
    ) -> list[PgTimeseries]:
        db = PgTimeseries._openDb(source)
        seriescodes = pd.read_sql(PgTimeseries.GetCodeSolutionTypeSql, db)
        db.close()
        series = []
        found = []
        stypes = solutiontype.split("+") if solutiontype else None
        for i in seriescodes.index:
            code, solntype = (seriescodes.code[i], seriescodes.solution_type[i])
            if stypes is None or solntype in stypes:
                solntype = solutiontype or solntype
                if (code, solntype) in found:
                    continue
                series.append(
                    PgTimeseries(
                        source,
                        code,
                        solntype,
                        after=after,
                        before=before,
                        normalize=normalize,
                    )
                )
                found.append((code, solntype))
        return series

    def __init__(
        self,
        source: str | None,
        code: str,
        solutiontype: str = "default",
        xyz0: np.ndarray | None = None,
        transform: TransformFunc | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
    ):
        Timeseries.__init__(
            self,
            code,
            solutiontype=solutiontype,
            xyz0=xyz0,
            transform=transform,
            after=after,
            before=before,
            normalize=normalize,
        )
        self._source = source

    def _loadData(self) -> None:
        solntype = self.solutiontype()
        code = self.code()
        logging.debug(
            f"CORS_Timeseries: Loading Postgres data from {self._source} for {code}:{solntype}"
        )
        db = self._openDb(self._source)
        if "+" not in solntype:
            sql = PgTimeseries.GetCoordinatesSql
        else:
            sql = PgTimeseries.GetMultiSolutionCoordinatesSql

        minepoch = None
        maxepoch = None
        if self._before is not None:
            maxepoch = self._before.strftime("%Y-%m-%d")
        if self._after is not None:
            minepoch = self._after.strftime("%Y-%m-%d")
        params = [solntype, code, minepoch, maxepoch]
        data = pd.read_sql(sql, db, params=params, index_col="epoch")
        db.close()
        data.set_index(pd.to_datetime(data.index), inplace=True)
        return data


class CoordApiTimeseries(Timeseries):
    """
    Timeseries provide by coordinate API from LINZ AWS linz-gnss-coord-analysis template/stack.

    The API is defined by a url which gives two endpoints.  The solutions_stats endpoint is used
    to get a list of strategies and codes available from the coordinates endpoint.  (These APIs
    match those provided by the GNSS processing daily processing on AWS).

      {baseurl}/solution_stats [strategy=xxx[+xxx]*] [code=xxxx[+xxxx]*]
      {baseurl}/coordinates?code=xxxx&strategy=xxx[+xxx]* [from_epoch=YYYY:DDD] [to_epoch=YYYY:DDD]
    """

    DefaultApiUrl = "http://results.prod.gnss-processing.awsint.linz.govt.nz/api/v1"
    ApiUrlEnvVar = "GNSS_COORD_API_URL"

    @staticmethod
    def _getApiUrl(apiurl: str | None) -> str:
        if apiurl is None or apiurl == "":
            url = os.environ.get(
                CoordApiTimeseries.ApiUrlEnvVar, CoordApiTimeseries.DefaultApiUrl
            )
            return url
        return apiurl

    @staticmethod
    def _getData(source, apiquery: str) -> Any:
        import requests

        if source.startswith("api:"):
            source = source[4:]
        url = source + "/" + apiquery
        try:
            logging.debug(f"CoordApiTimeseries: Requesting {url}")
            with requests.get(url) as response:
                if response.status_code != 200:
                    raise RuntimeError(f"Status code {response.status_code}")
                data = response.json()
            return data
        except Exception as ex:
            raise RuntimeError(f"Failed to load {url}: {ex}")

    @staticmethod
    def seriesList(
        solutiontype: str | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
        apiurl: str | None = None,
    ) -> list[CoordApiTimeseries]:
        apiurl = CoordApiTimeseries._getApiUrl(apiurl)
        series: list[CoordApiTimeseries] = []
        query = "solution_stats"
        if solutiontype:
            query = query + f"?strategy={solutiontype}"
        for strategy in CoordApiTimeseries._getData(apiurl, query)["strategies"]:
            for code in strategy["codes"]:
                series.append(
                    CoordApiTimeseries(
                        code["code"],
                        strategy["strategy"],
                        after=after,
                        before=before,
                        normalize=normalize,
                        apiurl=apiurl,
                    )
                )
        return series

    def __init__(
        self,
        code: str,
        solutiontype: str = "default",
        xyz0: np.ndarray | None = None,
        transform: TransformFunc | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
        apiurl: str | None = None,
    ):
        Timeseries.__init__(
            self,
            code,
            solutiontype=solutiontype,
            xyz0=xyz0,
            transform=transform,
            after=after,
            before=before,
            normalize=normalize,
        )
        self._apiurl = self._getApiUrl(apiurl)

    def _loadData(self) -> None:
        solntype = self.solutiontype()
        code = self.code()
        logging.debug(
            f"CORS_Timeseries: Loading Coordinate API data from {self._apiurl} for {code}:{solntype}"
        )
        query = f"coordinates?code={code}&strategy={solntype}"
        if self._after is not None:
            query = query + "&from_epoch=" + self._after.strftime("%Y:%j")
        if self._before is not None:
            query = query + "&to_epoch=" + self._before.strftime("%Y:%j")
        result = self._getData(self._apiurl, query)
        tsdata = result["data"]
        columns = ["solution_type" if c == "strategy" else c for c in result["fields"]]
        data = pd.DataFrame(tsdata, columns=columns)
        data.set_index(pd.to_datetime(data["sessionid"], format="%Y:%j"), inplace=True)
        return data


class FileTimeseries(Timeseries):

    """
    Defines a CORS time series from a file for analysis.

    Assumes the time series file has columns name, epoch, x, y, z and is
    a comma delimited CSV file, or if delim_whitespace is true then a whitespace delimited file.

    The xyz0 parameter can define a reference point used for calculating e, n, u
    components.  The default is to use the first coordinate.

    The transform argument can define a function which is applied to each coordinate
    as it is read.  The transform function is called for each coordinate as

    xyz=transform(xyz,epoch)

    where epoch is a datetime.datetime object

    Flags observations as excluded if they match the dates stored with the model.
    """

    _columns = ("name", "epoch", "x", "y", "z")
    _epoch_col = "epoch"
    _filename: str

    @staticmethod
    def seriesList(
        filepattern: str,
        solutiontype: str = "default",
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
        delim_whitespace: bool = False,
    ):
        """
        Get the potential time series files, any file matching the filepattern.  The
        pattern can include {code} in the filename to represent the code to use
        """
        path, file = os.path.split(filepattern)
        fileparts = file.split("{code}")
        if len(fileparts) != 2:
            raise RuntimeError(
                "File timeseries file pattern must include {code}: " + filepattern
            )
        filere = re.compile(
            re.escape(fileparts[0])
            + r"(?P<code>\w\w\w\w)"
            + re.escape(fileparts[1])
            + "$"
        )
        path = path or "."
        series = []
        for fn in os.listdir(path or "."):
            match = filere.match(fn)
            if not match:
                continue
            if path:
                fn = os.path.join(path, fn)
            code = match.group("code")
            series.append(
                FileTimeseries(
                    fn,
                    code,
                    solutiontype,
                    after=after,
                    before=before,
                    normalize=normalize,
                    delim_whitespace=delim_whitespace,
                )
            )
        return series

    def __init__(
        self,
        filename: str,
        code: str | None = None,
        solutiontype: str = "default",
        xyz0: np.ndarray | None = None,
        transform: TransformFunc | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
        delim_whitespace: bool = False,
    ):
        if code is None:
            code = ""
        Timeseries.__init__(
            self,
            code,
            solutiontype=solutiontype,
            xyz0=xyz0,
            transform=transform,
            after=after,
            before=before,
            normalize=normalize,
        )
        self._filename = filename
        self._delim_whitespace = delim_whitespace

    def _loadData(self) -> None:

        filename = self._filename
        code = self._code
        logging.debug(f"FileTimeseries: Loading data from {filename} for {code}")
        if not os.path.exists(filename):
            raise RuntimeError(
                "CORS time series file {0} does not exist".format(filename)
            )

        try:
            data = pd.read_csv(
                filename,
                delim_whitespace=self._delim_whitespace,
                parse_dates=True,
                index_col=self._epoch_col,
            )
        except:
            raise RuntimeError(
                "Cannot read time series data from {0}: {1}".format(
                    filename, sys.exc_info()[1]
                )
            )

        columns = data.columns
        missing_cols = [
            x for x in self._columns if x not in columns and x != self._epoch_col
        ]
        if missing_cols:
            raise RuntimeError(
                "Columns {0} missing from CORS time series file {1}".format(
                    " ".join(missing_cols), filename
                )
            )

        data.rename(columns={"name": "code"}, inplace=True)

        code = code if code else data.code.iloc[0]
        self.setName(code)
        data = data[data.code == code]
        return data

    def filename(self) -> str:
        return self._filename


TimeseriesCalcFunc = Callable[[dt.datetime], list[float] | np.ndarray]


class FunctionTimeseries(Timeseries):

    """
    Creates a CORS timeseries based on a function. The function must take a date
    as an input and return an XYZ coordinate.
    """

    # This defines _function as a method.  How to define a callable non-method function?
    # _function: TimeseriesCalcFunc

    def __init__(
        self,
        function: TimeseriesCalcFunc,
        code: str = "Function",
        solutiontype: str = "function",
        xyz0: np.ndarray | None = None,
        transform: TransformFunc | None = None,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        increment: int | float = 1,
        index: pd.Index | None = None,
        fillDays: bool = False,
    ):
        """
        Calculates a time series either each day from after to before, or for the days specified by index.
        """
        Timeseries.__init__(
            self,
            code,
            solutiontype=solutiontype,
            xyz0=xyz0,
            transform=transform,
            after=after,
            before=before,
        )
        if not callable(function):
            raise RuntimeError("FunctionTimeseries needs a callable function")
        self._function = function
        self._isfunction = True
        self.setDates(after=after, before=before, index=index, fillDays=fillDays)

    def setDates(
        self,
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        increment: int | float = 1,
        index: pd.Index | None = None,
        fillDays: bool = False,
        dates: list[dt.datetime] | None = None,
    ):
        """
        Set the dates at which the function timeseries will be calculated.
        Can initiallize with start and end dates (after=, before=), and existing pd.DateTimeIndex
        (index=..) optionally with fillDays to fill missing days in the time series, or with an array
        of dates (dates=[...]).
        """

        self.setDateRange(reset=True, after=after, before=before)
        if dates is not None:
            index = pd.DatetimeIndex(data=dates)

        if index is not None and fillDays:
            after = index[0]
            before = index[-1]
            index = None
        if index is None:
            if after is None:
                dtfrom = dt.date(2000, 1, 1)
            else:
                dtfrom = dt.date(after.year, after.month, after.day)
            if before is None:
                dtto = dt.date.today()
            else:
                dtto = dt.date(before.year, before.month, before.day)
            incdef = str(int(increment)) + "D"
            index = pd.DatetimeIndex(start=dtfrom, end=dtto, freq=incdef)
        self._index = index
        self._loaded = False

    def _loadData(self) -> None:
        index = self._index
        function = self._function
        xyzdata = [function(d) for d in index.to_pydatetime()]
        data = pd.DataFrame(data=xyzdata, index=index.copy(), columns=("x", "y", "z"))
        return data


class TimeseriesList(object):
    """
    Defines a source of coordinate time series for different codes.  Once the series is created then
    time series for specific codes can be retrieved from it.  Using a string source to defined the input.
    The source can be one of:

        sqlite:filename
        pg:connection_str
        file: filepath
        api:url

    (A file: source filepath must include "{code}" that is the part of the filename that is the mark code)

    """

    def __init__(
        self,
        source=None,
        solutiontype=None,
        after=None,
        before=None,
        normalize=False,
        series=None,
    ):
        self._get = None
        eparams = {}
        if after is not None:
            eparams["after"] = after
        if before is not None:
            eparams["before"] = before
        if source is not None and series is None:
            series = []
            if source.startswith("file:"):
                filesource = source[5:]
                if "{code}" not in filesource:
                    raise RuntimeError(
                        f'TimeseriesList: Specification of file time series {source} must include "{{code}}"'
                    )
                series.extend(
                    FileTimeseries.seriesList(
                        filesource,
                        solutiontype,
                        after=after,
                        before=before,
                        normalize=normalize,
                    )
                )
                self._get = lambda code, **params: FileTimeseries(
                    filename=source.replace("{code}", code),
                    code=code,
                    normalize=normalize,
                    **eparams,
                    **params,
                )
            if source.startswith("pg:"):
                series.extend(
                    PgTimeseries.seriesList(
                        source[3:],
                        solutiontype,
                        after=after,
                        before=before,
                        normalize=normalize,
                    )
                )
                self._get = lambda code, **params: PgTimeseries(
                    source=source, code=code, normalize=normalize, **eparams, **params
                )
            elif source.startswith("api:"):
                series.extend(
                    CoordApiTimeseries.seriesList(
                        solutiontype,
                        after=after,
                        before=before,
                        normalize=normalize,
                        apiurl=source[4:],
                    )
                )
                self._get = lambda code, **params: CoordApiTimeseries(
                    source=source, code=code, normalize=normalize, **eparams, **params
                )
            elif source.startswith("sqlite:"):
                series.extend(
                    SqliteTimeseries.seriesList(
                        source[7:],
                        solutiontype,
                        after=after,
                        before=before,
                        normalize=normalize,
                    )
                )
                self._get = lambda code, **params: SqliteTimeseries(
                    dbfile=source, code=code, normalize=normalize, **eparams, **params
                )
        self._series = series

    def codes(self):
        codes = {}
        for f in self._series:
            codes[f.code()] = 1
        return sorted(codes.keys())

    def solutiontypes(self):
        solutiontypes = set()
        for f in self._series:
            for s in f.solutiontype().split("+"):
                solutiontypes.add(s)
        return sorted(solutiontypes)

    def get(self, code, solutiontype=[], after=None, before=None):
        """
        Return the time series for the requested station codes.  Can
        specify a solution type wanted, or a list of solution types in
        order of preference.  If solutiontype is not specified then
        an error will be raised if there is more than one matching solution.

        The time series can be restricted to a date range by using before and after
        parameters.  These can be python datetime.datetime objects, or strings formatted
        as yyyy-mm-dd.
        """
        if isinstance(solutiontype, str):
            solutiontype = solutiontype.split("+")
        series_types = []
        baseseries = None
        for series in self._series:
            if series.code() != code:
                continue
            baseseries = series
            series_types.append(series.solutiontype())
        if len(solutiontype) > 0:
            series_types = [s for s in series_types if s in solutiontype]
        elif len(series_types) > 1:
            raise RuntimeError(
                "Ambiguous time series requested - multiple possible solution types"
            )
        if len(series_types) == 0:
            raise RuntimeError(
                "Invalid code/solution type requested - no matching time series found"
            )
        if self._get is not None and len(solutiontype) > 1:
            return self._get(
                code, solutiontype=series_types, after=after, before=before
            )
        return baseseries

    def compareWith(self, other, after=None, before=None, normalize=False):
        """
        Return a TimeseriesList.Comparison object which compares this list with other (results are in the
        sense self-other)
        """
        return TimeseriesList.Comparison(
            self, other, after=after, before=before, normalize=normalize
        )

    class Comparison(object):
        """
        Class representing a comparison between two time series lists.
        """

        ts1: TimeseriesList
        ts2: TimeseriesList
        before: dt.datetime | None
        after: dt.datetime | None
        normalize: bool

        def __init__(
            self,
            ts1: TimeseriesList,
            ts2: TimeseriesList,
            after: dt.datetime | None = None,
            before: dt.datetime | None = None,
            normalize: bool = False,
        ):
            self.ts1 = ts1
            self.ts2 = ts2
            codes1 = self.ts1.codes()
            codes2 = self.ts2.codes()
            self.before = before
            self.after = after
            self.normalize = normalize
            self.codes = [c for c in codes1 if c in codes2]
            self.codes.sort()

        def get(
            self,
            code: str,
            after: dt.datetime | None = None,
            before: dt.datetime | None = None,
            plot: bool = False,
            plotDiff: bool = False,
            stats: bool = False,
            detrend: bool = True,
        ) -> Timeseries.Comparison:
            """
            Returns the difference between two time series for a code.  Returns a Timeseries.Comparison object
            """

            after = after or self.after
            before = before or self.before
            ts1 = self.ts1.get(code, after=after, before=before)
            ts2 = self.ts2.get(code, after=after, before=before)
            comparison = Timeseries.Comparison(ts1, ts2, normalize=self.normalize)
            if plot:
                comparison.plot(detrend=detrend)
            if plotDiff:
                comparison.plotDiff()
            if stats:
                print(comparison.stats())
            return comparison

        def stats(self) -> pd.DataFrame:
            """
            Get summary statistics for common stations between the two time series.
            Returns a DataFrame indexed on code, with columns rse_sol1_[enu],
            rse_sol2_[enu] where sol1 and sol2 are the solution type,s and then
            diff_e, std_diff_e, diff_n, std_diff_n, diff_u, std_diff_u, where diff
            and

            """
            data = []
            stype1 = "ts1"
            stype2 = "ts2"
            for code in self.codes:
                cmp = self.get(code)
                stype1 = cmp.ts1.solutiontype()
                stype2 = cmp.ts2.solutiontype()
                row = [code]
                for ts in (cmp.ts1, cmp.ts2):
                    row.extend(ts.robustStandardError())
                diff = cmp.diff.getData()
                stats = diff.describe()
                stats2 = np.abs(diff).describe()
                row.extend(
                    (
                        stats.loc["mean", "e"],
                        stats.loc["std", "e"],
                        stats2.loc["max", "e"],
                        stats.loc["mean", "n"],
                        stats.loc["std", "n"],
                        stats2.loc["max", "n"],
                        stats.loc["mean", "u"],
                        stats.loc["std", "u"],
                        stats2.loc["max", "u"],
                    )
                )
                data.append(row)

            columns = ["code"]
            columns.extend((stype1 + "_rse_e", stype1 + "_rse_n", stype1 + "_rse_u"))
            columns.extend((stype2 + "_rse_e", stype2 + "_rse_n", stype2 + "_rse_u"))
            columns.extend(
                (
                    "diff_mean_e",
                    "diff_sdt_e",
                    "diff_maxabs_e",
                    "diff_mean_n",
                    "diff_std_n",
                    "diff_maxabs_n",
                    "diff_mean_u",
                    "diff_std_u",
                    "diff_maxabs_u",
                )
            )
            result = pd.DataFrame(data, columns=columns)
            result.set_index(result.code, inplace=True)
            result.dropna(inplace=True)
            return result


class SqliteTimeseriesList(TimeseriesList):
    def __init__(
        self,
        dbfile: str,
        solutiontype: str = "default",
        after: dt.datetime | None = None,
        before: dt.datetime | None = None,
        normalize: bool = False,
    ):
        if dbfile.startswith("sqlite:"):
            dbfile = dbfile[7:]
        if not os.path.isfile(dbfile):
            raise RuntimeError(
                "SqliteTimeseriesList source is not a file: " + str(dbfile)
            )
        self._dbfile = dbfile
        self._solutiontype = solutiontype
        TimeseriesList.__init__(
            self,
            source=f"sqlite:{dbfile}",
            solutiontype=solutiontype,
            after=after,
            before=before,
            normalize=normalize,
        )

    def solutionTypes(self):
        _sqlTypes = """
            select distinct solution_type from mark_coordinate m
            """
        db = SqliteTimeseries._openDb(self._dbfile)
        types = pd.read_sql(SqliteTimeseries._sqlTypes, db)
        db.close()
        return [types.solution_type[i] for i in types.index]

    def stationCoordinates(
        self,
        date: dt.datetime,
        solutiontype: str | None = None,
        codes: list[str] | None = None,
        selectday: bool = False,
    ):
        """
        Return a DataFrame of stations and coordinates that apply at a specific
        epoch (selectdate=False), or during (selectday=True) a (UTC) day.
        """

        sql = """
            select code as code, solution_type as solutiontype, epoch, X as x, Y as y, Z as z 
            from mark_coordinate m
            where solution_type=?
            {when}
            {codelist}
            order by code
            """
        db = SqliteTimeseries._openDb(self._dbfile)

        solutiontype = solutiontype or self._solutiontype
        if solutiontype is None:
            types = self.solutionTypes()
            if len(types) == 1:
                solutiontype = type
            else:
                raise RuntimeError(
                    "Ambiguous solution types in SqliteTimeseriesList.stationCoordinates"
                )
        elif "+" in solutiontype:
            raise RuntimeError(
                "Cannot use SqliteTimeseriesList.stationCoordinates with multiple solution types"
            )

        params = [solutiontype]
        when = ""
        dateformat = "%Y-%m-%d %H:%M:%S"
        if selectday:
            when = "and epoch between ? and ?"
            params.extend(
                [
                    dt.datetime(date.year, date.month, date.day, 0, 0, 0).strftime(
                        dateformat
                    ),
                    dt.datetime(date.year, date.month, date.day, 23, 59, 59).strftime(
                        dateformat
                    ),
                ]
            )
        else:
            when = "and epoch=?"
            params.append(date.strftime(dateformat))
        codelist = ""
        if codes is not None:
            codelist = (
                "and code in (" + (",".join(["'" + c + "'" for c in codes])) + ")"
            )
        sql = sql.replace("{when}", when)
        sql = sql.replace("{codelist}", codelist)
        data = pd.read_sql(sql, db, params=params, index_col="code")
        db.close()
        return data
