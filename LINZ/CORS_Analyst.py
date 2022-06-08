import sys
import os.path
import shutil
import urllib.request
import json
import re
import tempfile
import zipfile
import numpy as np
import pandas as pd
import yaml
from datetime import datetime
from collections import namedtuple
from pandas.tseries.offsets import DateOffset

from LINZ.Geodetic import GDB
from .StationCoordinateModel import Model as StationCoordinateModel

from . import CORS_Timeseries
from . import GDB_Timeseries

StationData = namedtuple(
    "StationData", "code timeseries stationCoordModel scmTimeseries gdbTimeseries"
)
StationOffsetStats = namedtuple("StationOffsetStats", "date gdbOffset scmOffset")
RefCoordinateDate = datetime(2018, 1, 1)


class CORS_Analyst(object):

    refconfigfile = os.path.splitext(__file__)[0] + ".yaml"
    configfile = os.path.basename(refconfigfile)
    httpre = re.compile(r"^https?\:\/\/")

    def __init__(
        self,
        verbose=False,
        write_plot_files=False,
        configfile=None,
        params=None,
    ):
        self._verbose = verbose
        self._writecsv = write_plot_files
        self._configfile = configfile or CORS_Analyst.configfile
        self._tempdirs = []

        if not os.path.exists(self._configfile):
            sys.stderr.write(
                "Cannot find configuration file {0}\n".format(self._configfile)
            )
            sys.exit()

        if self._verbose:
            print("Reading configuration from", self._configfile)

        cfg = {"test_stations": "", "test_stations_file": ""}
        filecfg = yaml.load(open(self._configfile), Loader=yaml.SafeLoader)
        cfg.update(filecfg)
        if params:
            cfg.update(params)
        self._cfg = cfg
        self._tsdb = self.getCfg("timeseries_database")
        self._tssolutiontype = self.getCfg("timeseries_solution_type")
        self._tslist = CORS_Timeseries.TimeseriesList(
            source=self._tsdb, solutiontype=self._tssolutiontype
        )
        self._scmpath = self.getCfg("station_coordinate_model_path")
        if self.httpre.match(self._scmpath):
            zipspec, scmfile = self._scmpath.rsplit("/", maxsplit=1)
            self._scmpath = self.downloadZipFile(zipspec) + "/" + scmfile
        self._offsetdays = self.getCfg("offset_calc_ndays", parser=int)
        self._teststat = self.getCfg("offset_test_stat")
        self._gdbwarning = self.loadWarningLevels("gdb_offset_warning")
        self._scmwarning = self.loadWarningLevels("scm_offset_warning")

        self.loadDeformationModel()

    def getCfg(self, cfgitem, default=None, parser=None):
        if cfgitem in self._cfg:
            value = self._cfg[cfgitem]
            if parser:
                try:
                    value = parser(value)
                except:
                    raise RuntimeError(
                        f"Invalid value {value} for configuration item {cfgitem}"
                    )
            return value
        elif default is not None:
            return default
        else:
            raise RuntimeError(f"Missing configuration item {cfgitem}")

    def loadWarningLevels(self, warning_type):
        try:
            levels = [float(x) for x in str(self.getCfg(warning_type)).split()]
            if len(levels) < 2:
                levels.append[levels[0]]
        except:
            raise RuntimeError(
                "Invalid value for " + warning_type + " configuration item"
            )
        return levels

    def downloadZipFile(self, zipfileurl):

        try:
            datapath = None
            if "+" in zipfileurl:
                zipfileurl, datapath = zipfileurl.split("+", maxsplit=1)
            tempdir = tempfile.TemporaryDirectory()
            data = urllib.request.urlopen(zipfileurl)
            dirname = tempdir.name
            zipfilename = os.path.join(dirname, "download.zip")
            with open(zipfilename, "wb") as zfh:
                shutil.copyfileobj(data, zfh)
            with zipfile.ZipFile(zipfilename, "r") as zf:
                zf.extractall(dirname)
            if datapath is not None:
                dirname = os.path.join(dirname, datapath)
            if os.path.isdir(dirname):
                # Ensure temp directory remains for duration of script
                self._tempdirs.append(tempdir)
                return dirname
            raise RuntimeError(
                "Directory {0} not found in zipfile {1}".format(datapath, zipfileurl)
            )
        except Exception as e:
            raise RuntimeError("Failed to download {0}: {1}".format(zipfileurl, e))

    def loadDeformationModel(self):
        path = self.getCfg("deformation_model_path")
        if self.httpre.match(path):
            path = self.downloadZipFile(path)
        self.gdbCalculator = GDB_Timeseries.GDB_Timeseries_Calculator(path)
        self.deformationModelVersion = self.gdbCalculator.deformationModelVersion()

    def gdbTimeseries(self, code, dates=None, xyz0=None, xyz0Date=None):
        if not GDB.get(code):
            raise RuntimeError("Cannot get GDB data for {0}".format(code))
        gdbts = self.gdbCalculator.get(code, xyz0=xyz0, index=dates, xyz0Date=xyz0Date)
        return gdbts

    def scmTimeseries(self, model, dates, xyz0=None):
        if model is None:
            return None
        spmxyz = model.calc(dates.to_pydatetime(), enu=False)
        code = model.station
        scmts = CORS_Timeseries.Timeseries(
            code, solutiontype="scm", dates=dates, xyz=spmxyz, xyz0=xyz0
        )
        return scmts

    def getStationData(self, code):
        # Load the time series

        ts = self._tslist.get(code)

        # Define a reference coordinate for the station to apply to the time series
        # (so that we get compatible xyz offsets)

        dates = ts.dates()
        if len(dates) == 0:
            raise RuntimeError("No timeseries data found for code {0}".format(code))

        # Create the GDB coordinate time series, and use it to set the reference
        # coordinate.

        gdbts = self.gdbTimeseries(code, dates, xyz0Date=RefCoordinateDate)
        xyz0 = gdbts.xyz0()
        ts.setXyzTransform(xyz0=xyz0)

        # Set up a time series of based on the station coordinate model
        spmts = None
        model = None
        if self._scmpath != "none":
            # Load the station coordinate model corresponding to the station
            scmfile = self._scmpath.replace("{code}", code)
            model = StationCoordinateModel(filename=scmfile)
            spmts = self.scmTimeseries(model, dates, xyz0)

        return StationData(code, ts, model, spmts, gdbts)

    def getStationOffsets(self, stnData, ndays):
        # Calculate mean offsets the last nDays data

        dates = stnData.timeseries.dates()

        refdate = dates[-1] + DateOffset(days=-ndays)
        calcdate = refdate + DateOffset(days=ndays / 2)
        select = dates > refdate
        obsts = stnData.timeseries.getData()[select]
        gdbts = stnData.gdbTimeseries.getData()[select]
        scmdiff = None
        if stnData.scmTimeseries is not None:
            scmts = stnData.scmTimeseries.getData()[select]
            scmdiff = (obsts - scmts).describe()
        return StationOffsetStats(calcdate, (obsts - gdbts).describe(), scmdiff)

    def writeTimeseriesCSV(self, code, tsdata, csv_cfg_type):
        results_dir = self.getCfg("result_path")
        results_file = self.getCfg(csv_cfg_type)
        csv_file = os.path.join(results_dir, results_file.replace("{code}", code))
        self.createFileDirectory(csv_file)
        tsdata.to_csv(csv_file, float_format="%.4f", index_label="epoch")

    def roundFloatsInStruct(self, ndp):
        def rounder(o):
            if isinstance(o, float):
                return round(o, ndp)
            if isinstance(o, dict):
                return {k: rounder(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [rounder(x) for x in o]
            return o

        return rounder

    def writeResultsJSON(self, results):
        results_dir = self.getCfg("result_path")
        results_file = self.getCfg("summary_file")
        json_file = os.path.join(results_dir, results_file)
        rounder = self.roundFloatsInStruct(4)
        self.createFileDirectory(json_file)
        with open(json_file, "w") as jf:
            jf.write(json.dumps(rounder(results), indent=2, sort_keys=True))

    def createFileDirectory(self, filename):
        dirname = os.path.dirname(filename)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

    def trimmedModelTimeseries(self, tsdata, precision):
        # Remove redundant data for plotting model time series based
        # Removes data points which are on line between included points
        # to a given precision.
        if tsdata is None:
            return None
        enu = np.vstack((tsdata.e, tsdata.n, tsdata.u)).T
        timevar = tsdata.index.to_pydatetime()
        reftime = timevar[0]
        timevar = np.array(
            [(t - reftime).total_seconds() / (60 * 60 * 24.0) for t in timevar]
        )
        usepoint = [False] * len(timevar)
        precision *= precision

        def addPoints(start, end, usepoint):
            usepoint[start] = True
            usepoint[end] = True
            if end <= start + 1:
                return

            timediff = timevar[end] - timevar[start]
            diff = (enu[end] - enu[start]).reshape(1, 3)
            offset = timevar[start : end + 1].reshape(end + 1 - start, 1).dot(diff)
            offset = enu[start : end + 1] - enu[start, :] - offset / timediff
            offset = np.square(offset).sum(axis=1)
            i = np.argmax(offset)
            if i > 0 and i < end - start and offset[i] > precision:
                addPoints(start, start + i, usepoint)
                addPoints(start + i, end, usepoint)

        addPoints(0, len(usepoint) - 1, usepoint)
        return tsdata[usepoint]

    def compileOffsetStats(self, code, stats, teststat, levels, testtype):
        results = {}
        fixnan = lambda f: 0.0 if np.isnan(f) else f
        for m in "count mean std median min max".split():
            mloc = "50%" if m == "median" else m
            mstats = stats.loc[mloc]
            results[m] = [fixnan(mstats.e), fixnan(mstats.n), fixnan(mstats.u)]
        results["count"] = results["count"][0]
        teststats = results[teststat]
        testh = np.hypot(teststats[0], teststats[1])
        testv = abs(teststats[2])
        status = "warn" if testh > levels[0] or testv > levels[1] else "good"
        message = None
        if status == "warn":
            message = ""
            if testh > levels[0]:
                message = (
                    message
                    + "{0} {1} {2} horizontal offset {3:.4f} m exceeds tolerance {4:.4f} m\n".format(
                        code, testtype, teststat, testh, levels[0]
                    )
                )
            if testv > levels[1]:
                message = (
                    message
                    + "{0} {1} {2} vertical offset {3:.4f} m exceeds tolerance {4:.4f} m \n".format(
                        code, testtype, teststat, testv, levels[1]
                    )
                )
            if self._verbose:
                print(
                    "  {0} coordinates are significantly offset ({1:.4f} {2:.4f} m)".format(
                        testtype, testh, testv
                    )
                )

        results["status"] = status
        results["status_message"] = message
        return results

    def compileStationData(self, code):
        if self._verbose:
            print("Compiling data for", code)
        stndata = self.getStationData(code)
        offsets = self.getStationOffsets(stndata, self._offsetdays)
        stnresults = {
            "code": code,
            "offset_test_date": offsets.date.strftime("%Y-%m-%d"),
            "gdb_offset": self.compileOffsetStats(
                code,
                offsets.gdbOffset,
                self._teststat,
                self._gdbwarning,
                "GDB coordinate",
            ),
        }
        if stndata.stationCoordModel is not None:
            stnresults["scm_offset"] = self.compileOffsetStats(
                code,
                offsets.scmOffset,
                self._teststat,
                self._scmwarning,
                "station coordinate model",
            )
            stnresults["scm_version"] = stndata.stationCoordModel.versiondate.strftime(
                "%Y-%m-%d %H:%M:%S"
            )

        if self._writecsv:
            # Write the CSV file options

            fillmodel = self.getCfg("fill_model_timeseries").lower() == "yes"
            trimmodel = self.getCfg("trim_model_timeseries").lower() == "yes"
            combinets = self.getCfg("combine_timeseries").lower() == "yes"
            precision = self.getCfg("trim_model_precision", parser=float)

            # Generate data for plotting gdb coord and stn pred model
            # Get the observed timeseries
            ts = stndata.timeseries
            tsdata = ts.getData()

            # Fill the model
            if fillmodel:
                # If filling the model, normalize to whole dates
                tsdata.set_index(tsdata.index.normalize(), inplace=True)
                plotmargin = self.getCfg("plot_margin", parser=int)
                startdate = tsdata.index[0] - DateOffset(days=plotmargin)
                enddate = tsdata.index[-1] + DateOffset(days=plotmargin)
                dates = pd.date_range(startdate, enddate, freq="d")
                xyz0 = ts.xyz0()
                gdbmodel = self.gdbTimeseries(code, dates, xyz0).getData()
                scmmodel = self.scmTimeseries(stndata.stationCoordModel, dates, xyz0)
            else:
                gdbmodel = stndata.gdbTimeseries.getData()
                scmmodel = stndata.scmTimeseries
            if scmmodel is not None:
                scmmodel = scmmodel.getData()

            # Calculate the time series

            if trimmodel:
                gdbmodel = self.trimmedModelTimeseries(gdbmodel, precision)
                scmmodel = self.trimmedModelTimeseries(scmmodel, precision)

            if combinets:
                gdbmodel.columns = ("gdb_e", "gdb_n", "gdb_u")
                tsdata = tsdata.join(gdbmodel, how="outer")
                if scmmodel is not None:
                    scmmodel.columns = ("scm_e", "scm_n", "scm_u")
                    tsdata = tsdata.join(scmmodel, how="outer")
                self.writeTimeseriesCSV(code, tsdata, "timeseries_csv")
            else:
                self.writeTimeseriesCSV(code, tsdata, "timeseries_csv")
                self.writeTimeseriesCSV(code, gdbmodel, "gdb_timeseries_csv")
                if scmmodel is not None:
                    self.writeTimeseriesCSV(code, scmmodel, "scm_timeseries_csv")

        return stnresults

    def compileAllStationData(self, codelist=None):
        if codelist is None:
            codelist = []
        if len(codelist) == 0:
            codeliststr = self.getCfg("test_stations").strip()
            if codeliststr != "":
                codelist.extend(codeliststr.split())

            codelistfile = self.getCfg("test_stations_file").strip()
            if codelistfile != "":
                with open(codelistfile) as codef:
                    for l in codef:
                        l = l.strip()
                        if not l.startswith("#"):
                            codelist.extend(l.split())

        results = {
            "calcdate": datetime.now().strftime("%Y-%m-%d: %H:%M:%S"),
            "deformation_model_version": self.deformationModelVersion,
            "number_of_days_tested": self._offsetdays,
            "test_statistic": self._teststat,
        }
        coderesults = []
        for c in sorted(codelist):
            try:
                coderesults.append(self.compileStationData(c))
            except KeyboardInterrupt:
                break
            except:
                print("Station " + c + ": " + str(sys.exc_info()[1]))
        results["station_summary"] = coderesults
        self.writeResultsJSON(results)
        return results


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyse PositioNZ CORS data to monitor datum integrity. "
        + "Results of analysis are written to files specified in the configuration file",
        usage="%(prog)s [options] config_file code code ...",
    )
    parser.add_argument("config_file", help="Name of the configuration file")
    parser.add_argument(
        "code_or_param",
        nargs="*",
        help="Specific station codes to analyse.  Default is as defined in the configuration file.  Or param=value to replace config default param",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "-p",
        "--plot-files",
        action="store_true",
        help="Generate CSV files for plotting",
    )
    parser.add_argument(
        "-c",
        "--dump-config-file",
        action="store_true",
        help="Print a sample configuration file",
    )
    args = parser.parse_args()

    if args.dump_config_file:
        with open(CORS_Analyst.refconfigfile) as cf:
            for l in cf:
                print(l.strip())
        sys.exit()

    codes = []
    params = {}
    for c in args.code_or_param:
        if "=" in c:
            k, v = c.split("=", 1)
            params[k] = v
        else:
            codes.append(c.upper())

    a = CORS_Analyst(
        configfile=args.config_file,
        verbose=args.verbose,
        write_plot_files=args.plot_files,
        params=params,
    )
    a.compileAllStationData(codes)


if __name__ == "__main__":
    main()
