# Script to calculate an ITRF time series for a geodetic mark based on the
# official coordinate from the geodetic database and the NZGD2000 deformation
# model.

from datetime import datetime, timedelta
from collections import namedtuple

from LINZ.Geodetic.Ellipsoid import GRS80
from LINZ.Geodetic import GDB
from LINZ.DeformationModel.Model import Model as DeformationModel
from LINZ.DeformationModel.ITRF_NZGD2000 import Transformation as ITRF_Transformation
from .StationCoordinateModel import Model as StationCoordinateModel
from . import CORS_Timeseries


class GDB_Timeseries_Calculator(object):
    """
    Class to calculate the time series of a geodetic database mark based on the official 
    coordinate and the deformation model
    """

    StationData = namedtuple("StationData", "code xyz markdata timeseries")

    def __init__(
        self,
        deformationModelDirectory,
        itrf="ITRF2008",
        defaultDates=None,
        version=None,
    ):
        """
        Initiallize the time series calculator with the deformation model (defined
        by the directory in which it is located).  Defines the ITRF to calculate for
        (default ITRF2008) and the default dates for calculating the time series.
        """
        self._deformationModelDirectory = deformationModelDirectory
        self._deformationModelVersion = version
        self.itrf = itrf
        self._defaultDates = None
        self._itrfTransformation = None
        self.loadDeformationModel()

    def loadDeformationModel(self, deformationModelDirectory=None, version=None):
        if deformationModelDirectory is not None:
            self._deformationModelDirectory = deformationModelDirectory
        if version is not None:
            self._deformationModelVersion = version
        if self._deformationModelDirectory is not None:
            mod = DeformationModel(
                self._deformationModelDirectory, version=self._deformationModelVersion
            )
            self._deformationModelVersion = mod.version()
            self.setupItrf()

    def setupItrf(self, itrf=None):
        if itrf is not None:
            self.itrf = itrf
        if self._deformationModelDirectory is not None:
            self._itrfTransformation = ITRF_Transformation(
                modeldir=self._deformationModelDirectory,
                version=self._deformationModelVersion,
                toNZGD2000=False,
                itrf=self.itrf,
            )

    def deformationModelVersion(self):
        """
        Returns the deformation model version
        """
        return self._deformationModelVersion

    def itrf(self):
        """
        Returns the ITRF for which the time series are calculated
        """
        return self.itrf

    def getMarkData(
        self,
        code,
        index=None,
        fillDays=False,
        after=None,
        before=None,
        increment=1.0,
        xyz0=None,
        xyz0Date=None,
    ):
        """
        Get the time series for a geodetic mark.  Retrieves the official coordinates, calculates the
        deformation model at the required dates, and constructs a timeseries from them.

        Can supply a list of dates as an index, a from date, to date, and increment, 
        or use the default dates.

        The reference coordinate can be explicitly supplied, or calculated for a reference date. 
        Otherwise the NZGD2000 XYZ value will be used.
        """

        if self._deformationModelDirectory is None:
            raise RuntimeError(
                "GDB timeseries deformation moel transformation is not defined"
            )
        if self._itrfTransformation is None:
            raise RuntimeError("GDB timeseries ITRF transformation is not defined")

        markdata = GDB.get(code)
        if markdata is None:
            raise RuntimeError(
                "GDB timeseries for " + code + " not available - mark not in GDB"
            )
        coord = markdata.coordinate
        lon, lat, hgt = coord.longitude, coord.latitude, coord.height
        markxyz = GRS80.xyz(lon, lat, hgt)
        function = lambda d: GRS80.xyz(self._itrfTransformation(lon, lat, hgt, d))
        if xyz0 is None and xyz0Date is not None:
            xyz0 = function(xyz0Date)
        if xyz0 is None:
            xyz0 = markxyz
        gdbts = CORS_Timeseries.FunctionTimeseries(
            function,
            code=code,
            solutiontype="gdb",
            index=index,
            fillDays=fillDays,
            after=after,
            before=before,
            xyz0=xyz0,
        )
        return GDB_Timeseries_Calculator.StationData(code, markxyz, markdata, gdbts)

    def get(
        self,
        code,
        index=None,
        fillDays=False,
        after=datetime(1999, 12, 31),
        before=None,
        increment=1.0,
        xyz0=None,
        xyz0Date=None,
    ):
        data = self.getMarkData(
            code,
            index=index,
            fillDays=fillDays,
            before=before,
            after=after,
            increment=increment,
            xyz0=xyz0,
            xyz0Date=xyz0Date,
        )
        return data.timeseries
