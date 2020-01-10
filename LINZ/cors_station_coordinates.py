"""
Script to export cors station coordinates for selected stations and selected epoch
"""

import sys
import math
import re
import argparse
import datetime as dt
import csv
import os
import pandas as pd

from .CORS_Timeseries import TimeseriesList
from LINZ.DeformationModel.Model import deformationModelArguments, loadDeformationModel
from LINZ.DeformationModel.Error import OutOfRangeError
from LINZ.DeformationModel import ITRF_NZGD2000
from LINZ.Geodetic import ITRF
from LINZ.Geodetic.Ellipsoid import GRS80


def main():

    epilog = """
    Available coordinate systems are NZGD2000, ITRFxxxx.  These are geographic coordinates.
    For geocentric coordinates add _XYZ.  Multiple coordinate systems can be entered, eg
    ITRF2008+NZGD2000

    The epoch can be formatted as YYYY-MM-DD or now-#### for #### days before now.
    A range of dates can be selected as date1:date2
    Appended to this can be criteria :W# to select a day of the week (0-6, 0=Monday) or
    :M# for the day of month.
    

    The output CSV file can contain strings {yyyy} {mm} and {dd} that are replaced 
    the year, month and day.
    """
    parser = argparse.ArgumentParser(
        description="Export cors station coordinates from time series database",
        epilog=epilog,
        parents=[deformationModelArguments()],
    )
    parser.add_argument(
        "timeseries_data", help="Data source for CORS timeseries coordinates"
    )
    parser.add_argument(
        "epoch", help="Epoch or range at which coordinates are required"
    )
    parser.add_argument("output_csv_file", help="File to write coordinates to")
    parser.add_argument(
        "-c",
        "--coordinate-systems",
        help="Coordinate systems to output (default NZGD2000)",
    )
    parser.add_argument(
        "-i",
        "--cors-itrf",
        default="ITRF2008",
        help="CORS timeseries ITRF reference frame",
    )
    parser.add_argument(
        "-s",
        "--stations",
        help="Specifies stations to export, separated by +, or @filename",
    )
    parser.add_argument(
        "-n",
        "--nz-only",
        action="store_true",
        help="Only list New Zealand stations (including Chathams)",
    )
    parser.add_argument(
        "-t",
        "--solution-type",
        help="Solution type to extract (required if data includes more than one type)",
    )

    args = parser.parse_args()

    cors_itrf = args.cors_itrf
    tslist = None
    outputfiles = {}
    try:
        tslist = TimeseriesList(
            args.timeseries_data, solutiontype=args.solution_type, normalize=True
        )

        try:
            parts = args.epoch.split(":")
            if len(parts) < 2:
                parts.append(parts[0])
            startenddate = []
            for p in parts[:2]:
                match = re.match(r"^now-(\d+)$", p, re.I)
                if match:
                    dtp = dt.date.today() - dt.timedelta(days=int(match.group(1)))
                    startenddate.append(dt.datetime.combine(dtp, dt.time(0, 0, 0)))
                else:
                    startenddate.append(dt.datetime.strptime(args.epoch, "%Y-%m-%d"))
            useday = lambda d: True
            if len(parts) > 2:
                match = re.match(r"^([MW])(\d\d?)$", parts[2].upper())
                if not match:
                    raise RuntimeError("Invalid epoch date selector " + parts[2])
                dayno = int(match.group(2))
                if match.group(1) == "M":
                    useday = lambda d: d.day == dayno
                else:
                    useday = lambda d: d.weekday() == dayno
            startdate, enddate = startenddate
            increment = dt.timedelta(days=1)
            while startdate <= enddate:
                calcdate = startdate
                startdate = startdate + increment
                if not useday(calcdate):
                    continue
                filename = args.output_csv_file
                filename = filename.replace("{yyyy}", "{0:04d}".format(calcdate.year))
                filename = filename.replace("{mm}", "{0:02d}".format(calcdate.month))
                filename = filename.replace("{dd}", "{0:02d}".format(calcdate.day))
                filename = filename.replace(
                    "{ddd}", "{0:03d}".format(calcdate.timetuple().tm_yday)
                )
                if filename not in outputfiles:
                    outputfiles[filename] = []
                outputfiles[filename].append(calcdate)
        except:
            raise RuntimeError(
                "Invalid calculation epoch - must be YYYY-MM-DD:" + args.epoch
            )
        if len(outputfiles) == 0:
            raise RuntimeError("No dates defined for station coordinates")

        itrf_nzgd2000 = None
        nzgd2000_version = ""
        conversions = []
        geodetic_suffix = "_lon _lat _ehgt".split()
        xyz_suffix = "_X _Y _Z".split()
        coord_cols = []

        for cs in (args.coordinate_systems or "NZGD2000").upper().split("+"):
            match = re.match(r"^(NZGD2000|ITRF(?:19|20)\d\d)(_XYZ)?$", cs)
            if not match:
                print("Invalid coordinate system requested:", cs)
                sys.exit()
            csbase = match.group(1)
            isgeodetic = match.group(2) != "_XYZ"
            transformation = None
            if csbase == "NZGD2000":
                if not itrf_nzgd2000:
                    defmodel = loadDeformationModel(args)
                    nzgd2000_version = defmodel.version()
                    itrf_nzgd2000 = ITRF_NZGD2000.Transformation(
                        itrf=cors_itrf, model=defmodel
                    )

                def transformation(xyz, date):
                    llh = GRS80.geodetic(xyz)
                    llh2k = itrf_nzgd2000.transform(llh[0], llh[1], llh[2], date)
                    return GRS80.xyz(llh2k)

            else:
                transformation = ITRF.Transformation(
                    from_itrf=cors_itrf, to_itrf=csbase
                ).transform
            conversions.append((transformation, isgeodetic))
            coord_cols.extend(
                (
                    csbase + suffix
                    for suffix in (geodetic_suffix if isgeodetic else xyz_suffix)
                )
            )

        check_code = lambda code: True
        if args.stations:
            stations = []
            for s in args.stations.split("+"):
                if s.startswith("@"):
                    try:
                        with open(s[1:]) as sf:
                            stations.extend(sf.read().upper().split())
                    except:
                        raise RuntimeError("Cannot open station list file " + s[1:])
                elif s != "":
                    stations.append(s.upper())
            if len(stations) == 0:
                raise RuntimeError("No stations specifed in " + args.stations)
            check_code = lambda code: code.upper() in stations

        check_xyz = lambda xyz: True
        if args.nz_only:
            nzxyz = GRS80.xyz(177.0, -42.0)
            nzdist = 1000000.0
            check_xyz = lambda xyz: (
                math.sqrt(
                    (xyz[0] - nzxyz[0]) ** 2
                    + (xyz[1] - nzxyz[1]) ** 2
                    + (xyz[2] - nzxyz[2]) ** 2
                )
                < nzdist
            )

        for output_csv_file in sorted(outputfiles):
            calcdates = outputfiles[output_csv_file]
            ncoord = 0
            buildfile = output_csv_file + ".temp"
            with open(buildfile, "w") as csvf:
                writer = csv.writer(csvf)
                header = "code epoch".split()
                if nzgd2000_version:
                    header.append("nzgd2000_version")
                header.extend(coord_cols)
                writer.writerow(header)
                for code in tslist.codes():
                    if not check_code(code):
                        continue
                    for calcdate in calcdates:
                        ts = tslist.get(code).getData(enu=False)
                        try:
                            crd = ts.ix[calcdate]
                        except KeyError:
                            continue

                        if type(crd) is pd.DataFrame:
                            print(
                                "Ambiguous coordinates {0} at date {1}".format(
                                    code, calcdate
                                )
                            )
                            crd = crd[-1:].ix[calcdate]
                        xyz_2008 = [crd.x, crd.y, crd.z]
                        if not check_xyz(xyz_2008):
                            continue

                        row = [code, calcdate.strftime("%Y-%m-%d")]
                        if nzgd2000_version:
                            row.append(nzgd2000_version)
                        for transformation, isgeodetic in conversions:
                            try:
                                xyz = transformation(xyz_2008, calcdate)
                                if isgeodetic:
                                    llh = GRS80.geodetic(xyz)
                                    row.extend(
                                        [
                                            "{0:.9f}".format(llh[0]),
                                            "{0:.9f}".format(llh[1]),
                                            "{0:.4f}".format(llh[2]),
                                        ]
                                    )
                                else:
                                    row.extend(
                                        [
                                            "{0:.4f}".format(xyz[0]),
                                            "{0:.4f}".format(xyz[1]),
                                            "{0:.4f}".format(xyz[2]),
                                        ]
                                    )
                            except OutOfRangeError:
                                row.extend(["", "", ""])
                        writer.writerow(row)
                        ncoord += 1
            if ncoord == 0:
                os.unlink(buildfile)
                print(output_csv_file + " not built as coordinates not available")
            else:
                os.rename(buildfile, output_csv_file)

    except RuntimeError as e:
        print(e.message)


if __name__ == "__main__":
    main()
