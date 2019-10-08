#!/usr/bin/python
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import os.path
import re
import argparse
from datetime import datetime,date

from .StationCoordinateModel import Model as spm

# Note: This code needs to be refactored! Structure has become very ugly!

_description='Script to compare NZGD2000 coordinates with station time series'
_epilog='''

NZGD2000 coordinates of the reference stations can either be specified
in a CSV data file with columns 'geodeticcode lat lon ellhgt', or 
columns 'code gdb_lon gdb_lat gdb_hgt'.  Alternatively if the filename
is 'gdb' then the current coordinates will be downloaded from the 
online geodetic database.

'''
model_itrf='ITRF2008'
default_gx_error=['0.5','1.0']


def main():
    calcDate=datetime(date.today().year,1,1)

    parser=argparse.ArgumentParser(_description,epilog=_epilog)
    parser.add_argument('gdb_file',default='gdb_coords.txt',help="Name of gdb file input file, or 'gdb' to use online database, or 'none' to not compare with gdb")
    parser.add_argument('check_file',nargs='?',default='coord_differences.csv',help="Output csv file of coordinate comparisons")
    parser.add_argument('-u','--update-csv-file',help="Name of gdb coordinate upload file (default none)")
    parser.add_argument('-x','--snap-gx-file',help="Name of snap data file (default none)")
    parser.add_argument('-c','--snap-crd-file',help="Name of snap coordinate file (default none)")
    parser.add_argument('-f','--def-dir',help="Deformation model base dir (contains tools and model)")
    parser.add_argument('-m','--model-dir',help="Directory containing deformation model")
    parser.add_argument('-s','--stn-dir',default='stations',help="Directory containing SPM xml definitions")
    parser.add_argument('-d','--calc-date',help="Calculation date for coordinates (default "+calcDate.strftime("%Y-%m-%d")+")")
    parser.add_argument('-a','--all-marks',action='store_true',help="Output coord information for all marks including not in gdb")
    parser.add_argument('-e','--extend-dates',action='store_true',help="Extrapolate models before first or after last dates")
    parser.add_argument('-i','--itrf-only',action='store_true',help="Calculate ITRF coordinates but not GDB")
    parser.add_argument('--snap-gx-error',nargs=2,type=float,
                        default=default_gx_error,help="Snap horizontal and vertical mm error")
    parser.add_argument('--snap-gx-refframe',default=model_itrf,help="Snap GX observation reference frame code")
    parser.add_argument('--snap-crd-use-itrf',action='store_true',help="Snap GX observation reference frame code")

    args=parser.parse_args()

    gdbfile=args.gdb_file
    chkfile=args.check_file
    updfile=args.update_csv_file
    snapgxfile=args.snap_gx_file
    snapgxrf=args.snap_gx_refframe
    snapcrdfile=args.snap_crd_file
    stndir=args.stn_dir
    defdir=args.def_dir
    allmarks=args.all_marks
    alldates=args.extend_dates
    itrfonly=args.itrf_only
    gx_error=args.snap_gx_error
    crd_itrf=args.snap_crd_use_itrf

    mdldir=args.model_dir or (defdir+'/model' if defdir else 'model')

    if args.calc_date is not None:
        calcDate=datetime.strptime(args.calc_date,"%Y-%m-%d")

    if args.def_dir:
        sys.path.append(defdir+'/tools/LINZ')
    from LINZ.DeformationModel import Model as DefModel
    from LINZ.Geodetic.Ellipsoid import GRS80
    from LINZ.Geodetic.ITRF import Transformation

    itrf_tfm=Transformation(from_itrf=model_itrf,to_itrf='ITRF96')


    calcgdb=not itrfonly
    if gdbfile == 'gdb':
        from LINZ.Geodetic import GDB
        GDB.setCached()
        def gdbcrds( code ):
            try:
                markdata=GDB.get(code)
                coord=markdata.coordinate
                return [coord.longitude,coord.latitude,coord.height]
            except:
                return None
    elif gdbfile == 'none': 
        gdbcrds=lambda code: None
        allmarks=True
        calcgdb=False
    else:
        markcrds={}
        with open(gdbfile,'r') as gdbf:
            l=gdbf.readline()
            l=l.lower().replace("\"","").replace(","," ")
            if l.split() == 'code gdb_lon gdb_lat gdb_hgt'.split():
                crdorder=(1,2,3)
            elif l.split() == 'geodeticcode lat lon ellhgt'.split():
                crdorder=(2,1,3)
            else:
                raise RuntimeError("Invalid fields in "+gdbfile)
            for l in gdbf:
                l=l.lower().replace("\"","").replace(","," ")
                try:
                    parts=l.split()
                    crds=[float(parts[x]) for x in crdorder]
                    markcrds[parts[0].upper()]=crds
                except:
                    pass
            gdbcrd=lambda code: markcrds.get(code)

    needgdb=False
    neednz2k=False
    csv=None
    if chkfile and chkfile != 'none':
        csv=open(chkfile,'w')
        itrfc=model_itrf.lower()
        columns=[
            'code',
            'scm_version',
            'deformation_version',
            'calc_date',
            itrfc+'_X',itrfc+'_Y',itrfc+'_Z',
            itrfc+'_lon',itrfc+'_lat',itrfc+'_hgt',
            ]
        if not itrfonly:
            neednz2k=True
            columns.extend([
                'itrf96_lon','itrf96_lat','itrf96_hgt',
                'nzgd2000_lon','nzgd2000_lat','nzgd2000_hgt',
                ])
        if calcgdb:
            needgdb=True
            columns.extend([
                'gdb_lon','gdb_lat','gdb_hgt',
                'e_diff','n_diff','h_diff',
                ])
        csv.write(','.join(columns))
        csv.write("\n")

    csvu=None
    logf=None
    if updfile:
        needgdb=True
        neednz2k=True
        csvu=open(updfile,'w')
        csvu.write("CODE,LAT,LON,EHGT,ROBG,COSY,DATE,COMM\n");
        updcomment="SPM version {0:%Y-%m-%d %H:%M:%S}, Deformation model version {1}, Calc date {2:%Y-%m-%d %H:%M:%S}"
        logfile=updfile+'.log'
        logf=open(logfile,'w')

    datf=None
    if snapgxfile:
        datf=open(snapgxfile,'w')
        datf.write("CORS reference stations coordinates\n")
        datf.write("\n")
        datf.write("#date {0:%Y-%m-%d}\n".format(calcDate))
        datf.write("#gx_enu_error {0} {0} {1} mm\n".format(*gx_error))
        datf.write("#reference_frame {0}\n".format(snapgxrf))
        datf.write("#classify gx source scm\n".format(model_itrf))
        datf.write("#classification scmversion\n".format(model_itrf))
        datf.write("#classification scmlastobs\n".format(model_itrf))
        datf.write("#data gx scmversion scmlastobs value no_heights\n\n".format(model_itrf))

    crdf=None
    if snapcrdfile:
        neednz2k=not crd_itrf
        crdf=open(snapcrdfile,'w')

    if neednz2k:
        defmodel=DefModel.Model(mdldir)

    if logf:
        logf.write("calc_refstation_coordinates\n\n")
        logf.write("Run time: {0:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()))
        logf.write("Deformation model: {0}\n".format(mdldir))
        logf.write("Deformation model version: {0}\n".format(defmodel.version()))
        logf.write("Station coordinate model directory: {0}\n".format(stndir))
        logf.write("Coordinate comparisons written to: {0}\n".format(chkfile))
        logf.write("GDB update file written to: {0}\n".format(updfile))
        logf.write("Coordinate calculation date: {0:%Y-%m-%d %H:%M:%S}\n".format(calcDate))
        logf.write("GDB coordinates: {0}\n".format(gdbfile))

    if crdf:
        crdf.write("PositioNZ coordinates calculated at {0:%Y-%m-%d %H:%M:%S}\n".format(calcDate))
        if crd_itrf:
            crdf.write("{0}@{1:%Y%m%d}\n".format(model_itrf,calcDate))
        else:
            crdf.write("NZGD2000_{0}\n".format(defmodel.version()))
        crdf.write("options no_geoid ellipsoidal_heights degrees c=scmversion c=scmlastobs\n\n")

    codes=[]
    for f in os.listdir(stndir):
        m=re.match(r'^(\w{4}).xml$',f)
        if not m:
            continue
        codes.append(m.group(1))

    for code in sorted(codes):
        f=code+'.xml'
        
        gcrds=gdbcrds(code)
        if gcrds is None and not allmarks:
            continue

        print("Processing",code)

        try:
            m=spm(filename=stndir+'/'+f)
            if logf is not None:
                logf.write("{0} model version: {1:%Y-%m-%d %H:%M:%S}\n".format(code,m.versiondate))
            # Don't want annual and semi-annual components...
            for c in m.components:
                if 'annual' in c.componentType():
                    c.setEnabled(False)

            if not alldates: 
                if m.startdate is not None and calcDate < m.startdate:
                   continue
                if m.enddate is not None and calcDate  > m.enddate:
                   continue

            xyz08=m.calc(calcDate,enu=False)
            llh08=GRS80.geodetic(xyz08)
            if neednz2k:
                llh96=itrf_tfm.transformLonLat(llh08[0],llh08[1],llh08[2],calcDate)
                if llh96[0] < 0:
                    llh96[0] += 360.0
                llhnz2k=defmodel.applyTo(llh96,date=calcDate,subtract=True)
                if llh96[0] > 180:
                    llh96[0] -= 360.0
                if llhnz2k[0] > 180:
                    llhnz2k[0] -= 360.0

            if csv:
                csv.write('"{0}"'.format(code))
                csv.write(',"{0:%Y-%m-%d %H:%M:%S}"'.format(m.versiondate))
                csv.write(',"{0}"'.format(defmodel.version()))
                csv.write(',"{0:%Y-%m-%d %H:%M:%S}"'.format(calcDate))
                csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(*xyz08))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh08))
                if not itrfonly:
                    csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh96))
                    csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llhnz2k))

                    ucode=code.upper()
                if calcgdb:
                    if gcrds is not None:
                        csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*gcrds))
                        dedln,dndlt=GRS80.metres_per_degree(*gcrds)
                        edif=(gcrds[0]-llhnz2k[0])*dedln
                        ndif=(gcrds[1]-llhnz2k[1])*dndlt
                        hdif=gcrds[2]-llhnz2k[2]
                        csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(edif,ndif,hdif))
                    else:
                        csv.write(',,,,,')
                csv.write("\n")
            if csvu:
                comment=updcomment.format(m.versiondate,defmodel.version(),calcDate).replace('"','""')
                csvu.write('"{0}",{2:.9f},{1:.9f},{3:.4f},"B10","NZGD2000","2000.01.01","{4}"\n'.format(
                    code.upper(),llhnz2k[0],llhnz2k[1],llhnz2k[2],comment))

            scmlastobs="{0:%Y%m%d}".format(m.enddate) if m.enddate is not None else "-"
            scmversion="SCM_{0:%Y%m%d}".format(m.versiondate)
            if datf:
                datf.write("{0} {1} {2} {3:.4f} {4:.4f} {5:.4f}\n".format(
                    code.upper(),scmversion,scmlastobs,xyz08[0],xyz08[1],xyz08[2]))

            if crdf:
                snapcrd=llh08 if crd_itrf else llhnz2k
                crdf.write("{0} {1:.10f} {2:.10f} {3:.4f} {4} {5} {0}\n".format(code.upper(),
                    snapcrd[1],snapcrd[0],snapcrd[2],scmversion,scmlastobs))
        except:
            print(sys.exc_info()[1])

    for f in (csv, csvu, logf, datf, crdf ):
        if f is not None:
            f.close()

if __name__=="__main__":
    main()
