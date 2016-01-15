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

_description='Script to compare NZGD2000 coordinates with station time series'
_epilog='''

NZGD2000 coordinates of the reference stations can either be specified
in a CSV data file with columns 'geodeticcode lat lon ellhgt', or 
columns 'code gdb_lon gdb_lat gdb_hgt'.  Alternatively if the filename
is 'gdb' then the current coordinates will be downloaded from the 
online geodetic database.

'''


def main():
    calcDate=datetime(date.today().year,1,1)

    parser=argparse.ArgumentParser(_description,epilog=_epilog)
    parser.add_argument('gdb_file',default='gdb_coords.txt',help="Name of gdb file input file, or 'gdb' to use online database")
    parser.add_argument('check_file',nargs='?',default='coord_differences.csv',help="Output csv file of coordinate comparisons")
    parser.add_argument('update_csv_file',nargs='?',help="Name of gdb coordinate upload file (default none)")
    parser.add_argument('-f','--def-dir',help="Deformation model base dir (contains tools and model)")
    parser.add_argument('-m','--model-dir',help="Directory containing deformation model")
    parser.add_argument('-s','--stn-dir',default='stations',help="Directory containing SPM xml definitions")
    parser.add_argument('-d','--calc-date',help="Calculation date for coordinates (default "+calcDate.strftime("%Y-%m-%d")+")")
    parser.add_argument('-a','--all-marks',help="Output coord information for all marks including not in gdb")

    args=parser.parse_args()

    gdbfile=args.gdb_file
    chkfile=args.check_file
    updfile=args.update_csv_file
    stndir=args.stn_dir
    defdir=args.def_dir
    allmarks=args.all_marks

    mdldir=args.model_dir or (defdir+'/model' if defdir else 'model')

    if args.calc_date != '':
        calcDate=datetime.strptime(args.calc_date,"%Y-%m-%d")

    if args.def_dir:
        sys.path.append(defdir+'/tools/LINZ')
    from LINZ.DeformationModel import Model as DefModel
    from LINZ.Geodetic.Ellipsoid import GRS80
    from LINZ.Geodetic.ITRF import Transformation

    itrf_tfm=Transformation(from_itrf='ITRF2008',to_itrf='ITRF96')
    defmodel=DefModel.Model(mdldir)

    logfile=None
    logf=None
    if updfile:
        logfile=updfile+'.log'
        logf=open(logfile,'w')
        logf.write("calc_refstation_coordinates\n\n")
        logf.write("Run time: {0:%Y-%m-%d %H:%M:%S}\n".format(datetime.now()))
        logf.write("Deformation model: {0}\n".format(mdldir))
        logf.write("Deformation model version: {0}\n".format(defmodel.version()))
        logf.write("Station coordinate model directory: {0}\n".format(stndir))
        logf.write("Coordinate comparisons written to: {0}\n".format(chkfile))
        logf.write("GDB update file written to: {0}\n".format(updfile))
        logf.write("Coordinate calculation date: {0:%Y-%m-%d %H:%M:%S}\n".format(calcDate))

    if gdbfile == 'gdb':
        from LINZ.Geodetic import GDB
        def gdbcrds( code ):
            try:
                markdata=GDB.get(code)
                return markdata.official_coordinate
            except:
                return None
        if logf:
            logf.write("GDB coordinates: from online geodetic database\n")
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
        if logf:
            logf.write("GDB coordinates: {0}\n".format(gdbfile))

    csvu=None
    if updfile:
        csvu=open(updfile,'w')

    with open(chkfile,'w') as csv:
        csv.write(','.join((
            'code',
            'scm_version',
            'deformation_version',
            'calc_date',
            'itrf2008_X','itrf2008_Y','itrf2008_Z',
            'itrf2008_lon','itrf2008_lat','itrf2008_hgt',
            'itrf96_lon','itrf96_lat','itrf96_hgt',
            'nzgd2000_lon','nzgd2000_lat','nzgd2000_hgt',
            'gdb_lon','gdb_lat','gdb_hgt',
            'e_diff','n_diff','h_diff',
        )))
        csv.write('\n');
        if csvu:
            csvu.write("CODE,LAT,LON,EHGT,ROBG,COSY,DATE,COMM\n");

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

                xyz08=m.calc(calcDate,enu=False)
                llh08=GRS80.geodetic(xyz08)
                llh96=itrf_tfm.transformLonLat(llh08[0],llh08[1],llh08[2],calcDate)
                if llh96[0] < 0:
                    llh96[0] += 360.0
                llhnz2k=defmodel.applyTo(llh96,date=calcDate,subtract=True)
                if llh96[0] > 180:
                    llh96[0] -= 360.0
                if llhnz2k[0] > 180:
                    llhnz2k[0] -= 360.0

                csv.write('"{0}"'.format(code))
                csv.write(',"{0:%Y-%m-%d %H:%M:%S}"'.format(m.versiondate))
                csv.write(',"{0}"'.format(defmodel.version()))
                csv.write(',"{0:%Y-%m-%d %H:%M:%S}"'.format(calcDate))
                csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(*xyz08))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh08))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh96))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llhnz2k))

                ucode=code.upper()
                if gcrds is not None:
                    csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*gcrds))
                    dedln,dndlt=GRS80.metres_per_degree(*gcrds)
                    edif=(gcrds[0]-llhnz2k[0])*dedln
                    ndif=(gcrds[1]-llhnz2k[1])*dndlt
                    hdif=gcrds[2]-llhnz2k[2]
                    csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(edif,ndif,hdif))
                    if csvu:
                        comment="SPM version {0:%Y-%m-%d %H:%M:%S}, Deformation model version {1}, Calc date {2:%Y-%m-%d %H:%M:%S}"
                        comment=comment.format(m.versiondate,defmodel.version(),calcDate).replace('"','""')
                        csvu.write('"{0}",{2:.9f},{1:.9f},{3:.4f},"B10","NZGD2000","2000.01.01","{4}"\n'.format(
                            code.upper(),llhnz2k[0],llhnz2k[1],llhnz2k[2],comment))
                else:
                    csv.write(',,,,,')
                csv.write("\n")
            except:
                print(sys.exc_info()[1])

        if logf is not None:
            logf.close()

if __name__=="__main__":
    main()
