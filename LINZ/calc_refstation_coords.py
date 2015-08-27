
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

def main():
    calcDate=datetime(date.today().year,1,1)

    parser=argparse.ArgumentParser('Script to compare NZGD2000 coordinates with station time series')
    parser.add_argument('gdb_file',default='gdb_coords.txt',help="Name of gdb file input file")
    parser.add_argument('check_file',nargs='?',default='coord_differences.csv',help="Output csv file of coordinate comparisons")
    parser.add_argument('update_csv_file',nargs='?',help="Name of gdb coordinate upload file (default none)")
    parser.add_argument('-f','--def-dir',help="Deformation model base dir (contains tools and model)")
    parser.add_argument('-m','--model-dir',help="Directory containing deformation model")
    parser.add_argument('-s','--stn-dir',default='stations',help="Directory containing SPM xml definitions")
    parser.add_argument('-d','--calc-date',help="Calculation date for coordinates (default "+calcDate.strftime("%Y-%m-%d")+")")

    args=parser.parse_args()

    gdbfile=args.gdb_file
    chkfile=args.check_file
    updfile=args.update_csv_file
    stndir=args.stn_dir
    defdir=args.def_dir

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

    gdbcrds={}
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
                gdbcrds[parts[0].upper()]=crds
            except:
                pass

    csvu=None
    if updfile:
        csvu=open(updfile,'w')

    with open(chkfile,'w') as csv:
        csv.write(','.join((
            'code',
            'itrf2008_X','itrf2008_Y','itrf2008_Z',
            'itrf2008_lon','itrf2008_lat','itrf2008_hgt',
            'itrf96_lon','itrf96_lat','itrf96_hgt',
            'nzgd2000_lon','nzgd2000_lat','nzgd2000_hgt',
            'gdb_lon','gdb_lat','gdb_hgt',
            'e_diff','n_diff','h_diff',
        )))
        csv.write('\n');
        if csvu:
            csvu.write("CODE,LAT,LON,EHGT,ROBG,COSY,DATE\n");


        codes=[]
        for f in os.listdir(stndir):
            m=re.match(r'^(\w{4}).xml$',f)
            if not m:
                continue
            codes.append(m.group(1))

        for code in sorted(codes):
            f=code+'.xml'
            if code not in gdbcrds:
                continue

            print("Processing",code)

            try:
                m=spm(filename=stndir+'/'+f)
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
                csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(*xyz08))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh08))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llh96))
                csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*llhnz2k))

                ucode=code.upper()
                if ucode in gdbcrds:
                    gcrds=gdbcrds[ucode] 
                    csv.write(',{0:.9f},{1:.9f},{2:.4f}'.format(*gcrds))
                    dedln,dndlt=GRS80.metres_per_degree(*gcrds)
                    edif=(gcrds[0]-llhnz2k[0])*dedln
                    ndif=(gcrds[1]-llhnz2k[1])*dndlt
                    hdif=gcrds[2]-llhnz2k[2]
                    csv.write(',{0:.4f},{1:.4f},{2:.4f}'.format(edif,ndif,hdif))
                else:
                    csv.write(',,,,,')
                csv.write("\n")
                if csvu and ucode in gdbcrds:
                    csvu.write('"{0}",{1:.9f},{2:.9f},{3:.4f},"B10","NZGD2000","2000.01.01"\n'.format(
                        code.upper(),*llhnz2k))
            except:
                print(sys.exc_info()[1])

if __name__=="__main__":
    main()
