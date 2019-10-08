
# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import csv
import os.path
import numpy as np
from LINZ.Geodetic import GDB
from LINZ.Geodetic.Ellipsoid import GRS80

class CoordinateUpdate( object ):

    def __init__( self, code, lon=None, lat=None, hgt=None, robg='B10', cosy='NZGD2000', rdate='2000.01.01', comment='' ):
        self.code=code
        markdata=GDB.get(code)
        if not markdata:
            raise RuntimeError(code+' is not in the geodetic database')
        gdbcoord=[coord.longitude,coord.latitude,coord.height]
        if lon is None or lat is None or hgt is None:
            lon,lat,hgt=gdbcoord

        self.lon=lon
        self.lat=lat
        self.hgt=hgt
        self.robg=robg
        self.cosy=cosy
        self.rdate=rdate
        self.comment=comment
        self.gdbCoord=gdbcoord
        self.xyz0=GRS80.xyz(self.gdbCoord)
        self.enu_axes=GRS80.enu_axes(self.gdbCoord[0],self.gdbCoord[1])
        pass

    def xyzOffset( self ):
        xyz1=GRS80.xyz(self.lon,self.lat,self.hgt)
        dxyz=xyz1-self.xyz0
        return self.enu_axes.dot(dxyz)

    def setXyzOffset( self, denu, comment=None ):
        denu=np.array(denu)
        xyz=self.xyz0+self.enu_axes.T.dot(denu)
        self.lon,self.lat,self.hgt=GRS80.geodetic(xyz)
        if comment is not None:
            self.comment=comment
    

class Updates( object ):
    '''
    CSV file representing updates to station coordinates ready for loading into
    crd_crs script. (Mainly intended for updating CORS station coordinates)
    '''

    columns='CODE LON LAT EGHT ROBG COSY DATE COMM'.split()

    def __init__( self, csvFile ):
        self._csvFile=csvFile
        self.updates={}
        self.read()

    def read( self ):
        if not os.path.isfile(self._csvFile):
            return
        self.updates={}
        with open(self._csvFile,'rb') as csvfh:
            reader=csv.DictReader(csvfh)
            fieldnames=reader.fieldnames
            if fieldnames is None:
                return
            fields={f.upper():f for f in fieldnames}
            for c in Updates.columns[:-1]:
                if c not in fields:
                    raise RuntimeError(self._csvFile+" doesn't have column"+c)

            nerrors=0
            for row in reader:
                code=row[fields['CODE']]
                if code == '':
                    continue
                lon=row[fields['LON']] 
                lat=row[fields['LAT']]
                hgt=row[fields['EGHT']]
                robg=row[fields['ROBG']]
                cosy=row[fields['COSY']]
                rdate=row[fields['DATE']]
                comment=row.get(fields['COMM'],'')
                lon=float(lon) if lon != '' else None
                lat=float(lat) if lat != '' else None
                hgt=float(hgt) if hgt != '' else None
                try:
                    self.updates[code.upper()]=CoordinateUpdate(code,lon,lat,hgt,robg,cosy,rdate,comment)
                except:
                    nerrors += 1

    def save( self, filename=None ):
        filename = filename if filename is not None else self._csvFile
        with open(filename,'wb') as csvfh:
            writer=csv.writer(csvfh)
            writer.writerow(Updates.columns)
            for update in [self.updates[k] for k in sorted(self.updates.keys())]:
                row=(
                    update.code,
                    "{0:.10f}".format(update.lon),
                    "{0:.10f}".format(update.lat),
                    "{0:.4f}".format(update.hgt),
                    format(update.robg),
                    format(update.cosy),
                    format(update.rdate),
                    format(update.comment),
                    )
                writer.writerow(row)

    def getUpdate( self, code, canCreate=False ):
        code=code.upper()
        if code not in self.updates and canCreate:
            self.updates[code]=CoordinateUpdate(code)
        return self.updates.get(code,None)
