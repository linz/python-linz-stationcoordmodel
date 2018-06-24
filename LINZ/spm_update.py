#!/usr/bin/python
"""
Recalculation of station coordinate model from time series
"""

# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import os.path
import argparse
import datetime as dt
from datetime import date, datetime, timedelta
import re

import numpy as np
import pandas as pd

from . import StationCoordinateModel as spm
from .CORS_Timeseries import TimeseriesList, robustStandardError
from LINZ.Geodetic.Ellipsoid import GRS80

default_model_file='stations/{code}.xml'
default_model_backup_file=None # stations/{code}.xml.{fdatetime}
default_timeseries_file='timeseries/{code}_igs08_xyz.dat'


class StationCoordModelUpdater( object ):

    summary_cols='''station longitude latitude
                   start_date end_date dh dv
                   end_de end_dn end_du 
                   end_mean_de end_mean_dn end_mean_du
                   residual_date res_de res_dn res_du'''.split()

    def __init__(self, cfgfile, options=None, testonly=None, verbose=False ):
        self.testonly=testonly
        self.verbose=verbose
        self.read_config(cfgfile,options)
        self.backedup=set()

    def read_config(self,cfgfile,options=None):
        config={x:None for x in '''
                model_file model_backup_file timeseries_file timeseries_type 
                before after
                update_availability
                robust_se_percentile 
                outlier_reject_level
                outlier_test_range
                set_model_start_date'''.split()}
        if cfgfile and os.path.exists(cfgfile):
            with open(cfgfile) as cfg:
                for l in cfg:
                    if re.match(r'^\s*(\#|$)',l):
                        continue
                    parts=re.split(r'\s+',l.strip(),1)
                    if len(parts) != 2:
                        parts.append(None)
                    if parts[0] not in config:
                        continue
                        # raise RuntimeError('Invalid configuration item: '+parts[0])
                    config[parts[0]]=parts[1]

        if options:
            config.update(options)
        config={k:v for k,v in config.iteritems() if v is not None}
        self.config=config

        self.model_file=config.get('model_file',default_model_file)
        self.model_backup_file=config.get('model_backup_file',default_model_backup_file)
        self.update_availability=config.get('update_availability','False').lower()=='true'
        self.set_start_date=config.get('set_model_start_date','False').lower()=='true'
        if '{code}' not in self.model_file:
            raise RuntimeError('Configuration item model_file must include "{code}"')
        if ( self.model_backup_file and 
             '{code}' not in self.model_backup_file and 
             '{model_file}' not in self.model_backup_file ):
            raise RuntimeError('Configuration item model_backup_file must include "{model_file}" or "{code}"')
        self.timeseries_file=config.get('timeseries_file',default_timeseries_file)
        self.robust_se_percentile=float(config.get('robust_se_percentile','95.0'))
        self.outlier_reject_level=float(config.get('outlier_reject_level','5.0'))
        self.outlier_test_range=float(config.get('outlier_test_range','10.0'))
        self.auto_reject_obs=str(config.get('auto_reject_obs','False')).lower()=='true'
        solutiontypes=config.get('timeseries_type','')
        self.solutiontypes=solutiontypes
        self.timeseries_list=TimeseriesList(self.timeseries_file,solutiontypes or None,
              after=config.get('after'),before=config.get('before'))

    def backupModel( self ):
        if self.model is None or self.model_backup_file is None:
            return
        code = self.model.station
        if code in self.backedup:
            return
        self.backedup.add(code)
        oldfile=self.model.filename
        if not os.path.exists(oldfile):
            return
        filedate=dt.datetime.fromtimestamp(os.path.getmtime(oldfile))
        ymd=filedate.strftime('%Y%m%d')
        ymdhms=filedate.strftime('%Y%m%d%H%M%S')
        backupfile = (self.model_backup_file.replace('{code}',code).
                      replace('{model_file}',oldfile).
                      replace('{fdate}',ymd).
                      replace('{fdatetime}',ymdhms))
        try:
            os.renames(oldfile,backupfile)
        except:
            shutil.copy(oldfile,backupfile)

    def saveModel( self ):
        if not self.model:
            return
        if self.model.changed():
            self.backupModel();
        self.model.save(updateAvailability=self.update_availability)

    def loadModel(self,code):
        filename=self.model_file.replace('{code}',code)
        try:
            self.model=spm.Model(station=code,filename=filename,loadfile=True)
            self.original_model=self.model.copy()
            self.timeseries=self.timeseries_list.get(code)
            self.model.loadTimeSeries(self.timeseries,setdate=self.set_start_date)
        except RuntimeError as ex:
            self.model=None
            raise

    def calculateStats( self ):
        self.obs_rse=[0,0,0]
        self.residual_rse=[0,0,0]
        self.obs_count=0
        dates,obs,useobs = self.model.getObs()
        if dates is None or obs is None or len(obs) == 0:
            return

        self.obs_count=len(dates)
        self.obs_rse=robustStandardError(obs)
        residuals=[0,0,0]
        calc=None
        if self.model:
            calc=self.model.calc(dates)
            diff=obs-calc
            # This code was to remove mean difference between observed and calced..
            # Should not be doing this!
            # offset=np.mean(diff,axis=0)
            # calc += offset
            # diff -= offset

            # Calculate a "robust residual" measure based on 95%ile
            residuals=np.percentile(np.abs(diff[useobs]),95.0,axis=0)/(1.96*np.sqrt(2.0))
        return residuals

    def modelFile(self,code):
        filename=self.model_file.replace('{code}',code)
        return filename

    def stationCodes(self):
        for code in self.timeseries_list.codes():
            if os.path.exists(self.modelFile(code)):
                yield code

    def autoRejectObs( self ):
        if self.model:
            self.model.autoRejectObs(
                ndays=self.outlier_test_range,
                tolerance=self.outlier_reject_level,
                percentile=self.robust_se_percentile)

    def fitAllLinear(self):
        ok, message = self.model.fitAllLinear()
        if not ok:
            raise RuntimeError('Error fitting model: '+message)

    def _calcChange( self, dates, excludeTypes=None ):
        return (self.model.calc(dates,excludeTypes=excludeTypes) 
                - self.original_model.calc(dates,excludeTypes=excludeTypes))

    def calcChange( self, ndays=10, start_date=None, end_date=None ):
        '''
        Calculates coordinate changes between original and current model.
        '''
        if start_date is None:
            start_date=self.model.refdate
        if start_date is None:
            start_date=spm.refdate
        if end_date is None:
            end_date=dt.datetime.combine(dt.date.today(),dt.time(0,0,0))
        calc_date=start_date
        dates=[calc_date]
        tdelta=dt.timedelta(days=ndays)
        while calc_date < end_date:
            calc_date += tdelta
            dates.append(calc_date)
        change=self._calcChange(dates)
        meanChange=self._calcChange(dates,excludeTypes=spm.cyclic)
        return dates, change, meanChange

    def changeSummary( self ):
        dates, change, meanChange=self.calcChange()
        lon,lat,h=GRS80.geodetic(self.model.calc(dates[-1],enu=False))
        if lon < 0:
            lon += 360.0
        return dict(
            station=self.model.station,
            longitude=lon,
            latitude=lat,
            start_date=dates[0],
            end_date=dates[-1],
            dh=np.max(np.hypot(change[:,0],change[:,1])),
            dv=np.max(np.abs(change[:,2])),
            end_de=change[-1,0],
            end_dn=change[-1,1],
            end_du=change[-1,2],
            end_mean_de=meanChange[-1,0],
            end_mean_dn=meanChange[-1,1],
            end_mean_du=meanChange[-1,2]
            )

    def residuals( self ):
        residual=self.timeseries.subtract(lambda d: self.model.calc(d,enu=False))
        residual=residual.withoutOutliers()
        return residual

    def residualSummary( self, ndays=30 ):
        '''
        Calculate the mean enu difference of the last n days of residuals
        '''
        residuals=self.residuals()
        dates,diff=residuals.getObs()
        ndays=int(min(len(dates),ndays))
        nmid=int((ndays+1)/2)
        residual_date=dates[-nmid]
        denu=np.mean(diff[-ndays:],axis=0)
        return dict(
            station=self.model.station,
            residual_date=residual_date,
            res_de=denu[0],
            res_dn=denu[1],
            res_du=denu[2]
            )

    def updateStation( self, code ):
        if self.verbose:
            print("Updating {0}".format(code))
        try:
            self.loadModel( code )
            if self.auto_reject_obs:
                self.autoRejectObs()
            self.fitAllLinear()
            if not self.testonly:
                self.saveModel()
        except RuntimeError as ex:
            print("Station {0}: {1}".format(code,ex.message))
            return False
        return True

    def update( self, codes=None, calcSummary=None ):
        if not codes:
            codes=list(self.stationCodes())
        summary=[]
        for code in codes:
            ok=self.updateStation(code)
            if ok and calcSummary:
                change=self.changeSummary()
                residuals=self.residualSummary()
                summary.append([change.get(c,residuals.get(c)) 
                                for c in StationCoordModelUpdater.summary_cols])
        if self.verbose and self.testonly:
            print("Station coordinate model updates NOT saved")
        result=None
        if calcSummary:
            result=pd.DataFrame(summary,columns=StationCoordModelUpdater.summary_cols)
        return result

def main():
    cfgfile='spm_editor.cfg'
    usercfg='~/.'+os.path.basename(cfgfile)
    if os.path.exists(usercfg):
        cfgfile=usercfg
    usercfg=os.path.basename(cfgfile)
    if os.path.exists(usercfg):
        cfgfile=usercfg
    parser=argparse.ArgumentParser(description='Refresh station coordinate models by linear fitting to data')
    parser.add_argument('config_file',nargs='?',help='Station prediction model configuration file')
    parser.add_argument('-c','--codes',nargs='+',help='Station codes to evaluate')
    parser.add_argument('-r','--auto-reject',action='store_true',help='Auto reject observations')
    parser.add_argument('-v','--verbose',action='store_true',help='Generate more output')
    parser.add_argument('-s','--summary-file',help='Name of change/residual summary CSV file')

    args=parser.parse_args()
    options=dict(
        auto_reject_obs=args.auto_reject
        )

    codes=[]
    if args.codes is not None:
        codes=[c.upper() for c in args.codes]
        for c in codes:
            if not re.match(r'^\w\w\w\w$',c):
                print("Invalid station code {0}".format(c))
                sys.exit()

    if args.config_file:
        cfgfile=args.config_file

    updater=StationCoordModelUpdater(cfgfile,options,verbose=args.verbose)
    summary=updater.update(codes,calcSummary=args.summary_file is not None)
    if summary is not None and args.summary_file:
        summary.to_csv(args.summary_file,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    
if __name__ == "__main__":
    main()
