import sys
import os.path
import numpy as np
import datetime as dt
import pandas as pd
from LINZ.geodetic.ellipsoid import grs80


class CORS_timeseries( object ):

    _columns=('name','epoch','x','y','z')
    _epoch_col='epoch'

    @staticmethod
    def robustStandardError(obs):
        '''
        Estimate the standard error for components of a time series.

        Standard errors are estimated using the differences between consecutive elements
        of the time series.  The 95 percentile of the absoluted differences is used as a
        robust estimator for the 95% cumulative probability (two tailed), so is scaled by 
        1.96 to get standard error.
        '''
        errors=[0]*3
        for axis in range(3):
            diffs=np.abs(obs[1:,axis]-obs[:-1,axis])
            # Note sqrt(2) accounts for the fact that these are differences
            # between two observations
            se=np.percentile(diffs,95.0)/(1.96*np.sqrt(2))
            errors[axis]=se
        return errors

    def __init__( self, filename, name=None, xyz0=None, transform=None ):
        '''
        Loads the CORS time series from a file for analsis.

        Assumes the time series file has columns name, epoch, x, y, z and is 
        whitespace delimited.  

        The xyz0 parameter can define a reference point used for calculating e, n, u 
        components.  The default is to use the first coordinate.

        The transform argument can define a function which is applied to each coordinate 
        as it is read.  The transform function is called for each coordinate as 

        xyz=transform(xyz,epoch)

        where epoch is a datetime.datetime object

        Flags observations as excluded if they match the dates stored with the model.
        '''

        if not os.path.exists(filename):
            raise RuntimeError("CORS time series file {0} does not exist".format(filename))

        columns=[]
        try:
            with open(filename) as f:
                columns=f.readline().split()
        except:
            raise RuntimeError('CORS time series file {0} missing or empty'.format(filename))

        missing_cols=[x for x in self._columns if x not in columns]
        if missing_cols:
            raise RuntimeError('Columns {0} missing from CORS time series file {1}'.format
                                   (' '.join(missing_cols),filename))
        try:
            data=pd.read_csv(filename,delim_whitespace=True,parse_dates=True,index_col=self._epoch_col)
        except:
            raise RuntimeError("Cannot read time series data from {0}: {1}".format
                               (filename,sys.exc_info()[1]))

        name=name if name else data.name.iloc[0]
        data=data[data.name==name]
        data=data.sort_index()

        xyz=np.vstack((data.x,data.y,data.z)).T
        if transform:
            xyz=[transform(x) for x in xyz]
            data['x']=xyz[:,0]
            data['y']=xyz[:,1]
            data['z']=xyz[:,2]

        xyz0=np.array(xyz0) if xyz0 is not None else xyz[0]
        lon,lat=grs80.geodetic(xyz0)[:2]
        enu_axes=grs80.enu_axes(lon,lat)

        diff=xyz-xyz0
        enu=(xyz-xyz0).dot(enu_axes.T)
        data['e']=enu[:,0]
        data['n']=enu[:,1]
        data['u']=enu[:,2]

        self._filename=filename
        self._name=name
        self._xyz0=xyz0
        self._data=data

        # Add instance version of class methods

        self.robustStandardError=lambda: CORS_timeseries.robustStandardError(self.getObs()[1])

    def filename( self ):
        return self._filename

    def xyz0( self ):
        return self._xyz0

    def name( self ):
        return self._name

    def getObs( self, enu=True, index=None ):
        '''
        Returns the time series as time and date arrays

        Parameters:
            enu      Select the ENU data rather than XYZ (default True)
            index    Pandas data frame index to extract subset of data

        Returns two values:
            dates    An array of dates for the time series
            enu      An numpy array of east,north,up values relative to the
                     model reference coordinate, or xyz coordinates if enu is False
        '''
        data=self._data
        if index:
            data=data[index]
        cols=('e','n','u') if enu else ('x','y','z')
        return data.index.to_pydatetime(),np.vstack((data[x] for x in cols)).T

    def getData( self, enu=True, index=None ):
        ''' 
        Returns the time series as a pandas DataFrame. 

        Parameters:
            enu      Select the ENU data rather than XYZ (default True)
            index    Pandas data frame index to extract subset of data
        '''
        data=self._data
        if index:
            data=data[index]
        return data[['e','n','u']] if enu else data[['x','y','z']]

