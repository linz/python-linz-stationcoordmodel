import sys
import os
import os.path
import re
import sqlite3

import numpy as np
import datetime as dt
import pandas as pd

from LINZ.Geodetic.Ellipsoid import GRS80

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

class Timeseries( object ):

        
    def __init__( self, code, solutiontype='default', 
                 data=None, 
                 dates=None, xyz=None,
                 xyz0=None, 
                 xyzenu=None,
                 transform=None,
                 after=None):
        '''
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
        after           The earliest date for data of interest
        '''

        self._loaded=False
        self._code=code
        self._solutiontype=solutiontype
        self._xyz0=xyz0
        self._xyzenu=xyzenu
        self._transform=transform
        self._after=after
        if xyz is not None and dates is not None:
            data=pd.DataFrame(xyz,columns=('x','y','z'))
            data.set_index(pd.to_datetime(dates),inplace=True)
        self._sourcedata=data

    def _loadData( self ):
        '''
        Abstract class, overridden in base classes
        '''
        return self._sourcedata.copy()

    def _load( self ):

        if self._loaded:
            return

        # Call provider specific data load
        # Should return a dataframe with x,y,z indexed on the date/time of the each point

        data=self._loadData()
        if data is None:
            raise RuntimeError('No data provided for time series')
        
        if self._after is not None:
            data=data[data.index > self._after]

        data.sort_index(inplace=True)
        xyz=np.vstack((data.x,data.y,data.z)).T
        if self._transform:
            xyz=[self._transform(x) for x in xyz]
            data['x']=xyz[:,0]
            data['y']=xyz[:,1]
            data['z']=xyz[:,2]

        xyz0=self._xyz0
        xyz0=np.array(xyz0) if xyz0 is not None else xyz[0]
        xyzenu=self._xyzenu
        xyzenu=np.array(xyzenu) if xyzenu is not None else xyz0
        lon,lat=GRS80.geodetic(xyzenu)[:2]
        enu_axes=GRS80.enu_axes(lon,lat)

        diff=xyz-xyz0
        enu=(xyz-xyz0).dot(enu_axes.T)
        data['e']=enu[:,0]
        data['n']=enu[:,1]
        data['u']=enu[:,2]

        self._xyz0=xyz0
        self._data=data

        # Add instance version of class methods

        self.robustStandardError=lambda: Timeseries.robustStandardError(self.getObs()[1])

    def setXyzTransform( self, xyz0=None, xyzenu=None, transform=None ):
        if xyz0 is not None and (
            self._xyz0 is None or
            (self._xyz0 != xyz0).any()):
            self._xyz0=xyz0
            self._loaded=False
        if xyzenu is not None and (
            self._xyzenu is None or
            (self._xyzenu != xyzenu).any()):
            self._xyzenu=xyzenu
            self._loaded=False
        if transform != self._transform:
            self._transform=transform
            self._loaded=False

    def xyz0( self ):
        return self._xyz0

    def xyzenu( self ):
        return self._xyzenu

    def code( self ):
        return self._code

    def solutiontype( self ):
        return self._solutiontype

    def setName( self, code ):
        self._code = code

    def trend( self, columns='e n u'.split() ):
        import matplotlib.dates as mdates
        self._load()
        columns=[columns] if isinstance(columns,basestring) else list(columns)
        series=self._data[columns]
        trends=[]
        days=mdates.date2num(self._data.index)
        pfit=np.polyfit(days,series,1)
        trends=[np.poly1d(pfit[:,i])(days) for i,c in enumerate(columns)]
        return pd.DataFrame(np.vstack(trends).T,index=self._data.index,columns=columns)

    def getObs( self, enu=True, detrend=True, index=None ):
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
        self._load()
        data=self._data
        if index:
            data=data[index]
        cols=('e','n','u') if enu else ('x','y','z')
        return data.index.to_pydatetime(),np.vstack((data[x] for x in cols)).T

    def getData( self, enu=True, index=None, normalize=False ):
        ''' 
        Returns the time series as a pandas DataFrame. 

        Parameters:
            enu      Select the ENU data rather than XYZ (default True)
            index    Pandas data frame index to extract subset of data
        '''
        self._load()
        data=self._data
        if index:
            data=data[index]
        result=data[['e','n','u']] if enu else data[['x','y','z']]
        if normalize:
            result.set_index(result.index.normalize(),inplace=True)
        return result

    def dates( self ):
        '''
        Returns the dates of the time series
        '''
        self._load()
        return self._data.index

    def plot( self, detrend=False, samescale=False, **kwds ):
        '''
        Plot the time series onto 3 separate graphs

           detrend=True to remove the trend from the plots
           samescale=True to force the X,Y,Z axes to share the same scale
           
           Additional keywords are passed to the pyplot.subplots function
           call.
        '''
        self._load()

        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates

        title="{0} {1} timeseries".format(self._code,self._solutiontype)
        if detrend:
            title=title+' detrended'

        settings={'sharex':True,'sharey':samescale,'figsize':(8,6),'dpi':100}
        settings.update(kwds)
        fig, plots=plt.subplots(3,1,**settings)
        fig.suptitle(title)
        data=self._data
        
        axis_labels=('East)','North','Up')
        for i,axis in enumerate(('e','n','u')):
            series=data[axis]
            ylabel=axis_labels[i]
            if detrend:
                days=mdates.date2num(data.index)
                trendp=np.poly1d(np.polyfit(days,series,1))
                trend=trendp(days)
                series=(series-trend)*1000
                ylabel=ylabel+' mm'
            plots[i].plot(data.index,series,'b+',label='Time series',picker=5)
            plots[i].set_ylabel(ylabel)
            plots[i].tick_params(labelsize=8)
            plots[i].format_xdata=mdates.DateFormatter('%d-%m-%Y')

    def compare( self, other, newcode=None, newtype=None ):
        '''
        Returns a time series comparing two others

        By default requires that both series have the same code.  Will
        return an error if not.  Over-ride by including a newcode parameter
        '''
        if newcode is None:
            if self._code != other._code:
                raise RuntimeError('Cannot compare two series with different codes')
            newcode=self._code
        if newtype is None:
            newtype=other._solutiontype+'-'+self._solutiontype
        d1=self.getData(enu=False)
        d2=other.getData(enu=False)
        join=d1.join(d2,rsuffix='2',how='inner')
        data=pd.DataFrame(data={'x':join.x2-join.x,'y':join.y2-join.y,'z':join.z2-join.z})
        xyz0=[0,0,0]
        xyzenu=self._xyzenu if self._xyzenu is not None else self._xyz0
        return Timeseries(newcode,newtype,data=data,xyz0=xyz0,xyzenu=xyzenu)

class SqliteTimeseries( Timeseries ):

    _sql='''
        select epoch, code as code, X as x, Y as y, Z as z 
        from mark_coordinate 
        where code=? and solution_type=?
        order by epoch
        '''

    _sqlMultiType='''
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
        order by m.epoch
        '''

    _sqlList='''
        select distinct code, solution_type
        from mark_coordinate
        order by code, solution_type
        '''

    @staticmethod
    def _openDb( dbfile ):
        if not os.path.exists(dbfile):
            raise RuntimeError(str(dbfile)+' does not exist')
        db=sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES)
        return db

    @staticmethod
    def seriesList( dbfile, solutiontype=None,after=None ):
        db=SqliteTimeseries._openDb( dbfile )
        seriescodes=pd.read_sql(SqliteTimeseries._sqlList, db )
        db.close()
        series=[]
        for i in seriescodes.index:
            code,solntype=(seriescodes.code[i],seriescodes.solution_type[i])
            if solutiontype is None or solutiontype == solntype:
                series.append(SqliteTimeseries(dbfile,code,solntype,after=after))
        return series

    def __init__( self, dbfile, code, solutiontype='default', xyz0=None, transform=None, after=None ):
        Timeseries.__init__( self, code, solutiontype=solutiontype, xyz0=xyz0, transform=transform, after=after )
        self._dbfile=dbfile

    def _loadData( self ):
        db=SqliteTimeseries._openDb( self._dbfile )
        solntype=self.solutiontype()
        code=self.code()
        if '+' not in solntype:
            data=pd.read_sql(SqliteTimeseries._sql,db,params=(code,solntype),index_col='epoch')
        else:
            types=solntype.split('+')
            casesql='CASE solution_type'
            for i,t in enumerate(types):
                casesql=casesql+' WHEN ? THEN '+str(i)
            casesql=casesql+' ELSE '+str(len(types))+' END'
            sql=SqliteTimeseries._sqlMultiType
            sql=sql.replace('{case}',casesql)
            params=[]
            params.extend(types)
            params.append(code)
            params.extend(types)
            params.append(len(types))
            params.append(code)
            params.extend(types)
            data=pd.read_sql(sql,db,params=params,index_col='epoch')
        db.close()
        #, parse_dates=['epoch'])
        data.set_index(pd.to_datetime(data.index),inplace=True)
        return data

class FileTimeseries( Timeseries ):

    '''
    Defines a CORS time series from a file for analysis.

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

    _columns=('name','epoch','x','y','z')
    _epoch_col='epoch'

    @staticmethod
    def seriesList( filepattern, solutiontype='default',after=None ):
        '''
        Get the potential time series files, any file matching the filepattern.  The 
        pattern can include {code} in the filename to represent the code to use
        '''
        path,file=os.path.split(filepattern)
        fileparts=file.split('{code}')
        if len(fileparts) != 2:
            raise RuntimeError('File timeseries file pattern must include {code}: '+filepattern)
        filere=re.compile(re.escape(fileparts[0])+r'(?P<code>\w\w\w\w)'+re.escape(fileparts[1])+'$')
        path=path or '.'
        series=[]
        for fn in os.listdir(path or '.'):
            match=filere.match(fn)
            if not match:
                continue
            if path:
                fn=os.path.join(path,fn)
            code=match.group('code')
            series.append(FileTimeseries(fn,code,solutiontype,after=after))
        return series


    def __init__( self, filename, code=None, solutiontype=None, xyz0=None, transform=None, after=None ):
        Timeseries.__init__(self,code,solutiontype=solutiontype,xyz0=xyz0,transform=transform,after=after)
        self._filename=filename

    def _loadData( self ):

        filename=self._filename
        code=self._code

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

        data.rename(columns={'name':'code'},inplace=True)

        code=code if code else data.code.iloc[0]
        self.setName(code)
        data=data[data.code==code]
        return data


    def filename( self ):
        return self._filename

class TimeseriesList( list ):

    def __init__( self, source=None, solutiontype=None, after=None ):
        list.__init__(self)
        if source is not None:
            if os.path.isfile(source):
                self.extend(SqliteTimeseries.seriesList(source, solutiontype, after ))
            else:
                self.extend(FileTimeseries.seriesList(source,solutiontype, after))

    def codes( self ):
        codes={}
        for f in self:
            codes[f.code()]=1
        return sorted(codes.keys())

    def solutiontypes( self ):
        solutiontypes={}
        for f in self:
            solutiontypes[f.solutiontype()]=1
        return sorted(solutiontypes.keys())

    def timeseries( self, code, solutiontype=[] ):
        '''
        Return the time series for the requested station codes.  Can
        specify a solution type wanted, or a list of solution types in
        order of preference.  If solutiontype is not specified then
        an error will be raised if there is more than one matching solution.
        '''
        if isinstance(solutiontype,basestring):
            solutiontype=[solutiontype] 
        potential=None
        priority=len(solutiontype)
        for series in self:
            if series.code() != code:
                continue
            if len(solutiontype) > 0:
                solntype=series.solutiontype()
                if solntype in solutiontype:
                    spriority=solutiontype.index(solntype)
                    if spriority < priority:
                        potential=series
                        priority=spriority
            elif potential is None:
                potential=series
            else:
                raise RuntimeError('Ambiguous time series requested - multiple solution types')
        if not potential:
            raise RuntimeError('No time series found for requested station '+code)
        return potential






