
function OrdinateTimeseries( ordinate, data, columnNames, onupdate )
{
    this.ordinate=ordinate;
    this.data=data;
    this.columnNames=columnNames;
    this.ftrend0=0;
    this.ftrend=0;
    this._onupdate=onupdate;
}

/** 
 * Add a column from another time series
 */

OrdinateTimeseries.prototype.addColumn=function( srcts,srccol, tgtcol )
{
    var srccolno=srcts.columnNames.indexOf(srccol);
    if( srccolno < 0 ) return;
    var tgtcolno=this.columnNames.indexOf(tgtcol);
    var nsrc=0;
    var maxsrc=srcts.data.length-1;
    var ndata=this.data.length;
    var data=this.data;
    var srcdata=srcts.data;
    var ft0=this.ftrend0-srcts.ftrend0;
    var ft=this.ftrend-srcts.ftrend;
    var tolerance=43200*1000; /* date matching tolerance in milliseconds */
    if( tgtcolno < 0 )
    {
        this.columnNames.push(tgtcol);
    }
    for( var i=0; i<ndata; i++ )
    {
        var epoch=data[i][0];
        while( nsrc <= maxsrc && (srcdata[nsrc][0]-epoch) < -tolerance ) nsrc++;
        var value=NaN;
        if( nsrc < maxsrc && (srcdata[nsrc][0]-epoch) < tolerance )
        {
            value=srcdata[nsrc][srccolno]-ft0-ft*(epoch-data[0][0]);
            nsrc++;
        }
        if( tgtcolno < 0 )
        {
            data[i].push(value);
        }
        else
        {
            data[i][tgtcolno]=value;
        }
    }
}

/**
 * Set a handler for updating of the timeseries
 */

OrdinateTimeseries.prototype.onUpdate=function( onupdate )
{
    this._onupdate=onupdate;
}

OrdinateTimeseries.prototype._removeTrend=function( ftrend0, ftrend )
{
    var data=this.data;
    var ndata=data.length;
    if( ndata < 1 ) return;
    if( ftrend0 == 0 && ftrend == 0 ) return;

    var nvar=data[0].length;
    for( var i=0; i<ndata; i++ )
    {
        for( var j=1; j < nvar; j++ )
        {
            data[i][j]=data[i][j]-ftrend0-ftrend*(data[i][0]-data[0][0]);
        }
    }
    this.ftrend0 += ftrend0;
    this.ftrend += ftrend;
    if( this._onupdate != null ){ this._onupdate(); }
}

/**
 *  Detrend the timeseries using the series specified
 *
 *  @param {integer} series the number of the series to use for detrending (1,...)
 *  @param {npoints} npoints the number of points from each end of the series to base the trend on
 *
 */

OrdinateTimeseries.prototype.detrend=function( series, npoints )
{
    var data=this.data;
    var total=data.length;
    var testpoints=[];

    for( var i=0; i<total; i++ )
    {
        if( ! isNaN(data[i][series] )) 
        {
            testpoints.push(i);
            if( testpoints.length >= npoints ) break;
        }
    }
    if( testpoints.length == 0 ) return;

    testpoints.sort(function(i,j){ return data[i][series]-data[j][series]; });
    var i0=testpoints[Math.round((testpoints.length-1)/2)]

    testpoints=[];
    for( var i=total; i-- >0; )
    {
        if( ! isNaN(data[i][series] )) 
        {
            testpoints.push(i);
            if( testpoints.length >= npoints ) break;
        }
    }

    testpoints.sort(function(i,j){ return data[i][series]-data[j][series]; });
    var i1=testpoints[Math.round((testpoints.length-1)/2)];

    if( i1 == i0 ) return;

    var ftrend=(data[i1][series]-data[i0][series])/(data[i1][0]-data[i0][0]);
    var ftrend0=(data[0][0]-data[i0][0])*ftrend+data[i0][series];

    this._removeTrend( ftrend0, ftrend );
}

/**
 *  Restore the original time series (ie cancel any detrending)
 *
 */

OrdinateTimeseries.prototype.retrend=function()
{
    if( 'ftrend0' in this )
    {
        this._removeTrend( -this.ftrend0, -this.ftrend );
    }
}

OrdinateTimeseries.prototype.applyTrend=function(ots )
{
    this.retrend();
    if( 'ftrend0' in ots )
    {
        this._removeTrend( ots.ftrend0, ots.ftrend );
    }
}


function CORS_timeseries(strategy)
{
    this.strategy=strategy.toUpperCase();
    this.ordinates=["east","north","up"];
    this.datacols=[]
    this.data={};
}


CORS_timeseries.prototype.parseTimeseriesCsv=function( data, prefix )
{
    var lines=data.split(/\n/);
    var header=lines.shift();
    var colnames=header.toLowerCase().split(/\,/).map(function(s){ return s.trim()});
    var datecol=colnames.indexOf('epoch');
    var ordinates=this.ordinates;
    var tsdata={};
    var ccols={};
    ordinates.forEach(function(ordinate){
        tsdata[ordinate]=[]
        c=ordinate.charAt(0);
        ccols[ordinate]=[c,"gdb_"+c,"scm_"+c].map(function(m)
                    {
                    return colnames.indexOf(m);
                    });
            });

    lines.forEach( function(line)
            {
                columns=line.split(/\,/);
                date=new Date(columns[datecol]);
                if( isNaN(date.getTime())) return;
                ordinates.forEach(function(ord)
                    {
                        cdata=[date].concat(
                            ccols[ord].map( function(ci)
                            {
                                if( columns[ci] == '' ) return NaN;
                                // return columns[ci];
                                return parseFloat(columns[ci]);
                            }));
                        tsdata[ord].push(cdata);
                    });

            });

    var corsdata=this;
    ordinates.forEach(function(ordinate){
        var labels=["Epoch",prefix+" "+ordinate,"GDB "+ordinate,"SCM "+ordinate];
        corsdata.data[ordinate]=new OrdinateTimeseries(ordinate,tsdata[ordinate],labels);
        });
    this.datacols=[prefix+" ","GDB ","SCM "];

    return tsdata;
}

CORS_timeseries.prototype.loadCsv=function( url,action )
{
    var corsdata=this;
    var prefix=this.strategy;
    $.ajax({
        url: url,
        success: function(data){
            tsdata=corsdata.parseTimeseriesCsv(data,prefix);
            action(corsdata);
            },
        dataType: "text"
        });
}

CORS_timeseries.prototype.addColumn = function( source, source_prefix, target_prefix )
{
    if( source.datacols.indexOf(source_prefix) < 0 )
    {
        return false;
    }
    if( this.datacols.indexOf(target_prefix) < 0 ) this.datacols.push(target_prefix);
    for( ord in this.data )
    {
        var srcts=source.data[ord];
        var tgtts=this.data[ord];
        var srccol=source_prefix+ord;
        var tgtcol=target_prefix+ord;
        tgtts.addColumn(srcts,srccol,tgtcol);
    }
}

CORS_timeseries.prototype.detrend=function( series, npoints )
{
    for( ord in this.data )
    {
        this.data[ord].detrend(series,npoints);
    }
}


CORS_timeseries.prototype.retrend=function()
{
    for( ord in this.data )
    {
        this.data[ord].retrend();
    }
}

CORS_timeseries.prototype.applyTrend=function( corsts )
{
    for( ord in this.data )
    {
        this.data[ord].applyTrend( cordts.data[ord] );
    }
}
