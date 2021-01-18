var strategies=['nzg'];
var procsum='Processing summary';
var procfiletemplates={
    nzg: {
        'Processing summary': '/results/nzg/{year}/{doy}/R2S{yy}{doy}0.PRC',
        'Bernese coordinate file': '/results/nzg/{year}/{doy}/F1_{yy}{doy}0.CRD',
        'Baselines summary file': '/results/nzg/{year}/{doy}/BSL{yy}{doy}0.BSL',
        'Final sinex file': '/results/nzg/{year}/{doy}/F1_{yy}{doy}0.SNX.gz'
    }
}

var resultdirtemplate='/analysis/{strategy}/';
var summaryfiletemplate=resultdirtemplate+'cors_summary.json';
var csvfiletemplate=resultdirtemplate+'{code}_enu_ts.csv';

// Web pages

var statusurl='index.html?strategy={strategy}';
var ploturl='plot.html?strategy={strategy}&code={code}';
var processingurl='processing.html?strategy={strategy}&year={year}&doy={doy}';
