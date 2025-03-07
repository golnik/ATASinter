var au2fs = 0.024188843265864002;
var au2eV = 27.211386245988;
var c = 137.;

var tmin = 0./au2fs;
var tmax = 100./au2fs;
var Nt = 100;
var dt = (tmax - tmin) / (Nt - 1);

var tgrid = new Array(Nt);

for(var it = 0; it < Nt; it++){
  tgrid[it] = tmin + dt * it;
}

var Emin = 21.4/au2eV;
var Emax = 21.75/au2eV;
var NE = 100;
var dE = (Emax - Emin) / (NE - 1);

var Egrid = new Array(NE);

for(var iE = 0; iE < NE; iE++){
  Egrid[iE] = Emin + dE * iE;
}

var NI = 2;
var EI = [
  0.00000000/au2eV,
  0.17749368/au2eV
];

var NF = 2;
var EF = [
  21.62407061/au2eV,
  21.67504277/au2eV
];

var G = 0.02/au2eV;

function compute_cs(t, w, D){
  var res = 0.;

  for(var I = 0; I < NI; I++){
    for(var J = 0; J < NI; J++){
      var i = math.complex(0.,1.);
      var phi = (EI[I] - EI[J]) * t;
      var xiIJ = math.exp( math.multiply(i,phi) );

      for(var F = 0; F < NF; F++){
        var dIF = D[I][F];
        var dFI = D[J][F];

        var EFI = math.divide(1.,math.complex(EF[F] - EI[I] - w,-0.5*G));
        var EFJ = math.divide(1.,math.complex(EF[F] - EI[J] + w,+0.5*G));

        var tmp = math.multiply(xiIJ,math.multiply(dIF*dFI,math.add(EFI,EFJ)));
        res += math.im(tmp);
      }
    }
  }

  return (4 * Math.PI * w / c) * res;
}

function compute_spec(D){
  //2d array of spectra
  var sigma = new Array(Nt);
  for(var it = 0; it < Nt; it++) {
    sigma[it] = new Array(NE);
    var t = tgrid[it];
    for(var iE = 0; iE < NE; iE++) {
      var w = Egrid[iE];
      sigma[it][iE] = compute_cs(t,w,D);
    }
  }
  return sigma;
}

function compute_lineouts(D){
  //2d array of spectra
  var sigma = new Array(NI*NF);
  var indx = 0;
  for(var i = 0; i < NI; i++) {
    for(var f = 0; f < NF; f++) {
      sigma[indx] = new Array(Nt);
      for(var it = 0; it < Nt; it++) {
        var t = tgrid[it];
        var w = EF[f] - EI[i];
        sigma[indx][it] = compute_cs(t,w,D);
      }
      indx+=1;
    }
  }

  //normalize to 1
  for(var indx = 0; indx < NI*NF; indx++) {
    var sigma_min = Math.min(...sigma[indx].flat());
    var sigma_max = Math.max(...sigma[indx].flat());
    for(var it = 0; it < Nt; it++) {
      sigma[indx][it] = (sigma[indx][it] - sigma_min) / (sigma_max - sigma_min) - 0.5;
    }
  }

  return sigma;
}

var tgrid_fs = tgrid.map(function(t){ return t * au2fs; });
var Egrid_eV = Egrid.map(function(E){ return E * au2eV; });

var global_margin = {
  l: 100,
  r: 50,
  b: 80,
  t: 25,
  pad: 0
}

var global_font = {
  family: '"Open Sans", verdana, arial, sans-serif',
  size: 20,
  color: '#444'
}


var layout = {
  showlegend: false,
  font: global_font,
  margin: global_margin,
  displayModeBar: false,
  responsive: true,
  xaxis: {
    title: 'Time [fs]',
  },
  yaxis: {
    title: {
      text: 'Energy [eV]',
      standoff: 30
    }
  },
};

var lineout_layout = {
  showlegend: false,
  font: global_font,
  margin: global_margin,
  displayModeBar: false,
  responsive: true,
  yaxis: {
    range: [Egrid_eV[0],Egrid_eV[NE-1]]
  },
  xaxis: {
    title: 'Time [fs]',
  }
};

var colors = ["rgb(255, 0, 0)","rgb(0, 0, 255)","rgb(0, 150, 0)","rgb(255, 0, 255)"];

function plot_spec(D1I1F,D1I2F,D2I1F,D2I2F) {
  var D = [
    [D1I1F, D1I2F],
    [D2I1F, D2I2F]
  ];

  var sigma    = compute_spec(D);
  var lineouts = compute_lineouts(D);

  var data = [ {
      z: sigma,
      x: tgrid_fs,
      y: Egrid_eV,
      type: 'contour',
      transpose: true,
      showlines: true,
      colorscale: 'Jet',
      ncontours: 100,
      line:{
        width:0
      }
    }
  ];
  
  Plotly.newPlot('ATAS', data, layout);

  var traces = [];
  var indx = 0;
  for(var i = 0; i<NI; i++){
    for(var f = 0; f<NF; f++){
      var factor = D[i][f];
      var w = (EF[f] - EI[i]) * au2eV;
      var lineout = lineouts[indx].map(function(v){ return 0.08*Math.abs(factor)*v + w; });

      var trace = {
        x: tgrid_fs,
        y: lineout,
        mode: 'lines',
        line: {
          color: colors[indx],
          width: 3
        }        
      };
      traces.push(trace);

      indx++;
    }
  }

  Plotly.newPlot('lineout', traces, lineout_layout);
};