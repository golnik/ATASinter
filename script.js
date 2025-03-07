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
      var phi = -(EI[I] - EI[J]) * t;
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

var tgrid_fs = tgrid.map(function(t){ return t * au2fs; });
var Egrid_eV = Egrid.map(function(E){ return E * au2eV; });

var global_margin = {
  l: 150,
  r: 15,
  b: 100,
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
    title: 'Energy [eV]',
  },
};

function plot_spec(D1I1F,D1I2F,D2I1F,D2I2F) {
  var D = [
    [D1I1F, D1I2F],
    [D2I1F, D2I2F]
  ];

  sigma = compute_spec(D);

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
  Plotly.newPlot('myDiv', data, layout);
};