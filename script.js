var colors = [
  "rgb(255, 0, 0)",
  "rgb(0, 0, 255)",  
  "rgb(0, 150, 0)",
  "rgb(255, 0, 255)"
];

document.getElementById("sD1").style.backgroundColor = colors[0];
document.getElementById("sD2").style.backgroundColor = colors[1];
document.getElementById("sD3").style.backgroundColor = colors[2];
document.getElementById("sD4").style.backgroundColor = colors[3];

//important, colors in the spectra appear in energy order
var colors_energy = [
  colors[2],
  colors[3],
  colors[0],
  colors[1]
];

var slider = document.getElementById('slider');

wI1F1 = (EF[0] - EI[0])*au2eV;
wI1F2 = (EF[1] - EI[0])*au2eV;
wI2F1 = (EF[0] - EI[1])*au2eV;
wI2F2 = (EF[1] - EI[1])*au2eV;

noUiSlider.create(slider, {
  range: {
    'min': Emin*au2eV,
    'max': Emax*au2eV
  },

  step: 0.001/au2eV,

  // Handles start at ...
  start: [wI2F1,wI2F2,wI1F1,wI1F2],
  tooltips: wNumb({decimals: 3}),

  // Put '0' at the bottom of the slider
  direction: 'rtl',
  orientation: 'vertical'
});

var i = math.complex(0.,1.);

function compute_cs(t, w, D, c1){
  var res = 0.;
  var CI = [Math.sqrt(c1),Math.sqrt(1.-c1)];

  for(var I = 0; I < NI; I++){
    var cI = CI[I];
    for(var J = 0; J < NI; J++){
      var cJ = CI[J];

      var phi = (EI[I] - EI[J]) * t;
      var cIJ = math.multiply(cI,cJ);
      var xiIJ = math.multiply(cIJ, math.exp( math.multiply(i,phi) ));

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

function compute_spec(D,cI){
  //2d array of spectra
  var sigma = new Array(Nt);
  for(var it = 0; it < Nt; it++) {
    sigma[it] = new Array(NE);
    var t = tgrid[it];
    for(var iE = 0; iE < NE; iE++) {
      var w = Egrid[iE];
      sigma[it][iE] = compute_cs(t,w,D,cI);
    }
  }

  global_sigma_min = Math.min(...sigma.flat());
  global_sigma_max = Math.max(...sigma.flat());

  return sigma;
}

function compute_lineouts(D,cI,energies){
  var size = energies.length;

  var sigma = new Array(size);
  for(var indx = 0; indx < size; indx++) {
    sigma[indx] = new Array(Nt);
    var w = energies[indx];
    for(var it = 0; it < Nt; it++) {
      var t = tgrid[it];
      sigma[indx][it] = compute_cs(t,w,D,cI);
    }
  }

  //normalize to 1
  for(var indx = 0; indx < size; indx++) {
    //var sigma_min = Math.min(...sigma[indx].flat());
    //var sigma_max = Math.max(...sigma[indx].flat());
    for(var it = 0; it < Nt; it++) {
      sigma[indx][it] = (sigma[indx][it] - global_sigma_min) / (global_sigma_max - global_sigma_min);
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
    title: {
      text: 'Time [fs]'
    }
  },
  yaxis: {
    title: {
      text: 'Energy [eV]',
      standoff: 30
    }
  },

  shapes: []   
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
    title: {
      text: 'Time [fs]'
    }
  },

  shapes: []
};

// Appending money-formatting
var Format = wNumb({
	decimals: 1,
});

function get_parameters() {
  var D1I1F = Format.from(document.getElementById("sD1").value);
  var D1I2F = Format.from(document.getElementById("sD2").value);
  var D2I1F = Format.from(document.getElementById("sD3").value);
  var D2I2F = Format.from(document.getElementById("sD4").value);
  var CI = Format.from(document.getElementById("sCI").value);

  document.getElementById("vD1").innerHTML = Format.to(D1I1F);
  document.getElementById("vD2").innerHTML = Format.to(D1I2F);
  document.getElementById("vD3").innerHTML = Format.to(D2I1F);
  document.getElementById("vD4").innerHTML = Format.to(D2I2F);
  document.getElementById("vCI").innerHTML = Format.to(CI);

  var D = [
    [D1I1F, D1I2F],
    [D2I1F, D2I2F]
  ];

  return [D,CI];
}

const zeros = (m, n) => [...Array(m)].map(e => Array(n).fill(0));

var data = [ {
    z: zeros(NI,NF),
    x: tgrid_fs,
    y: Egrid_eV,
    type: 'contour',
    transpose: true,
    showlines: true,
    colorscale: 'Jet',
    ncontours: 100,
    line:{
      width:0
    },
    colorbar: {
      orientation: 'v'
    },
    showscale: false,  //hide colorbar    
  }
];

Plotly.newPlot('ATAS', data, layout);
Plotly.newPlot('lineout', [], lineout_layout);

function plot_ATAS() {
  var [D,cI] = get_parameters();
  var sigma = compute_spec(D,cI);
  data[0].z = sigma;
  Plotly.react('ATAS', data, layout);
};

function plot_lineouts() {
  var [D,cI] = get_parameters();

  var lineout_energies_eV = slider.noUiSlider.get(true);
  var lineout_energies = lineout_energies_eV.map(function(E){ return E/au2eV; });
  var nlineouts = lineout_energies.length;

  var lineouts = compute_lineouts(D,cI,lineout_energies);

  var traces = [];
  var shapes = [];
  for(var indx = 0; indx<nlineouts; indx++){
    var w = lineout_energies_eV[indx];

    var min = Math.min(...lineouts[indx].flat());
    var max = Math.max(...lineouts[indx].flat());
    var shift = 0.5 * (max - min) + min;
    var lineout = lineouts[indx].map(function(v){ return 0.05*(v-shift) + w; });

    var trace = {
      x: tgrid_fs,
      y: lineout,
      mode: 'lines',
      line: {
        color: colors_energy[indx],
        width: 3
      }        
    };
    traces.push(trace);

    var shape = {
      type: 'line',
      x0: 80,
      y0: lineout_energies_eV[indx],
      x1: 100,
      y1: lineout_energies_eV[indx],
      line:{
        color: colors_energy[indx],
        width: 4,
      }
    };
    shapes.push(shape);
  }

  layout.shapes = shapes;
  //lineout_layout.shapes = shapes;

  Plotly.relayout('ATAS', layout);
  Plotly.react('lineout', traces, lineout_layout);
};

function plot() {
  plot_ATAS();
  plot_lineouts();
}

// Binding signature
slider.noUiSlider.on('update',plot_lineouts);