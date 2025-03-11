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

var global_sigma_min = 0;
var global_sigma_max = 0;