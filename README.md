# Attosecond Transient Absorption Spectroscopy (ATAS) Learning Tool

An interactive, in-browser visualization of attosecond transient absorption spectroscopy (ATAS), used to probe coherent electron dynamics in a quantum system prepared in a superposition of its quantum states. Built as a learning tool for exploring the technique through live, adjustable plots rather than static figures.

**[Live demo](https://atas.ngolubev.com/)**

## What it covers

The interface is organized around the model diagram and two linked plots:

- **System Energy Levels** - the two initial states, \(\Phi_1\) and \(\Phi_2\), and the two higher-lying final states, \(\Phi_{F_1}\) and \(\Phi_{F_2}\), reached upon absorption of the probe photon.
- **Transition Dipole Moments & Population** - sliders controlling the four transition dipole matrix elements \(\langle \Phi_I | \hat{\mu} | \Phi_F \rangle\) and the population weight \(c_1^2\) of the initial superposition.
- **Attosecond Transient Absorption Spectrum** - the two-dimensional spectrum \(\sigma(\omega,\tau)\) as a function of pump-probe delay and photon energy.
- **Energy-Resolved Lineouts** - time traces of the spectrum at four selectable photon energies, matched by color to the spectrum and the energy-level diagram.

Each plot has an accompanying info icon with a plain-language explanation of what's being shown.

## Theory & references

The underlying ATAS theory implemented here is developed in:

- R. Santra *et al.*, [Phys. Rev. A **83**, 033405 (2011)](https://doi.org/10.1103/PhysRevA.83.033405)
- N. Golubev *et al.*, [Phys. Rev. Lett. **127**, 123001 (2021)](https://doi.org/10.1103/PhysRevLett.127.123001)

## Running locally

This is a static site with no build step. Serve the directory with any static file server and open it in a browser, for example:

```bash
python3 -m http.server
```

Then visit `http://localhost:8000`.

## Tech stack

- Vanilla HTML/CSS/JS (no framework)
- [Plotly.js](https://plotly.com/javascript/) for plotting
- [MathJax](https://www.mathjax.org/) for math rendering
- [math.js](https://mathjs.org/) for complex-number arithmetic
- [noUiSlider](https://refreshless.com/nouislider/) / [wNumb](https://refreshless.com/wnumb/) for the energy-selection slider

## Developers

Nikolay Golubev

## License

[GPL-3.0](LICENSE)
