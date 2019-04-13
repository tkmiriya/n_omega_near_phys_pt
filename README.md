# n omega near phys pt

This is the repository of the published paper
"N-Omega dibaryon from lattice QCD near the physical point",
Phys. Lett. B792(2019)284-289, DOI:10.1016/j.physletb.20199.03.050
[arXiv:1810.03416 [hep-lat]](https://arxiv.org/abs/1810.03416).

* notes

  + N Omega potential and fit.ipynb

  plot NOmega(5S2) potential and fitting using gaussian + Yukawa-squared with form factor

  + R-correlator.ipynb

  plot time-dependence of the R-correlator

  + baryon masses.ipynb

  analyze effective masses

  + binding energy.ipynb

  calculate the binding energy and the root mean square distance from fitted parameter (Run "N Omega potential and fit.ipynb" before)

  + quark mass dep. of binding energies.ipynb

  estimate quark mass dependence of the binding energy (using physical baryon masses) (Run "N Omega potential and fit.ipynb" before)

  + scattering phase shift.ipynb

  calculate the scattering phase shift from fitted parameter (Run "N Omega potential and fit.ipynb" before)

* data

  + corr_jk_full_stat_20bins.pkl

  pickle file for the jack-knife samples of the single baryon correlator
  (See "baryon masses.ipynb")

  + n_omega_spinX_Rcorr_jk.pkl

  pickle file for the jack-knife samples of the R-correlator (See "R-correlator.ipynb")

  + n_omega_spinX_pot_av.pkl, n_omega_spinX_pot_jk.pkl

  pickle file for average and jack-knife samples of the potential (See "N Omega potential and fit.ipynb")

* src

  + n_omega_pot.py

  python code to calculate potentials from raw NBS files
