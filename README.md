# Xenon

Contents:

- 8B_flux: solar 8B neutrino flux from Bahcall
- core: core functions. Contains the differentntial cross-section, density matrix, and differential event rate functions
- FFM: computes the form factor for the differential cross section, from 2007.08529
- interpolators: just some basic first and second degree interpolators
- main: a nice and tidy main

The main function is the differential rate one, dNdT. To break it down a little:

double dNdT(double T, double norm, double flux_err, double flux_pull, double theta_13, double theta_12, double theta_23, double Dm2_21, double (*eps_e)[3], double (*eps_u)[3], double (*eps_d)[3])

It yields the recoil rate in keV^-1. You need to input the recoil energy T in keV, the normalization, in number of targets*seconds, the flux error (a 16% from Bahcall) and the pull (as the error enters in a normalization of the flux as (1.0 + error*pull)*flux), the oscillation parameters, and three NSI matrices, defined as in 2302.12846, equation 2. NSIs enter both in the density matrix (which would be the Sun part, NSIs in propagation though it) and the cross-sections (which would be the detection part via CEvNS).

