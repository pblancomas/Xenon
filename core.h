#ifndef __CORE_H
#define __CORE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double eps_D(double (*eps_input)[3], double theta_13, double theta_12, double theta_23);

double eps_N(double (*eps_input)[3], double theta_13, double theta_12, double theta_23);

double sol_density(double x0, int type);

double p(double x0, double (*eps_e)[3], double (*eps_p)[3], double (*eps_n)[3], double Enu, double theta_13, double theta_12, double theta_23, double Dm2_21);

double q(double x0, double (*eps_e)[3], double (*eps_p)[3], double (*eps_n)[3], double Enu, double theta_13, double theta_12, double theta_23, double Dm2_21);

double f_p_8B(double x0);

int compute_density_matrix(double (*output_matrix)[3], double (*eps_e)[3], double (*eps_u)[3], double (*eps_d)[3], double Enu, double theta_13, double theta_12, double theta_23, double Dm2_21);


int compute_dsigmadT(double (*out_matrix)[3], double T, double Enu, int Z, double mN, double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]);

double dNdT(double T, double norm, double flux_err, double flux_pull, double theta_13, double theta_12, double theta_23, double Dm2_21, double (*eps_e)[3], double (*eps_u)[3], double (*eps_d)[3]);

#endif
