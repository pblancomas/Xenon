#ifndef __FFM_H
#define __FFM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double FplusM(double q, int Z);

double FminusM(double q, int Z);

double Fplusphi(double q, int Z);

double Fminusphi(double q, int Z);

double FpM(double q, int Z);

double FnM(double q, int Z);

double Fpphi(double q, int Z);

double Fnphi(double q, int Z);

int CuV(double (*out_matrix)[3], double (*eps_uV)[3]);

int CdV(double (*out_matrix)[3], double (*eps_dV)[3]);

int CsV(double (*out_matrix)[3], double (*eps_sV)[3]);

int gvp(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]);

int gvn(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]);

int gvb(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]);

int gvpp(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN, double mN);

int gv2p(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double kp, double kn, double ksN);

int gvnp(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN, double mN);

int gv2n(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double kp, double kn, double ksN);


int QwFw_BSM(double (*out_matrix)[3], double T, int Z, double mN, double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN);



#endif
