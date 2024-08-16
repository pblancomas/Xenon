#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FFM.h"


/////////////////////////////////
// form factor from 2007.08529 //
/////////////////////////////////

double FplusM(double q, int Z) {
    //q in GeV
    
    if (Z == 55){  //Cs
    double coefs[6] = {133,-134.2,38.9577,-4.12938,0.151119,-0.00103353};  
    double b = 2.2976; 
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 6; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    else if (Z == 53){  //I
    double coefs[6] = {127,-125.164,35.3993,-3.62687,0.125083,-0.000670162};  
    double b = 2.2821; 
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 6; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    else if (Z == 54){  //Xe
    
    //Xenon relative abundance array
    double abund[7] = {0.019102, 0.264006, 0.040710, 0.212324, 0.269086, 0.104357, 0.088573};
    
    //b vector
    double b_vec[7] = {2.2847, 2.2873, 2.2899, 2.2925, 2.2950, 2.3001, 2.3051};
    
    // Array 2D 
    double coefs[6][7] = {
    {128.0,     129.0,   130.0,    131.0,   132.0,    134.0,    136.0},
    {-126.455, -128.09, -129.753, -131.26, -132.835, -135.861, -138.787},
    {35.82, 36.4367, 37.2381, 37.8232, 38.4665, 39.6872, 40.9048},
    {-3.66991, -3.75317, -3.82921, -3.97171, -4.06999, -4.24713, -4.41984},
    {0.125062, 0.129553, 0.139778, 0.142995, 0.149636, 0.159053, 0.165388},
    {-5.63731e-04, -6.55816e-04, -9.30032e-04, -9.12955e-04, -0.00111463, -0.00125724, -0.00109211}
    };
    
    double resultado = 0.;
    double u;
    for (int j = 0; j < 7; j++){  //Xe isotope loop
    u = q * q * (b_vec[j] * 1000 / 197.327) * (b_vec[j] * 1000 / 197.327) / 2;
    
    for (int i = 0; i < 6; i++) {
        resultado += abund[j]*exp(-u / 2)*coefs[i][j] * pow(u, i);
    }
    }
    
    return resultado;
}
    
    
}

double FminusM(double q, int Z) {

    if (Z == 55){  //Cs
    double coefs[6] = {-23,33.9495,-13.9502,2.04567,-0.102733,0.000944352};  //Cs
    double b = 2.2976; 
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 6; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    else if (Z == 53){  //I
    double coefs[6] = {-21,30.4307,-12.321,1.78131,-0.0870947,0.000697815};  
    double b = 2.2821; 
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 6; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    
    else if (Z == 54){  //Xe
    
    //Xenon relative abundance array
    double abund[7] = {0.019102, 0.264006, 0.040710, 0.212324, 0.269086, 0.104357, 0.088573};
    
    //b vector
    double b_vec[7] = {2.2847, 2.2873, 2.2899, 2.2925, 2.2950, 2.3001, 2.3051};
    
    // Array 2D 
    double coefs[6][7] = {
    {108.0 - 128.0, 108.0 - 129.0, 108.0 - 130.0, 108.0 - 131.0,108.0 - 132.0, 108.0 - 134.0,  108.0 - 136.0},
    {29.0588, 30.6854, 32.2019, 33.7021, 35.253, 38.2701, 41.2081},
    {-11.7104, -12.3687, -13.1152, -13.7433, -14.4437, -15.773, -17.0848},
    {1.68447, 1.77928, 1.90775, 2.00031, 2.11305, 2.32061, 2.52635},
    {-0.0820044, -0.0868754, -0.0948184, -0.0991364, -0.105689, -0.116557, -0.12686},
    {6.65781e-04, 7.39474e-04, 8.47975e-04, 8.60686e-04, 9.61344e-04, 0.00106693, 0.00110965}
};

    
    double resultado = 0.;
    double u;
    for (int j = 0; j < 7; j++){  //Xe isotope loop
    u = q * q * (b_vec[j] * 1000 / 197.327) * (b_vec[j] * 1000 / 197.327) / 2;
    
    for (int i = 0; i < 6; i++) {
        resultado += abund[j]*exp(-u / 2)*coefs[i][j] * pow(u, i);
    }
    }
    
    return resultado;
}
    
}

double Fplusphi(double q, int Z) {
    
    if (Z == 55){  //Cs
    double coefs[6] = {-28.2527,20.4868,-4.09303,0.275572,-0.0051254};  
    double b = 2.2976;
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 5; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
     
    }
    else if (Z == 53){  //I
    double coefs[6] = {-26.1218,18.1692,-3.50413,0.223523,-0.00360552};  
    double b = 2.2821; 
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 5; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    else if (Z == 54){  //Xe
    
    //Xenon relative abundance array
    double abund[7] = {0.019102, 0.264006, 0.040710, 0.212324, 0.269086, 0.104357, 0.088573};
    
    //b vector
    double b_vec[7] = {2.2847, 2.2873, 2.2899, 2.2925, 2.2950, 2.3001, 2.3051};
    
    // Array 2D 
    double coefs[5][7] = {
    {-25.211, -26.1264, -27.7106, -28.0443, -28.7972, -29.5095, -29.8571},
    {17.592, 18.4401, 19.7108, 20.0888, 20.7751, 21.5578, 22.0402},
    {-3.46466, -3.64669, -3.85805, -3.94934, -4.0995, -4.27308, -4.37033},
    {0.224722, 0.239379, 0.252667, 0.260624, 0.272865, 0.287393, 0.296134},
    {-0.00353316, -0.00399779, -0.00444209, -0.00468846, -0.00507527, -0.00555437, -0.0059684}
};


    
    double resultado = 0.;
    double u;
    for (int j = 0; j < 7; j++){  //Xe isotope loop
    u = q * q * (b_vec[j] * 1000 / 197.327) * (b_vec[j] * 1000 / 197.327) / 2;
    
    for (int i = 0; i < 5; i++) {
        resultado += abund[j]*exp(-u / 2)*coefs[i][j] * pow(u, i);
    }
    }
    
    return resultado;
}
}

double Fminusphi(double q, int Z) {
    
    if (Z == 55){  //Cs
    double coefs[5] = {8.98993,-8.67714,2.21868,-0.189453,0.00473947};  //Cs
    double b = 2.2976; //Cs
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 5; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    else if (Z == 53){  //I
    double coefs[5] = {3.58476,-4.58091,1.46191,-0.139708,0.0035109};  
    double b = 2.2821; 
    
    double u = q * q * (b * 1000 / 197.327) * (b * 1000 / 197.327) / 2;
    double resultado = 0;
    int i;
    for (i = 0; i < 5; i++) {
        resultado = resultado + coefs[i] * pow(u, i);
    }
    return exp(-u / 2) * resultado;
    
    }
    else if (Z == 54){  //Xe
    
    //Xenon relative abundance array
    double abund[7] = {0.019102, 0.264006, 0.040710, 0.212324, 0.269086, 0.104357, 0.088573};
    
    //b vector
    double b_vec[7] = {2.2847, 2.2873, 2.2899, 2.2925, 2.2950, 2.3001, 2.3051};
    
    // Array 2D 
    double coefs[5][7] = {
    {3.89629, 5.47022, 6.28519, 6.90542, 7.93145, 9.3351, 10.1433},
    {-4.73163, -5.96963, -6.63842, -7.17962, -8.01086, -9.20279, -9.96123},
    {1.48489, 1.7533, 1.85406, 1.97217, 2.12817, 2.35489, 2.48784},
    {-0.140203, -0.160094, -0.166079, -0.175248, -0.186148, -0.202364, -0.212062},
    {0.00344765, 0.00387983, 0.00413453, 0.00437613, 0.00469887, 0.00519463, 0.00559688}
    };

    
    double resultado = 0.;
    double u;
    for (int j = 0; j < 7; j++){  //Xe isotope loop
    u = q * q * (b_vec[j] * 1000 / 197.327) * (b_vec[j] * 1000 / 197.327) / 2;
    
    for (int i = 0; i < 5; i++) {
        resultado += abund[j]*exp(-u / 2)*coefs[i][j] * pow(u, i);
    }
    }
    
    return resultado;
}
}

double FpM(double q, int Z) {
    double result_plusM = (FplusM(q, Z) + FminusM(q, Z)) / 2;
    return result_plusM;
}

double FnM(double q, int Z) {
    double result_minusM = (FplusM(q, Z) - FminusM(q, Z)) / 2;
    return result_minusM;
}

double Fpphi(double q, int Z) {
    double result_plusphi = (Fplusphi(q, Z) + Fminusphi(q, Z)) / 2;
    return result_plusphi;
}

double Fnphi(double q, int Z) {
    double result_minusphi = (Fplusphi(q, Z) - Fminusphi(q, Z)) / 2;
    return result_minusphi;
}

int CuV(double (*out_matrix)[3], double (*eps_uV)[3]) {

    double GF = 1.1663787e-5; // Constante de Fermi en GeV⁻²
    double s2thw = 0.23122; // Sin cuadrado del ángulo de mezcla débil
    
    double delta;
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    delta = (i == j) ? 1.0 : 0.0;
    
    out_matrix[i][j] = -delta*GF*(1 - (8.0 / 3)*s2thw)/sqrt(2) - sqrt(2)*GF*eps_uV[i][j];
    
    }
    }
    
    return 0;
}

int CdV(double (*out_matrix)[3], double (*eps_dV)[3]) {
    double GF = 1.1663787e-5; // Constante de Fermi en GeV⁻²
    double s2thw = 0.23122; // Sin cuadrado del ángulo de mezcla débil
    
    double delta;
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    delta = (i == j) ? 1.0 : 0.0;
    
    out_matrix[i][j] = delta*GF*(1 - (4.0 / 3)*s2thw)/sqrt(2) - sqrt(2)*GF*eps_dV[i][j];
    
    //printf("eps_d[%i][%i]: %g\n", i, j, eps_dV[i][j] );
    }
    }
    
    
    return 0;
}

int CsV(double (*out_matrix)[3], double (*eps_sV)[3]) {
    double GF = 1.1663787e-5; // Constante de Fermi en GeV⁻²
    double s2thw = 0.23122; // Sin cuadrado del ángulo de mezcla débil
    
    double delta;
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    delta = (i == j) ? 1.0 : 0.0;
    
    out_matrix[i][j] = delta*GF*(1 - (4.0 / 3)*s2thw)/sqrt(2) - sqrt(2)*GF*eps_sV[i][j];
    
    }
    }
    
    return 0;
}


int gvp(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]){

    double mat_u[3][3], mat_d[3][3];
    int up = CuV(mat_u, eps_uV);
    int down = CdV(mat_d, eps_dV) ;
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = 2*mat_u[i][j] + mat_d[i][j];
    
    }
    }

    return 0;
}

int gvn(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]) {

    double mat_u[3][3], mat_d[3][3];
    int up = CuV(mat_u, eps_uV);
    int down = CdV(mat_d, eps_dV) ;
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = mat_u[i][j] + 2*mat_d[i][j];
    
    }
    }

    return 0;
}

int gvb(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]) {

    double mat_u[3][3], mat_d[3][3], mat_s[3][3];
    int up = CuV(mat_u, eps_uV);
    int down = CdV(mat_d, eps_dV) ;
    int strange = CsV(mat_s, eps_sV) ;
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = mat_u[i][j] + mat_d[i][j] + mat_s[i][j];
    
    }
    }

    return 0;
}

int gvpp(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN, double mN) {
    
    double gvp_mat[3][3], gvn_mat[3][3], gvb_mat[3][3];
    
    gvp(gvp_mat, eps_uV, eps_dV, eps_sV);
    gvn(gvn_mat, eps_uV, eps_dV, eps_sV);
    gvb(gvb_mat, eps_uV, eps_dV, eps_sV);
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = gvp_mat[i][j]*(r2Ep/6 - kp/(4*mN*mN)) + gvn_mat[i][j]*(r2En/6 - kn/(4*mN*mN)) + gvb_mat[i][j]*(r2EsN/6 - ksN/(4*mN*mN));
    
    }
    }
    
    return 0;
}

int gv2p(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double kp, double kn, double ksN) {

    double gvp_mat[3][3], gvn_mat[3][3], gvb_mat[3][3];
    
    gvp(gvp_mat, eps_uV, eps_dV, eps_sV);
    gvn(gvn_mat, eps_uV, eps_dV, eps_sV);
    gvb(gvb_mat, eps_uV, eps_dV, eps_sV);
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = gvp_mat[i][j]*kp + gvn_mat[i][j]*kn + gvb_mat[i][j]*ksN;
    
    }
    }
    
    return 0;
}

int gvnp(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN, double mN) {
    
    double gvp_mat[3][3], gvn_mat[3][3], gvb_mat[3][3];
    
    gvp(gvp_mat, eps_uV, eps_dV, eps_sV);
    gvn(gvn_mat, eps_uV, eps_dV, eps_sV);
    gvb(gvb_mat, eps_uV, eps_dV, eps_sV);
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = gvn_mat[i][j]*(r2Ep/6 - kp/(4*mN*mN)) + gvp_mat[i][j]*(r2En/6 - kn/(4*mN*mN)) + gvb_mat[i][j]*(r2EsN/6 - ksN/(4*mN*mN));
    
    }
    }
    
    return 0;
}

int gv2n(double (*out_matrix)[3], double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double kp, double kn, double ksN) {

    double gvp_mat[3][3], gvn_mat[3][3], gvb_mat[3][3];
    
    gvp(gvp_mat, eps_uV, eps_dV, eps_sV);
    gvn(gvn_mat, eps_uV, eps_dV, eps_sV);
    gvb(gvb_mat, eps_uV, eps_dV, eps_sV);
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    out_matrix[i][j] = gvn_mat[i][j]*kp + gvp_mat[i][j]*kn + gvb_mat[i][j]*ksN;
    
    }
    }
    
    return 0;
}


int QwFw_BSM(double (*out_matrix)[3], double T, int Z, double mN, double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3], double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN) {
//T: recoil energy (in GeV). All inputs in GeV related units
//output in GeV^-2! (comes from the g's)
    
    double q = sqrt(2*mN*T); //momentum transfer
    double t = q * q;
    
    double gvp_mat[3][3], gvn_mat[3][3], gvpp_mat[3][3], gv2p_mat[3][3], gvnp_mat[3][3], gv2n_mat[3][3];
    
    gvp(gvp_mat, eps_uV, eps_dV, eps_sV);
    gvn(gvn_mat, eps_uV, eps_dV, eps_sV);
    gvpp(gvpp_mat, eps_uV, eps_dV, eps_sV, r2Ep, r2En, r2EsN, kp, kn, ksN, mN);
    gv2p(gv2p_mat, eps_uV, eps_dV, eps_sV, kp, kn, ksN);
    gvnp(gvnp_mat, eps_uV, eps_dV, eps_sV, r2Ep, r2En, r2EsN, kp, kn, ksN, mN);
    gv2n(gv2n_mat, eps_uV, eps_dV, eps_sV, kp, kn, ksN);
    
    double p1, p2, p3, p4;
    
    for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
    
    p1 = gvp_mat[i][j] + gvpp_mat[i][j]* t + ((gvp_mat[i][j] + 2 * gv2p_mat[i][j]) * t / (8 * mN * mN));
    p2 = gvn_mat[i][j] + gvnp_mat[i][j] * t + ((gvn_mat[i][j] + 2 * gv2n_mat[i][j]) * t / (8 * mN * mN));
    p3 = (gvp_mat[i][j] + 2 * gv2p_mat[i][j]) * t / (4 * mN * mN);
    p4 = (gvn_mat[i][j] + 2 * gv2n_mat[i][j]) * t / (4 * mN * mN);
    
    out_matrix[i][j] = (p1 * FpM(q, Z) + p2 * FnM(q, Z) - p3 * Fpphi(q, Z) - p4 * Fnphi(q, Z));
    
    }
    }
    
    return 0;
}

/*
double Qw_BSM_fun(int Z, int N, double mN, double eps_uV, double eps_dV, double eps_sV, double r2Ep, double r2En, double r2EsN, double kp, double kn, double ksN) {
//all inputs in GeV related units
//output adimensional
    
    double GF = 1.16637*pow(10,-11);  //Fermi constant, in MeV⁻²
    double gvp_ = gvp(eps_uV, eps_dV, eps_sV);
    double gvpp_ = gvpp(eps_uV, eps_dV, eps_sV, r2Ep, r2En, r2EsN, kp, kn, ksN, mN);
    double gv2p_ = gv2p(eps_uV, eps_dV, eps_sV, kp, kn, ksN);
    double gvn_ = gvn(eps_uV, eps_dV, eps_sV);
    double gvnp_ = gvnp(eps_uV, eps_dV, eps_sV, r2Ep, r2En, r2EsN, kp, kn, ksN, mN);
    double gv2n_ = gv2n(eps_uV, eps_dV, eps_sV, kp, kn, ksN);
    
    double Qw_BSM = Z * gvp_ + N * gvn_;  //in GeV⁻²
    return sqrt(0.5)*Qw_BSM*pow(10,-6)/GF;  //factor is due to sigma normalization, see eq. (92) in 2007.08529 and Pilar`s paper
}
*/





