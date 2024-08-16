#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "8B_flux.h"
#include "interpolators.h"
#include "FFM.h"



#define elec 0
#define muon 1
#define tau 2
#define NA 6.022*1e23 //Avogadro's number
#define GF 1.16637*pow(10,-23) //eV^-2


double eps_D(double (*eps_input)[3], double theta_13, double theta_12, double theta_23){
//eq. (17) in 2302.12846v2 

    double s13 = sin(theta_13);
    double c13 = cos(theta_13);
    double s12 = sin(theta_12);
    double c12 = cos(theta_12);
    double s23 = sin(theta_23);
    double c23 = cos(theta_23);

    double p1, p2, p3, p4, epsilon_D;
    
    p1 = c13*s13*(s23*eps_input[elec][muon] + c23*eps_input[elec][tau]);
    p2 = (1 + s13*s13)*c23*s23*eps_input[muon][tau];
    p3 = 0.5*c13*c13*(eps_input[elec][elec] - eps_input[muon][muon]);
    p4 = 0.5*(s23*s23 - s13*s13*c23*c23)*(eps_input[tau][tau] - eps_input[muon][muon]);
    
    epsilon_D = p1 - p2 - p3 + p4;
    
return epsilon_D;
}


double eps_N(double (*eps_input)[3], double theta_13, double theta_12, double theta_23){
//eq. (18) in 2302.12846v2 

    double s13 = sin(theta_13);
    double c13 = cos(theta_13);
    double s12 = sin(theta_12);
    double c12 = cos(theta_12);
    double s23 = sin(theta_23);
    double c23 = cos(theta_23);

    double p1, p2, epsilon_N;
    
    p1 = c13*(c23*eps_input[elec][muon] - s23*eps_input[elec][tau]);
    p2 = s13*(s23*s23*eps_input[muon][tau] - c23*c23*eps_input[muon][tau] + c23*s23*(eps_input[tau][tau] - eps_input[muon][muon]));
    
    epsilon_N = p1 + p2;
    
return epsilon_N;
}


double sol_density(double x0, int type){
//returns electron (type = 0) or neutron (type = 1) number density (in eV^3) inside the Sun, at a fraccional solar radius x0. Data from astro-ph/0511337, Table 11, GS98 solar model

     double X[28] = {
        0.000, 0.025, 0.050, 0.075, 0.100,
        0.150, 0.200, 0.250, 0.300, 0.350,
        0.400, 0.450, 0.500, 0.550, 0.600,
        0.650, 0.700, 0.750, 0.800, 0.850,
        0.900, 0.940, 0.950, 0.960, 0.970,
        0.980, 0.990, 1.000
    };
    
    double y0;
    
    if (type == 0){ //electron
    
    double Y[28] = {
        2.0125, 1.9981, 1.9581, 1.8998, 1.8295,
        1.6648, 1.4711, 1.2509, 1.0129, 0.7678,
        0.5244, 0.2877, 0.0602, -0.1574, -0.3649,
        -0.5621, -0.7460, -0.9130, -1.1024, -1.3341,
        -1.6462, -2.0306, -2.1685, -2.3399, -2.5683,
        -2.9124, -3.5662, -6.8436
    };  //log10(ne/NA)
    
    y0 = quadraticInterpolation(X, Y, 28, x0);
    
    return NA*pow(10,y0)*(197.327*pow(10.0,-7.0))*(197.327*pow(10.0,-7.0))*(197.327*pow(10.0,-7.0)); //convert from cm^-3 to eV^3
    }
    else if (type == 1){ //neutron
    
    double Y[28] = {
        1.6990, 1.6665, 1.5770, 1.4495, 1.3037,
        1.0102, 0.7440, 0.4950, 0.2462, -0.0050,
        -0.2560, -0.4911, -0.7209, -0.9407, -1.1501,
        -1.3489, -1.5582, -1.7395, -1.9289, -2.1605,
        -2.4726, -2.8571, -2.9950, -3.1663, -3.3947,
        -3.7388, -4.3927, -7.6700
    };  //log10(nn/NA)
    
    y0 = quadraticInterpolation(X, Y, 28, x0);
    
    return NA*pow(10,y0)*(197.327*pow(10.0,-7.0))*(197.327*pow(10.0,-7.0))*(197.327*pow(10.0,-7.0)); //convert from cm^-3 to eV^3
    }
    else{
        printf("ERROR! SELECT A DENSITY TYPE PLEASE.\n");
        return 1.0;
    }
}

double p(double x0, double (*eps_e)[3], double (*eps_p)[3], double (*eps_n)[3], double Enu, double theta_13, double theta_12, double theta_23, double Dm2_21){
//eq. (21) in 2302.12846v2 . Dm2_21 in ev^2, enu in MeV

    double epse_N, epsp_N, epsn_N, Ne, Nn, Acc;
    
    epse_N = eps_N(eps_e, theta_13, theta_12, theta_23);
    epsp_N = eps_N(eps_p, theta_13, theta_12, theta_23);
    epsn_N = eps_N(eps_n, theta_13, theta_12, theta_23);
    Ne = sol_density(x0, 0);
    Nn = sol_density(x0, 1);
    //Enu = 10*1e6;                                                //~10 MeV in eV, see for example https://arxiv.org/pdf/nucl-th/9601044 Fig.5 
    Acc = 2.0*sqrt(2.0)*GF*Ne*Enu*1e6;
    
    return sin(2*theta_12) + 2*(epse_N + epsp_N + (Nn/Ne)*epsn_N)*Acc/Dm2_21;
}

double q(double x0, double (*eps_e)[3], double (*eps_p)[3], double (*eps_n)[3], double Enu, double theta_13, double theta_12, double theta_23, double Dm2_21){
//eq. (21) in 2302.12846v2 . Dm2_21 in ev^2, Enu in MeV

    double epse_D, epsp_D, epsn_D, Ne, Nn, Acc;
    
    epse_D = eps_D(eps_e, theta_13, theta_12, theta_23);
    epsp_D = eps_D(eps_p, theta_13, theta_12, theta_23);
    epsn_D = eps_D(eps_n, theta_13, theta_12, theta_23);
    Ne = sol_density(x0, 0);
    Nn = sol_density(x0, 1);
    //Enu = 10*1e6;                                                //~10 MeV in eV
    Acc = 2*sqrt(2)*GF*Ne*Enu*1e6;
    
    return cos(2*theta_12) + (2*(epse_D + epsp_D + (Nn/Ne)*epsn_D) - cos(theta_13)*cos(theta_13))*(Acc/Dm2_21);
}

double f_p_8B(double x0){
//8B spatial distribution function for the production of 8B neutrinos. From http://etheses.dur.ac.uk/14853/1/main.pdf?DDD25+
//it goes to 0 after x = 0.1605, and is defined for x (fractional radius) between 0 and 1
    
double X[93] = {
    0.0000, 0.0025, 0.0050, 0.0068, 0.0082, 0.0094,
    0.0104, 0.0118, 0.0126, 0.0134, 0.0142,
    0.0150, 0.0158, 0.0166, 0.0174, 0.0182,
    0.0190, 0.0198, 0.0206, 0.0214, 0.0222,
    0.0230, 0.0238, 0.0246, 0.0254, 0.0262,
    0.0270, 0.0278, 0.0283, 0.0291, 0.0299,
    0.0309, 0.0321, 0.0330, 0.0345, 0.0361,
    0.0377, 0.0400, 0.0436, 0.0476, 0.0504,
    0.0522, 0.0536, 0.0550, 0.0564, 0.0576,
    0.0588, 0.0599, 0.0611, 0.0623, 0.0637,
    0.0647, 0.0659, 0.0671, 0.0683, 0.0695,
    0.0707, 0.0719, 0.0731, 0.0743, 0.0755,
    0.0767, 0.0776, 0.0788, 0.0797, 0.0812,
    0.0828, 0.0843, 0.0859, 0.0872, 0.0890,
    0.0914, 0.0936, 0.0960, 0.0986, 0.1016,
    0.1053, 0.1097, 0.1141, 0.1185, 0.1229,
    0.1272, 0.1318, 0.1360, 0.1404, 0.1448,
    0.1488, 0.1535, 0.1577, 0.1605, 0.1606, 0.1607, 1.0
};

double Y[93] = {
    0.0000, 0.2518, 0.6276, 1.1019, 1.5619, 2.1044,
    2.6274, 3.1203, 3.6408, 4.0954, 4.5845,
    5.0862, 5.5753, 6.0695, 6.5630, 7.0772,
    7.5977, 8.1151, 8.6293, 9.1404, 9.6608,
    10.1782, 10.6893, 11.2035, 11.7334, 12.2183,
    12.6866, 13.1506, 13.5896, 14.0223, 14.4456,
    14.9190, 15.4677, 15.9459, 16.4988, 17.0159,
    17.5572, 18.0361, 18.3159, 18.1371, 17.7115,
    17.3051, 16.8567, 16.4209, 15.9412, 15.4928,
    15.0162, 14.5553, 14.1101, 13.6398, 13.0910,
    12.6170, 12.1316, 11.6174, 11.0875, 10.5921,
    10.0967, 9.5919, 9.0839, 8.6011, 8.1527,
    7.7068, 7.3869, 7.0975, 6.5724, 6.0843,
    5.6116, 5.2265, 4.6817, 4.3028, 3.8728,
    3.3874, 2.9384, 2.5089, 2.0990, 1.7023,
    1.2954, 0.9294, 0.6524, 0.4523, 0.3206,
    0.1974, 0.1177, 0.0880, 0.0782, 0.0658,
    0.0679, 0.0538, 0.0557, 0.0632, 0.0000, 0.0000, 0.0000
};

    double y0 = linearInterpolation(X, Y, 93, x0);
    
    return y0;
}

int compute_density_matrix(double (*output_matrix)[3], double (*eps_e)[3], double (*eps_u)[3], double (*eps_d)[3], double Enu, double theta_13, double theta_12, double theta_23, double Dm2_21){
//eqs. (35)-(40) in 2302.12846v2 . Delta m^2 in eV^2, Enu in MeV, writes the result in output_matrix.
    
    /////////////////////////////////////////
    // build epsilon matrices for nucleons //
    /////////////////////////////////////////

    double eps_p[3][3], eps_n[3][3];

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
        
        eps_p[i][j] = 2*eps_u[i][j] + eps_d[i][j];
        eps_n[i][j] = eps_u[i][j] + 2*eps_d[i][j];
        
        }
    }    

    ////////////////////////////////////////
    // compute averaged cos(2*theta_12^m) //
    ////////////////////////////////////////
    
    double cos2thm, Delta, Pee;
    
    //equation (43), integration
    
    double N_steps = 50.0;  //for 0.005 precision
    cos2thm = 0.0;
    
    double x, fp_x, cos2thm_x, pval, qval;
    
    for (double i = 0.0; i < N_steps; i += 1.0){
    
        x = i/(N_steps - 1.0); 
        pval = p(x, eps_e, eps_p, eps_n, Enu, theta_13, theta_12, theta_23, Dm2_21);
        qval = q(x, eps_e, eps_p, eps_n, Enu, theta_13, theta_12, theta_23, Dm2_21);
        fp_x = f_p_8B(x);
        
        cos2thm_x = qval/sqrt(pval*pval + qval*qval);
        
        cos2thm += fp_x*cos2thm_x*(1.0/N_steps);
    }
    
    ////////////////////////////////
    // compute rest of parameters //
    ////////////////////////////////
    
    double s13 = sin(theta_13);
    double c13 = cos(theta_13);
    double s12 = sin(theta_12);
    double c12 = cos(theta_12);
    double s23 = sin(theta_23);
    double c23 = cos(theta_23);
    
    Pee = 0.5*(1 + cos(2*theta_12)*cos2thm);
    Delta = 0.5*s13*sin(2*theta_12)*sin(2*theta_23)*cos2thm;
    
    /////////////////
    // fill matrix //
    /////////////////
    
    output_matrix[elec][elec] = pow(s13,4) + pow(c13,4)*Pee;
    output_matrix[muon][muon] = c13*c13*(c23*c23*(1.0 - Pee) + s13*s13*s23*s23*(1.0 + Pee) + Delta);
    output_matrix[tau][tau] = c13*c13*(s23*s23*(1.0 - Pee) + s13*s13*c23*c23*(1.0 + Pee) - Delta);
    output_matrix[elec][muon] = c13*s13*s13*s13*s23 - 0.5*pow(c13,3)*(2.0*s13*s23*Pee + c23*sin(2.0*theta_12)*cos2thm);
    output_matrix[elec][tau] = c13*s13*s13*s13*c23 - 0.5*pow(c13,3)*(2.0*s13*c23*Pee - s23*sin(2*theta_12)*cos2thm);
    output_matrix[muon][tau] = 0.5*pow(c13,2)*(sin(2*theta_23)*((1.0 + s13*s13)*Pee - c13*c13) + 2.0*(1.0/tan(2*theta_23))*Delta);
    output_matrix[muon][elec] = output_matrix[elec][muon];
    output_matrix[tau][elec] = output_matrix[elec][tau];
    output_matrix[tau][muon] = output_matrix[muon][tau];
    
return 0;
}


int compute_dsigmadT(double (*out_matrix)[3], double T, double Enu, int Z, double mN, double (*eps_uV)[3], double (*eps_dV)[3], double (*eps_sV)[3]){
//T: recoil energy (in GeV). All dimensional inputs in GeV related units, output in GeV^-3
    
    ////////////////////
    // compute Q^2F^2 //
    ////////////////////
    
    int i, j, k;
    
    double r2Ep = 0.7071*(1000/197.327)*(1000/197.327);  //from 2007.08529v2
    double r2En = -0.1161*(1000/197.327)*(1000/197.327);
    double r2EsN = -0.0048*(1000/197.327)*(1000/197.327);
    double kp = 1.79284734462;
    double kn = -1.91304273;
    double ksN = -0.017;

    double QF1[3][3], QF2[3][3], result[3][3];
    
    QwFw_BSM(QF1, T, Z, mN, eps_uV, eps_dV, eps_sV, r2Ep, r2En, r2EsN, kp, kn, ksN);
    QwFw_BSM(QF2, T, Z, mN, eps_uV, eps_dV, eps_sV, r2Ep, r2En, r2EsN, kp, kn, ksN);
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            result[i][j] = 0.0;
        }
    }

    // Multiply QF1 and QF2 and store the result in result
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                result[i][j] += QF1[i][k] * QF2[k][j];
            }
        }
    }
    
    ////////////////////////////////////////////////////
    // give differential cross section in matrix form //
    ////////////////////////////////////////////////////
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            out_matrix[i][j] = mN*(1 - (mN*T/(2*Enu*Enu)) - (T/Enu))*result[i][j]/(2*M_PI);
        }
    }

    return 0;
}

double dNdT(double T, double norm, double flux_err, double flux_pull, double theta_13, double theta_12, double theta_23, double Dm2_21, double (*eps_e)[3], double (*eps_u)[3], double (*eps_d)[3]){
// differential rate of events with recoil energy T.
// Input: T, in keV, normalizacion, number of targets*time in seconds, flux porcentual error and pull, NSI matrices for e, u, d. 
// Output: dN/dT in #events*keV^-1

    double T_GeV = T*0.000001;
    double eps_s[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}; 
    double Enu;  
    double result = 0.0;
    double sigma_matrix[3][3], rho_matrix[3][3], prod_matrix[3][3];
    
    // experiment data //
    
    int Z_Xe = 54;
    double m_Xe = 131.293*0.9314941 ;  //GeV

    
    //////////////////////////////////
    // integrate in neutrino energy //
    //////////////////////////////////
    
    int i, j, k;
    
    double Enu_min = 0.5*(T_GeV + sqrt(T_GeV*T_GeV + 2*m_Xe*T_GeV)); //minimum neutrino energy in GeV
    double Enu_max = 16.0/1000.0; //maximum neutrino energy in GeV
    double step = (Enu_max - Enu_min)/300.0;  //300 steps for 0.001 precision
    double flux, conv;
    conv = (197.327*pow(10.0,-13.0))*(197.327*pow(10.0,-13.0))*0.001;//MeV^-1 cm^-2 s^-1 to GeV/s
    
    for (Enu = Enu_min; Enu < Enu_max; Enu += step){
        
        // Initialize matrices to 0
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                prod_matrix[i][j] = 0.0;
                sigma_matrix[i][j] = 0.0;
                rho_matrix[i][j] = 0.0;
            }
        }
        
        flux = (1.0 + flux_err*flux_pull)*dphidEnu(Enu*1000.0)*conv;  //the input for our flux is in MeV, the output is in MeV^-1 cm^-2 s^-1, we want it in GeV/s
        
        compute_density_matrix(rho_matrix, eps_e, eps_u, eps_d, Enu*1000.0, theta_13, theta_12, theta_23, Dm2_21);
        
        compute_dsigmadT(sigma_matrix, T_GeV, Enu, Z_Xe, m_Xe, eps_u, eps_d, eps_s);
    
        // Multiply and store the result in result
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    prod_matrix[i][j] += rho_matrix[i][k] * sigma_matrix[k][j];
                }
            }
        }
    
        // traze of the matrix //
        
        double traze = 0.0;
    
        for (int i = 0; i < 3; i++) traze += prod_matrix[i][i];
        result += norm*traze*flux*step;  // s*(GeV^-3)*(GeV/s)*GeV = GeV^-1
    }

    return result*0.000001;  //output in keV^-1
}
