#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "core.h"

#define NA 6.022*1e23 //Avogadro's number


//////////////////
////// MAIN //////
//////////////////


int main(int argc, char *argv[]) 
{
    double eps_mat[3][3], rho_matrix[3][3];
    int i, j;
    
        // Initialize the result matrix to 0
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                eps_mat[i][j] = 0.0;  //no NSIs for tests
                rho_matrix[i][j] = 0.0;
            }
        }
    
    
    // from nufit //

    double theta_13 = 8.58*M_PI/180;
    double theta_12 = 33.67*M_PI/180;
    double theta_23 = 42.3*M_PI/180;
    double Dm2_21 = 7.41*pow(10,-5); //in eV^2
    
    
    double norm = 1.0*1000000.0*(1.0/131.293)*NA*24.0*3600.0; // 1 tonne.day
    double flux_error = 0.16;  //from Bahcall
    double flux_pull = 0.0; //no to add errors in the flux need now
    
    double N = 0.0;
    
    printf("T (keV) : dNdT (ton*day*keV)^-1\n");
    
    for (double T = 0.3; T < 3.0; T += 0.01){
    
    N = dNdT(T, norm, flux_error, flux_pull, theta_13, theta_12, theta_23, Dm2_21, eps_mat, eps_mat, eps_mat);
    printf("%g : %g\n", T, N); 
    
    }  
    

exit(0);
}

