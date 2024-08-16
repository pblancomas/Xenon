#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolators.h"

// Función para encontrar el índice del elemento más cercano en X
int findClosest(double X[], int size, double x0) {

    if (x0 < X[0] || x0 > X[size - 1]){
        printf("ERROR. x0 OUT OF BOUNDS\n");
        return -1;
    }

    int closestIndex = 0;
    double minDiff = fabs(X[0] - x0);
    for (int i = 1; i < size; i++) {
        double diff = fabs(X[i] - x0);
        if (diff < minDiff) {
            minDiff = diff;
            closestIndex = i;
        }
    }
    return closestIndex;
}

// Función para realizar la interpolación cuadrática
double quadraticInterpolation(double X[], double Y[], int size, double x0) {
    // Encuentra el índice del elemento más cercano a x0
    int idx = findClosest(X, size, x0);

    // Asegúrate de que tienes tres puntos para la interpolación
    int idx1 = (idx - 1 < 0) ? idx : idx - 1;
    int idx2 = idx;
    int idx3 = (idx + 1 >= size) ? idx : idx + 1;

    // Si idx está en el borde, ajusta los índices
    if (idx1 == idx) {
        idx2 = idx + 1;
        idx2 = idx + 2;
    }
    if (idx3 == idx) {
        idx1 = idx - 2;
        idx2 = idx - 1;
    }

    double x1 = X[idx1], x2 = X[idx2], x3 = X[idx3];
    double y1 = Y[idx1], y2 = Y[idx2], y3 = Y[idx3];

    // Calcula los coeficientes del polinomio de Lagrange
    double L1 = (x0 - x2)*(x0 - x3)/((x1 - x2)*(x1 - x3));
    double L2 = (x0 - x1)*(x0 - x3)/((x2 - x1)*(x2 - x3));
    double L3 = (x0 - x1)*(x0 - x2)/((x3 - x1)*(x3 - x2));

    // Evalúa el polinomio en x0
    double y0 = y1*L1 + y2*L2 + y3*L3;

    return y0;
}


// Función para realizar la interpolación lineal
double linearInterpolation(double X[], double Y[], int size, double x0) {
    
    if (x0 < X[0] || x0 > X[size - 1]){
        printf("ERROR. x0 OUT OF BOUNDS\n");
        return 1e20;
    }
    
    // Encuentra idx1, idx2 tal que x0 está entre ambos
    int idx1, idx2;
    for (int i = 0; i < size - 1; i++){
        if (X[i] <= x0 && x0 < X[i + 1]){
        idx1 = i;
        idx2 = i+1;
        break;
        }
    }
    
    //check if x0 is right edge
    if (x0 == X[size - 1]){
        idx1 = size - 2;
        idx2 = size - 1;
    }

    double x1 = X[idx1], x2 = X[idx2];
    double y1 = Y[idx1], y2 = Y[idx2];

    // usamos la fŕomula de interpolación lineal
    double y0 = y1 + (x0 - x1)*(y2 - y1)/(x2 - x1);
    
    return y0;
}
