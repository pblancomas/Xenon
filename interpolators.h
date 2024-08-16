#ifndef __INTERPOLATORS_H
#define __INTERPOLATORS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int findClosest(double X[], int size, double x0);
double quadraticInterpolation(double X[], double Y[], int size, double x0);
double linearInterpolation(double X[], double Y[], int size, double x0);


#endif
