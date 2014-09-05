//Author: Shelby Straight, Paesani Lab
//Date Last Edited: Sep 3, 2014
//
//Gaussian objected compute random numbers between 0 and 1
//This algorithm is an adaptation of code originally written
//in fortran 90 by F. Paesani, which was originally adapted
//from Numerical Recipes

#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include "gaussian.h"
#include <time.h>

long double gaussian::compute_gaussian(long double vel){
    
    srand(time(NULL)+vel);

    bool is_set=false;

    long double v1, v2, rsq;

    while (is_set==false){

	long double random_01_1=((long double)rand())/ ((long double)RAND_MAX);
	long double random_01_2=((long double)rand())/ ((long double)RAND_MAX);

	v1=(2.0*random_01_1) - 1.0;
	v2=(2.0*random_01_2) - 1.0;

	rsq = ((v1*v1)+(v2*v2));

	if (rsq>0.0 && rsq<1.0) is_set=true;

    }

    long double factor= std::sqrt((-2.0*std::log(rsq))/rsq);

    long double gset=factor*v1;
    long double gauss2=factor*v2;

    return gauss2;
}
