//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 3, Sep, 2014
//
//implementation of transition_matrix class. 

#include "transition_matrix.h"
#include <cmath>
#include <iomanip>
#include <iostream>

const long double round_off=1.0e-15;

void transition_matrix::deallocate(){
 
    delete nm_to_cart_;
    delete cart_to_nm_;
    delete diagonals_;
    delete eigenvalues_;
    delete freqs_;

}

void transition_matrix::convert_coordinates(long double spring_omega, int num_beads, long double phys_mass, long double gamma){

	number_of_beads_=num_beads;

	int matrix_size=number_of_beads_+1;

	cart_to_nm_=new long double*[matrix_size];
	nm_to_cart_=new long double*[matrix_size];
	diagonals_=new long double*[matrix_size];

	eigenvalues_=new long double[matrix_size];
	freqs_= new long double[matrix_size];

	omegaP_=spring_omega;

	for (int i=1; i<=number_of_beads_; ++i){

		cart_to_nm_[i]=new long double[matrix_size];
		nm_to_cart_[i]=new long double[matrix_size];
		diagonals_[i]=new long double[matrix_size];

	}

	this->normal_to_cart();

	this->cart_to_normal();

	this->check_unitary();

	this->compute_eigen_freqs(gamma);
}

//used for double checking the conversion matrices have been initialized correctly
long double transition_matrix::compute_metric(){

    long double sum=0.0;

    for (int i=2; i<=number_of_beads_; ++i){
	sum+=(eigenvalues_[i]*cart_to_nm_[i][i]);
    }

    return sum;
}


void transition_matrix::compute_eigen_freqs(long double gamma){

	eigenvalues_[1]=0.0;

	for (int i=1; i<=((number_of_beads_/2)-1); ++i){

		long double arg= 2* M_PI*((long double) i)/((long double) number_of_beads_);
		long double eigenvalue= 2.0* (1.0-std::cos(arg));
		
		eigenvalues_[(2*i)]=eigenvalue*number_of_beads_;
		eigenvalues_[(2*i)+1]=eigenvalue*number_of_beads_;

	}

	if (number_of_beads_>1) eigenvalues_[number_of_beads_]=4.0*number_of_beads_;

	for (int i=1; i<=number_of_beads_; ++i){
	    if (i>1) freqs_[i]=omegaP_/gamma;
	    else if (i==1) freqs_[i]=omegaP_; 
	}

}

//routine used to check if the matrix transformation is performed correctly
//to check, run check_unitary and then print the 2D diagonal matrix using the
//print 2d matrix routine in this class. The output of the diagonals_ matrix
//should be the identity matrix for the system under consideration.
void transition_matrix::check_unitary(){

	for (int i=1; i<=number_of_beads_; ++i){
		for (int j=1; j<=number_of_beads_;++j){
			long double diag_val=0.0;
			for (int k=1; k<=number_of_beads_; ++k){
				diag_val+=(nm_to_cart_[i][k]*cart_to_nm_[k][j]);
			}
			diagonals_[i][j]=diag_val;
			if ((diagonals_[i][j]*diagonals_[i][j])<(round_off*round_off)) diagonals_[i][j]=0.0;
		}
	}
}

//create the cartesian-to-normal coordinate conversion matrix 
void transition_matrix::cart_to_normal(){

	for (int i=1; i<=number_of_beads_; ++i){

		for (int k=1; k<=number_of_beads_; ++k){

			cart_to_nm_[i][k]=(nm_to_cart_[k][i]/((long double)number_of_beads_));
		}

	}
}

//create the normal-mode-to-cartesian coordinate conversion matrix
void transition_matrix::normal_to_cart(){

	for (int i=1; i<=number_of_beads_; ++i){
		nm_to_cart_[i][1]=1.0;
	}

	for (int i=1; i<=(number_of_beads_/2); ++i){
		nm_to_cart_[(2*i)][number_of_beads_]=-1.0;
		nm_to_cart_[(2*i)-1][number_of_beads_]=1.0;
	}

	for (int i=1; i<=number_of_beads_; ++i){
		
		long double arg= 2.0 * ((long double)(i-1)) * M_PI / ((long double)number_of_beads_);
		
		for (int j=1; j<=((number_of_beads_/2)-1); ++j){

			nm_to_cart_[i][(2*j)]= std::sqrt(2.0)*std::cos(arg*((long double)j));

			nm_to_cart_[i][(2*j)+1]= -1.0*std::sqrt(2.0)*std::sin(arg*((long double)j));
		}

	}
}

void transition_matrix::print_1D_matrix(long double * matrix){

	std::cout<<'\n';

	for (int i=1; i<=number_of_beads_; ++i){

		std::cout<<std::setw(18)<<matrix[i];

	}

	std::cout<<'\n';
}

void transition_matrix::print_2D_matrix(long double ** matrix){

	std::cout<<'\n';

	for (int i=1; i<=number_of_beads_; ++i){

		for (int j=1; j<=number_of_beads_; ++j){

			std::cout<<std::setw(18)<<matrix[i][j];
		}

		std::cout<<'\n';
	}
}

