#include "transition_matrix.h"
#include <cmath>
#include <iomanip>
#include <iostream>

const long double round_off=1.0e-15;

void transition_matrix::deallocate(){
 

    delete nm_to_cart;
    delete cart_to_nm;
    delete diagonals;
    delete eigenvalues;
    delete freqs_root;
    delete freqs;
    delete fictitious_masses;

}

void transition_matrix::convert_coordinates(long double spring_omega, int num_beads, long double phys_mass, long double gamma){

	number_of_beads=num_beads;

	int matrix_size=number_of_beads+1;

	cart_to_nm=new long double*[matrix_size];
	nm_to_cart=new long double*[matrix_size];
	diagonals=new long double*[matrix_size];

	eigenvalues=new long double[matrix_size];
	fictitious_masses=new long double[matrix_size];
	freqs= new long double[matrix_size];
	freqs_root=new long double[matrix_size];

	physical_mass=phys_mass;
	omegaP=spring_omega;

	for (int i=1; i<=number_of_beads; ++i){

		cart_to_nm[i]=new long double[matrix_size];
		nm_to_cart[i]=new long double[matrix_size];
		diagonals[i]=new long double[matrix_size];

	}

	this->normal_to_cart();

	this->cart_to_normal();

	this->check_unitary();

	this->compute_eigen_freqs_fictm(gamma);
}

long double transition_matrix::compute_metric(){

    long double sum=0.0;

    for (int i=2; i<=number_of_beads; ++i){
	sum+=(eigenvalues[i]*cart_to_nm[i][i]);
    }

    return sum;
}


void transition_matrix::compute_eigen_freqs_fictm(long double gamma){

	eigenvalues[1]=0.0;

	for (int i=1; i<=((number_of_beads/2)-1); ++i){

		long double arg= 2* M_PI*((long double) i)/((long double) number_of_beads);
		long double eigenvalue= 2.0* (1.0-std::cos(arg));
		
		eigenvalues[(2*i)]=eigenvalue*number_of_beads;
		eigenvalues[(2*i)+1]=eigenvalue*number_of_beads;

	}

	if (number_of_beads>1) eigenvalues[number_of_beads]=4.0*number_of_beads;

	for (int i=1; i<=number_of_beads; ++i){
	    if (i>1) freqs[i]=omegaP/gamma;
	    else if (i==1) freqs[i]=omegaP; 
    	    freqs_root[i]=std::sqrt(freqs[i]);
	    fictitious_masses[i]=physical_mass*freqs[i];
	}

}

void transition_matrix::check_unitary(){

	for (int i=1; i<=number_of_beads; ++i){
		for (int j=1; j<=number_of_beads;++j){
			long double diag_val=0.0;
			for (int k=1; k<=number_of_beads; ++k){
				diag_val+=(nm_to_cart[i][k]*cart_to_nm[k][j]);
			}
			diagonals[i][j]=diag_val;
			if ((diagonals[i][j]*diagonals[i][j])<(round_off*round_off)) diagonals[i][j]=0.0;
		}
	}
}

void transition_matrix::cart_to_normal(){

	for (int i=1; i<=number_of_beads; ++i){

		for (int k=1; k<=number_of_beads; ++k){

			cart_to_nm[i][k]=(nm_to_cart[k][i]/((long double)number_of_beads));
		}

	}
}

void transition_matrix::normal_to_cart(){

	for (int i=1; i<=number_of_beads; ++i){
		nm_to_cart[i][1]=1.0;
	}

	for (int i=1; i<=(number_of_beads/2); ++i){
		nm_to_cart[(2*i)][number_of_beads]=-1.0;
		nm_to_cart[(2*i)-1][number_of_beads]=1.0;
	}

	for (int i=1; i<=number_of_beads; ++i){
		
		long double arg= 2.0 * ((long double)(i-1)) * M_PI / ((long double)number_of_beads);
		
		for (int j=1; j<=((number_of_beads/2)-1); ++j){

			nm_to_cart[i][(2*j)]= std::sqrt(2.0)*std::cos(arg*((long double)j));
//			if ((nm_to_cart[i][(2*j)])*(nm_to_cart[i][(2*j)])<=(round_off*round_off)) nm_to_cart[i][(2*j)]=0;

			nm_to_cart[i][(2*j)+1]= -1.0*std::sqrt(2.0)*std::sin(arg*((long double)j));
//			if ((nm_to_cart[i][(2*j)+1]*nm_to_cart[i][(2*j)+1])<=(round_off*round_off)) nm_to_cart[i][(2*j)+1]=0;
		}

	}
}

void transition_matrix::print_1D_matrix(long double * matrix){

	std::cout<<'\n';

	for (int i=1; i<=number_of_beads; ++i){

		std::cout<<std::setw(18)<<matrix[i];

	}

	std::cout<<'\n';
}

void transition_matrix::print_2D_matrix(long double ** matrix){

	std::cout<<'\n';

	for (int i=1; i<=number_of_beads; ++i){

		for (int j=1; j<=number_of_beads; ++j){

			std::cout<<std::setw(18)<<matrix[i][j];
		}

		std::cout<<'\n';
	}
}

