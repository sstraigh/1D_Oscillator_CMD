#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

class transition_matrix{

public:

	long double ** nm_to_cart;
	long double ** cart_to_nm;
	long double ** diagonals;
	long double * eigenvalues;
	long double * fictitious_masses;

	long double * freqs;
	long double * freqs_root;


	long double omegaP;
	long double physical_mass;
	long double kT;
	long double Beta;
	long double hbar;
	
	int number_of_beads;

	void convert_coordinates(long double spring_omega, int num_beads, long double phys_mass, long double gamma);

	void normal_to_cart();
	void cart_to_normal();
	void check_unitary();
	void compute_eigen_freqs_fictm(long double gamma);

	void print_2D_matrix(long double ** matrix);
	void print_1D_matrix(long double * matrix);

	void deallocate();                                                         

	long double compute_metric();

};
#endif
