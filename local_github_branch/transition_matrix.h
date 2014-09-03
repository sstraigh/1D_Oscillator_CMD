//Author: Shelby Straight, Paesani Lab
//Date last edited: 3, Sep, 2014
//
//tranistion_matrix.h defines the transition_matrix class
//This is a class responsible for coordinate transformations
//between the normal and cartesian representations.
//
//The algorithm which computes the matrices responsible for
//converting coordinates was adapted from Fortran code written
//by F. Paesani. For further reference on the algorithm see:
//     ref 1.
//     ref 2.
//     ref 3.


#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

class transition_matrix{

 public:
  long double get_nm_to_cart(int row, int col){return nm_to_cart_[row][col];}
  long double get_cart_to_nm(int row, int col){return cart_to_nm_[row][col];}
  long double get_frequency(int index){return freqs_[index];}
  long double get_eigenvalue(int index){return eigenvalues_[index];}

  //creates the conversion matrices based on simulation parameters
  void convert_coordinates(long double spring_omega, int num_beads, long double phys_mass, long double gamma);

  //frees memory associated with dynamically allocated arrays
  void deallocate();

  long double * freqs_; //must be passed into another data structure, therefore it is publice (for now??)

 private:
  long double ** nm_to_cart_;
  long double ** cart_to_nm_;
  long double ** diagonals_; //is this used outside of check_unitary? 
                             //no, but it is important to maintain the 
			     //ability to check that the code is working
  long double * eigenvalues_;

  long double omegaP_;
	
  int number_of_beads_;

  //these four functions are called by convert_coordinates to create the matrices
  void normal_to_cart();
  void cart_to_normal();
  void check_unitary();
  void compute_eigen_freqs(long double gamma);

  //used to double check for correct initialization
  long double compute_metric();

  //prints matrices associated with the class in accordance with 
  //indexes starting at 1
  void print_2D_matrix(long double ** matrix);
  void print_1D_matrix(long double * matrix);

};
#endif
