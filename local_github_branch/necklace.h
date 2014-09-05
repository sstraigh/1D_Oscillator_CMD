//Author: Shelby Straight, UCSD, Paesani Lab summer research rotation
//Date last edited: Sep 4444014
//
//Necklace.h is the header file which defines the data structure
//responsible for the propagation of particle dynamics.
//A single particle's quantum effects are captured by modelling it
//as a collection of particles, each of which are connected to its
//neighbors by a harmonic spring. The spring frequency is dictated
//by the parameters of the simulation

#ifndef NECKLACE_H
#define NECKLACE_H

#include "particle.h"

class necklace{

 public:
  particle * particle_array;

  long double compute_total_hamiltonian();
  long double compute_NHC_hamiltonian();

  long double get_necklace_KE(){return kinetic_E;}
  long double get_necklace_PE(){return potential_E;}

  long double get_necklace_bead_potential(){return bead_potential;}
  long double get_necklace_classical_potential(){return classical_potential;}

  long double compute_bead_force(int bead_index);

  void normal_modes(long double hbar, int num_beads, long double * mass, 
		    long double * initial_location, 
		    long double * initial_velocity, long double * spring_const, 
		    int ndof, int nh_chainlength, long double timestep, 
		    long double kT, bool activate_CMD, int pot_choice);

  void create_necklace(long double hbar, int num_beads, long double mass, 
		       long double initial_location, 
		       long double initial_velocity, long double spring_const, 
		       int ndof, int nh_chainlength, long double timestep, 
		       long double kT, int pot_choice);

  void deallocate();

  int get_num_beads(){return number_of_beads;}
  long double compute_external_potential(int bead_index); 
  long double compute_bead_KE(int bead_index);

 private:
  long double physical_mass;
  long double spr_const;
  long double m_prime;
  int potential_choice;
  int number_of_beads;
  long double bead_freq;
  long double nhc_contribution;

  long double kinetic_E;
  long double potential_E;

  long double classical_potential;
  long double bead_potential;
 
  long double compute_bead_PE(int bead_index);
  long double compute_external_force(int bead_index);
};
#endif
