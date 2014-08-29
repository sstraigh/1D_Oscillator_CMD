#ifndef NECKLACE_H
#define NECKLACE_H

#include "particle.h"

class necklace{

public:

    particle * particle_array;

    int potential_choice;

    int number_of_beads;
    long double bead_freq;

    long double m_prime;

    long double spr_const;
    long double nhc_contribution;

    long double kinetic_E;
    long double potential_E;

    long double classical_potential;
    long double bead_potential; 

    long double physical_mass;

    necklace();

    long double compute_potential_energy();
    long double compute_total_hamiltonian();
    long double compute_NHC_hamiltonian();

    long double compute_bead_KE(int bead_index);
    long double compute_bead_PE(int bead_index);

    long double compute_external_force(int bead_index);
    long double compute_external_potential(int bead_index);

    long double compute_bead_force(int bead_index);

    void normal_modes(long double hbar, int num_beads, long double * mass, long double * initial_location, long double * initial_velocity, long double * spring_const, int ndof, int nh_chainlength, long double timestep, long double kT, bool activate_CMD, int pot_choice);
    void create_necklace(long double hbar, int num_beads, long double mass, long double initial_location, long double initial_velocity, long double spring_const, int ndof, int nh_chainlength, long double timestep, long double kT, int pot_choice);

    void deallocate();

};
#endif
