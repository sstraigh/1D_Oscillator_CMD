//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 4, Sep, 2014
//
//Particle_tracker is a class which records and updates the instantaneous and
//cumulative averages of each particle's kinetic and potential energy-for the
//sake of the tracker, the inter-bead interactions on the necklace are ignored.
//In addition, the particle's instantaneous position and velocity is also recorded
//
//This is useful to check simulation convergence by monitoring the centroid
//variable's average kinetic and potential energies
//
//This class also includes a print_xyz function, which enables simulation 
//visualization in a package such as VMD

#ifndef PARTICLE_TRACKER_H
#define PARTICLE_TRACKER_H

#include "necklace.h"

class particle_tracker{

public:
    void initialize(int number_of_particles);
    void update_tracker(necklace my_necklace, int weight_of_average);

    void deallocate();

    void print_xyz(int num_particles);

    long double get_avg_KE(int bead_index){return avg_kinetic_energies_[bead_index];}
    long double get_avg_PE(int bead_index){return avg_potential_energies_[bead_index];} 
    long double get_instant_KE(int bead_index){return instant_KE_[bead_index];}
    long double get_instant_PE(int bead_index){return instant_PE_[bead_index];}

private:
    long double * avg_kinetic_energies_;
    long double * avg_potential_energies_;
    long double * instant_PE_;
    long double * instant_KE_;

    long double * instant_positions_;
    long double * instant_velocities_;

};

#endif
