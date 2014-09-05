//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 4, Sep, 2014
//
//Simulation set-up is the "base" of the simulation-
//it is responsible for reading input (if continuing from a previous
//simulation) and for printing output, both to the screen and to the
//files which would dictate the starting configuration of the next
//simulation

#ifndef SIMULATION_SETUP_H
#define SIMULATION_SETUP_H

#include "necklace.h"
#include "transition_matrix.h"
#include "virial_estimator.h"
#include "barker_estimator.h"

class simulation_setup{

public:
    void read_starting_configuration(long double gamma, int detach_cthermo, int rstrtPIMD,
				     necklace normal_mode_ring, necklace cartesian_ring,
				     transition_matrix tm, int potential_choice);

    void print_captured_trajectory(long double time_stamp, necklace normal_mode_ring,
				   virial_estimator virial, barker_estimator barker,
				   long double conserved, long double analytical_solution);

    void print_final_configuration(necklace normal_mode_ring);

    void propagate_frame(necklace normal_mode_ring, necklace cartesian_ring,
			 transition_matrix tm, long double dt);

private:
    long double gamma_;
    int potential_choice_;
    int detach_centroid_thermostat_;
    int restart_PIMD_;
};
#endif
 
