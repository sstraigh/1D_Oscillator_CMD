#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>

#include "thermostat.h"
#include "particle.h"
#include "necklace.h"
#include "transition_matrix.h"
#include "virial_estimator.h"
#include "barker_estimator.h"
#include "simulation_setup.h"
#include "readin.h"
#include "particle_tracker.h"

//Morse potential expression
const long double alpha=1.1563;
const long double D=0.18748;
const long double r0=1.8325;
const long double frac=7.0/12.0;

int main(int argc, char *argv[]){

    long double plancks_h, hbar, spring_constant, particle_mass, v0, r_init, dt,
	 kT, Beta, gamma;

    int spatial_deg_of_freedom, nhc_length, number_of_steps, number_of_particles,
       	capture_rate, potential_choice;
    
    bool detach_cthermo=false;
    int c_thermo=-1;
    if (c_thermo==1) detach_cthermo=true;

    int rstrtPIMD=-1;

    readin reading_frame;
    reading_frame.read_parameters(number_of_steps, number_of_particles, nhc_length, kT, capture_rate,
				 gamma, potential_choice, c_thermo, rstrtPIMD, hbar, spring_constant,
				 particle_mass, v0, r_init, dt, spatial_deg_of_freedom);

    Beta=1.0/kT;

    long double * fictitious_masses=new long double[number_of_particles+1];
    long double * normal_mode_coord=new long double[number_of_particles+1]; 
    long double * nm_velocities=new long double[number_of_particles+1];

    //neck_beads is the data structure which holds the particles in its cartesian coordinates
    //it is responsible for computing the forces and potentials (both from the external-ie central-
    //spring and the inter-bead springs)
    
    necklace neck_beads;
    neck_beads.create_necklace(hbar, number_of_particles, particle_mass, r_init, v0, spring_constant, spatial_deg_of_freedom, nhc_length, dt, kT, potential_choice);

    long double root_beads=std::sqrt(number_of_particles);
    transition_matrix tm;
    tm.convert_coordinates((root_beads*kT), number_of_particles, particle_mass, gamma);

    //    fictitious masses are used for the propagation of dynamics in the PIMD/CMD representation
    fictitious_masses[1]=particle_mass;
    for (int i=2; i<=number_of_particles; ++i) fictitious_masses[i]=particle_mass*tm.get_eigenvalue(i)*gamma*gamma;

    //this routine uses the tm conversion matrices to initialize the NM coordinates and velocities
    for (int i=1; i<=number_of_particles; ++i){
	nm_velocities[i]=0.0L;
	normal_mode_coord[i]=0.0L;
	for (int k=1; k<=number_of_particles; ++k){
	    nm_velocities[i]+=tm.get_cart_to_nm(i,k)*neck_beads.particle_array[k-1].get_velocity();
	    normal_mode_coord[i]+=tm.get_cart_to_nm(i,k)*neck_beads.particle_array[k-1].get_location();
	}
    }

    //normal_mode_ring is the data structure which holds the positions of the NMs (important because of CMD)
    //it is also responsible for the computation of the kinetic energies associated with the fictitious Hamiltonian
    necklace normal_mode_ring;
    normal_mode_ring.normal_modes(hbar, number_of_particles, fictitious_masses, normal_mode_coord, nm_velocities, tm.freqs_, spatial_deg_of_freedom, nhc_length, dt, kT, detach_cthermo, potential_choice);

    simulation_setup simulation_base;
    simulation_base.read_starting_configuration(gamma, c_thermo, rstrtPIMD, normal_mode_ring, neck_beads, tm, potential_choice);

    //The analytical solution to the harmonic oscillator at temperature Beta=1/kT is calculated from the 
    //partition function of the harmonic oscillator
    long double analytical_solution=0.5*hbar*spring_constant;
    analytical_solution+=(hbar*spring_constant*std::exp(-1.0*Beta*hbar*spring_constant))/(1-std::exp(-1.0*Beta*hbar*spring_constant));

    virial_estimator virial_calculator;
    barker_estimator barker_calculator;

    if (potential_choice==1){
	virial_calculator.initialize_estimator(number_of_particles, spring_constant);
	barker_calculator.initialize_estimator(number_of_particles, spring_constant, particle_mass, kT, hbar);
    }

    particle_tracker nm_tracker;
    nm_tracker.initialize(number_of_particles);

    //This loop starts the propagation of the simulation
    for (int j=1; j<=number_of_steps; ++j){

	long double Total_E=0.0;
	long double Kin_E=0.0;
	long double conserved=0.0;
	long double NHam=0.0;

	simulation_base.propagate_frame(normal_mode_ring, neck_beads, tm, dt);

	//use the thermostat hamiltonian and kinetic energy contributions from the 
	//ring representation in the normal mode coordinates...
	NHam    = normal_mode_ring.compute_NHC_hamiltonian();
	Total_E = normal_mode_ring.compute_total_hamiltonian();
	Kin_E   = normal_mode_ring.get_necklace_KE();

	//...and the potential energy contributions (both internal necklace-spring
	//and external portions) to compute the conserved extended hamiltonian
	long double Total_E1 = neck_beads.compute_total_hamiltonian();
	long double bead_pot = neck_beads.get_necklace_bead_potential();
	long double ext_pot  = neck_beads.get_necklace_classical_potential();
	
	//calculation of the virial and barker energy estimators
	if (j%capture_rate==0){
	    nm_tracker.update_tracker(normal_mode_ring, (j/capture_rate));

	    if (potential_choice==1){
    		virial_calculator.update_energy(neck_beads, (j/capture_rate));
    		barker_calculator.update_energy(neck_beads, (j/capture_rate));
    	    }
	}
	//compute the extended hamiltonian regardless of NM vs. primitive
    	conserved=(Kin_E+bead_pot+ext_pot+NHam);

	simulation_base.print_captured_trajectory((j*dt), normal_mode_ring, virial_calculator, barker_calculator, conserved, analytical_solution);
    }

    simulation_base.print_final_configuration(normal_mode_ring);

    delete fictitious_masses;
    delete normal_mode_coord;
    delete nm_velocities;

    tm.deallocate();
    neck_beads.deallocate();
    normal_mode_ring.deallocate();

    nm_tracker.deallocate();

    return 0;
}
