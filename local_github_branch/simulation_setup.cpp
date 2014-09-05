#include "simulation_setup.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

void simulation_setup::propagate_frame(necklace normal_mode_ring, necklace cartesian_ring,
                              transition_matrix tm, long double dt){

    long double * nm_forces=new long double[normal_mode_ring.get_num_beads()+1];

    //compute initial forces in NM representation
    for (int i=1; i<=normal_mode_ring.get_num_beads(); ++i){
	nm_forces[i]=0.0L;
	for (int k=1; k<=normal_mode_ring.get_num_beads(); ++k){
 	    nm_forces[i]+=tm.get_cart_to_nm(i,k)*cartesian_ring.compute_bead_force(k-1)*normal_mode_ring.get_num_beads();
	}
    }
	   
    //1st part of Velocity Verlet algorithm

    //(note: for CMD-NVE, there is no thermostat, so the detach_cthermo flag ensures that no discontinuities are
    //fed into the calculations)

    //update forces and accelerations in the normal mode representations, scale the 
    //velocities according to the thermostats
    for (int i=0; i<normal_mode_ring.get_num_beads(); ++i){
	long double acc=nm_forces[i+1]/normal_mode_ring.particle_array[i].get_mass();
	normal_mode_ring.particle_array[i].set_acceleration(acc);
	if (detach_centroid_thermostat_==0 || i!=0) normal_mode_ring.particle_array[i].integrate_thermostat(dt/2.0L);
    }

    //integrate velocities to half-time, positions to full-time
    for (int i=0; i<normal_mode_ring.get_num_beads(); ++i){
	long double vel=normal_mode_ring.particle_array[i].get_velocity();
	vel += (dt/2.0L)*normal_mode_ring.particle_array[i].get_acceleration();
	normal_mode_ring.particle_array[i].set_velocity(vel);

	long double loc=normal_mode_ring.particle_array[i].get_location();
	loc+=dt*normal_mode_ring.particle_array[i].get_velocity();
	normal_mode_ring.particle_array[i].set_location(loc);
    }

    //update cartesian coordinate in order to recompute forces
    for (int i=1; i<=normal_mode_ring.get_num_beads(); ++i){
	cartesian_ring.particle_array[i-1].set_location(0.0L);
	long double loc=0.0L;
	for (int k=1; k<=normal_mode_ring.get_num_beads(); ++k){
	    loc+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_location();
	}
	cartesian_ring.particle_array[i-1].set_location(loc);
    }	    

    //compute forces in NM representation
    for (int i=1; i<=normal_mode_ring.get_num_beads(); ++i){
	nm_forces[i]=0.0L;
	for (int k=1; k<=normal_mode_ring.get_num_beads(); ++k){
	    nm_forces[i]+=(tm.get_cart_to_nm(i,k)*cartesian_ring.compute_bead_force(k-1))*normal_mode_ring.get_num_beads();
	}
    }

    //second portion of VV
    for (int i=0; i<normal_mode_ring.get_num_beads(); ++i){
	normal_mode_ring.particle_array[i].set_acceleration(nm_forces[i+1]/normal_mode_ring.particle_array[i].get_mass());
	long double vel=normal_mode_ring.particle_array[i].get_velocity();
	vel+=(dt/2.0L)*normal_mode_ring.particle_array[i].get_acceleration();
	normal_mode_ring.particle_array[i].set_velocity(vel);
	if (detach_centroid_thermostat_==0 || i!=0){
	    normal_mode_ring.particle_array[i].integrate_thermostat(dt/2.0L);
	}
    }

    //update cartesian coordinate-necessary to compute the conserved hamiltonian
    for (int i=1; i<=normal_mode_ring.get_num_beads(); ++i){
	cartesian_ring.particle_array[i-1].set_velocity(0.0L);
	cartesian_ring.particle_array[i-1].set_location(0.0L);
	long double vel=0.0L;
	long double loc=0.0L;
	for (int k=1; k<=normal_mode_ring.get_num_beads(); ++k){
	    vel+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_velocity();
	    loc+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_location();
	}
	cartesian_ring.particle_array[i-1].set_velocity(vel);
	cartesian_ring.particle_array[i-1].set_location(loc);
    }

    delete nm_forces;

}


void simulation_setup::read_starting_configuration(long double gamma, int detach_cthermo, 
						   int rstrtPIMD, necklace normal_mode_ring, 
						   necklace cartesian_ring, transition_matrix tm,
						   int potential_choice){

    gamma_=gamma;
    detach_centroid_thermostat_=detach_cthermo;
    restart_PIMD_=rstrtPIMD;

    potential_choice_=potential_choice;
    //if gamma (the adiabaticity parameter) is less than 1 and the thermostat
    //has NOT been detached, this is the CMD-NVT ensemble. The configuration
    //file is from a PIMD simulation
    //
    //if the restartPIMD flag is set, the user is specifying to continue
    //the previous PIMD simulation from its final coordinates.

    if ((gamma<1 && detach_cthermo==0) || rstrtPIMD==1){
 	std::ifstream input_file;
	input_file.open("PIMD.out");
	int count=0;
	while (count<normal_mode_ring.get_num_beads()){
	    int index=0;
	    long double pos=0.0L;
	    long double vel=0.0L;
	    input_file>>index;
	    input_file>>pos;
	    input_file>>vel;
	    normal_mode_ring.particle_array[index].set_location(pos);
	    normal_mode_ring.particle_array[index].set_velocity(vel);
	    ++count;
	}
    }				   

    //if gamma is less than one and the thermostat is detached from the centroid,
    //then this is a CMD-NVE simulation
    else if (gamma<1 && detach_cthermo==1){
	std::ifstream input_file;
	input_file.open("CMD.out");
	int count=0;
	while (count<normal_mode_ring.get_num_beads()){
	    int index=0;
	    long double pos=0.0L;
	    long double vel=0.0L;
	    input_file>>index;
	    input_file>>pos;
	    input_file>>vel;
	    normal_mode_ring.particle_array[index].set_location(pos);
	    normal_mode_ring.particle_array[index].set_velocity(vel);
	    ++count;
	}
    }

    //if either of the above conditions are met, the cartesian positions and
    //velocities will also need to be updated to match

    if (gamma<1.0 || rstrtPIMD==1){
	for (int i=1; i<=cartesian_ring.get_num_beads(); ++i){
	    cartesian_ring.particle_array[i-1].set_velocity(0.0L);
	    cartesian_ring.particle_array[i-1].set_location(0.0L);
	    long double temp_vel=0.0L;
	    long double temp_loc=0.0L;
	    for (int k=1; k<=normal_mode_ring.get_num_beads(); ++k){
		temp_vel+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_velocity();
		temp_loc+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_location();
	    }
	    cartesian_ring.particle_array[i-1].set_velocity(temp_vel);
	    cartesian_ring.particle_array[i-1].set_location(temp_loc);
	}
    }
}
void simulation_setup::print_captured_trajectory(long double time_stamp, necklace normal_mode_ring,
						 virial_estimator virial, barker_estimator barker, 
						 long double conserved, long double analytical_solution){

    std::cout
	<<std::setw(18)<<time_stamp
	//the centroid kinetic and potential energies are conserved throughout the simulation (as detaching it from the thermostat places it in the microcanonical ensemble)
	<<std::setw(18)<<normal_mode_ring.particle_array[0].get_location()
        <<std::setw(18)<<normal_mode_ring.particle_array[0].get_velocity()
	<<std::setw(18)<<conserved;
 
    if (potential_choice_==1){
	std::cout<<
	    std::setw(18)<<analytical_solution<<
	    std::setw(18)<<barker.get_avg_total_E()<<
            std::setw(18)<<virial.get_average_energy();
    }
    std::cout
	<<'\n';

}


void simulation_setup::print_final_configuration(necklace normal_mode_ring){

    if (gamma_==1 && detach_centroid_thermostat_==0){
    	std::ofstream output_file;
	output_file.open("PIMD.out");
	for (int i=0; i<normal_mode_ring.get_num_beads(); ++i){
	    output_file
		<<std::setw(18)<<i
		<<std::setw(18)<<normal_mode_ring.particle_array[i].get_location()
		<<std::setw(18)<<normal_mode_ring.particle_array[i].get_velocity()
		<<'\n';
	}
    }

    else if(gamma_<1 && detach_centroid_thermostat_==0){
	std::ofstream output_file;
	output_file.open("CMD.out");
	for (int i=0; i<normal_mode_ring.get_num_beads(); ++i){
	    output_file
		<<std::setw(18)<<i
		<<std::setw(18)<<normal_mode_ring.particle_array[i].get_location()
		<<std::setw(18)<<normal_mode_ring.particle_array[i].get_velocity()
		<<'\n';
	}
    }
}
