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

//Morse potential expression
const long double alpha=1.1563;
const long double D=0.18748;
const long double r0=1.8325;
const long double frac=7.0/12.0;

int main(int argc, char *argv[]){

    long double c_pot_avg=0.0L;
    long double c_kin_avg=0.0L;

    if (argc<10){
	std::cerr<<"Usage: ./1D_oscillator Number_of_steps Number_of_particles Length_of_Thermostat Temperature(in units of kT) Capture_rate Adiabaticity_parameter Potential_choice Thermostat&Centroid RstrtPIMD"<<std::endl;
	exit(1);
    }
//PARAMETERS FOR EXPERIMENTAL VARIATION:
//     k=spring constant
//     m=mass
//     v0=initial velocity
//     r0=initial displacement
//     dt=timestep=0.01
//     beta=1/kT
//     gamma=adiabadicity parameter
//     omega=normal mode vibrational frequency

    long double plancks_h = 6.626e-34;
    long double hbar= 1.0;

    long double spring_constant=1.0;
    long double particle_mass=1.0;
    long double v0=1.0;
    long double r_init=1.0;
    long double dt=0.0001;
    long double kT, Beta;
    const int spatial_deg_of_freedom=1;
    int nhc_length=0;
    long double avg_barker_kin=0.0;
    long double avg_barker_pot=0.0;

    long double avg_virial_estimator=0.0;

    //read in no. of steps from command line
    std::istringstream iss(argv[1]);
    int number_of_steps=0;
    iss>>number_of_steps;

    //read in no. particles from command line 
    std::istringstream ss(argv[2]);
    int number_of_particles=0;
    ss>>number_of_particles;

    //read in length of nose-hoover chain from command line 
    std::istringstream ssi(argv[3]);
    ssi>>nhc_length;

    //read in temperature from command line
    std::istringstream ssii(argv[4]);
    ssii>>kT;
    Beta=1.0/kT;

    //read in the rate at which the trajectory
    //is forwarded to the output from command line
    int capture_rate=0;
    std::istringstream ssiii(argv[5]);
    ssiii>>capture_rate;

    //read in the adiabaticity parameter from command line
    long double gamma=0.0L;
    std::istringstream sss(argv[6]);
    sss>>gamma;

    int potential_choice=0;
    std::istringstream sis(argv[7]);
    sis>>potential_choice;

    std::cerr<<'\n'<<potential_choice<<'\n';

    //use the adiabaticity parameter to determine
    //whether or not CMD or NMPIMD simulation will be run
    bool detach_cthermo=false;

    int c_thermo=-1;
    std::istringstream siss(argv[8]);
    siss>>c_thermo;

    if (c_thermo==1) detach_cthermo=true;

    int rstrtPIMD=-1;
    std::istringstream sisi(argv[9]);
    sisi>>rstrtPIMD;

    long double * fictitious_masses=new long double[number_of_particles+1];
    long double * normal_mode_coord=new long double[number_of_particles+1]; 
    long double * nm_velocities=new long double[number_of_particles+1];

    //neck_beads is the data structure which holds the particles in its cartesian coordinates
    //it is responsible for computing the forces and potentials (both from the external-ie central-
    //spring and the inter-bead springs)
    
    necklace neck_beads;
    neck_beads.create_necklace(hbar, number_of_particles, particle_mass, r_init, v0, spring_constant, spatial_deg_of_freedom, nhc_length, dt, kT, potential_choice);

    //if an NVT ensemble using PIMD has already been run, the final positions and velocities of the ensemble
    //will be stored in the file PIMD.out
    //
    //These positions and velocities must be used in a CMD simulation in order to guarantee a proper sampling
    //of the microcanonical
    //ensemble for the centroid variable\

    long double * avg_ext_pot=new long double[number_of_particles];
    long double * avg_kin=new long double[number_of_particles];

    std::fill(avg_ext_pot, avg_ext_pot+number_of_particles, 0.0L);
    std::fill(avg_kin, avg_kin+number_of_particles, 0.0L);

    long double root_beads=std::sqrt(number_of_particles);

    transition_matrix tm;
    tm.convert_coordinates((root_beads*kT), number_of_particles, particle_mass, gamma);

    //    fictitious masses are used for the propagation of dynamics in the PIMD/CMD representation
    fictitious_masses[1]=particle_mass;

    for (int i=2; i<=number_of_particles; ++i){
	fictitious_masses[i]=particle_mass*tm.get_eigenvalue(i)*gamma*gamma;
    }

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
    
    //if gamma<1, we are reading in avg NM pot_E and kin_E (without spring contributions)

    if ((gamma<1 && detach_cthermo==false) || rstrtPIMD==1){
	std::ifstream input_file;
	input_file.open("PIMD.out");

	int count=0;

	while (count<number_of_particles){
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

    else if (gamma<1 && detach_cthermo==true){
	std::ifstream input_file;
	input_file.open("CMD.out");

	int count=0;

	while (count<number_of_particles){

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

    //if doing a CMD simulation (or if restarting a PIMD simulation), cartesian positions need to be updated before propagation
    if (gamma<1.0 || rstrtPIMD==1){
	for (int i=1; i<=number_of_particles; ++i){
	    neck_beads.particle_array[i-1].set_velocity(0.0L);
	    neck_beads.particle_array[i-1].set_location(0.0L);

	    long double temp_vel=0.0L;
	    long double temp_loc=0.0L;

	    for (int k=1; k<=number_of_particles; ++k){
		temp_vel+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_velocity();
		temp_loc+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_location();
	    }

	    neck_beads.particle_array[i-1].set_velocity(temp_vel);
	    neck_beads.particle_array[i-1].set_location(temp_loc);
	}
    }

    //The analytical solution to the harmonic oscillator at temperature Beta=1/kT is calculated from the 
    //partition function of the harmonic oscillator
    long double analytical_solution=0.5*hbar*spring_constant;
    analytical_solution+=(hbar*spring_constant*std::exp(-1.0*Beta*hbar*spring_constant))/(1-std::exp(-1.0*Beta*hbar*spring_constant));

    //This loop starts the propagation of the simulation
    for (int j=1; j<=number_of_steps; ++j){

	long double Total_E=0.0;
	long double Kin_E=0.0;
	long double Pot_E=0.0;
	long double conserved=0.0;
	long double NHam=0.0;

	long double Npos=0.0L;
        long double Nkin=0.0L;

	long double * nm_forces=new long double[number_of_particles+1];

	//compute initial forces in NM representation
	for (int i=1; i<=number_of_particles; ++i){
	    nm_forces[i]=0.0L;
	    for (int k=1; k<=number_of_particles; ++k){
		nm_forces[i]+=tm.get_cart_to_nm(i,k)*neck_beads.compute_bead_force(k-1)*number_of_particles;
	    }
	}

	//1st part of Velocity Verlet algorithm
	//(note: for CMD, there is no thermostat, so the detach_cthermo flag ensures that no discontinuities are
	//fed into the calculations)
	for (int i=0; i<number_of_particles; ++i){
	    long double acc=nm_forces[i+1]/fictitious_masses[i+1];
	    normal_mode_ring.particle_array[i].set_acceleration(acc);
	    if (detach_cthermo==false || i!=0) normal_mode_ring.particle_array[i].integrate_thermostat(dt/2.0L);
	}

	//CAN THIS LOOP BE ROLLED INTO THE ONE ABOVE??
	for (int i=0; i<number_of_particles; ++i){
	    long double vel=normal_mode_ring.particle_array[i].get_velocity();
	    vel += (dt/2.0L)*normal_mode_ring.particle_array[i].get_acceleration();
	    normal_mode_ring.particle_array[i].set_velocity(vel);

	    long double loc=normal_mode_ring.particle_array[i].get_location();
	    loc+=dt*normal_mode_ring.particle_array[i].get_velocity();
	    normal_mode_ring.particle_array[i].set_location(loc);
	}

	//update cartesian positions in order to recompute forces
	for (int i=1; i<=number_of_particles; ++i){
	    neck_beads.particle_array[i-1].set_location(0.0L);

	    long double loc=0.0L;
	    for (int k=1; k<=number_of_particles; ++k){
		loc+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_location();
	    }

	    neck_beads.particle_array[i-1].set_location(loc);
	}

	//compute updated forces in NM representation
	//CAN THIS LOOP BE ROLLED INTO THE ONE ABOVE??
        for (int i=1; i<=number_of_particles; ++i){
	    nm_forces[i]=0.0L;
	    for (int k=1; k<=number_of_particles; ++k){
		nm_forces[i]+=(tm.get_cart_to_nm(i,k)*neck_beads.compute_bead_force(k-1))*number_of_particles;
	    }
	}	

	//continue with the second portion of VV
	for (int i=0; i<number_of_particles; ++i){
	    normal_mode_ring.particle_array[i].set_acceleration(nm_forces[i+1]/fictitious_masses[i+1]);

	    long double vel=normal_mode_ring.particle_array[i].get_velocity();
	    vel+=(dt/2.0L)*normal_mode_ring.particle_array[i].get_acceleration();
	    normal_mode_ring.particle_array[i].set_velocity(vel);

	    if (detach_cthermo==false || i!=0){
	
		normal_mode_ring.particle_array[i].integrate_thermostat(dt/2.0L);
	
	    }
	}
	//update cartesian coordinates
  	for (int i=1; i<=number_of_particles; ++i){
	    neck_beads.particle_array[i-1].set_velocity(0.0L);
	    neck_beads.particle_array[i-1].set_location(0.0L);

	    long double vel=0.0L;
	    long double loc=0.0L;
	    for (int k=1; k<=number_of_particles; ++k){
		vel+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_velocity();
		loc+=tm.get_nm_to_cart(i,k)*normal_mode_ring.particle_array[k-1].get_location();
	    }

	    neck_beads.particle_array[i-1].set_velocity(vel);
	    neck_beads.particle_array[i-1].set_location(loc);

	}
	//use the thermostat hamiltonian and kinetic energy contributions from the 
	//ring representation in the normal mode coordinates...
	NHam    = normal_mode_ring.compute_NHC_hamiltonian();
	Total_E = normal_mode_ring.compute_total_hamiltonian();
	Kin_E   = normal_mode_ring.get_necklace_KE();
	Pot_E   = normal_mode_ring.get_necklace_PE();

	//...and the potential energy contributions (both internal necklace-spring
	//and external portions) to compute the conserved extended hamiltonian
	long double Total_E1 = neck_beads.compute_total_hamiltonian();
	long double bead_pot = neck_beads.get_necklace_bead_potential();
	long double ext_pot  = neck_beads.get_necklace_classical_potential();


/*	
	//calculation of the barker energy estimator (for kinetic energies)
	//this algorithm computes the summation term (eq 3.4, Cao and Berne, 1989) using the primitive algorithm
	long double pre_factor=particle_mass*number_of_particles*kT*kT/(2*hbar*hbar);
	long double summation=0.0;
	for (int i=0; i<(number_of_particles-1); ++i) summation+=(neck_beads.particle_array[i].location-neck_beads.particle_array[i+1].location)*(neck_beads.particle_array[i].location-neck_beads.particle_array[i+1].location);
	summation+=(neck_beads.particle_array[number_of_particles-1].location-neck_beads.particle_array[0].location)*(neck_beads.particle_array[number_of_particles-1].location-neck_beads.particle_array[0].location);
	long double barker_kin_est=pre_factor*summation;
	barker_kin_est=(number_of_particles*spatial_deg_of_freedom*kT/2.0)-barker_kin_est;

	//compute running average of the kinetic term of the barker estimator
	avg_barker_kin*=(j-1);	
	avg_barker_kin+=barker_kin_est;
	avg_barker_kin/=j;

	//compute the potential term of the barker estimator
	long double barker_potential_estimator=0.0;
	for (int i=0; i<number_of_particles; ++i) barker_potential_estimator+=0.5*spring_constant*neck_beads.particle_array[i].location*neck_beads.particle_array[i].location/((long double) number_of_particles);

	//update the running average of the barker potential energy estimator
	avg_barker_pot*=(j-1);
	avg_barker_pot+=barker_potential_estimator;
	avg_barker_pot/=j;


	//calculation of the virial energy estimator
	long double multiplier=1.0/((long double)number_of_particles);
	long double summation_term=0.0;
	for (int i=0; i<number_of_particles; ++i){
	   summation_term+=spring_constant*neck_beads.particle_array[i].location*neck_beads.particle_array[i].location;
	}
	long double virial_estimator=multiplier*summation_term;

	//update the running average of the virial estimator
	avg_virial_estimator*=(j-1);
	avg_virial_estimator+=virial_estimator;
	avg_virial_estimator/=j;
*/
	//compute the extended hamiltonian regardless of NM vs. primitive
    	conserved=(Kin_E+bead_pot+ext_pot+NHam);

	
	/*
	long double c_kin, c_pot;

	if (j%capture_rate==0){

	c_kin=normal_mode_ring.compute_bead_KE(0);
	c_kin_avg*=((j/capture_rate)-1);
	c_kin_avg+=c_kin;
	c_kin_avg/=(j/capture_rate);

	c_pot=normal_mode_ring.compute_external_potential(0)*number_of_particles;
	c_pot_avg*=((j/capture_rate)-1);
	c_pot_avg+=c_pot;
        c_pot_avg/=(j/capture_rate);

	}
	//for some reason the following definitions of the potential applied to the centroid are not consistent (but they do converge?) throughout the course of the CMD simulaition
	//     originally these were written using the cartesian variables. Lets try these using the normal modes?

	for (int i=0; i<number_of_particles; ++i){

	    if (j%capture_rate==0){
	    avg_ext_pot[i]*=((j/capture_rate)-1);
	    long double pot_En=normal_mode_ring.compute_external_potential(i)*number_of_particles;
	    avg_ext_pot[i]+=pot_En;
	    avg_ext_pot[i]/=(j/capture_rate);

	    long double kin_En=normal_mode_ring.compute_bead_KE(i);
	    avg_kin[i]*=((j/capture_rate)-1);
	    avg_kin[i]+=kin_En;
	    avg_kin[i]/=(j/capture_rate);
	    }
	}
*/

	delete nm_forces;
	//this statement relays all new computed information to the screen at the given capture rate
	if (j%capture_rate==0){
    	    std::cout
    		<<std::setw(18)<<(j*dt)

		//the centroid kinetic and potential energies are conserved throughout the simulation (as detaching it from the thermostat places it in the microcanonical ensemble)
/*		<<std::setw(18)<<c_kin
		<<std::setw(18)<<c_kin_avg
       		<<std::setw(18)<<c_pot
		<<std::setw(18)<<c_pot_avg
*/
		<<std::setw(18)<<normal_mode_ring.particle_array[0].get_location()
	   	<<std::setw(18)<<conserved
		<<std::setw(18)<<normal_mode_ring.particle_array[0].get_velocity()
		<<std::setw(18)<<(avg_barker_pot+avg_barker_kin)
    		<<std::setw(18)<<avg_virial_estimator
    		<<std::setw(18)<<analytical_solution
    		<<'\n';
/*
	    for (int i=0; i<number_of_particles; ++i){

		std::cerr<<
		    std::setw(18)<<normal_mode_ring.particle_array[i].location;
	
   	    }

	    std::cerr<<std::endl;
*/		
	}
    }
/*
    if (gamma==1 && detach_cthermo==false){
	std::ofstream output_file;
	output_file.open("PIMD.out");
	
	for (int i=0; i<number_of_particles; ++i){
	    output_file
		<<std::setw(18)<<i
		<<std::setw(18)<<normal_mode_ring.particle_array[i].location
		<<std::setw(18)<<normal_mode_ring.particle_array[i].velocity
		<<'\n';
	}
    }

    else if(gamma<1 && detach_cthermo==false){
	std::ofstream output_file;
	output_file.open("CMD.out");
	for (int i=0; i<number_of_particles; ++i){
	    output_file
		<<std::setw(18)<<i
		<<std::setw(18)<<normal_mode_ring.particle_array[i].location
		<<std::setw(18)<<normal_mode_ring.particle_array[i].velocity
		<<'\n';
	}
    }
*/
    delete fictitious_masses;
    delete normal_mode_coord;
    delete nm_velocities;
    delete avg_ext_pot;
    delete avg_kin;

    tm.deallocate();
    neck_beads.deallocate();
    normal_mode_ring.deallocate();



    return 0;

}
