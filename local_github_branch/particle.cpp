//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 3 Sep, 2014


#include "particle.h"
#include <cmath>

//Particle class-responsible for instantiating a "bead", which is contained in the 
//"necklace" class. Particle also attached a Nose-Hoover chain to each bead's degree
//of freedom (1 in the 1D case without holonomic constraints); the integrate_thermostat
//routine ensures proper velocity scaling with the attached thermostat within the
//framework of the Velocity-verlet algorithm

void particle::deallocate(){
    nhc.deallocate();
}


void particle::initialize(long double mas, long double displacement, long double spr_constant, long double vel, int ndof, int chain_length, long double timestep, long double kT, long double spring_freq){
    omega_k=spring_freq;
    mass=mas;
    velocity=vel;
    location=displacement;
    Beta=1.0L/kT;
    nhc.initialize_thermostat(ndof, chain_length, timestep, kT, velocity, mass, spring_freq);
}

//Integration scheme adapted from source code written by F. Paesani
//For further reference on the Suzuki-Yoshida integration scheme, see
//     ref 1
//     ref 2
//     ref 3
//     ref 4


void particle::integrate_thermostat(long double timestep){

    long double velocity_scale_factor=1.0L;

    long double Kin_E_2 = mass*velocity*velocity;

    for (int i=0; i<5; ++i){

	long double eps = timestep * nhc.suzuki_coefficients[i];
	long double eps2= eps * 0.5L;
	long double eps4= eps * 0.25L;

	//Update force acting upon the first nose particle
	//using 2*kinetic_energy of the particle it is thermostatting
	nhc.nose_accel[1]= (Kin_E_2 - (1.0L/Beta))/nhc.Q_param[1];

	//update thermostat velocities for the first half of the integration scheme
	nhc.nose_velocity[nhc.chain_length]+=(eps2*nhc.nose_accel[nhc.chain_length]);
	for (int j=(nhc.chain_length-1); j>=1; --j){
	    long double val=std::exp(-1.0L*eps4*nhc.nose_velocity[j+1]);
	    nhc.nose_velocity[j]=val*(val*nhc.nose_velocity[j] + eps2*nhc.nose_accel[j]);
	}

	//update nose particle positions
	for (int j=1; j<=nhc.chain_length; ++j) nhc.nose_pos[j]+=(eps*nhc.nose_velocity[j]);

	//update the system (i.e. "bead") velocity scaling factor
	long double temp_val= std::exp(-1.0L*eps*nhc.nose_velocity[1]);
	velocity_scale_factor*=temp_val;
	Kin_E_2*=(temp_val*temp_val);

	//update the nose particle velocities for the second half of the integration
	//scheme

	//update first nose particle velocity before the others, as it is dependent on the system
	nhc.nose_accel[1] = (Kin_E_2 - (1.0L/Beta))/nhc.Q_param[1];
	long double temp=std::exp(-1.0L*eps4*nhc.nose_velocity[2]);
	nhc.nose_velocity[1]=temp*(temp*nhc.nose_velocity[1] + eps2*nhc.nose_accel[1]);

	//update middle nose particle velocities
	for (int j=2; j<nhc.chain_length; ++j){
	    nhc.nose_accel[j]=((nhc.Q_param[j-1]*nhc.nose_velocity[j-1]*nhc.nose_velocity[j-1] - (1/Beta))/nhc.Q_param[j]);
	    temp=std::exp(-1.0L*eps4*nhc.nose_velocity[j+1]);
	    nhc.nose_velocity[j]=temp*(temp*nhc.nose_velocity[j] + eps2*nhc.nose_accel[j]);
	}

	//update last nose particle velocity
	nhc.nose_accel[nhc.chain_length]=((nhc.Q_param[nhc.chain_length-1]*nhc.nose_velocity[nhc.chain_length-1]*nhc.nose_velocity[nhc.chain_length-1] - (1/Beta))/nhc.Q_param[nhc.chain_length]);
	nhc.nose_velocity[nhc.chain_length]+=(eps2*nhc.nose_accel[nhc.chain_length]);

    }
    //update system velocity
    velocity*=velocity_scale_factor;
}
