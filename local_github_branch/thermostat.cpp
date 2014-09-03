//Author: Shelby Straight, Paesani Lab
//Date last edited: 3, Sep 2014
//
//This is the implementation file for the thermostat object.
//It includes routines to:
//      1. initialize the thermostat
//      2. deallocate the memory assigned to the thermostat variables
//      3. calculate the contribution of the thermostat object to the 
//         extended hamiltonian

#include "thermostat.h"
#include <iomanip>
#include <iostream>
#include "gaussian.h"

void thermostat::deallocate(){

    delete nose_pos;
    delete nose_velocity;
    delete nose_accel;
    delete Q_param;

}

void thermostat::initialize_thermostat(int ndof, int num_pairs, long double timestep, long double kT, long double particle_velocity, long double particle_mass, long double spring_freq){

    chain_length=num_pairs;

    tau=1.0/spring_freq;

    Beta=1.0/kT;

    degrees_of_freedom=ndof;

    p_velocity=particle_velocity;
    p_mass=particle_mass;

    Q_param=new long double[chain_length+1];
    nose_pos=new long double[chain_length+2];
    nose_velocity=new long double[chain_length+2];
    nose_accel=new long double[chain_length+2];

    std::fill(nose_pos, nose_pos+chain_length+2, 0.0);
    std::fill(nose_velocity, nose_velocity+chain_length+2, 0.0);
    std::fill(nose_accel, nose_accel+chain_length+2, 0.0);  
   
    //initialize Q-vector, parameters responsible for adjusting the rate of change of
    //introduced (fake) particle momenta (fake particle masses)
    for (int i=1; i<=chain_length; ++i){
	Q_param[i]=kT*tau*tau;
    }

    //random initialization of velocities using the gaussian object
    //initialization of positions to zero
    gaussian my_gauss;
    for (int i=1; i<=chain_length; ++i){
	long double vel=my_gauss.compute_gaussian((long double) i);
	vel*=std::sqrt(kT/Q_param[i]);
    	nose_velocity[i]=vel;
	nose_pos[i]=0.0;
    }

    //initialize nose particle accelerations
    nose_accel[1]=((p_mass*p_velocity*p_velocity)-kT)/Q_param[1];
    for (int i=2; i<=chain_length; ++i){
	nose_accel[i]=((Q_param[i-1]*nose_velocity[i-1]*nose_velocity[i-1]) - kT ) /Q_param[i];
    }

    //initialize array of suzuki coefficients
    long double power= 1.0L/3.0L;
    long double val= 1.0/(4.0-(std::pow(4, power)));
    for (int i=0; i<5; ++i){
	if (i==2) suzuki_coefficients[i]=(1.0-(4.0*val));
	else suzuki_coefficients[i]=val;
    }

}

long double thermostat::calculate_Nhamiltonian(){

    Nhamiltonian=0.0;
    N_kin_H=0.0L;
    N_pos_H=0.0L;

    for (int i=1; i<=chain_length; ++i){

	N_kin_H+=(0.5*Q_param[i]*nose_velocity[i]*nose_velocity[i]);
	N_pos_H+=(nose_pos[i]/Beta);
    }

    if (degrees_of_freedom>1 && chain_length>=1){
	N_pos_H+=(nose_pos[1]*(degrees_of_freedom-1))/Beta;
    }

    Nhamiltonian = N_kin_H + N_pos_H;

    return Nhamiltonian;
}
