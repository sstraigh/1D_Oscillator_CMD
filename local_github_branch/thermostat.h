#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>

class thermostat{

public:

    thermostat();

    void initialize_thermostat(int ndof, int chain_length, long double timestep, long double kT, long double particle_velocity, long double particle_mass, long double spring_freq);

    void print_status();
    
    void deallocate();

    long double calculate_Nhamiltonian();
    
    int chain_length; //M in Tuckerman(2010), "size" of the chain
    int degrees_of_freedom;

    long double p_velocity;
    long double p_mass;
    long double tau;
    long double Beta;
    long double Nhamiltonian;

    long double N_pos_H;
    long double N_kin_H;

    long double * nose_accel;
    long double * nose_pos;
    long double * nose_velocity;
    long double * Q_param;

    long double suzuki_coefficients[5];
};

#endif
