#include "readin.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

void readin::read_parameters(int& number_of_steps, int& number_of_beads,
   			     int& thermostat_length, long double& kT,
      			     int& capture_rate, long double& gamma, int& potential_choice,
 			     int& detach_centroid_thermostat, int& continue_previous_sim,
			     long double& hbar, long double& spring_constant,
			     long double& particle_mass, long double& v0,
			     long double& r_init, long double& dt,
			     int& spatial_deg_of_freedom){

    std::ifstream input_file;

    input_file.open("parameters.dat");

    int num_parameters=16;
    int count=1;

    while(!input_file.eof() && count<=num_parameters){

	std::string str;

	double data;

	input_file>>str;
	input_file>>data;

	if (count==1) number_of_steps=data;
	else if (count==2) number_of_beads=data;
	else if (count==3) thermostat_length=data;
	else if (count==4) kT=data;
	else if (count==5) capture_rate=data;
	else if (count==6) gamma=data;
	else if (count==7) potential_choice=data;
	else if (count==8) detach_centroid_thermostat=data;
	else if (count==9) continue_previous_sim=data;
	else if (count==10) hbar=data;
	else if (count==11) spring_constant=data;
	else if (count==12) particle_mass=data;
	else if (count==13) v0=data;
	else if (count==14) r_init=data;
	else if (count==15) dt=data;
	else if (count==16) spatial_deg_of_freedom=data;
	++count;	
    }

}
