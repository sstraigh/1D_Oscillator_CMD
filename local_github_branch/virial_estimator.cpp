//Author: Shelby Straight, Paesani Lab
//Date last modified: 4 Sep, 2014
//

#include "virial_estimator.h"

void virial_estimator::initialize_estimator(int number_of_beads, 
					    long double spring_constant){
    avg_estimated_energy_=0.0L;
    instantaneous_estimated_energy_=0.0L;
    multiplier_=1.0L/((long double) number_of_beads);

    spring_constant_=spring_constant;
}

void virial_estimator::update_energy(necklace my_necklace, int weight_of_average){

    long double summation_term=0.0;
    
    for (int i=0; i<my_necklace.get_num_beads(); ++i){
	summation_term+=spring_constant_*my_necklace.particle_array[i].get_location()
	    *my_necklace.particle_array[i].get_location();
    }

    instantaneous_estimated_energy_=summation_term*multiplier_;

    avg_estimated_energy_*=(weight_of_average-1);
    avg_estimated_energy_+=instantaneous_estimated_energy_;
    avg_estimated_energy_/=weight_of_average;

}

