#include "barker_estimator.h"

void barker_estimator::initialize_estimator(int number_of_beads, long double spring_constant, long double particle_mass, long double kT, long double hbar){
    avg_kinetic_energy_=0.0L;
    avg_potential_energy_=0.0L;
    instant_kinetic_energy_=0.0L;
    instant_potential_energy_=0.0L;

    spring_constant_=spring_constant;
    kinetic_prefactor_=particle_mass*number_of_beads*kT*kT/(2*hbar*hbar);

}

void barker_estimator::update_energy(necklace my_necklace, int weight_of_average){

    long double summation_term=0.0L;

    for (int i=0; i<(my_necklace.get_num_beads()-1); ++i){
	summation_term+=((my_necklace.particle_array[i].get_location()-my_necklace.particle_array[i+1].get_location())
	    *(my_necklace.particle_array[i].get_location()-my_necklace.particle_array[i+1].get_location()));
    }

    summation_term+=(my_necklace.particle_array[my_necklace.get_num_beads()-1].get_location()-my_necklace.particle_array[0].get_location())
	*(my_necklace.particle_array[my_necklace.get_num_beads()-1].get_location()-my_necklace.particle_array[0].get_location());

    instant_kinetic_energy_=kinetic_prefactor_*summation_term;

    avg_kinetic_energy_*=(weight_of_average-1);
    avg_kinetic_energy_+=instant_kinetic_energy_;
    avg_kinetic_energy_/=weight_of_average;

    instant_potential_energy_=0.0L;
    for (int i=0; i<my_necklace.get_num_beads(); ++i){
	instant_potential_energy_+=0.5*spring_constant_*my_necklace.particle_array[i].get_location()*my_necklace.particle_array[i].get_location()/((long double) my_necklace.get_num_beads());
    }

    avg_potential_energy_*=(weight_of_average-1);
    avg_potential_energy_+=instant_potential_energy_;
    avg_potential_energy_/=weight_of_average;
}
