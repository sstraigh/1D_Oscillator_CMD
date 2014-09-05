#include "particle_tracker.h"

void particle_tracker::initialize(int number_of_particles){

    avg_kinetic_energies_=new long double[number_of_particles];
    avg_potential_energies_=new long double[number_of_particles]; 
    instant_PE_= new long double[number_of_particles];
    instant_KE_= new long double[number_of_particles];
    instant_positions_= new long double[number_of_particles];
    instant_velocities_ = new long double[number_of_particles];
}

void particle_tracker::deallocate(){

    delete avg_kinetic_energies_;
    delete avg_potential_energies_;
    delete instant_PE_;
    delete instant_positions_;
    delete instant_KE_;
    delete instant_velocities_;

}

void particle_tracker::print_xyz(int num_particles){

    char atom_symbol='C';

    std::cerr
	<<num_particles<<"\nComment_Line\n";

    for (int i=0; i<num_particles; ++i){
	std::cerr
	    <<atom_symbol
	    <<std::setw(18)<<instant_positions_[i]
	    <<std::setw(18)<<"0.000"
	    <<std::setw(18)<<"0.000\n";
    }
}

void particle_tracker::update_tracker(necklace my_necklace, int weight_of_average){

    for (int i=0; i<my_necklace.get_num_beads(); ++i){

	avg_kinetic_energies_[i]*=(weight_of_average-1);
	avg_potential_energies_[i]*=(weight_of_average-1);

	instant_PE_[i]=my_necklace.compute_external_potential(i)*my_necklace.get_num_beads();
	instant_KE_[i]=my_necklace.compute_bead_KE(i);

	instant_positions_[i]=my_necklace.particle_array[i].get_location();
	instant_velocities_[i]=my_necklace.particle_array[i].get_velocity();

	avg_kinetic_energies_[i]+=instant_KE_[i];
	avg_potential_energies_[i]+=instant_PE_[i];

	avg_kinetic_energies_[i]/=weight_of_average;
	avg_potential_energies_[i]/=weight_of_average;

    }

    print_xyz(my_necklace.get_num_beads());

}
			 
