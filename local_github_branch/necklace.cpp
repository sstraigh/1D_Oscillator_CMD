#include "necklace.h"
#include "gaussian.h"

void necklace::deallocate(){

    for (int i=0; i<number_of_beads; ++i){
	particle_array[i].deallocate();
    }

    delete particle_array;
}

long double necklace::compute_external_potential(int bead_index){

    long double classical_PE=0.0L;
    //harmonic potential
    if (potential_choice==1) classical_PE=particle_array[bead_index].get_location()*particle_array[bead_index].get_location()*spr_const*0.5/((long double)number_of_beads);
     
    //morse potential
    //build in this later once parameters are specified
    else if (potential_choice==2) classical_PE=0.0L;

    //Anharmonic-a
    if (potential_choice==3){
	classical_PE+=(0.5 *particle_array[bead_index].get_location()*particle_array[bead_index].get_location());
	classical_PE+=(0.1 *particle_array[bead_index].get_location()*particle_array[bead_index].get_location()*particle_array[bead_index].get_location());
	classical_PE+=(0.01*particle_array[bead_index].get_location()*particle_array[bead_index].get_location()*particle_array[bead_index].get_location()*particle_array[bead_index].get_location());
	classical_PE/=number_of_beads;
    } 
    //Anharmonic-b
    if (potential_choice==4){
	classical_PE=(0.25*particle_array[bead_index].get_location()*particle_array[bead_index].get_location()*particle_array[bead_index].get_location()*particle_array[bead_index].get_location())/((long double)number_of_beads);
    }

    return classical_PE;

}

long double necklace::compute_bead_KE(int bead_index){

    return (0.5*particle_array[bead_index].get_velocity()*particle_array[bead_index].get_velocity()*particle_array[bead_index].get_mass());

}

long double necklace::compute_bead_PE(int bead_index){

    int index_below= (bead_index-1);
    if (index_below<0) index_below=number_of_beads-1;

    int index_above=bead_index+1;
    if (index_above==number_of_beads) index_above=0;

    long double classical_PE=this->compute_external_potential(bead_index);

    //Quantum interaction potential energy referes to the potential added to the particle in addition to the classical potential from the potential connecting the necklace and its
    //equilibrium distance.
    
    //This potential takes into account the interaction of the bead with its nearest neighbors and the springs interconnecting them
    long double quantum_interaction_PE=0.0;

    long double pre_factor=0.5*bead_freq*bead_freq*particle_array[bead_index].get_mass();

    if (number_of_beads>1){
	quantum_interaction_PE = (pre_factor*(particle_array[index_above].get_location()-particle_array[bead_index].get_location())*(particle_array[index_above].get_location()-particle_array[bead_index].get_location()));
    }

    bead_potential+=quantum_interaction_PE;

    classical_potential+=classical_PE;

    return (classical_PE + quantum_interaction_PE);

}

long double necklace::compute_external_force(int bead_index){

    long double classical_force=0.0L;

    if (potential_choice==1) classical_force=(-1.0*spr_const*particle_array[bead_index].get_location())/number_of_beads;

    else if (potential_choice==2) classical_force=0.0L;

    else if (potential_choice==3){
	classical_force= (-1.0L*particle_array[bead_index].get_location());
	classical_force-=(0.3L*particle_array[bead_index].get_location() * particle_array[bead_index].get_location());
	classical_force-=(0.04L*particle_array[bead_index].get_location()* particle_array[bead_index].get_location() *particle_array[bead_index].get_location());
	classical_force/=number_of_beads;
    }

    else if (potential_choice==4){
	classical_force=(-1.0L *particle_array[bead_index].get_location()* particle_array[bead_index].get_location() *particle_array[bead_index].get_location())/number_of_beads;
    }
    return classical_force;
}

long double necklace::compute_bead_force(int bead_index){

    //The force acting on the bead at any time is the 1st derivative of the potential with respect to the position of the bead (since the potential is constant with respect to time)
    int index_below= (bead_index-1);
    if (index_below<0) index_below=number_of_beads-1;

    int index_above= (bead_index+1);
    if (index_above==number_of_beads) index_above=0;

    long double classical_force=this->compute_external_force(bead_index);

    long double spring_force=0.0;
    //Spring force takes into account the two nearest neighbor interactions between the beads-i.e., the X(k) bead feels the springs connecting it
    //to both the X(k+1) and X(k-1) beads. In the case of num_beads==1, this reduces to no net spring force, and the "classical" potential is 
    //recovered

    //bead freq is dependent on the index of the bead in question (it is either omega_P or the decoupled omega_P/gamma) however, forces
    //are computed in the cartesian coordinates, where the frequencies are all the same

    spring_force  =(1.0)*bead_freq*bead_freq*((particle_array[index_above].get_location()-particle_array[bead_index].get_location())*particle_array[bead_index].get_mass());
    spring_force+=(-1.0)*bead_freq*bead_freq*((particle_array[bead_index].get_location()-particle_array[index_below].get_location())*particle_array[bead_index].get_mass());

    long double total_force= (classical_force + spring_force);

    return total_force;
}

long double necklace::compute_total_hamiltonian(){

    //Computes the total energy associated with the necklace, accounting for all kinetic interactions
    //and forces arising from the bead_springs and the spring connecting the necklace to its center
    long double total_E=0.0;

    kinetic_E=0.0;
    potential_E=0.0;

    bead_potential=0.0L;
    classical_potential=0.0L;

    for (int i=0; i<number_of_beads; ++i){

	potential_E+=(this->compute_bead_PE(i));
	kinetic_E+=(this->compute_bead_KE(i));

    }

    return (potential_E+kinetic_E);

}

long double necklace::compute_NHC_hamiltonian(){

    //This function uses a subroutine in the thermostat class
    //to compute the contribution of the nose-hoover masses to
    //the overall conserved Nose Hamiltonian

    long double NHC_contributions=0.0;

    for (int i=0; i<number_of_beads; ++i){

	NHC_contributions += particle_array[i].nhc.calculate_Nhamiltonian();
    }

    nhc_contribution=NHC_contributions;

    return NHC_contributions;

}

void necklace::normal_modes(long double hbar, int num_beads, long double * masses, long double * initial_locations, long double * initial_velocities, long double * spring_const, int ndof, int nh_chainlength, long double timestep, long double kT, bool activate_CMD, int pot_choice){
  
    potential_choice=pot_choice;

    //necklace which uses arrays for initialization rather than scalars.
    //specifically designed for usage with the normal modes
    particle_array=new particle[num_beads];
    number_of_beads=num_beads;

    spr_const=1.0L;

    bead_freq=(std::sqrt(number_of_beads)*kT)/hbar;

    for (int i=0; i<num_beads; ++i){
	if (i==0 && activate_CMD==true){
	    particle_array[i].initialize(masses[i+1], initial_locations[i+1], spr_const, initial_velocities[i+1], ndof, 0, timestep, kT, spring_const[i+1]);
	}
	else{
	    particle_array[i].initialize(masses[i+1], initial_locations[i+1], spr_const, initial_velocities[i+1], ndof, nh_chainlength, timestep, kT, spring_const[i+1]); 
	}
    }
    
}

void necklace::create_necklace(long double hbar, int num_beads, long double mass, long double initial_location, long double initial_velocity, long double spring_const, int ndof, int nh_chainlength, long double timestep, long double kT, int pot_choice){

    potential_choice=pot_choice;

    //A new necklace (i.e., an array holding particle objects, each particle is coupled to a thermostat)
    //is instantiated. As a result, the thermostat and particle objects are also created within the array

    particle_array=new particle[num_beads];

    number_of_beads=num_beads;

    //overall potential [U(x)] felt by the particles is modeled by the external potential of the necklace to its center
    //and by the "bead_freq" which regulates the interaction of the particles within the necklace with each other
    spr_const=1.0;

    //m_prime does not effect the thermodynamics of the system and can be set as any value
    m_prime= (mass*number_of_beads)/(2*M_PI*hbar*2*M_PI*hbar);

    physical_mass=mass;

    bead_freq=std::sqrt(number_of_beads)*kT/hbar;

    gaussian my_distribution;

    for (int i=0; i<num_beads; ++i){
	long double vel=my_distribution.compute_gaussian((double)i);
	vel*=std::sqrt(kT/mass);
	long double pos=my_distribution.compute_gaussian(initial_location);
	particle_array[i].initialize(mass, initial_location, spring_const, vel, ndof, nh_chainlength, timestep, kT, bead_freq);
    }
}
