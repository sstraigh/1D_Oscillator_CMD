1D_Oscillator_CMD
=================

There are three executables in the software suite:

1D_Oscillator: 

     usage: (from command line, example)

          "./1D_Oscillator 1000 4 4 1.0 5 1.0 1 0 0 > STDOUT 2> STDERR"

               Runs a simulation of a 1D-particle attracted to a central potential. 1st, a NMPIMD simulation is run to equilibrate the system. Once
               an ergodic trajectory is confirmed (by monitoring the average Kinetic and Potential energy values of the system throughout the 
               simulation), a 2nd CMD-NVT simulation is run-this is done by changing the value of the adiabaticity parameter from 1.0 to a value less
	       than one, typically 1/P. Following this first Centroid molecular dynamics simulation, a second CMD simulation is run where the
	       thermostat has been detached from the centroid, which effectively places the Centroid in the microcanonical ensemble while leaving the 
	       other (now adiabatically decoupled) normal modes in the canonical ensemble.

	       The NMPIMD simulation produces an output file called PIMD.out-this is an instantenous snapshop of the system dynamics at the time the
	       simulation ended. This is the file that is used to either restart another NMPIMD simulation (if the Restart_Simulation flag is enabled,
	       see below), or that will be used to continue the simulation of the system as it progresses to the 1st CMD simulation. The CMD simulation
	       outputs a file called CMD.out, which, analogously to PIMD.out, is an instantaneous snapshot of the system at the end of the CMD simulation.
	       This can either be used to run another simulation of the system in the CMD-NVT ensemble (if Detach_centroid_thermostat remains 0 and 
	       Restart_simulation is set to 1) or it can be used to start a simulation of the CMD-NVE ensemble (if the above flags are set to 1 and 0,
	       respectively).

     Input Parameters from Command line: Number_of_Steps, Number_of_particles, Thermostat_Length, Temperature, Capture_rate,
                                         Adiabaticity_parameter, Potential_choice, Detach_centroid_thermostat, Restart_Simulation

		Parameter Description: Number_of_steps: integer, details how long the simulation will run
		                       
		                       Number_of_particles: integer (a power of two is preferred), the number of beads used to 
				                            model the particle in the central potential

	                               Thermostat_Length: the length of the Nose-Hoover thermostat chain-see ref >>> for more details

				       Temperature: temperature of the simulation in units of kT

				       Capture_rate: the rate at which trajectory snapshots will be taken and relayed to STDOUT

				       Adiabaticity_paramter: a double between >0 and <=1.0; if the parameter is taken to be 1.0, the
				                              simulation that will be run will correspond to a Normal-Mode Path Integral
							      Molecular Dynamics (NMPIMD) simulation. If the adiabaticity parameter is <1.0
							      then a Centroid Molecular Dynamics (CMD) simulation is run. The parameter
							      must always be greater than 0.0. see ref ()))(& for details

				       Potential_Choice: an integer option between 1-4 which species the type of central potential to
				                         which the particle is attracted. The following potential choices have been 
							 specified:
							      1- Harmonic potential (i.e., V(x)=0.5*k*x*x)
							      2- Morse potential (NOT YET PROGRAMMED IN!!)
							      3- Anharmonic potential, option one (V(x)=(0.5*x*x)+(0.1*x*x*x)+(0.01*x*x*x*x))
							      4- Anharmonic potential, option two (V(x)=0.25*x*x*x*x)

							            **room for improvement, read in an arbitrary potential and base dynamics on that?

				       Detach_Centroid_Thermostat: 0 or 1 (bool value), used to detach the thermostat from the centroid in the final
				                                   simulation if the value is 1; if zero, then the centroid is kept in the NVT
								   ensemble.

				       Restart_Simulation: 0 or 1 (bool value), used to restart a simulation in the same ensemble it finished in, if
				                           upon review of the simulation trajectory it appears that the sampling was insufficient.

Create-Position-Autocorrelation: 

     Usage: (from command-line, example)

          "./Create-Position-Autocorrelation 6 11 1000 < STDIN > STDOUT"

	       Uses the trajectory generated from a simulation to compute an autocorrelation function of an observable-in this case, the observable
	       is x^2(t), and the corresponding autocorrelation function generated is "C(t) = <x^2(0)*x^2(t)> ". The trajectory is redirected to
	       STDIN of the program, and the output can either be read from the screen or (more usefully) redirected to a file.

     Input Parameters from command-line: Column_containing_positions, total_number_of_columns, number_of_data_points

          Parameter Description: Column_containing_positions: integer, specifies which column number in the input fil corresponds to the column
	                                                      displaying the observable of interest at that point in the trajectory (position, in 
					                      this case)

			         Total_number_of_columns: integer, total number of columns contained in the input file
				                                ***this could be done dynamically, note for improvement

				 Number_of_data_points: integer, total number of data points from which an autocorrelation function will be computed.
				                        This is typically equal to the number of snapshots in the trajectory file-the STDOUT file
							specified during the execution of 1D_oscillator.

Position-Probability-Distribution:

     Usage (from command-line, example)

          "./Position-Probability-Distribution < STDIN > STDOUT "

	       Uses the trajectory file output from the simulation to create a probability distribution chart of the particle during the course of
	       simulation. The column containing the particle position (or the centroid position, in this case) and the total number of columns in
	       the STDIN are hardcoded into the program (create_pos_histogram.cpp). 
	            **Input parameters can be improved

Background:

The purpose of this program is to model a 1-dimensional particle attracted to a central potential.

The simulation utilizes the theorems of Centroid Molecular Dynamics (CMD) to capture Quantum Dyanamical effects on the system.
These effects are captured by the calculation of autocorrelation functions and probability densities.

CMD is an extension of the Path Integral formulation of Quantum Statistical Mechanics, originally proposed by Feynman (ca. ??)

The Suzuki integration scheme for the propagation of Nose-Hoover-Chain thermostats was originally provided by (ref 1, ...) and then
implemented by Paesani and Voth (ca. 2005). Previous editions of the software utilized integration schemes derived from a simpler 
velocity-verlet framework probided by Voth and Jang (1997), but it was found that the error int simpler schemes grew too much with
larger values of P and M (see Tuckerman, Statistical Mechanics: Theory and Simulation, for the full decomposition of the Liouville
Operator with respect to the aformentioned integrator).
