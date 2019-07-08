#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <yaml-cpp/yaml.h>

class Parameters{

	private:

		int n_updates   = 0; 
		int n_particles = 0; 
		long seed = 0; 

		double weight = 0; 

		bool rigidBC 	= 0; 
		bool periodicBC = 0; 

		bool hardDisk = 0;  
		bool lenJones = 0;  	
		bool WAC	  = 0; 

	public: 
		
		Parameters(std::string yamlFile); 

		int getUpdates(); 	
		int getNumParticles(); 
		long getSeed(); 

		double getWeight(); 

		bool getRigidBC(); 
		bool getPeriodicBC(); 

		bool getHardDisk(); 
		bool getLenJones(); 
		bool getWAC(); 

		// NOTE: THE SETTERS ARE NOT NECESSARY
		// GIVEN THAT THESE PARAMETERS ARE CONSTANT
		// THROUGHOUT THE SIMULATION

		// void setUpdates(int u);
		// void setNumParticles(int num); 
		// void setSeed(long s); 
		
		// void setWeight(double w); 

		// void setRigidBC(bool r); 
		// void setPeriodicBC(bool p); 

		// void setHardDisk(bool h); 
		// void setLenJones(bool l); 
		// void setWAC(bool w); 
}; 
#endif 