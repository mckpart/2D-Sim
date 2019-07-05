#ifndef INTERACTION_H
#define INTERACTION_H

#include <yaml-cpp/yaml.h>
#include "Particle.h"
#include "kiss.h"

class Interaction{ 
	
	private:

		double LJ_wellDepth = 0; 	

	public: 

		Interaction(std::string yamlFile); 

		double lennardJones		(Particle* particles, int index, int n_particles); 
		double WACpotential		(Particle* particles, int index, int n_particles);
		bool hardDisks			(Particle* particles, int index, int n_particles); 
		void initialPosition	(Particle* particles, int n_particles, KISSRNG randVal);

		////// GETTERS AND SETTERS //////

		double getWellDepth(); 

		void setWellDepth(double w); 	 
};
#endif