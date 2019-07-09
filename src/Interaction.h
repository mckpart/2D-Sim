#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>
#include <yaml-cpp/yaml.h>

#include "Particle.h"
#include "kiss.h"

class Interaction{ 
	
private:

   double LJ_wellDepth = 0; 	

public: 

   void initializeInteraction(std::string yamlFile); 

   double lennardJones(std::vector<Particle>* particles, int index, int n_particles); 
   double WCApotential(std::vector<Particle>* particles, int index, int n_particles);
   bool hardDisks     (std::vector<Particle>* particles, int index, int n_particles); 
   // void initialPosition	(std::vector<Particle>* particles, int n_particles, KISSRNG randVal);

};
#endif