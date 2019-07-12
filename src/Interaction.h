#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>
#include <yaml-cpp/yaml.h>

#include "Particle.h"
#include "kiss.h"

class Interaction{ 
	
private:
   double LJ_wellDepth = 0; // for LJ Potential
   double sigma = 0; 

   double restLength = 0;  // for crosslinking 
   double sprConstant = 0; // potential  	
  
   double truncDist = 0;
   double beta = 0; 
public: 

   void initializeInteraction(std::string yamlFile); 

   double lennardJones(std::vector<Particle>* particles, int index, int n_particles); 
   double WCApotential(std::vector<Particle>* particles, int index, int n_particles);
   double crosslinkers(std::vector<Particle>* particles, int index, int n_particles); 
   bool hardDisks     (std::vector<Particle>* particles, int index, int n_particles);


};
#endif
