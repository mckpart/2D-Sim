#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include <yaml-cpp/yaml.h>

#include "Particle.h"
#include "kiss.h"

class Boundary{

private:
   double boxLength = 0; 

public:

   void initializeBoundary(std::string yamlFile); 	

   void initialPosition(std::vector<Particle>* particles, int n_particles, KISSRNG randVal); 
   void periodicBoundary(std::vector<Particle>* particles, int index);
   bool rigidBoundary(std::vector<Particle>* particles, int index);
   double externalWell(std::vector<Particle>* particles, int index);  
};
#endif