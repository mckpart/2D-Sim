#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>

#include "Parameters.h"
#include "Particle.h"
#include "kiss.h"

class Boundary{

private:
   double boxLength = 0; 
   double sigma = 0; 
   int interact_type = 0; 
   int n_particles = 0; 

//   bool LJ = 0; 
public:

   void initializeBoundary(Parameters *p); 	

   void initialPosition (std::vector<Particle>* particles, KISSRNG randVal); 
   int initialHexagonal (std::vector<Particle>* particles); // not random
   int initialSquare    (std::vector<Particle>* particles); 

   void periodicBoundary(std::vector<Particle>* particles, int index);
   bool rigidBoundary   (std::vector<Particle>* particles, int index);
   double externalWell  (std::vector<Particle>* particles, int index);  
};
#endif
