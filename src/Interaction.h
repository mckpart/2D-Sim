#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>

#include "Parameters.h"
#include "Particle.h"
#include "kiss.h"

class Interaction{ 
	
private:
   // double LJ_wellDepth = 0; // for LJ Potential
   double sigma = 0; 
   double redDens = 0; 
   double tail_corr = 0; 

//   double restLength = 0;  // for crosslinking 
  
//   double sprConstant = 0; // potential  	
   double truncDist = 0;
   double truncShift = 0; 
//   double beta = 0;

   double boxLength = 0; 
   int n_particles = 0;  
public: 

   void initializeInteraction(Parameters* p); 
   void populateCellArray(double x, double y, 
		   std::vector<std::vector<double>>* cellPositions); 
   
   double lennardJones(std::vector<Particle>* particles, int index); 
   double WCApotential(std::vector<Particle>* particles, int index);
   double crosslinkers(std::vector<Particle>* particles, int index); 
   bool hardDisks     (std::vector<Particle>* particles, int index);


};
#endif
