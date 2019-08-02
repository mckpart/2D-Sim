#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>

#include "Parameters.h"
#include "Particle.h"
#include "kiss.h"

class Interaction{ 
	
private:
   double sigma = 0; 
   double redDens = 0; 
   double LJ_par = 0; 
   double LJ_antipar = 0; 

   int interact_type = 0; 

   double tail_corr = 0; 
   double truncDist = 0;
   double truncShift = 0; 

   double boxLength = 0; 
   int n_particles = 0;  
public: 

   void initializeInteraction(Parameters* p); 
   void populateCellArray(double x, double y, 
		          std::vector<std::vector<double>>* cellPositions); 

   double distance(double x1, double x2, double y1, double y2);
   
   double lenjones_energy(double r, double c); 
   double WCA_energy     (double r, double c); 

   double nonPeriodicInteraction(std::vector<Particle>* particles, int index); 
   double periodicInteraction   (std::vector<Particle>* particles, int index); 
   bool hardDisks               (std::vector<Particle>* particles, int index);


};
#endif
