#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <vector>
#include <cmath>
#include <yaml-cpp/yaml.h>
#include <fstream>

#include "Particle.h"

class Properties{

private:
   double f_energy = 0; 
   double Z_func = 0; 

   std::vector<double> sum_Fdot_r; 
   std::vector<double> sum_energy; 

   double sigma = 0; 
   double epsilon = 0; 
   double KbT = 0;
   double truncDist = 0; 

   double LJ = 0;

   double boxLength = 0; 
   int n_particles = 0; 

   double redDens = 0; 
   double redTemp = 0;  
public: 

   void initializeProperties(std::string yf); 

   double lenJonesEnergy(double r);  
   double lenJonesForce(double r);   
   
   void calcForces(std::vector<Particle>* particles); 
   void calcEnergy(double r_dist);
   
   double calcPressure();    

   void writeProperties(); 
}; 
#endif
