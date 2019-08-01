#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <vector>
#include <cmath>
#include <fstream>

#include "Parameters.h"
#include "Particle.h"

class Properties{

private:
   double f_energy = 0; 
   double f_r = 0; 

   std::vector<double> sum_Fdot_r; 
   std::vector<double> sum_energy; 
   std::vector<double> num_density; 

   double delta_r = 0; 
   double sigma = 0; 
   double truncDist = 0; 
   double truncShift = 0; 

   bool LJ = 0;

   double LJ_par = 0; 
   double LJ_antipar = 0; 

   double boxLength = 0; 
   int n_particles = 0; 

   double redDens = 0; 
   double redTemp = 0;  
public: 

   void initializeProperties(Parameters* p); 

   void populateCellArray(double x,double y, std::vector<std::vector<double>>* cellPositions); 
   double lenJonesEnergy(double r, double c);  
   double lenJonesForce(double r, double c);   

   void updateNumDensity(double r);

   void calcPeriodicProp(std::vector<Particle>* particles,std::ofstream* r_dist_file); 
   
   void calcEnergy(double r, double c);
   void calcVirial(double r, double c); 
   
   double radDistance(double x1, double x2, double y1, double y2); 

   double calcPressure();    
   double calcAvgEnergy(); 
   
   void writeProperties(); 
}; 
#endif
