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
   std::vector<double> par_num_density; 
   std::vector<double> antp_num_density; 
   
//   std::vector<double> x_num_density; // maybe make this into a number density as well? 
//   std::vector<double> y_num_density; 

//   std::vector<double> y_relate_x;
//   std::vector<double> x_relate_y; 

   std::vector<std::vector<double>> xy_num_density; // 2D arrays used for calculating 
   std::vector<std::vector<double>> par_xy_density; // the pair correlation functions
   std::vector<std::vector<double>> antp_xy_density; 
   
   double delta_r = 0;
   double cell_L = 0;  
   double sigma = 0; 
   double truncDist = 0; 
   double truncShift = 0; 

   double interact_type = 0; 

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

   void updateNumDensity(double r, int ID);
//   void update_xy_corr(double x, double y, int ID); 
   void calc_xy_dens(double x, double y, int ID); 

   void calcPeriodicProp(std::vector<Particle>* particles,std::ofstream* r_dist_file); 
   
   void calcEnergy(double r, double c);
   void calcVirial(double r, double c); 
   
   double radDistance(double x1, double x2, double y1, double y2); 

   double calcPressure();    
   double calcAvgEnergy(); 
   
   void writeProperties(); 
}; 
#endif
