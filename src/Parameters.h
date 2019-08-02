#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <yaml-cpp/yaml.h>
#include <cmath>

class Parameters{

private:

   double redDensity = 0; 
   double redTemp = 0; 
   double sigma = 0;
   double boxLength = 0; 

   double LJ_const_1 = 0; 
   double LJ_const_2 = 0; 

   int eq_sweep = 0; 
   int d_interval = 0; 

   int n_updates = 0; 
   int n_particles = 0; 
   long seed = 0; 

   double weight = 0; 
   double beta = 0; 

   int init_type = 0; 
   int interact_type = 0; 

   bool rigidBC = 0; 
   bool periodicBC = 0; 
   bool extWell = 0; 
 
   bool hardDisk = 0;  
   bool lenJones = 0;  			
   bool WCA = 0; 
   bool c_linkers = 0; 

public: 
		
   void initializeParameters(std::string yamlFile); // may need to add a default constructor

   int getUpdates(); 	
   int getNumParticles(); 
   long getSeed(); 

   double getWeight(); 
   double getBeta(); 

   int getInit_Type(); 
   int getInteract_Type(); 

   int getEq_sweep(); 
   int getData_interval(); 

   bool getRigidBC(); 		
   bool getPeriodicBC(); 		
   bool getExtWell(); 
//
//   bool getHardDisk(); 
//   bool getLenJones(); 
//   bool getWCA(); 

   double getLJ_const_1(); 
   double getLJ_const_2(); 

   double getRedDens(); 
   double getRedTemp();
   double getSigma();
   double getBoxLength(); 

   // NOTE: THE SETTERS ARE NOT NECESSARY
   // GIVEN THAT THESE PARAMETERS ARE CONSTANT
   // THROUGHOUT THE SIMULATION

}; 
#endif 
