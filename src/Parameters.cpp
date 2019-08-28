#include <iostream>
#include "Parameters.h"


void Parameters::initializeParameters(std::string yamlFile){

   YAML::Node node = YAML::LoadFile(yamlFile); 

   /* SET THE REDUCED PARAMETERS OF THE SYTEM */

   redDensity  = node["reducedDens"].as<double>();
   redTemp     = node["reducedTemp"].as<double>(); 
   sigma       = node["sigma"].as<double>(); 
   n_particles = node["totalParticles"].as<int>(); 
   k_spring    = node["springConstant"].as<double>(); 
  
   rest_L      = node["rest_length"].as<double>(); 
//   rest_L = 2.6; // this is = 65nm/25nm = c-c/diam
   
//   LJ_const_1  = node["LJ_constant_1"].as<double>(); 
//   LJ_const_2  = node["LJ_constant_2"].as<double>(); 
   
   a_ref       = node["reference_affinity"].as<double>(); 
   a_mult      = node["affinity_multiple"].as<double>(); 

   eq_sweep    = node["equilibriate_sweep"].as<int>(); 
   d_interval  = node["data_collect_interval"].as<int>(); 

   if(boxLength == 0){
      boxLength = sigma * sqrt(n_particles/redDensity); 
   }
   else if(sigma == 0){  // add a check such that a bunch of zeros can't be added
      sigma = boxLength * sqrt(redDensity/n_particles); 
   }
   else if(redDensity == 0){
      redDensity = n_particles * pow(sigma/boxLength,2);
   }
   std::cout << "reduced density: " << redDensity 
	     << "\nreduced temp: " << redTemp 
	     << "\nsigma: " << sigma 
	     << "\nbox length: " << boxLength 
	     << "\nrest length: " << rest_L << std::endl;

   seed        = node["seed"].as<long>(); 
   n_updates   = node["numberUpdates"].as<int>(); 
 
   init_type     = node["initializationType"].as<int>(); 
   interact_type = node["interactionType"].as<int>(); 
   bound_type    = node["boundaryType"].as<int>(); 
}

///// GETTERS ////////////////

int Parameters::getData_interval(){return d_interval;}
int Parameters::getEq_sweep()     {return eq_sweep;}
int Parameters::getUpdates()      {return n_updates;} 	
int Parameters::getNumParticles() {return n_particles;}
long Parameters::getSeed()        {return seed;}

int Parameters::getInit_Type()    {return init_type;}
int Parameters::getInteract_Type(){return interact_type;}
int Parameters::getBound_Type()   {return bound_type;}

//double Parameters::getLJ_const_1(){return LJ_const_1;}
//double Parameters::getLJ_const_2(){return LJ_const_2;}

double Parameters::getRefAffinity() {return a_ref;}
double Parameters::getAffinityMult(){return a_mult;}; 

double Parameters::getSigma()     {return sigma;}
//double Parameters::getSigmaPar(){  return sigma_par;}
double Parameters::getRedDens()   {return redDensity;}
double Parameters::getRedTemp()   {return redTemp;}
double Parameters::getBoxLength() {return boxLength;}

double Parameters::getRestLength(){return rest_L;}
double Parameters::getSprConst()  {return k_spring;}
