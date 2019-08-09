#include <iostream>
#include "Parameters.h"


void Parameters::initializeParameters(std::string yamlFile){

   YAML::Node node = YAML::LoadFile(yamlFile); 

   /* SET THE REDUCED PARAMETERS OF THE SYTEM */

   redDensity  = node["reducedDens"].as<double>();
   redTemp     = node["reducedTemp"].as<double>(); 
   sigma       = node["sigma"].as<double>(); 
   n_particles = node["totalParticles"].as<int>(); 
   sigma_par = sigma; // BAD, COME BACK TO FIX
   k_spring    = node["springConstant"].as<double>(); 
//   rest_L      = node["rest_length"].as<double>(); 
   rest_L = sigma; // this is temporary
   LJ_const_1  = node["LJ_constant_1"].as<double>(); 
   LJ_const_2  = node["LJ_constant_2"].as<double>(); 

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
//   sigma_par = sqrt(2)*sigma; 
   std::cout << "reduced density: " << redDensity 
	     << "\nredTemp: " << redTemp 
	     << "\nsigma: " << sigma 
	     << "\nboxLength: " << boxLength << std::endl;

   seed        = node["seed"].as<long>(); 
   n_updates   = node["numberUpdates"].as<int>(); 
 
   init_type     = node["initializationType"].as<int>(); 
   interact_type = node["interactionType"].as<int>(); 
   bound_type    = node["boundaryType"].as<int>(); 
}

///// GETTERS ////////////////

int Parameters::getData_interval(){return d_interval;}
int Parameters::getEq_sweep(){     return eq_sweep;}
int Parameters::getUpdates(){      return n_updates;} 	
int Parameters::getNumParticles(){ return n_particles;}
long Parameters::getSeed(){        return seed;}

double Parameters::getWeight(){    return weight;}
double Parameters::getBeta(){      return beta;} // come back to remove

int Parameters::getInit_Type(){    return init_type;}
int Parameters::getInteract_Type(){return interact_type;}
int Parameters::getBound_Type(){   return bound_type;}

double Parameters::getLJ_const_1(){return LJ_const_1;}
double Parameters::getLJ_const_2(){return LJ_const_2;}

double Parameters::getSigma(){     return sigma;}
double Parameters::getSigmaPar(){  return sigma_par;}
double Parameters::getRedDens(){   return redDensity;}
double Parameters::getRedTemp(){   return redTemp;}
double Parameters::getBoxLength(){ return boxLength;}

double Parameters::getRestLength(){return rest_L;}
double Parameters::getSprConst()  {return k_spring;}
