#include <iostream>
#include "Parameters.h"


void Parameters::initializeParameters(std::string yamlFile){

   YAML::Node node = YAML::LoadFile(yamlFile); 

   /* SET THE REDUCED PARAMETERS OF THE SYTEM */

   redDensity  = node["reducedDens"].as<double>();
   redTemp     = node["reducedTemp"].as<double>(); 
   sigma       = node["sigma"].as<double>(); 
   boxLength   = node["boxLength"].as<double>();
   n_particles = node["totalParticles"].as<int>();

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

   std::cout << "reduced density: " << redDensity 
	     << "\n redTemp: " << redTemp 
	     << "\n sigma: " << sigma 
	     << "\n boxLength: " << boxLength << std::endl;

   seed        = node["seed"].as<long>(); 
   n_updates   = node["numberUpdates"].as<int>(); 
 
   init_type   = node["initializationType"].as<int>(); 

   rigidBC     = node["rigidBoundary"].as<bool>(); 
   periodicBC  = node["periodicBoundary"].as<bool>(); 
   extWell     = node["externalWell"].as<bool>(); 

   hardDisk    = node["hardDisks"].as<bool>(); 
   lenJones    = node["lennardJones"].as<bool>(); 	
   WCA         = node["WCA"].as<bool>(); 
}

///// GETTERS ////////////////

int Parameters::getData_interval(){
   return d_interval; 
}
int Parameters::getEq_sweep(){
   return eq_sweep; 
}
int Parameters::getUpdates(){
   return n_updates; 
} 	
int Parameters::getNumParticles(){
   return n_particles; 
}
long Parameters::getSeed(){
   return seed; 
}

double Parameters::getWeight(){
   return weight; 
}
double Parameters::getBeta(){
   return beta; 
}

int Parameters::getInit_Type(){
   return init_type; 
}

bool Parameters::getRigidBC(){
   return rigidBC; 
}
bool Parameters::getPeriodicBC(){
   return periodicBC; 
}
bool Parameters::getExtWell(){
   return extWell; 
}

bool Parameters::getHardDisk(){
   return hardDisk; 
}
bool Parameters::getLenJones(){
   return lenJones; 
}
bool Parameters::getWCA(){
   return WCA; 
}

double Parameters::getLJ_const_1(){
   return LJ_const_1; 
}
double Parameters::getLJ_const_2(){
   return LJ_const_2; 
}

double Parameters::getSigma(){
   return sigma; 
}
double Parameters::getRedDens(){
   return redDensity; 
}
double Parameters::getRedTemp(){
   return redTemp; 
}
double Parameters::getBoxLength(){
   return boxLength; 
}

