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
   
   if(boxLength == 0){
      boxLength = sigma * sqrt(n_particles/redDensity); 
   }
   else if(sigma == 0){  // add a check such that a bunch of zeros can't be added
      sigma = boxLength * sqrt(redDensity/n_particles); 
   }

   std::cout << "reduced density: " << redDensity << " redTemp: " << redTemp 
	     << " sigma: " << sigma << " boxLength: " << boxLength << std::endl;


   seed        = node["seed"].as<long>(); 
   n_updates   = node["numberUpdates"].as<int>(); 
 
//   beta        = node["beta"].as<double>();

   init_type   = node["initializationType"].as<int>(); 

   rigidBC     = node["rigidBoundary"].as<bool>(); 
   periodicBC  = node["periodicBoundary"].as<bool>(); 
   extWell     = node["externalWell"].as<bool>(); 

   hardDisk    = node["hardDisks"].as<bool>(); 
   lenJones    = node["lennardJones"].as<bool>(); 	
   WCA         = node["WCA"].as<bool>(); 
//   c_linkers   = node["crosslinkers"].as<bool>(); 
}


///// GETTERS ////////////////

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

