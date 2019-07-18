#include <iostream>
#include "Parameters.h"


void Parameters::initializeParameters(std::string yamlFile){

   YAML::Node node = YAML::LoadFile(yamlFile); 

   seed        = node["seed"].as<long>(); 
   n_particles = node["totalParticles"].as<int>();
   n_updates   = node["numberUpdates"].as<int>(); 
 
   beta        = node["beta"].as<double>();

   init_type   = node["initializationType"].as<int>(); 

   rigidBC     = node["rigidBoundary"].as<bool>(); 
   periodicBC  = node["periodicBoundary"].as<bool>(); 
   extWell     = node["externalWell"].as<bool>(); 

   hardDisk    = node["hardDisks"].as<bool>(); 
   lenJones    = node["lennardJones"].as<bool>(); 	
   WCA         = node["WCA"].as<bool>(); 
   c_linkers   = node["crosslinkers"].as<bool>(); 
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
bool Parameters::getCrosslinkers(){
   return c_linkers; 
}

