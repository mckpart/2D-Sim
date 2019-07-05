#include <iostream>
#include "Parameters.h"

// Parameters::Parameters(){
// 	n_updates = 0; 
// 	n_particles = 0; 
// 	seed = 0; 
// 	weight = 0; 
// }


Parameters::Parameters(std::string yamlFile){

	YAML::Node node = YAML::LoadFile(yamlFile); 

	seed 		= node["seed"].as<long>(); 
	n_particles = node["totalParticles"].as<int>();
	n_updates 	= node["numberUpdates"].as<int>(); 

	rigidBC 	= node["rigidBoundary"].as<bool>(); 
	periodicBC	= node["periodicBoundary"].as<bool>(); 

	hardDisk 	= node["hardDisks"].as<bool>(); 
	lenJones 	= node["lennardJones"].as<bool>(); 	
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

bool Parameters::getRigidBC(){
	return rigidBC; 
}
bool Parameters::getPeriodicBC(){
	return periodicBC; 
}

bool Parameters::getHardDisk(){
	return hardDisk; 
}
bool Parameters::getLenJones(){
	return lenJones; 
}

///////// SETTERS /////////////

void Parameters::setUpdates(int u){
	n_updates = u; 
}
void Parameters::setNumParticles(int num){
	n_particles = num; 
}
void Parameters::setSeed(long s){
	seed = s; 
}

void Parameters::setWeight(double w){
	weight = w; 
}

void Parameters::setRigidBC(bool r){
	rigidBC = r; 
}
void Parameters::setPeriodicBC(bool p){
	periodicBC = p; 
}

void Parameters::setHardDisk(bool h){
	hardDisk = h; 
}
void Parameters::setLenJones(bool l){
	lenJones = l; 
}