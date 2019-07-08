#ifndef SIMULATION_H
#define SIMULATION_H

#include <yaml-cpp/yaml.h>
#include <string>
#include <fstream>
#include <vector>

#include "Particle.h"
#include "Parameters.h"
#include "Interaction.h"
#include "Boundary.h"

class Simulation{

private: 

	std::string yamlFile;

	Parameters param; 
	Interaction interact;
	Boundary bound;

	std::vector<Particle> particles; 	  

	int n_particles = 0; 

public: 

	Simulation(std::string yf); 

	void runSimulation(); 
	void setParticleParams(); 

}; 
#endif