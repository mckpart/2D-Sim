#include <iostream>
#include <yaml-cpp/yaml.h>
#include <string>
#include <fstream>

#include "Particle.h"
#include "Parameters.h"
#include "Interaction.h"
#include "Boundary.h"

Parameters setSystemParams(Parameters param, std::string yamlFile); 
void setParticleParams(Particle* particles,std::string yamlFile);

void simulation(Particle* particles, Parameters param){
	
	Boundary bound;
	Interaction interact; 
	Particle prt; 

	int n_particles = 0; 
	int n_updates 	= 0;

	double n_rejects = 0;  
	double perc_rej  = 0; 

	double x_trial = 0; 
	double y_trial = 0; 

	bool accept = 0; 

	std::ofstream pos_file; 
	pos_file.open("positions.txt"); 

	KISSRNG randVal; 
	randVal.InitCold(param.getSeed()); 

	n_particles = param.getNumParticles(); 
	n_updates 	= param.getUpdates(); 

	interact.initialPosition(particles,n_particles,randVal); 

	for(int n = 0; n < n_updates; n++){
		for(int k = 0; k < n_particles; k++){
			prt = particles[k]; 

			x_trial = prt.x_trial(randVal.RandomUniformDbl());  
			y_trial = prt.y_trial(randVal.RandomUniformDbl()); 

			if(param.getHardDisk() == 1){
				accept = interact.hardDisks(particles,k,n_particles,x_trial,y_trial);
				// std::cout << "hard disk accept: " << accept << std::endl; 
			}
			if(param.getRigidBC() == 1 && accept == 1){
				accept = bound.rigidBoundary(particles,k,n_particles,x_trial,y_trial); 
				std::cout << "rigid boundary accept: " << accept << std::endl; 
			}

			if(accept == 1){
				prt.setX_Position(x_trial);
				prt.setY_Position(y_trial); 
				particles[k] = prt; 

				pos_file << prt.getX_Position() << " "; 
				pos_file << prt.getY_Position() << " "; 
			}
			else{
				n_rejects++; 
				k = k - 1;   					// reject trial move
			}

		}

		pos_file << std::endl; 					// starts new line in position file
	}

	perc_rej = n_rejects / (n_rejects + (n_updates * n_particles)) * 100.0; 
	std::cout << perc_rej << "% of the moves were rejected." << std::endl; 
}


int main(int num, char *input[]){

	if(num > 1){	

		int n_particles = 0; 
		std::string yamlFile;

		yamlFile = input[1];

		Parameters param;
		param = setSystemParams(param,yamlFile); 

		Particle particles[param.getNumParticles()]; 
		setParticleParams(particles,yamlFile); 

		simulation(particles, param); 
	}
	else{
		std::cout << "ERROR: NO .YAML FILE" << std::endl; 
	}
}

Parameters setSystemParams(Parameters param, std::string yamlFile){		// change this to a pointer

	YAML::Node node = YAML::LoadFile(yamlFile); 

	param.setSeed(			node["seed"].as<long>()); 
	param.setNumParticles(	node["totalParticles"].as<int>());
	param.setUpdates(		node["numberUpdates"].as<int>()); 

	param.setRigidBC(		node["rigidBoundary"].as<bool>()); 
	param.setHardDisk(		node["hardDisks"].as<bool>()); 

	return param; 
}

void setParticleParams(Particle* particles,std::string yamlFile){

	int num_part1 = 0; 
	int num_part2 = 0; 

	double radius_1 = 0; 
	double radius_2 = 0; 

	YAML::Node node = YAML::LoadFile(yamlFile);

	num_part1 = node["type1_Particles"].as<int>(); 
	num_part2 = node["type2_Particles"].as<int>(); 

	radius_1 = node["particleRadius_1"].as<double>(); 
	radius_2 = node["particleRadius_2"].as<double>(); 

	Particle prt; 
	prt.setStepWeight(node["weight"].as<double>());

	prt.setRadius(radius_1);
	prt.setType(1); 

	for(int k = 0; k < num_part1; k++){
		prt.setIdentifier(k); 
		particles[k] = prt; 
	}

	prt.setRadius(radius_2); 
	prt.setType(2); 

	for(int k = 0; k < num_part2; k++ ){
		prt.setIdentifier(k + num_part1); 
		particles[k + num_part1] = prt; 
	}

}