#include <iostream>
#include "Interaction.h"
#include <cmath>
#include <array>


double distance(double x1,double x2,double y1, double y2){
	return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2)); 	
}

double Interaction::WACpotential(Particle* particles, int index, int n_particles){

	Particle current_prt; 
	Particle compare_prt; 						

	double x_temp	  = 0; 
	double y_temp 	  = 0; 
	double x_curr	  = 0; 
	double y_curr	  = 0; 
	double x_comp	  = 0; 
	double y_comp	  = 0; 

	double sigma = 0; 

	double delta_energy = 0; 
	double energy_curr   = 0; 
	double energy_temp   = 0; 

	double rad_curr   = 0; 
	double rad_comp   = 0; 
	double num 	  	  = 0; 

	double dist_curr  = 0; 
	double dist_temp  = 0; 	

	current_prt = particles[index];

	x_temp = current_prt.getX_TrialPos(); 
	y_temp = current_prt.getY_TrialPos(); 

	x_curr   = current_prt.getX_Position(); 
	y_curr   = current_prt.getY_Position(); 
	rad_curr = current_prt.getRadius(); 

	for(int k = 0; k < n_particles; k++){

		compare_prt = particles[k]; 

		if(current_prt.getIdentifier() != compare_prt.getIdentifier()){ 

			x_comp   = compare_prt.getX_Position(); 
			y_comp   = compare_prt.getY_Position(); 
			rad_comp = compare_prt.getRadius(); 

			dist_curr = distance(x_curr,x_comp,y_curr,y_comp); 
			dist_temp = distance(x_temp,x_comp,y_temp,y_comp); 

			sigma = (rad_curr + rad_comp); 				

			if(dist_curr > pow(2,1/6) * sigma){
				energy_curr = 0; 
			}	
			else{
				energy_curr = 4 * LJ_wellDepth * (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6) + .25); 
			}

			if(dist_temp > pow(2,1/6) * sigma){
				energy_temp = 0; 
			}
			else{
				energy_temp = 4 * LJ_wellDepth * (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6) + .25); 								
			}


			delta_energy = delta_energy + (energy_temp - energy_curr); 
		}
	}

	return delta_energy; 	
}

double Interaction::lennardJones(Particle* particles, int index, int n_particles){

	Particle current_prt; 
	Particle compare_prt; 						

	double x_temp	  = 0; 
	double y_temp 	  = 0; 
	double x_curr	  = 0; 
	double y_curr	  = 0; 
	double x_comp	  = 0; 
	double y_comp	  = 0; 

	double sigma = 0; 

	double delta_energy = 0; 
	double energy_curr   = 0; 
	double energy_temp   = 0; 

	double rad_curr   = 0; 
	double rad_comp   = 0; 
	double num 	  	  = 0; 

	double dist_curr  = 0; 
	double dist_temp  = 0; 	

	current_prt = particles[index];

	x_temp = current_prt.getX_TrialPos(); 
	y_temp = current_prt.getY_TrialPos(); 

	x_curr   = current_prt.getX_Position(); 
	y_curr   = current_prt.getY_Position(); 
	rad_curr = current_prt.getRadius(); 

	for(int k = 0; k < n_particles; k++){

		compare_prt = particles[k]; 

		if(current_prt.getIdentifier() != compare_prt.getIdentifier()){ 

			x_comp   = compare_prt.getX_Position(); 
			y_comp   = compare_prt.getY_Position(); 
			rad_comp = compare_prt.getRadius(); 

			dist_curr = distance(x_curr,x_comp,y_curr,y_comp); 
			dist_temp = distance(x_temp,x_comp,y_temp,y_comp); 

			sigma = (rad_curr + rad_comp); 

			if(current_prt.getType() == compare_prt.getType()){ // interaction betweeen like particles

				energy_curr = 30 * LJ_wellDepth * (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); 
				energy_temp = 30 * LJ_wellDepth * (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 				
			} 
			else{												// interaction between unlike particles

				energy_curr = 4 * LJ_wellDepth * (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); 
				energy_temp = 4 * LJ_wellDepth * (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 
			}	

			delta_energy = delta_energy + (energy_temp - energy_curr); 
		}
	}

	return delta_energy; 	
}

bool Interaction::hardDisks(Particle* particles, int index, int n_particles){ 
	
	Particle current_prt; 
	Particle compare_prt; 						

	double x_temp = 0; 
	double y_temp = 0; 
	double x_comp = 0; 
	double y_comp = 0; 

	double rad_temp   = 0; 
	double rad_comp   = 0; 
	double num 	  	  = 0; 

	bool accept = 0; 


	current_prt = particles[index];

	x_temp = current_prt.getX_TrialPos(); 
	y_temp = current_prt.getY_TrialPos(); 
	rad_temp = current_prt.getRadius(); 

	accept = 1; 

//////// CHECK FOR PARTICLE-PARTICLE COLLISION //////////
		
for(int k = 0; k < n_particles; k++){

	compare_prt = particles[k]; 		// compare all particles positions
										// to current particle's position	
	if(current_prt.getIdentifier() != compare_prt.getIdentifier()){
											
		x_comp = compare_prt.getX_Position();
		y_comp = compare_prt.getY_Position(); 
		rad_comp = compare_prt.getRadius(); 

		if(distance(x_temp,x_comp,y_temp,y_comp) 
					< rad_comp + rad_temp){
			accept = 0; 
			break; 
		}
	}	
}

	return accept; 							// returns 1 if trial move is accepted 
}


void Interaction::initialPosition(Particle* particles, int n_particles, 
								  KISSRNG randVal){
	Particle current_prt; 
	Particle compare_prt; 

	double x_temp = 0; 
	double y_temp = 0; 
	double x_comp = 0; 
	double y_comp = 0; 
	double x_wall = 0; 
	double y_wall = 0; 							// the locations of the nearest 'wall'

	double rad_temp = 0; 
	double rad_comp = 0; 
	double num 	  	= 0; 
	double wall_bound = 0; 						

	bool accept = 0; 

	for(int k = 0; k < n_particles; k++){

		current_prt = particles[k];
		rad_temp = current_prt.getRadius(); 
		wall_bound = 1;			 				// checks that distance is < 2 * radius

		x_temp = randVal.RandomUniformDbl(); 
		y_temp = randVal.RandomUniformDbl(); 

		////// DETERMINE THE 'QUADRANT' //////////////

		num = randVal.RandomUniformDbl(); 

		if(k % 2 == 0 && num <.5){				// the acceptance presented in 
			x_temp = -1 * x_temp; 				// 'check collision' is not implemented here
												// since the initial position cannot exceed 1
			x_wall = -1 * wall_bound;
			y_wall = wall_bound;  
		}
		else if(k % 2 == 0 && num >= .5){
			y_temp = -1 * y_temp; 

			x_wall = wall_bound; 
			y_wall = -1 * wall_bound; 
		}
		else if(k % 2 != 0 && num < .5){
			x_temp = -1 * x_temp; 
			y_temp = -1 * y_temp; 

			x_wall = -1 * wall_bound; 
			y_wall = -1 * wall_bound; 
		}
		else{
			x_wall = wall_bound; 
			y_wall = wall_bound; 
		}

		//////// DETERMINE THE ACCEPTANCE //////////

		accept = 1; 

		if(abs(x_temp - x_wall) < rad_temp){
			accept = 0;
		}
		else if(abs(y_temp - y_wall) < rad_temp){
			accept = 0; 
		}
		else if(k != 0){
			for(int n = 0; n < k; n++){
				compare_prt = particles[n]; 

				x_comp = compare_prt.getX_Position();
				y_comp = compare_prt.getY_Position(); 
				rad_comp = compare_prt.getRadius(); 

				if(distance(x_temp,x_comp,y_temp,y_comp) 
					        < rad_comp + rad_temp){
					accept = 0; 
					break; 
				}
			}
		}

		if(accept == 1){
			current_prt.setX_Position(x_temp);
			current_prt.setY_Position(y_temp); 

			particles[k] = current_prt; 
		}
		else{
			k = k - 1; 							// generate a new random position for the SAME particle
		}

	}
}

////// DEFAULT CONSTRUCTOR ////////////////

Interaction::Interaction(std::string yamlFile){

	YAML::Node node = YAML::LoadFile(yamlFile);

	LJ_wellDepth = node["wellDepth"].as<double>(); 	
}



/////// GETTERS AND SETTERS ///////////////

double Interaction::getWellDepth(){
	return LJ_wellDepth; 
}

void Interaction::setWellDepth(double w){
	LJ_wellDepth = w; 
}