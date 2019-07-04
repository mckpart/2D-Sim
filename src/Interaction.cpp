#include <iostream>
#include "Interaction.h"
#include <cmath>
#include <array>


double distance(double x1,double x2,double y1, double y2){
	return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2)); 	
}


bool Interaction::rigidCollisions(Particle* particles, int index, int n_particles, double x_temp, double y_temp){ // consider making x,y trials part of particle class
	
	Particle current_prt; 
	Particle compare_prt; 						

	double x_comp = 0; 
	double y_comp = 0; 
	double x_wall = 0; 
	double y_wall = 0; 							// the locations of the nearest 'wall'

	double rad_temp   = 0; 
	double rad_comp   = 0; 
	double num 	  	  = 0; 
	double wall_bound = 0; 						

	bool accept = 0; 


	current_prt = particles[index];
	rad_temp = current_prt.getRadius(); 
	wall_bound = 1; 	

	accept = 1; 

////// CHECK FOR PARTICLE-WALL COLLISION //////////////

	if(x_temp * -1 < 0){		
		x_wall = wall_bound; 		// x_temp > 0, closest x_wall > 0
		if(x_temp - x_wall > 0){
			accept = 0; 
		}
	}
	else{
		x_wall = -1 * wall_bound; 
		if(x_temp - x_wall < 0){
			accept = 0; 
		}
	}
	
	if(y_temp * -1 < 0){			
		y_wall = wall_bound; 		// y_temp > 0, closest y_wall > 0
		if(y_temp - y_wall > 0){
			accept = 0; 
		}
	}
	else{
		y_wall = -1 * wall_bound; 
		if(y_temp - y_wall < 0){
			accept = 0;
		}
	}

//////// CHECK FOR PARTICLE-PARTICLE COLLISION //////////

	if(accept != 0){
		if(abs(x_temp - x_wall) < rad_temp){
			accept = 0;
		}
		else if(abs(y_temp - y_wall) < 2 * rad_temp){
			accept = 0; 
		}
		else{
			for(int k = 0; k < n_particles; k++){
				if(k != index){
					compare_prt = particles[k]; 		// compare all particles positions
														// to current particle's position
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
		}
	}
		
	return accept; 										// returns 1 if trial move is accepted 
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