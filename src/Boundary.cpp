#include <iostream>
#include <cmath>

#include "Boundary.h"

bool Boundary::rigidBoundary(Particle* particles, int index, int n_particles, double x_temp, double y_temp){

	Particle current_prt; 
	Particle compare_prt; 	

	double x_wall = 0; 
	double y_wall = 0; 	

	double rad_temp   = 0; 
	double num 	  	  = 0; 
	double wall_bound = 0; 	

	bool accept = 0;	

	current_prt = particles[index];
	rad_temp = current_prt.getRadius(); 
	wall_bound = 1; 	

	accept = 1; 

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

	if(fabs(x_temp - x_wall) < rad_temp){
		accept = 0;
	}
	else if(fabs(y_temp - y_wall) < rad_temp){
		accept = 0; 
	}

	return accept; 
}