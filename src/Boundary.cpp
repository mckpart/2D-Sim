#include <iostream>
#include <cmath>

#include "Boundary.h"

double distance(double x1,double x2){
	return fabs(x2 - x1); 
}

bool Boundary::rigidBoundary(Particle* particles, int index, int n_particles){

	Particle current_prt; 

	double x_wall = 0; 
	double y_wall = 0; 	

	double x_temp	= 0; 
	double y_temp	= 0; 
	double rad_temp	= 0; 

	double num 	  	  = 0; 
	double wall_bound = 0; 	

	bool accept = 0;	

	current_prt = particles[index];

	x_temp = current_prt.getX_TrialPos(); 
	y_temp = current_prt.getY_TrialPos(); 
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

void Boundary::periodicBoundary(Particle* particles, int index, int n_particles){

	Particle current_prt; 

	double x_wall = 0; 
	double y_wall = 0; 	

	double x_curr 	= 0 ; 
	double y_curr 	= 0 ;
	double x_temp	= 0; 
	double y_temp	= 0; 
	double rad_temp = 0; 

	double num 	  	  = 0; 
	double wall_bound = 0; 	

	double dist_xtravel = 0; 
	double dist_ytravel = 0; 
	double dist_xwall 	= 0; 
	double dist_ywall 	= 0 ;

	int flag = 0; 

	// bool accept = 0;	

	current_prt = particles[index];

	x_curr = current_prt.getX_Position(); 
	y_curr = current_prt.getY_Position(); 
	x_temp = current_prt.getX_TrialPos(); 
	y_temp = current_prt.getY_TrialPos(); 	
	rad_temp = current_prt.getRadius(); 

	wall_bound = 1; 					// will be added as private variable

	dist_xtravel = distance(x_curr,x_temp); 
	dist_ytravel = distance(y_curr,y_temp); 

	if(x_curr * -1 < 0){
		x_wall = wall_bound - rad_temp; 

		if(x_temp - x_wall > 0){

			dist_xwall = distance(x_temp,x_wall); 
			x_temp = -1 * x_wall + dist_xwall; 
		}
	}
	else{
		x_wall = -1 * (wall_bound - rad_temp);

		if(x_temp - x_wall < 0){

			dist_xwall  = distance(x_temp,x_wall) ; 
			x_temp = -1 * x_wall - dist_xwall; 
		} 
	}


	if(y_curr * -1 < 0){
		y_wall = wall_bound - rad_temp; 

		if(y_temp - y_wall > 0){

			dist_ywall = distance(y_temp,y_wall); 
			y_temp = -1 * y_wall + dist_ywall; 
		}
	}
	else{
		y_wall = -1 * (wall_bound - rad_temp);

		if(y_temp - y_wall < 0){

			dist_ywall  = distance(y_temp,y_wall) ; 
			y_temp = -1 * y_wall - dist_ywall; 
		} 
	}

	current_prt.setX_TrialPos(x_temp); 	// puts the updated particle with updated 
	current_prt.setY_TrialPos(y_temp); 	// trial positions back into the array 

	particles[index] = current_prt; 

}
