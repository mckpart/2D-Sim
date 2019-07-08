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

	double x_temp = 0; 
	double y_temp = 0; 
	double x_curr = 0; 
	double y_curr = 0; 
	double x_comp = 0; 
	double y_comp = 0; 

	double sigma = 0; 

	double delta_energy	= 0; 
	double energy_curr	= 0; 
	double energy_temp	= 0; 

	double rad_curr	= 0; 
	double rad_comp	= 0; 
	double num		= 0; 

	double dist_curr = 0; 
	double dist_temp = 0; 	

	current_prt = particles[index];			// assign current particle

	x_temp = current_prt.getX_TrialPos(); 	// assign the current and trial
	y_temp = current_prt.getY_TrialPos(); 	// positions and the radius of 
											// the current particle
	x_curr   = current_prt.getX_Position(); 
	y_curr   = current_prt.getY_Position(); 
	rad_curr = current_prt.getRadius(); 

	for(int k = 0; k < n_particles; k++){	// the O(N^2)... think about 
											// implementing neighbor lists here
		compare_prt = particles[k]; 		// assign comparison particle

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
				energy_curr = 4 * LJ_wellDepth * 
				(pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6) + .25); 
			}

			if(dist_temp > pow(2,1/6) * sigma){
				energy_temp = 0; 
			}
			else{
				energy_temp = 4 * LJ_wellDepth * 
				(pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6) + .25); 								
			}

			delta_energy = delta_energy + (energy_temp - energy_curr); 
		}
	}

	return delta_energy; 	
}

double Interaction::lennardJones(Particle* particles, int index, int n_particles){

	Particle current_prt; 
	Particle compare_prt; 						

	double x_temp = 0; 
	double y_temp = 0; 
	double x_curr = 0; 
	double y_curr = 0; 
	double x_comp = 0; 
	double y_comp = 0; 

	double sigma = 0; 

	double delta_energy = 0; 
	double energy_curr	= 0; 
	double energy_temp	= 0; 

	double rad_curr	= 0; 
	double rad_comp	= 0; 
	double num		= 0; 

	double dist_curr = 0; 
	double dist_temp = 0; 	

	current_prt = particles[index];			// assign current particle

	x_temp = current_prt.getX_TrialPos(); 	// assign the current and trial
	y_temp = current_prt.getY_TrialPos(); 	// positions and the radius of 
											// the current particle
	x_curr   = current_prt.getX_Position(); 
	y_curr   = current_prt.getY_Position(); 
	rad_curr = current_prt.getRadius(); 

	for(int k = 0; k < n_particles; k++){

		compare_prt = particles[k]; 		// assign the comparison particle


		/* 	- CHECK THAT THE COMPARISON PARTICLE IS NOT THE CURRENT PARTICLE
			- SET THE COMPARISON POSITION 
			- COMPUTE THE CURRENT DISTANCE BETWEEN THE TWO PARTICLES 
			- COMPUTE THE TRIAL DISTANCE BETWEEN THE TWO PARTICLES
			- IF THE PARTICLES ARE OF THE SAME TYPE, HAVE A DEEPER
			  LENNARD JONES POTENTIAL
			- IF THE PARTICLES ARE NOT OF THE SAME TYPE, HAVE A WEAKER
		      LENNARD JONES POTENTIAL 
		    - CALCULATE THE SUMMATION OF TOTAL CHANGE IN ENERGY 
		      FOR THE PARTICLE IN RELATION TO PARTICLE INTERACTIONS  
		*/

		if(current_prt.getIdentifier() != compare_prt.getIdentifier()){ 

			x_comp   = compare_prt.getX_Position(); 
			y_comp   = compare_prt.getY_Position(); 
			rad_comp = compare_prt.getRadius(); 

			dist_curr = distance(x_curr,x_comp,y_curr,y_comp); 
			dist_temp = distance(x_temp,x_comp,y_temp,y_comp); 

			sigma = (rad_curr + rad_comp); 						// isigma = deal separation distance
																// between particles
			if(current_prt.getType() == compare_prt.getType()){ // interaction betweeen like particles

				energy_curr = 50 * LJ_wellDepth * 
				(pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); // 6-12 potential 

				energy_temp = 50 * LJ_wellDepth * 
				(pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 				
			} 
			else{												// interaction between unlike particles

				energy_curr = 4 * LJ_wellDepth * 
				(pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); // 6-12 potential 

				energy_temp = 4 * LJ_wellDepth * 
				(pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 
			}	

			delta_energy = delta_energy + (energy_temp - energy_curr); 	// running sum of total change
		}																// of current particle's energy 
	}

	return delta_energy; 	// returns the total change in energy 
}

bool Interaction::hardDisks(Particle* particles, int index, int n_particles){ 
	
	Particle current_prt; 
	Particle compare_prt; 						

	double x_temp = 0; 
	double y_temp = 0; 
	double x_comp = 0; 
	double y_comp = 0; 

	double rad_temp	= 0; 
	double rad_comp	= 0; 
	double num		= 0; 

	bool accept = 0; 


	current_prt = particles[index];			// assign the current particle

	x_temp = current_prt.getX_TrialPos(); 	// assign x,y trial position
	y_temp = current_prt.getY_TrialPos(); 	// and the radius of the current
	rad_temp = current_prt.getRadius(); 	// particle 

	accept = 1; 

//////// CHECK FOR PARTICLE-PARTICLE COLLISION //////////
		
for(int k = 0; k < n_particles; k++){

	compare_prt = particles[k]; 		// compare all particles positions
										// to current particle's position	
	if(current_prt.getIdentifier() != compare_prt.getIdentifier()){
											
		x_comp = compare_prt.getX_Position();
		y_comp = compare_prt.getY_Position(); 
		rad_comp = compare_prt.getRadius(); 

		if(distance(x_temp,x_comp,y_temp,y_comp) // if current particle's center
					< rad_comp + rad_temp){		// is closer than the radius of
			accept = 0; 						// current partice plus the radius
			break; 								// of comparison particle, reject move
		}
	}	
}

	return accept; 						// returns 1 if trial move is accepted 
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

		/*	- GENERATE A RANDOM NUMBER FROM [0,1)
			- PROVIDES DIFFERENT 'QUADRANTS' FOR THE PARTICLE TO BE GENERATED IN 
			- ASSIGNS THE NEAREST X,Y 'WALLS'
		*/

		num = randVal.RandomUniformDbl(); 

		if(k % 2 == 0 && num <.5){				// the acceptance presented in 
			x_temp = -1 * x_temp; 				// 'check collision' is not implemented here
												// since the initial position cannot exceed 1
			x_wall = -1 * wall_bound;
			y_wall = wall_bound;  				// randomly assigns a negative/positive
		}										// value to the generated position
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

		accept = 1; 		// the following statements can only change 
							// accept to 0
		if(abs(x_temp - x_wall) < rad_temp){	// if particle's center is closer than 
			accept = 0;							// one radius to the wall, reject move
		}
		else if(abs(y_temp - y_wall) < rad_temp){
			accept = 0; 
		}
		else if(k != 0){						
			for(int n = 0; n < k; n++){
				compare_prt = particles[n]; 			// assign comparison particle 

				x_comp = compare_prt.getX_Position();	// assign the comparison x,y position
				y_comp = compare_prt.getY_Position(); 	// and radius
				rad_comp = compare_prt.getRadius(); 

				if(distance(x_temp,x_comp,y_temp,y_comp) // if the distance between particles is 
					        < rad_comp + rad_temp){		// is closer than the sum of the radii, 
					accept = 0; 						// reject position
					break; 
				}
			}
		}

		if(accept == 1){						// if the position is accepted, assign the 
			current_prt.setX_Position(x_temp);	// x,y position to current particle
			current_prt.setY_Position(y_temp); 

			particles[k] = current_prt; 		// put initialized particle into array 
		}
		else{
			k = k - 1; 			// generate a new random position for the SAME particle
		}

	}
}

////// DEFAULT CONSTRUCTOR ////////////////

Interaction::Interaction(std::string yamlFile){

	YAML::Node node = YAML::LoadFile(yamlFile);

	LJ_wellDepth = node["wellDepth"].as<double>(); 	
}



/////// GETTERS AND SETTERS ///////////////

double Interaction::getWellDepth(){	// this is a purely 'internal' variable
	return LJ_wellDepth; 			// and thus a getter in unnecessary
}

// void Interaction::setWellDepth(double w){
// 	LJ_wellDepth = w; 
// }