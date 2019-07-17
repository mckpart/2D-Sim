#include <iostream>
#include "Interaction.h"
#include <cmath>


double distance(double x1,double x2,double y1, double y2){
   return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2)); 	
}
double dist1D(double x1,double x2){
   return fabs(x2 - x1); 
}

double Interaction::WCApotential(std::vector<Particle>* particles, int index, int n_particles){

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
   double energy_curr = 0; 
   double energy_temp = 0; 

   double rad_curr = 0; 
   double rad_comp = 0; 
   double num = 0; 

   double dist_curr = 0; 
   double dist_temp = 0; 	

   current_prt = (*particles)[index];  // assign current particle

   x_temp = current_prt.getX_TrialPos(); // assign the current and trial
   y_temp = current_prt.getY_TrialPos(); // positions and the radius of 
                                         // the current particle
   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 
   rad_curr = current_prt.getRadius(); 

   for(int k = 0; k < n_particles; k++){  // the O(N^2)... think about 
                                          // implementing neighbor lists here
      compare_prt = (*particles)[k];      // assign comparison particle


      /*   - CHECK THAT THE COMPARISON PARTICLE IS NOT THE CURRENT PARTICLE 
           - SET THE COMPARISON POSITION AND RADIUS
           - IF THE DISTANCE BETWEEN PARTICLES IS > SIGMA * 2^(1/6), THE CHANGE
             IN ENERGY IS ZERO
           - IF LESS THAN SIGMA * 2^(1/6), THE LENNARD JONES REPULSION IS USED
           - SUM THE TOTAL CHANGE IN ENERGY
      */	 

      if(current_prt.getIdentifier() != compare_prt.getIdentifier()){ 

         x_comp = compare_prt.getX_Position(); 
         y_comp = compare_prt.getY_Position(); 
         rad_comp = compare_prt.getRadius(); 

         dist_curr = distance(x_curr,x_comp,y_curr,y_comp); 
         dist_temp = distance(x_temp,x_comp,y_temp,y_comp); 

         sigma = (rad_curr + rad_comp); 				

         if(dist_curr > pow(2,1/6) * sigma){  // the potential does not extend past 
            energy_curr = 0;                  // this maximum distance - there is no
         }                                    // change in energy beyond this distance
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

   return delta_energy;  // returns the total change in energy associated with this move
}

double Interaction::lennardJones(std::vector<Particle>* particles, int index, int n_particles){

   Particle current_prt; 
   Particle compare_prt; 						

   double x_temp = 0; 
   double y_temp = 0; 
   double x_curr = 0; 
   double y_curr = 0; 
   double x_comp = 0; 
   double y_comp = 0; 

   double delta_energy = 0; 
   double energy_curr = 0; 
   double energy_temp = 0; 

   double rad_curr = 0; 
   double rad_comp = 0; 
   double num = 0; 

   double dist_curr_x = 0;
   double dist_curr_y = 0;  
   double dist_temp_x = 0; 	
   double dist_temp_y = 0; 

   double dist_curr_tot = 0; 
   double dist_temp_tot = 0; 

   double d_curr_wall_x = 0; 
   double d_curr_wall_y = 0; 
   double d_temp_wall_x = 0; 
   double d_temp_wall_y = 0; 
   double d_comp_wall_x = 0; 
   double d_comp_wall_y = 0; 

   double dist_trunc = 0; 

   current_prt = (*particles)[index];    // assign current particle
 
   x_temp = current_prt.getX_TrialPos(); // assign the current and trial
   y_temp = current_prt.getY_TrialPos(); // positions and the radius of 
                                         // the current particle
   d_temp_wall_x = boxLength - fabs(x_temp); // x,y distance from 
   d_temp_wall_y = boxLength - fabs(y_temp); // nearest walls

   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 
   rad_curr = current_prt.getRadius(); 

   d_curr_wall_x = boxLength - fabs(x_curr); // x,y distance from
   d_curr_wall_y = boxLength - fabs(y_curr); // nearets walls
   
   for(int k = 0; k < n_particles; k++){

      compare_prt = (*particles)[k];   // assign the comparison particle


      /* - CHECK THAT THE COMPARISON PARTICLE IS NOT THE CURRENT PARTICLE
         - SET THE COMPARISON POSITION AND RADIUS
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

         x_comp = compare_prt.getX_Position(); 
         y_comp = compare_prt.getY_Position(); 
         rad_comp = compare_prt.getRadius(); 
         
	 d_comp_wall_x = boxLength - fabs(x_comp);
         d_comp_wall_y = boxLength - fabs(y_comp);
         
	 dist_curr_x = dist1D(x_curr,x_comp);
	 dist_curr_y = dist1D(y_curr,y_comp); 
         dist_temp_x = dist1D(x_temp,x_comp); 
	 dist_temp_y = dist1D(y_temp,y_comp); 


         if(current_prt.getType() == compare_prt.getType()){ // interaction betweeen like particles
	    
            /* IF THE SUMMATION OF THE X DISTANCES FROM THE WALL IS WITHIN THE 
	     *    DISTANCE OF INTERACTION AND THE PARTICLES ARE NOT ON THE SAME SIDE
	     *    OF THE BOX, UPDATE THE X DISTANCE BETWEEN THE PARTICLES
	     * IF THE SUMMATION OF THE Y DISTANCES FROM THE WALL IS WITHIN THE
	     *    DISTANCE OF INTERACTION AND THE PARTICLES ARE NOT ON THE SAME SIDE
	     *    OF THE BOX, UPDATE THE Y DISTAANCE BETWEEN THE PARTICLES
	     * COMPUTE THE RADIAL DISTANCE BETWEEN THE CURRENT PARTICLE AND 
	     *    COMPARISON PARTICLE ONCE X,Y DISTANCES ARE UPDATED ACCORDINGLY
	     */
		 
	    if(d_curr_wall_x + d_comp_wall_x < truncDist && x_curr * x_comp < 0){    // truncation dist = sigma * 2.5
               dist_curr_x = 2 * boxLength - dist_curr_x;   
	    }
	    if(d_curr_wall_y + d_comp_wall_y < truncDist && y_curr * y_comp < 0){
	       dist_curr_y = 2 * boxLength - dist_curr_y; 
	    }

            if(d_temp_wall_x + d_comp_wall_x < truncDist && x_temp * x_comp < 0){
	       dist_temp_x = 2 * boxLength - dist_temp_x; 
	    }
	    if(d_temp_wall_y + d_comp_wall_y < truncDist && y_temp * y_comp < 0){
	       dist_temp_y = 2 * boxLength - dist_temp_y; 
	    }

	    dist_curr_tot = sqrt(pow(dist_curr_x,2) + pow(dist_curr_y,2)); 
            dist_temp_tot = sqrt(pow(dist_temp_x,2) + pow(dist_temp_y,2));  

	    if(dist_curr_tot < truncDist){
               energy_curr = 4 * LJ_wellDepth *        // 6-12 potential 
               (pow(sigma/dist_curr_tot,12) - pow(sigma/dist_curr_tot,6)); 
	    }
	    else{
	       energy_curr = 0;
	    }

	    if(dist_temp_tot < truncDist){
               energy_temp = 4 * LJ_wellDepth *        // 6-12 potential 
               (pow(sigma/dist_temp_tot,12) - pow(sigma/dist_temp_tot,6)); 				
	    }
	    else{
	       energy_temp = 0; 
	    }   
	 } 
//         else{                                             // interaction between unlike particles

//           energy_curr = 4 * LJ_wellDepth * 
//           (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); // 6-12 potential 

//            energy_temp = 4 * LJ_wellDepth * 
//            (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 
//         }	

         delta_energy = delta_energy + (energy_temp - energy_curr);  // running sum of total change
      }                                                              // of current particle's energy 
   }
   return delta_energy;      // returns the total change in energy 
}

bool Interaction::hardDisks(std::vector<Particle>* particles, int index, int n_particles){ 
	
   Particle current_prt; 
   Particle compare_prt; 						

   double x_temp = 0; 
   double y_temp = 0; 
   double x_comp = 0; 
   double y_comp = 0; 

   double rad_temp = 0; 
   double rad_comp = 0; 
   double num = 0; 

   bool accept = 0; 

   current_prt = (*particles)[index];   // assign the current particle

   x_temp = current_prt.getX_TrialPos(); // assign x,y trial position
   y_temp = current_prt.getY_TrialPos(); // and the radius of the current
   rad_temp = current_prt.getRadius();   // particle 

   accept = 1; 

//////// CHECK FOR PARTICLE-PARTICLE COLLISION //////////
		
   for(int k = 0; k < n_particles; k++){

      compare_prt = (*particles)[k]; // compare all particles positions
                                     // to current particle's position	
      if(current_prt.getIdentifier() != compare_prt.getIdentifier()){
										
         x_comp = compare_prt.getX_Position();
         y_comp = compare_prt.getY_Position(); 
         rad_comp = compare_prt.getRadius(); 

         if(distance(x_temp,x_comp,y_temp,y_comp) // if current particle's center
                     < rad_comp + rad_temp){      // is closer than the radius of
            accept = 0;                           // current partice plus the radius
            break;                                // of comparison particle, reject move
         }
      }	
   }
   return accept;     // returns 1 if trial move is accepted 
}

double Interaction::crosslinkers(std::vector<Particle>* particles, int index, int n_particles){
   
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
   double energy_curr = 0; 
   double energy_temp = 0; 

   double rad_curr = 0; 
   double rad_comp = 0; 
   double num = 0; 

   double dist_curr = 0; 
   double dist_temp = 0;   

   current_prt = (*particles)[index];    // assign current particle
 
   x_temp = current_prt.getX_TrialPos(); // assign the current and trial
   y_temp = current_prt.getY_TrialPos(); // positions of current particle 
                                         
   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 

   for(int k = 0; k < n_particles; k++){

      compare_prt = (*particles)[k]; // assign the comparison particle 
                                     // comparison particle cannot be current
      if(compare_prt.getIdentifier() != current_prt.getIdentifier()){

         x_comp = compare_prt.getX_Position();  // assign x,y comparison 
         y_comp = compare_prt.getY_Position();  // position 

         dist_curr = distance(x_curr,x_comp,y_curr,y_comp); // current distance 
                                                            // between particles
	 if(dist_curr < truncDist){                                             
   	    energy_curr = 1/beta * exp(-.5 * sprConstant * beta * 
		        pow(dist_curr - restLength,2)); 
         }
	 else{
	    energy_curr = 0; 
	 }

	 dist_temp = distance(x_temp,x_comp,y_temp,y_comp); // trial distance  
                                                            // between particles 
	if(dist_temp < truncDist){
	   energy_temp = 1/beta * exp(-.5 * sprConstant * beta * 
		         pow(dist_temp - restLength,2));  
	} 
        else{
	   energy_temp = 0; 
	}	

	delta_energy = delta_energy + (energy_temp - energy_curr); // running sum of
      }                                                             // total change in 
   }                                                                // energy 
   // std::cout << "the change in energy is: " << delta_energy << std::endl; 
   return delta_energy;  // returns total change in energy 
}

////// DEFAULT CONSTRUCTOR ////////////////

void Interaction::initializeInteraction(std::string yamlFile){

   YAML::Node node = YAML::LoadFile(yamlFile);

   LJ_wellDepth = node["wellDepth"].as<double>();
   sigma        = node["sigma"].as<double>();    
   restLength   = node["restLength"].as<double>();
   sprConstant  = node["springConstant"].as<double>(); 
//   truncDist    = node["truncationDist"].as<double>();
   beta         = node["beta"].as<double>();   
   boxLength    = node["boxLength"].as<double>(); 

   truncDist = 2.5 * sigma; 
}

