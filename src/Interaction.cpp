#include <iostream>
#include "Interaction.h"
#include <cmath>


double distance(double x1,double x2,double y1, double y2){
   return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2)); 	
}
double dist1D(double x1,double x2){
   return fabs(x2 - x1); 
}

double Interaction::WCApotential(std::vector<Particle>* particles, int index){

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
            energy_curr = 4 *  
            (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6) + .25); 
         }

         if(dist_temp > pow(2,1/6) * sigma){
            energy_temp = 0; 
         }
         else{
            energy_temp = 4 * 
            (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6) + .25); 								
         }

         delta_energy = delta_energy + (energy_temp - energy_curr); 
      }
   }

   return delta_energy;  // returns the total change in energy associated with this move
}

void Interaction::populateCellArray(double x,double y, std::vector<std::vector<double>>* cellPositions){
   
//   std::cout << "x: " << x << " y: " <<  y << std::endl; 

   (*cellPositions)[0][0] = x;             (*cellPositions)[0][1] = y; 
   (*cellPositions)[1][0] = x;             (*cellPositions)[1][1] = y + boxLength; 
   (*cellPositions)[2][0] = x;             (*cellPositions)[2][1] = y - boxLength; 
   (*cellPositions)[3][0] = x + boxLength; (*cellPositions)[3][1] = y; 
   (*cellPositions)[4][0] = x + boxLength; (*cellPositions)[4][1] = y + boxLength; 
   (*cellPositions)[5][0] = x + boxLength; (*cellPositions)[5][1] = y - boxLength; 
   (*cellPositions)[6][0] = x - boxLength; (*cellPositions)[6][1] = y; 
   (*cellPositions)[7][0] = x - boxLength; (*cellPositions)[7][1] = y + boxLength; 
   (*cellPositions)[8][0] = x - boxLength; (*cellPositions)[8][1] = y - boxLength; 

}

double Interaction::lennardJones(std::vector<Particle>* particles, int index){

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

   std::vector<std::vector<double>> cellPositions(9,std::vector<double>(2,0));  

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

   d_curr_wall_x = 0.5 * boxLength - fabs(x_curr); // x,y distance from
   d_curr_wall_y = 0.5 * boxLength - fabs(y_curr); // nearets walls
   
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

	 dist_curr_tot = distance(x_curr,x_comp,y_curr,y_comp); 
	 dist_temp_tot = distance(x_temp,x_comp,y_temp,y_comp); 
	   
//	 std::cout << "curr dist: " << dist_curr_tot << " temp dist: " << dist_temp_tot << std::endl;
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
	    if(dist_curr_tot > truncDist || dist_temp_tot > truncDist){
              
	       populateCellArray(x_comp,y_comp,&cellPositions);
               for(int z = 0; z < 9; z++){
               
	          x_comp = cellPositions[z][0]; 
	          y_comp = cellPositions[z][1];

	          dist_curr_tot = distance(x_curr,x_comp,y_curr,y_comp); 
                  dist_temp_tot = distance(x_temp,x_comp,y_temp,y_comp);  
                   
		  if(dist_curr_tot < truncDist){
                     energy_curr = 4 *         // 6-12 potential 
                     (pow(sigma/dist_curr_tot,12) - pow(sigma/dist_curr_tot,6) + truncShift); // this will calulate the reduced 
		  }
                  else{
		     energy_curr = 0; 
		  }
		  
		  if(dist_temp_tot < truncDist){
                     energy_temp = 4 *        // 6-12 potential 
                     (pow(sigma/dist_temp_tot,12) - pow(sigma/dist_temp_tot,6) + truncShift); 				
		  }
		  else{
		     energy_temp = 0; 
		  }

		  delta_energy = delta_energy + (energy_temp - energy_curr); // running summation of change in energy
	    
	       }
	    }
	    else{
	       energy_curr = 4*(pow(sigma/dist_curr_tot,12) - pow(sigma/dist_curr_tot,6) + truncShift); 
	       energy_temp = 4*(pow(sigma/dist_temp_tot,12) - pow(sigma/dist_temp_tot,6) + truncShift);
	    
	       delta_energy = delta_energy + (energy_temp - energy_curr); 
	    }
         } 
      }                                                              // of current particle's energy 
   }
   return delta_energy + tail_corr;      // returns the total change in energy 
}

bool Interaction::hardDisks(std::vector<Particle>* particles, int index){ 
	
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

////// DEFAULT CONSTRUCTOR ////////////////

void Interaction::initializeInteraction(Parameters* p){

   boxLength    = p->getBoxLength(); 
   n_particles  = p->getNumParticles(); 
   sigma        = p->getSigma(); 
   redDens      = p->getRedDens(); 

   truncDist = 2.5 * sigma;
   truncShift = -1 * (pow(sigma/truncDist,12) 
		    - pow(sigma/truncDist,6));

   tail_corr =  3.141592654 * redDens * (.4 * pow(sigma/truncDist,10) 
		                 - pow(sigma/truncDist,4)); 

   std::cout << "the shifted potential value is: " << truncShift << std::endl; 
}

