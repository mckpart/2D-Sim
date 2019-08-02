#include <iostream>

#include "Interaction.h"
#include <cmath>


double Interaction::distance(double x1,double x2,double y1, double y2){
   return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2))/sigma; // returns the characteristic length 	
}

double Interaction::lenjones_energy(double r, double c){
   return 4*c*(pow(1/r,12) - pow(1/r,6) + truncShift);
}

double Interaction::WCA_energy(double r, double c){
   return 4*c*(pow(1/r,12) - pow(1/r,6) + .25);
}

void Interaction::populateCellArray(double x,
		                    double y, 
				    std::vector<std::vector<double>>* cellPositions){
   // creates the cell images of the comparison particle... should be cleaned up
 
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

double Interaction::periodicInteraction(std::vector<Particle>* particles, 
		                        int index){
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

   double num = 0; 
   double c = 0; 

   double dist_curr_x = 0;
   double dist_curr_y = 0;  
   double dist_temp_x = 0; 	
   double dist_temp_y = 0; 

   double dist_curr_tot = 0; 
   double dist_temp_tot = 0; 

   std::vector<std::vector<double>> cellPositions(9,std::vector<double>(2,0));  

   current_prt = (*particles)[index];    // assign current particle
 
   x_temp = current_prt.getX_TrialPos(); // assign the current and trial
   y_temp = current_prt.getY_TrialPos(); // positions of the current
                                         // particle
   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 
   
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

         x_comp = compare_prt.getX_Position(); // set the comparison  
         y_comp = compare_prt.getY_Position(); // particles position

	 dist_curr_tot = distance(x_curr,x_comp,y_curr,y_comp); 
	 dist_temp_tot = distance(x_temp,x_comp,y_temp,y_comp); 
	   
	 if(current_prt.getType() == compare_prt.getType()){ // interaction betweeen like particles
	    c = LJ_par; 
	 }
         else if(current_prt.getType() != compare_prt.getType()){
	    c = LJ_antipar; 
	 }	 
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
               
	       x_comp = cellPositions[z][0];  // creates the 'phantom' particles 
	       y_comp = cellPositions[z][1];  // in the other cell images

	       dist_curr_tot = distance(x_curr,x_comp,y_curr,y_comp); 
               dist_temp_tot = distance(x_temp,x_comp,y_temp,y_comp);  
                   
     	       if(dist_curr_tot < truncDist){ // maybe have a different truncation dist for different 
		  switch(interact_type){      // types of interactions...
		     case 0: std::cout << "mehhh" << std::endl;	  // also different energy corrections should
                             break;                               // be included for the different interaction
		     case 1: energy_curr = lenjones_energy(dist_curr_tot,c); // types
			     break;
	             case 2: energy_curr = WCA_energy(dist_curr_tot,c);
			     break;  
	          }
	       }
               else{
                  energy_curr = 0; 
	       }
		  
	       if(dist_temp_tot < truncDist){
                  switch(interact_type){   
		     case 0: std::cout << "mehhh" << std::endl;	  
                             break;                               
		     case 1: energy_temp = lenjones_energy(dist_temp_tot,c); 
			     break;
	             case 2: energy_temp = WCA_energy(dist_temp_tot,c); 
			     break;
		  }
	       }
	       else{
		 energy_temp = 0; 
	       }
	       
	       delta_energy = delta_energy + (energy_temp - energy_curr); // running summation of change in energy
	    }
	 }
	 else{

            switch(interact_type){      // types of interactions...
	       case 0: std::cout << "mehhh" << std::endl;	  // also different energy corrections should
                       break;                               // be included for the different interaction
	       case 1: energy_curr = lenjones_energy(dist_curr_tot,c); // types
                       energy_temp = lenjones_energy(dist_temp_tot,c); 
		       break;
	       case 2: energy_curr = WCA_energy(dist_curr_tot,c);
		       energy_temp = WCA_energy(dist_temp_tot,c); 
		       break;
	    }
	    delta_energy = delta_energy + (energy_temp - energy_curr); 
	    
	 }
      }                                                                
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


void Interaction::initializeInteraction(Parameters* p){

   boxLength    = p->getBoxLength();           // assign all private variables
   n_particles  = p->getNumParticles();        // used in this class
   sigma        = p->getSigma(); 
   redDens      = p->getRedDens(); 
  
   interact_type = p->getInteract_Type();
   
   LJ_par       = p->getLJ_const_1(); 
   LJ_antipar   = p->getLJ_const_2(); 

   truncDist = 2.5;    // this is really 2.5 * sigma / sigma
   truncShift = -1 * (pow(1/truncDist,12) 
		    - pow(1/truncDist,6));

   tail_corr =  3.141592654 * redDens * (.4 * pow(1/truncDist,10) 
		                 - pow(1/truncDist,4)); 
}

