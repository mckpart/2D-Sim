#include <iostream>
#include "Interaction.h"
#include <cmath>


double distance(double x1,double x2,double y1, double y2){
   return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2)); 	
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
   y_temp = current_prt.getY_TrialPos(); // positions and the radius of 
                                         // the current particle
   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 
   rad_curr = current_prt.getRadius(); 

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

         dist_curr = distance(x_curr,x_comp,y_curr,y_comp); 
         dist_temp = distance(x_temp,x_comp,y_temp,y_comp); 

         sigma = (rad_curr + rad_comp) * pow(2,.5);                      // sigma = deal separation distance
                                                             // between particles
         if(current_prt.getType() == compare_prt.getType()){ // interaction betweeen like particles

            energy_curr = 10 * LJ_wellDepth * 
            (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); // 6-12 potential 

            energy_temp = 10 * LJ_wellDepth * 
            (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 				
         } 
         else{                                             // interaction between unlike particles

            energy_curr = 4 * LJ_wellDepth * 
            (pow(sigma/dist_curr,12) - pow(sigma/dist_curr,6)); // 6-12 potential 

            energy_temp = 4 * LJ_wellDepth * 
            (pow(sigma/dist_temp,12) - pow(sigma/dist_temp,6)); 
         }	

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
   y_temp = current_prt.getY_TrialPos(); // positions and the radius of 
                                         // the current particle
   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 
   // rad_curr = current_prt.getRadius(); 

   for(int k = 0; k < n_particles; k++){

      compare_prt = (*particles)[k]; // assign the comparison particle 
                                     // comparison particle cannot be current
      if(compare_prt.getIdentifier() != current_prt.getIdentifier()){

         x_comp = compare_prt.getX_Position();  // assign x,y comparison 
         y_comp = compare_prt.getY_Position();  // position 
         // rad_comp = compare_prt.getRadius(); 

         dist_curr = distance(x_curr,x_comp,y_curr,y_comp); // current distance 
                                                            // between particles
         if(dist_curr > truncDist){    // particles do not interact if current
            energy_curr = 0;           // distance is greater than the truncation
	 }                             // distance
         else{
	    energy_curr = 0.5 * sprConstant * pow(dist_curr - restLength,2);  
	 }

	 dist_temp = distance(x_temp,x_comp,y_temp,y_comp); // trial distance  
                                                            // between particles
         if(dist_temp > truncDist){ // particles do not interact if trial 
            energy_temp = 0;        // trial distance is greater than the 
	 }                          // truncation distance
         else{
	    energy_temp = 0.5 * sprConstant * pow(dist_temp - restLength,2); 
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
   restLength   = node["restLength"].as<double>();
   sprConstant  = node["springConstant"].as<double>(); 
   truncDist    = node["truncationDist"].as<double>();  
}

