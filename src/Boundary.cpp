#include <iostream>
#include <cmath>

#include "Boundary.h"

double distance(double x1,double x2){
   return fabs(x2 - x1); 
}

bool Boundary::rigidBoundary(std::vector<Particle>* particles, int index){

   Particle current_prt; 

   double x_wall = 0; 
   double y_wall = 0; 	

   double x_temp = 0; 
   double y_temp = 0; 
   double rad_temp = 0; 

   double num = 0; 
   double wall_bound = 0; 	

   bool accept = 0;	

   current_prt = (*particles)[index];

   x_temp = current_prt.getX_TrialPos(); // reads in the current x,y trial position
   y_temp = current_prt.getY_TrialPos(); // and the radius of the trial particle
   rad_temp = current_prt.getRadius(); 
   wall_bound = 1; 	

   accept = 1; 

   if(x_temp * -1 < 0){		
      x_wall = wall_bound;  // x_temp > 0, closest x_wall > 0

      if(x_temp - x_wall > 0){ // the particle moved past the wall 
         accept = 0;           // along the x-boundary
      }
   }
   else{
      x_wall = -1 * wall_bound; //x_temp < 0, closest wall < 0 

      if(x_temp - x_wall < 0){
         accept = 0; 
      }
   }

   if(y_temp * -1 < 0){			
      y_wall = wall_bound;   // y_temp > 0, closest y_wall > 0

      if(y_temp - y_wall > 0){  // the particle moved past the wall 
         accept = 0;            // along the y-boundary
      }
   }
   else{
      y_wall = -1 * wall_bound; 

      if(y_temp - y_wall < 0){
         accept = 0;
      }
   }

   if(fabs(x_temp - x_wall) < rad_temp){   // the particle's center of mass is less than 
      accept = 0;                          // a radius' distance from the x-boundary
   }	
   else if(fabs(y_temp - y_wall) < rad_temp){  // is less than a radius' distance from the 
      accept = 0;                              // y-boundary
   }

   return accept;   // returns 1 if trial move is accepted
}

void Boundary::periodicBoundary(std::vector<Particle>* particles, int index){

   Particle current_prt; 

   double x_wall = 0; 
   double y_wall = 0; 	

   double x_curr = 0; 
   double y_curr = 0;
   double x_temp = 0; 
   double y_temp = 0; 
   double rad_temp = 0; 

   double num = 0; 
   double wall_bound = 0; 	

   double dist_xtravel = 0; 
   double dist_ytravel = 0; 
   double dist_xwall = 0; 
   double dist_ywall = 0;

   current_prt = (*particles)[index];

   x_curr = current_prt.getX_Position(); // sets the current particle's x,y
   y_curr = current_prt.getY_Position(); // position
   x_temp = current_prt.getX_TrialPos(); // set the x,y trial position for 
   y_temp = current_prt.getY_TrialPos(); // current particle
   rad_temp = current_prt.getRadius(); 

   wall_bound = 1;  // will be added as private variable


   /*   - FINDS THE NEAREST X,Y WALLS
        - IF THE PARTICLE HAS MOVED PAST THE NEAREST WALL,
          COMPUTE THE DISTANCE TRAVELED PAST WALL - RADIUS
        - MOVE PARTICLE TO THE OPPOSITE WALL AND TRAVEL THE 
          REMAINING DISTANCE
   */

   if(x_curr * -1 < 0){               // x > 0, nearest x-wall > 0
      x_wall = wall_bound - rad_temp; // the 'wall' includes the radius
                                      // to make further computations simpler
      if(x_temp - x_wall > 0){

         dist_xwall = distance(x_temp,x_wall); 	// distance from x trial position to nearest x-wall
         x_temp = -1 * x_wall + dist_xwall;     // moves to opposite wall and travels	 
      }	                                        // the calculated distance
   }
   else{
      x_wall = -1 * (wall_bound - rad_temp);    // x < 0, nearest x-wall < 0

      if(x_temp - x_wall < 0){

         dist_xwall  = distance(x_temp,x_wall); 
         x_temp = -1 * x_wall - dist_xwall; 
      } 
   }


   if(y_curr * -1 < 0){	             // y > 0, nearest y-wall > 0
      y_wall = wall_bound - rad_temp; 

      if(y_temp - y_wall > 0){

         dist_ywall = distance(y_temp,y_wall); // distance from y trial position 
         y_temp = -1 * y_wall + dist_ywall;    // to the nearest y wall 
      }
   }
   else{
      y_wall = -1 * (wall_bound - rad_temp);   // y < 0, nearest y-wall < 0

      if(y_temp - y_wall < 0){

         dist_ywall  = distance(y_temp,y_wall);
         y_temp = -1 * y_wall - dist_ywall; 
      } 
   }

   current_prt.setX_TrialPos(x_temp); 	// puts the updated particle with updated 
   current_prt.setY_TrialPos(y_temp); 	// trial positions back into the array 

   (*particles)[index] = current_prt; 

}

double Boundary::externalWell(std::vector<Particle>* particles, int index){

   Particle current_prt; 

   double x_temp = 0; 
   double y_temp = 0; 
   double x_curr = 0; 
   double y_curr = 0; 

   double delta_energy = 0; 
   double energy_curr  = 0; 
   double energy_temp  = 0; 

   double num = 0; 

   double dist_curr = 0; 
   double dist_temp = 0; 	

   current_prt = (*particles)[index];   // assign current particle

   x_temp = current_prt.getX_TrialPos(); // assign the current and trial
   y_temp = current_prt.getY_TrialPos(); // positions and the radius of 
                                         // the current particle
   x_curr   = current_prt.getX_Position(); 
   y_curr   = current_prt.getY_Position(); 

   energy_curr = 10 * (pow(x_curr,2) + pow(y_curr,2)); // treat the constants as 1 for now.. change later
   energy_temp = 10 * (pow(x_temp,2) + pow(x_temp,2)); 

   delta_energy = energy_temp - energy_curr; 

   return delta_energy; 
}
