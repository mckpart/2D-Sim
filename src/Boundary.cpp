#include <iostream>
#include <cmath>

#include "Boundary.h"

double dist_wall(double x1,double x2){
   return fabs(x2 - x1); 
}

double dist_part(double x1,double x2,double y1, double y2){
   return sqrt(pow(x2 - x1,2) + pow(y2 - y1,2));   
}

bool Boundary::rigidBoundary(std::vector<Particle>* particles, int index){

   Particle current_prt; 

   double x_wall = 0; 
   double y_wall = 0; 	

   double x_temp = 0; 
   double y_temp = 0; 
   double rad_temp = 0; 

   double num = 0; 
   bool accept = 0;	

   current_prt = (*particles)[index];

   x_temp = current_prt.getX_TrialPos(); // reads in the current x,y trial position
   y_temp = current_prt.getY_TrialPos(); // and the radius of the trial particle
   rad_temp = current_prt.getRadius(); 

   accept = 1; 

   if(x_temp * -1 < 0){		
      x_wall = boxLength;  // x_temp > 0, closest x_wall > 0

      if(x_temp - x_wall > 0){ // the particle moved past the wall 
         accept = 0;           // along the x-boundary
      }
   }
   else{
      x_wall = -1 * boxLength; //x_temp < 0, closest wall < 0 

      if(x_temp - x_wall < 0){
         accept = 0; 
      }
   }

   if(y_temp * -1 < 0){			
      y_wall = boxLength;   // y_temp > 0, closest y_wall > 0

      if(y_temp - y_wall > 0){  // the particle moved past the wall 
         accept = 0;            // along the y-boundary
      }
   }
   else{
      y_wall = -1 * boxLength; 

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


   /*   - FINDS THE NEAREST X,Y WALLS
        - IF THE PARTICLE HAS MOVED PAST THE NEAREST WALL,
          COMPUTE THE DISTANCE TRAVELED PAST WALL - RADIUS
        - MOVE PARTICLE TO THE OPPOSITE WALL AND TRAVEL THE 
          REMAINING DISTANCE
   */

   if(x_curr * -1 < 0){               // x > 0, nearest x-wall > 0
      x_wall = boxLength - rad_temp;  // the 'wall' includes the radius
                                      // to make further computations simpler
      if(x_temp - x_wall > 0){

         dist_xwall = dist_wall(x_temp,x_wall); // distance from x trial position to nearest x-wall
         x_temp = -1 * x_wall + dist_xwall;     // moves to opposite wall and travels	 
      }	                                       // the calculated distance
   }
   else{
      x_wall = -1 * (boxLength - rad_temp);    // x < 0, nearest x-wall < 0

      if(x_temp - x_wall < 0){

         dist_xwall  = dist_wall(x_temp,x_wall); 
         x_temp = -1 * x_wall - dist_xwall; 
      } 
   }


   if(y_curr * -1 < 0){	             // y > 0, nearest y-wall > 0
      y_wall = boxLength - rad_temp; 

      if(y_temp - y_wall > 0){

         dist_ywall = dist_wall(y_temp,y_wall); // distance from y trial position 
         y_temp = -1 * y_wall + dist_ywall;     // to the nearest y wall 
      }
   }
   else{
      y_wall = -1 * (boxLength - rad_temp);   // y < 0, nearest y-wall < 0

      if(y_temp - y_wall < 0){

         dist_ywall  = dist_wall(y_temp,y_wall);
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
   x_curr = current_prt.getX_Position(); 
   y_curr = current_prt.getY_Position(); 

   energy_curr = 4 * (pow(x_curr,2) + pow(y_curr,2)); // treat the constants as 1 for now.. change later
   energy_temp = 4 * (pow(x_temp,2) + pow(x_temp,2)); 

   delta_energy = energy_temp - energy_curr; 

   return delta_energy; 
}

void Boundary::initialPosition(std::vector<Particle>* particles, int n_particles, 
                          KISSRNG randVal){
   Particle current_prt; 
   Particle compare_prt; 

   double x_temp = 0; 
   double y_temp = 0; 
   double x_comp = 0; 
   double y_comp = 0; 
   double x_wall = 0; 
   double y_wall = 0;   // the locations of the nearest 'wall'

   double rad_temp = 0; 
   double rad_comp = 0; 
   double num = 0; 

   bool accept = 0; 

   for(int k = 0; k < n_particles; k++){

      current_prt = (*particles)[k];
      rad_temp = current_prt.getRadius();  

      x_temp = randVal.RandomUniformDbl() * boxLength; 
      y_temp = randVal.RandomUniformDbl() * boxLength; 

      /*   - GENERATE A RANDOM NUMBER FROM [0,1)
           - PROVIDES DIFFERENT 'QUADRANTS' FOR THE PARTICLE TO BE GENERATED IN 
           - ASSIGNS THE NEAREST X,Y 'WALLS'
      */

      num = randVal.RandomUniformDbl(); 

      if(k % 2 == 0 && num <.5){      // the acceptance presented in 
         x_temp = -1 * x_temp;        // 'check collision' is not implemented here
                                      // since the initial position cannot exceed 1
         x_wall = -1 * boxLength;
         y_wall = boxLength;           // randomly assigns a negative/positive
      }                                // value to the generated position
      else if(k % 2 == 0 && num >= .5){
         y_temp = -1 * y_temp; 

         x_wall = boxLength; 
         y_wall = -1 * boxLength; 
      }
      else if(k % 2 != 0 && num < .5){
         x_temp = -1 * x_temp; 
         y_temp = -1 * y_temp; 

         x_wall = -1 * boxLength; 
         y_wall = -1 * boxLength; 
      }
      else{
         x_wall = boxLength; 
         y_wall = boxLength; 
      }

      accept = 1;                           // the following statements can only change 
                                            // accept to 0
      if(abs(x_temp - x_wall) < rad_temp){  // if particle's center is closer than 
         accept = 0;                        // one radius to the wall, reject move
      }
      else if(abs(y_temp - y_wall) < rad_temp){
         accept = 0; 
      }
      else if(k != 0){                 
          for(int n = 0; n < k; n++){
             compare_prt = (*particles)[n];           // assign comparison particle 

             x_comp = compare_prt.getX_Position();    // assign the comparison x,y position
             y_comp = compare_prt.getY_Position();    // and radius
             rad_comp = compare_prt.getRadius(); 

             if(dist_part(x_temp,x_comp,y_temp,y_comp) // if the distance between particles is 
                         < rad_comp + rad_temp){      // is closer than the sum of the radii, 
                accept = 0;                           // reject position
                break; 
             }
         }
      }

      if(accept == 1){                      // if the position is accepted, assign the 
         current_prt.setX_Position(x_temp); // x,y position to current particle
         current_prt.setY_Position(y_temp); 

         (*particles)[k] = current_prt;     // put initialized particle into array 
      }
      else{
         k = k - 1; // generate a new random position for the SAME particle
      }
   }
}

void Boundary::initializeBoundary(std::string yamlFile){
   
   YAML::Node node = YAML::LoadFile(yamlFile); 

   boxLength = node["boxLength"].as<double>(); 
}
