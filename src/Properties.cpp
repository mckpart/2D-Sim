#include <iostream>
#include "Properties.h" 

double radDistance(double deltaX, double deltaY){
   return sqrt(pow(deltaX,2) + pow(deltaY,2));  
}   
double Properties::lenJonesForce(double r){
   return 24 * epsilon/sigma * (2 * pow(sigma/r,13) - pow(sigma/r,7));  // come back to check this calculation
}
double Properties::lenJonesEnergy(double r){
   return 4 * epsilon * (pow(sigma/r,12) - pow(sigma/r,6)); 
}

void Properties::calcEnergy(double r_dist){
   f_energy = f_energy + lenJonesEnergy(r_dist); 
}

void Properties::calcForces(std::vector<Particle>* particles){ // this needs to be reorganized
   Particle curr_prt; 
   Particle comp_prt; 

   double x_curr = 0;
   double y_curr = 0; 
   double x_comp = 0; 
   double y_comp = 0; 

   double d_curr_wall_x = 0; 
   double d_curr_wall_y = 0; 
   double d_comp_wall_x = 0; 
   double d_comp_wall_y = 0;

   double dist_curr_x = 0; 
   double dist_curr_y = 0; 

   double x_force = 0; 
   double y_force = 0; 

   double r_dist = 0; 
   double force_tot = 0; 
   double f_r = 0; 

   f_energy = 0;   // make sure that the free energy previously calculated is reset
  
   for(int k = 0; k < n_particles; k++){
      curr_prt = (*particles)[k]; 
      
      x_curr = curr_prt.getX_Position(); 
      y_curr = curr_prt.getY_Position(); 
   
      d_curr_wall_x = boxLength - fabs(x_curr);
      d_curr_wall_y = boxLength - fabs(y_curr); 

      for(int n = k + 1; n < n_particles; n++){ // each particle-particle interaction
         comp_prt = (*particles)[n]; 

	 x_comp = comp_prt.getX_Position(); 
	 y_comp = comp_prt.getY_Position(); 

         d_comp_wall_x = boxLength - fabs(x_comp);
         d_comp_wall_y = boxLength - fabs(y_comp); 
	
         dist_curr_x = x_comp - x_curr; 
	 dist_curr_y = y_comp - y_curr;
	 
	 if(LJ == 1){
            
            if(d_curr_wall_x + d_comp_wall_x < truncDist && x_curr * x_comp < 0){    // truncation dist = sigma * 2.5
	       dist_curr_x = d_curr_wall_x + d_comp_wall_x; 
	    }

	    if(d_curr_wall_y + d_comp_wall_y < truncDist && y_curr * y_comp < 0){
	       dist_curr_y = d_curr_wall_y + d_comp_wall_y; 
	    }
	     
	    r_dist = radDistance(dist_curr_x, dist_curr_y);

	    if(r_dist < truncDist){
	       force_tot = lenJonesForce(r_dist); 	    
	    }
	    else{
	       force_tot = 0; 
	    }

	 }                       
         
	 calcEnergy(r_dist);  // this allows the energy to be calculated at the same time as
                              // the forces. If this method works, implement the other method
	 f_r = f_r + force_tot * r_dist; 
      }
   }
   sum_Fdot_r.push_back(f_r); 
   sum_energy.push_back(f_energy); 
}


double Properties::calcPressure(){
   double len = 0; 
   double avgEnergy = 0; 
   double redPressure = 0; 

   len = double(sum_Fdot_r.size());
   
   for(int k = 0; k < len; k++){
      avgEnergy = avgEnergy + sum_Fdot_r[k];  
   }
   avgEnergy = avgEnergy/len;

   redPressure = redDens * (redTemp + avgEnergy / (epsilon * 2 * n_particles)); 
 
   return redPressure; 
}

void Properties::writeProperties(){
   
   double len = 0; 

   std::ofstream force_file; 
   std::ofstream energy_file; 

   force_file.open("forces.txt");
   energy_file.open("energies.txt"); 

   len = double(sum_Fdot_r.size()); // the force and energy vector are the same 
                                    // size hence are put into one for-loop
   for(int k = 0; k < len; k++){
      force_file << sum_Fdot_r[k] << " ";
      energy_file << sum_energy[k] << " "; 
   }
}
void Properties::initializeProperties(std::string yf){
   YAML::Node node = YAML::LoadFile(yf);

   sigma   = node["sigma"].as<double>();
   KbT     = 1/ node["beta"].as<double>(); 
   epsilon = node["wellDepth"].as<double>();

   boxLength   = node["boxLength"].as<double>();
   n_particles = node["totalParticles"].as<double>(); 

   LJ      = node["lennardJones"].as<bool>();

   redDens = (n_particles/(pow(2 * boxLength,2)) * pow(sigma,2)); 
   redTemp = KbT / epsilon; 

   truncDist = 2.5 * sigma; 

//   std::cout << "the reduced density is: " << redDens 
//	     << " and red temp is: " << redTemp << std::endl; 
}
