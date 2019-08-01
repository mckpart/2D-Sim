#include <iostream>
#include "Properties.h" 

// returns the characteristic length between two particles
double Properties::radDistance(double x1, double x2, 
		               double y1, double y2){
   return sqrt(pow(x2-x1,2) + pow(y2-y1,2)) / sigma;  
}   

double Properties::lenJonesForce(double r, double c){
   return 24*c/sigma * (2 * pow(1/r,13) - pow(1/r,7));  // r is really r/sigma
}

double Properties::lenJonesEnergy(double r, double c){
   return 4*c * (pow(1/r,12) - pow(1/r,6) + truncShift); // r is really r/sigma 
}

void Properties::calcEnergy(double r, double c){
   f_energy = f_energy + lenJonesEnergy(r,c); // sums the energy of the current 
}                                           // configuration
void Properties::calcVirial(double r, double c){      // sums the virial of current
   f_r = f_r + r * lenJonesForce(r,c);        // configuration
}

void Properties::updateNumDensity(double r){
   int val = r/delta_r;
   int index = 0; 
    
   if(r < boxLength/2){
      if(r > (val + 0.5) * delta_r){
         index = val + 1; 
      }
      else{
         index = val; // come back to check this... 
      } 
      num_density[index] = num_density[index] + 1;
   }
}
void Properties::populateCellArray(double x,double y, std::vector<std::vector<double>>* cellPositions){
   
   // defines the 8 images of the comparison particle's position
   (*cellPositions)[0][0] = x;             (*cellPositions)[0][1] = y + boxLength; 
   (*cellPositions)[1][0] = x;             (*cellPositions)[1][1] = y - boxLength; 
   (*cellPositions)[2][0] = x + boxLength; (*cellPositions)[2][1] = y; 
   (*cellPositions)[3][0] = x + boxLength; (*cellPositions)[3][1] = y + boxLength; 
   (*cellPositions)[4][0] = x + boxLength; (*cellPositions)[4][1] = y - boxLength; 
   (*cellPositions)[5][0] = x - boxLength; (*cellPositions)[5][1] = y; 
   (*cellPositions)[6][0] = x - boxLength; (*cellPositions)[6][1] = y + boxLength; 
   (*cellPositions)[7][0] = x - boxLength; (*cellPositions)[7][1] = y - boxLength; 

}


void Properties::calcPeriodicProp(std::vector<Particle>* particles, 
		                  std::ofstream* r_dist_file){ // this needs to be reorganized
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
   
   double LJ_constant = 0; 
   double r_dist = 0; 
   double force_tot = 0; 
   
   std::vector<std::vector<double>> cellPositions(9,std::vector<double>(2,0));
   
   f_energy = 0;   // make sure that the free energy previously calculated is reset
   f_r = 0;        // the free energy is only the energy that comes from the positions  
                   // within the configuration 
   for(int k = 0; k < n_particles; k++){ 
      curr_prt = (*particles)[k]; 
      
      x_curr = curr_prt.getX_Position();    // set current x,y position
      y_curr = curr_prt.getY_Position(); 
   
      for(int n = 0; n < n_particles; n++){ // each particle-particle interaction
	      
	 comp_prt =(*particles)[n]; 
         x_comp = comp_prt.getX_Position(); // set comparison x,y position
	 y_comp = comp_prt.getY_Position(); 

         // the particle cannot interact with itself
	 if(curr_prt.getIdentifier() != comp_prt.getIdentifier()){  
	    r_dist = radDistance(x_curr,x_comp,y_curr,y_comp);
	    updateNumDensity(r_dist);   
	 }
         
         if(curr_prt.getType() == comp_prt.getType()){
	    LJ_constant = LJ_par; 
	 }
	 else if(curr_prt.getType() != comp_prt.getType()){
	    LJ_constant = LJ_antipar; 
	 }

	 if(n > k){
	    if(r_dist > truncDist){
	       populateCellArray(x_comp,y_comp,&cellPositions);
               for(int z = 0; z < 8; z++){
               
	          x_comp = cellPositions[z][0]; // this is very inefficient but works... 
	          y_comp = cellPositions[z][1]; // creates the 8 cells surrounding 
                                                // the reference cell 
		  r_dist = radDistance(x_curr, x_comp, y_curr, y_comp);
		  for(int j = 0; j < 2; j++){
                     updateNumDensity(r_dist); 
		  }
		  if(r_dist < truncDist){
	             for(int j = 0; j < 2; j ++){
		        calcEnergy(r_dist,LJ_constant);  
		        calcVirial(r_dist,LJ_constant); 
		     }
                  } 
               }
            }
            else{
	       calcEnergy(r_dist,LJ_constant);  
	       calcVirial(r_dist,LJ_constant); 
	    }
	 }                          
      }
   }
   sum_Fdot_r.push_back(f_r); 
   sum_energy.push_back(f_energy);
}


double Properties::calcPressure(){
   double len = 0; 
   double avgEnergy = 0; 
   double redPressure = 0; 

   len = double(sum_Fdot_r.size()); // the pressure correction is 
                                    // added in the analysis code
   for(int k = 0; k < len; k++){
      avgEnergy = avgEnergy + sum_Fdot_r[k];  
   }
   avgEnergy = avgEnergy/len;
   redPressure = redDens * (redTemp + avgEnergy / (2 * n_particles));  
   return redPressure; 
}

double Properties::calcAvgEnergy(){
   double len = 0; 
   double avgEnergy = 0; 

   len = double(sum_energy.size()); 

   for(int k = 0; k < len; k++){
      avgEnergy = avgEnergy + sum_energy[k]; // sum all terms in the   
   }                                         // energy vector to compute
   avgEnergy = avgEnergy/len;                // the average energy
   return avgEnergy;                         // of the whole config
}

void Properties::writeProperties(){
   double len = 0; 

   std::ofstream virial_file; 
   std::ofstream energy_file; 
   std::ofstream n_dens_file; 

   virial_file.open("forces.txt");  // open each file that will be written to
   energy_file.open("energies.txt"); 
   n_dens_file.open("numDensity.txt");

   len = double(sum_Fdot_r.size()); // the force and energy vector are the same                                    
   for(int k = 0; k < len; k++){    // size hence are put into one for-loop
      virial_file << sum_Fdot_r[k] << " ";
      energy_file << sum_energy[k] << " "; 
   }
   virial_file.close(); 
   energy_file.close(); 

   len = double(num_density.size());
   for(int k = 0; k < len; k++){
      n_dens_file << num_density[k] << " ";
   }
   n_dens_file.close();             // close all files written to 
}

void Properties::initializeProperties(Parameters* p){

   boxLength = p->getBoxLength();      // assign private variables used in 
   n_particles = p->getNumParticles(); // class
   sigma = p->getSigma(); 
   redDens = p->getRedDens(); 
   redTemp = p->getRedTemp(); 

   LJ = p->getLenJones(); 
  
   LJ_par = p->getLJ_const_1(); 
   LJ_antipar = p->getLJ_const_2(); 

   delta_r = sigma/20; // this might not be the best way to define delta_r
   
   truncDist = 2.5;                 // really is 2.5 * sigma / sigma 
   truncShift = -1 * lenJonesEnergy(truncDist,.25); // shifts the cutoff to zero 
   num_density.resize((boxLength/2)/delta_r + 1);
}
