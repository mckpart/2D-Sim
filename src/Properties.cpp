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

double Properties::lenJonesEnergy(double r, double a){  
   return 4*a * (pow(1/r,12) - pow(1/r,6) + truncShift); // r is really r/sigma 
}

// might be worth making a new class with different
// types of energies and their associated forces
// NOTE: THE BINDING AFFINITY SHOULD NOT BE ATTACHED TO THE 
// WCA POTENTIAL SINCE THE WCA POTENTIAL IS SERVING THE 
// PURPOSE OF A SOFT DISK INTERACTION
double Properties::WCA_force(double r){
   double val = 0; 
   if(r <= pow(2.0,1.0/6.0)){
      val = 24/sigma * (2*pow(1/r,13) - pow(1/r,7)); 
   }
   return val; // if the particles are further than sigma * 2^(1/6)
}              // apart then there is no force/potential energy

double Properties::WCA_energy(double r){ // a behaves as the
   double val = 0;                       // binding affinity
   if(r<= pow(2.0,1.0/6.0)){
      val = 4*(pow(1/r,12)-pow(1/r,6) + .25); // WCA is more of a piecewise
   }                                          // potential 
   return val; 
}

double Properties::simple_spring_force(double r, double a){
   return a*red_temp*k_spring*(r-rest_L)*
	  exp(-k_spring/2*pow(r-rest_L,2))*
          (k_spring/2*pow(r-rest_L,2)-1); 
}

double Properties::simple_spring_energy(double r, double a){
   return a*red_temp*k_spring/2*pow(r-rest_L,2)
	   *exp(-k_spring/2.0*pow(r-rest_L,2.0)); // add the spring constant later
}

void Properties::calcEnergy(double r, double a){
   double val = 0;
   switch(interact_type){
      case 1: val = lenJonesEnergy(r,a); 
	      break; 
      case 2: val = WCA_energy(r); 
	      break;
      case 3: val = WCA_energy(r) + 
	      simple_spring_energy(r,a); 
	      break;
   }
   f_energy = f_energy + val; // calculates energy of current 
}                             // configuration
void Properties::calcVirial(double r, double a){ // sums the virial of current
//   f_r = f_r + r * lenJonesForce(r,a);           // configuration
   double val = 0; 
   switch(interact_type){
      case 1: val = lenJonesForce(r,a);
	      break;
      case 2: val = WCA_force(r); 
	      break; 
      case 3: val = WCA_energy(r) + 
	      simple_spring_force(r,a); 
	      break;
   }
   f_r = f_r + r * val; 
}

void Properties::updateNumDensity(double r, int ID){ // this could always return the index to increase... 
   int val = r/delta_r;
   int index = 0; 
    
   if(r < boxLength/2){
      if(r > (val + 0.5) * delta_r){
         index = val + 1; 
      }
      else{
         index = val;
      } 

      switch(ID){
         case 0: num_density[index] = num_density[index] + 1; 
                 break;
	 case 1: par_num_density[index] = par_num_density[index] + 1; 
		 break;
	 case 2: antp_num_density[index] = antp_num_density[index] + 1; 
		 break;
      }
   }
}

void Properties::calc_xy_dens(double x, double y, int ID){
   int ind_1 = 0; 
   int ind_2 = 0; 

   ind_1 = (x + sqrt(2)*boxLength/2)/cell_L; // cell_L = delta x = delta y 
   ind_2 = (y + sqrt(2)*boxLength/2)/cell_L;
//   std::cout << "ind_1 " << ind_1 << " ind_2 " << ind_2 << std::endl;   
//   std::cout << "x " << x << " y " << y << std::endl;
   if(fabs(x) < sqrt(2)*boxLength/2 - cell_L && fabs(y) < sqrt(2)*boxLength/2 - cell_L){
      if(x > (ind_1 + 0.5) * cell_L - sqrt(2)*boxLength/2){
         ++ind_1; 
      }
      if(y > (ind_2 + 0.5) * cell_L - sqrt(2)*boxLength/2){
         ++ind_2; 
      }
      
      switch(ID){
         case 0: xy_num_density[ind_1][ind_2] = xy_num_density[ind_1][ind_2] + 1; // increment the value here. 
	         break;
	 case 1: par_xy_density[ind_1][ind_2] = par_xy_density[ind_1][ind_2] + 1; // increment the value here. 
	         break;
	 case 2: antp_xy_density[ind_1][ind_2] = antp_xy_density[ind_1][ind_2] + 1; // increment the value here. 
                 break;
      }
   }
}

void Properties::calcNonPerProp(std::vector<Particle>* particles){

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

   int temp = 0;  
   
//   std::vector<std::vector<double>> cellPositions(9,std::vector<double>(2,0));
   
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
	    updateNumDensity(r_dist,0);   // updates overall number density
	    calc_xy_dens(x_comp - x_curr, y_comp - y_curr,0); 

	    if(curr_prt.getType() == comp_prt.getType()){ // interaction of 
	       LJ_constant = LJ_par;                // parallel microtubules
	       updateNumDensity(r_dist,1);          // updates number density 
	       calc_xy_dens(x_comp-x_curr,y_comp-y_curr,1); //for parallel
	    }                                       // interactions
	    else if(curr_prt.getType() != comp_prt.getType()){ // interaction
	       LJ_constant = LJ_antipar;       // of antiparallel microtubules
	       updateNumDensity(r_dist,2);     // updates number density for
	       calc_xy_dens(x_comp-x_curr,y_comp-y_curr,2); // antiparallel
	    }                                               // interactions
	 }
	 if(n > k && r_dist < truncDist){
//            if(r_dist < truncDist){
	       calcEnergy(r_dist,LJ_constant);  
	       calcVirial(r_dist,LJ_constant); 
//	    }
	 }                          
      }
   }
   sum_Fdot_r.push_back(f_r); 
   sum_energy.push_back(f_energy);
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

void Properties::calcPeriodicProp(std::vector<Particle>* particles){ // this needs to be reorganized
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

   int temp = 0;  
   
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
	    updateNumDensity(r_dist,0);   // updates overall number density
	    calc_xy_dens(x_comp - x_curr, y_comp - y_curr,0); 

	    if(curr_prt.getType() == comp_prt.getType()){ // interaction of 
	       LJ_constant = LJ_par;                // parallel microtubules
	       updateNumDensity(r_dist,1);          // updates number density 
	       calc_xy_dens(x_comp-x_curr,y_comp-y_curr,1); //for parallel
	    }                                       // interactions
	    else if(curr_prt.getType() != comp_prt.getType()){ // interaction
	       LJ_constant = LJ_antipar;       // of antiparallel microtubules
	       updateNumDensity(r_dist,2);     // updates number density for
	       calc_xy_dens(x_comp-x_curr,y_comp-y_curr,2); // antiparallel
	    }                                               // interactions
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
                     updateNumDensity(r_dist,0); 
		     calc_xy_dens(x_comp-x_curr,y_comp-y_curr,0);
		     
		     if(curr_prt.getType() == comp_prt.getType()){
	                updateNumDensity(r_dist,1);  
	                calc_xy_dens(x_comp-x_curr,y_comp-y_curr,1); 
		     }
	             else if(curr_prt.getType() != comp_prt.getType()){
	                updateNumDensity(r_dist,2); 
	                calc_xy_dens(x_comp-x_curr,y_comp-y_curr,2); 
		     }	  
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

// this calculation has been moved to external analysis code
// and can probably be removed. 
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
   redPressure = redDens * (red_temp + avgEnergy / (2 * n_particles));  
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
   std::ofstream par_dens_file; 
   std::ofstream antp_dens_file;
   
   std::ofstream xy_dens_file;  
   std::ofstream par_xy_file; 
   std::ofstream antp_xy_file; 

   virial_file.open("forces.txt"); // open each file that will be written to
   energy_file.open("energies.txt"); 
   
   n_dens_file.open("numDensity.txt");
   par_dens_file.open("par_numDensity.txt");
   antp_dens_file.open("antp_numDensity.txt");
   
   xy_dens_file.open("xy_numDensity.txt"); 
   par_xy_file.open("par_xy_numDensity.txt");
   antp_xy_file.open("antp_xy_numDensity.txt"); 

   len = double(sum_Fdot_r.size()); // the force and energy vector are the same                                    
   for(int k = 0; k < len; k++){    // size hence are put into one for-loop
      virial_file << sum_Fdot_r[k] << " ";
      energy_file << sum_energy[k] << " "; 
   }
   virial_file.close(); 
   energy_file.close(); 

   len = double(num_density.size());
   for(int k = 0; k < len; k++){            // all of the number density 
      n_dens_file << num_density[k] << " "; // vectors are the same length 
      antp_dens_file << antp_num_density[k] << " "; 
      par_dens_file << par_num_density[k] << " "; 
   }
   n_dens_file.close();             // close all files written to 
   par_dens_file.close(); 
   antp_dens_file.close(); 
   len = double(xy_num_density.size()); 
   for(int k = 0; k < len; ++k){
      for(int n = 0; n < len; ++n){
         xy_dens_file << xy_num_density[k][n]  << " "; 
         par_xy_file  << par_xy_density[k][n]  << " "; 
         antp_xy_file << antp_xy_density[k][n] << " "; 
      }
   }
   xy_dens_file.close();
   par_xy_file.close(); 
   antp_xy_file.close();  
}

void Properties::truncation_dist(){
   switch(interact_type){
//      case 1:  
	       
//              break;
//      case 2: break;  
      case 3: truncDist = boxLength/2.0;
	      break;
      default:  truncDist = 2.5;
                truncShift = -1*(pow(1/truncDist,12)
			        -pow(1/truncDist,6));
                break;
   }
}

void Properties::initializeProperties(Parameters* p){ // maybe make this into a 
                                                      // constructor
   boxLength = p->getBoxLength();      // assign private variables used in 
   n_particles = p->getNumParticles(); // class
   sigma = p->getSigma(); 
   redDens = p->getRedDens(); 
   red_temp = p->getRedTemp(); 

   k_spring = p->getSprConst(); 

   interact_type = p->getInteract_Type(); 
   rest_L = p->getRestLength(); 

   LJ_par = p->getLJ_const_1(); 
   LJ_antipar = p->getLJ_const_2(); 

   delta_r = sigma/20; // this might not be the best way to define delta_r
   cell_L = sigma/20; 
//   truncDist = 2.5;                 // really is 2.5 * sigma / sigma 
//   truncShift = -1 * lenJonesEnergy(truncDist,.25); // shifts the cutoff to zero 
   truncation_dist();   
   // define the various RDF vectors (dependent upon r)
   num_density.resize((boxLength/2)/delta_r + 1);
   par_num_density.resize((boxLength/2)/delta_r + 1);
   antp_num_density.resize((boxLength/2)/delta_r + 1);
   
   int val = sqrt(2)*boxLength/cell_L + 1;
   std::cout << "the size of the vector is " << val << std::endl; 
   xy_num_density.resize(val,std::vector<double>(val,0)); 
   par_xy_density.resize(val,std::vector<double>(val,0)); 
   antp_xy_density.resize(val,std::vector<double>(val,0)); 
}
