#include <iostream>
#include "Simulation.h"


double boltzmannFactor(double energy){
   return exp(-1 * energy); 
}

Simulation::Simulation(std::string yf){
   yamlFile = yf; 

   param.initializeParameters(yamlFile);    // initialize the parameters
   interact.initializeInteraction(yamlFile);// for the simulation

   n_particles = param.getNumParticles();   // initialize vector of
   particles.resize(n_particles);           // particles and set particle
   setParticleParams();	                    // parameters
}

void Simulation::runSimulation(){
	
   Particle prt; 

   int n_updates = 0;

   double n_rejects = 0;  
   double perc_rej  = 0; 

   double x_trial = 0; 
   double y_trial = 0; 

   double delta_energy = 0; 
   double total_prob   = 0; 

   bool accept = 0; 

   std::ofstream pos_file; 
   pos_file.open("positions.txt"); 

   KISSRNG randVal; 
   randVal.InitCold(param.getSeed()); 	// warms the RNG

   n_updates 	= param.getUpdates(); 

   interact.initialPosition(&particles,n_particles,randVal); 

   for(int n = 0; n < n_updates; n++){
      for(int k = 0; k < n_particles; k++){
	     prt = particles[k]; 

         x_trial = prt.x_trial(randVal.RandomUniformDbl());  
         y_trial = prt.y_trial(randVal.RandomUniformDbl()); 

         prt.setX_TrialPos(x_trial); 	// generate and set the
         prt.setY_TrialPos(y_trial); 	// x,y trial position

         particles[k] = prt;  

         accept = 1; 		

         if(param.getRigidBC() == 1){        // run sim with hard boundaries
            accept = bound.rigidBoundary(&particles,k);  
         }
         else if(param.getPeriodicBC() == 1){     // run simulation with 
            bound.periodicBoundary(&particles,k); // periodic boundaries

            prt = particles[k]; 

            x_trial = prt.getX_TrialPos();   // updates trial position in function
            y_trial = prt.getY_TrialPos();   // then particle - particle 
         }                                   // interactions are checked
         else if(param.getExtWell() == 1){
            delta_energy = bound.externalWell(&particles,k); 
         }

      /*  - RUNS DIFFERENT TYPES OF PARTICLE-PARTICLE INTERACTIONS
          - HARD DISKS IS A 0 - 1 PROBABILITY THUS A DELTA ENERGY 
            IS NOT RETURNED
          - THE CHANGE IN ENERGY IS RETURNED FROM LENJONES AND WCA
            POTENTIAL 
          - THIS TOTAL CHANGE IS SENT INTO THE BOLTZMANN FACTOR 
            FUNCTION TO CALCULATE THE TOTAL PROBABILITY OF
            ACCEPTING THE TRIAL MOVE. IF THE TOTAL CHANGE < 0, THE
            MOVE IS ACCEPTED. ELSE A RANDOM NUMBER IS GENERATED TO 
            DETERMINE WHETHER THE MOVE IS TO BE ACCEPTED
      */		  
				  
         if(param.getHardDisk() == 1 && accept == 1){
            accept = interact.hardDisks(&particles,k,n_particles);
         }
         else if(param.getLenJones() == 1 && accept == 1){
            delta_energy = delta_energy +  
                           interact.lennardJones(&particles,k,n_particles); 
         }
         else if(param.getWCA() == 1 && accept == 1){
            delta_energy = delta_energy +  
                           interact.WCApotential(&particles,k,n_particles); 
         }

         if(accept == 1 && delta_energy > 0){
            total_prob = boltzmannFactor(delta_energy); // compute acceptance probability

            if(randVal.RandomUniformDbl() < total_prob){
               accept = 1; 
            }
            else{
               accept = 0; 
            }
         }

         if(accept == 1){
            prt.setX_Position(x_trial);	// if trial move is accepted, update the 
            prt.setY_Position(y_trial); // position of the current particle
            particles[k] = prt;         // and put back into array

            pos_file << prt.getX_Position() << " "; // read updated positions into  
            pos_file << prt.getY_Position() << " "; // file
         }
         else{
            n_rejects++;   // keeps count of total moves rejected
            k = k - 1;     // reject trial move - will generate new trial 
         }                 // position for the same particle 	
      }

      pos_file << std::endl;   // starts new line in position file
   }

   perc_rej = n_rejects / (n_rejects + (n_updates * n_particles)) * 100.0; 
   std::cout << perc_rej << "% of the moves were rejected." << std::endl; 
}

// THIS IS THE NEXT PIECE TO BE ALTERED ////

void Simulation::setParticleParams(){	

   int num_part1 = 0; 
   int num_part2 = 0; 

   double radius_1 = 0; 
   double radius_2 = 0; 

   YAML::Node node = YAML::LoadFile(yamlFile);

   num_part1 = node["type1_Particles"].as<int>(); 
   num_part2 = node["type2_Particles"].as<int>(); 

   radius_1 = node["particleRadius_1"].as<double>(); 
   radius_2 = node["particleRadius_2"].as<double>(); 

   Particle prt; 
   prt.setStepWeight(node["weight"].as<double>());

   prt.setRadius(radius_1);
   prt.setType(1); 

   for(int k = 0; k < num_part1; k++){	
      prt.setIdentifier(k);   // initialize all particles of type 1
      particles[k] = prt;     // with specific radius and assign 
   }                          // unique ID

   prt.setRadius(radius_2); 
   prt.setType(2); 

   for(int k = 0; k < num_part2; k++ ){  // initialize all particles
      prt.setIdentifier(k + num_part1);  // of type 2 with specific 
      particles[k + num_part1] = prt;    // radius and unique ID
   }
}
