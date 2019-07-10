#include <iostream>
#include <string>

#include "Simulation.h"


int main(int argc, char *argv[]){

   if(argc > 1){	
      std::string yamlFile;

      Simulation sim(argv[1]); 	// initialize the simulation
      sim.runSimulation();      // run the simulation
   }
   else{
      std::cout << "ERROR: NO .YAML FILE" << std::endl; 
   }
}
