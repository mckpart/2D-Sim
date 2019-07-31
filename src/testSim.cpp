#include <iostream>
#include "Simulation.h"
#include <fstream>

int main(){
   Properties prop; 
   std::vector<Particle> particles; 
   Particle prt; 

   std::ofstream rfile; 
   rfile.open("radialDistance.txt");
   if(rfile.is_open()){
      std::cout << "good" << std::endl; 
   } 

   double x = 0; 
   double y = 0;
   for(int k = 0; k < 5; k++){
      switch(k){
         case 0: 
            x = 1; y = 1; break;
	 case 1:
	    x = 1; y = -1; break; 
         case 2:
	    x = -1; y = 1; break;
	 case 3: 
	    x = -1; y = -1; break; 
         case 4: 
	    x = 0; y = 0; break; 
      }
      prt.setX_Position(x); 
      prt.setY_Position(y); 
      particles.push_back(prt); 
   }
   std::cout << "hmmm " << std::endl;
   prop.calcPeriodicProp(&particles,&rfile);  
   rfile.close();  
   return 0; 
}
