#ifndef PARTICLE_H
#define PARTICLE_H

class Particle{

private:
   int type = 0;	// may be worth changing type to a string later on
   int identifier = 0; 

   double radius = 0; 
   double x_position = 0; 
   double y_position = 0; 

   double x_trialPos = 0; 
   double y_trialPos = 0;  

   double stepWeight = 0; 

public: 

   // Particle(); 
   // // Particle(int t, double rad, double w); 

   int getType(); 
   int getIdentifier(); 

   double getRadius(); 
   double getX_Position(); 
   double getY_Position();
   double getX_TrialPos(); 
   double getY_TrialPos(); 

   double getStepWeight(); 	

   void setType(int t); 
   void setIdentifier(int id); 

   void setRadius(double rad); 
   void setX_Position(double x); 
   void setY_Position(double y); 
   void setX_TrialPos(double x); 
   void setY_TrialPos(double y); 

   void setStepWeight(double w); 

   double x_trial(double randVal); 
   double y_trial(double randVal); 

};
#endif
