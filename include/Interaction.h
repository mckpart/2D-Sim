#ifndef INTERACTION_H
#define INTERACTION_H

#include "Particle.h"
#include "kiss.h"

class Interaction{ 
	
	public: 

		bool rigidCollisions(Particle* particles, int index, int n_particles, double x_temp, double y_temp); 
		void initialPosition(Particle* particles, int n_particles, KISSRNG randVal);
};
#endif