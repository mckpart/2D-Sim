#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Particle.h"

class Boundary{

	public:

		void periodicBoundary(Particle* particles, int index, int n_particles);

		bool rigidBoundary(Particle* particles, int index, int n_particles);  
};
#endif