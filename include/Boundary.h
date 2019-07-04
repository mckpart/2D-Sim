#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Particle.h"

class Boundary{

	public:

		bool rigidBoundary(Particle* particles, int index, int n_particles, double x_temp, double y_temp);  
};
#endif