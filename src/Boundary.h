#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include "Particle.h"

class Boundary{

	public:

		void periodicBoundary(std::vector<Particle>* particles, int index, int n_particles);

		bool rigidBoundary(std::vector<Particle>* particles, int index, int n_particles);  
};
#endif