#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include "Particle.h"

class Boundary{

	public:

		void periodicBoundary(std::vector<Particle>* particles, int index);
		bool rigidBoundary(std::vector<Particle>* particles, int index);
		double externalWell(std::vector<Particle>* particles, int index);  
};
#endif