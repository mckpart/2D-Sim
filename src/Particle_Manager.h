#ifndef PARTICLE_MANAGER_H
#define PARTICLE_MANAGER_H

// can probably remove this
#include <vector>

#include "Parameters.h"
#include "Particle.h"
#include "Properties.h"
#include "kiss.h"

class Particle_Manager {
    // add the particle vector to this class - will need to reorganize sim
  private:
    Parameters param;
    Properties prop;

  public:
    Particle_Manager();
    Particle_Manager(std::string yaml_file);
    void init_particle_params(std::vector<Particle> *particles,
                              KISSRNG *randVal);
};
#endif

