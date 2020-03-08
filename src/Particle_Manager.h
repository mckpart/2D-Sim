#ifndef PARTICLE_MANAGER_H
#define PARTICLE_MANAGER_H

// can probably remove this
#include <iostream>
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
    std::ofstream force_file;

  public:
    Particle_Manager();
    Particle_Manager(std::string yaml_file);
    void init_particle_params(std::vector<Particle> *particles,
                              KISSRNG *randVal);

    void assign_forces(double r, double LJ_const, Particle *curr,
                       Particle *ref);
    void reset_forces(std::vector<Particle> *particles);
};
#endif

