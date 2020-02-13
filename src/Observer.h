#ifndef PARTICLE_OBS_H
#define PARTICLE_OBS_H

#include <vector>

#include "Properties.h"

class Particle_obs {
  private:
    // vector of strings that holds the different combinations for  .. this may
    // be a good unit test later...
    //
    // current total force (2 by 1 vector)
  public:
    // update total force on the particle
    void update_forces();
};
#endif
