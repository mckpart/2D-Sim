#ifndef SIMULATION_H
#define SIMULATION_H

#include <fstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "Boundary.h"
#include "Interaction.h"
#include "Parameters.h"
#include "Particle.h"
#include "Properties.h"
#include "System_Manager.h"

class Simulation {

  private:
    std::string yaml_file;

    Parameters param;
    Interaction interact;
    Boundary bound;
    Properties prop;
    System_Manager sim_manage;
    std::vector<Particle> particles;

    KISSRNG randVal;

    int n_particles = 0;
    double red_temp = 0;

    double delta_energy = 0;
    std::vector<double> trial_position;


  public:
    Simulation(std::string yf);

    double boltzmannFactor(double delta_energy);

    void runSimulation();
    // recall that this function is going to get moved to the particle manager
    // class
    void setParticleParams();
    void writePositions(std::ofstream *pos_file);
    void testSimulation();

    // private functions that I want to have access to the properties of the
    // class
    void periodic_pos_trial(Particle *p);
    bool nonperiodic_pos_trial(Particle *p, bool accept);
    void init_configuration();

    //    void calcNonPerProp(std::vector<Particle> particles);
    // have this commented out since the particles are currently part of the sim
    // class
    void calcNonPerProp();
};
#endif
