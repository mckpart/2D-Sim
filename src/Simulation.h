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

class Simulation {

  private:
    std::string yamlFile;

    Parameters param;
    Interaction interact;
    Boundary bound;
    Properties prop;

    std::vector<Particle> particles;

    KISSRNG randVal;

    int n_particles = 0;
    double red_temp = 0;

  public:
    Simulation(std::string yf);

    double boltzmannFactor(double delta_energy);

    void runSimulation();
    void setParticleParams();
    void writePositions(std::ofstream *pos_file);
    void testSimulation();
};
#endif
