#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <yaml-cpp/yaml.h>

class Parameters {

  private:
    double redDensity = 0;
    double redTemp = 0;
    double sigma = 0;
    double sigma_par = 0;
    double boxLength = 0;
    double rest_L = 0;
    double k_spring = 0;

    double a_ref = 0;
    double a_mult = 0;

    int eq_sweep = 0;
    int d_interval = 0;

    int n_updates = 0;
    int n_particles = 0;
    long seed = 0;

    double ext_well_d = 0;

    int init_type = 0;
    int interact_type = 0;
    int bound_type = 0;

  public:
    // may be worth adding a default constructor that reads in the yaml file
    void initializeParameters(std::string yamlFile);

    int getUpdates();
    int getNumParticles();
    long getSeed();

    int getInit_Type();
    int getInteract_Type();
    int getBound_Type();

    int getEq_sweep();
    int getData_interval();

    double getRefAffinity();
    double getAffinityMult();

    double getExtWellDepth();

    double getSprConst();
    double getRestLength();
    double getRedDens();
    double getRedTemp();
    double getSigma();
    double getBoxLength();

    // NOTE: THE SETTERS ARE NOT NECESSARY
    // GIVEN THAT THESE PARAMETERS ARE CONSTANT
    // THROUGHOUT THE SIMULATION
};
#endif
