#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <cmath>
#include <fstream>
#include <vector>

#include "Parameters.h"
#include "Particle.h"

class Properties {

  private:
    double f_energy = 0;
    double f_r = 0;

    int force_num = 0;
    std::vector<double> avg_force{std::vector<double>(2, 0)};
    std::ofstream avg_force_particle;

    std::vector<double> sum_Fdot_r;
    std::vector<double> sum_energy;

    std::vector<double> num_density;
    std::vector<double> par_num_density;
    std::vector<double> antp_num_density;

    // 2D arrays used for calculating the PCFs
    std::vector<std::vector<double>> xy_num_density;
    std::vector<std::vector<double>> par_xy_density;
    std::vector<std::vector<double>> antp_xy_density;

    // this is the average force per particle - holds x,y coordinate
    //    std::vector<std::vector<double>> avg_force;

    double delta_r = 0;
    double cell_L = 0;

    double sigma = 0;
    double truncDist = 0;
    double truncShift = 0;
    int interact_type = 0;

    double k_spring = 0;
    double rest_L = 0;

    double a_ref = 0;
    double a_mult = 0;

    double boxLength = 0;
    int n_particles = 0;

    double redDens = 0;
    double red_temp = 0;

  public:
    void initializeProperties(Parameters *p);
    double truncation_dist();

    void populateCellArray(double x, double y,
                           std::vector<std::vector<double>> *cellPositions);
    double lenJonesEnergy(double r, double a);
    double lenJonesForce(double r, double a);

    double WCA_energy(double r);
    double WCA_force(double r);
    double simple_spring_energy(double r, double a);
    double simple_spring_force(double r, double a);

    void updateNumDensity(double r, int ID);
    void calc_xy_dens(double x, double y, int ID);

    void calcPeriodicProp(std::vector<Particle> *particles);
    void calcNonPerProp(std::vector<Particle> *particles);
    double calcEnergy(double r, double c);
    double calcVirial(double x, double y, double r, double a);

    double radDistance(double x1, double x2, double y1, double y2);

    double calcPressure();
    double calcAvgEnergy();

    // remember, these are being moved to the particle manager class eventually
    void calc_force_vec(double x, double y, double r,
                        std::vector<double> *F_vec);
    void avg_force_vec(std::vector<std::vector<double>> *F);

    void writeProperties();
    void writeAvgForces();

    void open_files();
    void close_files();
};
#endif
