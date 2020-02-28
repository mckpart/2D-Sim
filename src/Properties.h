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

    // these are included in the sysmanager class... go through and update the
    // appropriate functions that use these
    std::vector<double> sum_Fdot_r;
    std::vector<double> sum_energy;

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

    double calcEnergy(double r, double c);
    double calcVirial(double x, double y, double r, double a);

    double radDistance(double x1, double x2, double y1, double y2);

    double calcPressure();
    double calcAvgEnergy();

    // remember, these are being moved to the particle manager class eventually
    void calc_force_vec(double x, double y, double r,
                        std::vector<double> *F_vec);
    void avg_force_vec(std::vector<std::vector<double>> *F);

    void writeProperties(std::vector<double> *sum_energy,
                         std::vector<double> *sum_virial,
                         std::vector<double> *nd, std::vector<double> *par_nd,
                         std::vector<double> *antp_nd,
                         std::vector<std::vector<double>> *xy,
                         std::vector<std::vector<double>> *par_xy,
                         std::vector<std::vector<double>> *antp_xy);
    void writeForces(std::vector<Particle> *particles);

    void open_files();
    void close_files();
};
#endif
