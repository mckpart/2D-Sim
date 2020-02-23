#ifndef SYSTEM_MANAGER_H
#define SYSTEM MANAGER_H

#include "Particle.h"
#include "Properties.h"

#include <vector>

// it may be better to move the 'main loop' for this into the simulation class.
// The observer class can then organize the data (i.e. the total energy config
// of the system, the current RDF/PCF, handle updating the forces for each
// particle, etc.

class System_Manager {
//  private:
//    // vector of strings that holds the different combinations for  .. this
//    may
//    // be a good unit test later...
//    //
//    //
//    // maybe keep private objects here
//    // current total force (2 by 1 vector)
//    Properties prop;
//    // just moved over from the properties class
//
//    // these are values used for initializing the rdf & pcf vectors
//    double delta_r = 0;
//    double cell_L = 0;
//    // arrays used for calculating the RDF
//    std::vector<double> num_density;
//    std::vector<double> par_num_density;
//    std::vector<double> antp_num_density;
//
//    // 2D arrays used for calculating the PCFs
//    std::vector<std::vector<double>> xy_num_density;
//    std::vector<std::vector<double>> par_xy_density;
//    std::vector<std::vector<double>> antp_xy_density;
//
//  public:
//    // here is the default constructor
//    System_Manager();
//    // update total force on the particle
//    // move the update forces function into the particle manager class
//    //     void update_forces();
//
//    void Properties::updateNumDensity(double r, int ID);
//    void Properties::calc_xy_dens(double x, double y, int ID);
//    // recall that this routine should really be moved to the sim class
//    void calcPerProp(std::vector<Particle> *particles);
};
#endif
