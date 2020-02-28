#include <iostream>
#include "Properties.h" 

// returns the characteristic length between two particles
double Properties::radDistance(double x1, double x2, 
		               double y1, double y2){
   return sqrt(pow(x2-x1,2) + pow(y2-y1,2)) / sigma;  
}

double Properties::lenJonesForce(double r, double c){
    // r is really r/sigma
    return 24 * c / sigma * (2 * pow(1 / r, 13) - pow(1 / r, 7));
}
// Note: r is the characteristic length r/sigma
double Properties::lenJonesEnergy(double r, double a) {
    return 4 * a * (pow(1 / r, 12) - pow(1 / r, 6) + truncShift);
}

// might be worth making a new class with different
// types of energies and their associated forces
// NOTE: THE BINDING AFFINITY SHOULD NOT BE ATTACHED TO THE 
// WCA POTENTIAL SINCE THE WCA POTENTIAL IS SERVING THE 
// PURPOSE OF A SOFT DISK INTERACTION
double Properties::WCA_force(double r) {
    double val = 0;
    // as usual... r = r/sigma
    // if r > 2^1/6 then there is no force/potential energy
    if (r <= pow(2.0, 1.0 / 6.0)) {
        val = 24 / sigma * (2 * pow(1 / r, 13) - pow(1 / r, 7));
    }
    return val; // if the particles are further than sigma * 2^(1/6)
}

// the WCA energy is a piecewise potential
double Properties::WCA_energy(double r) { // a behaves as the
    double val = 0;                       // binding affinity
    if (r <= pow(2.0, 1.0 / 6.0)) {
        val = 4 * (pow(1 / r, 12) - pow(1 / r, 6) + .25);
    }
    return val;
}

double Properties::simple_spring_force(double r, double a){
    return a * red_temp * k_spring * (r - rest_L) *
           exp(-.5 * k_spring * pow(r - rest_L, 2)) *
           (k_spring * .5 * pow(r - rest_L, 2) - 1);
}

double Properties::simple_spring_energy(double r, double a) {
    return a * .5 * red_temp * k_spring * pow(r - rest_L, 2) *
           exp(-.5 * k_spring * pow(r - rest_L, 2.0));
}

// calculates the total energy of current configuration
double Properties::calcEnergy(double r, double a) {
    double val = 0;
    switch (interact_type) {
    case 1:
        val = lenJonesEnergy(r, a);
        break;
    case 2:
        val = WCA_energy(r);
        break;
    case 3:
        val = WCA_energy(r) + simple_spring_energy(r, a);
        break;
    }
    return val;
}

// sums the total virial of the current configuration
double Properties::calcVirial(double x, double y, double r, double a) {
    double val = 0;
    switch (interact_type) {
    case 1:
        val = lenJonesForce(r, a);
        break;
    case 2:
        val = WCA_force(r);
        break;
    case 3:
        val = WCA_force(r) + simple_spring_force(r, a);
        break;
    }

    // I temporarily change this to compute the avg virial per particle
    avg_force[0] = (x * val + avg_force[0] * force_num) / (1 + force_num);
    avg_force[1] = (y * val + avg_force[1] * force_num) / (1 + force_num);

    return r * val;
}

// this is for computing the non periodic properties... could probably also be
// moved to the simulation class
void Properties::populateCellArray(
    double x, double y, std::vector<std::vector<double>> *cellPositions) {

    // defines the 8 images of the comparison particle's position
    (*cellPositions)[0][0] = x;
    (*cellPositions)[0][1] = y + boxLength;
    (*cellPositions)[1][0] = x;
    (*cellPositions)[1][1] = y - boxLength;
    (*cellPositions)[2][0] = x + boxLength;
    (*cellPositions)[2][1] = y;
    (*cellPositions)[3][0] = x + boxLength;
    (*cellPositions)[3][1] = y + boxLength;
    (*cellPositions)[4][0] = x + boxLength;
    (*cellPositions)[4][1] = y - boxLength;
    (*cellPositions)[5][0] = x - boxLength;
    (*cellPositions)[5][1] = y;
    (*cellPositions)[6][0] = x - boxLength;
    (*cellPositions)[6][1] = y + boxLength;
    (*cellPositions)[7][0] = x - boxLength;
    (*cellPositions)[7][1] = y - boxLength;
}

// try to find a way to combine this calculation with the virial calculation
//
// JUST INCORORATED THE VECTOR INTO THIS FUNCTION. IF BAD, REMOVE!!!!!!!!!
void Properties::calc_force_vec(double x, double y, double r,
                                std::vector<double> *F_vec) {
    // I should NOT be forcing the interaction type to be 1 here.. bad practice
    interact_type = 1;
    //    std::vector<double> f(2, 0);
    double FF = 0;
    switch (interact_type) {
    case 1:
        FF = lenJonesForce(r, 1);
        break;
    case 2:
        FF = WCA_force(r);
        break;
    case 3:
        FF = WCA_force(r) + simple_spring_force(r, 1);
        break;
    }
    // here x = (x2 - x1) and y = (y2 - y1). It represents the vector pointing
    // from the reference particle to the comparison particle... thus a negative
    // sign is needed to appropriately orient the force vector
    (*F_vec)[0] = -x / r * FF + (*F_vec)[0];
    (*F_vec)[1] = -y / r * FF + (*F_vec)[1];
    // this currently matches for LJ - compute
}

// this routine is for calculating the total
void Properties::avg_force_vec(std::vector<std::vector<double>> *F) {
    std::vector<double> f(2);
    double fx = 0;
    double fy = 0;
    for (int i = 0; i < F->size(); ++i) {
        fx = fx + (*F)[i][0];
        fy = fy + (*F)[i][1];
    }
    std::cout << "in properties... fx = " << fx << " fy = " << fy << std::endl;
}

// this calculation has been moved to external analysis code
// and can probably be removed. 
double Properties::calcPressure() {
    double avgEnergy = 0;
    double redPressure = 0;

    double len = double(sum_Fdot_r.size()); // the pressure correction is
                                            // added in the analysis code
    for (int k = 0; k < len; k++) {
        avgEnergy = avgEnergy + sum_Fdot_r[k];
    }
    avgEnergy = avgEnergy / len;
    redPressure = redDens * (red_temp + avgEnergy / (2 * n_particles));
    return redPressure;
}

double Properties::calcAvgEnergy() {
    double len = 0;
    double avgEnergy = 0;

    len = double(sum_energy.size());

    for (int k = 0; k < len; k++) {
        avgEnergy = avgEnergy + sum_energy[k]; // sum all terms in the
    }                                          // energy vector to compute
    avgEnergy = avgEnergy / len;               // the average energy
    return avgEnergy;                          // of the whole config
}

void Properties::writeAvgForces() {
    for (int k = 0; k < 2; ++k) {
        avg_force_particle << avg_force[k] << " ";
    }
}

// this is probably a fine place to leave this function, but give it thought
void Properties::writeProperties(std::vector<double> *sum_energy,
                                 std::vector<double> *sum_virial,
                                 std::vector<double> *nd,
                                 std::vector<double> *par_nd,
                                 std::vector<double> *antp_nd,
                                 std::vector<std::vector<double>> *xy,
                                 std::vector<std::vector<double>> *par_xy,
                                 std::vector<std::vector<double>> *antp_xy) {
    double len = 0;

    std::ofstream virial_file;
    std::ofstream energy_file;

    std::ofstream n_dens_file;
    std::ofstream par_dens_file;
    std::ofstream antp_dens_file;

    std::ofstream xy_dens_file;
    std::ofstream par_xy_file;
    std::ofstream antp_xy_file;

    virial_file.open("forces.txt"); // open each file that will be written to
    energy_file.open("energies.txt");

    n_dens_file.open("numDensity.txt");
    par_dens_file.open("par_numDensity.txt");
    antp_dens_file.open("antp_numDensity.txt");

    xy_dens_file.open("xy_numDensity.txt");
    par_xy_file.open("par_xy_numDensity.txt");
    antp_xy_file.open("antp_xy_numDensity.txt");

    // the force and energy vector are the same size hence are put into one
    // for-loop
    len = double(sum_energy->size());
    for (int k = 0; k < len; k++) {
        virial_file << (*sum_virial)[k] << " ";
        energy_file << (*sum_energy)[k] << " ";
    }

    virial_file.close();
    energy_file.close();

    // len should be calculated ahead of time to avoid pushbacks in the sum
    // energy vector
    len = double(nd->size());
    for (int k = 0; k < len; k++) {           // all of the number density
        n_dens_file << (*nd)[k] << " ";       // vectors are the same length
        antp_dens_file << (*antp_nd)[k] << " ";
        par_dens_file << (*par_nd)[k] << " ";
    }
    n_dens_file.close(); // close all files written to
    par_dens_file.close();
    antp_dens_file.close();

    len = double(xy->size());
    for (int k = 0; k < len; ++k) {
        for (int n = 0; n < len; ++n) {
            xy_dens_file << (*xy)[k][n] << " ";
            par_xy_file << (*par_xy)[k][n] << " ";
            antp_xy_file << (*antp_xy)[k][n] << " ";
        }
    }
    xy_dens_file.close();
    par_xy_file.close();
    antp_xy_file.close();

    close_files();
}

// recall that this function only applies to the LJ force
// and the WCA force - maybe find a correction for the other
// forces implemented
double Properties::truncation_dist() {
    switch (interact_type) {
    case 3:
        truncDist = .5 * boxLength;
        break;
    default:
        truncDist = 2.5;
        truncShift = -1 * (pow(1 / truncDist, 12) - pow(1 / truncDist, 6));
        break;
    }
    return truncDist;
}

// should I put more of the files I am writing to into these functions?
void Properties::open_files() {
    avg_force_particle.open("avgForcePerParticle.txt");
}
void Properties::close_files() { avg_force_particle.close(); }
// assign private variable used in class
void Properties::initializeProperties(Parameters *p) {

    boxLength = p->getBoxLength();
    n_particles = p->getNumParticles();
    sigma = p->getSigma();
    redDens = p->getRedDens();
    red_temp = p->getRedTemp();

    k_spring = p->getSprConst();

    interact_type = p->getInteract_Type();
    rest_L = p->getRestLength();

    a_ref = p->getRefAffinity();
    a_mult = p->getAffinityMult();

    // determines truncation distance
    truncation_dist();
    open_files();
}
