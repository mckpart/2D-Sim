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
void Properties::calcEnergy(double r, double a) {
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
    f_energy = f_energy + val;
}

// sums the total virial of the current configuration
void Properties::calcVirial(double x, double y, double r, double a) {
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

    f_r = f_r + r * val;
}

// void Properties::updateNumDensity(double r, int ID) {
//    int val = r / delta_r;
//    int index = 0;
//
//    if (r < 0.5 * boxLength) {
//        if (r > (val + 0.5) * delta_r) {
//            index = val + 1;
//        } else {
//            index = val;
//        }
//
//        switch (ID) {
//        case 0:
//            num_density[index] = num_density[index] + 1;
//            break;
//        case 1:
//            par_num_density[index] = par_num_density[index] + 1;
//            break;
//        case 2:
//            antp_num_density[index] = antp_num_density[index] + 1;
//            break;
//        }
//    }
//}

// void Properties::calc_xy_dens(double x, double y, int ID) {
//    double half_boxL = .5 * boxLength;
//    int ind_1 = (x + half_boxL) / cell_L; // cell_L = delta x = delta y
//    int ind_2 = (y + half_boxL) / cell_L;
//
//    if (fabs(x) < half_boxL - cell_L && fabs(y) < half_boxL - cell_L) {
//        if (x > (ind_1 + 0.5) * cell_L - half_boxL) {
//            ++ind_1;
//        }
//        if (y > (ind_2 + 0.5) * cell_L - half_boxL) {
//            ++ind_2;
//        }
//
//        // increment vales of appropriate number densities
//        switch (ID) {
//        case 0:
//            xy_num_density[ind_1][ind_2] = xy_num_density[ind_1][ind_2] + 1;
//            break;
//        case 1:
//            par_xy_density[ind_1][ind_2] = par_xy_density[ind_1][ind_2] + 1;
//            break;
//        case 2:
//            antp_xy_density[ind_1][ind_2] = antp_xy_density[ind_1][ind_2] + 1;
//            break;
//        }
//    }
//}

// for now leave this routine as is... and check the periodic properties. If
// periodic properties still works, then move this to sim class as weell
void Properties::calcNonPerProp(std::vector<Particle> *particles) {
    //    Particle curr_prt;
    //    Particle comp_prt;
    //
    //    // can probably remove the L
    //    double LJ_constant = 0;
    //    double r_dist = 0;
    //
    //    // make sure that the free energy previously calculated is reset the
    //    free
    //    // energy is only the energy that comes from the positions within the
    //    // configuration
    //    f_energy = 0;
    //    f_r = 0;
    //
    //    for (int k = 0; k < n_particles; k++) {
    //        avg_force[0] = 0;
    //        avg_force[1] = 0;
    //        force_num = 0;
    //
    //        curr_prt = (*particles)[k];
    //
    //        // set current x,y position
    //        double x_curr = curr_prt.getX_Position();
    //        double y_curr = curr_prt.getY_Position();
    //
    //        // each particle-particle interaction
    //        for (int n = 0; n < n_particles; n++) {
    //
    //            comp_prt = (*particles)[n];
    //
    //            // set comparison x,y position
    //            double x_comp = comp_prt.getX_Position();
    //            double y_comp = comp_prt.getY_Position();
    //
    //            // the particle cannot interact with itself
    //            if (curr_prt.getIdentifier() != comp_prt.getIdentifier()) {
    //                r_dist = radDistance(x_curr, x_comp, y_curr, y_comp);
    //
    //                // updates overall number density
    //                updateNumDensity(r_dist, 0);
    //                calc_xy_dens(x_comp - x_curr, y_comp - y_curr, 0);
    //
    //                // if type == type interaction of parallel microtubules
    //                updates
    //                // number density for parallel interactions if type !=
    //                type:
    //                // interaction of antiparallel microtubules updates number
    //                // density for antiparallel interactions
    //                if (curr_prt.getType() == comp_prt.getType()) {
    //                    LJ_constant = a_ref;
    //                    updateNumDensity(r_dist, 1);
    //                    calc_xy_dens(x_comp - x_curr, y_comp - y_curr, 1);
    //                } else if (curr_prt.getType() != comp_prt.getType()) {
    //                    LJ_constant = a_ref * a_mult;
    //                    updateNumDensity(r_dist, 2);
    //                    calc_xy_dens(x_comp - x_curr, y_comp - y_curr, 2);
    //                }
    //            }
    //            if (n > k && r_dist < truncDist) {
    //                calcEnergy(r_dist, LJ_constant);
    //                calcVirial(x_curr - x_comp, y_curr - y_comp, r_dist,
    //                           LJ_constant);
    //            }
    //            //            calc_average_force(x_curr - x_comp, y_curr -
    //            y_comp,
    //            //            r_dist);
    //            ++force_num;
    //        }
    //        // this should probably be write total forces. the tot forces
    //        could then
    //        // be averaged in a separate python analysis program.
    //        writeAvgForces();
    //    }
    //    avg_force_particle << "\n";
    //    //    close_files();
    //    sum_Fdot_r.push_back(f_r);
    //    sum_energy.push_back(f_energy);

    std::cout << "in properties, nonperiodic stuff" << std::endl;
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

void Properties::calcPeriodicProp(std::vector<Particle> *particles) {
    //    Particle curr_prt;
    //    Particle comp_prt;
    //
    //    double LJ_constant = 0;
    //    double r_dist = 0;
    //
    //    std::vector<std::vector<double>> cellPositions(9,
    //                                                   std::vector<double>(2,
    //                                                   0));
    //    // make sure that the free energy previously calculated is reset
    //    // the free energy is only the energy that comes from the positions
    //    // within the configuration
    //    f_energy = 0;
    //    f_r = 0;
    //
    //    for (int k = 0; k < n_particles; k++) {
    //        curr_prt = (*particles)[k];
    //
    //        // set current x,y position
    //        double x_curr = curr_prt.getX_Position();
    //        double y_curr = curr_prt.getY_Position();
    //
    //        // takes into account each particle-particle interaction
    //        for (int n = 0; n < n_particles; n++) {
    //
    //            comp_prt = (*particles)[n];
    //
    //            // set comparison x,y position
    //            double x_comp = comp_prt.getX_Position();
    //            double y_comp = comp_prt.getY_Position();
    //
    //            // the particle cannot interact with itself
    //            if (curr_prt.getIdentifier() != comp_prt.getIdentifier()) {
    //                r_dist = radDistance(x_curr, x_comp, y_curr, y_comp);
    //
    //                // updates total radial num density and x,y num density
    //                updateNumDensity(r_dist, 0);
    //                calc_xy_dens(x_comp - x_curr, y_comp - y_curr, 0);
    //
    //                // if type == type: interaction of parallel microtubules
    //                // and parallel num density is updated
    //                // if type != type: interaction of antiparallel
    //                microtubules
    //                // and antiparallel num density is updated
    //                if (curr_prt.getType() == comp_prt.getType()) {
    //                    LJ_constant = a_ref;
    //                    updateNumDensity(r_dist, 1);
    //                    calc_xy_dens(x_comp - x_curr, y_comp - y_curr, 1);
    //                } else if (curr_prt.getType() != comp_prt.getType()) {
    //                    LJ_constant = a_ref * a_mult;
    //                    updateNumDensity(r_dist, 2);
    //                    calc_xy_dens(x_comp - x_curr, y_comp - y_curr, 2);
    //                }
    //            }
    //
    //            if (n > k) {
    //                if (r_dist > truncDist) {
    //                    populateCellArray(x_comp, y_comp, &cellPositions);
    //                    for (int z = 0; z < 8; z++) {
    //                        // creates the 8 cell images
    //                        x_comp = cellPositions[z][0];
    //                        y_comp = cellPositions[z][1];
    //                        r_dist = radDistance(x_curr, x_comp, y_curr,
    //                        y_comp);
    //
    //                        for (int j = 0; j < 2; j++) {
    //                            updateNumDensity(r_dist, 0);
    //                            calc_xy_dens(x_comp - x_curr, y_comp - y_curr,
    //                            0);
    //
    //                            if (curr_prt.getType() == comp_prt.getType())
    //                            {
    //                                updateNumDensity(r_dist, 1);
    //                                calc_xy_dens(x_comp - x_curr, y_comp -
    //                                y_curr,
    //                                             1);
    //                            } else if (curr_prt.getType() !=
    //                                       comp_prt.getType()) {
    //                                updateNumDensity(r_dist, 2);
    //                                calc_xy_dens(x_comp - x_curr, y_comp -
    //                                y_curr,
    //                                             2);
    //                            }
    //                        }
    //                        if (r_dist < truncDist) {
    //                            for (int j = 0; j < 2; j++) {
    //                                calcEnergy(r_dist, LJ_constant);
    //                                calcVirial(x_curr - x_comp, y_curr -
    //                                y_comp,
    //                                           r_dist, LJ_constant);
    //                            }
    //                        }
    //                    }
    //                } else {
    //                    calcEnergy(r_dist, LJ_constant);
    //                    calcVirial(x_curr - x_comp, y_curr - y_comp, r_dist,
    //                               LJ_constant);
    //                }
    //            }
    //        }
    //    }
    //    sum_Fdot_r.push_back(f_r);
    //    sum_energy.push_back(f_energy);
    std::cout << "in properties, periodic stuff" << std::endl;
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
void Properties::writeProperties() {
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

    len = double(sum_Fdot_r.size()); // the force and energy vector are the same
    for (int k = 0; k < len; k++) {  // size hence are put into one for-loop
        virial_file << sum_Fdot_r[k] << " ";
        energy_file << sum_energy[k] << " ";
    }
    virial_file.close();
    energy_file.close();

    len = double(num_density.size());
    for (int k = 0; k < len; k++) {           // all of the number density
        n_dens_file << num_density[k] << " "; // vectors are the same length
        antp_dens_file << antp_num_density[k] << " ";
        par_dens_file << par_num_density[k] << " ";
    }
    n_dens_file.close(); // close all files written to
    par_dens_file.close();
    antp_dens_file.close();

    len = double(xy_num_density.size());
    for (int k = 0; k < len; ++k) {
        for (int n = 0; n < len; ++n) {
            xy_dens_file << xy_num_density[k][n] << " ";
            par_xy_file << par_xy_density[k][n] << " ";
            antp_xy_file << antp_xy_density[k][n] << " ";
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
void Properties::truncation_dist() {
    switch (interact_type) {
    case 3:
        truncDist = .5 * boxLength;
        break;
    default:
        truncDist = 2.5;
        truncShift = -1 * (pow(1 / truncDist, 12) - pow(1 / truncDist, 6));
        break;
    }
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
    //    std::cout << k_spring << "  " << a_ref << std::endl;

    // may be worth putting this chunk of code into a separate function
    delta_r = sigma / 20; // this might not be the best way to define delta_r
    cell_L = sigma / 20;

    // determines truncation distance
    truncation_dist();
    open_files();

    // define the various RDF vectors(dependent upon r)
    int arr_size = 0.5 * boxLength / delta_r + 1;
    num_density.resize(arr_size);
    par_num_density.resize(arr_size);
    antp_num_density.resize(arr_size);

    int val = boxLength / cell_L + 1;
    //    std::cout << "the size of the vector is " << val << std::endl;
    xy_num_density.resize(val, std::vector<double>(val, 0));
    par_xy_density.resize(val, std::vector<double>(val, 0));
    antp_xy_density.resize(val, std::vector<double>(val, 0));
}
