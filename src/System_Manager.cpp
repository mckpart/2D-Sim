#include <iostream>
#include "System_Manager.h"

System_Manager::System_Manager() {}

System_Manager::System_Manager(std::string yaml_file) {
    // ask Adam for a better method for choosing the intervals for the pcf rdf
    // boxes
    // also - it might be better to actually pass these objects into this
    // function

    param.initializeParameters(yaml_file); // initialize the parameters
    prop.initializeProperties(&param);
    delta_r = param.getSigma() /
              20; // this might not be the best way to define delta_r
    cell_L = param.getSigma() / 20;

    // define the various RDF vectors (dependent upon r)
    int arr_size = 0.5 * param.getBoxLength() / delta_r + 1;
    num_density.resize(arr_size);
    par_num_density.resize(arr_size);
    antp_num_density.resize(arr_size);

    int val = param.getBoxLength() / cell_L + 1;
    //    std::cout << "the size of the vector is " << val << std::endl;
    xy_num_density.resize(val, std::vector<double>(val, 0));
    par_xy_density.resize(val, std::vector<double>(val, 0));
    antp_xy_density.resize(val, std::vector<double>(val, 0));
}

void System_Manager::updateNumDensity(double r, int ID) {
    int val = r / delta_r;
    int index = 0;

    if (r < 0.5 * param.getBoxLength()) {
        if (r > (val + 0.5) * delta_r) {
            index = val + 1;
        } else {
            index = val;
        }

        switch (ID) {
        case 0:
            num_density[index] = num_density[index] + 1;
            break;
        case 1:
            par_num_density[index] = par_num_density[index] + 1;
            break;
        case 2:
            antp_num_density[index] = antp_num_density[index] + 1;
            break;
        }
    }
}

void System_Manager::calc_xy_dens(double x, double y, int ID) {
    double half_boxL = .5 * param.getBoxLength();
    int ind_1 = (x + half_boxL) / cell_L; // cell_L = delta x = delta y
    int ind_2 = (y + half_boxL) / cell_L;

    if (fabs(x) < half_boxL - cell_L && fabs(y) < half_boxL - cell_L) {
        if (x > (ind_1 + 0.5) * cell_L - half_boxL) {
            ++ind_1;
        }
        if (y > (ind_2 + 0.5) * cell_L - half_boxL) {
            ++ind_2;
        }

        // increment vales of appropriate number densities
        switch (ID) {
        case 0:
            xy_num_density[ind_1][ind_2] = xy_num_density[ind_1][ind_2] + 1;
            break;
        case 1:
            par_xy_density[ind_1][ind_2] = par_xy_density[ind_1][ind_2] + 1;
            break;
        case 2:
            antp_xy_density[ind_1][ind_2] = antp_xy_density[ind_1][ind_2] + 1;
            break;
        }
    }
}
// ahead of time, determine how many 'system updates' there will be to set the
// size of total virial and free energy vectors
void System_Manager::updateVirial(double v) { f_r += v; }
void System_Manager::setTotalVirial() { sum_virial.push_back(f_r); }
void System_Manager::resetVirial() { f_r = 0; }

void System_Manager::updateEnergy(double e) { f_energy += e; }
void System_Manager::setTotalEnergy() { sum_energy.push_back(f_energy); }
void System_Manager::resetEnergy() { f_energy = 0; }

