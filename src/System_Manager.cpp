#include <iostream>
#include "System_Manager.h"

// System_Manager::System_Manager() {
//    // ask Adam for a better method for choosing the intervals for the pcf rdf
//    // boxes
//    delta_r = sigma / 20; // this might not be the best way to define delta_r
//    cell_L = sigma / 20;
//
//    // define the various RDF vectors (dependent upon r)
//    int arr_size = 0.5 * boxLength / delta_r + 1;
//    num_density.resize(arr_size);
//    par_num_density.resize(arr_size);
//    antp_num_density.resize(arr_size);
//
//    int val = boxLength / cell_L + 1;
//    //    std::cout << "the size of the vector is " << val << std::endl;
//    xy_num_density.resize(val, std::vector<double>(val, 0));
//    par_xy_density.resize(val, std::vector<double>(val, 0));
//    antp_xy_density.resize(val, std::vector<double>(val, 0));
//}
//
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
//
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
