#include <iostream>

#include "Particle_Manager.h"

Particle_Manager::Particle_Manager() {}

Particle_Manager::Particle_Manager(std::string yaml_file) {
    // set params & properties object
    param.initializeParameters(yaml_file); // initialize the parameters
    prop.initializeProperties(&param);
}

void Particle_Manager::init_particle_params(std::vector<Particle> *particles,
                                            KISSRNG *randVal) {
    Particle prt;

    std::ofstream type_file;
    type_file.open("particle_type.txt");

    double sigma = param.getSigma();
    double boxLength = param.getBoxLength();

    // ratio of parallel to total
    double ratio = (double)param.getNumParallel() / param.getNumParticles();
    std::cout << "the ratio is: " << ratio << "\n";

    double weight = sigma * sqrt(1 / (4 * param.getRedDens()));
    std::cout << "the stepping weight is: " << weight << std::endl;
    prt.setStepWeight(weight);

    double radius = .5 * sigma;
    prt.setRadius(radius);

    double num_par = 0;
    double num_ant = 0;
    int type = 0;

    for (int k = 0; k < param.getNumParticles(); ++k) {
        if (num_par < param.getNumParallel() &&
            num_ant < param.getNumAntiparallel()) {

            double n = randVal->RandomUniformDbl();
            if (n < ratio) {
                type = 1;
                ++num_par;
            } else {
                type = 2;
                ++num_ant;
            }
        } else if (num_par < param.getNumParallel()) {
            type = 1;
            ++num_par;
            if (k == 1) {
                std::cout << "type is 1" << std::endl;
            }
        } else if (num_ant < param.getNumAntiparallel()) {
            type = 2;
            ++num_ant;
        } else {
            std::cout << "count again" << std::endl;
        }
        prt.setType(type);
        prt.setIdentifier(k);
        (*particles)[k] = prt;

        type_file << type << " ";
    }
    type_file.close();
}

