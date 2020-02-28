#include <iostream>

#include "Particle_Manager.h"

Particle_Manager::Particle_Manager() {}

// MOVE THE PARTICLE VECTOR INTO THIS CLASS, AND RESIZE IT IN THIS CONSTRUCTOR
Particle_Manager::Particle_Manager(std::string yaml_file) {
    // set params & properties object
    param.initializeParameters(yaml_file); // initialize the parameters
    prop.initializeProperties(&param);
}

// Particle_Manager:: function that assigns forces to appropriate particles
/* 1. Since each particle acts on every other particle equally and oppositely,
 * loop through all possible combinations
 *     - this can be achieved by doing outer loop < npart and inner loop < n <
 *       n_particles in order to avoid double counting interactions
 * 2. This should be done each time right before the system's state is captured
 *     - when the state is captured, we then write the total force vecs to some
 *       file
 *     - forces can then be averaged over the entire run in external analysis
 *       code
 * 3. Before computing the forces again, the force for each particle should be
 * reset
 */

// MAKE SURE THAT THIS IS CHANGING THE CORRECT PARTICLE IN THE VECTOR OF
// PARTICLES!!!!!! //
void Particle_Manager::assign_forces(double r, double LJ_const, Particle *curr,
                                     Particle *ref) {
    std::vector<double> f(2);
    std::vector<double> del = {ref->getX_Position() - curr->getX_Position(),
                               ref->getY_Position() - curr->getY_Position()};
    prop.calc_force_vec(del[0], del[1], r, &f);
    ref->addForce(f[0], f[1]);
    curr->addForce(-1 * f[0], -1 * f[1]);
}

void Particle_Manager::reset_forces(std::vector<Particle> *particles) {
    for (int k = 0; k < param.getNumParticles(); ++k) {
        (*particles)[k].resetForce();
    }
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

