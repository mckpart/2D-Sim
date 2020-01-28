#include <iostream>
#include "Simulation.h"

double Simulation::boltzmannFactor(double energy) {
    return exp(-1 * energy / red_temp); // energy here is already reduced
}

Simulation::Simulation(std::string yf) {
    yamlFile = yf;

    param.initializeParameters(yamlFile);   // initialize the parameters
    interact.initializeInteraction(&param); // for the simulation
    bound.initializeBoundary(&param);
    prop.initializeProperties(&param);

    randVal.InitCold(param.getSeed());

    n_particles = param.getNumParticles(); // initialize vector of
    particles.resize(n_particles);         // particles and set particle
    setParticleParams();                   // parameters

    Particle prt;

    red_temp = param.getRedTemp();
}

void Simulation::writePositions(std::ofstream *pos_file) {
    Particle prt;

    if (pos_file->is_open()) {

        for (int k = 0; k < n_particles; k++) {
            prt = particles[k];

            // writes updated positions into position file
            (*pos_file) << prt.getX_Position() << " ";
            (*pos_file) << prt.getY_Position() << " ";
        }
        (*pos_file) << std::endl;
    } else {
        std::cout << "ERROR: THE .TXT FILE COULD NOT OPEN" << std::endl;
    }
}

// eventually replace this test with some sort of Catch2

void Simulation::testSimulation() {

    n_particles = 4;
    particles.resize(n_particles);
    Particle prt;
    double x = 0;
    double y = 0;

    std::ofstream pos_file;
    pos_file.open("positions.txt");

    //    for (int k = 0; k < 4; k++) {
    //        if (k == 0) {
    //            x = .25;
    //            y = .25;
    //        } else if (k == 1) {
    //            x = .25;
    //            y = -.25;
    //        } else if (k == 2) {
    //            x = -.25;
    //            y = .25;
    //        } else if (k == 3) {
    //            x = -.25;
    //            y = -.25;
    //        }
    //
    //        prt.setX_Position(x);
    //        prt.setY_Position(y);
    //        particles[k] = prt;
    //    }
    //
    //    for (int k = 0; k < 1000; k++) {
    //        writePositions(&pos_file);
    //    }
    x = -.5;
    y = -.5;
    double r = .5 * sqrt(2);
    prop.calc_average_force(x, y, r);
}

void Simulation::runSimulation() {

    Particle prt;

    double n_rejects = 0;
    double perc_rej = 0;

    std::ofstream rad_dist_file;
    rad_dist_file.open("radialDistance.txt");

    std::ofstream pos_file;
    pos_file.open("positions.txt");

    double n_updates = param.getUpdates();
    double n_initial = n_particles;

    // initializing particle positions
    if (param.getInit_Type() == 0) {
        bound.initialPosition(&particles, randVal);
    } else if (param.getInit_Type() == 1) {
        n_initial = bound.initialHexagonal(&particles);
    } else if (param.getInit_Type() == 2) {
        n_initial = bound.initialSquare(&particles);
    }

    // print warning if the system has too many particles
    if (n_initial < n_particles) {
        std::cout << "ERROR: TOO MANY PARTICLES. INITIALIZED " << n_initial
                  << " PARTICLES" << std::endl;
    } else {
        std::cout << "SUCCESSFULLY INITIALIZED " << n_initial << " PARTICLES"
                  << std::endl;
    }
    for (int sweepNum = 0; sweepNum < n_updates; sweepNum++) {
        for (int k = 0; k < n_particles; k++) {

            // choose random particle
            int curr_index = int(randVal.RandomUniformDbl() * n_particles);
            prt = particles[curr_index];

            // generate and set the x,y trial position
            double x_trial = prt.x_trial(randVal.RandomUniformDbl());
            double y_trial = prt.y_trial(randVal.RandomUniformDbl());

            prt.setX_TrialPos(x_trial);
            prt.setY_TrialPos(y_trial);
            particles[curr_index] = prt;

            bool accept = 1;
            double delta_energy = 0; // sets change in energy to 0

            if (param.getBound_Type() == 1) {

                // run sim with periodic boundaries
                bound.periodicBoundary(&particles, curr_index);
                prt = particles[curr_index];

                // updates trial position in function then particle - particle
                // interactions
                x_trial = prt.getX_TrialPos();
                y_trial = prt.getY_TrialPos();

                if (param.getInteract_Type() != 0) {
                    delta_energy =
                        interact.periodicInteraction(&particles, curr_index);
                }
            } else {
                if (param.getBound_Type() == 0) {
                    accept = bound.rigidBoundary(&particles, curr_index);
                } else if (param.getBound_Type() == 2) {
                    delta_energy = bound.externalWell(&particles, curr_index);
                }

                if (param.getInteract_Type() != 0) {
                    delta_energy =
                        delta_energy +
                        interact.nonPeriodicInteraction(&particles, curr_index);
                }
            }

            /* RUNS DIFFERENT TYPES OF PARTICLE-PARTICLE INTERACTIONS
             * HARD DISKS IS A 0 - 1 PROBABILITY THUS A DELTA ENERGY IS NOT
             * RETURNED THE CHANGE IN ENERGY IS RETURNED FROM LENJONES AND WCA
             * POTENTIAL THIS TOTAL CHANGE IS SENT INTO THE BOLTZMANN FACTOR
             * FUNCTION TO CALCULATE THE TOTAL PROBABILITY OF ACCEPTING THE
             * TRIAL MOVE. IF THE TOTAL CHANGE < 0, THE MOVE IS ACCEPTED. ELSE A
             * RANDOM NUMBER IS GENERATED TO DETERMINE WHETHER THE MOVE IS TO BE
             * ACCEPTED
             */

            if (param.getInteract_Type() == 0 && accept == 1) {
                accept = interact.hardDisks(&particles, curr_index);
            }

            if (accept == 1 && delta_energy > 0) {
                // compute acceptance probability
                double total_prob = boltzmannFactor(delta_energy);

                if (randVal.RandomUniformDbl() < total_prob) {
                    accept = 1;
                } else {
                    accept = 0;
                }
            }

            // if trial move is accepted, update the position of current
            // particle and put back in particle vector
            if (accept == 1) {
                prt.setX_Position(x_trial);
                prt.setY_Position(y_trial);
                particles[curr_index] = prt;
            } else {
                n_rejects++; // keeps count of total moves rejected
            }
        }

        if (sweepNum > param.getEq_sweep() &&
            sweepNum % param.getData_interval() == 0) {
            std::cout << "current sweep: " << sweepNum << std::endl;
            writePositions(&pos_file);
            if (param.getBound_Type() == 1) {
                prop.calcPeriodicProp(&particles);
            } else {
                prop.calcNonPerProp(&particles);
            }
        }
    }
    prop.writeProperties();
    //   std::cout << "The average energy of the system is " <<
    //   prop.calcAvgEnergy() << std::endl; std::cout << "The pressure of the
    //   system is " << prop.c alcPressure() << std::endl;
    perc_rej = n_rejects / (n_updates * n_particles) * 100.0;
    std::cout << perc_rej << "% of the moves were rejected." << std::endl;
}

// THIS IS THE NEXT PIECE TO BE ALTERED ////

void Simulation::setParticleParams() {

    Particle prt;

    std::ofstream type_file;
    type_file.open("particle_type.txt");

    YAML::Node node = YAML::LoadFile(yamlFile);

    int num_part_1 = node["type1_Particles"].as<int>();
    int num_part_2 = node["type2_Particles"].as<int>();

    // both particles have the same radius
    double radius = node["particleRadius"].as<double>();

    double sigma = param.getSigma();
    double boxLength = param.getBoxLength();

    double ratio = (double)num_part_1 / n_particles;
    std::cout << "the ratio is: " << ratio << "\n";

    double weight = sigma * sqrt(1 / (4 * param.getRedDens()));
    std::cout << "the stepping weight is: " << weight << std::endl;
    prt.setStepWeight(weight);

    if (param.getInteract_Type() != 0) { // applies to the LJ and WCA potentials
        radius = .5 * sigma;
    }
    prt.setRadius(radius);

    double num_1 = 0;
    double num_2 = 0;
    int type = 0;

    for (int k = 0; k < n_particles; ++k) {
        if (num_1 < num_part_1 && num_2 < num_part_2) {

            double n = randVal.RandomUniformDbl();
            if (n < ratio) {
                type = 1;
                ++num_1;
            } else {
                type = 2;
                ++num_2;
            }
        } else if (num_1 < num_part_1) {
            type = 1;
            ++num_1;
        } else if (num_2 < num_part_2) {
            type = 2;
            ++num_2;
        } else {
            std::cout << "count again" << std::endl;
        }
        prt.setType(type);
        prt.setIdentifier(k);
        particles[k] = prt;

        type_file << type << " ";
    }
    type_file.close();
}
