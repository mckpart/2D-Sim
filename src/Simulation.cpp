#include <iostream>
#include "Simulation.h"

// having this 'reduced energy' made sense in the context of the Lennard Jones
// system, but is confusing for all other systems. so this would be worth
// changing back & then rerunning LJ to make sure nothing broke.
double Simulation::boltzmannFactor(double energy) {
    return exp(-1 * energy / red_temp); // energy here is already reduced
}

Simulation::Simulation(std::string yf) {
    yaml_file = yf;

    param.initializeParameters(yaml_file);  // initialize the parameters
    interact.initializeInteraction(&param); // for the simulation
    bound.initializeBoundary(&param);
    prop.initializeProperties(&param);

    randVal.InitCold(param.getSeed());

    sim_manage = System_Manager(yaml_file);
    part_manage = Particle_Manager(yaml_file);

    Particle prt;
    n_particles = param.getNumParticles(); // initialize vector of
    particles.resize(n_particles);         // particles and set particle
    part_manage.init_particle_params(&particles, &randVal); // parameters

    red_temp = param.getRedTemp();

    trial_position.resize(2);
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

    //    n_particles = 4;
    //    particles.resize(n_particles);
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
    for (int k = 0; k < 1000; k++) {
        writePositions(&pos_file);
    }
    x = -.5;
    y = -.5;
    double r = .5 * sqrt(2);
    //    prop.calc_average_force(x, y, r);
}

// for now leave this routine as is... and check the periodic properties. If
// periodic properties still works, then move this to sim class as weell
// the main loop that actually runs the sim should probably be put in to the
// main.cpp file
void Simulation::calcNonPerProp() {
    Particle curr_prt;
    Particle comp_prt;

    // can probably remove the LJ constant at some point..
    double LJ_constant = 0;
    double r_dist = 0;
    std::vector<double> delta_pos(2);
    // here we should note that these are still being initialized in the
    // properties class and should be here OR there
    double a_ref = param.getRefAffinity();
    double a_mult = param.getAffinityMult();
    double trunc_dist = prop.truncation_dist();

    // make sure the energy is being properly reset each time the properties are
    // updated
    sim_manage.resetEnergy();
    sim_manage.resetVirial();

    for (int k = 0; k < n_particles; k++) {
        // this avg_force should be tracked in the sim manager class
        //        avg_force[0] = 0;
        //        avg_force[1] = 0;
        // this force number should be entirely unnecessary - calculation of
        // forces per particle will be kept more in the particle manageer class
        //        force_num = 0;

        curr_prt = particles[k];

        // set current x,y position
        std::vector<double> curr_pos(2);
        curr_pos[0] = curr_prt.getX_Position();
        curr_pos[1] = curr_prt.getY_Position();

        // each particle-particle interaction
        for (int n = 0; n < n_particles; n++) {

            comp_prt = particles[n];

            // set comparison x,y position
            std::vector<double> comp_pos{comp_prt.getX_Position(),
                                         comp_prt.getY_Position()};

            // the particle cannot interact with itself
            if (curr_prt.getIdentifier() != comp_prt.getIdentifier()) {
                r_dist = prop.radDistance(curr_pos[0], comp_pos[0], curr_pos[1],
                                          comp_pos[1]);
                delta_pos = {comp_pos[0] - curr_pos[0],
                             comp_pos[1] - curr_pos[1]};
                // updates overall number density
                sim_manage.updateNumDensity(r_dist, 0);
                sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1], 0);

                // if type == type interaction of parallel microtubules
                // updates
                // number density for parallel interactions if type !=
                // type :
                // interaction of antiparallel microtubules updates number
                // density for antiparallel interactions
                if (curr_prt.getType() == comp_prt.getType()) {
                    LJ_constant = a_ref;
                    sim_manage.updateNumDensity(r_dist, 1);
                    sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1], 1);
                }
                else if (curr_prt.getType() != comp_prt.getType()) {
                    LJ_constant = a_ref * a_mult;
                    sim_manage.updateNumDensity(r_dist, 2);
                    sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1], 2);
                }
            }
            if (n > k && r_dist < trunc_dist) {
                sim_manage.updateEnergy(prop.calcEnergy(r_dist, LJ_constant));
                // here I have curr - comp whereas everywhere else I have comp -
                // curr - for now keep it the same but come back to check later
                sim_manage.updateVirial(prop.calcVirial(
                    -1 * delta_pos[0], -1 * delta_pos[1], r_dist, LJ_constant));
            }
            //            calc_average_force(x_curr - x_comp, y_curr -
            //            y_comp,
            //            r_dist);
            //            ++force_num;
        }
        // this should probably be write total forces. the tot forces
        //   could then
        // be averaged in a separate python analysis program.
        //        writeAvgForces();
    }
    //    avg_force_particle << "\n";
    //    close_files();
    sim_manage.setTotalVirial();
    sim_manage.setTotalEnergy();

    std::cout << "in sim class, nonperiodic stuff" << std::endl;
}

void Simulation::calcPeriodicProp() {
    Particle curr_prt;
    Particle comp_prt;

    // can probably remove the LJ constant at some point..
    double LJ_constant = 0;
    double r_dist = 0;
    std::vector<double> delta_pos(2);

    double a_ref = param.getRefAffinity();
    double a_mult = param.getAffinityMult();
    double trunc_dist = prop.truncation_dist();

    sim_manage.resetEnergy();
    sim_manage.resetVirial();

    std::vector<std::vector<double>> cellPositions(9,
                                                   std::vector<double>(2, 0));
    // what if in particle manager class, the object had a curr particle, ref
    // particle, etc, to help track the behvaior here
    for (int k = 0; k < n_particles; k++) {

        curr_prt = particles[k];

        // set current x,y position
        std::vector<double> curr_pos(2);
        curr_pos[0] = curr_prt.getX_Position();
        curr_pos[1] = curr_prt.getY_Position();
        // takes into account each particle-particle interaction
        for (int n = 0; n < n_particles; n++) {

            comp_prt = particles[n];

            // set comparison x,y position
            std::vector<double> comp_pos{comp_prt.getX_Position(),
                                         comp_prt.getY_Position()};

            // the particle cannot interact with itself
            if (curr_prt.getIdentifier() != comp_prt.getIdentifier()) {
                r_dist = prop.radDistance(curr_pos[0], comp_pos[0], curr_pos[1],
                                          comp_pos[1]);
                delta_pos = {comp_pos[0] - curr_pos[0],
                             comp_pos[1] - curr_pos[1]};
                // updates overall number density
                sim_manage.updateNumDensity(r_dist, 0);
                sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1], 0);

                // if type == type: interaction of parallel microtubules
                // and parallel num density is updated
                // if type != type: interaction of antiparallel
                //                microtubules
                // and antiparallel num density is updated
                if (curr_prt.getType() == comp_prt.getType()) {
                    LJ_constant = a_ref;
                    sim_manage.updateNumDensity(r_dist, 1);
                    sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1], 1);
                }
                else if (curr_prt.getType() != comp_prt.getType()) {
                    LJ_constant = a_ref * a_mult;
                    sim_manage.updateNumDensity(r_dist, 2);
                    sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1], 2);
                }
            }
            // this should certainly be broken up into smaller functions
            if (n > k) {
                if (r_dist > trunc_dist) {
                    prop.populateCellArray(comp_pos[0], comp_pos[1],
                                           &cellPositions);
                    for (int z = 0; z < 8; z++) {
                        // creates the 8 cell images
                        comp_pos = {cellPositions[z][0], cellPositions[z][1]};

                        r_dist = prop.radDistance(curr_pos[0], comp_pos[0],
                                                  curr_pos[1], comp_pos[1]);

                        for (int j = 0; j < 2; j++) {
                            sim_manage.updateNumDensity(r_dist, 0);
                            sim_manage.calc_xy_dens(delta_pos[0], delta_pos[1],
                                                    0);

                            if (curr_prt.getType() == comp_prt.getType()) {
                                sim_manage.updateNumDensity(r_dist, 1);
                                sim_manage.calc_xy_dens(delta_pos[0],
                                                        delta_pos[1], 1);
                            } else if (curr_prt.getType() !=
                                       comp_prt.getType()) {
                                sim_manage.updateNumDensity(r_dist, 2);
                                sim_manage.calc_xy_dens(delta_pos[0],
                                                        delta_pos[1], 2);
                            }
                        }
                        if (r_dist < trunc_dist) {
                            for (int j = 0; j < 2; j++) {
                                sim_manage.updateEnergy(
                                    prop.calcEnergy(r_dist, LJ_constant));
                                // see comments about virial calculation above
                                sim_manage.updateVirial(prop.calcVirial(
                                    -1 * delta_pos[0], -1 * delta_pos[1],
                                    r_dist, LJ_constant));
                            }
                        }
                    }
                } else {
                    sim_manage.updateEnergy(
                        prop.calcEnergy(r_dist, LJ_constant));
                    // see comments about virial calculation above
                    sim_manage.updateVirial(
                        prop.calcVirial(-1 * delta_pos[0], -1 * delta_pos[1],
                                        r_dist, LJ_constant));
                }
            }
        }
    }
    sim_manage.setTotalVirial();
    sim_manage.setTotalEnergy();

    std::cout << "in sim class, periodic stuff" << std::endl;
}
// end sim if this error is thrown
void Simulation::init_configuration() {
    double n_initial = n_particles;

    // IT MIGHT BE EVEN BETTER TO PUT THIS INTO A PARTICLE MANAGER CLASS
    // make this into a function - could probaby write some sort of catch2 test
    // that checks the various types of initializations of the system
    // initializing particle positions
    if (param.getInit_Type() == 0) {
        bound.initialPosition(&particles, randVal);
    } else if (param.getInit_Type() == 1) {
        n_initial = bound.initialHexagonal(&particles);
    } else if (param.getInit_Type() == 2) {
        n_initial = bound.initialSquare(&particles);
    }
    // this should be placed into the same function mentioned above.
    // print warning if the system has too many particles
    if (n_initial < n_particles) {
        std::cout << "ERROR: TOO MANY PARTICLES. INITIALIZED " << n_initial
                  << " PARTICLES" << std::endl;
    } else {
        std::cout << "SUCCESSFULLY INITIALIZED " << n_initial << " PARTICLES"
                  << std::endl;
    }
}

void Simulation::periodic_pos_trial(Particle *p) {

    // updates trial position in function then particle - particle
    // interactions
    // for now I am going to make x_trial and y_trial private variables of the
    // sim class but they should be made into properties of the particle manager
    // class later
    trial_position[0] = p->getX_TrialPos();
    trial_position[1] = p->getY_TrialPos();
    double curr_index = p->getIdentifier();
    // later, rather than doing curr_index, recall that each particle has an
    // identifier that relates it back to its index in the particle array
    if (param.getInteract_Type() != 0) {
        delta_energy = interact.periodicInteraction(&particles, curr_index);
    }
}

bool Simulation::nonperiodic_pos_trial(Particle *p, bool accept) {
    if (param.getBound_Type() == 0) {
        accept = bound.rigidBoundary(&particles, p->getIdentifier());
    } else if (param.getBound_Type() == 2) {
        delta_energy = bound.externalWell(&particles, p->getIdentifier());
    }

    if (param.getInteract_Type() != 0) {
        delta_energy = delta_energy + interact.nonPeriodicInteraction(
                                          &particles, p->getIdentifier());
    }
    return accept;
}
// once sim is working again, consider making these private functions since they
// don't need to be used outside of the sim class

void Simulation::runSimulation() {

    Particle prt;

    double n_rejects = 0;
    double perc_rej = 0;

    std::ofstream rad_dist_file;
    rad_dist_file.open("radialDistance.txt");

    std::ofstream pos_file;
    pos_file.open("positions.txt");

    double n_updates = param.getUpdates();

    init_configuration();
    // turn this into a function like sweep()
    for (int sweepNum = 0; sweepNum < n_updates; sweepNum++) {
        for (int k = 0; k < n_particles; k++) {

            // choose random particle
            int curr_index = int(randVal.RandomUniformDbl() * n_particles);
            prt = particles[curr_index];

            // generate and set the x,y trial position
            trial_position[0] = prt.x_trial(randVal.RandomUniformDbl());
            trial_position[1] = prt.y_trial(randVal.RandomUniformDbl());

            // might be worth just putting this into the create trial position
            // function in particle class
            prt.setX_TrialPos(trial_position[0]);
            prt.setY_TrialPos(trial_position[1]);
            particles[curr_index] = prt;

            bool accept = 1;
            delta_energy = 0;
            //            double delta_energy = 0; // sets change in energy to 0

            if (param.getBound_Type() == 1) {

                // run sim with periodic boundaries
                bound.periodicBoundary(&particles, curr_index);
                prt = particles[curr_index];
                periodic_pos_trial(&prt);
            } else {
                accept = nonperiodic_pos_trial(&prt, accept);
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
            // recall that for hard disks hte delta energy will always remain 0
            // & there's no need for redundancy
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
                prt.setX_Position(trial_position[0]);
                prt.setY_Position(trial_position[1]);
                particles[curr_index] = prt;
            } else {
                n_rejects++; // keeps count of total moves rejected
            }
        }

        if (sweepNum > param.getEq_sweep() &&
            sweepNum % param.getData_interval() == 0) {
            std::cout << "current sweep: " << sweepNum << std::endl;
            writePositions(&pos_file);
            // recall that these functions are going to be moved into this class
            if (param.getBound_Type() == 1) {
                calcPeriodicProp();
            } else {
                calcNonPerProp();
            }
        }
    }

    prop.writeProperties(
        sim_manage.getTotalEnergy(), sim_manage.getTotalVirial(),
        sim_manage.getTotalNumDensity(), sim_manage.getParNumDensity(),
        sim_manage.getAntiNumDensity(), sim_manage.getTotalXYDensity(),
        sim_manage.getParXYDensity(), sim_manage.getAntiXYDensity());
    //   std::cout << "The average energy of the system is " <<
    //   prop.calcAvgEnergy() << std::endl; std::cout << "The pressure of the
    //   system is " << prop.c alcPressure() << std::endl;
    perc_rej = n_rejects / (n_updates * n_particles) * 100.0;
    std::cout << perc_rej << "% of the moves were rejected." << std::endl;
}

