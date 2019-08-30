#include <iostream>
#include <cmath>

#include "Boundary.h"

double dist_wall(double x1, double x2) { return fabs(x2 - x1); }

double dist_part(double x1, double x2, double y1, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

bool Boundary::rigidBoundary(std::vector<Particle> *particles, int index) {

    Particle current_prt;

    double x_wall = 0;
    double y_wall = 0;
    double num = 0;

    current_prt = (*particles)[index];

    // reads in the current x,y trial position and the radius of the trial
    // particle
    double x_temp = current_prt.getX_TrialPos();
    double y_temp = current_prt.getY_TrialPos();
    double rad_temp = current_prt.getRadius();

    bool accept = 1;

    if (x_temp * -1 < 0) {
        x_wall = 0.5 * boxLength; // x_temp > 0, closest x_wall > 0

        if (x_temp - x_wall > 0) { // the particle moved past the wall
            accept = 0;            // along the x-boundary
        }
    } else {
        x_wall = -0.5 * boxLength; // x_temp < 0, closest wall < 0

        if (x_temp - x_wall < 0) {
            accept = 0;
        }
    }

    if (y_temp * -1 < 0) {
        y_wall = 0.5 * boxLength; // y_temp > 0, closest y_wall > 0

        if (y_temp - y_wall > 0) { // the particle moved past the wall
            accept = 0;            // along the y-boundary
        }
    } else {
        y_wall = -0.5 * boxLength;

        if (y_temp - y_wall < 0) {
            accept = 0;
        }
    }
    // the particle's center of mass is less than a radius' distance from the
    // x-boundary
    if (fabs(x_temp - x_wall) < rad_temp) {
        accept = 0;
        // is less than a radius' distance from the y-boundary
    } else if (fabs(y_temp - y_wall) < rad_temp) {
        accept = 0;
    }

    // returns 1 if trial move is accepted
    return accept;
}

void Boundary::periodicBoundary(std::vector<Particle> *particles, int index) {

    Particle current_prt;

    double x_wall = 0;
    double y_wall = 0;

    double dist_xtravel = 0;
    double dist_ytravel = 0;
    double dist_xwall = 0;
    double dist_ywall = 0;

    double wallBound = 0.5 * boxLength;

    current_prt = (*particles)[index];

    // sets the current particle's x,y sets the current particle's x,y position
    double x_curr = current_prt.getX_Position();
    double y_curr = current_prt.getY_Position();

    // set the x,y trial position for current particle
    double x_temp = current_prt.getX_TrialPos();
    double y_temp = current_prt.getY_TrialPos();
    double rad_temp = current_prt.getRadius();

    /* FINDS THE NEAREST X,Y WALLS
     * IF THE PARTICLE HAS MOVED PAST THE NEAREST WALL, COMPUTE THE DISTANCE
     * TRAVELED PAST WALL - RADIUS MOVE PARTICLE TO THE OPPOSITE WALL AND TRAVEL
     * THE REMAINING DISTANCE
     */

    // x > 0, nearest x-wall > 0 the 'wall' includes the radius to make further
    // computations simpler
    if (x_curr * -1 < 0) {
        x_wall = wallBound - rad_temp;
        if (x_temp - x_wall > 0) {

            // distance from x trial position to nearest x-wall
            dist_xwall = dist_wall(x_temp, x_wall);

            // moves to opposite wall and travels the calculated distance
            x_temp = -1 * x_wall + dist_xwall;
        }
    } else {
        // x < 0, nearest x-wall < 0
        x_wall = -1 * (wallBound - rad_temp);

        if (x_temp - x_wall < 0) {
            dist_xwall = dist_wall(x_temp, x_wall);
            x_temp = -1 * x_wall - dist_xwall;
        }
    }

    if (y_curr * -1 < 0) {
        // y > 0, nearest y-wall > 0
        y_wall = wallBound - rad_temp;

        if (y_temp - y_wall > 0) {

            // distance from y trial position to the nearest y wall
            dist_ywall = dist_wall(y_temp, y_wall);
            y_temp = -1 * y_wall + dist_ywall;
        }
    } else {
        // y < 0, nearest y-wall < 0
        y_wall = -1 * (wallBound - rad_temp);

        if (y_temp - y_wall < 0) {

            // distance from y trial position to the nearest y wall
            dist_ywall = dist_wall(y_temp, y_wall);
            y_temp = -1 * y_wall - dist_ywall;
        }
    }

    current_prt.setX_TrialPos(x_temp); // puts the updated particle with updated
    current_prt.setY_TrialPos(y_temp); // trial positions back into the array

    (*particles)[index] = current_prt;
}

double Boundary::externalWell(std::vector<Particle> *particles, int index) {

    Particle current_prt;

    current_prt = (*particles)[index]; // assign current particle

    double x_temp = current_prt.getX_TrialPos(); // assign the current and trial
    double y_temp = current_prt.getY_TrialPos(); // positions and the radius of
                                                 // the current particle
    double x_curr = current_prt.getX_Position();
    double y_curr = current_prt.getY_Position();

    double energy_curr = ext_well_d * (pow(x_curr, 2) + pow(y_curr, 2));
    double energy_temp = ext_well_d * (pow(x_temp, 2) + pow(y_temp, 2));

    double delta_energy = energy_temp - energy_curr;

    return delta_energy;
}

void Boundary::initialPosition(std::vector<Particle> *particles,
                               KISSRNG randVal) {

    Particle current_prt;
    Particle compare_prt;

    double x_wall = 0;
    double y_wall = 0; // the locations of the nearest 'wall'

    bool accept = 0;

    for (int k = 0; k < n_particles; k++) {

        current_prt = (*particles)[k];
        double rad_temp = current_prt.getRadius();

        double x_temp = randVal.RandomUniformDbl() * 0.5 * boxLength;
        double y_temp = randVal.RandomUniformDbl() * 0.5 * boxLength;

        /* generate a random number from [0,1)
         * provides different 'quadrants' for the particle to be generated in
         * assigns the nearest x,y walls
         */

        double wallbound = 0.5 * boxLength;
        double num = randVal.RandomUniformDbl();

        if (k % 2 == 0 && num < .5) { // the acceptance presented in
            x_temp = -1 * x_temp; // 'check collision' is not implemented here
                                  // since the initial position cannot exceed 1
            x_wall = -1 * wallbound;
            y_wall = wallbound; // randomly assigns a negative/positive
        }                       // value to the generated position
        else if (k % 2 == 0 && num >= .5) {
            y_temp = -1 * y_temp;

            x_wall = wallbound;
            y_wall = -1 * wallbound;
        } else if (k % 2 != 0 && num < .5) {
            x_temp = -1 * x_temp;
            y_temp = -1 * y_temp;

            x_wall = -1 * wallbound;
            y_wall = -1 * wallbound;
        } else {
            x_wall = wallbound;
            y_wall = wallbound;
        }
        accept = 1;
        // the following statements can only change accept to 0
        if (abs(x_temp - x_wall) <
            rad_temp) { // if particle's center is closer than
            accept = 0; // one radius to the wall, reject move
        } else if (abs(y_temp - y_wall) < rad_temp) {
            accept = 0;
        } else if (k != 0) {
            for (int n = 0; n < k; n++) {
                // assign comparison particle
                compare_prt = (*particles)[n];

                // assign the comparison x,y position and radius
                double x_comp = compare_prt.getX_Position();
                double y_comp = compare_prt.getY_Position();
                double rad_comp = compare_prt.getRadius();

                // if the distance between particles is less than the sum of the
                // radii, reject position
                if (dist_part(x_temp, x_comp, y_temp, y_comp) <
                    rad_comp + rad_temp) {
                    accept = 0;
                    break;
                }
            }
        }
        // if position is accepted, assign the x,y position to current particle
        if (accept == 1) {
            current_prt.setX_Position(x_temp);
            current_prt.setY_Position(y_temp);

            // put initialized particle into array
            (*particles)[k] = current_prt;
        } else {
            k = k - 1; // generate a new random position for the same particle
        }
    }
}

int Boundary::initialHexagonal(std::vector<Particle> *particles) {

    Particle prt;

    int k = 0;
    int flag = 0;

    double x_dist = 0;
    double x_init_dist = 0;

    // the radius parameter can probably be removed
    double radius = 0; // use sigma if LJ is turned on
    int curr_row = 0;

    int return_num = n_particles;
    double h = 0.8660254038; // sin(pi/3)

    if (interact_type != 0) { // this should be true for WCA and LJ
        x_dist = sigma;
    } else {
        prt = (*particles)[0];
        radius = prt.getRadius(); // if LJ != 1, use particle diameter
        x_dist = 2 * radius;
    }

    double row1_num = boxLength / x_dist; // may be worth adding x_dist as
    double row2_num = (boxLength - 0.5 * x_dist) / (x_dist); // parameter

    // y_dist = dist between rows, y_init = y coor of the first row
    double y_dist = x_dist * h;
    double y_init_dist = 0.5 * boxLength - 0.5 * x_dist;

    while (k < n_particles) {

        /* ALTERNATE THE TWO TYPES OF ROWS
         * FOR EACH ROW, HAVE AN INITAL STARTING X POSITION
         *    THE INITIAL Y POSITION WILL BE INDEPENDENT OF
         *    THE ROW TYPE
         * SET THE INITIAL POSITION OF CURRENT PARTICLE
         * INCREASE COUNTER K TO ACCESS NEXT PARTICLE
         * CHANGE FLAG TO SWITCH TO THE OTHER TYPE OF ROW
         * UPDATE THE CURRENT Y POSITION OF NEXT ROW
         * IF INITIALIZATION IS OUTSIDE OF THE 'BOX', STOP
         * RETURN THE NUMBER OF PARTICLES SUCCESSFULLY
         *    INITIALIZED
         */

        if (flag == 0) {
            curr_row = row1_num;
            x_init_dist = -.5 * boxLength + 0.5 * x_dist; // the first row
        }                                                 // begins at the wall
        else {                                            // plus .5 unit x-dist
            curr_row = row2_num;
            x_init_dist = -.5 * boxLength + x_dist; // second row begins
        }                                           // at wall plus unit x_dist

        for (int n = 0; n < curr_row; n++) {
            prt = (*particles)[k];
            prt.setX_Position(x_init_dist + n * x_dist); // set x pos
            prt.setY_Position(y_init_dist);              // set y pos

            (*particles)[k] = prt;
            k++;

            if (k == n_particles) {
                break;
            }
        }

        flag = (flag == 0); // switch the row type
        y_init_dist = y_init_dist - y_dist;

        if (y_init_dist < -0.5 * boxLength + 0.5 * x_dist && k < n_particles) {
            return_num = k;
            break;
        }
    }

    // returns the number of particles successfully initialized
    return return_num;
}

int Boundary::initialSquare(std::vector<Particle> *particles) {
    Particle prt;

    int k = 0;
    int flag = 0;

    double x_dist = 0;
    double y_dist = 0;
    double radius = 0; // use sigma if LJ is turned on

    if (interact_type != 0) { // should be true for LJ and WCA
        x_dist = sigma;
        y_dist = sigma;
    } else {
        prt = (*particles)[0];

        x_dist = 2 * prt.getRadius();
        y_dist = 2 * prt.getRadius();
    }

    // determines the num of particles that fit into a single row
    double curr_row = boxLength / x_dist;
    double return_num = n_particles;

    double x_init_dist = -0.5 * boxLength + 0.5 * x_dist;
    double y_init_dist = 0.5 * boxLength - 0.5 * y_dist; // initial x,y position

    while (k < n_particles) {

        for (int n = 0; n < curr_row; n++) {
            prt = (*particles)[k];

            prt.setX_Position(x_init_dist + n * x_dist); // y pos stays the same
            prt.setY_Position(y_init_dist);              // in for-loop while
                                                         // x pos updates
            (*particles)[k] = prt;
            k++;

            if (k == n_particles) {
                break;
            }
        }
        // update y pos to the y_pos of the next row
        y_init_dist = y_init_dist - y_dist;

        if (y_init_dist < -0.5 * boxLength + 0.5 * y_dist && k < n_particles) {
            return_num = k;
            break;
        }
    }
    // returns num of particles successfully initialized
    return return_num;
}

void Boundary::initializeBoundary(Parameters *p) {

    boxLength = p->getBoxLength();
    sigma = p->getSigma();
    interact_type = p->getInteract_Type();
    n_particles = p->getNumParticles();
    ext_well_d = p->getExtWellDepth();
}
