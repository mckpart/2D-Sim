#include <catch2/catch.hpp>
#include <yaml-cpp/yaml.h>

#include "../src/Parameters.h"
#include "../src/Particle.h"
#include "../src/Properties.h"

// come back and "un" hardcode the yaml file.. prob bad practice
Parameters *init_test_params() {
    static Parameters param;
    param.initializeParameters("catch_testing/test_params.yaml");
    return &param;
}

void prepare_particles(std::vector<Particle> *p) {

    Particle p1;
    Particle p2;
    Particle p3;
    Particle p4;

    p1.setX_Position(-.5);
    p1.setY_Position(.5);
    p2.setX_Position(.5);
    p2.setY_Position(.5);
    p3.setX_Position(.5);
    p3.setY_Position(-.5);
    p4.setX_Position(-.5);
    p4.setY_Position(-.5);

    *p = {p1, p2, p3, p4};
}

// the x and y here represent the distance between the particles, not the
// position of the reference particle in the simulation
void calc_force_vec(double r, double x, double y, double force,
                    std::vector<double> *F) {
    double Fx = -x / r * force;
    double Fy = -y / r * force;
    std::cout << "Fx = " << Fx << " Fy = " << Fy << std::endl;
    (*F)[0] = Fx + (*F)[0];
    (*F)[1] = Fy + (*F)[1];
}

void average_force_vec(std::vector<std::vector<double>> *Force) {
    std::vector<double> f(2);
    double fx = 0;
    double fy = 0;
    for (int i = 0; i < Force->size(); ++i) {
        fx += (*Force)[i][0];
        fy += (*Force)[i][1];
    }
    // for the LJ model, when the particle at (-.5,.5) is the refernce particle,
    // the total force vector should be 22.875(-1,1)
    std::cout << "final x and y: " << fx << " " << fy << std::endl;
}

void compute_f_vecs(std::vector<double> *tot_f, int k, Properties *prop) {

    std::vector<Particle> particles;
    prepare_particles(&particles);
    Particle p_reff = particles[k];

    for (int i = 0; i < 4; ++i) {
        if (i != k) {
            Particle p = particles[i];
            double r =
                prop->radDistance(p_reff.getX_Position(), p.getX_Position(),
                                  p_reff.getY_Position(), p.getY_Position());
            // not sure what the c = 1 was originally for...
            double F = prop->lenJonesForce(r, 1);
            prop->calc_force_vec(p.getX_Position() - p_reff.getX_Position(),
                                 p.getY_Position() - p_reff.getY_Position(), r,
                                 tot_f);
        }
    }
}

// VERY IMPORTANT: THIS TEST CASE WORKS WHEN SIGMA IS 1! MAYBE MAKE A SIMPLE
// PARAMETER FILE TO BE READ IN BY THIS TEST FUNCTIONS!!!!!!!!!!!
TEST_CASE("Lennard Jones Force Calculation") {
    Properties prop;
    Parameters *param;

    param = init_test_params();
    prop.initializeProperties(param);

    std::vector<std::vector<double>> force_vecs(4, std::vector<double>(2));
    std::vector<double> tot_force_vec(2);

    compute_f_vecs(&tot_force_vec, 0, &prop);
    REQUIRE(tot_force_vec[0] == -22.875);
    REQUIRE(tot_force_vec[1] == 22.875);
    std::cout << "LJ FORCE PASSED TEST" << std::endl;

    for (int k = 0; k < 4; ++k) {
        compute_f_vecs(&force_vecs[k], k, &prop);
    }
    for (int k = 0; k < 4; ++k) {
        std::cout << force_vecs[k][0] << " and " << force_vecs[k][1]
                  << std::endl;
    }
}
