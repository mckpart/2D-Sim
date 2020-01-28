#include <catch2/catch.hpp>
#include <yaml-cpp/yaml.h>

#include "../src/Parameters.h"
#include "../src/Particle.h"
#include "../src/Properties.h"

// come back and "un" hardcode the yaml file.. prob bad practice
Parameters *init_test_params() {
    static Parameters param;
    param.initializeParameters("params.yaml");
    return &param;
}

// write this force calculation function here first and then if it seems to be
// working, add it to the properties class

// the x and y here represent the distance between the particles, not the
// position of the reference particle in the simulation
void calc_force_vec(double r, double x, double y, double force,
                    std::vector<double> *F) {
    double Fx = -x / r * force;
    double Fy = -y / r * force;
    std::cout << "Fx = " << Fx << " Fy = " << Fy << std::endl;
    (*F)[0] = Fx;
    (*F)[1] = Fy;
}

void average_force_vec(std::vector<std::vector<double>> *Force) {
    std::vector<double> f(2);
    double fx = 0;
    double fy = 0;
    //    std::cout << Force->size() << std::endl;
    for (int i = 0; i < Force->size(); ++i) {
        fx = fx + (*Force)[i][0];
        fy = fy + (*Force)[i][1];
    }
    // for the LJ model, when the particle at (-.5,.5) is the refernce particle,
    // the total force vector should be 22.875(-1,1)
    std::cout << "final x and y: " << fx << " " << fy << std::endl;
}
// TEST_CASE("does func return 1") { REQUIRE(func_1() == 2); }

TEST_CASE("Lennard Jones Force Calculation") {
    Properties prop;
    Parameters *param;

    param = init_test_params();
    prop.initializeProperties(param);
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

    std::vector<std::vector<double>> force_vecs(3, std::vector<double>(2));
    std::vector<Particle> particles = {p1, p2, p3, p4};
    for (int i = 0; i < 3; ++i) {
        Particle p = particles[i + 1];
        //        std::cout << "x of p1: " << p1.getX_Position()
        //                  << " x of p: " << p.getX_Position() << std::endl;
        double r = prop.radDistance(p1.getX_Position(), p.getX_Position(),
                                    p1.getY_Position(), p.getY_Position());
        std::cout << "rad dist is: " << r << std::endl;
        // not sure what the c = 1 was originally for...
        double F = prop.lenJonesForce(r, 1);
        std::cout << "force is = " << F << std::endl;

        calc_force_vec(r, p.getX_Position() - p1.getX_Position(),
                       p.getY_Position() - p1.getY_Position(), F,
                       &force_vecs[i]);
    }
    for (int i = 0; i < 3; ++i) {

        std::cout << force_vecs[i][0] << "  " << force_vecs[i][1] << std::endl;
    }

    average_force_vec(&force_vecs);
    std::cout << "here it is" << param->getSigma() << std::endl;
}

