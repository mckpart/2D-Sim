#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>

class Particle {

  private:
    int type = 0; // may be worth changing type to a string later on
    int identifier = 0;

    double radius = 0;
    double x_position = 0;
    double y_position = 0;

    double x_trialPos = 0;
    double y_trialPos = 0;

    double stepWeight = 0;

    std::vector<double> force;

  public:
    // DEFAULT CONSTRUCTOR //
    Particle();

    // GETTERS //
    int getType();
    int getIdentifier();

    double getRadius();
    double getX_Position();
    double getY_Position();
    double getX_TrialPos();
    double getY_TrialPos();

    double getX_Force();
    double getY_Force();

    // SETTERS //
    double getStepWeight();

    void setType(int t);
    void setIdentifier(int id);

    void setRadius(double rad);
    void setX_Position(double x);
    void setY_Position(double y);
    void setX_TrialPos(double x);
    void setY_TrialPos(double y);

    void setStepWeight(double w);

    void addForce(double fx, double fy);
    void resetForce();
    std::vector<double> *getForce();

    double x_trial(double randVal);
    double y_trial(double randVal);

};
#endif
