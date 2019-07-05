#include <iostream>
#include <Particle.h>

Particle::Particle(int t, double rad, double w){
	type = t; 
	radius = rad; 
	stepWeight = w; 
}
Particle::Particle(){		// this may be unnecessary 
	type = 0; 
	radius = 0; 
	stepWeight = 0; 
}
/////////// GETTERS ////////////////////

int Particle::getType(){
	return type; 
} 
int Particle::getIdentifier(){
	return identifier; 
}

double Particle::getRadius(){
	return radius; 
}
double Particle::getX_Position(){
	return x_position; 
} 
double Particle::getY_Position(){
	return y_position; 
}
double Particle::getX_TrialPos(){
	return x_trialPos; 
}
double Particle::getY_TrialPos(){
	return y_trialPos; 
}

double Particle::getStepWeight(){
	return stepWeight; 
}

//////////// SETTERS //////////////////////

void Particle::setType(int t){
	type = t; 
}
void Particle::setIdentifier(int id){ //may be unnecessary 
	identifier = id; 
}

void Particle::setRadius(double rad){
	radius = rad; 
}
void Particle::setX_Position(double x){
	x_position = x; 
}
void Particle::setY_Position(double y){
	y_position = y; 
} 
void Particle::setX_TrialPos(double x){
	x_trialPos = x; 
}
void Particle::setY_TrialPos(double y){
	y_trialPos = y; 
}

void Particle::setStepWeight(double w){
	stepWeight = w; 
}

////////// OBJECT FUNCTIONS ////////////////

double Particle::x_trial(double randVal){
	return x_position + stepWeight * (randVal - 0.5); 
}
double Particle::y_trial(double randVal){
	return y_position + stepWeight * (randVal - 0.5); 
}
