#ifndef POTENTIAL_H
#define POTENTIAL_H
#include <string>

class potential{

private: 
	std::string type;
	std::string yamlFile;  

public: 
	potential(); 
	potential(std::string file); 

	double getProbability(); 	
	
	std::string getType(); 
	void setType(std::string x); 
};
#endif