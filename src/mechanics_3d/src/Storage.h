
#ifndef STORAGE_H_
#define STORAGE_H_

#include <fstream>	     
#include <memory>
#include <iomanip>

class System;
class SystemBuilder;

//During graph deformation, this file stores position and velocity of nodes at a given time step
class Storage {
	
	std::weak_ptr<System> system;
	std::weak_ptr<SystemBuilder> builder;
	//std::shared_ptr<ExternalForce> grip;
	std::ofstream output;
	std::ofstream statesOutput;
	
	std::ofstream statesOutputStrain;
	std::string bn;

	unsigned stepCounter = 0;
	unsigned stepInterval = 10;

	double forceStep = 0.0;	 //increments in which force will be increased. 
	double magnitudeForce = 0.0;  //how much force we are currently at.
	int currentAddedEdges = 0;
	int previousAddedEdges = 0;
	unsigned iteration = 0;
	std::string str_animation;
	std::string str_state;
	std::string str_params;

public: 
	Storage(std::weak_ptr<System> a_system,
		std::weak_ptr<SystemBuilder> b_system);

	void save_current_state(void);
	void save_params(void);

	void updateStorage(void);
	void print_VTK_File(void);
};
#endif