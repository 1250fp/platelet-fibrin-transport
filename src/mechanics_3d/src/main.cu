

#include <iomanip>
#include <string>
#include <memory>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <inttypes.h>
#include <cstddef>
#include <sstream> 

#include "pugixml.hpp"

#include "System.h"

#include "System_Builder.h"
#include "Storage.h"


void advanceSystem(std::string filename, std::shared_ptr<System> system){
	std::ifstream ifs(filename.c_str());
	std::string temp;
	std::string line;

	
	if(!ifs) {
		std::cout << filename << " is not available" << std::endl;
		return;
	}

	unsigned node_increment = 0;
	unsigned plt_increment = 0;
	unsigned temp_edge_counter = 0;
	unsigned len_increment = 0;
	unsigned gnbr_increment = 0;
	unsigned is_node_incr = 0;
	unsigned fil_node_id_incr = 0;
	unsigned fil_node_type_incr = 0;
	unsigned filopodia_count_per_plt = 0;
		
	while (std::getline(ifs,line)) {
		

		std::stringstream ss;
		std::string temp;
		
		ss.str(line);
		std::getline(ss,temp,' ');

		if (temp == "time") {
			std::getline(ss,temp,' ');
			double time = std::atof(temp.c_str());
			system->generalParams.currentTime = time;
			std::cout<<"resetting time: " << time << std::endl;
			std::getline(ss,temp,'\n');
		}
		else if (temp == "node") {
			std::getline(ss,temp,' ');
			double x = std::atof(temp.c_str());
			std::getline(ss,temp,' ');
			double y = std::atof(temp.c_str());
			std::getline(ss,temp,' ');
			double z = std::atof(temp.c_str());

			system->nodeInfoVecs.nodeLocX[node_increment] = x;
			system->nodeInfoVecs.nodeLocY[node_increment] = y;
			system->nodeInfoVecs.nodeLocZ[node_increment] = z;
			node_increment+=1;
		}
		else if (temp == "plt") {
			std::getline(ss,temp,' ');
			double x = std::atof(temp.c_str());
			std::getline(ss,temp,' ');
			double y = std::atof(temp.c_str());
			std::getline(ss,temp,' ');
			double z = std::atof(temp.c_str());

			system->pltInfoVecs.pltLocX[plt_increment] = x;
			system->pltInfoVecs.pltLocY[plt_increment] = y;
			system->pltInfoVecs.pltLocZ[plt_increment] = z;
			plt_increment+=1;
		}
		else if (temp == "edge_count") {
			std::getline(ss,temp,' ');
			unsigned edge_count = std::stoi(temp.c_str());
			system->generalParams.currentEdgeCount = edge_count;
			std::cout<<"RESET EDGE COUNT: " << edge_count << std::endl;
		}
		else if (temp == "edge_lr") {			
			std::getline(ss,temp,' ');
			unsigned left_edge = std::stoi(temp.c_str());
			std::getline(ss,temp,' ');
			unsigned right_edge = std::stoi(temp.c_str());

			system->nodeInfoVecs.hostEdgeLeft[temp_edge_counter] = left_edge;
			system->nodeInfoVecs.hostEdgeRight[temp_edge_counter] = right_edge;
			temp_edge_counter+=1;
		}
		//lengths must be overwritten because added edges have lenzero
		else if (temp == "length_zero") {
			std::getline(ss,temp,' ');
			double len_zero = std::atof(temp.c_str());
			system->wlcInfoVecs.lengthZero[len_increment] = len_zero;
			len_increment+=1;
		}
		else if (temp == "global_nbr") {
			std::getline(ss,temp,' ');
			unsigned nbr_index = std::stoi(temp.c_str());
			system->wlcInfoVecs.globalNeighbors[gnbr_increment] = nbr_index;
			
			//check corresponding length of edge, if not filled use 0.1
			double test_length = system->wlcInfoVecs.lengthZero[gnbr_increment];
			if (test_length == 0.0 ){
				system->wlcInfoVecs.lengthZero[gnbr_increment] = system->generalParams.fiberDiameter;
			}
			gnbr_increment+=1;

		}
		else if (temp == "is_node_in_plt") {
			std::getline(ss,temp,' ');
			bool is_in_plt = std::stoi(temp.c_str());
			system->nodeInfoVecs.isNodeInPltVol[is_node_incr] = is_in_plt;
			is_node_incr+=1;
		}
		else if (temp == "is_node_in_plt") {
			std::getline(ss,temp,' ');
			unsigned fil_node_id = std::stoi(temp.c_str());
			system->pltInfoVecs.tndrlNodeId[fil_node_id_incr] = fil_node_id;
			fil_node_id_incr+=1;
		}		
		else if (temp == "is_node_in_plt") {
			std::getline(ss,temp,' ');
			unsigned fil_node_type = std::stoi(temp.c_str());
			system->pltInfoVecs.tndrlNodeType[fil_node_type_incr] = fil_node_type;
			fil_node_type_incr+=1;
		}
		else if (temp == "filopodia_count") {
			std::getline(ss,temp,' ');
			unsigned num_filopodia = std::stoi(temp.c_str());
			system->pltInfoVecs.tndrlNodeNum[filopodia_count_per_plt] = num_filopodia;
			filopodia_count_per_plt+=1;
		}
		//else {std::cout<<"INPUT UNKNOWN"<< std::endl;}

	}
	std::cout<<"finished state copy" << std::endl;
};

std::shared_ptr<System> createSystem(const char* schemeFile, std::shared_ptr<SystemBuilder> builder)	{
	pugi::xml_document doc;
	pugi::xml_parse_result parseResult = doc.load_file(schemeFile);

	if (!parseResult) {
		std::cout << "parse error in createNodeSystem: " << parseResult.description() << std::endl;
		return nullptr;
	}
	pugi::xml_node root = doc.child("data");
	pugi::xml_node nodes = root.child("nodes");
	pugi::xml_node plts = root.child("plts");
	pugi::xml_node links = root.child("links");
	pugi::xml_node props = root.child("settings");

	//first, we'll input settings
	if (!(root && nodes && links)) {
		std::cout << "couldn't find nessesary data\n";
		//return false;
	}

	//default settings in NSB.h
	double base_tor = 0.00226259342068;
	double base_diam = 0.1;
	int base_num_mon=1100;
	if (auto p = props.child("link-diameter")){//base diameter 0.1 microns
		builder->defaultLinkDiameter = (p.text().as_double());
		double r_ratio= (p.text().as_double())/base_diam;
		builder->defaultTorsionSpringStiffness = pow(r_ratio,2.4)*base_tor;//E prop to r^-1.6 and I prop to r^4
		builder->defaultNumMonFiberArea = base_num_mon*pow(r_ratio,2);// N prop to r^2
	}
	else{
		builder->defaultLinkDiameter = base_diam;
		builder->defaultTorsionSpringStiffness = base_tor;
		builder->defaultNumMonFiberArea = base_num_mon;
	}

	if (auto p = props.child("resistance_fibrin"))
		builder->viscousDamp_Fibrin = (p.text().as_double());

	if (auto p = props.child("resistance_plt"))
		builder->viscousDamp_Plt = (p.text().as_double());

	if (auto p = props.child("spring-stiffness"))
		builder->defaultSpringStiffness = (p.text().as_double());

	/* if (auto p = props.child("torsion-stiffness"))
		builder->defaultTorsionSpringStiffness = (p.text().as_double());*/

	if (auto p = props.child("persistance-length"))
		builder->defaultPersistanceLength = (p.text().as_double());

	if (auto p = props.child("absolute-temperature"))
		builder->defaultTemperature = (p.text().as_double());

	if (auto p = props.child("contour-length-multiplier"))
		builder->defaultContourLengthMultiplier = (p.text().as_double());

	if (auto p = props.child("units-per-extra-node"))
		builder->defaultUnitsPerExtraNode = (p.text().as_double());

	if (auto p = props.child("extra-nodes-per-edge"))
		builder->defaultExtraNodesPerEdge = (p.text().as_uint());

	if (auto p = props.child("use-extra-nodes"))
		builder->useExtraNodes = (p.text().as_bool());

	if (auto p = props.child("constant-extra-nodes"))
		builder->useConstantNumberOfExtraNodes = (p.text().as_bool());

	if (auto p = props.child("worm-like"))
		builder->wormLikeEnabled = (p.text().as_bool());

	if (auto p = props.child("use-linking"))
		builder->linking = (p.text().as_bool());
////////////////////////////////////////////////////
//platelets parameters
	if (auto p = props.child("plt_mass"))
		builder->defaultPltMass = (p.text().as_double());

	if (auto p = props.child("plt_force"))
		builder->pltForce = (p.text().as_double());

	if (auto p = props.child("plt_other_intrct"))
		builder->plt_other_intrct = (p.text().as_uint());

	if (auto p = props.child("plt_tndrl_intrct"))
		builder->plt_tndrl_intrct = (p.text().as_uint());

	if (auto p = props.child("tndrl_M"))
		builder->tndrl_M = (p.text().as_double());
	if (auto p = props.child("tndrl_S"))
		builder->tndrl_S = (p.text().as_double());


	if (auto p = props.child("plt_r"))
		builder->pltR = (p.text().as_double());

	if (auto p = props.child("plt_r_force")) {
		builder->pltRForce = (p.text().as_double());
	}
	//fraction of platelet radius that adhesion is in.
	if (auto p = props.child("plt_r_adhesion")) {
		double RAdhesion=(p.text().as_double());
		if (RAdhesion>0.0 && RAdhesion<1.0){
			builder->pltRAdhesion = RAdhesion;
		}
		else{
			std::cout << "parse error: platelet adhesion radius fraction is not valid\n";
			return 0;
		}
	}


	if (auto p = props.child("plt_density")) {
		builder->pltDensity = (p.text().as_double());
		std::cout << "setting density: " << builder->pltDensity << std::endl;
	}

	if (auto p = props.child("use-pltforcefield")){
		builder->pltfrcfld = (p.text().as_bool());
		std::cout<<"frcFld: "<< builder->pltfrcfld<<std::endl;
	}
	if (auto p = props.child("use-plttndrl")){
		builder->plttndrl = (p.text().as_bool());
		std::cout<<"plttndrl: "<< builder->plttndrl<<std::endl;
	}
	if (auto p = props.child("use-pltrelease")){
		builder->pltrelease = (p.text().as_bool());
		std::cout<<"pltrelease: "<< builder->pltrelease<<std::endl;
	}
	if (auto p = props.child("use-plthandhand")){
		builder->plthandhand = (p.text().as_bool());
		std::cout<<"plthandhand: "<< builder->plthandhand<<std::endl;
	}
	if (auto p = props.child("strain_switch")){
		builder->strainswitch = (p.text().as_double());
		std::cout<<"strainswitch: "<< builder->strainswitch<<std::endl;
	}

	if (auto p = props.child("use-pltonplt")) {
		builder->pltonplt = (p.text().as_bool());
		std::cout << "plt_interact: " << builder->pltonplt << std::endl;
	}

	if (auto p = props.child("dynamic-plt")) {
		builder->use_dynamic_plt_force = (p.text().as_bool());
		std::cout << "dynamic_plt: " << builder->use_dynamic_plt_force << std::endl;
	}

	if (auto p = props.child("dynamic-plt-max-force")) {
		builder->max_dynamic_plt_force = (p.text().as_double());
		std::cout << "dynamic_plt_max_force: " << builder->max_dynamic_plt_force << std::endl;
	}

	if (auto p = props.child("use_nonlinear_dynamic_force")) {
		builder->use_nonlinear_dynamic_force = (p.text().as_bool());
		//std::cout << "use_nonlinear_dynamic_force: " << builder->use_nonlinear_dynamic_force << std::end;
	}
	if (auto p = props.child("distribute_plt_force")) {
		builder->distribute_plt_force = (p.text().as_bool());
	}
	if (auto p = props.child("tndrl_plt_pdf")) {
		builder->tndrl_plt_pdf = (p.text().as_bool());
	}

	std::cout << "builder ptr address: " << builder << std::endl;
	std::vector<unsigned> originNodes;
//buid nodes
	double mass;
	double x, y, z; //variables to be used reading in data.
	double defaultMass = nodes.attribute("default-mass").as_double(-1.0);
	builder->defaultMass = defaultMass;

	//debuggin printl
	//std::cout<<"default mass: " << defaultMass << std::endl;
	for (auto node = nodes.child("node"); node; node = node.next_sibling("node")) {
		mass = node.attribute("mass").as_double(defaultMass);
		if (mass < 0.0) {
			std::cout << "parse error: node mass is undefined\n";
			return 0;
		}
		const char* text = node.text().as_string();

		//std::cout<<"mass: " << mass << std::endl;
		if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
			std::cout << "parse node error\n";
			return 0;
		}
		int unused = builder->addNode(mass, glm::dvec3(x, y, z));
	}

	unsigned from, to;
	for (auto link = links.child("link"); link; link = link.next_sibling("link")) {
		if (2 != sscanf(link.text().as_string(""), "%u %u" , &from, &to)) {
			std::cout << "parse link error\n";
			return 0;
		}
		builder->putSpring(from, to); //adds edges into saved vectors

	}

	std::cout << "post springs" << std::endl;


	//replace false entries with true if node is fixed.
	pugi::xml_node fixedRoot = root.child("fixed");
	if (fixedRoot) {
		for (auto node = fixedRoot.child("node"); node; node = node.next_sibling("node"))
			builder->fixNode(node.text().as_uint());
	}

	std::cout << "post fixed" << std::endl;

//build platelets

	double pltmass;
	double pltx, plty, pltz; //variables to be used reading in data.

	//only use platelet input if density is zero
	for (auto plt = plts.child("plt"); plt; plt = plt.next_sibling("plt")) {

		const char* text = plt.text().as_string();

		//std::cout<<"mass: " << mass << std::endl;
		if (3 != sscanf(text, "%lf %lf %lf", &pltx, &plty, &pltz)) {
			std::cout << "parse plt error\n";
			return 0;
		}
		int unused = builder->addPlt(builder->defaultPltMass, glm::dvec3(pltx, plty, pltz));
	}




	//replace false entries with true if node is fixed.
	pugi::xml_node fixedPltRoot = root.child("pltfixed");
	if (fixedPltRoot) {
		for (auto plt = fixedPltRoot.child("plt"); plt; plt = plt.next_sibling("plt"))
			builder->fixPlt(plt.text().as_uint());
	}

	std::cout << "post fixed" << std::endl;

	//last, set and add non resistance and non external constraints.
	auto model = builder->create();

	std::cout << "model built" << "\n";

	return model;
}

void run(int argc, char** argv) {

	time_t t0,t1;
	t0 = time(0);

	double epsilon = 0.01;
	double timeStep = 0.001;

	for (int i = 0; i < argc; i++) {

		std::string arg = argv[i];
		int pos = arg.find('=');

		std::string key = arg.substr(0, pos);
		std::string val = arg.substr(pos + 1);

		std::cout<<"argc: "<< argc <<std::endl;
		std::cout<<"arg: "<< arg <<std::endl;
		std::cout<<"pos: "<< pos <<std::endl;
		std::cout<<"key: "<< key <<std::endl;
		std::cout<<"val: "<< val <<std::endl;

		if (key == "-dt") {
			timeStep = std::atof(val.c_str());

			std::cout<<"setting timestep: "<< timeStep << std::endl;
			continue;
		}
		if (key == "-eps") {
			epsilon = std::atof(val.c_str());
			std::cout<<"setting epsilon: "<< epsilon << std::endl;
			continue;
		}
	}


	auto builder = std::make_shared<SystemBuilder>(epsilon, timeStep);

		auto system = createSystem(argv[argc-1], builder);

	//once the system is set, we'll store the initial values via the ptr system.
	//Storage storage( system, outputFileName);
	auto storage = std::make_shared<Storage>(system, builder);


	std::cout << "assigning fdiagram in main" << std::endl;
	system->assignStorage(storage);
	std::cout << "post fdiagram in main" << std::endl;

	std::cout << "solving system in main" << std::endl;
	system->solveSystem();


	t1 = time(0);  //current time at the end of solving the system.
	int total,hours,min,sec;
	total = difftime(t1,t0);
	hours = total / 3600;
	min = (total % 3600) / 60;
	sec = (total % 3600) % 60;
	std::cout << "Total time hh: " << hours << " mm:" << min << " ss:" << sec <<"\n";

};

void run_time_advance(int argc, char** argv) {

	double epsilon = 0.01;
	double timeStep = 0.001;

	for (int i = 0; i < argc; i++) {

		std::string arg = argv[i];
		int pos = arg.find('=');

		std::string key = arg.substr(0, pos);
		std::string val = arg.substr(pos + 1);

		std::cout<<"argc: "<< argc <<std::endl;
		std::cout<<"arg: "<< arg <<std::endl;
		std::cout<<"pos: "<< pos <<std::endl;
		std::cout<<"key: "<< key <<std::endl;
		std::cout<<"val: "<< val <<std::endl;

		if (key == "-dt") {
			timeStep = std::atof(val.c_str());

			std::cout<<"setting timestep: "<< timeStep << std::endl;
			continue;
		}
		if (key == "-eps") {
			epsilon = std::atof(val.c_str());
			std::cout<<"setting epsilon: "<< epsilon << std::endl;
			continue;
		}
	}


	auto builder = std::make_shared<SystemBuilder>(epsilon, timeStep);
	auto system = createSystem(argv[argc-2], builder);
	auto storage = std::make_shared<Storage>(system, builder);


	std::cout << "assigning fdiagram in main" << std::endl;
	system->assignStorage(storage);
	std::cout << "resetting time" << std::endl;
	advanceSystem(argv[argc-1], system);
		
	std::cout << "time reset" << std::endl;
	std::cout << "solving system in main" << std::endl;
	system->solveSystem();
};


int main(int argc, char** argv)
{
	std::cout << argc << std::endl;
	std::string mode = argv[1];
	if (mode=="reset"){
		
		std::cout << "running reset mode" << std::endl;
		run_time_advance(argc,argv);
	}else{ run(argc, argv); }
	
	

	return 0;
}
