#include <sys/stat.h>

#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "System.h"
#include "System_Builder.h"
#include "SystemStructures.h"
#include "Storage.h"


Storage::Storage(std::weak_ptr<System> a_system,
	std::weak_ptr<SystemBuilder> b_system ) {
	//std::cout << "FDM constructor" << std::endl;

	system = a_system;
	builder = b_system;

	std::shared_ptr<System> sys = system.lock();
	if (sys) { 

		std::stringstream stream_min;
		std::stringstream stream_max;

		unsigned domain_size = ceil((sys->domainParams.maxX + 
			sys->domainParams.maxY + 
			sys->domainParams.maxZ) / 3.0);

		unsigned plt_count = sys->generalParams.maxPltCount;
		unsigned num_filo = sys->generalParams.plt_tndrl_intrct;
		double max_force = sys->generalParams.max_dynamic_plt_force;
		double min_force = sys->generalParams.pltForce;
		bool response_nonlin = sys->generalParams.use_nonlinear_dynamic_force;
		bool response = sys->generalParams.use_dynamic_plt_force;
		bool force_scale = sys->generalParams.distribute_plt_force;
		
		stream_min << std::fixed << std::setprecision(2) << min_force;
		std::string str_min_force = stream_min.str();
		
		stream_max << std::fixed << std::setprecision(2) << max_force;
		std::string str_max_force = stream_max.str();

		std::string str_force_scale = "_pltForceScale_";
		std::string str_nonlinear_response = "_nonlinDynamicResp_";
		std::string str_dynamic_response = "_dynamicResp_";
		std::string str_domain = "_domain_";
		std::string str_plt_count = "_plt_count_";
		std::string str_filo_count = "_filo_count_";
		std::string str_maxF = "_maxForce_";
		std::string str_minF = "_minForce_";

		std::string str_a = "Animation_";
		std::string str_p = "Params_";
		std::string str_s = "State_";

		str_animation = str_a +
			str_domain + std::to_string(domain_size) +
			str_plt_count + std::to_string(plt_count) + 
			str_filo_count + std::to_string(num_filo) + 
			str_maxF + str_max_force + 
			str_minF + str_min_force +
			str_force_scale + std::to_string(force_scale) +
			str_dynamic_response + std::to_string(response) +
			str_nonlinear_response + std::to_string(response_nonlin);

		std::cout<<"domain: " << domain_size<<
			" str_plt_count: " << plt_count<<
			" str_filo_count: " << num_filo<<
			" str_minF: " << min_force<< 
			" str_response: " << response << std::endl;


		const int dir_err_anim = mkdir(str_animation.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_anim)
		{
			printf("Error creating directory animation test!n");
			//exit(1);
		}
		else {
			printf("making folder!n");
			printf(str_animation.c_str());
		}

		str_params = str_p + 
			str_domain + std::to_string(domain_size) +
			str_plt_count + std::to_string(plt_count) + 
			str_filo_count + std::to_string(num_filo) + 
			str_maxF + str_max_force + 
			str_minF + str_min_force +
			str_force_scale + std::to_string(force_scale) +
			str_dynamic_response + std::to_string(response) +
			str_nonlinear_response + std::to_string(response_nonlin);

		const int dir_err_params = mkdir(str_params.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_params)
		{
			printf("Error creating directory params!n");
			//exit(1);
		}
		else {
			printf("making folder!n");
			printf(str_params.c_str());
		}

		str_state = str_s + 
			str_domain + std::to_string(domain_size) +
			str_plt_count + std::to_string(plt_count) + 
			str_filo_count + std::to_string(num_filo) + 
			str_maxF + str_max_force + 
			str_minF + str_min_force +
			str_force_scale + std::to_string(force_scale) +
			str_dynamic_response + std::to_string(response) +
			str_nonlinear_response + std::to_string(response_nonlin);

		const int dir_err_state = mkdir(str_state.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_state)
		{
			printf("Error creating directory state!n");
			//exit(1);
		}
		else {
			printf("making folder!n");
			printf(str_state.c_str());
		}
	}
};

void Storage::save_current_state(){
	std::shared_ptr<System> sys = system.lock();
	if (!sys) return;
	// Read-only export snapshots: avoid a separate GPU transfer for every value.
	thrust::host_vector<bool> io_snapshot_nodeInfoVecs_isNodeInPltVol = sys->nodeInfoVecs.isNodeInPltVol;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_nodeLocX = sys->nodeInfoVecs.nodeLocX;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_nodeLocY = sys->nodeInfoVecs.nodeLocY;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_nodeLocZ = sys->nodeInfoVecs.nodeLocZ;
	thrust::host_vector<double> io_snapshot_pltInfoVecs_pltLocX = sys->pltInfoVecs.pltLocX;
	thrust::host_vector<double> io_snapshot_pltInfoVecs_pltLocY = sys->pltInfoVecs.pltLocY;
	thrust::host_vector<double> io_snapshot_pltInfoVecs_pltLocZ = sys->pltInfoVecs.pltLocZ;
	thrust::host_vector<unsigned> io_snapshot_pltInfoVecs_tndrlNodeId = sys->pltInfoVecs.tndrlNodeId;
	thrust::host_vector<unsigned> io_snapshot_pltInfoVecs_tndrlNodeNum = sys->pltInfoVecs.tndrlNodeNum;
	thrust::host_vector<unsigned> io_snapshot_pltInfoVecs_tndrlNodeType = sys->pltInfoVecs.tndrlNodeType;
	thrust::host_vector<unsigned> io_snapshot_wlcInfoVecs_globalNeighbors = sys->wlcInfoVecs.globalNeighbors;
	thrust::host_vector<double> io_snapshot_wlcInfoVecs_lengthZero = sys->wlcInfoVecs.lengthZero;

	if (sys) {

		//first create a new file using the current network strain
		
		std::string format = ".sta";
		
		std::string strain =  std::to_string(sys->generalParams.currentTime);
		std::string initial = str_state+"/State_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());



		//unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		
		//place nodes
		ofs << std::setprecision(30) <<std::fixed<< "time " << sys->generalParams.currentTime<<'\n';
		
		//thrust::copy(sys->nodeInfoVecs.nodeLocX.begin(),
		//	sys->nodeInfoVecs.nodeLocX.end(), hostLocX.begin());
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			
			double x = io_snapshot_nodeInfoVecs_nodeLocX[i];
			double y = io_snapshot_nodeInfoVecs_nodeLocY[i];
			double z = io_snapshot_nodeInfoVecs_nodeLocZ[i];
			ofs << std::setprecision(30) <<std::fixed<< "node " << x << " " << y << " " << z <<'\n';
		}
		
		//place plts
		for (unsigned i = 0; i < sys->pltInfoVecs.pltLocX.size(); i++) {
			double x = io_snapshot_pltInfoVecs_pltLocX[i];
			double y = io_snapshot_pltInfoVecs_pltLocY[i];
			double z = io_snapshot_pltInfoVecs_pltLocZ[i];
			ofs << std::setprecision(30) <<std::fixed<< "plt " << x << " " << y << " " << z <<'\n';
		
		}
		
		//
		ofs << std::setprecision(30) <<std::fixed<< "edge_count " << sys->generalParams.currentEdgeCount<<'\n';
		for (unsigned i = 0; i < sys->generalParams.currentEdgeCount; i++) {
			unsigned left = sys->nodeInfoVecs.hostEdgeLeft[i];
			unsigned right = sys->nodeInfoVecs.hostEdgeRight[i];

			ofs << std::setprecision(30) <<std::fixed<< "edge_lr " << left << " " << right <<'\n';
		}

		/*
		ofs << std::setprecision(30) <<std::fixed<< "minX " << sys->domainParams.originMinX<<'\n';
		ofs << std::setprecision(30) <<std::fixed<< "maxX " << sys->domainParams.originMaxX<<'\n';
		ofs << std::setprecision(30) <<std::fixed<< "minY " << sys->domainParams.originMinY<<'\n';
		ofs << std::setprecision(30) <<std::fixed<< "maxY " << sys->domainParams.originMaxY<<'\n';
		ofs << std::setprecision(30) <<std::fixed<< "minZ " << sys->domainParams.originMinZ<<'\n';
		ofs << std::setprecision(30) <<std::fixed<< "maxZ " << sys->domainParams.originMaxZ<<'\n';
		*/
		for (unsigned i = 0; i < sys->wlcInfoVecs.globalNeighbors.size(); i++) {
			unsigned nbr = io_snapshot_wlcInfoVecs_globalNeighbors[i];

			ofs << std::setprecision(30) <<std::fixed<< "global_nbr " << nbr << '\n';
		}

		for (unsigned i = 0; i < sys->wlcInfoVecs.lengthZero.size(); i++) {
			double len_zero = io_snapshot_wlcInfoVecs_lengthZero[i];

			ofs << std::setprecision(30) <<std::fixed<< "length_zero " << len_zero << '\n';
		}

		for (unsigned i = 0; i < sys->nodeInfoVecs.isNodeInPltVol.size(); i++) {
			bool is_in_plt = io_snapshot_nodeInfoVecs_isNodeInPltVol[i];

			ofs << std::setprecision(3) <<std::fixed<< "is_node_in_plt " << is_in_plt << '\n';
		}
		for (unsigned i = 0; i < sys->pltInfoVecs.tndrlNodeId.size(); i++) {
			unsigned fil_node_id = io_snapshot_pltInfoVecs_tndrlNodeId[i];

			ofs << std::setprecision(3) <<std::fixed<< "file_node_id " << fil_node_id << '\n';
		}
		for (unsigned i = 0; i < sys->pltInfoVecs.tndrlNodeType.size(); i++) {
			unsigned fil_node_type = io_snapshot_pltInfoVecs_tndrlNodeType[i];
			ofs << std::setprecision(3) <<std::fixed<< "fil_node_type " << fil_node_type << '\n';
		}

		for (unsigned i = 0; i < sys->pltInfoVecs.tndrlNodeNum.size(); i++) {
			unsigned num_filopodia = io_snapshot_pltInfoVecs_tndrlNodeNum[i];
			ofs << std::setprecision(3) <<std::fixed<< "filopodia_count " << num_filopodia << '\n';
		}
		
	}
}

void Storage::save_params(void) {
	std::shared_ptr<System> sys = system.lock();
	if (!sys) return;
	// Read-only export snapshots: avoid a separate GPU transfer for every value.
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_discretizedEdgeAlignment = sys->nodeInfoVecs.discretizedEdgeAlignment;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_discretizedEdgeStrain = sys->nodeInfoVecs.discretizedEdgeStrain;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_nodeForceX = sys->nodeInfoVecs.nodeForceX;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_nodeForceY = sys->nodeInfoVecs.nodeForceY;
	thrust::host_vector<double> io_snapshot_nodeInfoVecs_nodeForceZ = sys->nodeInfoVecs.nodeForceZ;
	thrust::host_vector<unsigned> io_snapshot_pltInfoVecs_nodeImagingConnection = sys->pltInfoVecs.nodeImagingConnection;
	thrust::host_vector<unsigned> io_snapshot_pltInfoVecs_pltImagingConnection = sys->pltInfoVecs.pltImagingConnection;
	thrust::host_vector<unsigned> io_snapshot_wlcInfoVecs_currentNodeEdgeCountVector = sys->wlcInfoVecs.currentNodeEdgeCountVector;
	thrust::host_vector<unsigned> io_snapshot_wlcInfoVecs_numOriginalNeighborsNodeVector = sys->wlcInfoVecs.numOriginalNeighborsNodeVector;

	if (sys) {

		//first create a new file using the current network strain
		
		std::string format = ".sta";
		
		std::string strain =  std::to_string(sys->generalParams.currentTime);
		std::string initial = str_params+"/Param_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());



		//unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		//Now first place strain
		ofs << std::setprecision(5) <<std::fixed<< "time " << sys->generalParams.currentTime<<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "minX " << sys->domainParams.minX<<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "maxX " << sys->domainParams.maxX<<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "minY " << sys->domainParams.minY<<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "maxY " << sys->domainParams.maxY<<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "minZ " << sys->domainParams.minX<<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "maxZ " << sys->domainParams.maxX<<'\n';

		ofs << std::setprecision(5) <<std::fixed<< "total_applied_force_fibers " << sys->generalParams.totalAppliedForce_Fiber <<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "total_applied_force_plts " << sys->generalParams.totalAppliedForce_Plt <<'\n';

		ofs << std::setprecision(5) <<std::fixed<< "original_node_count " << originalNodeCount <<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "node_count_discretize " << maxNodeCount <<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "original_edge_count " << originalEdgeCount <<'\n';
		ofs << std::setprecision(5) <<std::fixed<< "edge_count_discretize " << edgeCountDiscretize <<'\n';
		
		//place nodes
		
		thrust::host_vector<double> hostNodeLocX = sys->nodeInfoVecs.nodeLocX;
		thrust::host_vector<double> hostNodeLocY = sys->nodeInfoVecs.nodeLocY;
		thrust::host_vector<double> hostNodeLocZ = sys->nodeInfoVecs.nodeLocZ;

		thrust::host_vector<double> hostPltLocX = sys->pltInfoVecs.pltLocX;
		thrust::host_vector<double> hostPltLocY = sys->pltInfoVecs.pltLocY;
		thrust::host_vector<double> hostPltLocZ = sys->pltInfoVecs.pltLocZ;

		//thrust::copy(sys->nodeInfoVecs.nodeLocX.begin(),hostEdgeRight
		//	sys->nodeInfoVecs.nodeLocX.end(), hostLocX.begin());
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			
			double x = hostNodeLocX[i];
			double y = hostNodeLocY[i];
			double z = hostNodeLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "node " << x << " " << y << " " << z <<'\n';
		}
		
		//place plts
		for (unsigned i = 0; i < sys->pltInfoVecs.pltLocX.size(); i++) {
			double x = hostPltLocX[i];
			double y = hostPltLocY[i];
			double z = hostPltLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "plt " << x << " " << y << " " << z <<'\n';
		
		}

		thrust::host_vector<double> hostPltForceX = sys->pltInfoVecs.pltForceX;
		thrust::host_vector<double> hostPltForceY = sys->pltInfoVecs.pltForceY;
		thrust::host_vector<double> hostPltForceZ = sys->pltInfoVecs.pltForceZ;
				
		//place plts force
		for (unsigned i = 0; i < sys->pltInfoVecs.pltForceX.size(); i++) {
			double x = hostPltForceX[i];
			double y = hostPltForceY[i];
			double z = hostPltForceZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "force_on_plt " <<
			x << " " << y << " " << z <<'\n';
		}
		
		thrust::host_vector<double> hostTotalExerted = sys->pltInfoVecs.totalExerted;

		//place force exerted
		for (unsigned i = 0; i < sys->pltInfoVecs.totalExerted.size(); i++) {
			double exerted = hostTotalExerted[i];
			ofs << std::setprecision(5) << std::fixed << "force_plt_exerted " << 
			exerted << '\n';
		}

		//now print plt and corresponding node it's connected to		
		unsigned maxPltCount = sys->pltInfoVecs.pltLocX.size();
		unsigned num_connections = sys->pltInfoVecs.numConnections;
		for (unsigned edge = 0; edge < num_connections; edge++ ){

			//because nodes are placed after platelets, their vtk file id is incremented. 
			//notice that this represents the vtk id, not the id within the c++ program
			unsigned node_id_vtk = maxPltCount + edge;
			unsigned node_id = io_snapshot_pltInfoVecs_nodeImagingConnection[edge];
			
			unsigned plt_id = io_snapshot_pltInfoVecs_pltImagingConnection[edge];
			
			ofs << std::setprecision(5) << std::fixed<< 
			"plt_id_node_id "<< plt_id << " "<< node_id_vtk <<'\n';

		}


		//place force node is experiencing (vector form)
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			double f_x = io_snapshot_nodeInfoVecs_nodeForceX[i];
			double f_y = io_snapshot_nodeInfoVecs_nodeForceY[i];
			double f_z = io_snapshot_nodeInfoVecs_nodeForceZ[i];
			
			ofs << std::setprecision(5) <<std::fixed<< "force_on_node " << 
			f_x <<" " << f_y <<" " << f_z << '\n';
		
		}

		//place original edges
		for (unsigned edge = 0; edge < sys->generalParams.originEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[edge];
			ofs <<"original_edge_discretized " <<idL <<" "<< idR <<'\n';
			
		}
				 
		//place added edges
		for (unsigned edge = sys->generalParams.originEdgeCount; edge < sys->generalParams.currentEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[edge];
			ofs <<"added_edge " <<idL <<" "<< idR <<'\n';
			
		}

		//original edge strain
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = io_snapshot_nodeInfoVecs_discretizedEdgeStrain[i];

			ofs << std::setprecision(5)<< std::fixed<<"original_edge_strain " << val <<'\n';
		}
				
		//original edge alignment
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = io_snapshot_nodeInfoVecs_discretizedEdgeAlignment[i];
			ofs << std::setprecision(5)<< std::fixed<<"original_edge_alignment " << val <<'\n';
		}

		//added edge strain
		for (unsigned i = sys->generalParams.originEdgeCount; i < sys->generalParams.currentEdgeCount; i++ ){
			double val = io_snapshot_nodeInfoVecs_discretizedEdgeStrain[i];
			ofs << std::setprecision(5)<< std::fixed<<"added_edge_strain " << val <<'\n';
		}
		
		//added links per node.
		for (unsigned i = 0; i < sys->generalParams.maxNodeCount; i++ ){
			unsigned val = io_snapshot_wlcInfoVecs_currentNodeEdgeCountVector[i] - 
				io_snapshot_wlcInfoVecs_numOriginalNeighborsNodeVector[i];
			ofs << std::setprecision(5)<< std::fixed<<"bind_sites_per_node " << val <<'\n';
		}



	}
};


void Storage::print_VTK_File() {

	std::shared_ptr<System> sys = system.lock();
	if (!sys) return;
	// Read-only export snapshots: avoid a separate GPU transfer for every value.
	thrust::host_vector<unsigned> io_snapshot_pltInfoVecs_nodeImagingConnection = sys->pltInfoVecs.nodeImagingConnection;
	thrust::host_vector<unsigned> io_snapshot_wlcInfoVecs_globalNeighbors = sys->wlcInfoVecs.globalNeighbors;
	thrust::host_vector<double> io_snapshot_wlcInfoVecs_lengthZero = sys->wlcInfoVecs.lengthZero;

	if (sys) {
		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation + "/FibrinNetwork_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;
		thrust::host_vector<double> hostNodeLocX = sys->nodeInfoVecs.nodeLocX;
		thrust::host_vector<double> hostNodeLocY = sys->nodeInfoVecs.nodeLocY;
		thrust::host_vector<double> hostNodeLocZ = sys->nodeInfoVecs.nodeLocZ;

		ofs.open(Filename.c_str());


		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned maxNeighborCount = (sys->generalParams).maxNeighborCount;

		unsigned numEdges = sys->generalParams.currentEdgeCount;

		ofs << "# vtk DataFile Version 3.0" << '\n';
		ofs << "Point representing Sub_cellular elem model" << '\n';
		ofs << "ASCII" << '\n' << '\n';
		ofs << "DATASET UNSTRUCTURED_GRID" << '\n';


		ofs << "POINTS " << maxNodeCount << " float" << '\n';
		for (unsigned i = 0; i< maxNodeCount; i++) { 
			double xPos = hostNodeLocX[i];
			double yPos = hostNodeLocY[i];
			double zPos = hostNodeLocZ[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles


		unsigned numCells = numEdges;
		unsigned numNumsInCells = 3 * numEdges;


		ofs << "CELLS " << numCells << " " << numNumsInCells << '\n';

		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[edge];

			ofs << 2 << " " << idL << " " << idR << '\n';
		}
	/*	for (unsigned idA = 0; idA < maxNodeCount; idA++ ){

			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = io_snapshot_wlcInfoVecs_globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX

				//counts only half
				if ((idA < idB) && (idB < maxNodeCount) ) {
					ofs<< 2 << " " << idA << " " << idB << '\n';
				}
			}
		}*/

		ofs << "CELL_TYPES " << numCells << '\n';
		for (unsigned i = 0; i< numEdges; i++) {
			ofs << 3 << '\n'; //edge joining two points
		}



		//
		ofs << "CELL_DATA " << numCells << '\n';
		ofs << "SCALARS Fiber_Strain double " << '\n';
		ofs << "LOOKUP_TABLE default "  << '\n';
		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idA = sys->nodeInfoVecs.hostEdgeLeft[edge];
			unsigned idB = sys->nodeInfoVecs.hostEdgeRight[edge];

			unsigned begin = idA * sys->generalParams.maxNeighborCount;
			unsigned end = begin + sys->generalParams.maxNeighborCount;
			double L0;
			for (unsigned i = begin; i < end; i++) {
				unsigned idTemp = io_snapshot_wlcInfoVecs_globalNeighbors[i];
				if (idTemp == idB){
					L0 = io_snapshot_wlcInfoVecs_lengthZero[i];
				}
			}
			double xL = hostNodeLocX[idA];
			double yL = hostNodeLocY[idA];
			double zL = hostNodeLocZ[idA];
			double xR = hostNodeLocX[idB];
			double yR = hostNodeLocY[idB];
			double zR = hostNodeLocZ[idB];

			double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
			double strain = (L1 - L0) / L0;
			ofs << std::fixed << strain   << '\n';

		}
		ofs.close();

	}

	//now print platelets 
	if ((sys)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation+"/Platelet";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned maxPltCount = sys->generalParams.maxPltCount;
		
		double xPos;
		double yPos;
		double zPos;

		thrust::host_vector<double> hostPltLocX = sys->pltInfoVecs.pltLocX;
		thrust::host_vector<double> hostPltLocY = sys->pltInfoVecs.pltLocY;
		thrust::host_vector<double> hostPltLocZ = sys->pltInfoVecs.pltLocZ;
		
		ofs << "# vtk DataFile Version 3.0" << '\n';
		ofs << "Point representing Sub_cellular elem model" << '\n';
		ofs << "ASCII" << '\n' << '\n';
		ofs << "DATASET UNSTRUCTURED_GRID" << '\n';
		
		 
		ofs << "POINTS " << maxPltCount  << " float" << '\n';
		for (unsigned i = 0; i < maxPltCount; i++) {
			xPos = hostPltLocX[i];
			yPos = hostPltLocY[i];
			zPos = hostPltLocZ[i];

			ofs << std::setprecision(5) << std::fixed << xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}

		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;

		unsigned numNumsInCells = 1 + maxPltCount;

		ofs << "CELLS " << numCells << " " << numNumsInCells << '\n';
		
		//place edges as cells of type 2. 
		ofs << maxPltCount;
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs << " " << point;
		}
		ofs <<" "<< '\n';

		ofs << "CELL_TYPES " << numCells << '\n';  
		//set edges and last set scattered points
				
		ofs << 2 << '\n';//scatter points for capsid

		/*thrust::host_vector<double> temp1 = sys->pltInfoVecs.totalExerted;
		thrust::host_vector<double> temp2 = sys->pltInfoVecs.filopodiaSaturation;
		thrust::host_vector<double> temp3 = sys->pltInfoVecs.numFinalInteractions;
		for (unsigned j = 0; j < temp1.size(); j++) { 
			std::cout << "Force: " << temp1[j] << " on platelet " << j << std::endl;
			std::cout << "Filo_sat: " << temp2[j] << " on platelet " << j << std::endl;
			std::cout << "Count: " << temp3[j] << " on platelet " << j << std::endl;
		}*/
		
		thrust::host_vector<double> hostTotalExerted = sys->pltInfoVecs.totalExerted;
		ofs << '\n' << "POINT_DATA " << maxPltCount << '\n';
		ofs << "SCALARS Platelet_Force double" << '\n';
		ofs << "LOOKUP_TABLE default "  << '\n';
		double plt_force_vtk;
		for (unsigned point = 0; point < maxPltCount; point++) {
			plt_force_vtk = hostTotalExerted[point];
			ofs << std::fixed << plt_force_vtk << '\n';
		}

		thrust::host_vector<double> hostFilopodiaSaturation = sys->pltInfoVecs.filopodiaSaturation;
		ofs << '\n' << "SCALARS Filopodia_Saturation double" << '\n';
		ofs << "LOOKUP_TABLE default "  << '\n';
		double plt_filo_sat;
		for (unsigned point = 0; point < maxPltCount; point++) {
			plt_filo_sat = hostFilopodiaSaturation[point];
			ofs << std::fixed << plt_filo_sat << '\n';
		}
		
		thrust::host_vector<double> hostNumFinalInteractions = sys->pltInfoVecs.numFinalInteractions;
		ofs << '\n' << "SCALARS Filopodia_Saturation double" << '\n';
		ofs << "LOOKUP_TABLE default "  << '\n';
		double num_final_interactions;
		for (unsigned point = 0; point < maxPltCount; point++) {
			num_final_interactions = hostNumFinalInteractions[point];
			ofs << std::fixed << num_final_interactions << '\n';
		}
		ofs.close();
	}
	
	//now print platelets with attatchments
	if ((sys)) {
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation+"/PlateletConn";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());
		
	
		unsigned maxPltCount = sys->pltInfoVecs.pltLocX.size();

		unsigned num_connections = sys->pltInfoVecs.numConnections;

		thrust::host_vector<double> hostPltLocX = sys->pltInfoVecs.pltLocX;
		thrust::host_vector<double> hostPltLocY = sys->pltInfoVecs.pltLocY;
		thrust::host_vector<double> hostPltLocZ = sys->pltInfoVecs.pltLocZ;
		
		double xPos;
		double yPos;
		double zPos;
		
		ofs << "# vtk DataFile Version 3.0" << '\n';
		ofs << "Point representing Sub_cellular elem model" << '\n';
		ofs << "ASCII" << '\n' << '\n';
		ofs << "DATASET UNSTRUCTURED_GRID" << '\n';
		
		
		
		ofs << "POINTS " << maxPltCount + num_connections << " float" << '\n';
		for (unsigned i = 0; i < maxPltCount; i++) {
			xPos = hostPltLocX[i];
			yPos = hostPltLocY[i];
			zPos = hostPltLocZ[i];

			ofs << std::setprecision(5) << std::fixed << xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}

		//set location for nodes that plt is connected to
		//ie  

		thrust::host_vector<double> hostNodeLocX = sys->nodeInfoVecs.nodeLocX;
		thrust::host_vector<double> hostNodeLocY = sys->nodeInfoVecs.nodeLocY;
		thrust::host_vector<double> hostNodeLocZ = sys->nodeInfoVecs.nodeLocZ;

		for (unsigned i = 0; i < num_connections; i++ ) {
			unsigned node_id = io_snapshot_pltInfoVecs_nodeImagingConnection[i];
			if (node_id < sys->generalParams.maxNodeCount) {
				xPos = hostNodeLocX[node_id];
				yPos = hostNodeLocY[node_id];
				zPos = hostNodeLocZ[node_id];
				
				ofs << std::setprecision(5) << std::fixed << xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
			} else {
				std::cout<<"imaging not working, node out of bounds" << std::endl;
			}
		}


		//std::cout<<'here1'<<std::flush;
		
		unsigned numCells = 1;
		numCells += num_connections;//add conections cells for edges

		unsigned numNumsInCells = 1 + maxPltCount;
		numNumsInCells += 3 * num_connections;//3 numbers per edge

		ofs << "CELLS " << numCells << " " << numNumsInCells << '\n';
		
		//place edges as cells of type 2. 
		ofs<< maxPltCount;
		for (unsigned point = 0; point < maxPltCount; point++ ){
			ofs<< " " << point;
		}
		ofs<<" "<< '\n';

		thrust::host_vector<unsigned> hostNodeImagingConnection = sys->pltInfoVecs.nodeImagingConnection;
		thrust::host_vector<unsigned> hostPltImagingConnection = sys->pltInfoVecs.pltImagingConnection;

		for (unsigned edge = 0; edge < num_connections; edge++) {

			//because nodes are placed after platelets, their vtk file id is incremented. 
			//notice that this represents the vtk id, not the id within the c++ program
			unsigned node_id_vtk = maxPltCount + edge;
			unsigned node_id = hostNodeImagingConnection[edge];
			
			unsigned plt_id = hostPltImagingConnection[edge];
				
			ofs << 2 << " " << node_id_vtk << " " << plt_id << '\n';

		}
		ofs << "CELL_TYPES " << numCells << '\n';  
		//set edges and last set scattered points
				
		ofs << 2 << '\n';//scatter points for capsid
		
	//	std::cout<<'here3'<<std::flush;
		for (unsigned edge = 0; edge< num_connections; edge++ ){
			ofs<< 3 <<'\n';
		}
		ofs.close();
	}

};
