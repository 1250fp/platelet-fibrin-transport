
#include "Storage.h"
#include "Link_Nodes.h"
#include "WLC_Force.h"
#include "Torsion_Force.h"
#include "Plt_Arm_Node_Force_Time.h"
#include "Plt_Arm_Plt_Force.h"
#include "Plt_Field_Node_Force.h"
#include "Plt_Field_Plt_Force.h"
#include "Plt_Vol_Exc_Force.h"

#include "Params_Calc.h"
#include "Advance_Positions_Fibrin.h"
#include "Advance_Positions_Plt.h"

#include "Bucket_Net.h"
#include "Bucket_Plt.h"
#include "functor_misc.h"
#include "System.h"
#include <random>
#include <cstdlib>
#include <cerrno>
#include <stdexcept>


void System::setBucketScheme() {

	init_dim_general(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	init_net_inct_bucket(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	build_net_inct_bucket(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	extend_net_inct_bucket(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);
		
	init_plt_inct_bucket(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	build_plt_inct_bucket(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	extend_plt_inct_bucket(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);
		
};

void System::solveForces() {
 
	//RESET FORCE TO ZERO AT BEGINNING/////////////////////////////////////////////////
	thrust::fill(nodeInfoVecs.nodeForceX.begin(),nodeInfoVecs.nodeForceX.end(),0);
	thrust::fill(nodeInfoVecs.nodeForceY.begin(),nodeInfoVecs.nodeForceY.end(),0); 
	thrust::fill(nodeInfoVecs.nodeForceZ.begin(),nodeInfoVecs.nodeForceZ.end(),0);
	
	thrust::fill(pltInfoVecs.pltForceX.begin(),pltInfoVecs.pltForceX.end(),0);
	thrust::fill(pltInfoVecs.pltForceY.begin(),pltInfoVecs.pltForceY.end(),0);
	thrust::fill(pltInfoVecs.pltForceZ.begin(),pltInfoVecs.pltForceZ.end(),0);

	
	if (generalParams.linking == true) {
		
		Link_Nodes(
			nodeInfoVecs,
			wlcInfoVecs,
			auxVecs,
			generalParams);
	}
	Torsion_Force(nodeInfoVecs, torsionInfoVecs, generalParams);

	//std::cout<<"prewlc"<<std::endl;
	WLC_Force(nodeInfoVecs, wlcInfoVecs, generalParams);

	//platetelet-node forces
	//RESETS PLATELET FORCES
	if (generalParams.pltfrcfld == true) {// note: this force-field includes both pulling and pushing
		Plt_Field_Node_Force(//platelet on node force field
			nodeInfoVecs,
			wlcInfoVecs,
			generalParams,
			pltInfoVecs,
			auxVecs);
		if (generalParams.pltonplt == true) {
			Plt_Field_Plt_Force(//platelet on platelet interaction through force field
				generalParams,
				pltInfoVecs,
				auxVecs);
		}

	}
	else if (generalParams.plttndrl == true) { //note for now force-field type has priority over tndrl-type

		// filopodia-node pulling
		//using time dependence function
		Plt_Arm_Node_Force_Time(
		  nodeInfoVecs,
		  wlcInfoVecs,
		  generalParams,
		  pltInfoVecs,
		  auxVecs,
		  randVecs);

		//Tndrl-Plt pulling
		if (generalParams.pltonplt == true) {
			/*Plt_Arm_Plt_Force(//platelet on platelet interaction through tndrl
				generalParams,
				pltInfoVecs,
				auxVecs);*/
		}

		Plt_Vol_Exc_Force(//push for volume exclusion
			nodeInfoVecs,
			wlcInfoVecs,
			generalParams,
			pltInfoVecs,
			auxVecs);

	}




};


void System::solveSystem() {

	//make sure iterationCounter is zero for padding init. 
	generalParams.iterationCounter = 0;
	//set initial bucket scheme
	setBucketScheme();

	//set initial epsilon
	generalParams.epsilon = (1.0) *
		sqrt(6.0*generalParams.kB * generalParams.temperature * generalParams.dtTemp / generalParams.viscousDamp_Fibrin);

	double final_time = 45.0 * 60.0;
    // Optional pilot limit. Without the environment variable the original limit is preserved.
    if (const char* text = std::getenv("MICHAEL_FINAL_TIME_SECONDS")) {
        char* end = nullptr;
        errno = 0;
        const double requested = std::strtod(text, &end);
        if (errno != 0 || end == text || *end != '\0' || !std::isfinite(requested) || requested <= 0.0) {
            throw std::runtime_error("MICHAEL_FINAL_TIME_SECONDS must be finite and positive");
        }
        final_time = requested;
        std::cout << "Pilot simulated-time limit: " << final_time << " seconds" << std::endl;
    }//convert seconds to minutes requires *60

	//Fill force with zero to double check
	thrust::fill(nodeInfoVecs.nodeForceX.begin(),nodeInfoVecs.nodeForceX.end(),0);
	thrust::fill(nodeInfoVecs.nodeForceY.begin(),nodeInfoVecs.nodeForceY.end(),0); 
	thrust::fill(nodeInfoVecs.nodeForceZ.begin(),nodeInfoVecs.nodeForceZ.end(),0);
	
	thrust::fill(pltInfoVecs.pltForceX.begin(),pltInfoVecs.pltForceX.end(),0);
	thrust::fill(pltInfoVecs.pltForceY.begin(),pltInfoVecs.pltForceY.end(),0);
	thrust::fill(pltInfoVecs.pltForceZ.begin(),pltInfoVecs.pltForceZ.end(),0);

	thrust::fill(pltInfoVecs.totalExerted.begin(),pltInfoVecs.totalExerted.end(),0);
	thrust::fill(pltInfoVecs.filopodiaSaturation.begin(),pltInfoVecs.filopodiaSaturation.end(),0);
	thrust::fill(pltInfoVecs.numFinalInteractions.begin(),pltInfoVecs.numFinalInteractions.end(),0);

	while (generalParams.runSim == true) {

		double time_iter = (generalParams.iterationCounter);

		//Simulations force quite after 45 minutes of real time running 
		if ((time_iter * generalParams.dtTemp) > final_time ) {
			generalParams.runSim = false;
		};
		
		generalParams.iterationCounter++;
		/*std::cout<<"iter: " << generalParams.iterationCounter << std::endl;

		thrust::device_vector<double>::iterator iter_x =
  			thrust::max_element(nodeInfoVecs.nodeForceX.begin(), nodeInfoVecs.nodeForceX.end());
		thrust::device_vector<double>::iterator iter_y =
  			thrust::max_element(nodeInfoVecs.nodeForceY.begin(), nodeInfoVecs.nodeForceY.end());
		thrust::device_vector<double>::iterator iter_z =
  			thrust::max_element(nodeInfoVecs.nodeForceZ.begin(), nodeInfoVecs.nodeForceZ.end());

		double max_val_x = *iter_x;
		double max_val_y = *iter_y;
		double max_val_z = *iter_z;

		std::cout << "Max Forces: " << max_val_x << " " << max_val_y << " "<< max_val_z<< std::endl;
		*/
		
		generalParams.totalAppliedForce_Fiber = thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.nodeForceX.begin(),
					nodeInfoVecs.nodeForceY.begin(),
					nodeInfoVecs.nodeForceZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.nodeForceX.end(),
					nodeInfoVecs.nodeForceY.end(),
					nodeInfoVecs.nodeForceZ.end())),
				functor_norm(), 0.0, thrust::plus<double>() );
	
		generalParams.totalAppliedForce_Plt = thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					pltInfoVecs.pltForceX.begin(),
					pltInfoVecs.pltForceY.begin(),
					pltInfoVecs.pltForceZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					pltInfoVecs.pltForceX.end(),
					pltInfoVecs.pltForceY.end(),
					pltInfoVecs.pltForceZ.end())),
				functor_norm(), 0.0, thrust::plus<double>() );
		//std::cout<<"total applied force: " << generalParams.totalAppliedForce_Plt << std::endl;
		
		generalParams.currentTime += generalParams.dtTemp;

		Advance_Positions_Fibrin(
			nodeInfoVecs,
			generalParams,
			randVecs);

		Advance_Positions_Plt(
			pltInfoVecs,
			generalParams,
			randVecs);
		
		if (generalParams.iterationCounter % 100 == 0) {
			setBucketScheme();
			//std::cout<<"time: " << generalParams.currentTime << std::endl;
		}

		
		solveForces(); //resets and solves forces for next time step
		double temp_current_iter = static_cast<double>(generalParams.iterationCounter);
		double temp_current_time = temp_current_iter * generalParams.dtTemp;
		if ( ( fmod(temp_current_time, 60.0) == 0.0) || (generalParams.iterationCounter == 10) ) {

			storage->save_current_state();
			storage->print_VTK_File();
			//store sum of all forces on each node. Used in stress calculations
			//store before upadting storage class.

			//WARNING BEFORE CALLING SAVE_PARAMS CALCULATE THEM FIRST
			Params_Calc(
    			wlcInfoVecs,
    			nodeInfoVecs,
    			generalParams,
    			pltInfoVecs);

			storage->save_params();

			generalParams.epsilon = (1.0) *
				sqrt(6.0 * generalParams.kB * generalParams.temperature * generalParams.dtTemp / generalParams.viscousDamp_Fibrin);

		}
	}

};

System::System()  {};

void System::assignStorage(std::shared_ptr<Storage> _storage) {
	storage = _storage;
}

//__host__ __device__
void System::initializeSystem(
	thrust::host_vector<bool>& hostIsNodeFixed,
	thrust::host_vector<double>& hostPosX,
	thrust::host_vector<double>& hostPosY,
	thrust::host_vector<double>& hostPosZ,
	thrust::host_vector<unsigned>& hostWLCEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCEdgeRight,
	thrust::host_vector<double>& hostWLCLenZero,

	thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
	thrust::host_vector<double>& hostWLCSubLenZero,
	thrust::host_vector<unsigned>& hostTorsionIndexLeft,
	thrust::host_vector<unsigned>& hostTorsionIndexCenter,
	thrust::host_vector<unsigned>& hostTorsionIndexRight,
	thrust::host_vector<double>& hostTorsionAngleZero,
	//platelets
	thrust::host_vector<bool>& hostIsPltFixed,
	thrust::host_vector<double>& hostPltPosX,
	thrust::host_vector<double>& hostPltPosY,
	thrust::host_vector<double>& hostPltPosZ) {

	std::cout<< "total Edge Count: "<< generalParams.originEdgeCount << std::endl;
	std::cout << "max num nodes: " << generalParams.maxNodeCount << std::endl;
	//platelets

	std::cout << "max num platelets in device: " << generalParams.maxPltCount << std::endl;



	setPltVecs(
		hostIsPltFixed,
		hostPltPosX,
		hostPltPosY,
		hostPltPosZ);
	std::cout << "done setting plt" << std::endl;

	setNodeVecs(//calls initDimensionBucketScheme
		hostIsNodeFixed,
		hostPosX,
		hostPosY,
		hostPosZ);

	std::cout << "done setting node" << std::endl;

	setTorsionVecs(
		hostTorsionIndexLeft,
		hostTorsionIndexCenter,
		hostTorsionIndexRight,
		hostTorsionAngleZero);
	
	std::cout << "done setting torsion" << std::endl;

	setWLCVecs(
		hostWLCEdgeLeft,
		hostWLCEdgeRight,
		hostWLCLenZero );

	std::cout << "done setting wlc" << std::endl;

	setRandVecs();
		
};

void System::setRandVecs(){
	randVecs.filopodia_detach.resize(pltInfoVecs.tndrlNodeId.size());
	randVecs.gaussianData.resize(generalParams.maxNodeCount);
	randVecs.gaussianPltData.resize(generalParams.maxPltCount); //
	randVecs.bucketPltStart.resize(generalParams.maxPltCount);
		
}

void System::setNodeVecs(
	thrust::host_vector<bool>& hostIsNodeFixed,
	thrust::host_vector<double>& hostPosX,
	thrust::host_vector<double>& hostPosY,
	thrust::host_vector<double>& hostPosZ) {


	nodeInfoVecs.sumForcesOnNode.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeVelocity.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeLocX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLocY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLocZ.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeForceX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeForceY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeForceZ.resize(generalParams.maxNodeCount);

	nodeInfoVecs.discretizedEdgeStrain.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.discretizedEdgeAlignment.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);

	//sized larger for input later
	nodeInfoVecs.hostEdgeLeft.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.hostEdgeRight.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.deviceEdgeLeft.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.deviceEdgeRight.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);


	thrust::fill(nodeInfoVecs.discretizedEdgeStrain.begin(), nodeInfoVecs.discretizedEdgeStrain.end(),0.0);
	thrust::fill(nodeInfoVecs.hostEdgeRight.begin(), nodeInfoVecs.hostEdgeRight.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.hostEdgeLeft.begin(), nodeInfoVecs.hostEdgeLeft.end(), 0);	//fill force and velocity with zeros for computation.

	thrust::fill(nodeInfoVecs.sumForcesOnNode.begin(), nodeInfoVecs.sumForcesOnNode.end(), 0);


	thrust::copy(hostPosX.begin(), hostPosX.end(), nodeInfoVecs.nodeLocX.begin());
	thrust::copy(hostPosY.begin(), hostPosY.end(), nodeInfoVecs.nodeLocY.begin());
	thrust::copy(hostPosZ.begin(), hostPosZ.end(), nodeInfoVecs.nodeLocZ.begin());


	//copy fixed positions
	nodeInfoVecs.isNodeFixed.resize(generalParams.maxNodeCount);
	thrust::copy(hostIsNodeFixed.begin(), hostIsNodeFixed.end(), nodeInfoVecs.isNodeFixed.begin());

	nodeInfoVecs.isNodeInPltVol.resize(generalParams.maxNodeCount);
	thrust::fill(nodeInfoVecs.isNodeInPltVol.begin(),nodeInfoVecs.isNodeInPltVol.end(),false);

	nodeInfoVecs.linksThreadMade.resize(generalParams.maxNodeCount);
	nodeInfoVecs.delinksThreadMade.resize(generalParams.maxNodeCount);
	nodeInfoVecs.idMadeTempLeft.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);
	nodeInfoVecs.idMadeTempRight.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);

	std::cout << "setting host and dim" <<std::endl;
	nodeInfoVecs.host_id_left.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);
	nodeInfoVecs.host_id_right.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);


	//at this point all nodes are filled, so we can generate domainParams
	init_dim_general(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);


	domainParams.originMinX = domainParams.minX;
	domainParams.originMaxX = domainParams.maxX;
	domainParams.originMinY = domainParams.minY;
	domainParams.originMaxY = domainParams.maxY;
	domainParams.originMinZ = domainParams.minZ;
	domainParams.originMaxZ = domainParams.maxZ;

	std::cout<< "node count : " <<nodeInfoVecs.nodeLocY.size()<< std::endl;


	auxVecs.id_bucket_net_intc.resize(generalParams.maxNodeCount);
	auxVecs.id_value_net_intc.resize(generalParams.maxNodeCount);
	auxVecs.id_bucket_expanded_net_intc.resize(27 * (generalParams.maxNodeCount));
	auxVecs.id_value_expanded_net_intc.resize(27 *( generalParams.maxNodeCount ));

	
	auxVecs.id_bucket_plt_intc.resize(generalParams.maxNodeCount);
	auxVecs.id_value_plt_intc.resize(generalParams.maxNodeCount);
	auxVecs.id_bucket_expanded_plt_intc.resize(27 * (generalParams.maxNodeCount));
	auxVecs.id_value_expanded_plt_intc.resize(27 *( generalParams.maxNodeCount ));

	
	auxVecs.idPlt_bucket.resize(generalParams.maxPltCount);
	auxVecs.idPlt_value.resize(generalParams.maxPltCount);
	auxVecs.idPlt_bucket_expanded.resize(27 * (generalParams.maxPltCount));
	auxVecs.idPlt_value_expanded.resize(27 *( generalParams.maxPltCount ));

};

//platelet
void System::setPltVecs(
	thrust::host_vector<bool>& hostIsPltFixed,
	thrust::host_vector<double>& hostPltPosX,
	thrust::host_vector<double>& hostPltPosY,
	thrust::host_vector<double>& hostPltPosZ) {


	pltInfoVecs.sumForcesOnPlt.resize(generalParams.maxPltCount);

	pltInfoVecs.pltVelocity.resize(generalParams.maxPltCount);

	pltInfoVecs.pltLocX.resize(generalParams.maxPltCount);
	pltInfoVecs.pltLocY.resize(generalParams.maxPltCount);
	pltInfoVecs.pltLocZ.resize(generalParams.maxPltCount);

	//pltInfoVecs.pltLocs.resize(generalParams.maxPltCount);
	//pltInfoVecs.pltLocsAddresses.resize(generalParams.maxPltCount);
	//for (unsigned i = 0; i < pltInfoVecs.pltLocs.size(); i++) {
	//	pltInfoVecs.pltLocsAddresses[i] = pltInfoVecs.pltLocs.data();
	//}

	pltInfoVecs.pltForceX.resize(generalParams.maxPltCount);
	pltInfoVecs.pltForceY.resize(generalParams.maxPltCount);
	pltInfoVecs.pltForceZ.resize(generalParams.maxPltCount);
	//pltInfoVecs.pltForce.resize(generalParams.maxPltCount);
	//pltInfoVecs.pltForceAddresses.resize(generalParams.maxPltCount);
	//for (unsigned i = 0; i < pltInfoVecs.pltForceAddresses.size(); i++) {
	//	pltInfoVecs.pltForceAddresses[i] = pltInfoVecs.pltForce.data();
	//}
	
	pltInfoVecs.totalExerted.resize(generalParams.maxPltCount);
	pltInfoVecs.filopodiaSaturation.resize(generalParams.maxPltCount);
	pltInfoVecs.numFinalInteractions.resize(generalParams.maxPltCount);

	pltInfoVecs.pltImagingConnection.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeImagingConnection.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);

	pltInfoVecs.nodeUnreducedId.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeUnreducedForceX.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeUnreducedForceY.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeUnreducedForceZ.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);

	pltInfoVecs.nodeReducedId.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeReducedForceX.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeReducedForceY.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);
	pltInfoVecs.nodeReducedForceZ.resize(generalParams.maxPltCount * generalParams.plt_other_intrct);

	thrust::fill(pltInfoVecs.sumForcesOnPlt.begin(), pltInfoVecs.sumForcesOnPlt.end(), 0);

	thrust::copy(hostPltPosX.begin(), hostPltPosX.end(), pltInfoVecs.pltLocX.begin());
	thrust::copy(hostPltPosY.begin(), hostPltPosY.end(), pltInfoVecs.pltLocY.begin());
	thrust::copy(hostPltPosZ.begin(), hostPltPosZ.end(), pltInfoVecs.pltLocZ.begin());


	std::cout<<"num platelets: "<< pltInfoVecs.pltLocX.size() << std::endl;
	std::cout<<"num platelets var: "<< generalParams.maxPltCount << std::endl;
	//copy fixed positions
	pltInfoVecs.isPltFixed.resize(generalParams.maxPltCount);
	thrust::fill(pltInfoVecs.isPltFixed.begin(), pltInfoVecs.isPltFixed.end(), false);
	//thrust::copy(hostIsPltFixed.begin(), hostIsPltFixed.end(), pltInfoVecs.isPltFixed.begin());


	auxVecs.idPlt_bucket.resize(generalParams.maxPltCount);
	auxVecs.idPlt_value.resize(generalParams.maxPltCount);
	auxVecs.idPlt_bucket_expanded.resize(27 *( generalParams.maxPltCount ));
	auxVecs.idPlt_value_expanded.resize(27 * (generalParams.maxPltCount));

    pltInfoVecs.tndrlNodeId.resize(generalParams.maxPltCount * generalParams.plt_tndrl_intrct);
	pltInfoVecs.tndrlNodeType.resize(generalParams.maxPltCount * generalParams.plt_tndrl_intrct);
	pltInfoVecs.tndrlNodeNum.resize(generalParams.maxPltCount);

	//fill with flag vales.	
	std::cout<<"maxIdFlag: "<< generalParams.maxIdCountFlag<<std::endl;

	thrust::fill(pltInfoVecs.tndrlNodeId.begin(),pltInfoVecs.tndrlNodeId.end(), generalParams.maxIdCountFlag);
	
	thrust::fill(pltInfoVecs.tndrlNodeType.begin(),pltInfoVecs.tndrlNodeType.end(), 0);
	if (generalParams.tndrl_plt_pdf){
		std::default_random_engine generator;
		std::lognormal_distribution<double> distribution(generalParams.tndrl_M,generalParams.tndrl_S);

		for (int i=0; i<generalParams.maxPltCount; ++i) {
			double number = -1;

			while ((number<1.0)||(number>generalParams.plt_tndrl_intrct+1)) {
				number = distribution(generator);
				std::cout<<number<<std::endl;
			}
			std::cout<<"plt number filo"<<std::endl;

			pltInfoVecs.tndrlNodeNum[i]=unsigned(number);
			std::cout<<pltInfoVecs.tndrlNodeNum[i]<<std::endl;
		}
	}
	else{
		thrust::fill(pltInfoVecs.tndrlNodeNum.begin(),pltInfoVecs.tndrlNodeNum.end(), generalParams.plt_tndrl_intrct);
	}
  
};

void System::setTorsionVecs(
	thrust::host_vector<unsigned>& hostTorsionIndexLeft,
	thrust::host_vector<unsigned>& hostTorsionIndexCenter,
	thrust::host_vector<unsigned>& hostTorsionIndexRight,
	thrust::host_vector<double>& hostTorsionAngleZero) {


	torsionInfoVecs.leftIndex.resize(generalParams.totalTorsionCount);
	torsionInfoVecs.centerIndex.resize(generalParams.totalTorsionCount);
	torsionInfoVecs.rightIndex.resize(generalParams.totalTorsionCount);
	torsionInfoVecs.angleZero.resize(generalParams.totalTorsionCount);

	thrust::copy(hostTorsionIndexLeft.begin(), hostTorsionIndexLeft.end(), torsionInfoVecs.leftIndex.begin());
	thrust::copy(hostTorsionIndexCenter.begin(), hostTorsionIndexCenter.end(), torsionInfoVecs.centerIndex.begin());
	thrust::copy(hostTorsionIndexRight.begin(), hostTorsionIndexRight.end(), torsionInfoVecs.rightIndex.begin());

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				torsionInfoVecs.leftIndex.begin(),
				torsionInfoVecs.centerIndex.begin(),
				torsionInfoVecs.rightIndex.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				torsionInfoVecs.leftIndex.begin(),
				torsionInfoVecs.centerIndex.begin(),
				torsionInfoVecs.rightIndex.begin())) + generalParams.totalTorsionCount,
			torsionInfoVecs.angleZero.begin(),//save vector
		TorsionAngleFunctor(
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data())));

	torsionInfoVecs.forceX.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.forceY.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.forceZ.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceX.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceY.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceZ.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);

	thrust::fill(torsionInfoVecs.forceX.begin(), torsionInfoVecs.forceX.end(), 0.0);
	thrust::fill(torsionInfoVecs.forceY.begin(), torsionInfoVecs.forceY.end(), 0.0);
	thrust::fill(torsionInfoVecs.forceZ.begin(), torsionInfoVecs.forceZ.end(), 0.0);

	torsionInfoVecs.tempTorIndices.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.reducedIds.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
};

void System::setWLCVecs(
	thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
	thrust::host_vector<double>& hostWLCSubLenZero ) {

	wlcInfoVecs.globalNeighbors.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	wlcInfoVecs.currentNodeEdgeCountVector.resize(generalParams.maxNodeCount);

	wlcInfoVecs.lengthZero.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	wlcInfoVecs.numOriginalNeighborsNodeVector.resize(generalParams.maxNodeCount);

	//default value is maxNodeCount
	thrust::fill(wlcInfoVecs.globalNeighbors.begin(), wlcInfoVecs.globalNeighbors.end(), generalParams.maxNodeCount);
	thrust::fill(wlcInfoVecs.currentNodeEdgeCountVector.begin(), wlcInfoVecs.currentNodeEdgeCountVector.end(),0);
	thrust::fill(wlcInfoVecs.lengthZero.begin(), wlcInfoVecs.lengthZero.end(), 0.0);

	nodeInfoVecs.hostEdgeLeft = hostWLCSubEdgeLeft;
	nodeInfoVecs.hostEdgeRight = hostWLCSubEdgeRight;

	//scan through hostAdj and put in device.
	for (unsigned id = 0; id < hostWLCSubLenZero.size(); id++) {

		unsigned idL = hostWLCSubEdgeLeft[id];
		unsigned idR = hostWLCSubEdgeRight[id];

		double edgeLen = hostWLCSubLenZero[id];
		//we use the lengthZero vector to identify edges as well.
		//node id is row, column node is connected to row node.

		//add edge for left node
		unsigned edgeNumL = wlcInfoVecs.currentNodeEdgeCountVector[idL]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexL = idL*generalParams.maxNeighborCount + edgeNumL;
		wlcInfoVecs.lengthZero[indexL] = edgeLen;
		wlcInfoVecs.globalNeighbors[indexL] = idR;
		(wlcInfoVecs.currentNodeEdgeCountVector[idL])++; //right connects to left

		//add edge for right node
		unsigned edgeNumR = wlcInfoVecs.currentNodeEdgeCountVector[idR]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexR = idR*generalParams.maxNeighborCount + edgeNumR;
		wlcInfoVecs.lengthZero[indexR] = edgeLen;
		wlcInfoVecs.globalNeighbors[indexR] = idL;
		(wlcInfoVecs.currentNodeEdgeCountVector[idR])++; //left connects to right

		generalParams.currentEdgeCount++;

	}

	//at this point currentNodeEdgeCountVector holds the number of edges, copy this to
	thrust::copy(wlcInfoVecs.currentNodeEdgeCountVector.begin(), wlcInfoVecs.currentNodeEdgeCountVector.end(), wlcInfoVecs.numOriginalNeighborsNodeVector.begin());
};


