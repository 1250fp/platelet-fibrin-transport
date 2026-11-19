#ifndef FUNCTOR_PLT_ARM_NODE_TIME_H_
#define FUNCTOR_PLT_ARM_NODE_TIME_H_

#include "SystemStructures.h"

#include <curand.h>
#include <curand_kernel.h>

struct functor_plt_arm_node_time : public thrust::unary_function< U3DAD2D, double>  {
	double current_time;
	double CLM;
	double wlc_factor; // NKBT/P
	bool distribute_plt_force;

  	unsigned plt_tndrl_intrct;
  	double pltRForce;
  	double pltForce;
  	double pltR;

  	unsigned maxPltCount;
  	double fiberDiameter;
  	unsigned maxNodeCount;
  	unsigned maxIdCountFlag;
  	unsigned maxNeighborCount;
	bool pltrelease;
	bool plthandhand;
	double strain_switch;

  	double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
	double* nodeUForceXAddr;
	double* nodeUForceYAddr;
	double* nodeUForceZAddr;

  	unsigned* nodeUId;
  	unsigned* pltUId;

	unsigned* id_value_expanded;
	unsigned* keyBegin;
	unsigned* keyEnd;

  	unsigned* tndrlNodeId;
  	unsigned* tndrlNodeType;
	bool* isNodeInPltVolAddr;
  	unsigned* glblNghbrsId;
	double* lengthZero;
	unsigned* numOriginalNeighborsNodeVector;

  	double* pltLocXAddr;
	double* pltLocYAddr;
	double* pltLocZAddr;


   __host__ __device__
   //
       functor_plt_arm_node_time(
			double& _current_time,
			double& _CLM,
			double& _wlc_factor,
			bool& _distribute_plt_force,

            unsigned& _plt_tndrl_intrct,
            double& _pltRForce,
            double& _pltForce,
            double& _pltR,

            unsigned& _maxPltCount,
            double& _fiberDiameter,
            unsigned& _maxNodeCount,
            unsigned& _maxIdCountFlag,
            unsigned& _maxNeighborCount,
			bool& _pltrelease,
			bool& _plthandhand,
			double& _strain_switch,

            double* _nodeLocXAddr,
            double* _nodeLocYAddr,
            double* _nodeLocZAddr,
            double* _nodeUForceXAddr,
            double* _nodeUForceYAddr,
            double* _nodeUForceZAddr,

      		unsigned* _nodeUId,
            unsigned* _pltUId,

      		unsigned* _id_value_expanded,
      		unsigned* _keyBegin,
      		unsigned* _keyEnd,

            unsigned* _tndrlNodeId,
            unsigned* _tndrlNodeType,
			bool* _isNodeInPltVolAddr,
            unsigned* _glblNghbrsId,

			double* _lengthZero,
			unsigned* _numOriginalNeighborsNodeVector,

            double* _pltLocXAddr,
            double* _pltLocYAddr,
            double* _pltLocZAddr) :

	current_time(_current_time),
	CLM(_CLM),
	wlc_factor(_wlc_factor),
	distribute_plt_force(_distribute_plt_force),

    plt_tndrl_intrct(_plt_tndrl_intrct),
    pltRForce(_pltRForce),
    pltForce(_pltForce),
    pltR(_pltR),

    maxPltCount(_maxPltCount),
    fiberDiameter(_fiberDiameter),
    maxNodeCount(_maxNodeCount),
    maxIdCountFlag(_maxIdCountFlag),
    maxNeighborCount(_maxNeighborCount),
	pltrelease(_pltrelease),
	plthandhand(_plthandhand),
	strain_switch(_strain_switch),

    nodeLocXAddr(_nodeLocXAddr),
	nodeLocYAddr(_nodeLocYAddr),
	nodeLocZAddr(_nodeLocZAddr),
	nodeUForceXAddr(_nodeUForceXAddr),
	nodeUForceYAddr(_nodeUForceYAddr),
	nodeUForceZAddr(_nodeUForceZAddr),

    nodeUId(_nodeUId),
    pltUId(_pltUId),

	id_value_expanded(_id_value_expanded),
	keyBegin(_keyBegin),
	keyEnd(_keyEnd),

    tndrlNodeId(_tndrlNodeId),
    tndrlNodeType(_tndrlNodeType),
	isNodeInPltVolAddr(_isNodeInPltVolAddr),
    glblNghbrsId(_glblNghbrsId),

	lengthZero(_lengthZero),
	numOriginalNeighborsNodeVector(_numOriginalNeighborsNodeVector),

    pltLocXAddr(_pltLocXAddr),
	pltLocYAddr(_pltLocYAddr),
	pltLocZAddr(_pltLocZAddr){}

   __device__
 		CVec3UCVec2 operator()(const U3CVec7 &u3d6, const CVec2 & cv2) {
		
		unsigned pltId = thrust::get<0>(u3d6);//platelet id should be the counting iterator
		unsigned num_plt_filo = thrust::get<1>(u3d6);		
        unsigned bucketId = thrust::get<2>(u3d6);

		double unif_rand = thrust::get<3>(u3d6);//used for generator unif

        unsigned storageLocation = pltId * plt_tndrl_intrct;

        double pltLocX = thrust::get<4>(u3d6);
        double pltLocY = thrust::get<5>(u3d6);
        double pltLocZ = thrust::get<6>(u3d6);

        //use for return.
        double pltCurrentForceX = thrust::get<7>(u3d6);
        double pltCurrentForceY = thrust::get<8>(u3d6);
        double pltCurrentForceZ = thrust::get<9>(u3d6);

        double sumPltForceX = pltCurrentForceX;
        double sumPltForceY = pltCurrentForceY;
        double sumPltForceZ = pltCurrentForceZ;
		double sumTotalExerted = thrust::get<0>(cv2);
		double saturation = thrust::get<1>(cv2);

		//set random generator

        //beginning and end of attempted interaction network nodes.
		unsigned beginIndex = keyBegin[bucketId];
		unsigned endIndex = keyEnd[bucketId];
		double len = static_cast<double>(endIndex-beginIndex);

		double seed_d = round(len*unif_rand);
		unsigned seed = static_cast<unsigned>(seed_d);
		curandState state;
		curand_init ( seed, pltId, 0, &state);

        //Loop through the number of available filopodia
		unsigned final_interaction_count = 0;
        for(unsigned interactionCounter = 0; interactionCounter < num_plt_filo; interactionCounter++) {

            unsigned pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
            unsigned pullNode_type = tndrlNodeType[storageLocation + interactionCounter];

			//CHECK IF PLT PULLS NODE
			//WARNING: RESETS NBRs
            if (( pullNode_id != maxIdCountFlag) &&
               ( pullNode_type == 0)) {
			    //then we have a node connected to the plt.
                //TYPE 0 is network

                //Calculate distance from plt to node.
                //unsigned pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
                //Get position of node
                double vecN_PX = pltLocX - nodeLocXAddr[pullNode_id];
                double vecN_PY = pltLocY - nodeLocYAddr[pullNode_id];
                double vecN_PZ = pltLocZ - nodeLocZAddr[pullNode_id];

                double dist = sqrt(
                    (vecN_PX) * (vecN_PX) +
                    (vecN_PY) * (vecN_PY) +
                    (vecN_PZ) * (vecN_PZ));

                //check if the node is not in pulling range anymore.
                if ((dist >= pltRForce) || (dist <= (pltR + fiberDiameter / 2.0 ) ) ) {

                  	//then node is out of range, so we empty tendril if pltrelease is true

					if (pltrelease == true){
						tndrlNodeId[storageLocation + interactionCounter] = maxIdCountFlag;//reset
					}
                  	//try to find a new node to pull within connections of previous node if plthandhand is true

					//HOW TO MAKE THE PLT CHOOSE ANOTHER
					//POINT WHEN TWO PLT's MEET?
					if (plthandhand == true) {
						unsigned startIndex_1 = maxNeighborCount * pullNode_id;
						unsigned endIndex_1 = startIndex_1 + maxNeighborCount;

						for (unsigned nbr_loc = startIndex_1; nbr_loc < endIndex_1; nbr_loc++) {

								unsigned newPullNode_id = glblNghbrsId[ nbr_loc ];
								//check tentative node is not already connected
								for (unsigned checkId = 0; checkId < num_plt_filo; checkId++) {
									if (newPullNode_id != tndrlNodeId[storageLocation + checkId]) {

										bool isNodeInPltVol = false;
										if (pltrelease) {
											isNodeInPltVol = isNodeInPltVolAddr[newPullNode_id];
										}
										//then newPullNode_id isn't yet pulled.
										//We can break out of this for statement
										//and check if it is close enough

										//This ensures that pltrelease controls the node in plt volume
										//variables.
										if (isNodeInPltVol == false){
											break;
										}
									}
								}
							//check neighbor not empty
							if (newPullNode_id < maxNodeCount){//maxNodeCount is default neighbor value.
								 vecN_PX = pltLocX - nodeLocXAddr[newPullNode_id];
								 vecN_PY = pltLocY - nodeLocYAddr[newPullNode_id];
								 vecN_PZ = pltLocZ - nodeLocZAddr[newPullNode_id];
								//Calculate distance from plt to node.
								 dist = sqrt(
									(vecN_PX) * (vecN_PX) +
									(vecN_PY) * (vecN_PY) +
									(vecN_PZ) * (vecN_PZ));

								//check if new node is in interaction range and fill tendril with new node than break neighbors loop
								if ((dist < pltRForce) && (dist > (pltR + fiberDiameter / 2.0) ) ) {//pull this node
									tndrlNodeId[storageLocation + interactionCounter] = newPullNode_id;//bucketNbrsExp[i];
									tndrlNodeType[storageLocation + interactionCounter] = 0;//assign type
									break;
								}
							}
						}
					}
					//end hand-hand
                }
            }

        	//CHECK IF PLT PULLS PLT INSTEAD OF NETWORK
        	else if ( (pullNode_id != maxIdCountFlag) &&
        	    ( pullNode_type == 1) ) {//note this happens only if plt-plt interaction is on
        	  	//then we have a plt connected to the plt.
				//TYPE 1 is plt

        	    //Calc range
        	    unsigned pullPlt_id = tndrlNodeId[storageLocation + interactionCounter];//bucketNbrsExp[i];
        	    //Get position of plt
        	    double vecN_PX = pltLocX - pltLocXAddr[pullPlt_id];
        	    double vecN_PY = pltLocY - pltLocYAddr[pullPlt_id];
        	    double vecN_PZ = pltLocZ - pltLocZAddr[pullPlt_id];
        	    double dist = sqrt(
        	        (vecN_PX) * (vecN_PX) +
        	        (vecN_PY) * (vecN_PY) +
        	        (vecN_PZ) * (vecN_PZ));

        	    //check if the plt is not pulled  anymore and pltrelease is true
        	    if (((dist >= 2.0 * pltRForce) || (dist <= 2.0 * pltR)) && pltrelease==true ){
        	        //then plt is out of range so we disconnect


					tndrlNodeId[storageLocation + interactionCounter] = maxIdCountFlag;
        	    }
        	}

			//after check, re generate choice of node and type.
            pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
            pullNode_type = tndrlNodeType[storageLocation + interactionCounter];

        	//CHECK IF PLT HAS NOTHING
			//THEN SEARCH FOR NODE
			if (pullNode_id == maxIdCountFlag) {
				//then we have nothing to pull.
        	    //try to find a node to pull by searching.

				//ISSUE HERE: we need a random permutation of nodes.
				//
				double bucket_len = static_cast<double>(endIndex-beginIndex);
        	    for (unsigned newpull_index = beginIndex; newpull_index < endIndex; newpull_index++){

					float randomf = curand_uniform( &state );//uniform set above.
					unsigned rand_choice = static_cast<unsigned>(floor(randomf * bucket_len)) + beginIndex;
        	        unsigned newPullNode_id = id_value_expanded[ rand_choice];//could be newpull_index

					bool isNodeInPltVol = false;
					if (pltrelease) {
						isNodeInPltVol = isNodeInPltVolAddr[newPullNode_id];
					}

					bool node_is_new = true;
					//check tentative node is not already connected
				    for (unsigned checkId = 0; checkId < num_plt_filo; checkId++){
        	        	if (newPullNode_id == tndrlNodeId[storageLocation + checkId]){
							node_is_new = false;
        	        		break;
        	        	}
        	      	}

					//only pull on new nodes that(if pltrelease is on) are not in plt volume
					if ( (node_is_new) && (isNodeInPltVol == false) ){

						double vecN_PX = pltLocX - nodeLocXAddr[newPullNode_id];
						double vecN_PY = pltLocY - nodeLocYAddr[newPullNode_id];
						double vecN_PZ = pltLocZ - nodeLocZAddr[newPullNode_id];
						//Calculate distance from plt to node.
						double dist = sqrt(
						(vecN_PX) * (vecN_PX) +
						(vecN_PY) * (vecN_PY) +
						(vecN_PZ) * (vecN_PZ));

						//check if new node is in interaction range and fill tenril with new unique node
						//then break searching loop
						if ((dist < pltRForce) &&
							(dist > (pltR + fiberDiameter / 2.0) ) ) {
							//pull this node

							tndrlNodeId[storageLocation + interactionCounter] = newPullNode_id;//bucketNbrsExp[i];
							tndrlNodeType[storageLocation + interactionCounter] = 0;//assign type
							break;
						}
					}
        	    }

        	}


        	//WARNING: ONLY PULLING NODES
			//UP TO HERE WE HAVE CHECKED FOR NODES
        	//check if tendril has been filled with a node and apply pulling forces. Note if filled direction and distence of forces are already calculated

			//after last check, re generate choice of node and type.
            pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
            pullNode_type = tndrlNodeType[storageLocation + interactionCounter];

			if ( (pullNode_id != maxIdCountFlag)
        	   && ( pullNode_type == 0) ) {
        	    //then we have a post-search node we can pull.
				//Add force to it and the current platelet.

				pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
				double vecN_PX = pltLocX - nodeLocXAddr[pullNode_id];
        	    double vecN_PY = pltLocY - nodeLocYAddr[pullNode_id];
        	    double vecN_PZ = pltLocZ - nodeLocZAddr[pullNode_id];
        	    //Calculate distance from plt to node.
        	    double dist = sqrt(
        	       (vecN_PX) * (vecN_PX) +
        	       (vecN_PY) * (vecN_PY) +
        	       (vecN_PZ) * (vecN_PZ));
				    
                    if ((dist < pltRForce) && (dist > (pltR + fiberDiameter / 2.0))){
                        double forceNodeX;
                        double forceNodeY;
                        double forceNodeZ;
                        double mag_force;

                        double strain_count = 0;
                        double sum_strain=0;
                        //first calculate the strain of neighbors of pullNode_id
                        unsigned begin_index = pullNode_id * maxNeighborCount;
                        unsigned num_original_connections = numOriginalNeighborsNodeVector[pullNode_id];

                        //can change added value to maxNeighborCount
                        unsigned end_index = begin_index + num_original_connections;
                        for (unsigned pt_index = begin_index; pt_index < end_index; pt_index++) {
                            unsigned pt = glblNghbrsId[pt_index];
                            double dist_0 = lengthZero[pt_index];
                            if ((pt < maxNodeCount) && (pt != pullNode_id)) {
                                //then we have a neighbor and we can calculate the strain
                                double vecN_PX = nodeLocXAddr[pullNode_id] - nodeLocXAddr[pt];
                                double vecN_PY = nodeLocYAddr[pullNode_id] - nodeLocYAddr[pt];
                                double vecN_PZ = nodeLocZAddr[pullNode_id] - nodeLocZAddr[pt];
                                //Calculate distance from plt to node.
                                double dist = sqrt(
                                    (vecN_PX) * (vecN_PX)+
                                    (vecN_PY) * (vecN_PY)+
                                    (vecN_PZ) * (vecN_PZ));
                                sum_strain += fabsf((dist - dist_0) / dist_0);
                                strain_count+=1.0;
                            }
                        }
                        double ave_strain=0;
                        if (strain_count > 0) {
                            ave_strain = sum_strain / (strain_count);
                        }
                        double term1 = 1.0 / CLM;
                        double term2 = CLM*CLM / (2.0*(CLM - ave_strain));
                        double stiffness = wlc_factor * (term1  + term2);

                        //double tau = 636.6;
                        double a = 61.0;
                        double b = 3.582;
                        double c = 131.8;
                        //double d = 1.417;
                        double f_0 = 2.233;//
						double force_reduction_factor = 1.0; // 
                        double fitting_factor = (a * stiffness) / (b * stiffness + c);
                        //mag_force = f_0 + (1.0 - d*exp(-(current_time / tau))) * fitting_factor;
                        mag_force = force_reduction_factor * (f_0 + fitting_factor) / 2.0; // SI eqn 13

                        //at high strain for certain fitting parameters. 
                        if (mag_force < 0.0) {
                            mag_force = 0.0;
                        }

						sumTotalExerted += mag_force;
                        //Determine direction of force based on positions and multiply magnitude force
                        forceNodeX = (vecN_PX / dist) * (mag_force);
                        forceNodeY = (vecN_PY / dist) * (mag_force);
                        forceNodeZ = (vecN_PZ / dist) * (mag_force);

                        //count force for self plt.
                        sumPltForceX += (-1.0) * forceNodeX;
                        sumPltForceY += (-1.0) * forceNodeY;
                        sumPltForceZ += (-1.0) * forceNodeZ;

                        //store force in temporary vector if a node is pulled. Call reduction later.
                        //Note: this is the only force storage, we need a different increment (final_interaction_count)
                        //so that we can rescale only those forces stored. Interaction counter might not apply force, so it could be different.
                        nodeUForceXAddr[storageLocation + final_interaction_count] = forceNodeX;
                        nodeUForceYAddr[storageLocation + final_interaction_count] = forceNodeY;
                        nodeUForceZAddr[storageLocation + final_interaction_count] = forceNodeZ;
                        nodeUId[storageLocation + final_interaction_count] = pullNode_id;
                        pltUId[storageLocation + final_interaction_count] = pltId;

                        final_interaction_count += 1;

				}

        	}

        }

	//The last step is to scale the forces we stored if distribute_plt_force is true.
	double divisor = __ll2double_ru(final_interaction_count);
	if (distribute_plt_force == true) {
		if (divisor > 0.0) {
			sumPltForceX = sumPltForceX/divisor;
			sumPltForceY = sumPltForceY/divisor;
			sumPltForceZ = sumPltForceZ/divisor;
		}
		for (unsigned arm = 0; arm < final_interaction_count; arm++) {
			if (divisor > 0.0) {
				double forceNodeX = nodeUForceXAddr[storageLocation + arm];
				double forceNodeY = nodeUForceYAddr[storageLocation + arm];
				double forceNodeZ = nodeUForceZAddr[storageLocation + arm];
				//store force in temporary vector if a node is pulled. Call reduction later.
				nodeUForceXAddr[storageLocation + arm] = forceNodeX/divisor;
				nodeUForceYAddr[storageLocation + arm] = forceNodeY/divisor;
				nodeUForceZAddr[storageLocation + arm] = forceNodeZ/divisor;
			}
		}
	}
    //return platelet forces
	if (__ll2double_ru(num_plt_filo) > 0.0) { 
		saturation = divisor / __ll2double_ru(num_plt_filo) ;
	} else {
		saturation = -1;
	}
    return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ, final_interaction_count, sumTotalExerted, saturation);

   }
};

#endif
