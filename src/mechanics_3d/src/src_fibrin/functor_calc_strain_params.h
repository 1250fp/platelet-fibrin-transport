#ifndef FUNCTOR_CALC_STRAIN_PARAMS_H_
#define FUNCTOR_CALC_STRAIN_PARAMS_H_

#include "SystemStructures.h"
#include <climits>
struct functor_calc_strain_params : public thrust::unary_function<Tuu, Tdd> {
	unsigned originLinkCount;
	unsigned originEdgeCount;
	unsigned originNodeCount;
	unsigned maxNodeCount;
	unsigned maxNeighborCount;

	double* locXAddr;
	double* locYAddr;
	double* locZAddr;

	unsigned* originNbrVec;
	unsigned* currentNbrVec;
	unsigned* globalNeighbors;
	double* lenZero;

	__host__ __device__
	functor_calc_strain_params(
		unsigned& _originLinkCount,
		unsigned& _originEdgeCount,
		unsigned& _originNodeCount,
		unsigned& _maxNodeCount,
		unsigned& _maxNeighborCount,	

		double* _locXAddr, 
		double* _locYAddr, 
		double* _locZAddr,

		unsigned* _originNbrVec,
		unsigned* _currentNbrVec,
		unsigned* _globalNeighbors,
		double* _lenZero) :
		originLinkCount(_originLinkCount),
		originEdgeCount(_originEdgeCount),
		originNodeCount(_originNodeCount),
		maxNodeCount(_maxNodeCount),
		maxNeighborCount(_maxNeighborCount),

		locXAddr(_locXAddr), 
		locYAddr(_locYAddr), 
		locZAddr(_locZAddr),
		
		originNbrVec(_originNbrVec),
		currentNbrVec(_currentNbrVec),
		globalNeighbors(_globalNeighbors),
		lenZero(_lenZero){}
 
	__device__
	Tdd operator() (const Tuu& u2) {
		unsigned idL = thrust::get<0>(u2);//use IDL as row
		unsigned idR = thrust::get<1>(u2);//use IDR as col
		double invalid = nan("");
		double edgeStrain = 0.0;
		double edgeAlignment = 0.0;

		// The table allocation uses unsigned maxNodeCount * maxNeighborCount.
		if (idL >= maxNodeCount || idR >= maxNodeCount || maxNeighborCount == 0 ||
			maxNodeCount > UINT_MAX / maxNeighborCount) {
			return thrust::make_tuple(invalid, invalid);
		}
		unsigned idMat = UINT_MAX;
		for (unsigned slot = 0; slot < maxNeighborCount; slot++) {
			unsigned index = idL * maxNeighborCount + slot;
			if (idR == globalNeighbors[index]) {
				idMat = index;
				break;
			}
		}
		if (idMat == UINT_MAX) {
			return thrust::make_tuple(invalid, invalid);
		}
		double lengthZero = lenZero[idMat];

		if (lengthZero > 0.0) {
			double posXA_XB = locXAddr[idL] - locXAddr[idR];
			double posYA_YB = locYAddr[idL] - locYAddr[idR];
			double posZA_ZB = locZAddr[idL] - locZAddr[idR];
			double currentLength = sqrt(
				(posXA_XB) * (posXA_XB) +
				(posYA_YB) * (posYA_YB) +
				(posZA_ZB) * (posZA_ZB));
			edgeStrain = ((currentLength - lengthZero) / lengthZero);
			edgeAlignment = currentLength > 0.0 ? abs(1.0 * posZA_ZB / currentLength) : invalid;
		}
		return thrust::make_tuple(edgeStrain, edgeAlignment);
	}
};

#endif