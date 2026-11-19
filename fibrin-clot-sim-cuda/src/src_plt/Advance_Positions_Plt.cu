#include "SystemStructures.h"
#include "functor_advance_pos.h"
#include "System.h"
#include "Advance_Positions_Plt.h"


void Advance_Positions_Plt(
	PltInfoVecs& pltInfoVecs,
	GeneralParams& generalParams,
	RandVecs& randVecs) {

	unsigned _seedplt = rand();

	thrust::counting_iterator<unsigned> index_sequence_plt_begin(_seedplt);


	
 	thrust::transform(thrust::device, index_sequence_plt_begin, index_sequence_plt_begin + (generalParams.maxPltCount),
		randVecs.gaussianPltData.begin(), psrunifgen(-1.0, 1.0));

thrust::counting_iterator<unsigned> pltIndexBegin(0);

thrust::transform(
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 pltIndexBegin,
		 pltInfoVecs.pltLocX.begin(),
		 pltInfoVecs.pltLocY.begin(),
		 pltInfoVecs.pltLocZ.begin())),
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 pltIndexBegin,
		 pltInfoVecs.pltLocX.begin(),
		 pltInfoVecs.pltLocY.begin(),
		 pltInfoVecs.pltLocZ.begin())) + generalParams.maxPltCount,
 //second vector begin
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 randVecs.gaussianPltData.begin(),
		 pltInfoVecs.pltForceX.begin(),
		 pltInfoVecs.pltForceY.begin(),
		 pltInfoVecs.pltForceZ.begin())),
 //save result in third vector to test values
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 pltInfoVecs.pltLocX.begin(),
		 pltInfoVecs.pltLocY.begin(),
		 pltInfoVecs.pltLocZ.begin(),
		 pltInfoVecs.pltVelocity.begin())),
 functor_advance_pos(generalParams.dtTemp,
	 generalParams.viscousDamp_Plt,
	 generalParams.temperature,
	 generalParams.kB,
	 generalParams.pltMass,
	 generalParams.maxPltCount,
	 thrust::raw_pointer_cast(pltInfoVecs.isPltFixed.data())));


};
