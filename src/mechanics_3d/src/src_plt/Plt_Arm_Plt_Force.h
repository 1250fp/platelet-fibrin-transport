#ifndef PLT_ARM_PLT_FORCE_H_
#define PLT_ARM_PLT_FORCE_H_


struct GeneralParams;
struct PltInfoVecs;
struct AuxVecs;

void Plt_Arm_Plt_Force(
	GeneralParams& generalParams,
	PltInfoVecs& pltInfoVecs,
	AuxVecs& auxVecs);


#endif