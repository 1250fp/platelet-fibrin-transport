
#ifndef PLT_ARM_NODE_FORCE_TIME_H_
#define PLT_ARM_NODE_FORCE_TIME_H_


struct NodeInfoVecs;
struct WLCInfoVecs;
struct GeneralParams;
struct PltInfoVecs;
struct AuxVecs;
struct RandVecs;
    //filopodia-like force
void Plt_Arm_Node_Force_Time(
  	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
  	PltInfoVecs& pltInfoVecs,
  	AuxVecs& auxVecs,
	RandVecs& randVecs);

#endif