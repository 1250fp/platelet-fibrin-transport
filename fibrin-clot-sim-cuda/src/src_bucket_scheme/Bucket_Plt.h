#ifndef BUCKET_PLT_H_
#define BUCKET_PLT_H_

struct NodeInfoVecs;
struct PltInfoVecs;
struct GeneralParams;
struct DomainParams;
struct WLCInfoVecs;
struct TorsionInfoVecs;
struct AuxVecs;


void init_plt_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);


void build_plt_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);


void extend_plt_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

#endif