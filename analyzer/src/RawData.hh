/*
  RawData.hh

  2012/1/24
*/

#ifndef RawData_h
#define RawData_h

#include "DetectorID.hh"
#include <vector>

#ifdef MemoryLeak
#include "DebugCounter.hh"
#endif

class HodoRawHit;
class DCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector <DCRawHit *> DCRHitContainer;

class RawData
{

private:
  HodoRHitContainer GCRawHC;
  HodoRHitContainer BH1RawHC;
  HodoRHitContainer BH2RawHC;
  HodoRHitContainer BACRawHC;
  HodoRHitContainer BFTRawHC[NumOfPlaneBFT];
  HodoRHitContainer SFTRawHC[NumOfPlaneSFT];

  HodoRHitContainer CFTRawHC[NumOfPlaneCFT];
  HodoRHitContainer BGORawHC;
  HodoRHitContainer T54CounterRawHC[NumOfPlaneT54Counter];
  HodoRHitContainer PiVRawHC;

  HodoRHitContainer TOFRawHC;
  HodoRHitContainer CHRawHC;
  HodoRHitContainer ACRawHC;
  HodoRHitContainer T0RawHC;
  HodoRHitContainer AC1RawHC;
  HodoRHitContainer LCRawHC;

  DCRHitContainer BcInRawHC[NumOfLayersBcIn+1],   BcOutRawHC[NumOfLayersBcOut+1];
  DCRHitContainer SdcInRawHC[NumOfLayersSdcIn+1], SdcOutRawHC[NumOfLayersSdcOut+1];
  DCRHitContainer SsdRawHC[NumOfLayersSsd+1];

  HodoRHitContainer ScalerRawHC;
  HodoRHitContainer MiscRawHC;
  HodoRHitContainer MatrixRawHC;

#ifdef MemoryLeak
  static debug::Counter sm_counter;
#endif

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeHits();

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  bool AddHodoRawHit(HodoRHitContainer& cont,
 		     int DetId,
 		     int Plane,
 		     int Seg,
 		     int AorT,
 		     int UorD,
 		     int Data);


  bool AddDCRawHit( DCRHitContainer &cont, 
		    int Plane, 
		    int Wire, 
		    int Tdc,
		    int type=0); 

public:
  bool AddCFTHodoRawHit(int DetId,
			int Plane,
			int Seg,
			int AorT,
			int UorD,
			double Data );

  bool AddBGOHodoRawHit(int DetId,
			int Plane,
			int Seg,
			int AorT,
			int UorD,
			double Data );

  bool AddTOFHodoRawHit(int DetId,
			int Plane,
			int Seg,
			int AorT,
			int UorD,
			double Data );

  bool AddLCHodoRawHit(int DetId,
			int Plane,
			int Seg,
			int AorT,
			int UorD,
			double Data );

  bool AddCHHodoRawHit(int DetId,
		       int Plane,
		       int Seg,
		       int AorT,
		       int UorD,
		       double Data );

  bool AddACHodoRawHit(int DetId,
		       int Plane,
		       int Seg,
		       int AorT,
		       int UorD,
		       double Data );

  bool AddT0HodoRawHit(int DetId,
		       int Plane,
		       int Seg,
		       int AorT,
		       int UorD,
		       double Data );

  bool AddPiVHodoRawHit(int DetId,
			int Plane,
			int Seg,
			int AorT,
			int UorD,
			double Data );
  
  bool AddSdcInRawHit( int Plane,
		       int Wire,
		       double dTime);

  bool AddSdcOutRawHit( int Plane,
		       int Wire,
		       double dTime);

  bool AddBcOutRawHit( int Plane,
		       int Wire,
		       double dTime);

  const HodoRHitContainer& GetGCRawHC() const;
  const HodoRHitContainer& GetBH1RawHC() const;
  const HodoRHitContainer& GetBH2RawHC() const;
  const HodoRHitContainer& GetBACRawHC() const;
  const HodoRHitContainer & GetBFTRawHC( int plane ) const;
  const HodoRHitContainer & GetSFTRawHC( int plane ) const;

  const HodoRHitContainer& GetCHRawHC() const;
  const HodoRHitContainer& GetTOFRawHC() const;
  const HodoRHitContainer& GetACRawHC() const;
  const HodoRHitContainer& GetT0RawHC() const;
  const HodoRHitContainer& GetAC1RawHC() const;
  const HodoRHitContainer& GetLCRawHC() const;

  const DCRHitContainer & GetBcInRawHC( int layer ) const;
  const DCRHitContainer & GetBcOutRawHC( int layer ) const;
  const DCRHitContainer & GetSsdRawHC( int layer ) const;
  const DCRHitContainer & GetSdcInRawHC( int layer ) const;
  const DCRHitContainer & GetSdcOutRawHC( int layer ) const;

  const HodoRHitContainer& GetScalerRawHC() const;
  const HodoRHitContainer& GetMiscRawHC() const;
  const HodoRHitContainer& GetMatrixRawHC() const;

  const HodoRHitContainer & GetCFTRawHC( int plane ) const;
  const HodoRHitContainer & GetBGORawHC() const;
  const HodoRHitContainer & GetPiVRawHC() const;
  const HodoRHitContainer & GetT54CounterRawHC( int plane ) const;
};

#endif
