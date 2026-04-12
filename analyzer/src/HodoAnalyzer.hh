/*
  HodoAnalyzer.hh

  2012/1/24
*/

#ifndef HodoAnalyzer_h
#define HodoAnalyzer_h

#include <vector>
#include <cstddef>

#include "DetectorID.hh"
#include "RawData.hh"


class RawData;
class Hodo2Hit;
class Hodo1Hit;
class BH2Hit;
class FiberHit;
class CFTFiberHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;
class CFTFiberCluster;

typedef std::vector <Hodo2Hit*> Hodo2HitContainer;
typedef std::vector <Hodo1Hit*> Hodo1HitContainer;
typedef std::vector <BH2Hit*>   BH2HitContainer;
typedef std::vector <FiberHit*> FiberHitContainer;
typedef std::vector <CFTFiberHit*> CFTFiberHitContainer;

typedef std::vector <HodoCluster*>  HodoClusterContainer;
typedef std::vector <BH2Cluster*>   BH2ClusterContainer;
typedef std::vector <FiberCluster*> FiberClusterContainer;
typedef std::vector <CFTFiberCluster*> CFTFiberClusterContainer;

class HodoAnalyzer
{
public:
  HodoAnalyzer();
  ~HodoAnalyzer();

  static HodoAnalyzer& GetInstance();
private:
  HodoAnalyzer(const HodoAnalyzer &);
  HodoAnalyzer & operator =(const HodoAnalyzer &);

private:
  Hodo1HitContainer GCCont;
  Hodo2HitContainer BH1Cont;
  BH2HitContainer   BH2Cont;
  FiberHitContainer BFTCont;
  FiberHitContainer SFTCont[NumOfLayersSFT];
  CFTFiberHitContainer CFTCont[NumOfLayersCFT];
  Hodo2HitContainer    BGOCont;
  Hodo2HitContainer    PiVCont;
  Hodo1HitContainer BACCont;
  Hodo1HitContainer TgtCont;
  Hodo2HitContainer CHCont;
  Hodo2HitContainer TOFCont;
  Hodo1HitContainer ACCont;
  Hodo1HitContainer T0Cont;
  Hodo2HitContainer LCCont;

  HodoClusterContainer  BH1ClCont;
  BH2ClusterContainer   BH2ClCont;
  FiberClusterContainer BFTClCont;
  FiberClusterContainer SFTClCont[NumOfLayersSFT];
  CFTFiberClusterContainer CFTClCont[NumOfLayersCFT];
  HodoClusterContainer  BGOClCont;
  HodoClusterContainer  TOFClCont;
  HodoClusterContainer  LCClCont;

public:
  bool DecodeRawHits(RawData* rawData);

  bool DecodeGCHits(RawData* rawData);  
  bool DecodeBH1Hits(RawData* rawData);
  bool DecodeBH2Hits(RawData* rawData);
  bool DecodeBFTHits(RawData* rawData);
  bool DecodeSFTHits(RawData* rawData);
  bool DecodeCFTHits(RawData* rawData);
  bool DecodeBGOHits(RawData* rawData);
  bool DecodePiVHits(RawData* rawData);
  bool DecodeBACHits(RawData* rawData);
  bool DecodeTGTHits(RawData* rawData);
  bool DecodeCHHits(RawData* rawData);
  bool DecodeTOFHits(RawData* rawData);
  bool DecodeACHits(RawData* rawData);
  bool DecodeT0Hits(RawData* rawData);
  bool DecodeLCHits(RawData* rawData);

  int GetNHitsGC( void ) const { return GCCont.size(); };
  int GetNHitsBH1( void ) const { return BH1Cont.size(); };
  int GetNHitsBH2( void ) const { return BH2Cont.size(); };
  int GetNHitsBFT( void ) const { return BFTCont.size(); };
  int GetNHitsSFT( int layer ) const { return SFTCont[layer].size(); };
  int GetNHitsCFT( int layer ) const { return CFTCont[layer].size(); };
  int GetNHitsBGO( void ) const { return BGOCont.size(); };
  int GetNHitsPiV( void ) const { return PiVCont.size(); };
  int GetNHitsBAC( void ) const { return BACCont.size(); };
  int GetNHitsCH( void ) const { return CHCont.size(); };
  int GetNHitsTOF( void ) const { return TOFCont.size(); };
  int GetNHitsAC( void ) const { return ACCont.size(); };
  int GetNHitsT0( void ) const { return T0Cont.size(); };
  int GetNHitsLC( void ) const { return LCCont.size(); };

  inline Hodo1Hit * GetHitGC( int i ) const;
  inline Hodo2Hit * GetHitBH1( int i ) const;
  inline BH2Hit   * GetHitBH2( int i ) const;
  inline FiberHit * GetHitBFT(int seg) const;
  inline FiberHit * GetHitSFT(int lyaer, int seg) const;
  inline CFTFiberHit * GetHitCFT(int lyaer, int seg) const;
  inline Hodo2Hit * GetHitBGO( int i ) const;
  inline Hodo2Hit * GetHitPiV( int i ) const;
  inline Hodo1Hit * GetHitBAC( int i ) const;
  inline Hodo2Hit * GetHitCH( int i ) const;
  inline Hodo2Hit * GetHitTOF( int i ) const;
  inline Hodo1Hit * GetHitAC( int i ) const;
  inline Hodo1Hit * GetHitT0( int i ) const;
  inline Hodo2Hit * GetHitLC( int i ) const;

  int GetNClustersBH1( void ) const { return BH1ClCont.size(); };
  int GetNClustersBH2( void ) const { return BH2ClCont.size(); };
  int GetNClustersBFT( void ) const { return BFTClCont.size(); };
  int GetNClustersSFT( int layer ) const { return SFTClCont[layer].size(); };
  int GetNClustersCFT( int layer ) const { return CFTClCont[layer].size(); };
  int GetNClustersBGO( void ) const { return BGOClCont.size(); }
  int GetNClustersTOF( void ) const { return TOFClCont.size(); }
  int GetNClustersLC( void )  const { return LCClCont.size(); }

  inline HodoCluster  * GetClusterBH1( int i ) const;
  inline BH2Cluster   * GetClusterBH2( int i ) const;
  inline FiberCluster * GetClusterBFT( int i ) const;
  inline FiberCluster * GetClusterSFT( int layer, int i ) const;
  inline CFTFiberCluster * GetClusterCFT( int layer, int i ) const;
  inline HodoCluster  * GetClusterBGO( int i ) const;
  inline HodoCluster  * GetClusterTOF( int i ) const;
  inline HodoCluster  * GetClusterLC( int i )  const;

  bool ReCalcGCHits( bool applyRecursively=false );
  bool ReCalcBH1Hits( bool applyRecursively=false );
  bool ReCalcBH2Hits( bool applyRecursively=false );
  bool ReCalcBACHits( bool applyRecursively=false );
  bool ReCalcTOFHits( bool applyRecursively=false );
  bool ReCalcACHits( bool applyRecursively=false );
  bool ReCalcT0Hits( bool applyRecursively=false );
  bool ReCalcLCHits( bool applyRecursively=false );
  
  bool ReCalcBH1Clusters( bool applyRecursively=false );
  bool ReCalcBH2Clusters( bool applyRecursively=false );
  bool ReCalcTOFClusters( bool applyRecursively=false );
  bool ReCalcLCClusters( bool applyRecursively=false );

  bool ReCalcAll( void );
  
  void TimeCutBH1(double tmin, double tmax);
  void TimeCutBH2(double tmin, double tmax);
  void TimeCutBFT(double tmin, double tmax);
  void TimeCutSFT(int layer, double tmin, double tmax);
  void TimeCutCFT(int layer, double tmin, double tmax);

private:
  static HodoAnalyzer *MySelf_;

  void clearGCHits();
  void clearBH1Hits();
  void clearBH2Hits();
  void clearBFTHits();
  void clearSFTHits();
  void clearCFTHits();
  void clearBGOHits();
  void clearPiVHits();
  void clearBACHits();
  void clearCHHits( void );
  void clearTOFHits( void );
  void clearACHits( void );
  void clearT0Hits( void );
  void clearLCHits( void );

  bool fl_decoded_SFT;
  bool fl_decoded_CFT;

  void DecodeSFT(RawData* rawData,
		 int Begin, int End,
		 int layer,
		 char* Name);

  void DecodeCFT(RawData* rawData,
		 int Begin, int End,
		 int layer,
		 char* Name);

  template<typename TypeCluster>
  void TimeCut(std::vector<TypeCluster>& cont, double tmin, double tmax);

  static int MakeUpClusters( const Hodo2HitContainer & HitCont,
			     HodoClusterContainer & ClusterCont,
			     double maxTimeDif );

  static int MakeUpClustersBGO( const Hodo2HitContainer & HitCont,
				HodoClusterContainer & ClusterCont, 
				double maxTimeDif );

  static int MakeUpClusters( const BH2HitContainer & HitCont,
			     BH2ClusterContainer & ClusterCont,
			     double maxTimeDif );

  static int MakeUpClusters( const FiberHitContainer & cont,
			     FiberClusterContainer& ClusterCont,
			     double maxTimeDif);

  static int MakeUpClusters( const CFTFiberHitContainer & cont,
			     CFTFiberClusterContainer& ClusterCont,
			     double maxTimeDif);
};

inline HodoCluster * HodoAnalyzer::GetClusterBH1( int i ) const
{
  if( i>=0 && i<BH1ClCont.size() )
    return BH1ClCont[i];
  else
    return 0;
}

inline BH2Cluster * HodoAnalyzer::GetClusterBH2( int i ) const
{
  if( i>=0 && i<BH2ClCont.size() )
    return BH2ClCont[i];
  else
    return 0;
}

inline FiberCluster * HodoAnalyzer::GetClusterBFT( int i ) const
{
  if( i>=0 && i<BFTClCont.size() )
    return BFTClCont[i];
  else
    return 0;
}

inline FiberCluster * HodoAnalyzer::GetClusterSFT( int layer, int i ) const
{
  if( i>=0 && i<SFTClCont[layer].size() )
    return SFTClCont[layer][i];
  else
    return 0;
}

inline CFTFiberCluster * HodoAnalyzer::GetClusterCFT( int layer, int i ) const
{
  if( i>=0 && i<CFTClCont[layer].size() )
    return CFTClCont[layer][i];
  else
    return 0;
}

inline HodoCluster * HodoAnalyzer::GetClusterBGO( int i ) const
{
  if( i>=0 && i<BGOClCont.size() )
    return BGOClCont[i];
  else
    return 0;
}

inline HodoCluster * HodoAnalyzer::GetClusterTOF( int i ) const
{
  if( i>=0 && i<TOFClCont.size() )
    return TOFClCont[i];
  else
    return 0;
}

inline HodoCluster * HodoAnalyzer::GetClusterLC( int i ) const
{
  if( i>=0 && i<LCClCont.size() )
    return LCClCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitBH1( int i ) const
{
  if( i>=0 && i<BH1Cont.size() )
    return BH1Cont[i];
  else
    return 0;
}

inline BH2Hit * HodoAnalyzer::GetHitBH2( int i ) const
{
  if( i>=0 && i<BH2Cont.size() )
    return BH2Cont[i];
  else
    return 0;
}

inline FiberHit* HodoAnalyzer::GetHitBFT(int seg) const
{
  if( seg>=0 && seg<BFTCont.size() )
    return BFTCont[seg];
  else
    return NULL;
}

inline FiberHit* HodoAnalyzer::GetHitSFT(int layer, int seg) const
{
  if( seg>=0 && seg<SFTCont[layer].size() )
    return SFTCont[layer][seg];
  else
    return NULL;
}

inline CFTFiberHit* HodoAnalyzer::GetHitCFT(int layer, int seg) const
{
  if( seg>=0 && seg<CFTCont[layer].size() )
    return CFTCont[layer][seg];
  else
    return NULL;
}

inline Hodo2Hit * HodoAnalyzer::GetHitBGO( int i ) const
{
  if( i>=0 && i<BGOCont.size() )
    return BGOCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitPiV( int i ) const
{
  if( i>=0 && i<PiVCont.size() )
    return PiVCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitCH( int i ) const
{
  if( i>=0 && i<CHCont.size() )
    return CHCont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitTOF( int i ) const
{
  if( i>=0 && i<TOFCont.size() )
    return TOFCont[i];
  else
    return 0;
}

inline Hodo1Hit * HodoAnalyzer::GetHitAC( int i ) const
{
  if( i>=0 && i<ACCont.size() )
    return ACCont[i];
  else
    return 0;
}

inline Hodo1Hit * HodoAnalyzer::GetHitT0( int i ) const
{
  if( i>=0 && i<T0Cont.size() )
    return T0Cont[i];
  else
    return 0;
}

inline Hodo2Hit * HodoAnalyzer::GetHitLC( int i ) const
{
  if( i>=0 && i<LCCont.size() )
    return LCCont[i];
  else
    return 0;
}

inline HodoAnalyzer& HodoAnalyzer::GetInstance(){
  return *MySelf_;
}

#endif
