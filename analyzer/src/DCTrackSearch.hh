/*
  DCTrackSearch.hh

  2012/1/24
*/

#ifndef DCTrackSearch_h
#define DCTrackSearch_h 1

#include "DCAnalyzer.hh"

#include <vector>

struct DCPairPlaneInfo;
class DCPairHitCluster;
class DCLocalTrack;
class DCLTrackHit;
class MWPCCluster;
class CFTLocalTrack;

//BC1&2
int MWPCLocalTrackSearch( const DCHitContainer * HC,
			  std::vector <DCLocalTrack *> & TrackCont);

int MWPCLocalTrackSearch( const std::vector<std::vector<DCHitContainer> >& hcList,
			  std::vector <DCLocalTrack *> & trackCont);

int findUnusedHits(int n,
		   const DCHitContainer* src, 
		   std::vector<DCHitContainer>& dest);

//BC3&4 SDC1&2 : Pair Plane Tracking
int LocalTrackSearch( const DCHitContainer * HC,  
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> & trackCont,
		      int MinNumOfHits=6 );
/*
int LocalTrackSearchBcOut( const DCHitContainer * HC,  
			   const DCPairPlaneInfo * PpInfo,
			   int npp, std::vector <DCLocalTrack *> & trackCont,
			   int MinNumOfHits=6 );
*/
int LocalTrackSearchBcOutBh2T0( const DCHitContainer * HC,  
				const DCPairPlaneInfo * PpInfo,
				int npp, std::vector <DCLocalTrack *> & trackCont,
				int MinNumOfHits=6, int T0Seg=-1 );

int LocalTrackSearch( const std::vector<std::vector<DCHitContainer> > &hcAssemble,  
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> & TrackCont,
		      int MinNumOfHits=6 );
/*
int LocalTrackSearchCFT( const CFTFiberClusterContainer * HC,
			 int NPhiPlane, int NUPlane,
			 std::vector <CFTLocalTrack *> &TrackCont,
			 int MinNumOfHits=3 );
*/

int LocalTrackSearchCFT( const CFTFiberClusterContainer * HC,
			 int NPlane,
			 std::vector <CFTLocalTrack *> &TrackCont,
			 int MinNumOfHitsPhi=3, int MinNumOfHitsU=3 );

int LocalTrackSearchCFT_2nd( const CFTFiberClusterContainer * HC,
			     std::vector <ThreeVector> CFTVtxCont,
			     int NPlane,
			     std::vector <CFTLocalTrack *> &TrackCont,
			     int MinNumOfHitsPhi=2,  int MinNumOfHitsU=2);

//Added by Miwa
int LocalTrackSearchBh2T0( const std::vector<std::vector<DCHitContainer> > &hcAssemble,  
		      const DCPairPlaneInfo * PpInfo,
		      int npp, std::vector <DCLocalTrack *> & TrackCont,
		      int MinNumOfHits=6, int T0Seg=-1 );

int LocalTrackSearchSdcInFiber( const DCHitContainer * HC,  
				const DCPairPlaneInfo * PpInfo,
				int npp, std::vector <DCLocalTrack *> & trackCont,
				int MinNumOfHits=6 );

bool checkHitCombSdcIn(int  *hitWire);

//BC3&4 SDC1&2 : XUV Tracking 
int LocalTrackSearchVUX( const DCHitContainer * HC,  
			 const DCPairPlaneInfo * PpInfo,
			 int npp, std::vector <DCLocalTrack *> & TrackCont,
			 int MinNumOfHits=6 );

//SDC3&4 Old routine (Only used for study)
//int SdcOutLocalTrackSearch( const DCHitContainer * HC,
//			    std::vector <DCLocalTrack *> &TrackCont,
//			    std::vector <DCLocalTrack *> &TrackCont1,
//			    std::vector <DCLocalTrack *> &TrackCont2);

//SDC3&4 (**Calculate all hit combinations -> Too slow !!**)
int SdcOutLocalTrackSearch( const DCHitContainer * HC,
			    std::vector <DCLocalTrack *> &TrackCont);

int LocalTrackSearchBcOutSdcIn( const DCHitContainer * BcHC,  
                                 const DCPairPlaneInfo * BcPpInfo,
                                 const DCHitContainer * SdcHC,  
                                 const DCPairPlaneInfo * SdcPpInfo,
                                 int BcNpp, int SdcNpp,
                                 std::vector <DCLocalTrack *> &TrackCont,
                                 int MinNumOfHits=18 );

bool MakePairPlaneHitCluster( const DCHitContainer & HC1,
			      const DCHitContainer & HC2,
			      double CellSize,
			      std::vector <DCPairHitCluster *> & Cont );

bool MakeUnPairPlaneHitCluster( const DCHitContainer & HC,
				std::vector <DCPairHitCluster *> & Cont );

bool MakeMWPCPairPlaneHitCluster( const DCHitContainer & HC,  
				  std::vector <DCPairHitCluster *> & Cont );

DCLocalTrack *MakeTrack( const std::vector < std::vector <DCPairHitCluster *> > &CandCont,
			 const int *combination );

CFTLocalTrack * MakeTrack(  const CFTFiberClusterContainer * HC,
			   const int *combination, int NPhiPlane );

CFTLocalTrack * MakeTrack_2nd(  const CFTFiberClusterContainer*  HC,
				const std::vector <ThreeVector> CFTVtxCont,
				const int *combination, int NPlane );

DCLocalTrack *MakeTrackVUX( const std::vector < std::vector <DCPairHitCluster *> > &CandCont,
			    const int *combination, const DCLocalTrack *tp1  );


std::vector< std::vector<int> > makeindex( int ndim, const int *index1 ); 

std::vector< std::vector<int> > makeindex_SdcOut( int ndim_org, int minimumHit, int ndim, const int *index1 ); 

std::vector< std::vector<int> > makeindex_BcIn( int ndim_org, int minimumHit, int ndim, const int *index1 ); 

//XUV Tracking 
bool MakePairPlaneHitClusterVUX( const DCHitContainer & HC1,
				 const DCHitContainer & HC2,
				 double CellSize,
				 std::vector <DCPairHitCluster *> & Cont );

std::vector< std::vector<int> > makeindex_VXU( int ndim, int maximumHit, const int *index1 ); 
std::vector< std::vector<int> > makeindex_SdcOut_below( int ndim_org, int maximumHit, int ndim, const int *index1 ); 

#endif
