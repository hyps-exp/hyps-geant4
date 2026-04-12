/*
  CFTLocalTrack.hh

  2012/1/24
*/

#ifndef CFTLocalTrack_h
#define CFTLocalTrack_h 1

#include <vector>
#include <functional>
#include "ThreeVector.hh"
#include "CFTFiberCluster.hh"

class CFTFiberCluster;
class DCAnalyzer;

class CFTLocalTrack
{
public:
  explicit CFTLocalTrack();
  ~CFTLocalTrack();
private:
  CFTLocalTrack( const CFTLocalTrack & );
  CFTLocalTrack & operator = ( const CFTLocalTrack & );

private:
  std::vector <CFTFiberCluster *> hitArray;
  std::vector <CFTFiberCluster *> hitArrayU;
  ThreeVector CFTVtx_;

  int    xyFitFlag_;
  double Axy_, Bxy_;
  double Az_, Bz_;
  double chisqrXY_;
  double chisqrZ_;

  // three dementional track
  int    zTrackFlag_;
  double u0_, v0_, x0_, y0_;
  double ur0_, vr0_, f0_, z0_;
  ThreeVector Dir_, Pos0_;

  // dE informaiton
  double TotalDEHiGain_,TotalDELowGain_;
  double MaxDEHiGain_,MaxDELowGain_;
  double pathlength_;

public:
  void AddHit( CFTFiberCluster *hitp ) { hitArray.push_back( hitp ); }
  void AddCFTVtx(ThreeVector vtx) {CFTVtx_=vtx;}
  bool DoFitXY( void );
  void AddHitU( CFTFiberCluster *hitp ) { hitArrayU.push_back( hitp ); }
  bool DoFitZTrack();
  bool DoFitXY_wVtx( void );
  bool CheckPhi( void );
  bool CheckPhi_1st( void );
  bool CheckEdep( void );
  bool DoFitZTrack_wVtx();
  bool DoFit( void );
  bool FindZTrack();
  bool SetCalculatedValue();

  std::size_t GetNHit( void ) const {return hitArray.size(); }
  std::size_t GetNHitU( void ) const {return hitArrayU.size(); }
  CFTFiberCluster * GetHit( std::size_t nth ) const;
  CFTFiberCluster * GetHitOfLayerNumber( int lnum ) const;
  CFTFiberCluster * GetHitU( std::size_t nth ) const;
  int    GetFirstLayerPhi( void ) const;
  int    GetFirstLayerUV( void ) const;
  double GetChiSquareXY( void ) const { return chisqrXY_; }
  double GetChiSquareZ( void ) const { return chisqrZ_; }

  /*
  void SetAv( double Av) { Av_=Av; }
  void SetAx( double Ax) { Ax_=Ax; }
  void SetAu( double Au) { Au_=Au; }
  void SetChiv( double Chiv) { Chiv_=Chiv; }
  void SetChix( double Chix) { Chix_=Chix; }
  void SetChiu( double Chiu) { Chiu_=Chiu; }
  */
  double GetX0( void ) const { return x0_; }
  double GetY0( void ) const { return y0_; }
  double GetU0( void ) const { return u0_; }
  double GetV0( void ) const { return v0_; }


  double GetChiSquare( void ) const { return chisqr_; }
  double GetX( double z ) const { return x0_+u0_*z; } 
  double GetY( double z ) const { return y0_+v0_*z; } 
  int    GetXYFitFlag( void ) const { return xyFitFlag_; }
  double GetAxy( void ) const { return Axy_; }
  double GetBxy( void ) const { return Bxy_; }

  ThreeVector GetPos0( void ) const { return Pos0_; }
  ThreeVector GetDir( void )  const { return Dir_; }

  double GetPhiFromR( double r ) const { return f0_+ur0_*r; } 
  double GetZFromR( double r ) const { return z0_+vr0_*r; } 

  int   GetZTrackFlag( void ) const { return zTrackFlag_; }

  bool   GetPathLengthInFiber(double r, double phi, double *len);
  bool   GetCrossPointR(double r, double *phi1, double *phi2);
  bool   GetCrossPointR(double r, double *phi);
  double calcPhi(double x, double y);
  double CalculateZpos(double phi, CFTFiberCluster *cl);
  double CalculateUVpos(double phi, CFTFiberCluster *cl);

  bool   CalcNormalizedDE( void ); 

  double    TotalDEHiGain() const { return TotalDEHiGain_;};
  double    MaxDEHiGain() const { return MaxDEHiGain_;};
  double    TotalDELowGain() const { return TotalDELowGain_;};
  double    MaxDELowGain() const { return MaxDELowGain_;};
  
  double    GetTotalPathLength() const  { return pathlength_;};
  double    NormalizedTotalDEHiGain (void) const {return TotalDEHiGain_/pathlength_;};
  double    NormalizedTotalDELowGain (void) const {return TotalDELowGain_/pathlength_;};
  double    NormalizedMaxDEHiGain (void) const {return MaxDEHiGain_/pathlength_;};
  double    NormalizedMaxDELowGain (void) const {return MaxDELowGain_/pathlength_;};

  bool GetStatus( void ) const { return status_; } 
  bool GoodForTracking( void ) const { return gftstatus_; }
  bool GoodForTracking( bool status )
  { bool ret=gftstatus_; gftstatus_=status; return ret; } 
  //bool ReCalc( bool ApplyRecursively=false );  

private:
  bool status_;
  //double x0_, y0_, u0_, v0_;
  double a_,b_;
  double chisqr_;
  bool gftstatus_;
};

#if 0
struct CFTLTrackComp 
  : public std::binary_function <CFTLocalTrack *, CFTLocalTrack *, bool>
{
  bool operator()( const CFTLocalTrack * const p1, 
		   const CFTLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquareXY(),chi2=p2->GetChiSquareXY();
    if( (n1>n2+1) ){
      return true;
    }
    else if( (n2>n1+1)  ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }
  }
};
#endif

struct CFTLTrackComp 
  : public std::binary_function <CFTLocalTrack *, CFTLocalTrack *, bool>
{
  bool operator()( const CFTLocalTrack * const p1, 
		   const CFTLocalTrack * const p2 ) const
  {
    double thr_chisqr=5.;

    int n1=p1->GetNHit()+p1->GetNHitU();
    int n2=p2->GetNHit()+p2->GetNHitU();

    double chi1=p1->GetChiSquareXY()+p1->GetChiSquareZ();
    double chi2=p2->GetChiSquareXY()+p2->GetChiSquareZ();

    if( (n1>n2) ){
      if (chi1<thr_chisqr)
	return true;
      else
	return false;
    }
    else if( (n2>n1)  ){
      if (chi2<thr_chisqr)
	return false;
      else 
	return true;
    }
    else{
      return (chi1<=chi2);
    }
  }
};

#endif
