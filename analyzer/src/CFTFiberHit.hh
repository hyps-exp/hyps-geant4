/*
  CFTFiberHit.hh

  2012/1/24
*/

#ifndef CFTFIBER_HIT_H 
#define CFTFIBER_HIT_H
 
#include "HodoRawHit.hh"
#include "ThreeVector.hh"
#include <cstddef>

class RawData;
class CFTFLHit;

class CFTFiberHit
{
public:
  explicit CFTFiberHit( HodoRawHit *rhit, const char* name);
  virtual ~CFTFiberHit();

private:
  CFTFiberHit( const CFTFiberHit & );
  CFTFiberHit & operator = ( const CFTFiberHit & );

protected:
  HodoRawHit *raw_;
  bool Status_;
  double a1_, a2_, pe_;

  int multi_hit;
  std::vector <double> t1_, t2_, ct1_, ct2_, w_;

protected: // From FiberHit
  std::string DetectorName_;
  //double      position_;
  //double      offset_;
  int         pair_id_;

  double      phi_, r_;
  double      x_, y_;
  double      z0_, slope_;

  int         trackLayer_;

  std::vector<bool>   flJoin_;
  std::vector<CFTFLHit*> PtrCont_;

public:
  HodoRawHit * GetRawHit( void ) { return raw_; }
  bool calculate( void );

  bool status( void ) const { return Status_; } 

  double GetPhotonNum( void )  const { return pe_; }
  double GetAHiGain( void )  const { return a1_; }
  double GetALowGain( void ) const { return a2_; }
  double GetDEHiGain( void )  const { return a1_; }
  double GetDELowGain( void ) const { return a2_; }

  int    GetNumOfHit(void ) const {return raw_->SizeTdc1();}

  double GetT1( void )    const { return t1_.at(0); }
  double GetT1( int n )    const { return t1_.at(n); }
  double GetT2( void )  const { return t2_.at(0); }
  double GetT2( int n )  const { return t2_.at(n); }

  double GetCT1( void )    const { return ct1_.at(0); }
  double GetCT1( int n )    const { return ct1_.at(n); }
  double GetCT2( void )  const { return ct2_.at(0); }
  double GetCT2( int n )  const { return ct2_.at(n); }

  int DetectorId( void ) const { return raw_->DetectorId(); }
  int PlaneId( void ) const { return raw_->PlaneId(); }
  int SegmentId( void ) const { return raw_->SegmentId(); }

  virtual bool ReCalc( bool applyRecursively=false )
  { return calculate(); }

  void   SetDetectorName(const char* name){DetectorName_ = name;};

  // leading and trailing are raw value (not 0)
  double GetLeading() const {return raw_->GetTdc1(0);}
  double GetLeading(int n) const {return raw_->GetTdc1(n);}
  double GetTrailing() const {return raw_->GetTdc2(0);}
  double GetTrailing(int n) const {return raw_->GetTdc2(n);}

  // t_ is time for PID (should be 0 for pion)
  double GetTime() const {return t1_.at(0);}
  double GetTime(int n) const {return t1_.at(n);}

  // ct_ is corrected time after PHC
  double GetCTime() const {return ct1_.at(0);}
  double GetCTime(int n) const {return ct1_.at(n);}

  // In this class, a_ is width = leading - trailing
  // If there is no trailing data, -1 is inserted into width
  double GetWidth() const {return w_.at(0);}
  double GetWidth(int n) const {return w_.at(n);}

  // Treat BFT & SFT as 1 dimentional detector
  int    PairId() const {return pair_id_;}

  // make flag_[joined]
  void SetJoined(int m){flJoin_.at(m) = true;}
  bool Joined(int m) const {return flJoin_.at(m);}

  // Memorize the pointer of FLHit
  void RegisterHits(CFTFLHit* ptr){PtrCont_.push_back(ptr);}

  
  double GetPhi( void )  const { return phi_; }
  double GetR( void )  const { return r_; }
  double GetX( void )  const { return x_; }
  double GetY( void )  const { return y_; }
  double GetZ0( void )  const { return z0_; }
  double GetSlopeU( void )  const { return slope_; }

  int    GetTrackingLayer( void )  const { return trackLayer_; }


  static bool CompFiberHit(const CFTFiberHit* rLeft,
			   const CFTFiberHit* rRight);

private:
  void CleanUp();


};

inline bool CFTFiberHit::CompFiberHit(const CFTFiberHit* rLeft,
				      const CFTFiberHit* rRight)
{
  return rLeft->PairId() < rRight->PairId();
}

#endif
