//----------------------------------------------------------------------------
// CFTFLHit
//  This class corresponds to the DCLTrackHit
//----------------------------------------------------------------------------

#ifndef CFTFLHit_H
#define CFTFLHit_H

#include"CFTFiberHit.hh"

class CFTFLHit{
private:
  CFTFiberHit* Hit_;
  int       nth_Hit_;

public:
  CFTFLHit(CFTFiberHit* ptr, int index):
    Hit_(ptr), nth_Hit_(index)
  {
    Hit_->RegisterHits(this);
    Hit_->SetJoined(index);
  }

  // photon number of High Gain
  double GetPhotonNum( )  const { return Hit_->GetPhotonNum();}

  // dE of High Gain
  double GetDEHiGain( void )  const { return Hit_->GetDEHiGain();}

  // dE of Low Gain
  double GetDELowGain( void ) const { return Hit_->GetDELowGain();}

  // leading and trailing are raw value (not 0)
  double GetLeading() const {return Hit_->GetLeading(nth_Hit_);}
  double GetTrailing() const {return Hit_->GetTrailing(nth_Hit_);}

  // t_ is time for PID (should be 0 for pion)
  double GetTime() const {return Hit_->GetTime(nth_Hit_);}

  // ct_ is corrected time after PHC
  double GetCTime() const {return Hit_->GetCTime(nth_Hit_);}

  // In this class, a_ is width = leading - trailing
  // If there is no trailing data, -1 is inserted into width
  double GetWidth() const {return Hit_->GetWidth(nth_Hit_);}

  // Return local position.
  // Note : BFT is a position detector, there is no z information.
  //        Z position should be written in the geometry file.
  //double GetPosition() const {return Hit_->GetPosition();}

  // Treat BFT as 1 dimentional detector
  int    PairId() const {return Hit_->PairId();}

  // Join in Cluster
  bool   Joined() const {return Hit_->Joined(nth_Hit_);}


  double GetPhi( void )  const { return Hit_->GetPhi();}
  double GetR( void )  const { return Hit_->GetR();}
  double GetX( void )  const { return Hit_->GetX();}
  double GetY( void )  const { return Hit_->GetY();}
  double GetZ0( void )  const { return Hit_->GetZ0();}
  double GetSlopeU( void )  const { return Hit_->GetSlopeU();}

  int    GetTrackingLayer( void )  const { return Hit_->GetTrackingLayer();}

  // Dump
  void   Dump(){Debug();}

  friend class CFTFiberHit;

private:
  ~CFTFLHit(){;};

  CFTFLHit();
  CFTFLHit(const CFTFLHit& object);
  CFTFLHit& operator =(const CFTFLHit& object);

  void Debug(){
    std::cout << "plid " << Hit_->PairId() << " ";
    //std::cout << "Pos " << GetPosition() << " ";
    std::cout << "Time " << GetCTime() << std::endl;
  }
};

#endif
