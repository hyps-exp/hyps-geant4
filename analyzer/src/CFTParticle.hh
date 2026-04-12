#ifndef CFTParticle_h

#define CFTParticle_h 1

#include <vector>
#include "CFTLocalTrack.hh"
#include "ThreeVector.hh"

class CFTLocalTrack;
class Hodo2Hit;


class CFTParticle
{
public:
  CFTParticle(CFTLocalTrack *track);
  ~CFTParticle();

private:
  CFTLocalTrack* Track_;
  std::vector <Hodo2Hit *> BGOCont_;
  std::vector <Hodo2Hit *> PiVCont_;
  double  BGO_E_;
  double  PiV_E_;
  double  FiberTotal_E_;
  double  FiberMax_E_;
  double  NormalizedFiberTotal_E_;
  double  NormalizedFiberMax_E_;
  double  TotalE_;
  double  PathLength_;
  double  Mass_;
  double  BGO_Ereso; //1 MeV at 80 MeV
  double  BGO_Ereso2; // Constant term
  double  Fiber_Ereso; //0.4 MeV at 4 MeV
  ThreeVector CFTVtx_;
public:
  void AddBGOHit(Hodo2Hit* hit) {BGOCont_.push_back(hit);}
  void AddPiVHit(Hodo2Hit* hit) {PiVCont_.push_back(hit);}
  CFTLocalTrack * GetTrack() {return Track_;}
  int NHitBGO() { return BGOCont_.size();}
  Hodo2Hit * GetBGOHit(int i); 
  int NHitPiV() { return PiVCont_.size();}
  Hodo2Hit * GetPiVHit(int i); 
  bool Calculate();
  ThreeVector GetPos0 () { return Track_->GetPos0(); }
  ThreeVector GetDir () { return Track_->GetDir(); }
  double GetMass() { return Mass_; }
  double GetTotalE() { return TotalE_; }
  double GetBGO_E() { return BGO_E_; }
  double GetPiV_E() { return PiV_E_; }
  double GetFiberTotal_E() { return FiberTotal_E_; }
  double GetNormFiberTotal_E() { return NormalizedFiberTotal_E_; }
  bool checkProton( double BGO_E, double *delta);
  bool checkPi( double BGO_E);
  void   SetCFTVtx(ThreeVector vtx) {CFTVtx_=vtx;}
  ThreeVector GetCFTVtx() {return CFTVtx_;}
};


#endif
