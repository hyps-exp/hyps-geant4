//-----------------------------------------------------------------------------
// Class for CFT cluster
//-----------------------------------------------------------------------------

#ifndef CFTFIBERCLUSTER
#define CFTFIBERCLUSTER

#include<vector>

class CFTFLHit;

class CFTFiberCluster{
protected:
  typedef std::vector<CFTFLHit*> HitContainer;
  
  HitContainer HitAssemblage_;
  int          csize_;

  double       TotalPhotonNum_;
  double       MaxPhotonNum_;

  double       TotalDEHiGain_;
  double       MaxDEHiGain_;

  double       TotalDELowGain_;
  double       MaxDELowGain_;

  double       MeanTime_;
  double       bestWidth_;
  //double       MeanPos_;
  enum         FlagsFiber{Initialized, gfastatus,
			  sizeFlagsFiber};
  bool         flag_[sizeFlagsFiber];

  double       MeanPairId_;

  double       MeanX_, MeanY_;
  double       MeanPhi_, MeanR_;
  double       MeanZ0_, Slope_;

  double       xcal_, ycal_, zcal_, phical_;
  double       pathlength_;

  int       trackLayer_;
  bool      belongTrack_;

  int       xyFitFlag_;

public:
  bool      calculate();
  void      push_back(CFTFLHit* ptr){ HitAssemblage_.push_back(ptr); };

  // Return # of used hit to create cluster
  int       ClusterSize() const {return csize_;};
  
  // Return total photon num
  double    TotalPhotonNum() const { return TotalPhotonNum_;};

  // Return max photon num
  double    MaxPhotonNum() const { return MaxPhotonNum_;};

  // Return total dE high gain
  double    TotalDEHiGain() const { return TotalDEHiGain_;};

  // Return max dE high gain
  double    MaxDEHiGain() const { return MaxDEHiGain_;};

  // Return total dE low gain
  double    TotalDELowGain() const { return TotalDELowGain_;};

  // Return max dE low gain
  double    MaxDELowGain() const { return MaxDELowGain_;};

  // Return most closest value to 0
  double    CMeanTime() const {return MeanTime_;};

  // Return best width
  double    Width() const {return bestWidth_;};

  // Return mean local position
  //double    MeanPosition() const {return MeanPos_;};

  // Return mean pair Id
  double    MeanPairId() const {return MeanPairId_;};
  
  CFTFLHit*    GetHit(int i) const;
  bool      GoodForAnalysis() const {return flag_[gfastatus];};
  bool      GoodForAnalysis(bool status)
  {bool ret = flag_[gfastatus]; flag_[gfastatus]=status; return ret;};
  
  bool      ReCalc( bool applyRecusively=false );

  int       GetTrackingLayer() const {return trackLayer_;};
  void  setFlags() { belongTrack_=true; }
  void  clearFlags() { belongTrack_=false; }
  bool  showFlags() const { return belongTrack_; }

  double    GetX() const {return MeanX_;};
  double    GetY() const {return MeanY_;};
  double    GetPhi() const {return MeanPhi_;};
  double    GetR() const {return MeanR_;};
  double    GetZ0() const {return MeanZ0_;};
  double    GetSlopeU() const {return Slope_;};

  void      SetCalPosition(double x, double y) {xcal_=x; ycal_=y;};
  void      SetCalZ(double z) {zcal_=z;};
  void      SetCalPhi(double phi) {phical_=phi;};
  void      SetPathLength(double len) {pathlength_ = len;};

  void      SetXYFitFlag(int flag) {xyFitFlag_ = flag;};
  double    GetXcal (void) const {return xcal_;};
  double    GetYcal (void) const {return ycal_;};
  double    GetZcal (void) const {return zcal_;};
  double    GetCalPhi (void) const {return phical_;};
  double    GetPathLength (void) const {return pathlength_;};

  double    NormalizedTotalDEHiGain (void) const {return TotalDEHiGain_/pathlength_;};
  double    NormalizedTotalDELowGain (void) const {return TotalDELowGain_/pathlength_;};
  double    NormalizedMaxDEHiGain (void) const {return MaxDEHiGain_/pathlength_;};
  double    NormalizedMaxDELowGain (void) const {return MaxDELowGain_/pathlength_;};

  inline double  GetResidual (void);

  CFTFiberCluster();
  virtual ~CFTFiberCluster();

private:
  void calc_wo_width();
  void calc_w_width();

  void Initializer();
  void CleanUp();
  void Debug();

  CFTFiberCluster(const CFTFiberCluster& object);
  CFTFiberCluster& operator =(const CFTFiberCluster& object);

};

inline double CFTFiberCluster::GetResidual()
{
  if (xyFitFlag_==0) {
    return MeanY_ - ycal_;
  } else if (xyFitFlag_==1) {
    return MeanX_ - xcal_;
  } else {
    std::cout << "CFTFiberCluster::GetResidual : invalid xyFitFlag=" << xyFitFlag_ << std::endl;
    return -999.;

  }
}

#endif
