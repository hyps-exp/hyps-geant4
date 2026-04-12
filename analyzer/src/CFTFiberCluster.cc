#include"CFTFiberHit.hh"
#include"CFTFLHit.hh"
#include"CFTFiberCluster.hh"

#include<string>
#include<math.h>

#define DEBUG 0

static const std::string MyName = "CFTFiberCluster";

// Constructor
CFTFiberCluster::CFTFiberCluster()
  : csize_(0), xyFitFlag_(-1), phical_(-1)
{
  Initializer();
}

// Destructor
CFTFiberCluster::~CFTFiberCluster(){
  CleanUp();
}

// Initializer
void CFTFiberCluster::Initializer(){
  for(int i = 0; i<sizeFlagsFiber; ++i){
    flag_[i] = false;
  }
}

// CleanUp
void CFTFiberCluster::CleanUp(){
  HitAssemblage_.clear();
}

// GetHit
CFTFLHit *CFTFiberCluster::GetHit( int i ) const
{
  return HitAssemblage_.at(i);
}

// Calculate
bool CFTFiberCluster::calculate(){
  static const std::string MyFunc = "::calculate ";
  if(flag_[Initialized]){
    std::cerr << "#E " << MyName << MyFunc
	      << "Already initialied" << std::endl;
    return false;
  }

  if(0 == (csize_ = HitAssemblage_.size())){
    std::cerr << "#E " << MyName << MyFunc
	      << "No CFTFiberHit in local container" << std::endl;
    return false;
  }

  int width = (int)HitAssemblage_.at(0)->GetWidth();
  if(-1 == width){
    calc_wo_width();
  }else{
    calc_w_width();
  }

  double total_pe = -10., max_pe = -10.;
  double total_dE1 = -10., max_dE1 = -10.;
  double total_dE2 = -10., max_dE2 = -10.;
  double mean_pair_id = -1.;

  double mean_x, mean_y, mean_phi, mean_r, mean_z0;

  total_pe = HitAssemblage_.at(0)->GetPhotonNum();
  max_pe = HitAssemblage_.at(0)->GetPhotonNum();

  total_dE1 = HitAssemblage_.at(0)->GetDEHiGain();
  max_dE1 = HitAssemblage_.at(0)->GetDEHiGain();

  total_dE2 = HitAssemblage_.at(0)->GetDELowGain();
  max_dE2 = HitAssemblage_.at(0)->GetDELowGain();

  mean_pair_id = (double)HitAssemblage_.at(0)->PairId();

  MeanX_ = HitAssemblage_.at(0)->GetX();
  MeanY_ = HitAssemblage_.at(0)->GetY();
  MeanPhi_ = HitAssemblage_.at(0)->GetPhi();
  MeanR_ = HitAssemblage_.at(0)->GetR();
  MeanZ0_ = HitAssemblage_.at(0)->GetZ0();
  Slope_ =  HitAssemblage_.at(0)->GetSlopeU();
  trackLayer_ =  HitAssemblage_.at(0)->GetTrackingLayer();

  for(int i = 1; i<csize_; ++i){
    double pe  = HitAssemblage_.at(i)->GetPhotonNum();
    double dE1 = HitAssemblage_.at(i)->GetDEHiGain();
    double dE2 = HitAssemblage_.at(i)->GetDELowGain();

    total_pe  += pe;
    total_dE1 += dE1;
    total_dE2 += dE2;

    MeanX_   += HitAssemblage_.at(i)->GetX();
    MeanY_   += HitAssemblage_.at(i)->GetY();
    MeanPhi_ += HitAssemblage_.at(i)->GetPhi();
    MeanR_   += HitAssemblage_.at(i)->GetR();
    MeanZ0_  += HitAssemblage_.at(i)->GetZ0();

    mean_pair_id += (double)HitAssemblage_.at(i)->PairId();

    if (pe > max_pe)
      max_pe = pe;
    if (dE1 > max_dE1)
      max_dE1 = dE1;
    if (dE2 > max_dE2)
      max_dE2 = dE2;
    
  }

  TotalPhotonNum_ = total_pe;
  MaxPhotonNum_   = max_pe;
  TotalDEHiGain_  = total_dE1;
  MaxDEHiGain_    = max_dE1;
  TotalDELowGain_  = total_dE2;
  MaxDELowGain_    = max_dE2;

  mean_pair_id /= double(csize_);
  MeanPairId_ = mean_pair_id;

  MeanX_    /= double(csize_);
  MeanY_    /= double(csize_);
  MeanPhi_  /= double(csize_);
  MeanR_    /= double(csize_);
  MeanZ0_   /= double(csize_);

  flag_[Initialized] = true;

#if DEBUG
  Debug();
#endif

  return true;
}

// calc_wo_width
// Time closest to 0 is selected as mean time
void CFTFiberCluster::calc_wo_width(){
  double mx = 0., bt = 0.;

  //mx = HitAssemblage_.at(0)->GetPosition();
  bt = HitAssemblage_.at(0)->GetCTime();

  for(int i = 1; i<csize_; ++i){
    //mx += HitAssemblage_.at(i)->GetPosition();
    double t = HitAssemblage_.at(i)->GetCTime();

    if(fabs(t) < fabs(bt)){
      bt = t;
    }    

  }

  //mx /= double(csize_);
  
  MeanTime_=bt; bestWidth_=-1; //MeanPos_=mx;

}

// calc_w_width
// Time closest to 0 is selected as mean time
void CFTFiberCluster::calc_w_width(){
  double mx = 0., mt = 0., bw = 0.;

  //mx = HitAssemblage_.at(0)->GetPosition();
  mt = HitAssemblage_.at(0)->GetCTime();
  bw = HitAssemblage_.at(0)->GetWidth();

  for(int i = 1; i<csize_; ++i){
    //mx += HitAssemblage_.at(i)->GetPosition();
    double t = HitAssemblage_.at(i)->GetCTime();

    if(fabs(t) < fabs(mt)){
      mt = t;
    }    

    if(bw < HitAssemblage_.at(i)->GetWidth()){
      bw = HitAssemblage_.at(i)->GetWidth();
    }

  }

  //mx /= double(csize_);
  
  MeanTime_=mt; bestWidth_=bw; //MeanPos_=mx;
}

// Debug
void CFTFiberCluster::Debug(){
  std::cout << "Used hit\n";
  for(int i = 0; i<csize_; ++i){
    HitAssemblage_.at(i)->Dump();
  }
}
