/*
  CFTFiberHit.cc

  2012/1/24
*/

#include "CFTFiberHit.hh"
#include"CFTFLHit.hh"

#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include"DCGeomMan.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);


CFTFiberHit::CFTFiberHit( HodoRawHit *rhit, const char* name)
  : raw_(rhit), Status_(false), DetectorName_(name), pair_id_(0)
{
}

CFTFiberHit::~CFTFiberHit()
{
  CleanUp();
}

// CleanUp
void CFTFiberHit::CleanUp(){
  int NofPtr = PtrCont_.size();
  for(int i = 0; i<NofPtr; ++i){
    delete PtrCont_.at(i);
    PtrCont_.at(i) = NULL;
  }
  PtrCont_.clear();
}

bool CFTFiberHit::calculate( void )
{
  static const std::string funcname = "[CFTFiberHit::calculate]";

  if(Status_){
    std::cerr << "#E " << funcname 
	      << " Already calculated" << std::endl;
    return false;
  }

  if("" == DetectorName_){
    std::cerr << "#E " << funcname 
	      << " DetectorName is NULL" << std::endl;
    return Status_;
  }

  // Detector information
  int cid  = raw_->DetectorId();
  int plid = raw_->PlaneId();
  int seg  = raw_->SegmentId();

  pair_id_ = raw_->SegmentId();

  // For simulation
  /*
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cerr << "#E " << funcname
	      << "Cannot get confManager"
	      << std::endl; 
    return Status_;
  }
  HodoParamMan *hodoMan = confMan->GetHodoParamManager();
  HodoPHCMan   *phcMan  = confMan->GetHodoPHCManager();
  if(!hodoMan){
    std::cerr << "#E " << funcname
	      << "Cannot get HodoParamManager"
	      << std::endl; 
    return Status_;
  }
  */

  // hit information
  multi_hit = raw_->SizeTdc1();
  int UorD=0; // BFT is allways U
  int fl_EASIROC = (-1 != raw_->GetTdc2());

  for(int m = 0; m<multi_hit; ++m){
    int leading  = raw_->GetTdc1(m); // leading
    int width    = -1;
    if(fl_EASIROC){
      int trailing = raw_->GetTdc2(m); // trailing
      width = leading - trailing; // width
    }
    w_.push_back(width);

    // Time calibration
    if( leading<0 ){
      std::cerr << "#E " <<  funcname 
		<< "No valid TDC\n" 
		<< "TDC = " << leading
		<< std::endl;
      return Status_;
    }
  
    // time is overwriten
    double time = 0;
    /*
    if( !hodoMan->GetTime(cid, plid, seg, UorD, leading, time) ){
      std::cerr << "#E " << funcname
		<< "Something wrong at GetTime"
		<< std::endl; 
      std::cerr << cid << " "
		<< plid << " "
		<< seg << " "
		<< UorD << " "
		<< leading << " "
		<< time << std::endl;
      
      return Status_;
    }
    */
    time = (double)leading/1000.; //for simulation
    t1_.push_back(time);

    // ctime is overwriten
    double ctime = time;
    /*
    if( phcMan ){
      phcMan->doCorrection(cid, plid, seg, UorD, time, width, ctime );
    }
    */
    ct1_.push_back(ctime);

    flJoin_.push_back(false);
  }

  int adc1=raw_->GetAdc1(), adc2=raw_->GetAdc2();

  if( adc1>=0 ){
    /*
    if( !hodoMan->GetPhotonNum(cid,plid,seg,0,adc1,pe_) ) return Status_;
    if( !hodoMan->GetDe(cid,plid,seg,2,adc1,a1_) ) return Status_;
    */
    pe_ = (double)adc1/1000.;
    a1_ = (double)adc1/1000.;
  }
  else
    a1_=-999.;

  if( adc2>=0 ){
    /*
    if( !hodoMan->GetDe(cid,plid,seg,1,adc2,a2_) ) return Status_;
    */
    a2_ = (double)adc2/1000.;
  }
  else
    a2_=-999.;
  

  DCGeomMan& geom = DCGeomMan::GetInstance();
  int lnum = geom.GetDetectorId(DetectorName_);

  trackLayer_ = lnum;
  r_ = geom.GetLocalZ(lnum);
  if (pair_id_%2 == 0)
    r_ -= 0.4;
  else
    r_ += 0.4;

  if (plid == CFT_PHI1 || plid == CFT_PHI2 || plid == CFT_PHI3 || plid == CFT_PHI4) {
    phi_ = geom.calcWirePosition(lnum, pair_id_);
    x_ = r_ * cos(phi_*Deg2Rad);
    y_ = r_ * sin(phi_*Deg2Rad);
  } else if (plid == CFT_U1 || plid == CFT_V2 || plid == CFT_U3 || plid == CFT_V4) {
    z0_ =  geom.calcWirePosition(lnum, pair_id_);
    slope_ =  geom.GetTiltAngle(lnum);
  }

  return Status_=true;
} 
