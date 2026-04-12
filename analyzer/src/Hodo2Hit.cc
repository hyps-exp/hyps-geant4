/*
  Hodo2Hit.cc

  2012/1/24
*/

#include "Hodo2Hit.hh"
#include "ConfMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

Hodo2Hit::Hodo2Hit( HodoRawHit *rhit, int index)
  : raw_(rhit), Status_(false), index_(index), pid_(0)
{
}

Hodo2Hit::~Hodo2Hit()
{
}

double Hodo2Hit::DeltaE( void ) const
{
  return sqrt(fabs(a1_*a2_));
}

bool Hodo2Hit::calculate( void )
{
  static const std::string funcname = "[Hodo2Hit::calculate]";

  Status_=false;
  if( raw_->GetNumOfTdcHits()!=2 ) return Status_;

  int tdc1=raw_->GetTdc1(index_), tdc2=raw_->GetTdc2(index_);
  if( tdc1<0 || tdc2<0 ) return Status_;

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cerr << funcname << ": cannot get confManager" << std::endl;
 
    return Status_;
  }
  HodoParamMan *hodoMan = confMan->GetHodoParamManager();
  HodoPHCMan   *phcMan  = confMan->GetHodoPHCManager();
  if(!hodoMan){
    std::cerr << funcname << ": cannot get HodoParamManager" << std::endl; 
    return Status_;
  }

  int cid=raw_->DetectorId(), plid=raw_->PlaneId(),
    seg=raw_->SegmentId();
  int adc1=raw_->GetAdc1(), adc2=raw_->GetAdc2();

  if( !hodoMan->GetTime(cid,plid,seg,0,tdc1,t1_) ||
      !hodoMan->GetTime(cid,plid,seg,1,tdc2,t2_) ) return Status_;
  
  if( adc1>=0 ){
    if( !hodoMan->GetDe(cid,plid,seg,0,adc1,a1_) ) return Status_;
  }
  else
    a1_=0.;

  if( adc2>=0 ){
    if( !hodoMan->GetDe(cid,plid,seg,1,adc2,a2_) ) return Status_;
  }
  else
    a2_=0.;
  
  ct1_=t1_; ct2_=t2_;
  
  if( phcMan ){
    phcMan->doCorrection(cid,plid,seg,0,t1_,a1_,ct1_ );
    phcMan->doCorrection(cid,plid,seg,1,t2_,a2_,ct2_ );
  }

  return Status_=true;
} 

bool Hodo2Hit::calculateSimu( void )
{
  static const std::string funcname = "[Hodo2Hit::calculate]";

  Status_=false;
  if( raw_->GetNumOfTdcHits()!=2 ) return Status_;

  int tdc1=raw_->GetTdc1(index_), tdc2=raw_->GetTdc2(index_);
  if( tdc1<0 || tdc2<0 ) return Status_;

  int cid=raw_->DetectorId(), plid=raw_->PlaneId(),
    seg=raw_->SegmentId();
  int adc1=raw_->GetAdc1(), adc2=raw_->GetAdc2();

  t1_ = (double)tdc1/1000.;
  t2_ = (double)tdc2/1000.;

  if( adc1>=0 ){
    a1_ = (double)adc1/1000.;
  }
  else
    a1_=0.;

  if( adc2>=0 ){
    a2_ = (double)adc2/1000.;
  }
  else
    a2_=0.;
  
  ct1_=t1_; ct2_=t2_;
  
  return Status_=true;
} 


void Hodo2Hit::BGOPos(double *x, double *y) const
{
  int seg=raw_->SegmentId();
  int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  /* orignal BGO position*/
  /*
  double theta = 22.5+(double)UnitNum*45.;
  double x0 = RadiusOfBGOSurface+BGO_Y/2;
  double y0 = (double)(SegInUnit-1)*BGO_X;
  */

  /* new BGO position 2014/12/03*/

  double theta;
  double x0;
  double y0;

  if (SegInUnit==0 || SegInUnit==1 ) {
    theta = (double)UnitNum*45.;
    x0 = RadiusOfBGOSurface+BGO_Y/2;
    y0 = (double)((double)SegInUnit-0.5)*BGO_X;
  } else {
    theta = 22.5+(double)UnitNum*45.;
    x0 = RadiusOfBGOSurface2+BGO_Y/2;
    y0 = 0.;

  }


  *x = x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad);
  *y = x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad);
}

void Hodo2Hit::PiVPos(int layer, double *x, double *y) const
{
  int seg=raw_->SegmentId();

  if (seg>=0 && seg<=23) {
    int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
    int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

    double PiV_X = (RadiusOfPiVSurface + layer*PiV_Y)/RadiusOfPiVSurface*PiV_X0;
    
    //double theta = 22.5+(double)UnitNum*45.;
    /* change PIV pos with BGO 2014/12/03 */
    double theta = (double)UnitNum*45.;
    double x0 = RadiusOfPiVSurface + layer*PiV_Y +PiV_Y/2;
    double y0 = (double)(SegInUnit-1)*PiV_X;

    *x = x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad);
    *y = x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad);
  } else if (seg>=24 && seg<=31) {
    int UnitNum = seg-24;

    //double theta = (double)UnitNum*45.;
    /* change PIV pos with BGO 2014/12/03 */
    double PiV2_X = (RadiusOfPiV2Surface + layer*PiV2_Y)/RadiusOfPiV2Surface*PiV2_X0;

    double theta = 22.5+(double)UnitNum*45.;
    double x0 = RadiusOfPiV2Surface + layer*PiV2_Y + PiV2_Y/2;
    double y0 = 0.;

    *x = x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad);
    *y = x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad);

  } else {
    std::cerr << "Hodo2Hit::PiVPos No such segment " << seg << std::endl;
    *x = -999.;
    *y = -999.;
  }

}
