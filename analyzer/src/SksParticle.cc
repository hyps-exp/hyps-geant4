/*
  SksParticle.cc

*/

#include "SksParticle.hh"
#include "SksTrack.hh"
#include "HodoCluster.hh"
#include "Kinematics.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"

#include <TRandom.h>

ThreeVector SksParticle::Momentum( void ) const
{
  return Track_->PrimaryMomentum();
}

double SksParticle::Polarity( void ) const
{
  return Track_->polarity();
}

ThreeVector SksParticle::Position( void ) const
{
  return Track_->PrimaryPosition();
}

double SksParticle::PathLengthToTOF( void ) const
{
  return Track_->PathLengthToTOF();
}

/*
double SksParticle::PathLengthTotal( void ) const
{
  return Track_->PathLengthTotal();
}
*/

double SksParticle::MassSquare( ) const
{
  double ttof=Tof_->CMeanTime();

  double time_reso = 0.150; //ns 
  double sigma_time = gRandom->Gaus(0, time_reso);
  ttof += sigma_time;
return ::MassSquare( Momentum().mag(), PathLengthToTOF(), 
		       ttof);
}

double SksParticle::yTof() const
{
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );
  return Out_->GetY(zTof);
}

double SksParticle::xTof() const
{
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );
  return Out_->GetX(zTof);
}
