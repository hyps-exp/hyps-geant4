/*
  CFTTrackingEffMan.cc

  2012/1/24
*/

#include "CFTTrackingEffMan.hh"

#include <TRandom.h>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>

const int MaxChar = 200;

CFTTrackingEffMan *CFTTrackingEffMan::CftEffMan_=0;

CFTTrackingEffMan::CFTTrackingEffMan()
{


}

CFTTrackingEffMan::~CFTTrackingEffMan()
{}


CFTTrackingEffMan & CFTTrackingEffMan::GetInstance( void )
{
  if( !CftEffMan_ ){
    CftEffMan_ = new CFTTrackingEffMan();
  }
  return *CftEffMan_;
}


bool CFTTrackingEffMan::Initialize( void )
{
  static const std::string funcname = "[CFTTrackingEffMan::Initialize]";
  char str[MaxChar];
  char cname[MaxChar];

  FILE *fp;

  if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      int index_theta, index_vtz;
      double P0, P1, P2, P0_simu, P1_simu, P2_simu;

      if( sscanf(str, "%d %d %lf %lf %lf %lf %lf %lf", 
		 &index_theta, &index_vtz, &P0_simu, &P1_simu, &P2_simu, &P0, &P1, &P2)	== 8 ){
	p0[index_vtz][index_theta] = P0;
	p1[index_vtz][index_theta] = P1;
	p2[index_vtz][index_theta] = P2;
	p0_simu[index_vtz][index_theta] = P0_simu;
	p1_simu[index_vtz][index_theta] = P1_simu;
	p2_simu[index_vtz][index_theta] = P2_simu;
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << funcname << " Initialization finished." << std::endl;

  return true;
}


bool CFTTrackingEffMan::CheckTracking_P(double vtz, double theta, double mom)
{
  int interval = 300/NdivVertex;
  int index_vtx = (int)(vtz+150)/interval;
  if (index_vtx<0)
    index_vtx=0;
  else if (index_vtx>=NdivVertex)
    index_vtx=NdivVertex-1;

  int index = (int)theta;
  if (index>=89)
    index=88;
  else if (index<0)
    index=0;

  double eff = p0[index_vtx][index]*(1.0 - 1./(1+exp((mom-p1[index_vtx][index])/p2[index_vtx][index])));
  double eff_simu = p0_simu[index_vtx][index]*(1.0 - 1./(1+exp((mom-p1_simu[index_vtx][index])/p2_simu[index_vtx][index])));

  //std::cout << eff << ", " << eff_simu << std::endl;

  if (eff_simu < 0.0001)
    return false;
  else if (eff > eff_simu)
    return true;
  else {
    double ratio = eff/eff_simu;
    double val = gRandom->Uniform();
    if (val<ratio)
      return true;
    else
      return false;
  }

  return false;


}

bool CFTTrackingEffMan::CheckTracking_Pi(double vtz, double theta, double mom)
{
  const double PionMass        = 0.1395701;
  const double ProtonMass      = 0.93827200;

  double mom_pi = mom;
  double E_pi = sqrt(mom_pi*mom_pi+PionMass*PionMass);
  double beta = mom_pi/E_pi;
  double gamma=1./sqrt(1-beta*beta);
  double mom_pi_to_p = gamma*ProtonMass*beta;

  //std::cout << "mom_pi = " << mom_pi << ", beta = " << beta << ", gamma = " << gamma                           
  //<< ", mom_pi_to_p = " << mom_pi_to_p << std::endl;                                                           

  bool flagTracking_Pi=false;

  if (theta<=90) {
    flagTracking_Pi = CheckTracking_P(vtz, theta, mom_pi_to_p);
  } else {
    flagTracking_Pi = CheckTracking_P(-vtz, 180-theta, mom_pi_to_p);
  }

  return flagTracking_Pi;

}
