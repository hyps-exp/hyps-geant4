/*
  CFTTrackingEffMan.hh

  2012/1/24
*/

#ifndef CFTTrackingEffMan_h
#define CFTTrackingEffMan_h 1

#include "ThreeVector.hh"
#include "DetectorID.hh"
#include <string>
#include <vector>
#include <map>

const int NdivVertex=10;
const int NdivThetaEff=90;


class CFTTrackingEffMan 
{
private:
  CFTTrackingEffMan();
public:
  ~CFTTrackingEffMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; bool ret=Initialize(); return ret;}
  bool Initialize( const std::string &filename )
  { filename_=filename; bool ret=Initialize(); return ret;}

  static CFTTrackingEffMan & GetInstance( void );

  bool CheckTracking_P(double vtz, double theta, double mom);
  bool CheckTracking_Pi(double vtz, double theta, double mom);
private:
  static CFTTrackingEffMan *CftEffMan_;
  std::string filename_;

  double p0[NdivVertex][NdivThetaEff];
  double p1[NdivVertex][NdivThetaEff];
  double p2[NdivVertex][NdivThetaEff];
  double p0_simu[NdivVertex][NdivThetaEff];
  double p1_simu[NdivVertex][NdivThetaEff];
  double p2_simu[NdivVertex][NdivThetaEff];

};


#endif
