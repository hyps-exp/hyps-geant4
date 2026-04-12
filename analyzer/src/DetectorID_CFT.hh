/*
  DetectorID.hh

  2012/1/24
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//Detector ID
const int DetIdCFT  =113;
const int DetIdBGO  =114;
const int DetIdT54Counter  =115;

// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 8;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 534;
const int NumOfSegCFT_PHI2  = 692;
const int NumOfSegCFT_PHI3  = 800;
const int NumOfSegCFT_PHI4  = 910;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 534, 472, 692, 510, 800, 538, 910};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80., 84.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75};
const double CFT_Z = 400.;

const int    NumOfBGOUnit = 8;
const int    NumOfBGOInOneUnit = 3;
const double BGO_X = 30.;
const double BGO_Y = 25.;
const double BGO_Z = 400.;
const double RadiusOfBGOSurface = 110.;

const double Target_Length = 300.;


#endif
