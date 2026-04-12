/*
  DetectorID.hh

  2012/1/24
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//Detector ID


/* 0.75 mm setup */

// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 8;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 584;
const int NumOfSegCFT_PHI2  = 692;
const int NumOfSegCFT_PHI3  = 800;
const int NumOfSegCFT_PHI4  = 910;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 584, 472, 692, 510, 800, 538, 910};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80., 84.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75};
//const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
const double CFT_Z = 400.;


/* phi 1 mm setup (spacing 0.6)*/
/*
// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 7;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 564;
const int NumOfSegCFT_PHI2  = 670;
const int NumOfSegCFT_PHI3  = 774;
const int NumOfSegCFT_PHI4  = 878;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 564, 472, 670, 510, 774, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 1.0, 0.75, 1.0, 0.75, 1.00, 0.75};
const double CFT_Z = 400.;
*/

/* phi 1 mm setup (spacing 0.65)*/
/*
// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 7;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 522;
const int NumOfSegCFT_PHI2  = 618;
const int NumOfSegCFT_PHI3  = 714;
const int NumOfSegCFT_PHI4  = 812;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 522, 472, 618, 510, 714, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 1.0, 0.75, 1.0, 0.75, 1.00, 0.75};
const double CFT_Z = 400.;
*/

/* phi 1 mm setup (spacing 0.68)*/
/*
// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 7;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 498;
const int NumOfSegCFT_PHI2  = 590;
const int NumOfSegCFT_PHI3  = 682;
const int NumOfSegCFT_PHI4  = 776;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 498, 472, 590, 510, 682, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 1.0, 0.75, 1.0, 0.75, 1.00, 0.75};
const double CFT_Z = 400.;
*/

/* phi layer1: 0.75 mm(0.55 spacing ), layer2,3: 1 mm setup (spacing 0.65)*/
/*
// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 7;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 616;
const int NumOfSegCFT_PHI2  = 618;
const int NumOfSegCFT_PHI3  = 714;
const int NumOfSegCFT_PHI4  = 812;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 616, 472, 618, 510, 714, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 0.75, 0.75, 1.0, 0.75, 1.00, 0.75};
const double CFT_Z = 400.;
*/

/* phi layer 1,2: 1 mm setup (spacing 0.65), layer3: 0.75 (spacing 0.55) */
/*
// CFT config
enum CFT_PLANE{CFT_U1, CFT_PHI1, CFT_V2, CFT_PHI2, CFT_U3, CFT_PHI3, CFT_V4, CFT_PHI4};
const int NumOfPlaneCFT  = 7;
const int NumOfSegCFT_U1  = 426;
const int NumOfSegCFT_V2  = 472;
const int NumOfSegCFT_U3  = 510;
const int NumOfSegCFT_V4  = 538;
const int NumOfSegCFT_PHI1  = 522;
const int NumOfSegCFT_PHI2  = 618;
const int NumOfSegCFT_PHI3  = 844;
const int NumOfSegCFT_PHI4  = 812;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 522, 472, 618, 510, 844, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 1.0, 0.75, 1.0, 0.75, 0.75, 0.75};
const double CFT_Z = 400.;

*/

/*
const int    NumOfBGOUnit = 8;
const int    NumOfBGOInOneUnit = 3;
const double BGO_X = 30.;
const double BGO_Y = 25.;
const double BGO_Z = 400.;
const double RadiusOfBGOSurface = 110.;
*/

const double BGO_X = 30.;
const double BGO_Y = 25.;
const double BGO_Z = 400.;
const int    NumOfBGOUnit = 8;
const int    NumOfBGOInOneUnit = 2;
const double RadiusOfBGOSurface = 100.;
const int    NumOfBGOInOneUnit2 = 1;
const double RadiusOfBGOSurface2 = 120.;

const double BGO_PMT_X = 30.;
const double BGO_PMT_Y = 30.;
const double BGO_PMT_Z = 32.5;

const double BGO_MagShld_R = 89.0;
const double BGO_MagShld_Thick = 10.0;
const double BGO_MagShld_Z = 20.;

const int    NumOfPiVUnit = 8;
const int    NumOfPiVInOneUnit = 3;
const double PiV_X = 30.;
const double PiV_Y = 10.;
const double PiV_Z = 400.;
const double RadiusOfPiVSurface = 140.;

const double PiV2_X = 40.;
const double PiV2_Y = 10.;
const double PiV2_Z = 400.;
const double RadiusOfPiV2Surface = 160.;


// LH2 target
const double Target_Length = 300.;
const G4double LH2TgtR = 40.; //tyokkei
//const G4double LH2TgtR = 70.; //tyokkei
const G4double TargetVessThickness = 0.25;

const G4double VaccumChamWinR = 80.;
const G4double VaccumChamThickness = 1.;
const G4double VaccumChamWinZ = 500.;

// LH2 Holder part
// G10
const G4double TgtG10HolderZ = 20.;
const G4double TgtG10HolderThickness = 2.;
// SUS
const G4double TgtSUSHolderZ = 40.;
const G4double TgtSUSHolderThickness = 0.5;

// BC3
const int NumOfPlaneBC3  = 6;

// BC4
const int NumOfPlaneBC4  = 6;

// SDC2
const int NumOfPlaneSDC2  = 6;

// SDC3
const int NumOfPlaneSDC3  = 6;

// SDC4
const int NumOfPlaneSDC4  = 6;

// CH
//const int NumOfSegmentCH = 24;
const int NumOfSegmentCH = 64;

// SH1
const int NumOfPlaneSH1  = 4;
const int NumOfSegmentSH1 = 48;

// SH2
const int NumOfPlaneSH2  = 4;
const int NumOfSegmentSH2 = 48;

// FTOF
//const int NumOfSegmentFTOF = 24;
const int NumOfSegmentFTOF = 40;

// AC
//const int NumOfSegmentAC = 10;
const int NumOfSegmentAC = 1;

//Counter
const int NumOfCounter = 5; // FTOF, CH, PIV, AC, BGO_VP
//const int CounterSegmentMax = NumOfSegmentFTOF;
const int CounterSegmentMax = 32; // PiV segment number

#endif
