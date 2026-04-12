/*
  DetectorID.hh

  2012/1/24
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//Detector ID
//Base chambers
const int DetIdBC1  =101;
const int DetIdBC2  =102;
const int DetIdBC3  =103;
const int DetIdBC4  =104;
const int DetIdSDC1 =105;
const int DetIdSDC2 =106;
const int DetIdSDC3 =107;
const int DetIdSDC4 =108;
const int DetIdK6BDC=109;
const int DetIdBFT  =110;
const int DetIdSFT  =111;
const int DetIdSsd  =112;
const int DetIdCFT  =113;
const int DetIdBGO  =114;
const int DetIdT54Counter  =115;
const int DetIdPiV  =116;

//Base counters
const int DetIdBH1 =  3;
const int DetIdBH2 =  4;
const int DetIdTOF =  5;
const int DetIdAC  =  6;
const int DetIdLC  =  7;
const int DetIdGC  =  8;

//Test counters
const int DetIdBAC =  9;
const int DetIdTGT =  10;

//Misc
const int DetIdMisc  =  11;
const int DetIdMatrix=  12;

// CH
const int DetIdCH =  13;

//For SksMinus
const int DetIdSP0 =  20;

//For E27
const int DetIdRC  = 21;
const int DetIdAC1  = 22;
const int DetIdBVH  = 23;

//Scaler
const int DetIdScaler = 210;


//DC Number of Plane
const int PlMinBcIn  =  1;
const int PlMaxBcIn  = 12;
const int PlMinBcOut = 13;
//const int PlMaxBcOut = 24;
const int PlMaxBcOut = 28;

const int PlMinSdcIn  =  1;
//const int PlMaxSdcIn  =  9;
const int PlMaxSdcIn  =  10; // for HYPS

//const int PlMinSdcOut = 31;
//const int PlMaxSdcOut = 42;
const int PlMinSdcOut = 31; // for HYPS
const int PlMaxSdcOut = 40; // for HYPS

const int PlOffsBc = 100;
const int PlOffsSdcOut = 30;
const int PlOffsSsd = 130;

const int NumOfLayersBc     = 6;
const int NumOfLayersSdc    = 6;
const int NumOfLayersSFT    = 3;
const int NumOfLayersSdc2   = 6;
const int NumOfLayersSsd    = 2;
const int NumOfLayersSH1    = 4;
const int NumOfLayersSH2    = 4;
/*phi 0.75 setup*/
const int NumOfLayersCFT    = 8;
/*phi 1.00 setup*/
//const int NumOfLayersCFT    = 7;

const int NumOfLayersT54Counter   = 5;

const int NumOfLayersBcIn  = PlMaxBcIn  - PlMinBcIn  + 1;
const int NumOfLayersBcOut = PlMaxBcOut - PlMinBcOut + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;

//BFT config
const int NumOfPlaneBFT  = 2;
const int NumOfSegBFT    = 160;

// SFT config
// SFT X layer has U D plane.
// SFT UV layers have only U plnane.
enum SFT_PLANE{SFT_X1, SFT_X2, SFT_V, SFT_U};
const int NumOfPlaneSFT  = 4;
const int NumOfSegSFT_X  = 256;
const int NumOfSegSFT_UV = 320;

// Target
// Cylindrical target
const double CyLH2TgtR = 40.; //mm //tyokkei
const double CyLH2TgtZ = 300.; // mm
const double TargetVessThickness = 0.25; // mm

const double VaccumChamWinR = 80.;
const double VaccumChamThickness = 1.;
const double VaccumChamWinZ = 500.;


// CFT config
/*phi 0.75 mm setup*/

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
//const int NumOfSegCFT_PHI4  = 910;

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
//const int NumOfSegCFT_PHI4  = 812;

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
//const int NumOfSegCFT_PHI4  = 812;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 498, 472, 590, 510, 682, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 1.0, 0.75, 1.0, 0.75, 1.00, 0.75};
const double CFT_Z = 400.;
*/

/* phi layer1 : 0.75mm (spacing 0.55), layer2,3 : 1 mm setup (spacing 0.65)*/
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
//const int NumOfSegCFT_PHI4  = 812;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 616, 472, 618, 510, 714, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 0.75, 0.75, 1.0, 0.75, 1.00, 0.75};
const double CFT_Z = 400.;
*/

/* phi layer1, 2 : 1 mm setup (spacing 0.65), layer3: 0.75mm (spacing 0.55)*/
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
//const int NumOfSegCFT_PHI4  = 812;

const int NumOfSegCFT[NumOfPlaneCFT]  = {426, 522, 472, 618, 510, 844, 538};
const double RadiusOfCFT[NumOfPlaneCFT] = {50., 54., 60., 64, 70., 74., 80.};
const double DiameterOfCFTFiber[NumOfPlaneCFT] = {0.75, 1.0, 0.75, 1.0, 0.75, 0.75, 0.75};
const double CFT_Z = 400.;
*/



// BGO 
const double BGO_X = 30.;
const double BGO_Y = 25.;
const double BGO_Z = 400.;
const int    NumOfBGOUnit = 8;
//const int    NumOfBGOInOneUnit = 3;
//const double RadiusOfBGOSurface = 110.;
const int    NumOfBGOInOneUnit = 2;
const double RadiusOfBGOSurface = 100.;
const int    NumOfBGOInOneUnit2 = 1;
const double RadiusOfBGOSurface2 = 120.;
//const int NumOfSegBGO = NumOfBGOUnit*NumOfBGOInOneUnit;
const int NumOfSegBGO = NumOfBGOUnit*(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);


//PiV
const int NumOfSegPiV = 32;
const int NumOfPiVUnit = 8;
const int NumOfPiVInOneUnit = 3;
const int    NumOfPiVLayer = 40;
const double PiV_X0 = 30.;
const double PiV_Y = 15.;
const double PiV_Z = 400.;
//const double RadiusOfPiVSurface = 140.;
const double RadiusOfPiVSurface = 160.;

//const double PiV2_X = 40.;
const double PiV2_X0 = 38.;
const double PiV2_Y = 15.;
const double PiV2_Z = 400.;
//const double RadiusOfPiV2Surface = 160.;
const double RadiusOfPiV2Surface = 165.;


// T54 counter
const int NumOfPlaneT54Counter  = NumOfLayersT54Counter;

//SP0 Number of Plane
const int PlMinSP0 =  1;
const int PlMaxSP0 =  8;
const int NumOfLayersSP0 = PlMaxSP0 - PlMinSP0 + 1;

const int IdK18Target = 130;
const int IdBH1       = 201;
const int IdBFT       = 203;

//Hodo Segments
const int NumOfSegGC  =   1;
const int NumOfSegBH1 =  11;
const int NumOfSegBH2 =   8;
const int NumOfSegTOF =  32;
const int NumOfSegAC  =  30;
const int NumOfSegAC1  =  18;
const int NumOfSegBVH  =  4;
const int NumOfSegLC  =  28;
//const int NumOfSegCH =  24;
const int NumOfSegCH =  64;

const int NumOfSegBAC =  10;
const int NumOfSegTGT =   3;

const int NumOfMisc    =  20;
const int NumOfMatrix  =  10;

const int NumOfSegSP0  =  5;

//for E27
const int NumOfSegRC = 42;
const int NumOfLayersRC = 7; 
const int PlaneId_RC = 6;  

const int NumOfScaler = 32;



//Number of Wires
const int MaxWireBC1  =  256;
const int MaxWireBC2  =  256;
const int MaxWireBC3  =  64;
//const int MaxWireBC4  =  48;
const int MaxWireBC4  =  64;

const int MaxWireSsd  = 768;

const int MaxWireSDC1  =  64;
const int MaxWireSDC2  =  96;

const int MaxWireSDC3   =  120;
const int MaxWireSDC3X  =  108;
const int MaxWireSDC3U  =  120;
const int MaxWireSDC3V  =  120;

const int MaxWireSDC4   =  120;
const int MaxWireSDC4X  =  108;
const int MaxWireSDC4U  =  120;
const int MaxWireSDC4V  =  120;



#endif
