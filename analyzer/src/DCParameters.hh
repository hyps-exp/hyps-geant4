/*
  DCParameters.hh

  2012/1/24
*/

#ifndef DCParameters_h 
#define DCParameters_h 1 

struct DCPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

extern const DCPairPlaneInfo PPInfoBcOut[], PPInfoSdcIn[], PPInfoSdcOut[];
extern const int NPPInfoBcOut, NPPInfoSdcIn, NPPInfoSdcOut;

const DCPairPlaneInfo PPInfoSdcInHyps[] = {
  { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, // SDC0
  { false,  5,  0,  12.0 },  { false,  6,  0,  12.0 }, // SDC1
  { true,  7,  8,  12.0 }, { true,  9,  10,  12.0 } // SDC1
};

const DCPairPlaneInfo PPInfoSdcOutHyps[] = {
  { false,  1,  0,  20.0 }, { true,  2,  3,  20.0 }, { true,  4,  5,  20.0 }, // SDC2
  { false,  6,  0,  20.0 }, { true,  7,  8,  20.0 }, { true,  9,  10,  20.0 } // SDC3
};

const int NPPInfoSdcInHyps  = sizeof(PPInfoSdcInHyps) /sizeof(DCPairPlaneInfo);
const int NPPInfoSdcOutHyps  = sizeof(PPInfoSdcOutHyps) /sizeof(DCPairPlaneInfo);

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBcOut[] = {
  { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, { true,  5,  6,  3.0 },
  { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }, { true, 11, 12,  5.0 },
  { false, 13, 13,  5.0 }, { false, 14, 14,  5.0 },{ false, 15, 15,  5.0 },{ false, 16, 16,  5.0 }
};


const DCPairPlaneInfo PPInfoSdcIn[] = {
  { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, { true,  5,  6,  6.5 },
  { true,  7,  8,  6.5 }, { false,  9,  0,  6.5 }, { false,  10, 0,  6.5 }
};
/*
const DCPairPlaneInfo PPInfoSdcIn[] = {
  { false,  1,  1,  3.0 }, { false,  2,  2,  3.0 }, {false,  3,  3,  3.0 },{false,  4,  4,  3.0 },
  { false,  5,  5,  3.0 }, { false,  6,  6,  3.0 }, {false,  7,  7,  3.0 },
  { true,  8,  9,  5.0 }, { true,  10,  11,  5.0 }, { true,  12, 13,  5.0 }
};
*/
const DCPairPlaneInfo PPInfoSdcOut[] = {
  { false,  1,  1,  3.0 }, { false,  2,  2,  3.0 }, {false,  3,  3,  3.0 }, {false,  4,  4,  3.0 },
  { false,  5,  0,  3.0 }, { false,  6,  0,  3.0 }, {false,  7,  0,  3.0 }, 
  { false,  8,  0,  3.0 }, { false,  9,  0,  3.0 }, {false,  10,  0,  3.0 }, 
  { false,  11,  0,  3.0 }, { false,  12,  0,  3.0 }, {false,  13,  0,  3.0 }, 
  { false,  14,  0,  3.0 }, { false,  15,  0,  3.0 }, {false,  16,  0,  3.0 }, 
  { false,  17,  17,  3.0 }, { false,  18,  18,  3.0 }, {false,  19,  19,  3.0 }, {false,  20,  20,  3.0 },
};

const int NPPInfoBcOut = sizeof(PPInfoBcOut)/sizeof(DCPairPlaneInfo);
const int NPPInfoSdcIn  = sizeof(PPInfoSdcIn) /sizeof(DCPairPlaneInfo);
const int NPPInfoSdcOut  = sizeof(PPInfoSdcOut) /sizeof(DCPairPlaneInfo);
#endif

const int MinNumOfHitsBcIn   = 8;
const int MinNumOfHitsBcOut  = 8;
const int MinNumOfHitsSdcIn  = 10; // for HYPS
const int MinNumOfHitsSdcOut = 10; // for HYPS

//DL Ranges (BC1&2 for Time range -5 ns <[Time gate]<75 ns)
const double MinDLBc[29] = {
   0.0,
   // BC1
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC2
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
   // BC4
   -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
   // VFT1
   -0.5, -0.5, -0.5, -0.5
};

const double MaxDLBc[29] = {
  0.0,
  // BC1
  75.0, 75.0, 75.0, 75.0, 75.0, 75.0,
  //  100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
  // BC2
  75.0, 75.0, 75.0, 75.0, 75.0, 75.0,
  //  100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
  // BC3
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8,
  // BC4
  3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
  // VFT1
  3.0, 3.0, 3.0, 3.0
};

const double MinDLSdc[47] = {
  0.0,
  // SDC1
  -0.5, -0.5, -0.5, -0.5,
  // SDC2
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // Dummy Id=11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SH1
  // -0.5, -0.5, -0.5, -0.5,
  // SDC3
  -0.5, -0.5, -0.5, -0.5, -0.5,
  // SDC4
  -0.5, -0.5, -0.5, -0.5, -0.5,
  // SH2
  // -0.5, -0.5, -0.5, -0.5
}; 

const double MaxDLSdc[47] = {
  0.0,
  // SFT
  // 3.0, 3.0, 3.0, 3.0,
  // SDC0 for Hyps
  3.0, 3.0, 3.0, 3.0,  
  // SDC1 for Hyps
  6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 
  // Dummy Id=10-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SH1
  // 10.0, 10.0, 10.0, 10.0,
  // SDC3
  10.5, 10.5, 10.5, 10.5, 10.5, 
  // SDC4
  10.5, 10.5, 10.5, 10.5, 10.5, 
  // SH2
  // 10.0, 10.0, 10.0, 10.0
};

/*
const double MinDLSdc[47] = {
  0.0,
  // VFT2
  -5.0, -5.0, -5.0, -5.0
  // SFT
  -5.0, -5.0, -5.0,
  // SDC2
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // Dummy Id=10-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SH1
  -0.5, -0.5, -0.5, -0.5,
  // SDC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // SDC4
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // SH2
  -0.5, -0.5, -0.5, -0.5
}; 

const double MaxDLSdc[47] = {
  0.0,
  // VFT2
  75.0, 75.0, 75.0, 75.0,
  // SFT
  75.0, 75.0, 75.0,
  // SDC2
  6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 
  // Dummy Id=10-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  // SH1
  10.0, 10.0, 10.0, 10.0,
  // SDC 3
  13.0, 13.0, 13.0, 13.0, 13.0, 13.0,
  // SDC4
  13.0, 13.0, 13.0, 13.0, 13.0, 13.0,
  // SH2
  10.0, 10.0, 10.0, 10.0
};
*/
#endif
