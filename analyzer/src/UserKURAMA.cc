/*
  SksTracking.cc
*/

#include "ThreeVector.hh"
#include "LorentzVector.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "HodoAnalyzer.hh"
#include "CFTFiberHit.hh"
#include "CFTFiberCluster.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "DCAnalyzer.hh"
#include "DCLocalTrack.hh"
#include "SksTrack.hh"
#include "TrackHit.hh"
#include "Kinematics.hh"
#include "EnergyCorrection.hh"
#include "SksParticle.hh"
#include "TemplateLib.hh"

#include "HistHelper.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FieldMan.hh"
#include "MatrixTrigMan.hh"

#include <TRandom.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <signal.h>


//#define Interactive 1

////////////////Physics Parameters///////////////////
const double AtomicMassUnit  = 0.93149432;
const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double ProtonMass      = 0.93827200;
const double NeutronMass     = 0.93956563;
const double LambdaMass      = 1.115648;
const double SigmaMinusMass  = 1.197449;
////////////////////////////////////////////////////

const double MinOfMassSquareK = 0.1;
const double MaxOfMassSquareK = 0.4;

bool ProcessOneEvent( std::ifstream & );
void DefineHistograms( const char * );
void InitializeEvent( void );
double calcKMomFromTheta(double bmom, double cost, double p0);

const double MaxChiSqr = 100.;
const double PIni = 1.0;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double EdepThrTOF = 1.; //MeV
const double EdepThrCH  = 0.1; //MeV

#ifndef MaxHits 
#define MaxHits 30
#endif
#ifndef MaxHits2 
#define MaxHits2 60
#endif

struct Event{
  //Primary for Yd scattering
  double thetaMeson;
  double phiMeson;
  double thetaMesonCM;
  double phiMesonCM;
  double thetaScatHypCM;
  double phiScatHypCM;
  double thetaScatPLab;
  double thetaScatHypLab;
  double beammom;
  double momVectorScatMeson[3];
  double momVectorHypBeam[3];
  double momVectorHypScat[3];
  double momVectorProtonScat[3];
  double momVectorDecayPi[3];
  double momVectorDecayNucleon[3];
  double momScatMeson;
  double momHypBeam;
  double momHypScat;
  double momProtonScat;
  double momDecayPi;
  double momDecayNucleon;
  double primaryVertex[3];
  double scatPos0[3];
  double NNscatPos[3];
  double PiNscatPos[3];
  double decayPos[3];
  int    decayFlag;
  int    scatFlag;
  int    scatTarget;
  int    NNscatFlag;
  int    NNscatTarget;
  int    PiNscatFlag;
  int    PiNscatTarget;


  //Fiber1-4
  int    FiberHits[NumOfPlaneCFT];
  double    FiberSeg [NumOfPlaneCFT][MaxHits2];
  double FiberTime[NumOfPlaneCFT][MaxHits2];
  double FiberEdep[NumOfPlaneCFT][MaxHits2];
  int    FiberPID [NumOfPlaneCFT][MaxHits2];

  double fClEdepForScat[NumOfPlaneCFT];
  double fClEdepNormForScat[NumOfPlaneCFT];
  double fClPathLengthForScat[NumOfPlaneCFT];

  double fClEdepForPi[NumOfPlaneCFT];
  double fClEdepNormForPi[NumOfPlaneCFT];
  double fClPathLengthForPi[NumOfPlaneCFT];

  //Crystal
  int    CrystalHits;
  double    CrystalSeg [MaxHits2];
  double CrystalTime[MaxHits2];
  double CrystalEdep[MaxHits2];
  int    CrystalPID [MaxHits2];

  double cClEdepForScat;
  double cClEdepForPi;
  double cRcEdepForScat;
  double cRcEdepForPi;

  // TOF
  int nhTof;
  int    TofSeg[MaxHits];
  double TofTime[MaxHits];
  double TofEdep[MaxHits];

  // LC
  int nhLc;
  int    LcSeg[MaxHits];
  double LcTime[MaxHits];
  double LcEdep[MaxHits];

  // CH
  int nhCh;
  int    ChSeg[MaxHits];
  double ChTime[MaxHits];
  double ChEdep[MaxHits];

  // AC
  int nhAc;
  int    AcSeg[MaxHits];
  double AcTime[MaxHits];
  double AcEdep[MaxHits];

  // SFTX
  int nhSftx;
  double SftxSeg[MaxHits];

  // SdcIn
  int ntSdcIn;
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];

  // SdcOut
  int ntSdcOut;
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double xtofSdcOut[MaxHits];
  double ytofSdcOut[MaxHits];

  // Tracking
  int ntSks;
  double p[MaxHits];

  int ntSksPart;
  double pSksPart[MaxHits];
  double m2[MaxHits];

  int nK;
  double pK[MaxHits];
  double pKCor[MaxHits];
  double pKCal[MaxHits];
  double theta[MaxHits];
  double MissMass[MaxHits];
  double MissMassCor[MaxHits];
  double pSigmaCor[MaxHits];
  double pSigmaCal[MaxHits];
  double Vertex_x[MaxHits];
  double Vertex_y[MaxHits];
  double Vertex_z[MaxHits];

  // Matrix trigger
  int flagMatrixTrig;
  int flagSftTrig;
  int flagMassTrig;
};
static Event event;

const int MinHitForPi = 6;
const int MinHitForP = 4;

void closeFile(int sig)
{
  if(gFile)
    {
      gFile->Write();
      gFile->Close();
    }
}

enum eArg
  {
    kArgProcessName,
    kArgConfFile,
    kArgInFile,
    kArgOutRootFile,
    kArgc
  };


TRandom *gRand;
const double time_reso = 0.25;
const double time_reso_m2 = 0.15;

int main( int argc, char **argv )
{
  if( argc != kArgc ) {
    std::cerr << "Usage: " << argv[kArgProcessName] 
	      << " [confFile]  [dataFile] [RootFile]"
	      << std::endl;
    exit(-1);
  }

  gRand = new TRandom();

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  const std::string& inFile   = arg[kArgInFile];
  const std::string& rootFile = arg[kArgOutRootFile];

  ConfMan* gconfManager = new ConfMan(confFile);
  if (!gconfManager->Initialize()) {
    std::cerr << "#E gconfManager->Initialize failed ! in Main()" << std::endl;
  }
  
  DefineHistograms(rootFile.c_str());

  signal(SIGINT,closeFile);
  
  std::ifstream InputData( inFile.c_str() );
  std::cout<<"***************************************************"<<std::endl;
  
  while( ProcessOneEvent( InputData ) ){
#ifdef Interactive 
    int it; 
    std::cout << "# ";
    std::cin >> it;
    if( it<0 ) break;
#endif
  }

  std::cout<<"*************************************"<<std::endl;
  gFile->Write();
  gFile->Close();
  std::cout<<"*************************************"<<std::endl;

  return 0;
}

bool ProcessOneEvent( std::ifstream &In )
{
  const std::string funcname = "ProcessEvent";

  const int maxMultiSdcIn = 3;
  const int maxMultiSdcOut = 3;
  const int maxMultiBcOut = 3;

  static int events=0;
  events++;
  if (events%100==0)
    std::cout<<"Events:"<< events<<std::endl;


  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  if( !In ) return false;

  char line[256];
  // Event header
  In.getline(line, sizeof(line));
  // Primary Info header
  In.getline(line, sizeof(line));
  // primary Info data
  In.getline(line, sizeof(line));
  sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	 &event.thetaMeson,&event.phiMeson,&event.thetaMesonCM,
	 &event.phiMesonCM,&event.thetaScatHypCM,&event.phiScatHypCM, 
	 &event.thetaScatPLab, &event.thetaScatHypLab, &event.beammom);
  In.getline(line, sizeof(line));
  sscanf(line,"momVectorScatMeson %lf %lf %lf",
	 &event.momVectorScatMeson[0], &event.momVectorScatMeson[1], &event.momVectorScatMeson[2]);
  In.getline(line, sizeof(line));
  sscanf(line,"momVectorHypBeam %lf %lf %lf",
	 &event.momVectorHypBeam[0], &event.momVectorHypBeam[1], &event.momVectorHypBeam[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"momVectorHypScat %lf %lf %lf",
	 &event.momVectorHypScat[0], &event.momVectorHypScat[1], &event.momVectorHypScat[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"momVectorProtonScat %lf %lf %lf",
	 &event.momVectorProtonScat[0], &event.momVectorProtonScat[1], &event.momVectorProtonScat[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"momVectorDecayPi %lf %lf %lf",
	 &event.momVectorDecayPi[0], &event.momVectorDecayPi[1], &event.momVectorDecayPi[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"momVectorDecayNucleon %lf %lf %lf",
	 &event.momVectorDecayNucleon[0], &event.momVectorDecayNucleon[1], &event.momVectorDecayNucleon[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"primaryVertex %lf %lf %lf",
	 &event.primaryVertex[0], &event.primaryVertex[1], &event.primaryVertex[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"scatPos %lf %lf %lf",
	 &event.scatPos0[0], &event.scatPos0[1], &event.scatPos0[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"NNscatPos %lf %lf %lf",
	 &event.NNscatPos[0], &event.NNscatPos[1], &event.NNscatPos[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"PiNscatPos %lf %lf %lf",
	 &event.PiNscatPos[0], &event.PiNscatPos[1], &event.PiNscatPos[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"decayPos %lf %lf %lf",
	 &event.decayPos[0], &event.decayPos[1], &event.decayPos[2]);

  In.getline(line, sizeof(line));
  sscanf(line,"Flags %d %d %d %d %d %d %d",
	 &event.decayFlag, &event.scatFlag, &event.scatTarget, &event.NNscatFlag,
	 &event.NNscatTarget, &event.PiNscatFlag, &event.PiNscatTarget);

  event.momScatMeson = 
    sqrt(event.momVectorScatMeson[0]*event.momVectorScatMeson[0] + 
         event.momVectorScatMeson[1]*event.momVectorScatMeson[1] + 
         event.momVectorScatMeson[2]*event.momVectorScatMeson[2]);
  event.momHypBeam = 
    sqrt(event.momVectorHypBeam[0]*event.momVectorHypBeam[0] + 
         event.momVectorHypBeam[1]*event.momVectorHypBeam[1] + 
         event.momVectorHypBeam[2]*event.momVectorHypBeam[2]);
  event.momHypScat = 
    sqrt(event.momVectorHypScat[0]*event.momVectorHypScat[0] + 
         event.momVectorHypScat[1]*event.momVectorHypScat[1] + 
         event.momVectorHypScat[2]*event.momVectorHypScat[2]);
  event.momProtonScat = 
    sqrt(event.momVectorProtonScat[0]*event.momVectorProtonScat[0] + 
         event.momVectorProtonScat[1]*event.momVectorProtonScat[1] + 
         event.momVectorProtonScat[2]*event.momVectorProtonScat[2]);
  event.momDecayPi = 
    sqrt(event.momVectorDecayPi[0]*event.momVectorDecayPi[0] + 
         event.momVectorDecayPi[1]*event.momVectorDecayPi[1] + 
         event.momVectorDecayPi[2]*event.momVectorDecayPi[2]);
  event.momDecayNucleon = 
    sqrt(event.momVectorDecayNucleon[0]*event.momVectorDecayNucleon[0] + 
         event.momVectorDecayNucleon[1]*event.momVectorDecayNucleon[1] + 
         event.momVectorDecayNucleon[2]*event.momVectorDecayNucleon[2]);

#if 0
  std::cout << "thetaMeson = " << event.thetaMeson 
	    << " phiMeson = " << event.phiMeson 
	    << " thetaMesonCM = " << event.thetaMesonCM 
	    << " phiMesonCM = "   << event.phiMesonCM 
	    << " thetaScatHypCM = " << event.thetaScatHypCM 
	    << " phiScatHypCM = " << event.phiScatHypCM << std::endl;
  
  std::cout << "momVectorScatmeson = ( " 
	    << event.momVectorScatMeson[0] << ", " 
	    << event.momVectorScatMeson[1] << ", " 
	    << event.momVectorScatMeson[2] << ")" << std::endl ;
  std::cout << "momVectorHypBeam = ( " 
	    << event.momVectorHypBeam[0] << ",  " 
	    << event.momVectorHypBeam[1] << ", " 
	    << event.momVectorHypBeam[2] << ")"  << std::endl;
  std::cout << "momVectorHypScat = ( " 
	    << event.momVectorHypScat[0] << ", " 
	    << event.momVectorHypScat[1] << ", " 
	    << event.momVectorHypScat[2] << ")"  << std::endl;
  std::cout << "momVectorProtonScat = ( " 
	    << event.momVectorProtonScat[0] << ", " 
	    << event.momVectorProtonScat[1] << ", " 
	    << event.momVectorProtonScat[2] << ")" << std::endl; ;
  std::cout << "momVectorDecayPi = ( " 
	    << event.momVectorDecayPi[0] << ", " 
	    << event.momVectorDecayPi[1] << ", " 
	    << event.momVectorDecayPi[2] << ")" << std::endl ;
  std::cout << "momVectorDecayNucleon = ( " 
	    << event.momVectorDecayNucleon[0] << ", " 
	    << event.momVectorDecayNucleon[1] << ", " 
	    << event.momVectorDecayNucleon[2] << ")"  << std::endl;
  std::cout << "primaryVertex = ( " 
	    << event.primaryVertex[0] << ", " 
	    << event.primaryVertex[1] << ", " 
	    << event.primaryVertex[2] << ")"  << std::endl;
  std::cout << "scatPos = ( " 
	    << event.scatPos0[0] << ", " 
	    << event.scatPos0[1] << ", " 
	    << event.scatPos0[2] << ")"  << std::endl;
  std::cout << "NNscatPos = ( " 
	    << event.NNscatPos[0] << ", " 
	    << event.NNscatPos[1] << ", " 
	    << event.NNscatPos[2] << ")"  << std::endl;
  std::cout << "decayPos = ( " 
	    << event.decayPos[0] << ", " 
	    << event.decayPos[1] << ", " 
	    << event.decayPos[2] << ")"  << std::endl;
  
  std::cout << "decayFlag = " << event.decayFlag 
	    << ", scatFlag =  " << event.scatFlag 
	    << ", scatTarget =  " << event.scatTarget 
	    << ", NNscatFlag =  " << event.NNscatFlag 
	    << ", NNscatTarget = " << event.NNscatTarget << std::endl;
#endif

  /*  
  double gx,gy,gz,gp,gt,gf,gpb,gub,gvb;
  In >> gx >> gy >> gz >> gp >> gt >> gf >> gpb >> gub >> gvb;
  */
#if 0
  std::cout<< "*********************************************" <<std::endl;
  std::cout<<" "<< gx 
	   <<" "<< gy 
	   <<" "<< gz 
	   <<" "<< gp 
	   <<" "<< gt 
	   <<" "<< gf
	   <<" "<< gpb
	   <<" "<< gub
	   <<" "<< gvb
	   <<std::endl;
  std::cout<< "*********************************************" <<std::endl;
#endif

  /*
  DCAnalyzer DCAna;
  CylindricalAna CylAna;
  CylindricalAna CylAnaPi;
  */

  RawData *rawData = new RawData;

  // Fiber
  In.getline(line, sizeof(line));

  int nhCFT=0;
  sscanf(line, " --- CFT Fiber --- Total # of Hits : %d", &nhCFT);

  for (int i=0; i<nhCFT; i++) {
    int layer, segment, pid;
    double edep, time;
    In.getline(line, sizeof(line));
    sscanf(line,"%d %d %lf [ns] %lf [MeV]  %d",
	   &layer,&segment,&time,&edep, &pid);
    
    if( layer >= 61 && layer <= 68) {
      int flayer = layer-61;
      rawData->AddCFTHodoRawHit(DetIdCFT, flayer, segment, 0, 0, edep);
      rawData->AddCFTHodoRawHit(DetIdCFT, flayer, segment, 0, 1, edep);
      rawData->AddCFTHodoRawHit(DetIdCFT, flayer, segment, 1, 0, time);
      rawData->AddCFTHodoRawHit(DetIdCFT, flayer, segment, 1, 1, time);
      
    }
  }
  // BGO
  In.getline(line, sizeof(line));
  int nhCrystal=0;

  sscanf(line, " --- Crystal --- Total # of Hits : %d", &nhCrystal);

  for (int i=0; i<nhCrystal; i++) {
    int layer, segment, pid;
    double edep, time;
    In.getline(line, sizeof(line));
    sscanf(line,"%d %d %lf [ns] %lf [MeV]  %d",
	   &layer,&segment,&time,&edep, &pid);

    // Crystal
    if( layer >= 81 && layer <= 81){
      int clayer = layer-81;

      rawData->AddBGOHodoRawHit(DetIdBGO, clayer, segment, 0, 0, edep);
      rawData->AddBGOHodoRawHit(DetIdBGO, clayer, segment, 1, 0, time);
      rawData->AddBGOHodoRawHit(DetIdBGO, clayer, segment, 0, 1, edep);
      rawData->AddBGOHodoRawHit(DetIdBGO, clayer, segment, 1, 1, time);

      if (event.CrystalHits < MaxHits2) {
	event.CrystalSeg [event.CrystalHits] = segment;
	event.CrystalTime[event.CrystalHits] = time;
	event.CrystalEdep[event.CrystalHits] = edep;
	event.CrystalPID [event.CrystalHits] = pid;
	event.CrystalHits++;

      } else {
	std::cerr << "Too many crystal hits : " << event.CrystalHits
		  << std::endl;
      }
    }
  }

  // DC
  In.getline(line, sizeof(line));
  int nhDC=0;

  sscanf(line, " --- DC --- Total # of Hits : %d", &nhDC);

  for (int i=0; i<nhDC; i++) {
    int layer, segment;
    double dtime;
    In.getline(line, sizeof(line));
    sscanf(line,"%d  %d  %lf", &layer,&segment,&dtime);
    if (layer >= PlMinSdcIn && layer <= PlMaxSdcIn)
      rawData->AddSdcInRawHit(layer, segment, dtime);
    else if (layer >= PlMinSdcOut && layer <= PlMaxSdcOut)
      rawData->AddSdcOutRawHit(layer, segment, dtime);
    else if (layer >= PlOffsBc+PlMinBcOut && layer <= PlOffsBc+PlMaxBcOut) 
      rawData->AddBcOutRawHit(layer, segment, dtime);
    else 
      std::cout << "No such layer in DC " << layer << std::endl;
  }

  // Counter
  In.getline(line, sizeof(line));
  int nhCounter=0;

  sscanf(line, " --- Counter --- Total # of Hits : %d", &nhCounter);

  for (int i=0; i<nhCounter; i++) {
    int layer, segment;
    double time, edep;
    In.getline(line, sizeof(line));
    sscanf(line,"%d  %d  %lf %lf", &layer,&segment, &time, &edep);

    if (layer == 51) {
      rawData->AddTOFHodoRawHit(DetIdTOF, layer, segment, 0, 0, edep);
      rawData->AddTOFHodoRawHit(DetIdTOF, layer, segment, 1, 0, time);
      rawData->AddTOFHodoRawHit(DetIdTOF, layer, segment, 0, 1, edep);
      rawData->AddTOFHodoRawHit(DetIdTOF, layer, segment, 1, 1, time);
    } else if  (layer == 54) {
      rawData->AddLCHodoRawHit(DetIdLC, layer, segment, 0, 0, edep);
      rawData->AddLCHodoRawHit(DetIdLC, layer, segment, 1, 0, time);
      rawData->AddLCHodoRawHit(DetIdLC, layer, segment, 0, 1, edep);
      rawData->AddLCHodoRawHit(DetIdLC, layer, segment, 1, 1, time);
    } else if  (layer == 10) {
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 0, 0, edep);
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 1, 0, time);
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 0, 1, edep);
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 1, 1, time);
    } else if  (layer == 57) {
      rawData->AddACHodoRawHit(DetIdAC, layer, segment, 0, 0, edep);
      rawData->AddACHodoRawHit(DetIdAC, layer, segment, 1, 0, time);
      //rawData->AddACHodoRawHit(DetIdCH, layer, segment, 0, 1, edep);
      //rawData->AddACHodoRawHit(DetIdCH, layer, segment, 1, 1, time);
    }
  }

  HodoAnalyzer *hodoAna = new HodoAnalyzer;
  hodoAna->DecodeCFTHits(rawData);
  hodoAna->DecodeBGOHits(rawData);
  hodoAna->DecodeCHHits(rawData);
  hodoAna->DecodeTOFHits(rawData);
  hodoAna->DecodeLCHits(rawData);
  hodoAna->DecodeACHits(rawData);

  int nhTof = hodoAna->GetNHitsTOF();
  event.nhTof = nhTof;

  HF1(1001, nhTof);
  for(int i = 0; i<nhTof; ++i){
    const Hodo2Hit* hit = hodoAna->GetHitTOF(i);
    int seg = hit->SegmentId();
    double ctime = hit->CMeanTime();
    double dE = hit->DeltaE();

    HF1(1002, seg);
    HF1(1003, ctime);
    HF1(1004, dE);

    if (i<MaxHits) {
      event.TofSeg[i] = seg;
      event.TofTime[i] = ctime;
      event.TofEdep[i] = dE;
    }
  }

  int nhLc = hodoAna->GetNHitsLC();
  event.nhLc = nhLc;

  //HF1(1001, nhTof);
  for(int i = 0; i<nhLc; ++i){
    const Hodo2Hit* hit = hodoAna->GetHitLC(i);
    int seg = hit->SegmentId();
    double ctime = hit->CMeanTime();
    double dE = hit->DeltaE();

    //HF1(1002, seg);
    //HF1(1003, ctime);
    //HF1(1004, dE);

    if (i<MaxHits) {
      event.LcSeg[i] = seg;
      event.LcTime[i] = ctime;
      event.LcEdep[i] = dE;
    }
  }

  int nhCh = hodoAna->GetNHitsCH();
  HF1(1011, nhCh);
  event.nhCh = nhCh;

  for(int i = 0; i<nhCh; ++i){
    const Hodo2Hit* hit = hodoAna->GetHitCH(i);
    int seg = hit->SegmentId();
    double ctime = hit->CMeanTime();
    double dE = hit->DeltaE();

    HF1(1012, seg);
    HF1(1013, ctime);
    HF1(1014, dE);

    if (i<MaxHits) {
      event.ChSeg[i] = seg;
      event.ChTime[i] = ctime;
      event.ChEdep[i] = dE;
    }
  }

  int nhAc = hodoAna->GetNHitsAC();
  event.nhAc = nhAc;

  for(int i = 0; i<nhAc; ++i){
    const Hodo1Hit* hit = hodoAna->GetHitAC(i);
    int seg = hit->SegmentId();
    double ctime = hit->GetCT();
    double dE = hit->DeltaE();

    if (i<MaxHits) {
      event.AcSeg[i] = seg;
      event.AcTime[i] = ctime;
      event.AcEdep[i] = dE;
    }
  }


  int ncTof=hodoAna->GetNClustersTOF();

  DCAnalyzer *DCAna = new DCAnalyzer;
  DCAna->DecodeSdcInHits(rawData);
  DCAna->DecodeSdcOutHits(rawData);
  DCAna->DecodeBcOutHits(rawData);

  DCHitContainer SFTX_cont = DCAna->GetSdcInHC(3);
  int nhSftx = SFTX_cont.size();
  event.nhSftx = nhSftx;
  for (int ih=0; ih<nhSftx; ih++) {
    DCHit *hit = SFTX_cont[ih];
    double segment = hit->GetWire();

    if (ih<MaxHits)
      event.SftxSeg[ih] = segment;
  }

  /* Trigger selection */
  bool flagMatrixTrigger = false;
  bool flagSftTrigger = false;
  bool flagMassTrigger = false;


  double beta = event.beammom/sqrt(event.beammom*event.beammom+0.139*0.139);
  double time_offset = (210.+event.primaryVertex[2])/(300.*beta);
  double time_offset_tgt0 = (210.)/(300.*beta);

  for(int i = 0; i<nhCh; ++i){
    if (event.ChEdep[i] > EdepThrCH) {
      for(int j = 0; j<nhTof; ++j){
	if (event.TofEdep[j] > EdepThrTOF) {
	  if (MatrixTrigMan::GetInstance().MatrixTrigger(event.ChSeg[i], event.TofSeg[j])) {
	    flagMatrixTrigger = true;

	    //double time = event.TofTime[j] + gRand->Gaus(0, time_reso) + time_offset;
	    //event.TofTime[j] = time;
	    event.TofTime[j] += gRand->Gaus(0, time_reso);
	    double time = event.TofTime[j];
	    if (MatrixTrigMan::GetInstance().MassTrigger(event.ChSeg[i], event.TofSeg[j], time))
	      flagMassTrigger = true;


	    for (int k=0; k<nhSftx; k++) {
	      if (MatrixTrigMan::GetInstance().SftTrigger(event.ChSeg[i], event.TofSeg[j], event.SftxSeg[k])) {
		flagSftTrigger = true;
	      }
	    }
	  }
	}
      }
    }
  }

  if (flagMatrixTrigger)
    event.flagMatrixTrig = 1;
  if (flagSftTrigger)
    event.flagSftTrig = 1;
  if (flagMassTrigger)
    event.flagMassTrig = 1;


  DCAna->TrackSearchSdcInFiber();
  int ntSdcIn = DCAna->GetNtracksSdcIn();
  event.ntSdcIn = ntSdcIn;
  HF1(10, ntSdcIn);
  //std::cout << "ntSdcIn : " << ntSdcIn << std::endl;

  for( int it=0; it<ntSdcIn; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcIn(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();

    HF1( 11, double(nh) );
    HF1( 12, chisqr );
    HF1( 14, x0 ); HF1( 15, y0 );
    HF1( 16, u0 ); HF1( 17, v0 );
    HF2( 18, x0, u0 ); HF2( 19, y0, v0 );
    HF2( 20, x0, y0 );

    if (it<MaxHits) {
      event.u0SdcIn[it] = u0;
      event.v0SdcIn[it] = v0;
      event.x0SdcIn[it] = x0;
      event.y0SdcIn[it] = y0;
    }

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer();
      HF1( 13, layerId );
      int wire=(int)hit->GetWire();
      double dl=hit->GetDriftLength();
      double pos = hit->GetLocalHitPos();

      HF1( 100+10*layerId+1, wire-0.5 );
      HF1( 100+10*layerId+2, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double res=hit->GetResidual();
      HF1( 100+10*layerId+3, pos );
      HF1( 100+10*layerId+4, res );
      HF2( 100+10*layerId+5, pos, res );
      HF2( 100+10*layerId+6, xcal, ycal);
    }
  }

  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );

  DCAna->TrackSearchSdcOut();
  int ntSdcOut = DCAna->GetNtracksSdcOut();
  event.ntSdcOut = ntSdcOut;
  //std::cout << "ntSdcOut : " << ntSdcOut << std::endl;

  HF1( 30, double(ntSdcOut) );
  for( int it=0; it<ntSdcOut; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double x0=tp->GetX(zTof), y0=tp->GetY(zTof);

    HF1( 31, double(nh) );
    HF1( 32, chisqr );
    HF1( 34, x0 ); HF1( 35, y0 );
    HF1( 36, u0 ); HF1( 37, v0 );
    HF2( 38, x0, u0 ); HF2( 39, y0, v0 );
    HF2( 40, x0, y0 );

    if (it<MaxHits) {
      event.u0SdcOut[it] = u0;
      event.v0SdcOut[it] = v0;
      event.xtofSdcOut[it] = x0;
      event.ytofSdcOut[it] = y0;
    }

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      //int layerId=hit->GetLayer()-30;
      int layerId=hit->GetLayer()-26;
      HF1( 33, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 200+10*layerId+1, wire-0.5 );
      HF1( 200+10*layerId+2, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 200+10*layerId+3, pos );
      HF1( 200+10*layerId+4, res );
      HF2( 200+10*layerId+5, pos, res );
      HF2( 200+10*layerId+6, xcal, ycal);
    }
  }

  //DCAna->TrackSearchSks();
  DCAna->TrackSearchSksTmp(event.momScatMeson);
  int ntSks = DCAna->GetNTracksSks();
  event.ntSks = ntSks;

  std::vector <ThreeVector> KmomCont, KposCont;
  std::vector <SksParticle *> SksPartCont;

  //std::cout << "ntSks : " << ntSks << std::endl;
  HF1( 50, double(ntSks) );
  for( int i=0; i<ntSks; ++i ){
    SksTrack *tp=DCAna->GetSksTrack(i);
    if(!tp) continue;
    int nh=tp->GetNHits();
    double chisqr=tp->chisqr();
    ThreeVector Ppos=tp->PrimaryPosition();
    ThreeVector Pmom=tp->PrimaryMomentum();
    double pathL=tp->PathLengthToTOF();
    double xt=Ppos.x(), yt=Ppos.y();
    double p=Pmom.mag();
    double ut=Pmom.x()/p, vt=Pmom.y()/p;
    double cost = 1./sqrt(1.+ut*ut+vt*vt);
    double theta = acos(cost)*Rad2Deg;
    double phi   = atan2( ut, vt );

    event.p[i] = p;

    HF1( 51, double(nh) );
    HF1( 52, chisqr );
    HF1( 54, xt ); HF1( 55, yt ); HF1( 56, ut ); HF1( 57,vt );
    HF2( 58, xt, ut ); HF2( 59, yt, vt ); HF2( 60, xt, yt );
    HF1( 61, p ); HF1( 62, pathL );
    HF1(63, theta);

    for( int j=0; j<nh; ++j ){
      TrackHit *hit=tp->GetHit(j);
      if(!hit) continue;
      int layerId=hit->GetLayer();
      HF1( 53, layerId );

      double wire=hit->GetHit()->GetWire();
      double dt=hit->GetHit()->GetDriftTime();
      double dl=hit->GetHit()->GetDriftLength();
      double pos=hit->GetCalLPos(), res=hit->GetResidual();
      DCLTrackHit *lhit=hit->GetHit();
      double xcal=lhit->GetXcal(), ycal=lhit->GetYcal();
      //std::cout << "layerId : " << layerId << ", " << res << std::endl;

      if (layerId<15) {
	HF1( 400+10*layerId+1, double(wire)-0.5 );
	HF1( 400+10*layerId+2, dl );
	HF1( 400+10*layerId+3, pos );
	HF1( 400+10*layerId+4, res );
	HF2( 400+10*layerId+5, pos, res );
	HF2( 400+10*layerId+6, xcal, ycal );
      } else if (layerId>26){
	layerId -= 26;
	HF1( 500+10*layerId+1, double(wire)-0.5 );
	HF1( 500+10*layerId+2, dl );
	HF1( 500+10*layerId+3, pos ); 
	HF1( 500+10*layerId+4, res );
	HF2( 500+10*layerId+5, pos, res );
	HF2( 500+10*layerId+6, xcal, ycal );
      }
    }

    DCLocalTrack *tout = tp->GetLocalTrackOut();
    double xtof = tout->GetFTOFX();
    for (int j=0; j<ncTof; ++j) {
      HodoCluster *clTof=hodoAna->GetClusterTOF(j);
      if (!clTof || !clTof->GoodForAnalysis() ) continue;

      double meanSeg = clTof->MeanSeg();
      double tofSegPos = DCGeomMan::GetInstance().calcWirePosition(IdTof, meanSeg);

      double diffX = xtof-tofSegPos;
      HF1(64, diffX);
      //if (fabs(diffX)<60.) {
      if (fabs(diffX)<6000.) {
	if (chisqr<20) {
	  SksParticle *SksPart = new SksParticle(tp, clTof, tout);
	  SksPartCont.push_back(SksPart);
	}

	double m2 = MassSquare(p, pathL, clTof->CMeanTime()+time_offset-time_offset_tgt0+gRand->Gaus(0, time_reso_m2));
	HF1(65, m2);

	if (m2>MinOfMassSquareK && m2 < MaxOfMassSquareK) {
	  KmomCont.push_back(Pmom);
	  KposCont.push_back(Ppos);
	}

	if (m2>0)
	  HF1(66, sqrt(m2));
	
      }
    }
  }

  int ntSksPart = SksPartCont.size();
  event.ntSksPart = ntSksPart;
  for (int i=0; i<ntSksPart; i++) {
    event.pSksPart[i] = SksPartCont[i]->Momentum().mag();
    double tof = SksPartCont[i]->GetTof()->CMeanTime();
    double pathL = SksPartCont[i]->PathLengthToTOF();
    double m2 = MassSquare(event.pSksPart[i], pathL, tof+time_offset-time_offset_tgt0+gRand->Gaus(0, time_reso_m2));
    event.m2[i] = m2;
    //event.m2[i] = SksPartCont[i]->MassSquare();
  }



  DCAna->TrackSearchBcOut();
  int ntBcOut = DCAna->GetNtracksBcOut();

  std::vector <ThreeVector> BmomCont, BposCont;

  int layerId_SftU = DCGeomMan::GetInstance().GetDetectorId("SFT-u-1");
  int layerId_SftV = DCGeomMan::GetInstance().GetDetectorId("SFT-v-1");
  double localZ_SftU = DCGeomMan::GetInstance().GetLocalZ(layerId_SftU);
  double localZ_SftV = DCGeomMan::GetInstance().GetLocalZ(layerId_SftV);
  double tiltAngle_SftU = DCGeomMan::GetInstance().GetTiltAngle(layerId_SftU);
  double tiltAngle_SftV = DCGeomMan::GetInstance().GetTiltAngle(layerId_SftV);
  DCHitContainer SFTU_cont = DCAna->GetSdcInHC(1);
  DCHitContainer SFTV_cont = DCAna->GetSdcInHC(2);

  HF1( 70, double(ntBcOut) );
  for( int it=0; it<ntBcOut; ++it ){
    DCLocalTrack *tp=DCAna->GetTrackBcOut(it);
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double x0=tp->GetX0(), y0=tp->GetY0();

    HF1( 71, double(nh) );
    HF1( 72, chisqr );
    HF1( 74, x0 ); HF1( 75, y0 );
    HF1( 76, u0 ); HF1( 77, v0 );
    HF2( 78, x0, u0 ); HF2( 39, y0, v0 );
    HF2( 80, x0, y0 );

    ThreeVector bmom(event.beammom*u0/sqrt(1+u0*u0+v0*v0),
		     event.beammom*v0/sqrt(1+u0*u0+v0*v0), 
		     event.beammom/sqrt(1+u0*u0+v0*v0));
    ThreeVector bpos(x0, y0, 0);
    BmomCont.push_back(bmom);
    BposCont.push_back(bpos);

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer()-112;
      HF1( 73, layerId );
      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1( 700+10*layerId+1, wire-0.5 );
      HF1( 700+10*layerId+2, dl );
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 700+10*layerId+3, pos );
      HF1( 700+10*layerId+4, res );
      HF2( 700+10*layerId+5, pos, res );
      HF2( 700+10*layerId+6, xcal, ycal);
    }

    // BcOut at SFT-U
    double  x1 = u0*localZ_SftU + x0;
    double  y1 = v0*localZ_SftU + y0;
    double  xp1 = x1*cos(-tiltAngle_SftU*Deg2Rad)-y1*sin(-tiltAngle_SftU*Deg2Rad);
    int nhSftU = SFTU_cont.size();
    double resU = 9999.;
    for (int ih=0; ih<nhSftU; ih++) {
      DCHit *hit = SFTU_cont[ih];
      double localPos = hit->GetWirePosition();
      double res = xp1 - localPos;
      if (fabs(res) < resU)
	resU = res;
    }
    HF1(81, resU);

    // BcOut at SFT-V
    double  x2 = u0*localZ_SftV + x0;
    double  y2 = v0*localZ_SftV + y0;
    double  xp2 = x2*cos(-tiltAngle_SftV*Deg2Rad)-y2*sin(-tiltAngle_SftV*Deg2Rad);

    int nhSftV = SFTV_cont.size();
    double resV = 9999.;
    for (int ih=0; ih<nhSftV; ih++) {
      DCHit *hit = SFTV_cont[ih];
      double localPos = hit->GetWirePosition();
      double res = xp2 - localPos;
      if (fabs(res) < resV)
	resV = res;
    }
    HF1(82, resV);    

  }

  int nBeam=BmomCont.size();
  int nK=KmomCont.size();
  event.nK = nK;

  if (nBeam == 1) {
    ThreeVector bMom = BmomCont[0];
    ThreeVector bPos = BposCont[0];

    LorentzVector LvBeam(  bMom, sqrt(PionMass*PionMass+bMom.mag2()) );
    LorentzVector LvTgt( 0., 0., 0., ProtonMass );
    //std::cout << "beammom : " << event.beammom << std::endl;
    ThreeVector vtx0(event.primaryVertex[0], event.primaryVertex[1], event.primaryVertex[2]);
    
    for (int i=0; i<nK; i++) {
      ThreeVector priMom=KmomCont[i];
      ThreeVector kPos=KposCont[i];
      event.pK[i] = priMom.mag();

      ThreeVector vtx = VertexPoint(bPos, kPos, bMom, priMom);
      event.Vertex_x[i] = vtx.x();
      event.Vertex_y[i] = vtx.y();
      event.Vertex_z[i] = vtx.z();

      LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
      LorentzVector LvRc = LvBeam+LvTgt-LvScat;
      double mismass = LvRc.mag();

      HF1(67, mismass);
      event.MissMass[i] = mismass;
      
      double p0 = priMom.mag();
      double p_Cor, E_Cor;    
      CorrElossOut(&p_Cor, &E_Cor, p0, KAON, priMom/p0, vtx0);
      priMom *= p_Cor/p0;
      event.pKCor[i] = priMom.mag();
      
      LorentzVector LvScatCor( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
      LorentzVector LvRcCor = LvBeam+LvTgt-LvScatCor;
      double mismassCor = LvRcCor.mag();
      HF1(68, mismassCor);
      event.MissMassCor[i] = mismassCor;
      ThreeVector pSigmaCor = bMom-priMom;
      event.pSigmaCor[i] = pSigmaCor.mag();
      
      double cost = bMom*priMom/(bMom.mag()*priMom.mag());
      event.theta[i] = acos(cost)*Rad2Deg;
      
      double calMom = calcKMomFromTheta(event.beammom, cost, p_Cor);
      event.pKCal[i] = calMom;
      ThreeVector KaonMomCal = priMom*(calMom/priMom.mag());
      ThreeVector pSigmaCal = bMom-KaonMomCal;
      event.pSigmaCal[i] = pSigmaCal.mag();
      
    }
  }

  int nSksPart=SksPartCont.size();
  HF1(90, nSksPart);
  for (int i=0; i<nSksPart; i++) {
    SksParticle *SksPart = SksPartCont[i];
    double m2 = SksPart->MassSquare();
    HF1(91, m2);
  }
  if (nSksPart==2) {
    double m2_1 = SksPartCont[0]->MassSquare();
    double m2_2 = SksPartCont[1]->MassSquare();

    if (m2_1 > m2_2)
      HF2(92, m2_1, m2_2);
    else
      HF2(92, m2_2, m2_1);
  }

  tree->Fill();

  for_each(SksPartCont.begin(), SksPartCont.end(), DeleteObject());

  delete hodoAna;
  delete DCAna;
  delete rawData;

  return true;
}

double calcKMomFromTheta(double bmom, double cost, double p0)
{
  double E1 = sqrt(PionMass*PionMass + bmom*bmom);
  double A = PionMass*PionMass + ProtonMass*ProtonMass
    + KaonMass*KaonMass - SigmaMinusMass*SigmaMinusMass
    + 2*E1*ProtonMass;

  double hanbetsu = bmom*bmom*A*A*cost*cost-
    ((ProtonMass+E1)*(ProtonMass+E1)-bmom*bmom*cost*cost)*
    (4.*(ProtonMass+E1)*(ProtonMass+E1)*KaonMass*KaonMass-A*A);

  if (hanbetsu<0)
    return 0.;

  double p1, p2;

  p1 = (bmom*A*cost+sqrt(hanbetsu))/(2.*((ProtonMass+E1)*(ProtonMass+E1)-bmom*bmom*cost*cost));
  p2 = (bmom*A*cost-sqrt(hanbetsu))/(2.*((ProtonMass+E1)*(ProtonMass+E1)-bmom*bmom*cost*cost));

  //std::cout << "p1 : " << p1 << ", p2 : " << p2 << std::endl;

  if (fabs(p0-p1)<fabs(p0-p2))
    return p1;
  else
    return p2;
}

void InitializeEvent( void )
{
  //Primary for Yd scattering
  event.thetaMeson = -999.9;
  event.phiMeson = -999.9;
  event.thetaMesonCM = -999.9;
  event.phiMesonCM = -999.9;
  event.thetaScatHypCM = -999.9;
  event.phiScatHypCM = -999.9;
  event.thetaScatPLab = -999.9;
  event.thetaScatHypLab = -999.9;
  event.momVectorScatMeson[0] = -999.9;
  event.momVectorScatMeson[1] = -999.9;
  event.momVectorScatMeson[2] = -999.9;
  event.momVectorHypBeam[0] = -999.9;
  event.momVectorHypBeam[1] = -999.9;
  event.momVectorHypBeam[2] = -999.9;
  event.momVectorHypScat[0] = -999.9;
  event.momVectorHypScat[1] = -999.9;
  event.momVectorHypScat[2] = -999.9;
  event.momVectorProtonScat[0] = -999.9;
  event.momVectorProtonScat[1] = -999.9;
  event.momVectorProtonScat[2] = -999.9;
  event.momVectorDecayPi[0] = -999.9;
  event.momVectorDecayPi[1] = -999.9;
  event.momVectorDecayPi[2] = -999.9;
  event.momVectorDecayNucleon[0] = -999.9;
  event.momVectorDecayNucleon[1] = -999.9;
  event.momVectorDecayNucleon[2] = -999.9;
  event.momScatMeson = -999.9;
  event.momHypBeam = -999.9;
  event.momHypScat = -999.9;
  event.momProtonScat = -999.9;
  event.momDecayPi = -999.9;
  event.momDecayNucleon = -999.9;
  event.primaryVertex[0] = -999.9;
  event.primaryVertex[1] = -999.9;
  event.primaryVertex[2] = -999.9;
  event.scatPos0[0] = -999.9;
  event.scatPos0[1] = -999.9;
  event.scatPos0[2] = -999.9;
  event.NNscatPos[0] = -999.9;
  event.NNscatPos[1] = -999.9;
  event.NNscatPos[2] = -999.9;
  event.PiNscatPos[0] = -999.9;
  event.PiNscatPos[1] = -999.9;
  event.PiNscatPos[2] = -999.9;
  event.decayPos[0] = -999.9;
  event.decayPos[1] = -999.9;
  event.decayPos[2] = -999.9;
  event.decayFlag = -999;
  event.scatFlag = -999;
  event.scatTarget = -999;
  event.NNscatFlag = -999;
  event.NNscatTarget = -999;
  event.PiNscatFlag = -999;
  event.PiNscatTarget = -999;

  for (int i=0; i<NumOfPlaneCFT; i++) 
    event.FiberHits[i] = 0;
  for (int i=0; i<NumOfPlaneCFT; i++) {
    for (int j=0; j<MaxHits2; j++) {
      event.FiberSeg[i][j] = -999.0;
      event.FiberTime[i][j] = -999.0;
      event.FiberEdep[i][j] = -999.0;
      event.FiberPID[i][j] = -1;
    }
  }

  event.CrystalHits = 0;
  for (int j=0; j<MaxHits2; j++) {
    event.CrystalSeg [j] = -999.0;
    event.CrystalTime[j] = -999.0;
    event.CrystalEdep[j] = -999.0;
    event.CrystalPID [j] = -1;
  }


  for (int i=0; i<NumOfPlaneCFT; i++) {
    event.fClEdepForScat[i] = 0.;
    event.fClEdepNormForScat[i] = 0.;
    event.fClPathLengthForScat[i] = -999.0;

    event.fClEdepForPi[i] = 0.;
    event.fClEdepNormForPi[i] = 0.;
    event.fClPathLengthForPi[i] = -999.0;
  }

  event.cClEdepForScat = 0.0;
  event.cClEdepForPi = 0.0;
  event.cRcEdepForScat = 0.0;
  event.cRcEdepForPi = 0.0;

  event.nhTof = 0;
  for (int i=0; i<MaxHits; i++) {
    event.TofSeg[i] = -999;
    event.TofTime[i] = -999.;
    event.TofEdep[i] = -999.;
  }

  event.nhLc = 0;
  for (int i=0; i<MaxHits; i++) {
    event.LcSeg[i] = -999;
    event.LcTime[i] = -999.;
    event.LcEdep[i] = -999.;
  }

  event.nhCh = 0;
  for (int i=0; i<MaxHits; i++) {
    event.ChSeg[i] = -999;
    event.ChTime[i] = -999.;
    event.ChEdep[i] = -999.;
  }

  event.nhAc = 0;
  for (int i=0; i<MaxHits; i++) {
    event.AcSeg[i] = -999;
    event.AcTime[i] = -999.;
    event.AcEdep[i] = -999.;
  }

  event.nhSftx = 0;
  for (int i=0; i<MaxHits; i++) {
    event.SftxSeg[i] = -999.;
  }

  event.ntSdcIn = 0;
  for (int i=0; i<MaxHits; i++) {
    event.u0SdcIn[i] = -999.;
    event.v0SdcIn[i] = -999.;
    event.x0SdcIn[i] = -999.;
    event.y0SdcIn[i] = -999.;
  }

  event.ntSdcOut = 0;
  for (int i=0; i<MaxHits; i++) {
    event.u0SdcOut[i] = -999.;
    event.v0SdcOut[i] = -999.;
    event.xtofSdcOut[i] = -999.;
    event.ytofSdcOut[i] = -999.;
  }

  event.ntSks=0;
  event.ntSksPart=0;
  event.nK=0;
  for (int i=0; i<MaxHits; i++) {
    event.p[i]=-999.0;
    event.pK[i]=-999.0;
    event.pKCor[i]=-999.0;
    event.pKCal[i]=-999.0;
    event.theta[i]=-999.0;
    event.MissMass[i]=-999.0;
    event.MissMassCor[i]=-999.0;
    event.pSigmaCor[i]=-999.0;
    event.pSigmaCal[i]=-999.0;
    event.Vertex_x[i]=-999.0;
    event.Vertex_y[i]=-999.0;
    event.Vertex_z[i]=-999.0;
  }

  event.flagMatrixTrig = 0;
  event.flagSftTrig = 0;
  event.flagMassTrig = 0;
}


void DefineHistograms( const char *filename )
{
  new TFile( filename, "recreate" );

  char buf[100];

  HB1(10, "N track SdcIn", 10, 0, 10);
  HB1(11, "N Hit SdcIn tracking", 10, 0, 10);
  HB1(12, "Chisquare SdcIn tracking", 100, 0, 50);
  HB1(13, "Hit Layer  SdcIn tracking", 10, 0, 10);
  HB1(14, "X0  SdcIn tracking", 200, -50, 50);
  HB1(15, "Y0  SdcIn tracking", 200, -50, 50);
  HB1(16, "U0  SdcIn tracking", 200, -0.5, 0.5);
  HB1(17, "V0  SdcIn tracking", 200, -0.5, 0.5);
  HB2(18, "X0%U0  SdcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(19, "Y0%V0  SdcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(20, "X0%Y0  SdcIn tracking", 200, -50, 50, 200, -50, 50);

  for (int i=1; i<=NumOfLayersSdcIn; i++) {
    int hid;
    sprintf(buf, "Hit wire SdcIn Layer  %d", i);
    hid = 100 + i*10+1;
    int wireNum;
    if (i==1)
      wireNum=512;
    else if (i==2 || i==3)
      wireNum=320;
    else if (i>=4 && i<=9)
      wireNum=50;

    HB1(hid, buf, wireNum, 0, wireNum);

    hid = 100 + i*10+2;
    sprintf(buf, "Drift length SdcIn Layer  %d", i);    
    HB1(hid, buf, 100, -2, 6);

    hid = 100 + i*10+3;
    sprintf(buf, "Local Pos  SdcIn Layer  %d", i);    
    HB1(hid, buf, 1000, -300, 300);

    hid = 100 + i*10+4;
    sprintf(buf, "Residual  SdcIn Layer  %d", i);    
    HB1(hid, buf, 400, -2, 2);

    hid = 100 + i*10+5;
    sprintf(buf, "Pos%Residual  SdcIn Layer  %d", i);    
    HB2(hid, buf, 300, -300, 300, 400, -2, 2);

    hid = 100 + i*10+6;
    sprintf(buf, "Ycal%Xcal  SdcIn Layer  %d", i);    
    HB2(hid, buf, 300, -300, 300, 300, -300, 300);

  }

  HB1(30, "N track SdcOut", 10, 0, 10);
  HB1(31, "N Hit SdcOut tracking", 12, 0, 12);
  HB1(32, "Chisquare SdcOut tracking", 100, 0, 50);
  HB1(33, "Hit Layer  SdcOut tracking", 20, 0, 20);
  HB1(34, "X0 (TOF)  SdcOut tracking", 200, -600, 600);
  HB1(35, "Y0 (TOF)  SdcOut tracking", 200, -600, 600);
  HB1(36, "U0  SdcOut tracking", 200, -1., 1.);
  HB1(37, "V0  SdcOut tracking", 200, -1., 1.);
  HB2(38, "X0%U0  SdcOut tracking", 200, -600, 600, 200, -1, 1);
  HB2(39, "Y0%V0  SdcOut tracking", 200, -600, 600, 200, -1, 1);
  HB2(40, "X0%Y0  SdcOut tracking", 200, -600, 600, 200, -600, 600);

  for (int i=1; i<=NumOfLayersSdcOut; i++) {
    int hid;
    sprintf(buf, "Hit wire SdcOut Layer  %d", i);
    hid = 200 + i*10+1;
    HB1(hid, buf, 120, 0, 120);

    hid = 200 + i*10+2;
    sprintf(buf, "Drift length SdcOut Layer  %d", i);    
    HB1(hid, buf, 100, -2, 12);

    hid = 200 + i*10+3;
    sprintf(buf, "Local Pos  SdcOut Layer  %d", i);    
    HB1(hid, buf, 1000, -600, 600);

    hid = 200 + i*10+4;
    sprintf(buf, "Residual  SdcOut Layer  %d", i);    
    HB1(hid, buf, 400, -2, 2);

    hid = 200 + i*10+5;
    sprintf(buf, "Pos%Residual  SdcOut Layer  %d", i);    
    HB2(hid, buf, 300, -600, 600, 400, -2, 2);

    hid = 200 + i*10+6;
    sprintf(buf, "Ycal%Xcal  SdcOut Layer  %d", i);    
    HB2(hid, buf, 300, -600, 600, 300, -600, 600);

  }

  HB1(50, "N track Sks", 10, 0, 10);
  HB1(51, "N Hit Sks tracking", 40, 0, 40);
  HB1(52, "Chisquare Sks tracking", 400, 0, 100);
  HB1(53, "Hit Layer  Sks tracking", 40, 0, 40);
  HB1(54, "X0 (Target)  Sks tracking", 200, -50, 50);
  HB1(55, "Y0 (Target)  Sks tracking", 200, -50, 50);
  HB1(56, "U0  Sks tracking", 200, -0.5, 0.5);
  HB1(57, "V0  Sks tracking", 200, -0.5, 0.5);
  HB2(58, "X0%U0  Sks tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(59, "Y0%V0  Sks tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(60, "X0%Y0  Sks tracking", 200, -50, 50, 200, -50, 50);
  HB1(61, "Momentum", 200, 0.4, 1.0);
  HB1(62, "Path Length", 1000, 0., 5000.);
  HB1(63, "Theta", 300, 0., 30.);
  HB1(64, "X Diff at TOF", 200, -100., 100.);
  HB1(65, "Mass square", 400, -0.5, 1.5);
  HB1(66, "Mass", 400, -0.5, 1.5);
  HB1(67, "MissMass", 400, 0.5, 1.5);
  HB1(68, "MissMass (Eloss Corr)", 400, 0.5, 1.5);

  for (int i=1; i<=NumOfLayersSdcIn; i++) {
    int hid;
    sprintf(buf, "Hit wire SdcIn (SKS) Layer  %d", i);
    hid = 400 + i*10+1;
    int wireNum;
    if (i==1)
      wireNum=512;
    else if (i==2 || i==3)
      wireNum=320;
    else if (i>=4 && i<=9)
      wireNum=50;

    HB1(hid, buf, wireNum, 0, wireNum);

    hid = 400 + i*10+2;
    sprintf(buf, "Drift length SdcIn (SKS) Layer  %d", i);    
    HB1(hid, buf, 100, -2, 6);

    hid = 400 + i*10+3;
    sprintf(buf, "Local Pos  SdcIn (SKS) Layer  %d", i);    
    HB1(hid, buf, 1000, -300, 300);

    hid = 400 + i*10+4;
    sprintf(buf, "Residual  SdcIn (SKS) Layer  %d", i);    
    HB1(hid, buf, 400, -2, 2);

    hid = 400 + i*10+5;
    sprintf(buf, "Pos%Residual  SdcIn (SKS) Layer  %d", i);    
    HB2(hid, buf, 300, -300, 300, 400, -2, 2);

    hid = 400 + i*10+6;
    sprintf(buf, "Ycal%Xcal  SdcIn (SKS) Layer  %d", i);    
    HB2(hid, buf, 300, -300, 300, 300, -300, 300);

  }

  for (int i=1; i<=NumOfLayersSdcOut; i++) {
    int hid;
    sprintf(buf, "Hit wire SdcOut (SKS) Layer  %d", i);
    hid = 500 + i*10+1;
    HB1(hid, buf, 120, 0, 120);

    hid = 500 + i*10+2;
    sprintf(buf, "Drift length SdcOut (SKS) Layer  %d", i);    
    HB1(hid, buf, 100, -2, 12);

    hid = 500 + i*10+3;
    sprintf(buf, "Local Pos  SdcOut (SKS) Layer  %d", i);    
    HB1(hid, buf, 1000, -600, 600);

    hid = 500 + i*10+4;
    sprintf(buf, "Residual  SdcOut (SKS) Layer  %d", i);    
    HB1(hid, buf, 400, -2, 2);

    hid = 500 + i*10+5;
    sprintf(buf, "Pos%Residual  SdcOut (SKS) Layer  %d", i);    
    HB2(hid, buf, 300, -600, 600, 400, -2, 2);

    hid = 500 + i*10+6;
    sprintf(buf, "Ycal%Xcal  SdcOut (SKS) Layer  %d", i);    
    HB2(hid, buf, 300, -600, 600, 300, -600, 600);

  }


  HB1(70, "N track BcOut", 10, 0, 10);
  HB1(71, "N Hit BcOut tracking", 10, 0, 10);
  HB1(72, "Chisquare BcOut tracking", 100, 0, 50);
  HB1(73, "Hit Layer  BcOut tracking", 10, 0, 10);
  HB1(74, "X0  BcOut tracking", 200, -50, 50);
  HB1(75, "Y0  BcOut tracking", 200, -50, 50);
  HB1(76, "U0  BcOut tracking", 200, -0.5, 0.5);
  HB1(77, "V0  BcOut tracking", 200, -0.5, 0.5);
  HB2(78, "X0%U0  BcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(79, "Y0%V0  BcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(80, "X0%Y0  BcIn tracking", 200, -50, 50, 200, -50, 50);
  HB1(81, "Residual SFT-U  BcOut tracking", 400, -50, 50);
  HB1(82, "Residual SFT-V  BcOut tracking", 400, -50, 50);

  for (int i=1; i<=NumOfLayersBcOut; i++) {
    int hid;
    sprintf(buf, "Hit wire BcOut Layer  %d", i);
    hid = 700 + i*10+1;
    HB1(hid, buf, 120, 0, 120);

    hid = 700 + i*10+2;
    sprintf(buf, "Drift length BcOut Layer  %d", i);    
    HB1(hid, buf, 100, -2, 12);

    hid = 700 + i*10+3;
    sprintf(buf, "Local Pos  BcOut Layer  %d", i);    
    HB1(hid, buf, 1000, -600, 600);

    hid = 700 + i*10+4;
    sprintf(buf, "Residual  BcOut Layer  %d", i);    
    HB1(hid, buf, 400, -2, 2);

    hid = 700 + i*10+5;
    sprintf(buf, "Pos%Residual  BcOut Layer  %d", i);    
    HB2(hid, buf, 300, -600, 600, 400, -2, 2);

    hid = 700 + i*10+6;
    sprintf(buf, "Ycal%Xcal  BcOut Layer  %d", i);    
    HB2(hid, buf, 300, -600, 600, 300, -600, 600);

  }

  HB1(90, "N SksParticle", 10, 0, 10);
  HB1(91, "SksParticle MassSquare", 200, -0.2, 1.2);
  HB2(92, "SksParticle MassSquare (2 track)", 200, -0.2, 1.2, 200, -0.2, 1.2);

  HB1(1001, "Tof Nhit", 10, 0, 10);
  HB1(1002, "Tof Hit Segment", 24, 0, 24);
  HB1(1003, "Tof Time", 300, 0, 30);
  HB1(1004, "Tof dE", 300, 0, 30);

  HB1(1011, "CH Nhit", 10, 0, 10);
  HB1(1012, "CH Hit Segment", 24, 0, 24);
  HB1(1013, "CH Time", 200, 0, 20);
  HB1(1014, "CH dE", 200, 0, 10);


  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Primary for Yd scattering
  tree->Branch("thetaMeson",   &event.thetaMeson,  "thetaMeson/D");
  tree->Branch("phiMeson",   &event.phiMeson,  "phiMeson/D");
  tree->Branch("thetaMesonCM",   &event.thetaMesonCM,  "thetaMesonCM/D");
  tree->Branch("phiMesonCM",   &event.phiMesonCM,  "phiMesonCM/D");
  tree->Branch("thetaScatHypCM",   &event.thetaScatHypCM,  "thetaScatHypCM/D");
  tree->Branch("phiScatHypCM",   &event.phiScatHypCM,  "phiScatHypCM/D");
  tree->Branch("thetaScatPLab",   &event.thetaScatPLab,  "thetaScatPLab/D");
  tree->Branch("thetaScatHypLab",   &event.thetaScatHypLab,  "thetaScatHypLab/D");
  tree->Branch("momVectorScatMeson",   event.momVectorScatMeson,  "momVectorScatMeson[3]/D");
  tree->Branch("momVectorHypBeam",   event.momVectorHypBeam,  "momVectorHypBeam[3]/D");
  tree->Branch("momVectorHypScat",   event.momVectorHypScat,  "momVectorHypScat[3]/D");
  tree->Branch("momVectorProtonScat",   event.momVectorProtonScat,  "momVectorProtonScat[3]/D");
  tree->Branch("momVectorDecayPi",   event.momVectorDecayPi,  "momVectorDecayPi[3]/D");
  tree->Branch("momVectorDecayNucleon",   event.momVectorDecayNucleon,  "momVectorDecayNucleon[3]/D");
  tree->Branch("momScatMeson",   &event.momScatMeson,  "momScatMeson/D");
  tree->Branch("momHypBeam",   &event.momHypBeam,  "momHypBeam/D");
  tree->Branch("momHypScat",   &event.momHypScat,  "momHypScat/D");
  tree->Branch("momProtonScat",   &event.momProtonScat,  "momProtonScat/D");
  tree->Branch("momDecayPi",   &event.momDecayPi,  "momDecayPi/D");
  tree->Branch("momDecayNucleon",   &event.momDecayNucleon,  "momDecayNucleon/D");
  tree->Branch("primaryVertex",   event.primaryVertex,  "primaryVertex[3]/D");
  tree->Branch("scatPos0",   event.scatPos0,  "scatPos0[3]/D");
  tree->Branch("NNscatPos",   event.NNscatPos,  "NNscatPos[3]/D");
  tree->Branch("PiNscatPos",   event.PiNscatPos,  "PiNscatPos[3]/D");
  tree->Branch("decayPos",   event.decayPos,  "decayPos[3]/D");
  tree->Branch("decayFlag",   &event.decayFlag,  "decayFlag/I");
  tree->Branch("scatFlag",   &event.scatFlag,  "scatFlag/I");
  tree->Branch("scatTarget",   &event.scatTarget,  "scatTarget/I");
  tree->Branch("NNscatFlag",   &event.NNscatFlag,  "NNscatFlag/I");
  tree->Branch("NNscatTarget",   &event.NNscatTarget,  "NNscatTarget/I");
  tree->Branch("PiNscatFlag",   &event.PiNscatFlag,  "PiNscatFlag/I");
  tree->Branch("PiNscatTarget",   &event.PiNscatTarget,  "PiNscatTarget/I");

  for (int i=0; i<NumOfPlaneCFT; i++) {
    char buf1[100], buf2[100];
    sprintf(buf1, "FiberHits%d", i);
    sprintf(buf2, "FiberHits%d/I", i);
    tree->Branch(buf1,   &event.FiberHits[i],  buf2);

    sprintf(buf1, "FiberSeg%d", i);
    sprintf(buf2, "FiberSeg%d[FiberHits%d]/D", i, i);
    tree->Branch(buf1,   event.FiberSeg[i],  buf2);

    sprintf(buf1, "FiberTime%d", i);
    sprintf(buf2, "FiberTime%d[FiberHits%d]/D", i, i);
    tree->Branch(buf1,   event.FiberTime[i],  buf2);

    sprintf(buf1, "FiberEdep%d", i);
    sprintf(buf2, "FiberEdep%d[FiberHits%d]/D", i, i);
    tree->Branch(buf1,   event.FiberEdep[i],  buf2);

  }

  tree->Branch("CrystalHits",   &event.CrystalHits,  "CrystalHits/I");
  tree->Branch("CrystalSeg",   &event.CrystalSeg,  "CrystalSeg[CrystalHits]/D");
  tree->Branch("CrystalTime",   &event.CrystalTime,  "CrystalTime[CrystalHits]/D");
  tree->Branch("CrystalEdep",   &event.CrystalEdep,  "CrystalEdep[CrystalHits]/D");
  tree->Branch("CrystalPID",   &event.CrystalPID,  "CrystalPID[CrystalHits]/I");

  tree->Branch("nhTof",   &event.nhTof,  "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/I");
  tree->Branch("TofTime",   event.TofTime,  "TofTime[nhTof]/D");
  tree->Branch("TofEdep",   event.TofEdep,  "TofEdep[nhTof]/D");

  tree->Branch("nhLc",   &event.nhLc,  "nhLc/I");
  tree->Branch("LcSeg",   event.LcSeg,  "LcSeg[nhLc]/I");
  tree->Branch("LcTime",   event.LcTime,  "LcTime[nhLc]/D");
  tree->Branch("LcEdep",   event.LcEdep,  "LcEdep[nhLc]/D");

  tree->Branch("nhCh",   &event.nhCh,  "nhCh/I");
  tree->Branch("ChSeg",   event.ChSeg,  "ChSeg[nhCh]/I");
  tree->Branch("ChTime",   event.ChTime,  "ChTime[nhCh]/D");
  tree->Branch("ChEdep",   event.ChEdep,  "ChEdep[nhCh]/D");

  tree->Branch("nhAc",   &event.nhAc,  "nhAc/I");
  tree->Branch("AcSeg",   event.AcSeg,  "AcSeg[nhAc]/I");
  tree->Branch("AcTime",   event.AcTime,  "AcTime[nhAc]/D");
  tree->Branch("AcEdep",   event.AcEdep,  "AcEdep[nhAc]/D");

  tree->Branch("nhSftx",   &event.nhSftx,  "nhSftx/I");
  tree->Branch("SftxSeg",   event.SftxSeg,  "SftxSeg[nhSftx]/D");

  tree->Branch("ntSdcIn",   &event.ntSdcIn,  "ntSdcIn/I");
  tree->Branch("u0SdcIn",   event.u0SdcIn,  "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",   event.v0SdcIn,  "v0SdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",   event.x0SdcIn,  "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",   event.y0SdcIn,  "y0SdcIn[ntSdcIn]/D");

  tree->Branch("ntSdcOut",   &event.ntSdcOut,  "ntSdcOut/I");
  tree->Branch("u0SdcOut",   event.u0SdcOut,  "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",   event.v0SdcOut,  "v0SdcOut[ntSdcOut]/D");
  tree->Branch("xtofSdcOut",   event.xtofSdcOut,  "xtofSdcOut[ntSdcOut]/D");
  tree->Branch("ytofSdcOut",   event.ytofSdcOut,  "ytofSdcOut[ntSdcOut]/D");

  tree->Branch("ntSks",   &event.ntSks,  "ntSks/I");
  tree->Branch("p",   &event.p,  "p[ntSks]/D");

  tree->Branch("ntSksPart",   &event.ntSksPart,  "ntSksPart/I");
  tree->Branch("pSksPart",   &event.pSksPart,  "pSksPart[ntSksPart]/D");
  tree->Branch("m2",   &event.m2,  "m2[ntSksPart]/D");

  tree->Branch("ntK",   &event.nK,  "nK/I");
  tree->Branch("pK",   &event.pK,  "pK[nK]/D");
  tree->Branch("pKCor",   &event.pKCor,  "pKCor[nK]/D");
  tree->Branch("pKCal",   &event.pKCal,  "pKCal[nK]/D");
  tree->Branch("theta",   &event.theta,  "theta[nK]/D");
  tree->Branch("MissMass",   &event.MissMass,  "MissMass[nK]/D");
  tree->Branch("MissMassCor",   &event.MissMassCor,  "MissMassCor[nK]/D");
  tree->Branch("pSigmaCor",   &event.pSigmaCor,  "pSigmaCor[nK]/D");
  tree->Branch("pSigmaCal",   &event.pSigmaCal,  "pSigmaCal[nK]/D");
  tree->Branch("Vertex_x",   &event.Vertex_x,  "Vertex_x[nK]/D");
  tree->Branch("Vertex_y",   &event.Vertex_y,  "Vertex_y[nK]/D");
  tree->Branch("Vertex_z",   &event.Vertex_z,  "Vertex_z[nK]/D");

  tree->Branch("flagMatrixTrig",   &event.flagMatrixTrig,  "flagMatrixTrig/I");
  tree->Branch("flagSftTrig",   &event.flagSftTrig,  "flagSftTrig/I");
  tree->Branch("flagMassTrig",   &event.flagMassTrig,  "flagMassTrig/I");
}

bool ConfMan::InitializeParameterFiles(void)
{

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( DCGeomFileName_!="" )
    DCGeomManager_->Initialize(DCGeomFileName_);
  else
    DCGeomManager_->Initialize();

  if( FieldMapFileName_!="" ) {
    FieldMan::GetInstance().Initialize( FieldMapFileName_ );
    std::cout << "MapScale : " << MapScale_ << std::endl;
    FieldMan::GetInstance().SetMapScale(MapScale_);
  }


  if( mtTrigFileName_!="" )
    MatrixTrigMan::GetInstance().Initialize( mtTrigFileName_ );

  if( sftTrigFileName_!="" )
    MatrixTrigMan::GetInstance().InitializeSftTrig( sftTrigFileName_ );

  if( massTrigFileName_!="" )
    MatrixTrigMan::GetInstance().InitializeMassTrig( massTrigFileName_ );


  return true;
}
