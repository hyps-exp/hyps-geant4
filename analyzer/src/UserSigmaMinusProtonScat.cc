/*
  SksTracking.cc
*/

#include "ThreeVector.hh"
#include "LorentzVector.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "HodoAnalyzer.hh"
#include "CFTFiberHit.hh"
#include "CFTFLHit.hh"
#include "CFTFiberCluster.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoCluster.hh"
#include "DCAnalyzer.hh"
#include "DCLocalTrack.hh"
#include "CFTLocalTrack.hh"
#include "SksTrack.hh"
#include "TrackHit.hh"
#include "Kinematics.hh"
#include "EnergyCorrection.hh"
#include "SksParticle.hh"
#include "CFTParticle.hh"
#include "TemplateLib.hh"

#include "HistHelper.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FieldMan.hh"
#include "MatrixTrigMan.hh"

#include "EvDispCFT.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <signal.h>

#include <TRandom.h>
#include <time.h>

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

const double MinOfMassSigma = 1.15;
const double MaxOfMassSigma = 1.25;

//const double MinOfDeltaU_CFT = 0.1;
//const double MaxOfDeltaU_CFT = 0.18;
const double MinOfDeltaU_CFT = -0.04;
const double MaxOfDeltaU_CFT = 0.04;
//const double MinOfDeltaV_CFT = -0.03;
//const double MaxOfDeltaV_CFT = 0.03;
const double MinOfDeltaV_CFT = -0.05;
const double MaxOfDeltaV_CFT = 0.05;

bool ProcessOneEvent( std::ifstream & );
void DefineHistograms( const char * );
void InitializeEvent( void );
bool calcDecayPiMom(double p1, double m1, double m2, double m3, double cost, double *momCal);
bool calcBeamMomFrom2BodyKinema(double M1, double M2, double theta, double pScat,
				double *beamMomCal1, double *beamMomCal2);
bool calc2BodyKinema(double M1, double p1, double M2, double phi,
		     double *scatMomCal, double *scatEkinCal, double *scatThetaCM);
bool calcPiMomFromLambda(double M0, double M1, double p1, double M2, double cost, double *pCal);
bool calc2BodyInelastic(double m1, double p1, double m2, double m3, double m4, 
			double theta, double *pCal1, double *pCal2, double *ThetaCM);

double calcKMomFromTheta(double bmom, double cost, double p0);
ThreeVector gaiseki(ThreeVector vec1, ThreeVector vec2);
void calcThetaPhi(ThreeVector vec, double *theta, double *phi);

const double MaxChiSqr = 100.;
const double MaxChiSqrSks = 20.;
const double PIni = 1.0;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

#ifndef MaxHits 
#define MaxHits 30
#endif
#ifndef MaxHits2 
#define MaxHits2 60
#endif


struct PartVec {
  ThreeVector pos0;
  ThreeVector mom;
};

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

  // Tracking
  int ntSks;
  double p[MaxHits];

  int ntSksPart;
  double pSksPart[MaxHits];
  double m2[MaxHits];

  int nK;
  double pK[MaxHits];
  double u0SdcIn;
  double v0SdcIn;

  double resSftU;
  double resSftV;

  double MissMass[MaxHits];

  int nSigma;
  double pKCal[MaxHits];
  double theta[MaxHits];
  double pSigmaCal[MaxHits];
  double Vertex_x[MaxHits];
  double Vertex_y[MaxHits];
  double Vertex_z[MaxHits];

  int nPi_CFT;
  int nK_CFT;
  int nP_CFT;
  double EkinP;
  double EkinCorP;

  double vertexDecayPi[3];
  double cdistDecayPi;
  double thetaDecayPi;
  double decayPiMomCal;	     
  double decayPiMomVecCal[3];
  double decayNMomCal;	     
  double decayNMomVecCal[3];


  double vertexScat[3];
  double hypBeamVec[3];
  double cdistScat;
  double thetaScat;
  double scatMomCal;
  double scatEkinCal;
  double thetaScatCM;

  double thetaDecayPi2;
  double decayPiMomCal2;	     

  double vertexNpScat[3];
  double cdistNpScat;
  double thetaNpScat;
  double scatNpMomCal;
  double scatNpEkinCal;
  double thetaNpScatCM;

  double cdistLambdaDecay;
  double vertexLambdaDecay[3];
  double thetaLambdaDecay;
  double momPiFromLambda;
  double momProtonFromLambda;
  double momLambda;
  double momVecPiFromLambda[3];
  double momVecProtonFromLambda[3];
  double momVecLambda[3];
  double cdistLambdaNConv;
  double thetaLambdaNConv;
  double vertexLambdaNConv[3];
  double momCalLambda;
  double momCalLambda2;
  double thetaCMLambdaNConv;

  double vertexPiPScat[3];
  double cdistPiPScat;

  double piBeamMom;
  double thetaPiPScat;
  double scatPiMomCal;
  double scatPiMomCalFromE;

  int ntCFT;
  double xDirCFT[MaxHits2];
  double yDirCFT[MaxHits2];
  double zDirCFT[MaxHits2];
  double xPos0CFT[MaxHits2];
  double yPos0CFT[MaxHits2];
  double zPos0CFT[MaxHits2];
  double zBGOCFT[MaxHits2];

  double BGO_Edep[MaxHits2];
  double TotalEdep[MaxHits2];
  double CFT_TotalEdep[MaxHits2];
  double CFT_NormTotalEdep[MaxHits2];  
  double PiV_Edep[MaxHits2];

  double xDirCFT_P[MaxHits2];
  double yDirCFT_P[MaxHits2];
  double zDirCFT_P[MaxHits2];
  double xPos0CFT_P[MaxHits2];
  double yPos0CFT_P[MaxHits2];
  double zPos0CFT_P[MaxHits2];
  double zBGOCFT_P[MaxHits2];
  double CFTVtx_x_P[MaxHits2];
  double CFTVtx_y_P[MaxHits2];
  double CFTVtx_z_P[MaxHits2];

  double BGO_Edep_P[MaxHits2];
  double TotalEdep_P[MaxHits2];
  double CFT_TotalEdep_P[MaxHits2];
  double CFT_NormTotalEdep_P[MaxHits2];  
  double PiV_Edep_P[MaxHits2];

  double xDirCFT_Pi[MaxHits2];
  double yDirCFT_Pi[MaxHits2];
  double zDirCFT_Pi[MaxHits2];
  double xPos0CFT_Pi[MaxHits2];
  double yPos0CFT_Pi[MaxHits2];
  double zPos0CFT_Pi[MaxHits2];
  double zBGOCFT_Pi[MaxHits2];
  double CFTVtx_x_Pi[MaxHits2];
  double CFTVtx_y_Pi[MaxHits2];
  double CFTVtx_z_Pi[MaxHits2];

  double u0CFT_Pi[MaxHits2];
  double v0CFT_Pi[MaxHits2];

  double BGO_Edep_Pi[MaxHits2];
  double TotalEdep_Pi[MaxHits2];
  double CFT_TotalEdep_Pi[MaxHits2];
  double CFT_NormTotalEdep_Pi[MaxHits2];  
  double PiV_Edep_Pi[MaxHits2];

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


int main( int argc, char **argv )
{
  if( argc != kArgc ) {
    std::cerr << "Usage: " << argv[kArgProcessName] 
	      << " [confFile]  [dataFile] [RootFile]"
	      << std::endl;
    exit(-1);
  }

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  const std::string& inFile   = arg[kArgInFile];
  const std::string& rootFile = arg[kArgOutRootFile];

  ConfMan* gconfManager = new ConfMan(confFile);
  if (!gconfManager->Initialize()) {
    std::cerr << "#E gconfManager->Initialize failed ! in Main()" << std::endl;
  }

  if (gconfManager->GetEvDispCftFlag()) {
    /*
    char RunNum_char[10];
    unsigned int loc = inFile.find("0");
    inFile.copy(RunNum_char, 5, loc);
    int RunNum=atoi(RunNum_char);
    */
    int RunNum=1;
    std::cout << "RunNumber : " << RunNum << std::endl;
    std::cout << "Create Event Display" << std::endl;
    gconfManager->InitializeEvDispCFT(RunNum);
  }

  
  // Random seed
  gRandom->SetSeed(time(NULL));

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

  bool FlagEvDisp=false;
  FlagEvDisp=ConfMan::GetConfManager()->GetEvDispCftFlag();


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
  if (FlagEvDisp) {
    std::cout << "decayFlag = " << event.decayFlag 
	      << ", scatFlag =  " << event.scatFlag 
	      << ", NNscatFlag =  " << event.NNscatFlag 
	      << ", PiNscatFlag = " << event.PiNscatFlag << std::endl;

    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    evDisp.DrawVertex0(event.primaryVertex[0], event.primaryVertex[1], event.primaryVertex[2]);

  std::cout << "primaryVertex = ( " 
	    << event.primaryVertex[0] << ", " 
	    << event.primaryVertex[1] << ", " 
	    << event.primaryVertex[2] << ")"  << std::endl;
  }

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

      if (FlagEvDisp) {
	std::cout << "seg : " << segment << ", time : " << time << ", dE : "
		  << edep << ", pid : " << pid << std::endl;
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
    } else if  (layer == 10) {
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 0, 0, edep);
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 1, 0, time);
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 0, 1, edep);
      rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 1, 1, time);
    } else if (layer == 70) {
      int Layer = segment/NumOfSegPiV;
      int Segment = segment%NumOfSegPiV;

      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 0, 0, edep);
      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 1, 0, time);
      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 0, 1, edep);
      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 1, 1, time);

      if (FlagEvDisp) {
	std::cout << "PiV seg : " << segment << ", time : " << time << ", dE : "
		  << edep  << std::endl;
      }

    }
  }

  HodoAnalyzer *hodoAna = new HodoAnalyzer;
  hodoAna->DecodeCFTHits(rawData);

  if (FlagEvDisp) {
    for(int l = 0; l<NumOfLayersCFT; ++l){
      int nhit = hodoAna->GetNHitsCFT(l);

      for(int i = 0; i<nhit; ++i){
	const CFTFiberHit* hit = hodoAna->GetHitCFT(l, i);
	int seg = hit->SegmentId();
	double pe_hi = hit->GetAHiGain();

	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.ShowHitFiber(l, seg, pe_hi*100);
      }
    }
  }

  hodoAna->DecodeBGOHits(rawData);

  if (FlagEvDisp) {
    int nhit = hodoAna->GetNHitsBGO();
    for (int n=0; n<nhit; n++) {
      const Hodo2Hit* hit = hodoAna->GetHitBGO(n);
      int segment = hit->SegmentId();
      double dE   = hit->DeltaE();

      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      evDisp.ShowHitBGO(segment, dE);
    }
  }

  hodoAna->DecodeTOFHits(rawData);
  hodoAna->DecodePiVHits(rawData);

  if (FlagEvDisp) {
    int nhit = hodoAna->GetNHitsPiV();
    for (int n=0; n<nhit; n++) {
      const Hodo2Hit* hit = hodoAna->GetHitPiV(n);
      int segment = hit->SegmentId();
      int layer   = hit->PlaneId();
      double dE   = hit->DeltaE();
      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      evDisp.ShowHitPiV(layer, segment, dE);
    }
  }


  int ncTof=hodoAna->GetNClustersTOF();

  DCAnalyzer *DCAna = new DCAnalyzer;
  DCAna->DecodeSdcInHits(rawData);
  DCAna->DecodeSdcOutHits(rawData);
  DCAna->DecodeBcOutHits(rawData);

  DCAna->TrackSearchSdcInFiber();
  int ntSdcIn = DCAna->GetNtracksSdcIn();
  HF1(10, ntSdcIn);

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
  }

  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );

  DCAna->TrackSearchSdcOut();
  int ntSdcOut = DCAna->GetNtracksSdcOut();

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

  }

  //DCAna->TrackSearchSks();
  DCAna->TrackSearchSksTmp(event.momScatMeson);
  int ntSks = DCAna->GetNTracksSks();
  event.ntSks = ntSks;

  std::vector <SksParticle *> SksPartCont;

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

    if (chisqr>MaxChiSqrSks)
      continue;

    DCLocalTrack *tout = tp->GetLocalTrackOut();
    double xtof = tout->GetFTOFX();
    for (int j=0; j<ncTof; ++j) {
      HodoCluster *clTof=hodoAna->GetClusterTOF(j);
      if (!clTof || !clTof->GoodForAnalysis() ) continue;

      double meanSeg = clTof->MeanSeg();
      double tofSegPos = DCGeomMan::GetInstance().calcWirePosition(IdTof, meanSeg);

      double diffX = xtof-tofSegPos;
      HF1(64, diffX);
      if (fabs(diffX)<60.) {
	SksParticle *SksPart = new SksParticle(tp, clTof);
	SksPartCont.push_back(SksPart);

	double m2 = MassSquare(p, pathL, clTof->CMeanTime());
	HF1(65, m2);

	if (m2>0)
	  HF1(66, sqrt(m2));
	
      }
    }
  }

  int ntSksPart = SksPartCont.size();
  event.ntSksPart = ntSksPart;
  for (int i=0; i<ntSksPart; i++) {
    event.pSksPart[i] = SksPartCont[i]->Momentum().mag();
    event.m2[i] = SksPartCont[i]->MassSquare();
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
    event.resSftU = resU;

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
    event.resSftV = resV;
  }

  int nBeam=BmomCont.size();

  if (nBeam != 1) {
    tree->Fill();

    for_each(SksPartCont.begin(), SksPartCont.end(), DeleteObject());
    
    delete hodoAna;
    delete DCAna;
    delete rawData;
    
    if (FlagEvDisp) {
      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      /*
      //if (1) {
      if (event.scatFlag==1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1&&nPi_CFT==1) {
	evDisp.UpdateCanvas();
	evDisp.get_command();
      }
      */
      evDisp.EndOfEvent();
    }

    return true;
  }

  ThreeVector bMom = BmomCont[0];
  ThreeVector bPos = BposCont[0];

  std::vector <PartVec> SigmaVecCont;
  std::vector <SksParticle *> KaonPartCont;

  int nK=0;

  // Missing Mass Analysis
  int nSksPart = SksPartCont.size();
  int index=0;
  for (int i=0; i<nSksPart; i++) {
    SksParticle *SksPart = SksPartCont[i];
    double m2 = SksPart->MassSquare();
    if (m2>MinOfMassSquareK && m2<MaxOfMassSquareK) {
      ThreeVector priMom = (-1.)*SksPart->Momentum();

      LorentzVector LvBeam(  bMom, sqrt(PionMass*PionMass+bMom.mag2()) );
      LorentzVector LvTgt( 0., 0., 0., ProtonMass );
      LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
      LorentzVector LvRc = LvBeam+LvTgt-LvScat;
      double mismass = LvRc.mag();
      HF1(67, mismass);
      event.MissMass[nK] = mismass;
      event.pK[nK] = priMom.mag();

      nK++;

      if (mismass>MinOfMassSigma && mismass<MaxOfMassSigma) {
	KaonPartCont.push_back(SksPart);

	ThreeVector pPos=SksPart->Position();
	ThreeVector vtx = VertexPoint(bPos, pPos, bMom, priMom);
	event.Vertex_x[index] = vtx.x();
	event.Vertex_y[index] = vtx.y();
	event.Vertex_z[index] = vtx.z();

	double cost = bMom*priMom/(bMom.mag()*priMom.mag());
	event.theta[index] = acos(cost)*Rad2Deg;


	double p0 = priMom.mag();
	double p_Cor, E_Cor;    
	CorrElossOut(&p_Cor, &E_Cor, p0, KAON, priMom/p0, vtx);
	ThreeVector KaonMomCor = priMom* p_Cor/p0;
	ThreeVector pSigmaCor = bMom-KaonMomCor;

	double calMom = calcKMomFromTheta(event.beammom, cost, priMom.mag());
	event.pKCal[index] = calMom;
	ThreeVector KaonMomCal = priMom*(calMom/priMom.mag());
	ThreeVector pSigmaCal = bMom-KaonMomCal;
	event.pSigmaCal[index] = pSigmaCal.mag();

	PartVec sigmaVec;
	sigmaVec.pos0 = vtx;
	sigmaVec.mom  = pSigmaCal; // original
	//sigmaVec.mom  = pSigmaCal*event.momHypBeam/pSigmaCal.mag();
	//sigmaVec.mom  = ThreeVector(event.momVectorHypBeam[0],event.momVectorHypBeam[1],event.momVectorHypBeam[2]);
	//sigmaVec.mom  = pSigmaCor;
	SigmaVecCont.push_back(sigmaVec);
	
	index++;
      }
    }
  }

  event.nK = nK;

  int nSigma = SigmaVecCont.size();
  event.nSigma = nSigma;

  // CFT Analysis
  // detailed analysis is later
  DCAna->DecodeCFTHits(rawData);
  DCAna->TrackSearchCFT();
  std::vector <CFTParticle *> CFTPartCont;
  std::vector <CFTParticle *> CFTProtonCont;
  std::vector <CFTParticle *> CFTPionCont;
  std::vector <CFTParticle *> CFTKaonCont;

  std::vector <ThreeVector> CFTVtxCont;

  int nhitBGO = hodoAna->GetNHitsBGO();
  int nhitPiV = hodoAna->GetNHitsPiV();

  int ntCFT = DCAna->GetNtracksCFT();
  event.ntCFT = ntCFT;
  HF1(400, ntCFT);

  if (nSigma != 1) {
    // CFT analysis when Sigma is not identified.
    for (int it=0; it<ntCFT; it++) {
      CFTLocalTrack *tp = DCAna->GetTrackCFT(it);

      int xyFitFlag = tp->GetXYFitFlag();
      double Axy = tp->GetAxy();
      double Bxy = tp->GetBxy();
      ThreeVector Pos0 = tp->GetPos0();
      ThreeVector Dir = tp->GetDir();

      CFTParticle * CFTPart = new CFTParticle(tp);
      double dE_BGO=0.;

      for (int n=0; n<nhitBGO; n++) {
	Hodo2Hit* hit = hodoAna->GetHitBGO(n);
	int segment = hit->SegmentId();
	double time = hit->CMeanTime();
	double dE   = hit->DeltaE();

	double x, y;
	hit->BGOPos(&x, &y);
	double dist=1000.;
	if (xyFitFlag==0) {
	  dist = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);
	} else if (xyFitFlag==1) {
	  dist = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);
	}

	double u=Dir.x(), v=Dir.y();
	double x0=Pos0.x(), y0=Pos0.y();
	double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);

	if (t>=0) {      
	  HF1(421, dist);

	  if (fabs(dist)<25 && time>0 && time<5) {
	    CFTPart->AddBGOHit(hit);
	  }
	}
      }    
      
      double dE_PiV=0.;
      for (int n=0; n<nhitPiV; n++) {
	Hodo2Hit* hit = hodoAna->GetHitPiV(n);
	int segment = hit->SegmentId();
	int layer   = hit->PlaneId();	
	double time = hit->CMeanTime();
	double dE   = hit->DeltaE();

	double x, y;
	hit->PiVPos(layer, &x, &y);
	double dist=1000.;
	if (xyFitFlag==0) {
	  dist = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);
	} else if (xyFitFlag==1) {
	  dist = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);
	}
	
	HF1(450, dist);
	
	double u=Dir.x(), v=Dir.y();
	double x0=Pos0.x(), y0=Pos0.y();
	double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);
	
	if (t>=0) {
	  if (fabs(dist)<30 && time>0 && time<5) {
	    dE_PiV += dE;
	    CFTPart->AddPiVHit(hit);
	  }
	}
      }    
      
      CFTPart->Calculate();
      CFTPartCont.push_back(CFTPart);

      if (it<MaxHits2) {
	event.xDirCFT[it] = Dir.x();
	event.yDirCFT[it] = Dir.y();
	event.xDirCFT[it] = Dir.z();
	event.xPos0CFT[it] = Pos0.x();
	event.yPos0CFT[it] = Pos0.y();
	event.zPos0CFT[it] = Pos0.z();
      }


      if (CFTPart->GetMass() > 0.9)
	CFTProtonCont.push_back(CFTPart);
      else if (CFTPart->GetMass() > 0.0 && CFTPart->GetMass() < 0.2)
	CFTPionCont.push_back(CFTPart);

    }
    event.nP_CFT = CFTProtonCont.size();
    event.nPi_CFT = CFTPionCont.size();

    tree->Fill();

    for_each(SksPartCont.begin(), SksPartCont.end(), DeleteObject());
    //for_each(KaonPartCont.begin(), KaonPartCont.end(), DeleteObject());
    for_each(CFTPartCont.begin(), CFTPartCont.end(), DeleteObject());
    
    delete hodoAna;
    delete DCAna;
    delete rawData;
    
    if (FlagEvDisp) {
      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      /*
      if (1) {
	//if (event.scatFlag==-1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1&&nPi_CFT==1) {
	evDisp.UpdateCanvas();
	evDisp.get_command();
      }
      */
      evDisp.EndOfEvent();
    }

    return true;
  }

  if (KaonPartCont.size() != 1) {
    std::cout << "KaonPart size is NOT 1" << std::endl;

    tree->Fill();

    for_each(SksPartCont.begin(), SksPartCont.end(), DeleteObject());
    //for_each(KaonPartCont.begin(), KaonPartCont.end(), DeleteObject());
    
    delete hodoAna;
    delete DCAna;
    delete rawData;

    return true;
  }

  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    evDisp.DrawVertex(event.Vertex_x[0], event.Vertex_y[0], event.Vertex_z[0]);

    SksParticle *KaonPart = KaonPartCont[0];
    double u0 = KaonPart->GetTrack()->GetLocalTrackIn()->GetU0();
    double v0 = KaonPart->GetTrack()->GetLocalTrackIn()->GetV0();
    double x0 = KaonPart->GetTrack()->GetLocalTrackIn()->GetX0();
    double y0 = KaonPart->GetTrack()->GetLocalTrackIn()->GetY0();
    double z0 = 0;

    double z1 = 350.;
    double x1 = u0*z1+x0;
    double y1 = v0*z1+y0;
    bool flagK = true;

    evDisp.DrawTrackInXYPlane(x0, y0, x1, y1, flagK);
    evDisp.DrawTrackInZXPlane(z0, x0, z1, x1, flagK);
    evDisp.DrawTrackInZYPlane(z0, y0, z1, y1, flagK);
  }


  SksParticle *KaonPart = KaonPartCont[0];
  event.u0SdcIn = KaonPart->GetTrack()->GetLocalTrackIn()->GetU0();
  event.v0SdcIn = KaonPart->GetTrack()->GetLocalTrackIn()->GetV0();

  PartVec SigmaVec = SigmaVecCont[0];
  event.hypBeamVec[0] = SigmaVec.mom.x();
  event.hypBeamVec[1] = SigmaVec.mom.y();
  event.hypBeamVec[2] = SigmaVec.mom.z();

  // CFT Analysis
  for (int it=0; it<ntCFT; it++) {
    CFTLocalTrack *tp = DCAna->GetTrackCFT(it);
    int nhXY = tp->GetNHit();
    int nhZ = tp->GetNHitU();
    double chisqrXY = tp->GetChiSquareXY();
    double chisqrZ = tp->GetChiSquareZ();
    double x0 = tp->GetX0(), u0 = tp->GetU0();
    double y0 = tp->GetY0(), v0 = tp->GetV0();

    HF1(402, nhXY);
    HF1(403, nhZ);
    HF1(404, chisqrXY);
    HF1(405, chisqrZ);
    HF1(406, x0);    HF1(407, y0);
    HF1(408, u0);    HF1(409, v0);
    HF2(410, x0, y0);HF1(411, u0, v0);

    double totalDE = tp->TotalDEHiGain();
    double maxDE = tp->MaxDEHiGain();
    double pathlength = tp->GetTotalPathLength();
    double normalizedTotalDE = tp->NormalizedTotalDEHiGain();
    double normalizedMaxDE = tp->NormalizedMaxDEHiGain();
    
    HF1(412, totalDE);
    HF1(413, maxDE);
    HF1(414, pathlength);
    HF1(415, normalizedTotalDE);
    HF1(416, normalizedMaxDE);
    HF2(417, pathlength, totalDE);
    HF2(418, pathlength, maxDE);
    HF2(419, pathlength, normalizedTotalDE);
    HF2(420, pathlength, normalizedMaxDE);

    int xyFitFlag = tp->GetXYFitFlag();
    double Axy = tp->GetAxy();
    double Bxy = tp->GetBxy();
    ThreeVector Pos0 = tp->GetPos0();
    ThreeVector Dir = tp->GetDir();

    double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());
    double B=(Dir.x()*Pos0.x()+Dir.y()*Pos0.y());
    double C=Pos0.x()*Pos0.x()+Pos0.y()*Pos0.y()-RadiusOfBGOSurface*RadiusOfBGOSurface;
    double t1, t2;
    double zBGO=-999.;
    double trackPosBGO[3]={-999., -999., -999};
    if (B*B-A*C>=0) {
      t1 = (-B+sqrt(B*B-A*C))/A;
      t2 = (-B-sqrt(B*B-A*C))/A;

      if (t1>=0 && t2<=0) {
	zBGO = Dir.z()*t1 + Pos0.z();
	trackPosBGO[0] = Dir.x()*t1 + Pos0.x();
	trackPosBGO[1] = Dir.y()*t1 + Pos0.y();
	trackPosBGO[2] = Dir.z()*t1 + Pos0.z();
      } else if (t2>=0 && t1<=0) {
	zBGO = Dir.z()*t2 + Pos0.z();
	trackPosBGO[0] = Dir.x()*t2 + Pos0.x();
	trackPosBGO[1] = Dir.y()*t2 + Pos0.y();
	trackPosBGO[2] = Dir.z()*t2 + Pos0.z();
      } else {
	std::cout << "t1 : " << t1 << ", t2 : " << t2 << std::endl;
      }
	
    }

    if (xyFitFlag==0) {
      double x1=-100., x2=100.;
      double y1 = Axy*x1+Bxy;
      double y2 = Axy*x2+Bxy;
      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.DrawTrackInXYPlane(Pos0.x(), Pos0.y(), trackPosBGO[0], trackPosBGO[1]);
      }
    } else if (xyFitFlag==1) {
      double y1=-100., y2=100.;
      double x1 = Axy*y1+Bxy;
      double x2 = Axy*y2+Bxy;
      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.DrawTrackInXYPlane(Pos0.x(), Pos0.y(), trackPosBGO[0], trackPosBGO[1]);
      }
    }

    if (tp->GetZTrackFlag()==0 || tp->GetZTrackFlag()==1) {
      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.DrawTrackInZXPlane(Pos0.z(), Pos0.x(), trackPosBGO[2], trackPosBGO[0]);
	evDisp.DrawTrackInZYPlane(Pos0.z(), Pos0.y(), trackPosBGO[2], trackPosBGO[1]);
      }
    }

    CFTParticle * CFTPart = new CFTParticle(tp);

    double dE_BGO=0.;

    for (int n=0; n<nhitBGO; n++) {
      Hodo2Hit* hit = hodoAna->GetHitBGO(n);
      int segment = hit->SegmentId();
      double time = hit->CMeanTime();
      double dE   = hit->DeltaE();

      double x, y;
      hit->BGOPos(&x, &y);
      double dist=1000.;
      if (xyFitFlag==0) {
	dist = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);
      } else if (xyFitFlag==1) {
	dist = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);
      }

      double u=Dir.x(), v=Dir.y();
      double x0=Pos0.x(), y0=Pos0.y();
      double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);

      if (t>=0) {      
	HF1(421, dist);

	if (fabs(dist)<25 && time>0 && time<5) {
	  //dE_BGO += dE;
	  CFTPart->AddBGOHit(hit);
	}
      }
    }    

    double dE_PiV=0.;
    for (int n=0; n<nhitPiV; n++) {
      Hodo2Hit* hit = hodoAna->GetHitPiV(n);
      int segment = hit->SegmentId();
      int layer   = hit->PlaneId();
      double time = hit->CMeanTime();
      double dE   = hit->DeltaE();

      double x, y;
      hit->PiVPos(layer, &x, &y);
      double dist=1000.;
      if (xyFitFlag==0) {
	dist = (Axy*x-y+Bxy)/sqrt(Axy*Axy+1.*1.);
      } else if (xyFitFlag==1) {
	dist = (Axy*y-x+Bxy)/sqrt(Axy*Axy+1.*1.);
      }

      HF1(450, dist);

      double u=Dir.x(), v=Dir.y();
      double x0=Pos0.x(), y0=Pos0.y();
      double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);

      if (t>=0) {
	if (fabs(dist)<30 && time>0 && time<5) {
	  dE_PiV += dE;
	  CFTPart->AddPiVHit(hit);
	}
      }
    }    

    CFTPart->Calculate();
    dE_BGO = CFTPart->GetBGO_E();

    ThreeVector pos0 = CFTPart->GetPos0();
    ThreeVector dir  = CFTPart->GetDir();
    /*
    std::cout << "Pos0 ( " << pos0.x() << ", " << pos0.y() << ", " << pos0.z() << "), "
	      << " Dir ( " << dir.x()/dir.mag() << ", " << dir.y()/dir.mag() 
	      << ", " << dir.z()/dir.mag() << ") " << std::endl;
    */

    CFTPartCont.push_back(CFTPart);

    HF2(422, dE_BGO, totalDE);
    HF2(423, dE_BGO, maxDE);
    HF2(424, dE_BGO, normalizedTotalDE);
    HF2(425, dE_BGO, normalizedMaxDE);

    HF2(426, dE_BGO+totalDE, totalDE);
    HF2(427, dE_BGO+maxDE, maxDE);
    HF2(428, dE_BGO+totalDE, normalizedTotalDE);
    HF2(429, dE_BGO+maxDE, normalizedMaxDE);

    if (it<MaxHits2) {
      event.xDirCFT[it] = Dir.x();
      event.yDirCFT[it] = Dir.y();
      event.xDirCFT[it] = Dir.z();
      event.xPos0CFT[it] = Pos0.x();
      event.yPos0CFT[it] = Pos0.y();
      event.zPos0CFT[it] = Pos0.z();

      event.zBGOCFT[it] = zBGO;

      event.CFT_TotalEdep[it] = totalDE;
      event.CFT_NormTotalEdep[it] = normalizedTotalDE;
      event.BGO_Edep[it] = dE_BGO;
      event.PiV_Edep[it] = dE_PiV;
      event.TotalEdep[it] = dE_BGO+totalDE;
    }

    // Vertex search for next tracking
    double cdist;
    ThreeVector CFTVtx = VertexPoint3D(bPos,Pos0, bMom, Dir, cdist);
    CFTVtxCont.push_back(CFTVtx);
    CFTPart->SetCFTVtx(CFTVtx);


    DCLocalTrack *In = KaonPart->GetTrack()->GetLocalTrackIn();
    double kaonU0 = In->GetU0();
    double kaonV0 = In->GetV0();

    HF1(430, u0-kaonU0);
    HF1(431, v0-kaonV0);
    HF2(432, u0-kaonU0, v0-kaonV0);

    ThreeVector kaonPos0(In->GetX0(), In->GetY0(), 0);
    ThreeVector kaonDir(kaonU0, kaonV0, 1);
    double cdistKaonVtx;
    ThreeVector vtxKaon = VertexPoint3D(Pos0, kaonPos0, Dir, kaonDir, cdistKaonVtx);
    HF1(434, cdistKaonVtx);
    HF1(435, vtxKaon.z());
    HF2(436, vtxKaon.x(), vtxKaon.y());

    if (u0-kaonU0>MinOfDeltaU_CFT && u0-kaonU0<MaxOfDeltaU_CFT
	&& v0-kaonV0>MinOfDeltaV_CFT && v0-kaonV0<MaxOfDeltaV_CFT)
      CFTKaonCont.push_back(CFTPart);
    else if  (u0-kaonU0>-0.20 && u0-kaonU0<0.20
	      && v0-kaonV0>-0.20 && v0-kaonV0<0.20 && cdistKaonVtx < 5 && vtxKaon.z()>150 && vtxKaon.z()<370 && vtxKaon.x()>70) 
      CFTKaonCont.push_back(CFTPart);
    else if (CFTPart->GetMass() > 0.9)
      CFTProtonCont.push_back(CFTPart);
    else if (CFTPart->GetMass() > 0.0 && CFTPart->GetMass() < 0.2)
      CFTPionCont.push_back(CFTPart);


  }

  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    int nCFTVtx = CFTVtxCont.size();
    for (int iv=0; iv<nCFTVtx; iv++) {
      ThreeVector CFTVtx = CFTVtxCont[iv];
      evDisp.DrawCFTVertex(CFTVtx.x(), CFTVtx.y(), CFTVtx.z());
    }
  }
  
  int ntCFT_2nd = 0;

  if (CFTProtonCont.size() == 0 || CFTPionCont.size() == 0 ) {
    if (CFTVtxCont.size()>0) {
      DCAna->TrackSearchCFT_2nd(CFTVtxCont);
      ntCFT_2nd = DCAna->GetNtracksCFT2nd();
      //std::cout << "ntCFT_2nd : " << ntCFT_2nd << std::endl;
      event.ntCFT = ntCFT+ntCFT_2nd;
      for (int it=0; it<ntCFT_2nd; it++) {
	CFTLocalTrack *tp = DCAna->GetTrackCFT2nd(it);
	int nhXY = tp->GetNHit();
	int nhZ = tp->GetNHitU();
	double chisqrXY = tp->GetChiSquareXY();
	double chisqrZ = tp->GetChiSquareZ();
	double x0 = tp->GetX0(), u0 = tp->GetU0();
	double y0 = tp->GetY0(), v0 = tp->GetV0();
	
	HF1(502, nhXY);
	HF1(503, nhZ);
	HF1(504, chisqrXY);
	HF1(505, chisqrZ);
	HF1(506, x0);    HF1(507, y0);
	HF1(508, u0);    HF1(509, v0);
	HF2(510, x0, y0);HF1(511, u0, v0);
	
	double totalDE = tp->TotalDEHiGain();
	double maxDE = tp->MaxDEHiGain();
	double pathlength = tp->GetTotalPathLength();
	double normalizedTotalDE = tp->NormalizedTotalDEHiGain();
	double normalizedMaxDE = tp->NormalizedMaxDEHiGain();
	
	HF1(512, totalDE);
	HF1(513, maxDE);
	HF1(514, pathlength);
	HF1(515, normalizedTotalDE);
	HF1(516, normalizedMaxDE);
	HF2(517, pathlength, totalDE);
	HF2(518, pathlength, maxDE);
	HF2(519, pathlength, normalizedTotalDE);
	HF2(520, pathlength, normalizedMaxDE);
	
	int xyFitFlag = tp->GetXYFitFlag();
	double Axy = tp->GetAxy();
	double Bxy = tp->GetBxy();
	ThreeVector Pos0 = tp->GetPos0();
	ThreeVector Dir = tp->GetDir();
	
	double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());
	double B=(Dir.x()*Pos0.x()+Dir.y()*Pos0.y());
	double C=Pos0.x()*Pos0.x()+Pos0.y()*Pos0.y()-RadiusOfBGOSurface*RadiusOfBGOSurface;
	double t1, t2;
	double zBGO=-999.;
	double trackPosBGO[3]={-999., -999., -999};
	if (B*B-A*C>=0) {
	  t1 = (-B+sqrt(B*B-A*C))/A;
	  t2 = (-B-sqrt(B*B-A*C))/A;
	  
	  if (t1>=0 && t2<=0) {
	    zBGO = Dir.z()*t1 + Pos0.z();
	    trackPosBGO[0] = Dir.x()*t1 + Pos0.x();
	    trackPosBGO[1] = Dir.y()*t1 + Pos0.y();
	    trackPosBGO[2] = Dir.z()*t1 + Pos0.z();
	  } else if (t2>=0 && t1<=0) {
	    zBGO = Dir.z()*t2 + Pos0.z();
	    trackPosBGO[0] = Dir.x()*t2 + Pos0.x();
	    trackPosBGO[1] = Dir.y()*t2 + Pos0.y();
	    trackPosBGO[2] = Dir.z()*t2 + Pos0.z();
	  } else {
	    std::cout << "t1 : " << t1 << ", t2 : " << t2 << std::endl;
	  }
	  
	}
	
	if (xyFitFlag==0) {
	  double x1=-100., x2=100.;
	  double y1 = Axy*x1+Bxy;
	  double y2 = Axy*x2+Bxy;
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.DrawTrackInXYPlane(Pos0.x(), Pos0.y(), trackPosBGO[0], trackPosBGO[1]);
	  }
	} else if (xyFitFlag==1) {
	  double y1=-100., y2=100.;
	  double x1 = Axy*y1+Bxy;
	  double x2 = Axy*y2+Bxy;
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.DrawTrackInXYPlane(Pos0.x(), Pos0.y(), trackPosBGO[0], trackPosBGO[1]);
	  }
	}
	
	if (tp->GetZTrackFlag()==0 || tp->GetZTrackFlag()==1) {
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.DrawTrackInZXPlane(Pos0.z(), Pos0.x(), trackPosBGO[2], trackPosBGO[0]);
	    evDisp.DrawTrackInZYPlane(Pos0.z(), Pos0.y(), trackPosBGO[2], trackPosBGO[1]);
	  }
	}
	
	CFTParticle * CFTPart = new CFTParticle(tp);
	CFTPart->Calculate();
	
	ThreeVector pos0 = CFTPart->GetPos0();
	ThreeVector dir  = CFTPart->GetDir();
	/*
	  std::cout << "Pos0 ( " << pos0.x() << ", " << pos0.y() << ", " << pos0.z() << "), "
	  << " Dir ( " << dir.x()/dir.mag() << ", " << dir.y()/dir.mag() 
	  << ", " << dir.z()/dir.mag() << ") " << std::endl;
	*/
	
	CFTPartCont.push_back(CFTPart);
	
	/*
	  if (it<MaxHits2) {
	  event.xDirCFT[it] = Dir.x();
	  event.yDirCFT[it] = Dir.y();
	  event.xDirCFT[it] = Dir.z();
	  event.xPos0CFT[it] = Pos0.x();
	  event.yPos0CFT[it] = Pos0.y();
	  event.zPos0CFT[it] = Pos0.z();
	  
	  event.zBGOCFT[it] = zBGO;
	
	  event.CFT_TotalEdep[it] = totalDE;
	  event.CFT_NormTotalEdep[it] = normalizedTotalDE;
	  event.BGO_Edep[it] = dE_BGO;
	  event.PiV_Edep[it] = dE_PiV;
	  event.TotalEdep[it] = dE_BGO+totalDE;
	  }
	*/

	double cdist;
	ThreeVector CFTVtx = VertexPoint3D(bPos,Pos0, bMom, Dir, cdist);
	CFTVtxCont.push_back(CFTVtx);
	CFTPart->SetCFTVtx(CFTVtx);
	
	
	DCLocalTrack *In = KaonPart->GetTrack()->GetLocalTrackIn();
	double kaonU0 = In->GetU0();
	double kaonV0 = In->GetV0();
      
	HF1(530, u0-kaonU0);
	HF1(531, v0-kaonV0);
	HF2(532, u0-kaonU0, v0-kaonV0);
	
	ThreeVector kaonPos0(In->GetX0(), In->GetY0(), 0);
	ThreeVector kaonDir(kaonU0, kaonV0, 1);
	double cdistKaonVtx;
	ThreeVector vtxKaon = VertexPoint3D(Pos0, kaonPos0, Dir, kaonDir, cdistKaonVtx);
	HF1(534, cdistKaonVtx);
	HF1(535, vtxKaon.z());
	HF2(536, vtxKaon.x(), vtxKaon.y());

	if (u0-kaonU0>MinOfDeltaU_CFT && u0-kaonU0<MaxOfDeltaU_CFT
	    && v0-kaonV0>MinOfDeltaV_CFT && v0-kaonV0<MaxOfDeltaV_CFT)
	  CFTKaonCont.push_back(CFTPart);
	else if  (u0-kaonU0>-0.20 && u0-kaonU0<0.20
		  && v0-kaonV0>-0.20 && v0-kaonV0<0.20 && cdistKaonVtx < 5 && vtxKaon.z()>150 && vtxKaon.z()<370 && vtxKaon.x()>70) 
	  CFTKaonCont.push_back(CFTPart);
	else if (CFTPart->GetMass() > 0.9 && CFTProtonCont.size() == 0)
	  CFTProtonCont.push_back(CFTPart);
	else if (CFTPart->GetMass() > 0.0 && CFTPart->GetMass() < 0.2 && CFTPionCont.size() == 0)
	  CFTPionCont.push_back(CFTPart);
      }
    }
  }
  int nK_CFT = CFTKaonCont.size();
  int nP_CFT = CFTProtonCont.size();
  int nPi_CFT = CFTPionCont.size();
  event.nK_CFT = nK_CFT;
  event.nP_CFT = nP_CFT;
  event.nPi_CFT = nPi_CFT;

  for (int it=0; it<nP_CFT; it++) {
    if (it<MaxHits2) {
      CFTParticle *ProtonCFT = CFTProtonCont[it];
      CFTLocalTrack *tp = ProtonCFT->GetTrack();
      ThreeVector Pos0 = ProtonCFT->GetPos0();
      ThreeVector Dir = ProtonCFT->GetDir();
      
      double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());
      double B=(Dir.x()*Pos0.x()+Dir.y()*Pos0.y());
      double C=Pos0.x()*Pos0.x()+Pos0.y()*Pos0.y()-RadiusOfBGOSurface*RadiusOfBGOSurface;
      double t1, t2;
      double zBGO=-999.;
      if (B*B-A*C>=0) {
	t1 = (-B+sqrt(B*B-A*C))/A;
	t2 = (-B-sqrt(B*B-A*C))/A;
	
	if (t1>=0 && t2<=0)
	  zBGO = Dir.z()*t1 + Pos0.z();
	else if (t2>=0 && t1<=0)
	  zBGO = Dir.z()*t2 + Pos0.z();
	else {
	  std::cout << "t1 : " << t1 << ", t2 : " << t2 << std::endl;
	}
      }

      event.xDirCFT_P[it] = Dir.x();
      event.yDirCFT_P[it] = Dir.y();
      event.xDirCFT_P[it] = Dir.z();
      event.xPos0CFT_P[it] = Pos0.x();
      event.yPos0CFT_P[it] = Pos0.y();
      event.zPos0CFT_P[it] = Pos0.z();

      event.zBGOCFT_P[it] = zBGO;

      event.CFT_TotalEdep_P[it] = CFTProtonCont[it]->GetFiberTotal_E();
      event.CFT_NormTotalEdep_P[it] = CFTProtonCont[it]->GetNormFiberTotal_E();
      event.BGO_Edep_P[it] = CFTProtonCont[it]->GetBGO_E();
      event.PiV_Edep_P[it] = CFTProtonCont[it]->GetPiV_E();
      event.TotalEdep_P[it] = CFTProtonCont[it]->GetTotalE();

      ThreeVector CFTVtx = ProtonCFT->GetCFTVtx();
      event.CFTVtx_x_P[it] = CFTVtx.x();
      event.CFTVtx_y_P[it] = CFTVtx.y();
      event.CFTVtx_z_P[it] = CFTVtx.z();

      int NHit = tp->GetNHit();
      for (int i=0; i<NHit; i++) {
	CFTFiberCluster *fcl = tp->GetHit(i);
	double dE = fcl->TotalPhotonNum();
	int clSize = fcl->ClusterSize();
	int l = fcl->GetTrackingLayer() -61;

	HF1(440+l, dE);
	for (int j=0; j<clSize; j++) {
	  CFTFLHit *flhit = fcl->GetHit(j);
	  int seg = flhit->PairId();
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.ShowHitFiberProton(l, seg);
	  }
	}
      }
      int NHitU = tp->GetNHitU();
      for (int i=0; i<NHitU; i++) {
	CFTFiberCluster *fcl = tp->GetHitU(i);
	int clSize = fcl->ClusterSize();
	int l = fcl->GetTrackingLayer()-61;
	double dE = fcl->TotalPhotonNum();

	HF1(440+l, dE);

	for (int j=0; j<clSize; j++) {
	  CFTFLHit *flhit = fcl->GetHit(j);
	  int seg = flhit->PairId();
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.ShowHitFiberProton(l, seg);
	  }
	}
      }
      int NHitBGO = ProtonCFT->NHitBGO();
      for (int i=0; i<NHitBGO; i++) {
	Hodo2Hit *hit = ProtonCFT->GetBGOHit(i);
	int segment = hit->SegmentId();
	if (FlagEvDisp) {
	  const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	  evDisp.ShowHitBGO_Proton(segment);
	}
      }
      int NHitPiV = ProtonCFT->NHitPiV();
      for (int i=0; i<NHitPiV; i++) {
	Hodo2Hit *hit = ProtonCFT->GetPiVHit(i);
	int segment = hit->SegmentId();
	int layer   = hit->PlaneId();
	if (FlagEvDisp) {
	  const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	  evDisp.ShowHitPiV_Proton(layer, segment);
	}
      }
    }
  }

  for (int it=0; it<nPi_CFT; it++) {
    if (it<MaxHits2) {
      CFTParticle *PiCFT = CFTPionCont[it];
      CFTLocalTrack *tp = PiCFT->GetTrack();
      ThreeVector Pos0 = PiCFT->GetPos0();
      ThreeVector Dir = PiCFT->GetDir();
      
      double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());
      double B=(Dir.x()*Pos0.x()+Dir.y()*Pos0.y());
      double C=Pos0.x()*Pos0.x()+Pos0.y()*Pos0.y()-RadiusOfBGOSurface*RadiusOfBGOSurface;
      double t1, t2;
      double zBGO=-999.;
      if (B*B-A*C>=0) {
	t1 = (-B+sqrt(B*B-A*C))/A;
	t2 = (-B-sqrt(B*B-A*C))/A;
	
	if (t1>=0 && t2<=0)
	  zBGO = Dir.z()*t1 + Pos0.z();
	else if (t2>=0 && t1<=0)
	  zBGO = Dir.z()*t2 + Pos0.z();
	else {
	  std::cout << "t1 : " << t1 << ", t2 : " << t2 << std::endl;
	}
      }
      event.xDirCFT_Pi[it] = Dir.x();
      event.yDirCFT_Pi[it] = Dir.y();
      event.xDirCFT_Pi[it] = Dir.z();
      event.xPos0CFT_Pi[it] = Pos0.x();
      event.yPos0CFT_Pi[it] = Pos0.y();
      event.zPos0CFT_Pi[it] = Pos0.z();

      event.zBGOCFT_Pi[it] = zBGO;

      event.CFT_TotalEdep_Pi[it] = CFTPionCont[it]->GetFiberTotal_E();
      event.CFT_NormTotalEdep_Pi[it] = CFTPionCont[it]->GetNormFiberTotal_E();
      event.BGO_Edep_Pi[it] = CFTPionCont[it]->GetBGO_E();
      event.PiV_Edep_Pi[it] = CFTPionCont[it]->GetPiV_E();
      event.TotalEdep_Pi[it] = CFTPionCont[it]->GetTotalE();

      ThreeVector CFTVtx = PiCFT->GetCFTVtx();
      event.CFTVtx_x_Pi[it] = CFTVtx.x();
      event.CFTVtx_y_Pi[it] = CFTVtx.y();
      event.CFTVtx_z_Pi[it] = CFTVtx.z();

      event.u0CFT_Pi[it] = tp->GetU0();
      event.v0CFT_Pi[it] = tp->GetV0();


      int NHit = tp->GetNHit();
      for (int i=0; i<NHit; i++) {
	CFTFiberCluster *fcl = tp->GetHit(i);
	int clSize = fcl->ClusterSize();
	int l = fcl->GetTrackingLayer() -61;

	double dE = fcl->TotalPhotonNum();

	HF1(450+l, dE);

	for (int j=0; j<clSize; j++) {
	  CFTFLHit *flhit = fcl->GetHit(j);
	  int seg = flhit->PairId();
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.ShowHitFiberPi(l, seg);
	  }
	}
      }
      int NHitU = tp->GetNHitU();
      for (int i=0; i<NHitU; i++) {
	CFTFiberCluster *fcl = tp->GetHitU(i);
	int clSize = fcl->ClusterSize();
	int l = fcl->GetTrackingLayer()-61;
	double dE = fcl->TotalPhotonNum();

	HF1(450+l, dE);

	for (int j=0; j<clSize; j++) {
	  CFTFLHit *flhit = fcl->GetHit(j);
	  int seg = flhit->PairId();
	  if (FlagEvDisp) {
	    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	    evDisp.ShowHitFiberPi(l, seg);
	  }
	}
      }
      int NHitBGO = PiCFT->NHitBGO();
      for (int i=0; i<NHitBGO; i++) {
	Hodo2Hit *hit = PiCFT->GetBGOHit(i);
	int segment = hit->SegmentId();
	if (FlagEvDisp) {	
	  const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	  evDisp.ShowHitBGO_Pi(segment);
	}
      }
      int NHitPiV = PiCFT->NHitPiV();
      for (int i=0; i<NHitPiV; i++) {
	Hodo2Hit *hit = PiCFT->GetPiVHit(i);
	int segment = hit->SegmentId();
	int layer   = hit->PlaneId();
	if (FlagEvDisp) {
	  const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	  evDisp.ShowHitPiV_Pi(layer, segment);
	}
      }
    }
  }


  for(int l = 0; l<NumOfLayersCFT; ++l){
    int nhit = hodoAna->GetNClustersCFT(l);
    for(int n = 0; n<nhit; ++n){
      CFTFiberCluster* cl = hodoAna->GetClusterCFT(l, n);
      if ( ! cl->showFlags()) {
	double total_pe   = cl->TotalPhotonNum();

	HF1(460+l, total_pe);
      }
    }
  }



  bool    flagDecayCal = false;
  PartVec PiDecayVec;
  PartVec DecayNeutronVec;

  PartVec SigmaVecOrg;

  // Sigma decay assumption
  if (nPi_CFT==1) {
    CFTParticle *PiCFT = CFTPionCont[0];


    PiDecayVec.pos0 = PiCFT->GetPos0();
    PiDecayVec.mom = PiCFT->GetDir();

    double cdistDecayPi;
    ThreeVector VertSigmaDecay = VertexPoint3D(SigmaVec.pos0, PiDecayVec.pos0,
					       SigmaVec.mom, PiDecayVec.mom, 
					       cdistDecayPi);
    double costDecayPi=SigmaVec.mom*PiDecayVec.mom/(SigmaVec.mom.mag()*PiDecayVec.mom.mag());
    double thetaDecayPi=acos(costDecayPi)*Rad2Deg;

    event.vertexDecayPi[0] = VertSigmaDecay.x();
    event.vertexDecayPi[1] = VertSigmaDecay.y();
    event.vertexDecayPi[2] = VertSigmaDecay.z();

    event.cdistDecayPi = cdistDecayPi;
    event.thetaDecayPi = thetaDecayPi;

    double decayPiMomCal;

    flagDecayCal = calcDecayPiMom( SigmaVec.mom.mag(), 
				   SigmaMinusMass, 
				   NeutronMass,
				   PionMass,
				   costDecayPi,
				   &decayPiMomCal);

    if (flagDecayCal) {
      event.decayPiMomCal = decayPiMomCal;      
      event.decayPiMomVecCal[0] = decayPiMomCal*PiDecayVec.mom.x()/PiDecayVec.mom.mag();
      event.decayPiMomVecCal[1] = decayPiMomCal*PiDecayVec.mom.y()/PiDecayVec.mom.mag();
      event.decayPiMomVecCal[2] = decayPiMomCal*PiDecayVec.mom.z()/PiDecayVec.mom.mag();

      PiDecayVec.mom *= (decayPiMomCal/PiDecayVec.mom.mag());

      DecayNeutronVec.mom = SigmaVec.mom-PiDecayVec.mom;
      DecayNeutronVec.pos0 = VertSigmaDecay;

      event.decayNMomCal = DecayNeutronVec.mom.mag();
      event.decayNMomVecCal[0] = DecayNeutronVec.mom.x();
      event.decayNMomVecCal[1] = DecayNeutronVec.mom.y();
      event.decayNMomVecCal[2] = DecayNeutronVec.mom.z();
    }
  }

  bool flagSigmaPScat = false;
  PartVec ProtonScatVec;
  bool flagDecayCal2 = false;
  bool flagNpScat = false;

  bool flagLambdaInv = false;
  bool flagLambdaNConv = false;
  PartVec LambdaVec;

  // Proton analysis
  if (nP_CFT==1) {
    CFTParticle *ProtonCFT = CFTProtonCont[0];
    event.EkinP = ProtonCFT->GetTotalE();
    ProtonScatVec.pos0 = ProtonCFT->GetPos0();
    ProtonScatVec.mom = ProtonCFT->GetDir();


    /* Sigma p scattering assumption */

    double cdist;

    //ThreeVector vtx0(event.primaryVertex[0], event.primaryVertex[1], event.primaryVertex[2]);
    ThreeVector vtx(event.Vertex_x[0], event.Vertex_y[0], event.Vertex_z[0]);

    ThreeVector VertScat = VertexPoint3D( SigmaVec.pos0, ProtonScatVec.pos0,
					  SigmaVec.mom, ProtonScatVec.mom, cdist);

    /* Sigma beam energy correction */
    SigmaVecOrg = SigmaVec;
    ThreeVector HypBeamLength = VertScat-vtx;
    double pSigmaBeam0 = SigmaVec.mom.mag();
    double pSigmaBeamCor, eSigmaBeamCor;
    int LH2Tag=0;

    caldE(pSigmaBeam0,  // GeV/c
	  SigmaMinusMass, // GeV
	  HypBeamLength.mag(), // mm
	  &pSigmaBeamCor, // GeV/c
	  &eSigmaBeamCor, // GeV
	  LH2Tag);

    SigmaVec.mom = SigmaVec.mom*pSigmaBeamCor/SigmaVec.mom.mag();

    event.hypBeamVec[0] = SigmaVec.mom.x();
    event.hypBeamVec[1] = SigmaVec.mom.y();
    event.hypBeamVec[2] = SigmaVec.mom.z();


    double costScat=SigmaVec.mom*ProtonScatVec.mom/(SigmaVec.mom.mag()*ProtonScatVec.mom.mag());
    double thetaScat=acos(costScat)*Rad2Deg;

    event.vertexScat[0] = VertScat.x();
    event.vertexScat[1] = VertScat.y();
    event.vertexScat[2] = VertScat.z();

    event.cdistScat = cdist;
    event.thetaScat = thetaScat;

    if (thetaScat>0. && thetaScat<90.) {
      double scatMomCal;
      double scatEkinCal;
      double thetaScatCM;

      flagSigmaPScat = calc2BodyKinema(SigmaMinusMass, SigmaVec.mom.mag(),  
				       ProtonMass, thetaScat,
				       &scatMomCal, &scatEkinCal, &thetaScatCM);

      if (flagSigmaPScat) {
	event.scatMomCal = scatMomCal;
	event.scatEkinCal = scatEkinCal*1000;
	event.thetaScatCM = thetaScatCM;
	
	double totE_p = ProtonCFT->GetTotalE();
	double p_p = sqrt(totE_p*totE_p+2.*totE_p*ProtonMass*1000.);
	
	p_p /= 1000.; // GeV/c
	double p_cor, e_cor;

	/*
	CorrElossOut(&p_cor, &e_cor, p_p, PROTON, 
		     ProtonScatVec.mom/ProtonScatVec.mom.mag(), VertScat);
	*/

	CorrElossOutWithCFRP(&p_cor, &e_cor, p_p, PROTON, 
			     ProtonScatVec.mom/ProtonScatVec.mom.mag(), VertScat,
			     ProtonScatVec.pos0);

	double Ekin_cor = (e_cor-ProtonMass)*1000.;
	totE_p = Ekin_cor;

	event.EkinCorP = totE_p;
	
	// Decay pi analysis when we measure the momntum of pi
	if (nPi_CFT == 1) {
	  
	  ThreeVector ScatProtonMom = ProtonScatVec.mom*p_cor/ProtonScatVec.mom.mag();
	  ThreeVector ScatHypMom = SigmaVec.mom-ScatProtonMom;
	  
	  ThreeVector DecayPiPos = PiDecayVec.pos0;
	  ThreeVector DecayPiVec = PiDecayVec.mom;
	  
	  double costDecayPi2=ScatHypMom*DecayPiVec/(ScatHypMom.mag()*DecayPiVec.mag());
	  double thetaDecayPi2=acos(costDecayPi2)*Rad2Deg;
	  event.thetaDecayPi2 = thetaDecayPi2;
	  
	  double decayPiMomCal2;
	  flagDecayCal2 = calcDecayPiMom(ScatHypMom.mag(), SigmaMinusMass, NeutronMass, PionMass,
					 costDecayPi2, &decayPiMomCal2);
	  if (flagDecayCal2)
	    event.decayPiMomCal2 = decayPiMomCal2;	  
	}
      }
    }
    SigmaVec = SigmaVecOrg;
    /* End of Sigma p scattering assumption */

    /* np scattering assumption */
    if (flagDecayCal) {

      double cdistNpScat;

      ThreeVector VertNpScat = VertexPoint3D( DecayNeutronVec.pos0, ProtonScatVec.pos0, 
					      DecayNeutronVec.mom, ProtonScatVec.mom, 
					      cdistNpScat);
      double costNpScat=DecayNeutronVec.mom*ProtonScatVec.mom/
	(DecayNeutronVec.mom.mag()*ProtonScatVec.mom.mag());
      double thetaNpScat=acos(costNpScat)*Rad2Deg;
	
      event.vertexNpScat[0] = VertNpScat.x();
      event.vertexNpScat[1] = VertNpScat.y();
      event.vertexNpScat[2] = VertNpScat.z();
      event.cdistNpScat = cdistNpScat;
      event.thetaNpScat = thetaNpScat;

      if (thetaNpScat>0.&&thetaNpScat<90.) {
	double scatNpMomCal;
	double scatNpEkinCal;
	double thetaNpScatCM;

	flagNpScat = calc2BodyKinema(NeutronMass, DecayNeutronVec.mom.mag(),  
				     ProtonMass, thetaNpScat,
				     &scatNpMomCal, &scatNpEkinCal, &thetaNpScatCM);
	
	if (flagNpScat) {
	  event.scatNpMomCal  = scatNpMomCal;
	  event.scatNpEkinCal = scatNpEkinCal*1000;
	  event.thetaNpScatCM = thetaNpScatCM;
	}
      }
    }
    /* End of np scattering assumption */

    /* LambdaN conversion assumption */
    if (flagDecayCal) {

	double cdistLambdaDecay;
	ThreeVector VertLambdaDecay = 
	  VertexPoint3D( PiDecayVec.pos0, ProtonScatVec.pos0, 
			 PiDecayVec.mom, ProtonScatVec.mom, cdistLambdaDecay);
	double costLambdaDecay=PiDecayVec.mom*ProtonScatVec.mom/
	  (PiDecayVec.mom.mag()*ProtonScatVec.mom.mag());
	double thetaLambdaDecay=acos(costLambdaDecay)*Rad2Deg;

	event.cdistLambdaDecay = cdistLambdaDecay;
	event.thetaLambdaDecay = thetaLambdaDecay;
	event.vertexLambdaDecay[0] = VertLambdaDecay.x();
	event.vertexLambdaDecay[1] = VertLambdaDecay.y();
	event.vertexLambdaDecay[2] = VertLambdaDecay.z();

	double Ekin1 = ProtonCFT->GetTotalE();

	Ekin1 *= 0.001; // GeV
	double M1 = ProtonMass;
	double p1 = sqrt(Ekin1*Ekin1+2.*Ekin1*M1) ;
	double p_cor, e_cor;

	/*
	CorrElossOut(&p_cor, &e_cor, p1, PROTON, ProtonScatVec.mom, VertLambdaDecay);
	*/

	CorrElossOutWithCFRP(&p_cor, &e_cor, p1, PROTON, 
			     ProtonScatVec.mom/ProtonScatVec.mom.mag(), VertLambdaDecay,
			     ProtonScatVec.pos0);

	p1 = p_cor;

	double momPiFromLambda;

	flagLambdaInv = calcPiMomFromLambda(LambdaMass, ProtonMass, p1,
					    PionMass, costLambdaDecay, &momPiFromLambda);

	if (flagLambdaInv) {
	  ThreeVector momVecPi     = PiDecayVec.mom*momPiFromLambda/PiDecayVec.mom.mag();
	  ThreeVector momVecProton = ProtonScatVec.mom*p1/ProtonScatVec.mom.mag();
	  ThreeVector momVecLambda = momVecPi+momVecProton;

	  LambdaVec.mom = momVecLambda;
	  LambdaVec.pos0 = VertLambdaDecay;

	  event.momPiFromLambda = momVecPi.mag();
	  event.momProtonFromLambda = momVecProton.mag();
	  event.momLambda = momVecLambda.mag();
	  event.momVecPiFromLambda[0] = momVecPi.x();
	  event.momVecPiFromLambda[1] = momVecPi.y();
	  event.momVecPiFromLambda[2] = momVecPi.z();
	  event.momVecProtonFromLambda[0] = momVecProton.x();
	  event.momVecProtonFromLambda[1] = momVecProton.y();
	  event.momVecProtonFromLambda[2] = momVecProton.z();
	  event.momVecLambda[0] = momVecLambda.x();
	  event.momVecLambda[1] = momVecLambda.y();
	  event.momVecLambda[2] = momVecLambda.z();


	  /*
	    std::cout << "P_pi = "<< momVecPi.mag() << ", P_proton = " << momVecProton.mag() 
	    << ", P_lambda = " << momVecLambda.mag() << std::endl;
	  */



	  double cdistLambdaNConv;
	  ThreeVector VertLambdaNConv = 
	    VertexPoint3D( SigmaVec.pos0, LambdaVec.pos0,
			   SigmaVec.mom, LambdaVec.mom, cdistLambdaNConv);

	  double costLambdaNConv=SigmaVec.mom*LambdaVec.mom/
	    (SigmaVec.mom.mag()*LambdaVec.mom.mag());
	  double thetaLambdaNConv=acos(costLambdaNConv)*Rad2Deg;


	  event.cdistLambdaNConv = cdistLambdaNConv;
	  event.thetaLambdaNConv = thetaLambdaNConv;
	  event.vertexLambdaNConv[0] = VertLambdaNConv.x();
	  event.vertexLambdaNConv[1] = VertLambdaNConv.y();
	  event.vertexLambdaNConv[2] = VertLambdaNConv.z();

	  
	  double momCalLambda =-999.;
	  double momCalLambda2 =-999.;
	  double thetaCMLambdaNConv =-999.;


	  flagLambdaNConv = calc2BodyInelastic(SigmaMinusMass, SigmaVec.mom.mag(),
					       ProtonMass, LambdaMass, NeutronMass,
					       thetaLambdaNConv,
					       &momCalLambda, &momCalLambda2,
					       &thetaCMLambdaNConv);
	  /*
	  std::cout << "SigmaMinusMass : "  << SigmaMinusMass
		    << ", Sigma mom : " << SigmaVec.mom.mag()
		    << ", ProtonMass : "  << ProtonMass 
		    << ", LambdaMass : "  << LambdaMass 
		    << ", NeutronMass : "  << NeutronMass 
		    << ", thetaLambdaNConv : "  << thetaLambdaNConv
		    << ", momCalLambda : "  << momCalLambda
		    << ", momCalLambda2 : "  << momCalLambda2
		    << ", thetaCMLambdaN "<< thetaCMLambdaNConv << std::endl;
	  */
	  if (flagLambdaNConv) {
	    event.momCalLambda = momCalLambda;
	    event.momCalLambda2 = momCalLambda2;
	    event.thetaCMLambdaNConv = thetaCMLambdaNConv;
	  }
	}
    }
    /* End of LambdaN conversion assumption */

    /* pip scattering assumption */

    if (nPi_CFT == 1) {
      //std::cout << "PiNscatFlag : " << event.PiNscatFlag << std::endl;

      double cdistPiP;

      CFTParticle *PiCFT = CFTPionCont[0];
      PartVec PiVec;      
      PiVec.pos0 = PiCFT->GetPos0();
      PiVec.mom = PiCFT->GetDir();

      ThreeVector VertPiPScat = VertexPoint3D( ProtonScatVec.pos0, PiVec.pos0,
					       ProtonScatVec.mom, PiVec.mom, cdistPiP);
      event.vertexPiPScat[0] = VertPiPScat.x();
      event.vertexPiPScat[1] = VertPiPScat.y();
      event.vertexPiPScat[2] = VertPiPScat.z();
      event.cdistPiPScat = cdistPiP;

      double costPiP=(ProtonScatVec.mom*PiVec.mom)/(ProtonScatVec.mom.mag()*PiVec.mom.mag());
      double thetaPiP = acos(costPiP)*Rad2Deg;

      ThreeVector NormalVector = gaiseki(ProtonScatVec.mom, PiVec.mom);

      double totE = ProtonCFT->GetTotalE();
      double mom_p = sqrt(totE*totE+2.*totE*ProtonMass*1000.);
      mom_p /= 1000.; // GeV/c
      double mom_cor, E_cor;

      CorrElossOutWithCFRP(&mom_cor, &E_cor, mom_p, PROTON, 
			   ProtonScatVec.mom/ProtonScatVec.mom.mag(), VertPiPScat,
			   ProtonScatVec.pos0);
      mom_p = mom_cor;

      std::vector <double> ThetaCont, PiBeamCont1, PiBeamCont2, PiBeamDiffCont1, PiBeamDiffCont2;

      double thetaNormVec, phiNormVec;
      calcThetaPhi(NormalVector, &thetaNormVec, &phiNormVec);

      for (double theta=0.; theta<=thetaPiP; theta+=0.1 ) {
	ThreeVector PiBeam = ProtonScatVec.mom;
	PiBeam.rotate(theta*Deg2Rad, NormalVector/NormalVector.mag());

	/*
	ThreeVector Norm(0, 0, 1);
	ThreeVector Vec(1, 0, 0);
	ThreeVector TmpBeam = Vec.rotate(theta*Deg2Rad, Norm);
	std::cout << "theta : " << theta << "(x, y, z) = (" << TmpBeam.x()<< ", " << TmpBeam.y() << ", " << TmpBeam.z()
		  << ")" << std::endl;
	std::cout << "theta : " << theta << "(x, y, z) = (" << Vec.x()<< ", " << Vec.y() << ", " << Vec.z()
		  << ")" << std::endl;
	*/

	double beamMomCal1, beamMomCal2;
	bool flagPiBeam = calcBeamMomFrom2BodyKinema(PionMass, ProtonMass, theta, mom_p, &beamMomCal1, &beamMomCal2) ;
	
	//std::cout << theta << " ";

	if (beamMomCal1>0) {
	  double E1 = sqrt(beamMomCal1*beamMomCal1+PionMass*PionMass);
	  double E3 = E_cor;
	  double p3 = sqrt(E3*E3 - ProtonMass*ProtonMass);
	  double E4 = E1 + ProtonMass - E3;
	  double p4 = sqrt(E4*E4 - PionMass*PionMass);

	  double momconv_x = beamMomCal1 - (p3*cos(theta*Deg2Rad)+p4*cos((thetaPiP-theta)*Deg2Rad));
	  double momconv_y = (p3*sin(theta*Deg2Rad)-p4*sin((thetaPiP-theta)*Deg2Rad));

	  //std::cout << momconv_x << " " << momconv_y << " " << fabs(momconv_x)+fabs(momconv_y) << " ";
	  //std::cout << "1 momconv_x : " << momconv_x << ", momconv_y : " << momconv_y << std::endl;
	} else {
	  //std::cout << -999 << " " << -999 << " "<< -999 << " ";
	}

	if (beamMomCal2>0) {
	  double E1 = sqrt(beamMomCal2*beamMomCal2+PionMass*PionMass);
	  double E3 = E_cor;
	  double p3 = sqrt(E3*E3 - ProtonMass*ProtonMass);
	  double EscatPi = E1 + ProtonMass - E3;
	  double E4 = sqrt(EscatPi*EscatPi - PionMass*PionMass);
	  double p4 = sqrt(E4*E4 - PionMass*PionMass);

	  double momconv_x = beamMomCal2 - (p3*cos(theta*Deg2Rad)+p4*cos((thetaPiP-theta)*Deg2Rad));
	  double momconv_y = (p3*sin(theta*Deg2Rad)-p4*sin((thetaPiP-theta)*Deg2Rad));

	  //std::cout << momconv_x << " " << momconv_y << " " << fabs(momconv_x)+fabs(momconv_y) << " ";

	  //std::cout << "1 momconv_x : " << momconv_x << ", momconv_y : " << momconv_y << std::endl;
	} else {
	  //std::cout << -999 << " " << -999 << " "<< -999 << " ";
	}


	if (flagPiBeam) {
	  double costPiDecay = (PiBeam*SigmaVec.mom)/(PiBeam.mag()*SigmaVec.mom.mag());
	  double beamMomCalFromSigmaDecay;
	  bool flagPiBeam2 = calcDecayPiMom(SigmaVec.mom.mag(), SigmaMinusMass, NeutronMass, PionMass, costPiDecay, &beamMomCalFromSigmaDecay);

	  if (flagPiBeam2) {
	    //std::cout << "Theta : " << theta << ", beamMomCal1 : " << beamMomCal1 << ", beamMomCal2 " <<  beamMomCal2
	    //<< ", Diff1 : " << fabs(beamMomCal1-beamMomCalFromSigmaDecay) 
	    //<< ", Diff2 : " << fabs(beamMomCal2-beamMomCalFromSigmaDecay) << std::endl;
	    /*
	    std::cout << costPiDecay << " " << beamMomCal1 << "  " <<  beamMomCal2 << " " << beamMomCalFromSigmaDecay 
		      << " " << fabs(beamMomCal1-beamMomCalFromSigmaDecay) 
		      << " " << fabs(beamMomCal2-beamMomCalFromSigmaDecay) << std::endl;
	    */
	    ThetaCont.push_back(theta);
	    PiBeamCont1.push_back(beamMomCal1);
	    PiBeamCont2.push_back(beamMomCal2);
	    PiBeamDiffCont1.push_back(fabs(beamMomCal1-beamMomCalFromSigmaDecay));
	    PiBeamDiffCont2.push_back(fabs(beamMomCal2-beamMomCalFromSigmaDecay));
	  }
	}
      }

      int nSize=ThetaCont.size();
      int index1=-1;
      double min = 10000000.;
      for (int i=0; i<nSize; i++) {
	if (PiBeamDiffCont1[i]<=min) {
	  min = PiBeamDiffCont1[i];
	  index1 = i;
	}
      }

      double scatPiMomCal1,scatPiMomCal2;
      double scatPiMomCalFromE1,scatPiMomCalFromE2;

      if (index1>=0 && index1 <nSize) {
	double piBeamMom = PiBeamCont1[index1];
	double Theta = ThetaCont[index1];
	double ThetaScatPi = thetaPiP - Theta;

	double mom1, mom2, thetaPiCM;
	bool flagPiBeam3 = calc2BodyInelastic(PionMass, piBeamMom, ProtonMass, PionMass, ProtonMass,
					      ThetaScatPi, &mom1, &mom2, &thetaPiCM);

	double E1 = sqrt(piBeamMom*piBeamMom+PionMass*PionMass);
	double E3 = E_cor;
	double EscatPi = E1 + ProtonMass - E3;
	double scatPiMomCalFromE = sqrt(EscatPi*EscatPi - PionMass*PionMass);

	scatPiMomCal1 = mom1;
	scatPiMomCalFromE1 = scatPiMomCalFromE;

	//std::cout << "1 scatPiMomCal1 : " << scatPiMomCal1 <<  ", scatPiMomCalFromE : " << scatPiMomCalFromE 
	//<< "Delta : " << scatPiMomCal1 - scatPiMomCalFromE  << std::endl;
      }

      int index2=-1;
      min = 10000000.;
      for (int i=0; i<nSize; i++) {
	if (PiBeamDiffCont2[i]<=min) {
	  min = PiBeamDiffCont2[i];
	  index2 = i;
	}
      }
      
      if (index2>=0 && index2 <nSize) {
	double piBeamMom = PiBeamCont2[index2];
	double Theta = ThetaCont[index2];
	double ThetaScatPi = thetaPiP - Theta;

	double mom1, mom2, thetaPiCM;
	bool flagPiBeam3 = calc2BodyInelastic(PionMass, piBeamMom, ProtonMass, PionMass, ProtonMass,
					      ThetaScatPi, &mom1, &mom2, &thetaPiCM);

	double E1 = sqrt(piBeamMom*piBeamMom+PionMass*PionMass);
	double E3 = E_cor;
	double EscatPi = E1 + ProtonMass - E3;
	double scatPiMomCalFromE = sqrt(EscatPi*EscatPi - PionMass*PionMass);

	scatPiMomCal2 = mom1;
	scatPiMomCalFromE2 = scatPiMomCalFromE;
	
	//std::cout << "2 scatPiMomCal1 : " << scatPiMomCal1 <<  ", scatPiMomCalFromE : " << scatPiMomCalFromE 
	//<< ", Delta : " << scatPiMomCal1 - scatPiMomCalFromE  << std::endl;

      }

      if (fabs(scatPiMomCal1-scatPiMomCalFromE1) <= fabs(scatPiMomCal2-scatPiMomCalFromE2)) {
	event.piBeamMom = PiBeamCont1[index1];
	event.thetaPiPScat = ThetaCont[index1];
	event.scatPiMomCal = scatPiMomCal1;
	event.scatPiMomCalFromE = scatPiMomCalFromE1;
      } else {
	event.piBeamMom = PiBeamCont2[index2];
	event.thetaPiPScat = ThetaCont[index2];
	event.scatPiMomCal = scatPiMomCal2;
	event.scatPiMomCalFromE = scatPiMomCalFromE2;
      }
    }

  }

  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    //if (1) {
    //if (ntCFT_2nd>0) {
    //if (nPi_CFT>=2 && nP_CFT==1) {
    //if (event.scatFlag==-1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1&&nPi_CFT==1) {
    //if (event.scatFlag==-1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1) {
    //if (nP_CFT==1&&nPi_CFT==1) {
    if (event.scatFlag==1) {
      std::cout << "HypBeamMom : " << event.momHypBeam << ", cost : " << cos(event.thetaScatHypCM*3.14/180) << std::endl;
      std::cout << "Delta U : " << event.u0SdcIn - event.u0CFT_Pi[0] << ", Delta V : " << event.v0SdcIn - event.v0CFT_Pi[0] << std::endl;
      evDisp.UpdateCanvas();
      evDisp.get_command();
    }
    evDisp.EndOfEvent();
  }


  tree->Fill();

  for_each(SksPartCont.begin(), SksPartCont.end(), DeleteObject());
  //for_each(KaonPartCont.begin(), KaonPartCont.end(), DeleteObject());
  for_each(CFTPartCont.begin(), CFTPartCont.end(), DeleteObject());


  delete hodoAna;
  delete DCAna;
  delete rawData;

  return true;
}


bool calcDecayPiMom(double p1, double m1, double m2, double m3, double cost, double *momCal)
/*
  p1 : Sigma mom
  m1 : Sigma mass
  m2 : Neutron mass
  m3 : Pion mass
*/
{
  double E1 = sqrt(p1*p1+m1*m1);
  double A  = (m1*m1+m3*m3-m2*m2)/2.;
  double hanbetu = (A*p1*cost)*(A*p1*cost)-(E1*E1-p1*p1*cost*cost)*(E1*E1*m3*m3-A*A);
  if (hanbetu >= 0) {
    double ans1 = (A*p1*cost+sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    double ans2 = (A*p1*cost-sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    if (ans1>=0 && ans2<0) {
      *momCal = ans1;	  
      return true;
    } else if (ans2>=0 && ans1<0) {
      *momCal = ans2;	  
      return true;
    } else if (ans1>=0 && ans2>=0) {
      std::cout << "decayPiMomCal two answers: " << ans1 << ", " << ans2
		<< std::endl;
      return false;
    } else if (ans1<0 && ans2<0) {
      std::cout << "decayPiMomCal two negative answers: " 
		<< ans1 << ", " << ans2
		<< std::endl;
      return false;
    }
  }

  return false;
}

bool calcBeamMomFrom2BodyKinema(double M1, double M2, double theta, double pScat,
				double *beamMomCal1, double *beamMomCal2)
/*
  M1 : pi mass
  M2 : proton mass
  theta : scat proton angle
  pScat : scat proton momentum
*/
{
  double p3 = pScat;
  double E3 = sqrt(M2*M2+p3*p3);

  double A = (M2-E3)*(M2-E3)-p3*p3*cos(theta*Deg2Rad)*cos(theta*Deg2Rad);
  double B = (E3*M2-M2*M2)*p3*cos(theta*Deg2Rad);
  double C = (M2-E3)*(M2-E3)*M1*M1-(E3*M2-M2*M2)*(E3*M2-M2*M2);

  double hanbetsu = B*B-A*C;
  if (hanbetsu>=0) {
    double p1 = (-B+sqrt(hanbetsu))/A;
    double p2 = (-B-sqrt(hanbetsu))/A;

    //std::cout << "p1 : " << p1 << ", p2 : " << p2 << std::endl;

    *beamMomCal1 = p1;
    *beamMomCal2 = p2;
    return true;
  }

  return false;
}


bool calc2BodyKinema(double M1, double p1, double M2, double phi,
		     double *scatMomCal, double *scatEkinCal, double *scatThetaCM)
/*
  M1 : sigma mass
  p1 : sigma mom
  M2 : proton mass
  phi : scat proton angle

*/
{
  double M = M2;
  double E1 = sqrt(M1*M1+p1*p1);
  double p4 = (2.*M*(E1+M)*p1*cos(phi*Deg2Rad))/
	((E1+M)*(E1+M)-p1*p1*cos(phi*Deg2Rad)*cos(phi*Deg2Rad));
  double E4 = sqrt(M*M+p4*p4);
  double Ekin4 = E4 - M;
	    
  double beta = p1/(E1+M);
  double gamma = 1./sqrt(1.-beta*beta);

  double M3 = M1;
  double Ecm = sqrt((E1+M)*(E1+M)-p1*p1);
  double Pcm = 1./(2.*Ecm)*
    sqrt((Ecm-M3+M)*(Ecm-M3-M)*(Ecm+M3-M)*(Ecm+M3+M));
  double sintCM = p4*sin(phi*Deg2Rad)/Pcm;
  double thetaCM = asin(sintCM)*Rad2Deg;

  double P4cm_para = -beta*gamma*E4 + 
    gamma*p4*cos(phi*Deg2Rad);
  if (P4cm_para>0)
    thetaCM = 180. - thetaCM;

  *scatMomCal = p4;
  *scatEkinCal = Ekin4;
  *scatThetaCM = thetaCM;

  return true;
}


bool calcPiMomFromLambda(double M0, double M1, double p1, double M2, double cost, double *pCal)
/*
  M0 : Lambda mass
  M1 : Proton Mass
  p1 : Proton mom
  M2 : Pion Mass
  cost : cos(theta) of opening angle between pi and proton
*/
{
  // Calculate momentum of pi assume Lambda decay
  double E1 = sqrt(p1*p1+M1*M1);
	      
  double A = (M0*M0-(M1*M1+M2*M2))/2.;
  //std::cout << "p1 = " << p1 << ", Ekin1 = " << Ekin1 << std::endl;

  double hanbetu = (A*p1*cost)*(A*p1*cost)-(E1*E1-p1*p1*cost*cost)*(E1*E1*M2*M2-A*A);
  if (hanbetu>=0) {
    double ans1 = (A*p1*cost+sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    double ans2 = (A*p1*cost-sqrt(hanbetu))/(E1*E1-p1*p1*cost*cost);
    /*		
      std::cout << "--pi momentum" << std::endl;
      if (ans1>=0 && ans2>= 0) {
      std::cout << "1 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else if (ans1>=0) {
      std::cout << "2 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else if (ans2>=0) {
      std::cout << "3 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else {
      std::cout << "4 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      }
    */
    
    if (ans1>=0) {
      double p2 = ans1;
      *pCal = p2;

      return true;
    }
  }

  return false;
}

bool calc2BodyInelastic(double m1, double p1, double m2, double m3, double m4, 
			double theta, double *pCal1, double *pCal2, double *ThetaCM)
/*
  m1 = SigmaMinusMass;
  m2 = ProtonMass;
  m3 = LambdaMass;
  m4 = NeutronMass;
  theta = thetaLambdaNConv;
*/
{

  double E1 = sqrt(p1*p1+m1*m1);
  double A = m1*m1+m2*m2+m3*m3-m4*m4+2.*E1*m2;

  double cost = cos(theta*Deg2Rad);

  double hanbetu = 4.*m3*m3*(p1*p1*cost*cost-(E1+m2)*(E1+m2))+A*A;
  if (hanbetu>=0) {
    double ans1 = (A*p1*cost+(E1+m2)*sqrt(hanbetu))/(2.*((E1+m2)*(E1+m2)-p1*p1*cost*cost));
    double ans2 = (A*p1*cost-(E1+m2)*sqrt(hanbetu))/(2.*((E1+m2)*(E1+m2)-p1*p1*cost*cost));
    *pCal1 = ans1;
    *pCal2 = ans2;

    /*
      std::cout << "---Calculated Lambda mom" << std::endl;
      if (ans1>=0 && ans2>= 0) {
      std::cout << "1 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else if (ans1>=0) {
      std::cout << "2 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else if (ans2>=0) {
      std::cout << "3 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      } else {
      std::cout << "4 : ans1 = " << ans1 << ", ans2 = " << ans2 << std::endl;
      }
    */		    
    double beta = p1/(E1+m2);
    double gamma = 1./sqrt(1.-beta*beta);
		    
    double Ecm = sqrt((E1+m2)*(E1+m2)-p1*p1);
    double Pcm = 1./(2.*Ecm)*
      sqrt((Ecm-m3+m4)*(Ecm-m3-m4)*(Ecm+m3-m4)*(Ecm+m3+m4));
    //double p3 = event.momCalLambda;
    double p3 = ans1;
    double E3 = sqrt(p3*p3+m3*m3);
    double sintCM = p3*sin(theta*Deg2Rad)/Pcm;
    double thetaCM = asin(sintCM)*Rad2Deg;
		    
    double P3cm_para = -beta*gamma*E3 + 
      gamma*p3*cos(theta*Deg2Rad);
    if (P3cm_para<0)
      thetaCM = 180. - thetaCM;
    
    *ThetaCM = thetaCM;
    
    return true;
  }

  return false;
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

ThreeVector gaiseki(ThreeVector vec1, ThreeVector vec2)
{
  double vec_x = vec1.y()*vec2.z()-vec2.y()*vec1.z();
  double vec_y = vec1.z()*vec2.x()-vec2.z()*vec1.x();
  double vec_z = vec1.x()*vec2.y()-vec2.x()*vec1.y();

  return ThreeVector(vec_x, vec_y, vec_z);
}

void calcThetaPhi(ThreeVector vec, double *theta, double *phi)
{
  *theta = acos(vec.z()/vec.mag())*Rad2Deg;

  if (*theta<0.)
    *theta *= -1.;


  if (vec.x()>=0. && vec.y()>=0. ) {
    //std::cout << "1: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = atan(vec.y()/vec.x())*Rad2Deg;
  } else if (vec.x()<0. && vec.y()>=0. ) {
    //std::cout << "2: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = 180.+atan(vec.y()/vec.x())*Rad2Deg;
  } else if (vec.x()<0. && vec.y()<0. ) {
    //std::cout << "3: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = 360.-acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = 180.+atan(vec.y()/vec.x())*Rad2Deg;
  } else if (vec.x()>=0. && vec.y()<0. ) {
    //std::cout << "4: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = 360.-acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = 360.+atan(vec.y()/vec.x())*Rad2Deg;
  } else {
    fprintf(stderr, "PrimaryGeneratorAction::calcThetaPhi strange vector (%f, %f, %f)", vec.x(), vec.y(), vec.z());
    *phi = -1.;
  }
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

  event.ntSks=0;
  event.ntSksPart=0;
  event.nK=0;
  event.nSigma=0;
  for (int i=0; i<MaxHits; i++) {
    event.p[i]=-999.0;
    event.pK[i]=-999.0;
    event.pKCal[i]=-999.0;
    event.theta[i]=-999.0;
    event.MissMass[i]=-999.0;
    event.pSigmaCal[i]=-999.0;
    event.Vertex_x[i]=-999.0;
    event.Vertex_y[i]=-999.0;
    event.Vertex_z[i]=-999.0;

    event.pSksPart[i]=-999.0;
    event.m2[i]=-999.0;
  }
  event.u0SdcIn = -999.;
  event.v0SdcIn = -999.;

  event.resSftU = -999.;
  event.resSftV = -999.;


  event.nPi_CFT=0;
  event.nK_CFT=0;
  event.nP_CFT=0;
  event.EkinP = -999.0;
  event.EkinCorP = -999.0;

  event.vertexDecayPi[0] = -999.9;
  event.vertexDecayPi[1] = -999.9;
  event.vertexDecayPi[2] = -999.9;
  event.cdistDecayPi = -999.9;
  event.thetaDecayPi = -999.9;
  event.decayPiMomCal = -999.9;
  event.decayPiMomVecCal[0] = -999.9;
  event.decayPiMomVecCal[1] = -999.9;
  event.decayPiMomVecCal[2] = -999.9;
  event.decayNMomCal = -999.9;
  event.decayNMomVecCal[0] = -999.9;
  event.decayNMomVecCal[1] = -999.9;
  event.decayNMomVecCal[2] = -999.9;

  event.vertexScat[0] = -999.9;
  event.vertexScat[1] = -999.9;
  event.vertexScat[2] = -999.9;

  event.hypBeamVec[0] = -999.9;
  event.hypBeamVec[1] = -999.9;
  event.hypBeamVec[2] = -999.9;
  
  event.cdistScat = -999.9;
  event.thetaScat = -999.9;
  event.scatMomCal = -999.9;
  event.scatEkinCal = -999.9;
  event.thetaScatCM = -999.9;

  event.thetaDecayPi2 = -999.9;
  event.decayPiMomCal2 = -999.9;

  event.vertexNpScat[0] = -999.9;
  event.vertexNpScat[1] = -999.9;
  event.vertexNpScat[2] = -999.9;

  event.cdistNpScat = -999.9;
  event.thetaNpScat = -999.9;
  event.scatNpMomCal = -999.9;
  event.scatNpEkinCal = -999.9;
  event.thetaNpScatCM = -999.9;

  event.cdistLambdaDecay = -999.9;
  event.thetaLambdaDecay = -999.9;
  event.vertexLambdaDecay[0] = -999.9;
  event.vertexLambdaDecay[1] = -999.9;
  event.vertexLambdaDecay[2] = -999.9;

  event.momPiFromLambda = -999.9;
  event.momProtonFromLambda = -999.9;
  event.momLambda = -999.9;
  event.momVecPiFromLambda[0] = -999.9;
  event.momVecPiFromLambda[1] = -999.9;
  event.momVecPiFromLambda[2] = -999.9;
  event.momVecProtonFromLambda[0] = -999.9;
  event.momVecProtonFromLambda[1] = -999.9;
  event.momVecProtonFromLambda[2] = -999.9;
  event.momVecLambda[0] = -999.9;
  event.momVecLambda[1] = -999.9;
  event.momVecLambda[2] = -999.9;

  event.cdistLambdaNConv = -999.9;
  event.thetaLambdaNConv = -999.9;
  event.vertexLambdaNConv[0] = -999.9;
  event.vertexLambdaNConv[1] = -999.9;
  event.vertexLambdaNConv[2] = -999.9;

  event.momCalLambda = -999.9;
  event.momCalLambda2 = -999.9;
  event.thetaCMLambdaNConv = -999.9;

  event.vertexPiPScat[0] = -999.9;
  event.vertexPiPScat[1] = -999.9;
  event.vertexPiPScat[2] = -999.9;
  event.cdistPiPScat = -999.9;

  event.piBeamMom = -999.9;
  event.thetaPiPScat = -999.9;
  event.scatPiMomCal = -999.9;
  event.scatPiMomCalFromE = -999.9;


  event.ntCFT = 0;
  for (int i=0; i<MaxHits2; i++) {
    event.xDirCFT[i]=-999.;
    event.yDirCFT[i]=-999.;
    event.zDirCFT[i]=-999.;
    event.xPos0CFT[i]=-999.;
    event.yPos0CFT[i]=-999.;
    event.zPos0CFT[i]=-999.;
    event.zBGOCFT[i]=-999.;

    event.BGO_Edep[i]=-999.;
    event.TotalEdep[i]=-999.;
    event.CFT_TotalEdep[i]=-999.;
    event.CFT_NormTotalEdep[i]=-999.;
    event.PiV_Edep[i]=-999.;

    event.xDirCFT_P[i]=-999.;
    event.yDirCFT_P[i]=-999.;
    event.zDirCFT_P[i]=-999.;
    event.xPos0CFT_P[i]=-999.;
    event.yPos0CFT_P[i]=-999.;
    event.zPos0CFT_P[i]=-999.;
    event.zBGOCFT_P[i]=-999.;

    event.CFTVtx_x_P[i]=-999.;
    event.CFTVtx_y_P[i]=-999.;
    event.CFTVtx_z_P[i]=-999.;

    event.BGO_Edep_P[i]=-999.;
    event.TotalEdep_P[i]=-999.;
    event.CFT_TotalEdep_P[i]=-999.;
    event.CFT_NormTotalEdep_P[i]=-999.;
    event.PiV_Edep_P[i]=-999.;

    event.xDirCFT_Pi[i]=-999.;
    event.yDirCFT_Pi[i]=-999.;
    event.zDirCFT_Pi[i]=-999.;
    event.xPos0CFT_Pi[i]=-999.;
    event.yPos0CFT_Pi[i]=-999.;
    event.zPos0CFT_Pi[i]=-999.;
    event.zBGOCFT_Pi[i]=-999.;

    event.CFTVtx_x_Pi[i]=-999.;
    event.CFTVtx_y_Pi[i]=-999.;
    event.CFTVtx_z_Pi[i]=-999.;

    event.u0CFT_Pi[i]=-999.;
    event.v0CFT_Pi[i]=-999.;

    event.BGO_Edep_Pi[i]=-999.;
    event.TotalEdep_Pi[i]=-999.;
    event.CFT_TotalEdep_Pi[i]=-999.;
    event.CFT_NormTotalEdep_Pi[i]=-999.;
    event.PiV_Edep_Pi[i]=-999.;
  }


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


  HB1(30, "N track SdcOut", 10, 0, 10);
  HB1(31, "N Hit SdcOut tracking", 12, 0, 12);
  HB1(32, "Chisquare SdcOut tracking", 100, 0, 50);
  HB1(33, "Hit Layer  SdcOut tracking", 12, 0, 12);
  HB1(34, "X0 (TOF)  SdcOut tracking", 200, -600, 600);
  HB1(35, "Y0 (TOF)  SdcOut tracking", 200, -600, 600);
  HB1(36, "U0  SdcOut tracking", 200, -1., 1.);
  HB1(37, "V0  SdcOut tracking", 200, -1., 1.);
  HB2(38, "X0%U0  SdcOut tracking", 200, -600, 600, 200, -1, 1);
  HB2(39, "Y0%V0  SdcOut tracking", 200, -600, 600, 200, -1, 1);
  HB2(40, "X0%Y0  SdcOut tracking", 200, -600, 600, 200, -600, 600);


  HB1(50, "N track Sks", 10, 0, 10);
  HB1(51, "N Hit Sks tracking", 12, 0, 12);
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

  HB1(401, "N track in CFT tracking", 10, 0, 10);
  HB1(402, "N Hit in CFT tracking (XY)", 10, 0, 10);
  HB1(403, "N Hit in CFT tracking (Z)", 10, 0, 10);
  HB1(404, "Chisqr in CFT tracking (XY)", 100, 0, 50);
  HB1(405, "Chisqr in CFT tracking (Z)", 100, 0, 50);
  HB1(406, "X0 in CFT tracking ", 100, -100, 100);
  HB1(407, "Y0 in CFT tracking ", 100, -100, 100);
  HB1(408, "U0 in CFT tracking ", 100, -2, 2);
  HB1(409, "V0 in CFT tracking ", 100, -2, 2);
  HB2(410, "Y0%X0 in CFT tracking ", 100, -100, 100, 100, -100, 100);
  HB2(411, "V0%U0 in CFT tracking ", 100, -2, 2, 100, -2, 2);
  HB1(412, "Total DE", 100, 0, 20);
  HB1(413, "Total Maximum DE", 100, 0, 20);
  HB1(414, "Path Length in CFT", 200, 0, 20);
  HB1(415, "Normalized Total DE", 200, 0, 5);
  HB1(416, "Normalized Total Maximum DE", 200, 0, 5);
  HB2(417, "Total DE%Pathlength",200, 0, 20,  100, 0, 20);
  HB2(418, "Total Maximum DE%Pathlength", 200, 0, 20, 100, 0, 20);
  HB2(419, "Normalized Total DE%Pathlength",200, 0, 20,  100, 0, 5);
  HB2(420, "Normalized Total Maximum DE%Pathlength", 200, 0, 20, 100, 0, 5);
  HB1(421, "Delta Pos at BGO", 400, -100, 100);
  HB2(422, "Total DE%dE BGO",150, 0, 150,  100, 0, 50);
  HB2(423, "Total Maximum DE%dE BGO", 150, 0, 150, 100, 0, 50);
  HB2(424, "Normalized Total DE%dE BGO",150, 0, 150,  100, 0, 10);
  HB2(425, "Normalized Total Maximum DE%dE BGO", 150, 0, 150, 100, 0, 10);
  HB2(426, "Total DE%dE BGO+Fiber(total)",150, 0, 150,  100, 0, 50);
  HB2(427, "Total Maximum DE%dE BGO+Fiber(max)", 150, 0, 150, 100, 0, 50);
  HB2(428, "Normalized Total DE%dE BGO+Fiber(total)",150, 0, 150,  100, 0, 5);
  HB2(429, "Normalized Total Maximum DE%dE BGO+Fiber(max)", 150, 0, 150, 100, 0, 5);

  HB1(430, "#Delta U (CFT-Sp)", 100, -0.5, 0.5);
  HB1(431, "#Delta V (CFT-Sp)", 100, -0.5, 0.5);
  HB2(432, "#Delta V%#Delta U (CFT-Sp)", 100, -0.5, 0.5, 100, -0.5, 0.5);
  HB1(433, "#Delta dy/dx (CFT-Sp)", 100, -0.5, 0.5);
  HB1(434, "cdist (CFT-Sp)", 200, 0, 100);
  HB1(435, "Vertex Z (CFT-Sp)", 200, -300, 500);
  HB2(436, "Vertex XY (CFT-Sp)", 200, -200, 200, 200, -200, 200);

  HB1(450, "Delta Pos at PiV", 400, -100, 100);

  for (int l=0; l<NumOfLayersCFT; l++) {
    sprintf(buf, "#DeltaE (p) layer = %d", l);
    HB1(440+l, buf, 200, 0, 10);
  }

  for (int l=0; l<NumOfLayersCFT; l++) {
    sprintf(buf, "#DeltaE (pi) layer = %d", l);
    HB1(450+l, buf, 200, 0, 10);
  }

  for (int l=0; l<NumOfLayersCFT; l++) {
    sprintf(buf, "#DeltaE (Not Tracked) layer = %d", l);
    HB1(460+l, buf, 200, 0, 10);
  }

  HB1(501, "N track in CFT tracking (2nd)", 10, 0, 10);
  HB1(502, "N Hit in CFT tracking (XY) (2nd)", 10, 0, 10);
  HB1(503, "N Hit in CFT tracking (Z) (2nd)", 10, 0, 10);
  HB1(504, "Chisqr in CFT tracking (XY) (2nd)", 100, 0, 50);
  HB1(505, "Chisqr in CFT tracking (Z) (2nd)", 100, 0, 50);
  HB1(506, "X0 in CFT tracking  (2nd)", 100, -100, 100);
  HB1(507, "Y0 in CFT tracking  (2nd)", 100, -100, 100);
  HB1(508, "U0 in CFT tracking  (2nd)", 100, -2, 2);
  HB1(509, "V0 in CFT tracking  (2nd)", 100, -2, 2);
  HB2(510, "Y0%X0 in CFT tracking  (2nd)", 100, -100, 100, 100, -100, 100);
  HB2(511, "V0%U0 in CFT tracking  (2nd)", 100, -2, 2, 100, -2, 2);
  HB1(512, "Total DE (2nd)", 100, 0, 40);
  HB1(513, "Total Maximum DE (2nd)", 100, 0, 40);
  HB1(514, "Path Length in CFT (2nd)", 200, 0, 20);
  HB1(515, "Normalized Total DE (2nd)", 200, 0, 10);
  HB1(516, "Normalized Total Maximum DE (2nd)", 200, 0, 10);
  HB2(517, "Total DE%Pathlength (2nd)",200, 0, 20,  100, 0, 20);
  HB2(518, "Total Maximum DE%Pathlength (2nd)", 200, 0, 20, 100, 0, 20);
  HB2(519, "Normalized Total DE%Pathlength (2nd)",200, 0, 20,  100, 0, 5);
  HB2(520, "Normalized Total Maximum DE%Pathlength (2nd)", 200, 0, 20, 100, 0, 5);
  HB1(521, "Delta Pos at BGO (2nd)", 400, -100, 100);
  HB2(522, "Total DE%dE BGO (2nd)",150, 0, 150,  100, 0, 50);
  HB2(523, "Total Maximum DE%dE BGO (2nd)", 150, 0, 150, 100, 0, 50);
  HB2(524, "Normalized Total DE%dE BGO (2nd)",150, 0, 150,  100, 0, 10);
  HB2(525, "Normalized Total Maximum DE%dE BGO (2nd)", 150, 0, 150, 100, 0, 10);
  HB2(526, "Total DE%dE BGO+Fiber(total) (2nd)",150, 0, 150,  100, 0, 50);
  HB2(527, "Total Maximum DE%dE BGO+Fiber(max) (2nd)", 150, 0, 150, 100, 0, 50);
  HB2(528, "Normalized Total DE%dE BGO+Fiber(total) (2nd)",150, 0, 150,  100, 0, 5);
  HB2(529, "Normalized Total Maximum DE%dE BGO+Fiber(max) (2nd)", 150, 0, 150, 100, 0, 5);

  HB1(530, "#Delta U (CFT-Sp) (2nd)", 100, -0.5, 0.5);
  HB1(531, "#Delta V (CFT-Sp) (2nd)", 100, -0.5, 0.5);
  HB2(532, "#Delta V%#Delta U (CFT-Sp) (2nd)", 100, -0.5, 0.5, 100, -0.5, 0.5);
  HB1(533, "#Delta dy/dx (CFT-Sp) (2nd)", 100, -0.5, 0.5);
  HB1(534, "cdist (CFT-Sp) (2nd)", 200, 0, 100);
  HB1(535, "Vertex Z (CFT-Sp) (2nd)", 200, -300, 500);
  HB2(536, "Vertex XY (CFT-Sp) (2nd)", 200, -200, 200, 200, -200, 200);


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

  tree->Branch("ntSks",   &event.ntSks,  "ntSks/I");
  tree->Branch("p",   &event.p,  "p[ntSks]/D");

  tree->Branch("ntSksPart",   &event.ntSksPart,  "ntSksPart/I");
  tree->Branch("pSksPart",   event.pSksPart,  "pSksPart[ntSksPart]/D");
  tree->Branch("m2",   event.m2,  "m2[ntSksPart]/D");

  tree->Branch("nK",   &event.nK,  "nK/I");
  tree->Branch("pK",   event.pK,  "pK[nK]/D");
  tree->Branch("MissMass",   event.MissMass,  "MissMass[nK]/D");
  tree->Branch("u0SdcIn",   &event.u0SdcIn,  "u0SdcIn/D");
  tree->Branch("v0SdcIn",   &event.v0SdcIn,  "v0SdcIn/D");

  tree->Branch("resSftU",   &event.resSftU,  "resSftU/D");
  tree->Branch("resSftV",   &event.resSftV,  "resSftV/D");

  tree->Branch("nSigma",   &event.nSigma,  "nSigma/I");
  tree->Branch("pKCal",   event.pKCal,  "pKCal[nSigma]/D");
  tree->Branch("theta",   event.theta,  "theta[nSigma]/D");
  tree->Branch("pSigmaCal",   event.pSigmaCal,  "pSigmaCal[nSigma]/D");
  tree->Branch("Vertex_x",   event.Vertex_x,  "Vertex_x[nSigma]/D");
  tree->Branch("Vertex_y",   event.Vertex_y,  "Vertex_y[nSigma]/D");
  tree->Branch("Vertex_z",   event.Vertex_z,  "Vertex_z[nSigma]/D");

  tree->Branch("nPi_CFT",   &event.nPi_CFT,  "nPi_CFT/I");
  tree->Branch("nK_CFT",   &event.nK_CFT,  "nK_CFT/I");
  tree->Branch("nP_CFT",   &event.nP_CFT,  "nP_CFT/I");
  tree->Branch("EkinP",   &event.EkinP,  "EkinP/D");
  tree->Branch("EkinCorP",   &event.EkinCorP,  "EkinCorP/D");

  tree->Branch("vertexDecayPi",   event.vertexDecayPi,  "vertexDecayPi[3]/D");
  tree->Branch("cdistDecayPi",   &event.cdistDecayPi,  "cdistDecayPi/D");
  tree->Branch("thetaDecayPi",   &event.thetaDecayPi,  "thetaDecayPi/D");
  tree->Branch("decayPiMomCal",   &event.decayPiMomCal,  "decayPiMomCal/D");
  tree->Branch("decayPiMomVecCal",   event.decayPiMomVecCal,  "decayPiMomVecCal/D");
  tree->Branch("decayNMomCal",   &event.decayNMomCal,  "decayNMomCal/D");
  tree->Branch("decayNMomVecCal",   event.decayNMomVecCal,  "decayNMomVecCal[3]/D");

  tree->Branch("vertexScat",   event.vertexScat,  "vertexScat[3]/D");
  tree->Branch("hypBeamVec",   event.hypBeamVec,  "hypBeamVec[3]/D");
  tree->Branch("cdistScat",   &event.cdistScat,  "cdistScat/D");
  tree->Branch("thetaScat",   &event.thetaScat,  "thetaScat/D");
  tree->Branch("scatMomCal",   &event.scatMomCal,  "scatMomCal/D");
  tree->Branch("scatEkinCal",   &event.scatEkinCal,  "scatEkinCal/D");
  tree->Branch("thetaScatCM",   &event.thetaScatCM,  "thetaScatCM/D");

  tree->Branch("thetaDecayPi2",   &event.thetaDecayPi2,  "thetaDecayPi2/D");
  tree->Branch("decayPiMomCal2",   &event.decayPiMomCal2,  "decayPiMomCal2/D");

  tree->Branch("vertexNpScat",   event.vertexNpScat,  "vertexNpScat[3]/D");
  tree->Branch("cdistNpScat",   &event.cdistNpScat,  "cdistNpScat/D");
  tree->Branch("thetaNpScat",   &event.thetaNpScat,  "thetaNpScat/D");
  tree->Branch("scatNpMomCal",   &event.scatNpMomCal,  "scatNpMomCal/D");
  tree->Branch("scatNpEkinCal",   &event.scatNpEkinCal,  "scatNpEkinCal/D");
  tree->Branch("thetaNpScatCM",   &event.thetaNpScatCM,  "thetaNpScatCM/D");

  tree->Branch("cdistLambdaDecay",   &event.cdistLambdaDecay,  "cdistLambdaDecay/D");
  tree->Branch("thetaLambdaDecay",   &event.thetaLambdaDecay,  "thetaLambdaDecay/D");
  tree->Branch("vertexLambdaDecay",   event.vertexLambdaDecay,  "vertexLambdaDecay[3]/D");
  tree->Branch("momPiFromLambda",   &event.momPiFromLambda,  "momPiFromLambda/D");
  tree->Branch("momProtonFromLambda",   &event.momProtonFromLambda,  "momProtonFromLambda/D");
  tree->Branch("momLambda",   &event.momLambda,  "momLambda/D");
  tree->Branch("momVecPiFromLambda",   event.momVecPiFromLambda,  "momVecPiFromLambda/D");
  tree->Branch("momVecProtonFromLambda",   event.momVecProtonFromLambda,  "momVecProtonFromLambda/D");
  tree->Branch("momVecLambda",   event.momVecLambda,  "momVecLambda/D");
  tree->Branch("cdistLambdaNConv",   &event.cdistLambdaNConv,  "cdistLambdaNConv/D");
  tree->Branch("thetaLambdaNConv",   &event.thetaLambdaNConv,  "thetaLambdaNConv/D");
  tree->Branch("vertexLambdaNConv",   event.vertexLambdaNConv,  "vertexLambdaNConv[3]/D");
  tree->Branch("momCalLambda",   &event.momCalLambda,  "momCalLambda/D");
  tree->Branch("momCalLambda2",   &event.momCalLambda2,  "momCalLambda2/D");
  tree->Branch("thetaCMLambdaNConv",   &event.thetaCMLambdaNConv,  "thetaCMLambdaNConv/D");

  tree->Branch("vertexPiPScat",   event.vertexPiPScat,  "vertexPiPScat[3]/D");
  tree->Branch("cdistPiPScat",   &event.cdistPiPScat,  "cdistPiPScat/D");

  tree->Branch("piBeamMom",   &event.piBeamMom,  "piBeamMom/D");
  tree->Branch("thetaPiPScat",   &event.thetaPiPScat,  "thetaPiPScat/D");
  tree->Branch("scatPiMomCal",   &event.scatPiMomCal,  "scatPiMomCal/D");
  tree->Branch("scatPiMomCalFromE",   &event.scatPiMomCalFromE,  "scatPiMomCalFromE/D");

  tree->Branch("ntCFT",   &event.ntCFT,  "ntCFT/I");

  tree->Branch("xDirCFT",   event.xDirCFT,  "xDirCFT[ntCFT]/D");
  tree->Branch("yDirCFT",   event.yDirCFT,  "yDirCFT[ntCFT]/D");
  tree->Branch("zDirCFT",   event.zDirCFT,  "zDirCFT[ntCFT]/D");
  tree->Branch("xPos0CFT",   event.xPos0CFT,  "xPos0CFT[ntCFT]/D");
  tree->Branch("yPos0CFT",   event.yPos0CFT,  "yPos0CFT[ntCFT]/D");
  tree->Branch("zPos0CFT",   event.zPos0CFT,  "zPos0CFT[ntCFT]/D");
  tree->Branch("zBGOCFT",   event.zBGOCFT,  "zBGOCFT[ntCFT]/D");

  tree->Branch("BGO_Edep",   event.BGO_Edep,  "BGO_Edep[ntCFT]/D");
  tree->Branch("TotalEdep",   event.TotalEdep,  "TotalEdep[ntCFT]/D");
  tree->Branch("CFT_TotalEdep",   event.CFT_TotalEdep,  "CFT_TotalEdep[ntCFT]/D");
  tree->Branch("CFT_NormTotalEdep",   event.CFT_NormTotalEdep,  "CFT_NormTotalEdep[ntCFT]/D");
  tree->Branch("PiV_Edep",   event.PiV_Edep,  "PiV_Edep[ntCFT]/D");

  tree->Branch("xDirCFT_P",   event.xDirCFT_P,  "xDirCFT_P[nP_CFT]/D");
  tree->Branch("yDirCFT_P",   event.yDirCFT_P,  "yDirCFT_P[nP_CFT]/D");
  tree->Branch("zDirCFT_P",   event.zDirCFT_P,  "zDirCFT_P[nP_CFT]/D");
  tree->Branch("xPos0CFT_P",   event.xPos0CFT_P,  "xPos0CFT_P[nP_CFT]/D");
  tree->Branch("yPos0CFT_P",   event.yPos0CFT_P,  "yPos0CFT_P[nP_CFT]/D");
  tree->Branch("zPos0CFT_P",   event.zPos0CFT_P,  "zPos0CFT_P[nP_CFT]/D");
  tree->Branch("zBGOCFT_P",   event.zBGOCFT_P,  "zBGOCFT_P[nP_CFT]/D");

  tree->Branch("CFTVtx_x_P",   event.CFTVtx_x_P,  "CFTVtx_x_P[nP_CFT]/D");
  tree->Branch("CFTVtx_y_P",   event.CFTVtx_y_P,  "CFTVtx_y_P[nP_CFT]/D");
  tree->Branch("CFTVtx_z_P",   event.CFTVtx_z_P,  "CFTVtx_z_P[nP_CFT]/D");

  tree->Branch("BGO_Edep_P",   event.BGO_Edep_P,  "BGO_Edep_P[nP_CFT]/D");
  tree->Branch("TotalEdep_P",   event.TotalEdep_P,  "TotalEdep_P[nP_CFT]/D");
  tree->Branch("CFT_TotalEdep_P",   event.CFT_TotalEdep_P,  "CFT_TotalEdep_P[nP_CFT]/D");
  tree->Branch("CFT_NormTotalEdep_P",   event.CFT_NormTotalEdep_P,  "CFT_NormTotalEdep_P[nP_CFT]/D");
  tree->Branch("PiV_Edep_P",   event.PiV_Edep_P,  "PiV_Edep_P[nP_CFT]/D");

  tree->Branch("xDirCFT_Pi",   event.xDirCFT_Pi,  "xDirCFT_Pi[nPi_CFT]/D");
  tree->Branch("yDirCFT_Pi",   event.yDirCFT_Pi,  "yDirCFT_Pi[nPi_CFT]/D");
  tree->Branch("zDirCFT_Pi",   event.zDirCFT_Pi,  "zDirCFT_Pi[nPi_CFT]/D");
  tree->Branch("xPos0CFT_Pi",   event.xPos0CFT_Pi,  "xPos0CFT_Pi[nPi_CFT]/D");
  tree->Branch("yPos0CFT_Pi",   event.yPos0CFT_Pi,  "yPos0CFT_Pi[nPi_CFT]/D");
  tree->Branch("zPos0CFT_Pi",   event.zPos0CFT_Pi,  "zPos0CFT_Pi[nPi_CFT]/D");
  tree->Branch("zBGOCFT_Pi",   event.zBGOCFT_Pi,  "zBGOCFT_Pi[nPi_CFT]/D");

  tree->Branch("CFTVtx_x_Pi",   event.CFTVtx_x_Pi,  "CFTVtx_x_Pi[nPi_CFT]/D");
  tree->Branch("CFTVtx_y_Pi",   event.CFTVtx_y_Pi,  "CFTVtx_y_Pi[nPi_CFT]/D");
  tree->Branch("CFTVtx_z_Pi",   event.CFTVtx_z_Pi,  "CFTVtx_z_Pi[nPi_CFT]/D");

  tree->Branch("u0CFT_Pi",   event.u0CFT_Pi,  "u0CFT_Pi[nPi_CFT]/D");
  tree->Branch("v0CFT_Pi",   event.v0CFT_Pi,  "v0CFT_Pi[nPi_CFT]/D");

  tree->Branch("BGO_Edep_Pi",   event.BGO_Edep_Pi,  "BGO_Edep_Pi[nPi_CFT]/D");
  tree->Branch("TotalEdep_Pi",   event.TotalEdep_Pi,  "TotalEdep_Pi[nPi_CFT]/D");
  tree->Branch("CFT_TotalEdep_Pi",   event.CFT_TotalEdep_Pi,  "CFT_TotalEdep_Pi[nPi_CFT]/D");
  tree->Branch("CFT_NormTotalEdep_Pi",   event.CFT_NormTotalEdep_Pi,  "CFT_NormTotalEdep_Pi[nPi_CFT]/D");
  tree->Branch("PiV_Edep_Pi",   event.PiV_Edep_Pi,  "PiV_Edep_Pi[nPi_CFT]/D");


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
