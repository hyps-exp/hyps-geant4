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
#include <utility>
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
const double XiMinusMass     = 1.32131;
////////////////////////////////////////////////////

const double MinOfMassSquareK = 0.1;
const double MaxOfMassSquareK = 0.4;

const double MinOfMassSigma = 1.15;
const double MaxOfMassSigma = 1.25;

const double MinOfMassXiMinus = 1.290;
const double MaxOfMassXiMinus = 1.370;

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

  int nXi;
  double pKCal[MaxHits];
  double theta[MaxHits];
  double pXiCal[MaxHits];
  double Vertex_x[MaxHits];
  double Vertex_y[MaxHits];
  double Vertex_z[MaxHits];

  int nPi_CFT;
  int nK_CFT;
  int nP_CFT;
  double EkinP;
  double EkinCorP;
  double EkinP_2nd;
  double EkinCorP_2nd;

  double vertexDecayPi[3];
  double cdistDecayPi;
  double thetaDecayPi;
  double decayPiMomCal;	     
  double decayPiMomVecCal[3];
  double decayLambdaMomCal;	     
  double decayLambdaMomVecCal[3];


  double vertexScat[3];
  double hypBeamVec[3];
  double cdistScat;
  double thetaScat;
  double scatMomCal;
  double scatEkinCal;
  double thetaScatCM;

  double thetaDecayPi_Xi1;
  double decayPiMomCal_Xi1;	     
  double cdistDecayPi_Xi1;
  double vertexDecayPi_Xi1[3];

  double thetaDecayPi_Xi2;
  double decayPiMomCal_Xi2;	     
  double cdistDecayPi_Xi2;
  double vertexDecayPi_Xi2[3];

  double vertexScat_2nd[3];
  double hypBeamVec_2nd[3];
  double cdistScat_2nd;
  double thetaScat_2nd;
  double scatMomCal_2nd;
  double scatEkinCal_2nd;
  double thetaScatCM_2nd;

  double thetaDecayPi_2nd_Xi1;
  double decayPiMomCal_2nd_Xi1;	     
  double cdistDecayPi_2nd_Xi1;
  double vertexDecayPi_2nd_Xi1[3];

  double thetaDecayPi_2nd_Xi2;
  double decayPiMomCal_2nd_Xi2;	     
  double cdistDecayPi_2nd_Xi2;
  double vertexDecayPi_2nd_Xi2[3];


  double thetaDecayPi2;
  double decayPiMomCal2;	     
  double cdistDecayPi2;
  double vertexDecayPi2[3];

  double vertexNpScat[3];
  double cdistNpScat;
  double thetaNpScat;
  double scatNpMomCal;
  double scatNpEkinCal;
  double thetaNpScatCM;

  int    nLambda;
  double cdistLambdaDecay;
  double vertexLambdaDecay[3];
  double thetaLambdaDecay;
  double momPiFromLambda;
  double momProtonFromLambda;
  double momLambda;
  double momVecPiFromLambda[3];
  double momVecProtonFromLambda[3];
  double momVecLambda[3];
  double cdistLambdaLambdaConv;
  double thetaLambdaLambdaConv;
  double vertexLambdaLambdaConv[3];
  double momCalLambda;
  double momCalLambda2;
  double thetaCMLambdaLambdaConv;
  double invMass_PiP;

  double cdistLambdaDecay_2nd;
  double vertexLambdaDecay_2nd[3];
  double thetaLambdaDecay_2nd;
  double momPiFromLambda_2nd;
  double momProtonFromLambda_2nd;
  double momLambda_2nd;
  double momVecPiFromLambda_2nd[3];
  double momVecProtonFromLambda_2nd[3];
  double momVecLambda_2nd[3];
  double cdistLambdaLambdaConv_2nd;
  double thetaLambdaLambdaConv_2nd;
  double vertexLambdaLambdaConv_2nd[3];
  double momCalLambda_2nd;
  double momCalLambda2_2nd;
  double thetaCMLambdaLambdaConv_2nd;
  double invMass_PiP_2nd;

  double momPiFromXi;
  double momXi;
  double momVecPiFromXi[3];
  double momVecXi[3];
  double cdistLambdaPi;
  double thetaLambdaPi;
  double vertexLambdaPi[3];

  double cdistXiPScat_Xi;
  double vertexXiPScat_Xi[3];
  double costXiPScat_Xi;
  double thetaXiPScat_Xi;
  double thetaCMXiPScat_Xi;
  double momCalXi1;
  double momCalXi2;
  //double momCalRecoilProton;
  //double momRecoilProton;
  //double deltaThetaRecoilProton;

  double EkinP_Lp;			
  double EkinCorP_Lp;			
  double EkinP_2nd_Lp;			
  double EkinCorP_2nd_Lp;		
  double invMass_PiP_Lp;		
  double momLambda1_Lp;			
  double momLambda2_Lp;			
  double vertexLambdaDecay_Lp[3];	
  double cdistLambdaDecay_Lp;		
  double cdistLambdaPScat_Lambda;	
  double vertexLambdaPScat_Lambda[3];	
  double costLambdaPScat_Lambda;	
  double thetaLambdaPScat_Lambda;	
  double thetaCMLambdaPScat_Lambda;	
  double momCalLambda1_Lp;		
  double momCalLambda2_Lp;		
  double momCalRecoilProton_Lp;		
  double momRecoilProton_Lp;		
  double deltaThetaRecoilProton_Lp;	
  double vertexLambdaPScat_p[3];	
  double cdistLambdaPScat_p;		
  double thetaScatLambdaPScat_p;	
  double scatMomCalLambdaPScat_p;	
  double scatEkinCalLambdaPScat_p;	
  double thetaScatCMLambdaPScat_p;	
  double costScatOpenAngle_Lp;		
  double thetaScatOpenAngle_Lp;		
  double cdistDecayPi2_Lp;		
  double vertexDecayPi2_Lp[3];		
  double cdistDecayProton1_Lp_2p;	
  double vertexDecayProton1_Lp_2p[3];	
  double MissMassPi_1_Lp_2p;            

  double momLambda1_Lp_1st;		     
  double momLambda2_Lp_1st;		     
  double vertexLambdaDecay_Lp_1st[3];	     
  double cdistLambdaDecay_Lp_1st;	     
  double cdistLambdaPScat_Lambda_1st;	     
  double vertexLambdaPScat_Lambda_1st[3];    
  double costLambdaPScat_Lambda_1st;	     
  double thetaLambdaPScat_Lambda_1st;	     
  double thetaCMLambdaPScat_Lambda_1st;	     
  double momCalLambda1_Lp_1st;		     
  double momCalLambda2_Lp_1st;		     
  double momCalRecoilProton_Lp_1st;	     
  double momRecoilProton_Lp_1st;	     
  double deltaThetaRecoilProton_Lp_1st;	     
	 					     
  double momLambda1_Lp_2nd;		     
  double momLambda2_Lp_2nd;		     
  double vertexLambdaDecay_Lp_2nd[3];	     
  double cdistLambdaDecay_Lp_2nd;	     
  double cdistLambdaPScat_Lambda_2nd;	     
  double vertexLambdaPScat_Lambda_2nd[3];    
  double costLambdaPScat_Lambda_2nd;	     
  double thetaLambdaPScat_Lambda_2nd;	     
  double thetaCMLambdaPScat_Lambda_2nd;	     
  double momCalLambda1_Lp_2nd;		     
  double momCalLambda2_Lp_2nd;		     
  double momCalRecoilProton_Lp_2nd;	     
  double momRecoilProton_Lp_2nd;	     
  double deltaThetaRecoilProton_Lp_2nd;      


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

  std::vector <PartVec> XiVecCont;
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

      LorentzVector LvBeam(  bMom, sqrt(KaonMass*KaonMass+bMom.mag2()) );
      LorentzVector LvTgt( 0., 0., 0., ProtonMass );
      LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
      LorentzVector LvRc = LvBeam+LvTgt-LvScat;
      double mismass = LvRc.mag();
      HF1(67, mismass);
      event.MissMass[nK] = mismass;
      event.pK[nK] = priMom.mag();

      nK++;

      if (mismass>MinOfMassXiMinus && mismass<MaxOfMassXiMinus) {
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
	ThreeVector pXiCor = bMom-KaonMomCor;

	double calMom = calcKMomFromTheta(event.beammom, cost, priMom.mag());
	event.pKCal[index] = calMom;
	ThreeVector KaonMomCal = priMom*(calMom/priMom.mag());
	ThreeVector pXiCal = bMom-KaonMomCal;
	event.pXiCal[index] = pXiCal.mag();

	PartVec xiVec;
	xiVec.pos0 = vtx;
	xiVec.mom  = pXiCal; // original
	//xiVec.mom  = pXiCal*event.momHypBeam/pXiCal.mag();
	//xiVec.mom  = ThreeVector(event.momVectorHypBeam[0],event.momVectorHypBeam[1],event.momVectorHypBeam[2]);
	//xiVec.mom  = pXiCor;
	XiVecCont.push_back(xiVec);
	
	index++;
      }
    }
  }

  event.nK = nK;

  int nXi = XiVecCont.size();
  event.nXi = nXi;

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

  if (nXi != 1) {
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

  PartVec XiVec = XiVecCont[0];
  event.hypBeamVec[0] = XiVec.mom.x();
  event.hypBeamVec[1] = XiVec.mom.y();
  event.hypBeamVec[2] = XiVec.mom.z();

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
  PartVec DecayLambdaVec;

  PartVec XiVecOrg;
  int indexXiDecayPi = -1;

  // Xi decay assumption
  if (nPi_CFT==2 || nPi_CFT==1) {
    for (int i1=0; i1<nPi_CFT; i1++) {
      CFTParticle *PiCFT = CFTPionCont[i1];


      PiDecayVec.pos0 = PiCFT->GetPos0();
      PiDecayVec.mom = PiCFT->GetDir();
      
      double cdistDecayPi;
      ThreeVector VertXiDecay = VertexPoint3D(XiVec.pos0, PiDecayVec.pos0,
					      XiVec.mom, PiDecayVec.mom, 
					      cdistDecayPi);
      double costDecayPi=XiVec.mom*PiDecayVec.mom/(XiVec.mom.mag()*PiDecayVec.mom.mag());
      double thetaDecayPi=acos(costDecayPi)*Rad2Deg;

      double decayPiMomCal;
      
      flagDecayCal = calcDecayPiMom( XiVec.mom.mag(), 
				     XiMinusMass, 
				     LambdaMass,
				     PionMass,
				     costDecayPi,
				     &decayPiMomCal);
      
      if (flagDecayCal) {
	if (i1 == 0) {
	  indexXiDecayPi =i1;
	  event.vertexDecayPi[0] = VertXiDecay.x();
	  event.vertexDecayPi[1] = VertXiDecay.y();
	  event.vertexDecayPi[2] = VertXiDecay.z();
	  
	  event.cdistDecayPi = cdistDecayPi;
	  event.thetaDecayPi = thetaDecayPi;
	  
	  event.decayPiMomCal = decayPiMomCal;      
	  event.decayPiMomVecCal[0] = decayPiMomCal*PiDecayVec.mom.x()/PiDecayVec.mom.mag();
	  event.decayPiMomVecCal[1] = decayPiMomCal*PiDecayVec.mom.y()/PiDecayVec.mom.mag();
	  event.decayPiMomVecCal[2] = decayPiMomCal*PiDecayVec.mom.z()/PiDecayVec.mom.mag();
	  
	  PiDecayVec.mom *= (decayPiMomCal/PiDecayVec.mom.mag());
	  
	  DecayLambdaVec.mom = XiVec.mom-PiDecayVec.mom;
	  DecayLambdaVec.pos0 = VertXiDecay;
	  
	  event.decayLambdaMomCal = DecayLambdaVec.mom.mag();
	  event.decayLambdaMomVecCal[0] = DecayLambdaVec.mom.x();
	  event.decayLambdaMomVecCal[1] = DecayLambdaVec.mom.y();
	  event.decayLambdaMomVecCal[2] = DecayLambdaVec.mom.z();
	} else if (i1==1) {
	  if (std::abs(VertXiDecay.z() - event.Vertex_z[0]) < std::abs(event.vertexDecayPi[2] - event.Vertex_z[0])) {
	    indexXiDecayPi =i1;
	    event.vertexDecayPi[0] = VertXiDecay.x();
	    event.vertexDecayPi[1] = VertXiDecay.y();
	    event.vertexDecayPi[2] = VertXiDecay.z();
	  
	    event.cdistDecayPi = cdistDecayPi;
	    event.thetaDecayPi = thetaDecayPi;
	  
	    event.decayPiMomCal = decayPiMomCal;      
	    event.decayPiMomVecCal[0] = decayPiMomCal*PiDecayVec.mom.x()/PiDecayVec.mom.mag();
	    event.decayPiMomVecCal[1] = decayPiMomCal*PiDecayVec.mom.y()/PiDecayVec.mom.mag();
	    event.decayPiMomVecCal[2] = decayPiMomCal*PiDecayVec.mom.z()/PiDecayVec.mom.mag();
	  
	    PiDecayVec.mom *= (decayPiMomCal/PiDecayVec.mom.mag());
	    
	    DecayLambdaVec.mom = XiVec.mom-PiDecayVec.mom;
	    DecayLambdaVec.pos0 = VertXiDecay;
	    
	    event.decayLambdaMomCal = DecayLambdaVec.mom.mag();
	    event.decayLambdaMomVecCal[0] = DecayLambdaVec.mom.x();
	    event.decayLambdaMomVecCal[1] = DecayLambdaVec.mom.y();
	    event.decayLambdaMomVecCal[2] = DecayLambdaVec.mom.z();
	  }
	}
      }
    }
  }

  bool flagXiPScat = false;
  PartVec ProtonScatVec;
  bool flagDecayCal2 = false;
  bool flagNpScat = false;

  bool flagLambdaInv = false;
  bool flagLambdaLambdaConv = false;
  PartVec LambdaVec;

  std::vector <double> cdistPiPCont;
  typedef std::pair<int, int> pair;
  std::vector <pair> indexPiPCont;
  
  if ((nP_CFT==2 || nP_CFT==1) && (nPi_CFT==2 || nPi_CFT==1)) {
    for (int i1=0; i1<nP_CFT; i1++) {
      for (int i2=0; i2<nPi_CFT; i2++) {

        CFTParticle *ProtonCFT = CFTProtonCont[i1];
        ThreeVector DirProton = ProtonCFT->GetDir();
        ThreeVector PosProton = ProtonCFT->GetPos0();

        CFTParticle *PiCFT = CFTPionCont[i2];
        ThreeVector DirPi = PiCFT->GetDir();
        ThreeVector PosPi = PiCFT->GetPos0();

	double cdist;
	ThreeVector vertex = VertexPoint3D( PosProton, PosPi, DirProton, DirPi, cdist );

	if (cdist>5)
	  continue;


	indexPiPCont.push_back(pair(i1, i2));
	cdistPiPCont.push_back(cdist);

      }
    }
    

    std::vector <pair> indexPiPContOK;
    if (indexPiPCont.size()>=2) {
      // 1st vertex
      double minCdist=999.;
      int index1=-1;
      for (int i=0; i<cdistPiPCont.size(); i++) {
	if (cdistPiPCont[i] < minCdist) {
	  minCdist = cdistPiPCont[i];
	  index1 = i;
	}
      }

      // 2nd vertex
      minCdist=999.;
      int index2=-1;
      for (int i=0; i<cdistPiPCont.size(); i++) {
	if (i==index1)
	  continue;

	if (cdistPiPCont[i] < minCdist &&
	    indexPiPCont[i].first != indexPiPCont[index1].first &&
	    indexPiPCont[i].second != indexPiPCont[index1].second
	    ) {
	  minCdist = cdistPiPCont[i];
	  index2 = i;
	}
      }

      if (index1>=0 && index2>=0) {
	// 2 pip vertex
	indexPiPContOK.push_back(indexPiPCont[index1]);
	indexPiPContOK.push_back(indexPiPCont[index2]);
      } else if (index1>=0) {
	// 1 pip vertex
	indexPiPContOK.push_back(indexPiPCont[index1]);
      }
    } else if (indexPiPCont.size()==1) {
      // 1pip vertex
      indexPiPContOK.push_back(indexPiPCont[0]);
    }

    event.nLambda = indexPiPContOK.size();
    
    for (int i=0; i<indexPiPContOK.size(); i++) {
      int indexP  = indexPiPContOK[i].first;
      int indexPi = indexPiPContOK[i].second;

      CFTParticle *PiCFT = CFTPionCont[indexPi];
      ThreeVector DirPi = PiCFT->GetDir();
      ThreeVector PosPi = PiCFT->GetPos0();

      CFTParticle *ProtonCFT = CFTProtonCont[indexP];
      ThreeVector DirProton = ProtonCFT->GetDir();
      ThreeVector PosProton = ProtonCFT->GetPos0();

      double cdistLambdaDecay;
      ThreeVector VertLambdaDecay = VertexPoint3D( PosPi, PosProton, DirPi, DirProton, cdistLambdaDecay );

      double costLambdaDecay=DirPi*DirProton/(DirPi.mag()*DirProton.mag());
      double thetaLambdaDecay=acos(costLambdaDecay)*Rad2Deg;

      double totE_pi = PiCFT->GetTotalE();
      double p_pi = sqrt(totE_pi*totE_pi+2.*totE_pi*PionMass*1000.);
      p_pi /= 1000.; // GeV/c  
      ThreeVector MomPi1 = DirPi*(p_pi/DirPi.mag());

      
      double Ekin1 = ProtonCFT->GetTotalE();

      Ekin1 *= 0.001; // GeV
      double M1 = ProtonMass;
      double p1 = sqrt(Ekin1*Ekin1+2.*Ekin1*M1) ;
      double p_cor, e_cor;

      /*
	CorrElossOut(&p_cor, &e_cor, p1, PROTON, ProtonScatVec.mom, VertLambdaDecay);
      */

      CorrElossOutWithCFRP(&p_cor, &e_cor, p1, PROTON, 
			   DirProton/DirProton.mag(), VertLambdaDecay,
			   PosProton);

      p1 = p_cor;
      ThreeVector momVecProton = p1*DirProton/DirProton.mag();

      LorentzVector LvPi( MomPi1, 
			  std::sqrt( PionMass*PionMass+MomPi1.mag2()) );
      LorentzVector LvP( momVecProton, 
			 std::sqrt( ProtonMass*ProtonMass+momVecProton.mag2()) );

      LorentzVector LvPiP = LvPi + LvP;
      double InvMass = LvPiP.mag();

      if (i==0) {
	event.cdistLambdaDecay = cdistLambdaDecay;
	event.thetaLambdaDecay = thetaLambdaDecay;
	event.vertexLambdaDecay[0] = VertLambdaDecay.x();
	event.vertexLambdaDecay[1] = VertLambdaDecay.y();
	event.vertexLambdaDecay[2] = VertLambdaDecay.z();
	event.invMass_PiP = InvMass;
      } else if (i==1) {
	event.cdistLambdaDecay_2nd = cdistLambdaDecay;
	event.thetaLambdaDecay_2nd = thetaLambdaDecay;
	event.vertexLambdaDecay_2nd[0] = VertLambdaDecay.x();
	event.vertexLambdaDecay_2nd[1] = VertLambdaDecay.y();
	event.vertexLambdaDecay_2nd[2] = VertLambdaDecay.z();
	event.invMass_PiP_2nd = InvMass;
      }


      double momPiFromLambda;

      flagLambdaInv = calcPiMomFromLambda(LambdaMass, ProtonMass, p1,
					  PionMass, costLambdaDecay, &momPiFromLambda);

      if (flagLambdaInv) {
	ThreeVector momVecPi     = momPiFromLambda*DirPi/DirPi.mag();
	ThreeVector momVecLambda = momVecPi+momVecProton;

	LambdaVec.mom = momVecLambda;
	LambdaVec.pos0 = VertLambdaDecay;

	/*
	  std::cout << "P_pi = "<< momVecPi.mag() << ", P_proton = " << momVecProton.mag() 
	  << ", P_lambda = " << momVecLambda.mag() << std::endl;
	*/



	double cdistLambdaLambdaConv;
	ThreeVector VertLambdaLambdaConv = 
	  VertexPoint3D( XiVec.pos0, LambdaVec.pos0,
			 XiVec.mom, LambdaVec.mom, cdistLambdaLambdaConv);
	
	double costLambdaLambdaConv=XiVec.mom*LambdaVec.mom/
	  (XiVec.mom.mag()*LambdaVec.mom.mag());
	double thetaLambdaLambdaConv=acos(costLambdaLambdaConv)*Rad2Deg;

	if (i==0) {
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
	  
	  event.cdistLambdaLambdaConv = cdistLambdaLambdaConv;
	  event.thetaLambdaLambdaConv = thetaLambdaLambdaConv;
	  event.vertexLambdaLambdaConv[0] = VertLambdaLambdaConv.x();
	  event.vertexLambdaLambdaConv[1] = VertLambdaLambdaConv.y();
	  event.vertexLambdaLambdaConv[2] = VertLambdaLambdaConv.z();
  
	} else if (i==1) {
	  event.momPiFromLambda_2nd = momVecPi.mag();
	  event.momProtonFromLambda_2nd = momVecProton.mag();
	  event.momLambda_2nd = momVecLambda.mag();
	  event.momVecPiFromLambda_2nd[0] = momVecPi.x();
	  event.momVecPiFromLambda_2nd[1] = momVecPi.y();
	  event.momVecPiFromLambda_2nd[2] = momVecPi.z();
	  event.momVecProtonFromLambda_2nd[0] = momVecProton.x();
	  event.momVecProtonFromLambda_2nd[1] = momVecProton.y();
	  event.momVecProtonFromLambda_2nd[2] = momVecProton.z();
	  event.momVecLambda_2nd[0] = momVecLambda.x();
	  event.momVecLambda_2nd[1] = momVecLambda.y();
	  event.momVecLambda_2nd[2] = momVecLambda.z();
	  
	  event.cdistLambdaLambdaConv_2nd = cdistLambdaLambdaConv;
	  event.thetaLambdaLambdaConv_2nd = thetaLambdaLambdaConv;
	  event.vertexLambdaLambdaConv_2nd[0] = VertLambdaLambdaConv.x();
	  event.vertexLambdaLambdaConv_2nd[1] = VertLambdaLambdaConv.y();
	  event.vertexLambdaLambdaConv_2nd[2] = VertLambdaLambdaConv.z();
	}


	// Xi- decay (Xi-p scattering)
	if (indexXiDecayPi>=0 && indexPi != indexXiDecayPi) {
	  // Xi- decay check
	  double cdistLambdaPi;
	  CFTParticle *PiCFT2 = CFTPionCont[indexXiDecayPi];
	  ThreeVector DirPi2 = PiCFT2->GetDir();
	  ThreeVector PosPi2 = PiCFT2->GetPos0();
	  
	  ThreeVector VertLambdaPi = 
	    VertexPoint3D( PosPi2, LambdaVec.pos0,
			   DirPi2, LambdaVec.mom, cdistLambdaPi);
	  
	  double costLambdaPi=DirPi2*LambdaVec.mom/(DirPi2.mag()*LambdaVec.mom.mag());
	  double thetaLambdaPi=acos(costLambdaPi)*Rad2Deg;
	  
	  double momPiFromXi;
	  bool flagXiInv = calcPiMomFromLambda(XiMinusMass, LambdaMass, LambdaVec.mom.mag(),
					       PionMass, costLambdaPi, &momPiFromXi);
	  if (flagXiInv) {
	    ThreeVector momVecPi2 = momPiFromXi*DirPi2/DirPi2.mag();
	    ThreeVector momVecXi  = momVecPi2+LambdaVec.mom;
	    
	    event.momPiFromXi = momVecPi2.mag();
	    event.momXi       = momVecXi.mag();
	    event.momVecPiFromXi[0] = momVecPi2.x();
	    event.momVecPiFromXi[1] = momVecPi2.y();
	    event.momVecPiFromXi[2] = momVecPi2.z();
	    event.momVecXi[0] = momVecXi.x();
	    event.momVecXi[1] = momVecXi.y();
	    event.momVecXi[2] = momVecXi.z();
	    
	    event.cdistLambdaPi = cdistLambdaPi;
	    event.thetaLambdaPi = thetaLambdaPi;
	    event.vertexLambdaPi[0] = VertLambdaPi.x();
	    event.vertexLambdaPi[1] = VertLambdaPi.y();
	    event.vertexLambdaPi[2] = VertLambdaPi.z();


	    // Xi-p scattering check

	    bool flagXiPScat_Xi = false;	    
            double cdistXiPScat_Xi;
            ThreeVector VertXiPScat_Xi =
              VertexPoint3D( XiVec.pos0, VertLambdaPi,
			     XiVec.mom, momVecXi, cdistXiPScat_Xi);

            double costXiPScat_Xi=XiVec.mom*momVecXi/
              (XiVec.mom.mag()*momVecXi.mag());
            double thetaXiPScat_Xi=acos(costXiPScat_Xi)*Rad2Deg;

            double momCalXi1 =-999.;
            double momCalXi2 =-999.;
            double thetaCMXiPScat_Xi =-999.;

            flagXiPScat_Xi =
              calc2BodyInelastic(XiMinusMass, XiVec.mom.mag(),
				 ProtonMass, XiMinusMass, ProtonMass,
				 thetaXiPScat_Xi,
				 &momCalXi1, &momCalXi2,
				 &thetaCMXiPScat_Xi);

            /* 
	       std::cout << "Xi mom : " << XiVec.mom.mag()
	       << ", momCalXi : "  << momCalXi1
	       << ", momCalXi2 : "  << momCalXi2
	       << ", thetaCMXiPScat_Xi "<< thetaCMXiPScat_Xi << std::endl;
	    */

            event.cdistXiPScat_Xi = cdistXiPScat_Xi;
            event.vertexXiPScat_Xi[0] = VertXiPScat_Xi.x();
            event.vertexXiPScat_Xi[1] = VertXiPScat_Xi.y();
            event.vertexXiPScat_Xi[2] = VertXiPScat_Xi.z();

            event.costXiPScat_Xi = costXiPScat_Xi;
            event.thetaXiPScat_Xi = thetaXiPScat_Xi;
            event.thetaCMXiPScat_Xi = thetaCMXiPScat_Xi;
            event.momCalXi1 = momCalXi1;
            event.momCalXi2 = momCalXi2;

	  }
	}

	  
	double momCalLambda =-999.;
	double momCalLambda2 =-999.;
	double thetaCMLambdaLambdaConv =-999.;


	flagLambdaLambdaConv = calc2BodyInelastic(XiMinusMass, XiVec.mom.mag(),
					     ProtonMass, LambdaMass, LambdaMass,
					     thetaLambdaLambdaConv,
					     &momCalLambda, &momCalLambda2,
					     &thetaCMLambdaLambdaConv);
	/*
	  std::cout << "SigmaMinusMass : "  << SigmaMinusMass
	  << ", Sigma mom : " << SigmaVec.mom.mag()
	  << ", ProtonMass : "  << ProtonMass 
	  << ", LambdaMass : "  << LambdaMass 
	  << ", NeutronMass : "  << NeutronMass 
	  << ", thetaLambdaLambdaConv : "  << thetaLambdaLambdaConv
	  << ", momCalLambda : "  << momCalLambda
	  << ", momCalLambda2 : "  << momCalLambda2
	  << ", thetaCMLambdaLambda "<< thetaCMLambdaLambdaConv << std::endl;
	*/
	if (flagLambdaLambdaConv) {
	  if (i==0) {
	    event.momCalLambda = momCalLambda;
	    event.momCalLambda2 = momCalLambda2;
	    event.thetaCMLambdaLambdaConv = thetaCMLambdaLambdaConv;
	  } else if (i==1) {
	    event.momCalLambda_2nd = momCalLambda;
	    event.momCalLambda2_2nd = momCalLambda2;
	    event.thetaCMLambdaLambdaConv_2nd = thetaCMLambdaLambdaConv;
	  }
	}
      }


    }

  }

  if ((nP_CFT==2 || nP_CFT==1) && (nPi_CFT==2 || nPi_CFT==1)) {
    for (int i1=0; i1<nP_CFT; i1++) {
      CFTParticle *ProtonCFT = CFTProtonCont[i1];
      ThreeVector DirProton = ProtonCFT->GetDir();
      ThreeVector PosProton = ProtonCFT->GetPos0();

      /* Xi- p scattering assumption */
      double cdist1;

    //ThreeVector vtx0(event.primaryVertex[0], event.primaryVertex[1], event.primaryVertex[2]);
      ThreeVector vtx(event.Vertex_x[0], event.Vertex_y[0], event.Vertex_z[0]);

      ThreeVector VertScat1 = VertexPoint3D( XiVec.pos0, PosProton,
					     XiVec.mom, DirProton, cdist1);


      /* Xi beam energy correction */
      PartVec XiVec1 = XiVec;
      ThreeVector HypBeamLength1 = VertScat1-vtx;
      double pXiBeam0 = XiVec1.mom.mag();
      double pXiBeamCor, eXiBeamCor;
      int LH2Tag=0;

      caldE(pXiBeam0,  // GeV/c
	    XiMinusMass, // GeV
	    HypBeamLength1.mag(), // mm
	    &pXiBeamCor, // GeV/c
	    &eXiBeamCor, // GeV
	    LH2Tag);

      XiVec1.mom = XiVec1.mom*pXiBeamCor/XiVec1.mom.mag();

      double hypBeamVec1[3];
      hypBeamVec1[0] = XiVec1.mom.x();
      hypBeamVec1[1] = XiVec1.mom.y();
      hypBeamVec1[2] = XiVec1.mom.z();


      double costScat1=XiVec1.mom*DirProton/(XiVec1.mom.mag()*DirProton.mag());
      double thetaScat1=acos(costScat1)*Rad2Deg;

      double vertexScat1[3];
      vertexScat1[0] = VertScat1.x();
      vertexScat1[1] = VertScat1.y();
      vertexScat1[2] = VertScat1.z();

      double scatMomCal1 = -999.;
      double scatEkinCal1 = -999.;
      double thetaScatCM1 = -999.;
      double EkinCorP1 = -999.;

      
      if (thetaScat1>0. && thetaScat1<90.) {

	flagXiPScat = calc2BodyKinema(XiMinusMass, XiVec1.mom.mag(),  
				      ProtonMass, thetaScat1,
				      &scatMomCal1, &scatEkinCal1, &thetaScatCM1);

	if (flagXiPScat) {
	  scatEkinCal1 *= 1000;
	  /*
	    event.scatMomCal = scatMomCal;
	    event.scatEkinCal = scatEkinCal*1000;
	    event.thetaScatCM = thetaScatCM;
	  */
	
	  double totE_p = ProtonCFT->GetTotalE();
	  double p_p = sqrt(totE_p*totE_p+2.*totE_p*ProtonMass*1000.);
	
	  p_p /= 1000.; // GeV/c
	  double p_cor, e_cor;

	  /*
	    CorrElossOut(&p_cor, &e_cor, p_p, PROTON, 
	    ProtonScatVec.mom/ProtonScatVec.mom.mag(), VertScat);
	  */
	  CorrElossOutWithCFRP(&p_cor, &e_cor, p_p, PROTON, 
			       DirProton/DirProton.mag(), VertScat1,
			       PosProton);

	  EkinCorP1 = (e_cor-ProtonMass)*1000.;
	  

	  for (int i2=0; i2<nPi_CFT; i2++) {
	    CFTParticle *PiCFT = CFTPionCont[i2];
	    ThreeVector DirPi = PiCFT->GetDir();
	    ThreeVector PosPi = PiCFT->GetPos0();

	    ThreeVector ScatProtonMom = DirProton*p_cor/DirProton.mag();
	    ThreeVector ScatHypMom = XiVec1.mom-ScatProtonMom;
	  
	    ThreeVector DecayPiPos = PosPi;
	    ThreeVector DecayPiVec = DirPi;
	  
	    double costDecayPi2=ScatHypMom*DecayPiVec/(ScatHypMom.mag()*DecayPiVec.mag());
	    double thetaDecayPi2=acos(costDecayPi2)*Rad2Deg;
	    event.thetaDecayPi2 = thetaDecayPi2;
	  
	    double cdistDecayPi2;
	    ThreeVector vertexDecayPi2   = VertexPoint3D( VertScat1, DecayPiPos, ScatHypMom, DecayPiVec, cdistDecayPi2 );

	    double decayPiMomCal2;
	    flagDecayCal2 = calcDecayPiMom(ScatHypMom.mag(), 
					   XiMinusMass, 
					   LambdaMass, PionMass,
					   costDecayPi2, &decayPiMomCal2);



	    if (i1==0) {
	      if (i2==0) {
		event.cdistDecayPi_Xi1 = cdistDecayPi2;
		event.vertexDecayPi_Xi1[0] = vertexDecayPi2.x();
		event.vertexDecayPi_Xi1[1] = vertexDecayPi2.y();
		event.vertexDecayPi_Xi1[2] = vertexDecayPi2.z();
		if (flagDecayCal2)
		  event.decayPiMomCal_Xi1 = decayPiMomCal2;	  
	      } else if (i2==1) {
		event.cdistDecayPi_Xi2 = cdistDecayPi2;
		event.vertexDecayPi_Xi2[0] = vertexDecayPi2.x();
		event.vertexDecayPi_Xi2[1] = vertexDecayPi2.y();
		event.vertexDecayPi_Xi2[2] = vertexDecayPi2.z();
		if (flagDecayCal2)
		  event.decayPiMomCal_Xi2 = decayPiMomCal2;	  
		
		if (event.cdistDecayPi_Xi2>0 && 
		    event.cdistDecayPi_Xi2 < event.cdistDecayPi_Xi1 ) {
		  double cdist = event.cdistDecayPi_Xi1;
		  double  vertexDecayPi[3] = {event.vertexDecayPi_Xi1[0],
					      event.vertexDecayPi_Xi1[1],
					      event.vertexDecayPi_Xi1[2]};
		  double decayPiMomCal  = event.decayPiMomCal_Xi1;

		  event.cdistDecayPi_Xi1     = event.cdistDecayPi_Xi2    ;
		  event.vertexDecayPi_Xi1[0] = event.vertexDecayPi_Xi2[0];
		  event.vertexDecayPi_Xi1[1] = event.vertexDecayPi_Xi2[1];
		  event.vertexDecayPi_Xi1[2] = event.vertexDecayPi_Xi2[2];
		  event.decayPiMomCal_Xi1    = event.decayPiMomCal_Xi2   ;

		  event.cdistDecayPi_Xi2     = cdist;
		  event.vertexDecayPi_Xi2[0] = vertexDecayPi[0];
		  event.vertexDecayPi_Xi2[1] = vertexDecayPi[1];
		  event.vertexDecayPi_Xi2[2] = vertexDecayPi[2];
		  event.decayPiMomCal_Xi2    = decayPiMomCal;
		}
	      }
	    } else if (i1==1) {
	      if (i2==0) {
		event.cdistDecayPi_2nd_Xi1 = cdistDecayPi2;
		event.vertexDecayPi_2nd_Xi1[0] = vertexDecayPi2.x();
		event.vertexDecayPi_2nd_Xi1[1] = vertexDecayPi2.y();
		event.vertexDecayPi_2nd_Xi1[2] = vertexDecayPi2.z();
		if (flagDecayCal2)
		  event.decayPiMomCal_2nd_Xi1 = decayPiMomCal2;	  
	      } else if (i2==1) {
		event.cdistDecayPi_2nd_Xi2 = cdistDecayPi2;
		event.vertexDecayPi_2nd_Xi2[0] = vertexDecayPi2.x();
		event.vertexDecayPi_2nd_Xi2[1] = vertexDecayPi2.y();
		event.vertexDecayPi_2nd_Xi2[2] = vertexDecayPi2.z();
		if (flagDecayCal2)
		  event.decayPiMomCal_2nd_Xi2 = decayPiMomCal2;	  

		if (event.cdistDecayPi_2nd_Xi2>0 && 
		    event.cdistDecayPi_2nd_Xi2 < event.cdistDecayPi_2nd_Xi1 ) {
		  double cdist = event.cdistDecayPi_2nd_Xi1;
		  double  vertexDecayPi[3] = {event.vertexDecayPi_2nd_Xi1[0],
					      event.vertexDecayPi_2nd_Xi1[1],
					      event.vertexDecayPi_2nd_Xi1[2]};
		  double decayPiMomCal  = event.decayPiMomCal_2nd_Xi1;

		  event.cdistDecayPi_2nd_Xi1     = event.cdistDecayPi_2nd_Xi2    ;
		  event.vertexDecayPi_2nd_Xi1[0] = event.vertexDecayPi_2nd_Xi2[0];
		  event.vertexDecayPi_2nd_Xi1[1] = event.vertexDecayPi_2nd_Xi2[1];
		  event.vertexDecayPi_2nd_Xi1[2] = event.vertexDecayPi_2nd_Xi2[2];
		  event.decayPiMomCal_2nd_Xi1    = event.decayPiMomCal_2nd_Xi2   ;

		  event.cdistDecayPi_2nd_Xi2     = cdist;
		  event.vertexDecayPi_2nd_Xi2[0] = vertexDecayPi[0];
		  event.vertexDecayPi_2nd_Xi2[1] = vertexDecayPi[1];
		  event.vertexDecayPi_2nd_Xi2[2] = vertexDecayPi[2];
		  event.decayPiMomCal_2nd_Xi2    = decayPiMomCal;
		}
	      }
	    }


	  }

	}
      }


      if (i1==0) {
	event.EkinP = ProtonCFT->GetTotalE();
	event.vertexScat[0] = VertScat1.x();
	event.vertexScat[1] = VertScat1.y();
	event.vertexScat[2] = VertScat1.z();
	event.cdistScat = cdist1;
	event.thetaScat = thetaScat1;
	event.scatMomCal = scatMomCal1;
	event.scatEkinCal = scatEkinCal1;
	event.thetaScatCM = thetaScatCM1;
	event.EkinCorP = EkinCorP1;
	event.hypBeamVec[0] = XiVec1.mom.x();
	event.hypBeamVec[1] = XiVec1.mom.y();
	event.hypBeamVec[2] = XiVec1.mom.z();
      } else if (i1==1) {
	event.EkinP_2nd = ProtonCFT->GetTotalE();
	event.vertexScat_2nd[0] = VertScat1.x();
	event.vertexScat_2nd[1] = VertScat1.y();
	event.vertexScat_2nd[2] = VertScat1.z();
	event.cdistScat_2nd = cdist1;
	event.thetaScat_2nd = thetaScat1;
	event.scatMomCal_2nd = scatMomCal1;
	event.scatEkinCal_2nd = scatEkinCal1;
	event.thetaScatCM_2nd = thetaScatCM1;
	event.EkinCorP_2nd= EkinCorP1;
	event.hypBeamVec_2nd[0] = XiVec1.mom.x();
	event.hypBeamVec_2nd[1] = XiVec1.mom.y();
	event.hypBeamVec_2nd[2] = XiVec1.mom.z();

	double dE1 = -9999.;
	double dE2 = -9999.;

	if (event.EkinCorP>0 && event.scatEkinCal>0)
	  dE1 = event.EkinCorP - event.scatEkinCal;
	if (event.EkinCorP_2nd>0 && event.scatEkinCal_2nd>0)
	  dE2 = event.EkinCorP_2nd - event.scatEkinCal_2nd;

	if (std::abs(dE2) < std::abs(dE1)) {

	  double EkinP = event.EkinP;
	  double vertexScat[3] = {event.vertexScat[0],
				  event.vertexScat[1],
				  event.vertexScat[2]};
	  double cdistScat = event.cdistScat;
	  double thetaScat = event.thetaScat;
	  double scatMomCal = event.scatMomCal;
	  double scatEkinCal = event.scatEkinCal;
	  double thetaScatCM = event.thetaScatCM;
	  double EkinCorP = event.EkinCorP;
	  double hypBeamVec[3] = {event.hypBeamVec[0],
				  event.hypBeamVec[1],
				  event.hypBeamVec[2]};

	  double cdistDecayPi1 = event.cdistDecayPi_Xi1;
	  double vertexDecayPi1[3] = {event.vertexDecayPi_Xi1[0],
				      event.vertexDecayPi_Xi1[1],
				      event.vertexDecayPi_Xi1[2]};
	  double decayPiMomCal1  = event.decayPiMomCal_Xi1;
	  double cdistDecayPi2 = event.cdistDecayPi_Xi2;
	  double vertexDecayPi2[3] = {event.vertexDecayPi_Xi2[0],
				       event.vertexDecayPi_Xi2[1],
				       event.vertexDecayPi_Xi2[2]};
	  double decayPiMomCal2  = event.decayPiMomCal_Xi2;


	  event.EkinP         = event.EkinP_2nd        ;
	  event.vertexScat[0] = event.vertexScat_2nd[0];
	  event.vertexScat[1] = event.vertexScat_2nd[1];
	  event.vertexScat[2] = event.vertexScat_2nd[2];
	  event.cdistScat     = event.cdistScat_2nd    ;
	  event.thetaScat     = event.thetaScat_2nd    ;
	  event.scatMomCal    = event.scatMomCal_2nd   ;
	  event.scatEkinCal   = event.scatEkinCal_2nd  ;
	  event.thetaScatCM   = event.thetaScatCM_2nd  ;
	  event.EkinCorP      = event.EkinCorP_2nd    ;
	  event.hypBeamVec[0] = event.hypBeamVec_2nd[0];
	  event.hypBeamVec[1] = event.hypBeamVec_2nd[1];
	  event.hypBeamVec[2] = event.hypBeamVec_2nd[2];

	  event.cdistDecayPi_Xi1     = event.cdistDecayPi_2nd_Xi1    ;
	  event.vertexDecayPi_Xi1[0] = event.vertexDecayPi_2nd_Xi1[0];	  
	  event.vertexDecayPi_Xi1[1] = event.vertexDecayPi_2nd_Xi1[1];	  
	  event.vertexDecayPi_Xi1[2] = event.vertexDecayPi_2nd_Xi1[2];	  
	  event.decayPiMomCal_Xi1    = event.decayPiMomCal_2nd_Xi1   ;	  
				       				       
	  event.cdistDecayPi_Xi2     = event.cdistDecayPi_2nd_Xi2    ;	  
	  event.vertexDecayPi_Xi2[0] = event.vertexDecayPi_2nd_Xi2[0];	  
	  event.vertexDecayPi_Xi2[1] = event.vertexDecayPi_2nd_Xi2[1];	  
	  event.vertexDecayPi_Xi2[2] = event.vertexDecayPi_2nd_Xi2[2];	  
	  event.decayPiMomCal_Xi2    = event.decayPiMomCal_2nd_Xi2   ;	  

	  event.EkinP_2nd         = EkinP;
	  event.vertexScat_2nd[0] = vertexScat[0];
	  event.vertexScat_2nd[1] = vertexScat[1];
	  event.vertexScat_2nd[2] = vertexScat[2];
	  event.cdistScat_2nd     = cdistScat;
	  event.thetaScat_2nd     = thetaScat;
	  event.scatMomCal_2nd    = scatMomCal;
	  event.scatEkinCal_2nd   = scatEkinCal;
	  event.thetaScatCM_2nd   = thetaScatCM;
	  event.EkinCorP_2nd     = EkinCorP;
	  event.hypBeamVec_2nd[0] = hypBeamVec[0];
	  event.hypBeamVec_2nd[1] = hypBeamVec[1];
	  event.hypBeamVec_2nd[2] = hypBeamVec[2];

	  event.cdistDecayPi_2nd_Xi1     = cdistDecayPi1;
	  event.vertexDecayPi_2nd_Xi1[0] = vertexDecayPi1[0];	  
	  event.vertexDecayPi_2nd_Xi1[1] = vertexDecayPi1[1];	  
	  event.vertexDecayPi_2nd_Xi1[2] = vertexDecayPi1[2];	  
	  event.decayPiMomCal_2nd_Xi1    = decayPiMomCal1;	  
	  				          
	  event.cdistDecayPi_2nd_Xi2     = cdistDecayPi2;	  
	  event.vertexDecayPi_2nd_Xi2[0] = vertexDecayPi2[0];	  
	  event.vertexDecayPi_2nd_Xi2[1] = vertexDecayPi2[1];	  
	  event.vertexDecayPi_2nd_Xi2[2] = vertexDecayPi2[2];	  
	  event.decayPiMomCal_2nd_Xi2    = decayPiMomCal2;	  

	}
      }
    }
  }


  if (indexXiDecayPi != -1 ) {
    if (nPi_CFT ==2 && nP_CFT ==1) {
      for (int i1 = 0; i1 < nPi_CFT; i1++ ) {
        if (i1 == indexXiDecayPi)
          continue;

        CFTParticle *PiCFT1 = CFTPionCont[i1];

        ThreeVector Dir1 = PiCFT1->GetDir();
        ThreeVector Pos1 = PiCFT1->GetPos0();

        for (int i2 = 0; i2 < nP_CFT; i2++ ) {
          CFTParticle *ProtonCFT2 = CFTProtonCont[i2];

          ThreeVector Dir2 = ProtonCFT2->GetDir();
          ThreeVector Pos2 = ProtonCFT2->GetPos0();

          double cdist;
          ThreeVector vertex = VertexPoint3D( Pos1, Pos2, Dir1, Dir2, cdist );

          if (cdist>5)
            continue;

          double cost = Dir1*Dir2/(Dir1.mag()*Dir2.mag());
          double theta = std::acos(cost)*Rad2Deg;

          double totE_pi = PiCFT1->GetTotalE();
          double p_pi = sqrt(totE_pi*totE_pi+2.*totE_pi*PionMass*1000.);
          p_pi /= 1000.; // GeV/c  
          ThreeVector Mom1 = Dir1*(p_pi/Dir1.mag());

          double totE_proton = ProtonCFT2->GetTotalE();
          double p_proton = sqrt(totE_proton*totE_proton+
				 2.*totE_proton*ProtonMass*1000.);
          p_proton /= 1000.; // GeV/c 
	  double p_cor, e_cor;
	  //
	  CorrElossOutWithCFRP(&p_cor, &e_cor, p_proton, PROTON,
			       Dir2*(1./Dir2.mag()), vertex,
			       Pos2);
	  //
	  double EkinCorP1 = (e_cor-ProtonMass)*1000.;
          ThreeVector Mom2 = Dir2*(p_cor/Dir2.mag());

          LorentzVector LvPi( Mom1, 
			      std::sqrt( PionMass*PionMass+Mom1.mag2()) );
          LorentzVector LvP( Mom2, 
			     std::sqrt( ProtonMass*ProtonMass+Mom2.mag2()) );

          LorentzVector LvPiP = LvPi + LvP;
          double InvMass = LvPiP.mag();
          event.invMass_PiP_Lp = InvMass;

          double momPiFromLambda;

          bool flagLambdaInv2 =
            calcPiMomFromLambda(LambdaMass, ProtonMass, p_proton,
				PionMass, cost, &momPiFromLambda);

          ThreeVector momVecLambda(0, 0, 0);
          if (flagLambdaInv2) {
            ThreeVector momVecPi     = Mom1*(momPiFromLambda/Mom1.mag());
            ThreeVector momVecProton = Mom2;
            momVecLambda = momVecPi+momVecProton;

          }
          event.momLambda1_Lp = (Mom1+Mom2).mag();
          event.momLambda2_Lp = momVecLambda.mag();
          event.vertexLambdaDecay_Lp[0] = vertex.x();
          event.vertexLambdaDecay_Lp[1] = vertex.y();
          event.vertexLambdaDecay_Lp[2] = vertex.z();
          event.cdistLambdaDecay_Lp     = cdist;

          /*
	    std::cout << "Vtx CATCH pi & p( " << vertex.x()
	    << ", " << vertex.y() << ", " << vertex.z()
	    << "), cdistBeam = " << cdistBeanAndVertex << std::endl;
	    std::cout << "InvMass = " << InvMass << std::endl; 
	    std::cout << "Mom (Mom1+Mom2) : " << (Mom1+Mom2).mag()
            << ", MissMom = " << LambdaMissMom.mag() << std::endl;
	    std::cout << "Mom (Lambda decay assumption) : " << momVecLambda.mag()
	    << ", MissMom = " << LambdaMissMom.mag() << std::endl;
	  */

	  bool flagLambdaPScat_Lambda = false;
          if (flagLambdaInv2) {

            double cdistLambdaPScat_Lambda;
            ThreeVector VertLambdaPScat_Lambda =
              VertexPoint3D( DecayLambdaVec.pos0, vertex,
			     DecayLambdaVec.mom, momVecLambda, cdistLambdaPScat_Lambda);

            double costLambdaPScat_Lambda=DecayLambdaVec.mom*momVecLambda/
              (DecayLambdaVec.mom.mag()*momVecLambda.mag());
            double thetaLambdaPScat_Lambda=acos(costLambdaPScat_Lambda)*Rad2Deg;

            double momCalLambda1 =-999.;
            double momCalLambda2 =-999.;
            double thetaCMLambdaPScat_Lambda =-999.;

            flagLambdaPScat_Lambda =
              calc2BodyInelastic(LambdaMass, DecayLambdaVec.mom.mag(),
				 ProtonMass, LambdaMass, ProtonMass,
				 thetaLambdaPScat_Lambda,
				 &momCalLambda1, &momCalLambda2,
				 &thetaCMLambdaPScat_Lambda);

            /* 
	       std::cout << "Lambda mom : " << LambdaMissMom.mag()
	       << ", momCalLambda : "  << momCalLambda1
	       << ", momCalLambda2 : "  << momCalLambda2
	       << ", thetaCMLambdaPScat_Lambda "<< thetaCMLambdaPScat_Lambda << std::endl;
	    */

            event.cdistLambdaPScat_Lambda = cdistLambdaPScat_Lambda;
            event.vertexLambdaPScat_Lambda[0] = VertLambdaPScat_Lambda.x();
            event.vertexLambdaPScat_Lambda[1] = VertLambdaPScat_Lambda.y();
            event.vertexLambdaPScat_Lambda[2] = VertLambdaPScat_Lambda.z();

            event.costLambdaPScat_Lambda = costLambdaPScat_Lambda;
            event.thetaLambdaPScat_Lambda = thetaLambdaPScat_Lambda;
            event.thetaCMLambdaPScat_Lambda = thetaCMLambdaPScat_Lambda;
            event.momCalLambda1_Lp = momCalLambda1;
            event.momCalLambda2_Lp = momCalLambda2;

            if (flagLambdaPScat_Lambda) {
              double dp = momCalLambda1 - momVecLambda.mag();
              //std::cout << "Lambda p scat : dp = " << dp << std::endl;
            }
          }

	  // recoil proton assumption
	  {
	    double cdistLambdaPScat_p;
	    ThreeVector VertLambdaPScat_p =
	      VertexPoint3D( DecayLambdaVec.pos0, Pos2,
			     DecayLambdaVec.mom, Dir2, cdistLambdaPScat_p);

	    double totE_proton = ProtonCFT2->GetTotalE();
	    double p_proton = sqrt(totE_proton*totE_proton+
				   2.*totE_proton*ProtonMass*1000.);
	    p_proton /= 1000.; // GeV/c 
	    double p_cor, e_cor;
	    //
	    CorrElossOutWithCFRP(&p_cor, &e_cor, p_proton, PROTON,
				 Dir2*(1./Dir2.mag()), VertLambdaPScat_p,
				 Pos2);
	    //
	    double EkinCorP1 = (e_cor-ProtonMass)*1000.;
	    ThreeVector Mom2 = Dir2*(p_cor/Dir2.mag());
	    event.EkinP_Lp    = totE_proton;
	    event.EkinCorP_Lp = EkinCorP1;

	    double costScat1=DecayLambdaVec.mom*Mom2/(DecayLambdaVec.mom.mag()*Mom2.mag());
	    double thetaScat1=acos(costScat1)*Rad2Deg;

	    //std::cout << "Lambda p scat (1) : theta = " << thetaScat1 << std::endl;
      
	    event.vertexLambdaPScat_p[0] = VertLambdaPScat_p.x();
	    event.vertexLambdaPScat_p[1] = VertLambdaPScat_p.y();
	    event.vertexLambdaPScat_p[2] = VertLambdaPScat_p.z();
	    event.cdistLambdaPScat_p     = cdistLambdaPScat_p;

	    double scatMomCal1 = -999.;
	    double scatEkinCal1 = -999.;
	    double thetaScatCM1 = -999.;

	    if (thetaScat1>0. && thetaScat1<90.) {

	      bool flagLambdaPScat_p =
		calc2BodyKinema(LambdaMass, DecayLambdaVec.mom.mag(),
				ProtonMass, thetaScat1,
				&scatMomCal1, &scatEkinCal1, &thetaScatCM1);


	      if (flagLambdaPScat_p) {
		scatEkinCal1 *= 1000;
		event.scatMomCalLambdaPScat_p  = scatMomCal1;
		event.scatEkinCalLambdaPScat_p = scatEkinCal1;
		event.thetaScatLambdaPScat_p   = thetaScat1;
		event.thetaScatCMLambdaPScat_p = thetaScatCM1;

		double dE_LambdaPScat = EkinCorP1 - scatEkinCal1;
		//std::cout << "Lambda p scat (1) : dE = " << dE_LambdaPScat  << std::endl;
	  
		ThreeVector ScatProtonMom = Mom2;
		ThreeVector ScatHypMom = DecayLambdaVec.mom-ScatProtonMom;

		/*
		std::cout << "Lambda Mom Vec ( " << DecayLambdaVec.mom.x()
			  <<  ", " << DecayLambdaVec.mom.y() 
			  <<  ", " << DecayLambdaVec.mom.z()  << " )"  << std::endl;
		
		std::cout << "Proton Mom Vec ( " << ScatProtonMom.x()
			  <<  ", " << ScatProtonMom.y()
			  <<  ", " << ScatProtonMom.z()  << " )"  << std::endl;
		
		std::cout << "Recoil Sigma Mom Vec ( " << ScatHypMom.x()
			  <<  ", " << ScatHypMom.y()
			  <<  ", " << ScatHypMom.z()  << " )"  << std::endl;
		*/

		double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
		double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;

		event.costScatOpenAngle_Lp  = costScatOpenAngle;
		event.thetaScatOpenAngle_Lp = thetaScatOpenAngle;

		double cdistDecayPi;
		ThreeVector vertexDecayPi
		  = VertexPoint3D( VertLambdaPScat_p, Pos1, 
				   ScatHypMom, Dir1, 
				   cdistDecayPi );

		event.vertexDecayPi2_Lp[0] = vertexDecayPi.x();
		event.vertexDecayPi2_Lp[1] = vertexDecayPi.y();
		event.vertexDecayPi2_Lp[2] = vertexDecayPi.z();
		event.cdistDecayPi2_Lp     = cdistDecayPi;

		/*
		std::cout << "Decay Point of Scattered Lambda (1) ( " << vertexDecayPi.x()
			  << ", " << vertexDecayPi.y()
			  << ", " << vertexDecayPi.z()
			  << " ), cdist = " << cdistDecayPi << std::endl;
		*/
	      }
	    }
	  }
        }
      }
    } else if (nPi_CFT ==2 && nP_CFT ==2) {
      for (int i1 = 0; i1 < nPi_CFT; i1++ ) {
        if (i1 == indexXiDecayPi)
          continue;

        CFTParticle *PiCFT1 = CFTPionCont[i1];

        ThreeVector Dir1 = PiCFT1->GetDir();
        ThreeVector Pos1 = PiCFT1->GetPos0();

        for (int i2 = 0; i2 < nP_CFT; i2++ ) {
          CFTParticle *ProtonCFT2 = CFTProtonCont[i2];

          ThreeVector Dir2 = ProtonCFT2->GetDir();
          ThreeVector Pos2 = ProtonCFT2->GetPos0();

          double cdist;
          ThreeVector vertex = VertexPoint3D( Pos1, Pos2, Dir1, Dir2, cdist );

          if (cdist>5)
            continue;

          double cost = Dir1*Dir2/(Dir1.mag()*Dir2.mag());
          double theta = std::acos(cost)*Rad2Deg;

          double totE_pi = PiCFT1->GetTotalE();
          double p_pi = sqrt(totE_pi*totE_pi+2.*totE_pi*PionMass*1000.);
          p_pi /= 1000.; // GeV/c  
          ThreeVector Mom1 = Dir1*(p_pi/Dir1.mag());

          double totE_proton = ProtonCFT2->GetTotalE();
          double p_proton = sqrt(totE_proton*totE_proton+
				 2.*totE_proton*ProtonMass*1000.);
          p_proton /= 1000.; // GeV/c 

	  double p_cor, e_cor;
	  //
	  CorrElossOutWithCFRP(&p_cor, &e_cor, p_proton, PROTON,
			       Dir2*(1./Dir2.mag()), vertex,
			       Pos2);
	  //
	  ThreeVector Mom2 = Dir2*(p_cor/Dir2.mag());

          LorentzVector LvPi( Mom1, 
			      std::sqrt( PionMass*PionMass+Mom1.mag2()) );
          LorentzVector LvP( Mom2, 
			     std::sqrt( ProtonMass*ProtonMass+Mom2.mag2()) );

          LorentzVector LvPiP = LvPi + LvP;
          double InvMass = LvPiP.mag();
          event.invMass_PiP_Lp = InvMass;

          double momPiFromLambda;

          bool flagLambdaInv2 =
            calcPiMomFromLambda(LambdaMass, ProtonMass, Mom2.mag(),
				PionMass, cost, &momPiFromLambda);

          ThreeVector momVecLambda(0, 0, 0);
          if (flagLambdaInv2) {
            ThreeVector momVecPi     = Mom1*(momPiFromLambda/Mom1.mag());
            ThreeVector momVecProton = Mom2;
            momVecLambda = momVecPi+momVecProton;

          }
          event.momLambda1_Lp = (Mom1+Mom2).mag();
          event.momLambda2_Lp = momVecLambda.mag();
          event.vertexLambdaDecay_Lp[0] = vertex.x();
          event.vertexLambdaDecay_Lp[1] = vertex.y();
          event.vertexLambdaDecay_Lp[2] = vertex.z();
          event.cdistLambdaDecay_Lp     = cdist;

	  if (i2==0) {
	    event.momLambda1_Lp_1st = (Mom1+Mom2).mag();
	    event.momLambda2_Lp_1st = momVecLambda.mag();
	    event.vertexLambdaDecay_Lp_1st[0] = vertex.x();
	    event.vertexLambdaDecay_Lp_1st[1] = vertex.y();
	    event.vertexLambdaDecay_Lp_1st[2] = vertex.z();
	    event.cdistLambdaDecay_Lp_1st     = cdist;
	  } else if (i2==1) {
	    event.momLambda1_Lp_2nd = (Mom1+Mom2).mag();
	    event.momLambda2_Lp_2nd = momVecLambda.mag();
	    event.vertexLambdaDecay_Lp_2nd[0] = vertex.x();
	    event.vertexLambdaDecay_Lp_2nd[1] = vertex.y();
	    event.vertexLambdaDecay_Lp_2nd[2] = vertex.z();
	    event.cdistLambdaDecay_Lp_2nd     = cdist;
	  }

          /*
	    std::cout << "Vtx CATCH pi & p( " << vertex.x()
	    << ", " << vertex.y() << ", " << vertex.z()
	    << "), cdistBeam = " << cdistBeanAndVertex << std::endl;
	    std::cout << "InvMass = " << InvMass << std::endl; 
	    std::cout << "Mom (Mom1+Mom2) : " << (Mom1+Mom2).mag()
            << ", MissMom = " << LambdaMissMom.mag() << std::endl;
	    std::cout << "Mom (Lambda decay assumption) : " << momVecLambda.mag()
	    << ", MissMom = " << LambdaMissMom.mag() << std::endl;
	  */


	  bool flagLambdaPScat_Lambda = false;
          if (flagLambdaInv2) {

            double cdistLambdaPScat_Lambda;
            ThreeVector VertLambdaPScat_Lambda =
              VertexPoint3D( DecayLambdaVec.pos0, vertex,
			     DecayLambdaVec.mom, momVecLambda, cdistLambdaPScat_Lambda);

            double costLambdaPScat_Lambda=DecayLambdaVec.mom*momVecLambda/
              (DecayLambdaVec.mom.mag()*momVecLambda.mag());
            double thetaLambdaPScat_Lambda=acos(costLambdaPScat_Lambda)*Rad2Deg;

            double momCalLambda1 =-999.;
            double momCalLambda2 =-999.;
            double thetaCMLambdaPScat_Lambda =-999.;

            flagLambdaPScat_Lambda =
              calc2BodyInelastic(LambdaMass, DecayLambdaVec.mom.mag(),
				 ProtonMass, LambdaMass, ProtonMass,
				 thetaLambdaPScat_Lambda,
				 &momCalLambda1, &momCalLambda2,
				 &thetaCMLambdaPScat_Lambda);

            /* 
	       std::cout << "Lambda mom : " << DecayLambdaVec.mom.mag()
	       << ", momCalLambda : "  << momCalLambda1
	       << ", momCalLambda2 : "  << momCalLambda2
	       << ", thetaCMLambdaPScat_Lambda "<< thetaCMLambdaPScat_Lambda << std::endl;
	    */

            event.cdistLambdaPScat_Lambda = cdistLambdaPScat_Lambda;
            event.vertexLambdaPScat_Lambda[0] = VertLambdaPScat_Lambda.x();
            event.vertexLambdaPScat_Lambda[1] = VertLambdaPScat_Lambda.y();
            event.vertexLambdaPScat_Lambda[2] = VertLambdaPScat_Lambda.z();

            event.costLambdaPScat_Lambda = costLambdaPScat_Lambda;
            event.thetaLambdaPScat_Lambda = thetaLambdaPScat_Lambda;
            event.thetaCMLambdaPScat_Lambda = thetaCMLambdaPScat_Lambda;
            event.momCalLambda1_Lp = momCalLambda1;
            event.momCalLambda2_Lp = momCalLambda2;

	    if (i2==0) {
	      event.cdistLambdaPScat_Lambda_1st = cdistLambdaPScat_Lambda;
	      event.vertexLambdaPScat_Lambda_1st[0] = VertLambdaPScat_Lambda.x();
	      event.vertexLambdaPScat_Lambda_1st[1] = VertLambdaPScat_Lambda.y();
	      event.vertexLambdaPScat_Lambda_1st[2] = VertLambdaPScat_Lambda.z();
	      
	      event.costLambdaPScat_Lambda_1st = costLambdaPScat_Lambda;
	      event.thetaLambdaPScat_Lambda_1st = thetaLambdaPScat_Lambda;
	      event.thetaCMLambdaPScat_Lambda_1st = thetaCMLambdaPScat_Lambda;
	      event.momCalLambda1_Lp_1st = momCalLambda1;
	      event.momCalLambda2_Lp_1st = momCalLambda2;
	    } else if (i2==1) {
	      event.cdistLambdaPScat_Lambda_2nd = cdistLambdaPScat_Lambda;
	      event.vertexLambdaPScat_Lambda_2nd[0] = VertLambdaPScat_Lambda.x();
	      event.vertexLambdaPScat_Lambda_2nd[1] = VertLambdaPScat_Lambda.y();
	      event.vertexLambdaPScat_Lambda_2nd[2] = VertLambdaPScat_Lambda.z();
	      
	      event.costLambdaPScat_Lambda_2nd = costLambdaPScat_Lambda;
	      event.thetaLambdaPScat_Lambda_2nd = thetaLambdaPScat_Lambda;
	      event.thetaCMLambdaPScat_Lambda_2nd = thetaCMLambdaPScat_Lambda;
	      event.momCalLambda1_Lp_2nd = momCalLambda1;
	      event.momCalLambda2_Lp_2nd = momCalLambda2;
	    }

            if (flagLambdaPScat_Lambda) {
              double dp = momCalLambda1 - momVecLambda.mag();
              //std::cout << "Lambda p scat : dp = " << dp << std::endl;
            }


	    // recoil proton assumption
	    ThreeVector momVecCalRecoilProton = DecayLambdaVec.mom-momVecLambda;
	    
	    CFTParticle *RecoilProtonCFT = CFTProtonCont[0];
	    if (i2==0)
	      RecoilProtonCFT = CFTProtonCont[1];
	    
	    ThreeVector DirRP = RecoilProtonCFT->GetDir();
	    double totE_RP = RecoilProtonCFT->GetTotalE();
	    double p_RP    = sqrt(totE_RP*totE_RP+
				  2.*totE_RP*ProtonMass*1000.);
	    p_RP /= 1000.; // GeV/c 
	    
	    double p_RP_cor, e_RP_cor;
	    //
	    CorrElossOutWithCFRP(&p_RP_cor, &e_RP_cor, p_RP, PROTON,
				 Dir2*(1./Dir2.mag()), VertLambdaPScat_Lambda,
				 Pos2);
	    //
	    ThreeVector MomRP = DirRP*(p_RP_cor/DirRP.mag());
	    
	    double cost_diff = DirRP*momVecCalRecoilProton/
	      (DirRP.mag()*momVecCalRecoilProton.mag());
	    double theta_diff = std::acos(cost_diff)*Rad2Deg;
	    double p_RP_diff = momVecCalRecoilProton.mag() - MomRP.mag();
	    
	    if (i2==0) {
	      event.momCalRecoilProton_Lp_1st = momVecCalRecoilProton.mag();
	      event.momRecoilProton_Lp_1st = MomRP.mag();
	      event.deltaThetaRecoilProton_Lp_1st = theta_diff;
	    } else if (i2==1) {
	      event.momCalRecoilProton_Lp_2nd = momVecCalRecoilProton.mag();
	      event.momRecoilProton_Lp_2nd = MomRP.mag();
	      event.deltaThetaRecoilProton_Lp_2nd = theta_diff;
	    }
          }
	}
      }


      if (event.momCalLambda1_Lp_1st>=0 && event.momCalLambda1_Lp_2nd<0) {
	event.momLambda1_Lp           = event.momLambda1_Lp_1st;          
	event.momLambda2_Lp           = event.momLambda2_Lp_1st;          
	event.vertexLambdaDecay_Lp[0] = event.vertexLambdaDecay_Lp_1st[0];
	event.vertexLambdaDecay_Lp[1] = event.vertexLambdaDecay_Lp_1st[1];
	event.vertexLambdaDecay_Lp[2] = event.vertexLambdaDecay_Lp_1st[2];
	event.cdistLambdaDecay_Lp     = event.cdistLambdaDecay_Lp_1st;    

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_1st;    
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_1st[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_1st[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_1st[2];
					    				 
	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_1st;     
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_1st;    
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_1st;  
	event.momCalLambda1_Lp            = event.momCalLambda1_Lp_1st;              
	event.momCalLambda2_Lp            = event.momCalLambda2_Lp_1st;              

	event.momCalRecoilProton_Lp     = event.momCalRecoilProton_Lp_1st;    
	event.momRecoilProton_Lp        = event.momRecoilProton_Lp_1st;       
	event.deltaThetaRecoilProton_Lp = event.deltaThetaRecoilProton_Lp_1st;
      } else if (event.momCalLambda1_Lp_1st<0 && event.momCalLambda1_Lp_2nd>=0) {
	event.momLambda1_Lp           = event.momLambda1_Lp_2nd;          
	event.momLambda2_Lp           = event.momLambda2_Lp_2nd;          
	event.vertexLambdaDecay_Lp[0] = event.vertexLambdaDecay_Lp_2nd[0];
	event.vertexLambdaDecay_Lp[1] = event.vertexLambdaDecay_Lp_2nd[1];
	event.vertexLambdaDecay_Lp[2] = event.vertexLambdaDecay_Lp_2nd[2];
	event.cdistLambdaDecay_Lp     = event.cdistLambdaDecay_Lp_2nd;    

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_2nd;    
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_2nd[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_2nd[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_2nd[2];
					    				 
	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_2nd;     
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_2nd;    
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_2nd;  
	event.momCalLambda1_Lp            = event.momCalLambda1_Lp_2nd;              
	event.momCalLambda2_Lp            = event.momCalLambda2_Lp_2nd;              

	event.momCalRecoilProton_Lp     = event.momCalRecoilProton_Lp_2nd;    
	event.momRecoilProton_Lp        = event.momRecoilProton_Lp_2nd;       
	event.deltaThetaRecoilProton_Lp = event.deltaThetaRecoilProton_Lp_2nd;
      } else if (event.momCalLambda1_Lp_1st>=0 && event.momCalLambda1_Lp_2nd>=0
		 && std::abs(event.momCalLambda1_Lp_1st-event.momLambda2_Lp_1st) <= 
		 std::abs(event.momCalLambda1_Lp_2nd-event.momLambda2_Lp_2nd)) {
	event.momLambda1_Lp           = event.momLambda1_Lp_1st;          
	event.momLambda2_Lp           = event.momLambda2_Lp_1st;          
	event.vertexLambdaDecay_Lp[0] = event.vertexLambdaDecay_Lp_1st[0];
	event.vertexLambdaDecay_Lp[1] = event.vertexLambdaDecay_Lp_1st[1];
	event.vertexLambdaDecay_Lp[2] = event.vertexLambdaDecay_Lp_1st[2];
	event.cdistLambdaDecay_Lp     = event.cdistLambdaDecay_Lp_1st;    

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_1st;    
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_1st[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_1st[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_1st[2];
					    				 
	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_1st;     
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_1st;    
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_1st;  
	event.momCalLambda1_Lp            = event.momCalLambda1_Lp_1st;              
	event.momCalLambda2_Lp            = event.momCalLambda2_Lp_1st;              

	event.momCalRecoilProton_Lp     = event.momCalRecoilProton_Lp_1st;    
	event.momRecoilProton_Lp        = event.momRecoilProton_Lp_1st;       
	event.deltaThetaRecoilProton_Lp = event.deltaThetaRecoilProton_Lp_1st;
      }  else if (event.momCalLambda1_Lp_1st>=0 && event.momCalLambda1_Lp_2nd>=0
		 && std::abs(event.momCalLambda1_Lp_1st-event.momLambda2_Lp_1st) >= 
		 std::abs(event.momCalLambda1_Lp_2nd-event.momLambda2_Lp_2nd)) {
	event.momLambda1_Lp           = event.momLambda1_Lp_2nd;          
	event.momLambda2_Lp           = event.momLambda2_Lp_2nd;          
	event.vertexLambdaDecay_Lp[0] = event.vertexLambdaDecay_Lp_2nd[0];
	event.vertexLambdaDecay_Lp[1] = event.vertexLambdaDecay_Lp_2nd[1];
	event.vertexLambdaDecay_Lp[2] = event.vertexLambdaDecay_Lp_2nd[2];
	event.cdistLambdaDecay_Lp     = event.cdistLambdaDecay_Lp_2nd;    

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_2nd;    
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_2nd[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_2nd[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_2nd[2];
					    				 
	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_2nd;     
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_2nd;    
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_2nd;  
	event.momCalLambda1_Lp            = event.momCalLambda1_Lp_2nd;              
	event.momCalLambda2_Lp            = event.momCalLambda2_Lp_2nd;              

	event.momCalRecoilProton_Lp     = event.momCalRecoilProton_Lp_2nd;    
	event.momRecoilProton_Lp        = event.momRecoilProton_Lp_2nd;       
	event.deltaThetaRecoilProton_Lp = event.deltaThetaRecoilProton_Lp_2nd;
      }


    } 

    if ( nP_CFT ==2) {
      CFTParticle *ProtonCFT1 = CFTProtonCont[0];
      double EkinP1 = ProtonCFT1->GetTotalE();
      PartVec ProtonVec1;
      ProtonVec1.pos0 = ProtonCFT1->GetPos0();
      ProtonVec1.mom = ProtonCFT1->GetDir();

      /* Lambda p scattering assumption */
      double cdist1;

      ThreeVector vtx = DecayLambdaVec.pos0;

      ThreeVector VertScat1 =
        VertexPoint3D( vtx, ProtonVec1.pos0,
		       DecayLambdaVec.mom, ProtonVec1.mom, cdist1);

      double EkinCorP1 = -999.;
      {
        double totE_p = ProtonCFT1->GetTotalE();
        double p_p = sqrt(totE_p*totE_p+2.*totE_p*ProtonMass*1000.);

        p_p /= 1000.; // GeV/c 
        double p_cor, e_cor;
        //
	CorrElossOutWithCFRP(&p_cor, &e_cor, p_p, PROTON,
			     ProtonVec1.mom*(1./ProtonVec1.mom.mag()), VertScat1,
			     ProtonVec1.pos0);
        //
        EkinCorP1 = (e_cor-ProtonMass)*1000.;
        ProtonVec1.mom *= (p_p/ProtonVec1.mom.mag());
      }
      /*  
	  std::cout << "Lambda p vertex (1): ( " << VertScat1.x() << ", "
	  << VertScat1.y() << ", " << VertScat1.z() << " )"
	  << std::endl;
      */

      CFTParticle *ProtonCFT2 = CFTProtonCont[1];
      double EkinP2 = ProtonCFT2->GetTotalE();
      PartVec ProtonVec2;
      ProtonVec2.pos0 = ProtonCFT2->GetPos0();
      ProtonVec2.mom = ProtonCFT2->GetDir();

      /* Sigma p scattering assumption */
      double cdist2;

      ThreeVector VertScat2 =
        VertexPoint3D( vtx, ProtonVec2.pos0,
		       DecayLambdaVec.mom, ProtonVec2.mom, cdist2);

      /*  
	  std::cout << "Lambda p vertex (2): ( " << VertScat2.x() << ", "
	  << VertScat2.y() << ", " << VertScat2.z() << " )"
	  << std::endl;
      */

      double EkinCorP2 = -999.;
      {
        double totE_p = ProtonCFT2->GetTotalE();
        double p_p = sqrt(totE_p*totE_p+2.*totE_p*ProtonMass*1000.);

        p_p /= 1000.; // GeV/c  
        double p_cor, e_cor;

        //
	CorrElossOutWithCFRP(&p_cor, &e_cor, p_p, PROTON,
			     ProtonVec2.mom*(1./ProtonVec2.mom.mag()), VertScat2,
			     ProtonVec2.pos0);
        //
        EkinCorP2 = (e_cor-ProtonMass)*1000.;
        ProtonVec2.mom *= (p_p/ProtonVec2.mom.mag());
      }

      /* Sigma beam energy correction */
      double costScat1=DecayLambdaVec.mom*ProtonVec1.mom/(DecayLambdaVec.mom.mag()*ProtonVec1.mom.mag());
      double thetaScat1=acos(costScat1)*Rad2Deg;

      //std::cout << "Sigma p scat (1) : theta = " << thetaScat1 << std::endl;
      
      double vertexScat1[3];
      vertexScat1[0] = VertScat1.x();
      vertexScat1[1] = VertScat1.y();
      vertexScat1[2] = VertScat1.z();

      double scatMomCal1 = -999.;
      double scatEkinCal1 = -999.;
      double thetaScatCM1 = -999.;

      double vertexDecayProton1_2p[3] = {-999., -999., -999.};
      double cdistDecayProton1_2p = -999.;

      double MissMassPi_1_2p = -999.;

      if (thetaScat1>0. && thetaScat1<90.) {

        bool flagLambdaPScat =
          calc2BodyKinema(LambdaMass, DecayLambdaVec.mom.mag(),
			  ProtonMass, thetaScat1,
			  &scatMomCal1, &scatEkinCal1, &thetaScatCM1);

        if (flagLambdaPScat) {
          scatEkinCal1 *= 1000;
          /*
            event.scatMomCal = scatMomCal;
            event.scatEkinCal = scatEkinCal*1000;
            event.thetaScatCM = thetaScatCM;
	  */

          double dE_LambdaPScat = EkinCorP1 - scatEkinCal1;
          //std::cout << "Lambda p scat (1) : dE = " << dE_LambdaPScat  << std::endl;
	  
          ThreeVector ScatProtonMom = ProtonVec1.mom;
          ThreeVector ScatHypMom = DecayLambdaVec.mom-ScatProtonMom;

          /*
	    std::cout << "Lambda Mom Vec ( " << DecayLambdaVec.mom.x()
	    <<  ", " << DecayLambdaVec.mom.y() 
	    <<  ", " << DecayLambdaVec.mom.z()  << " )"  << std::endl;
	    
	    std::cout << "Proton Mom Vec ( " << ScatProtonMom.x()
	    <<  ", " << ScatProtonMom.y()
	    <<  ", " << ScatProtonMom.z()  << " )"  << std::endl;
	    
	    std::cout << "Recoil Sigma Mom Vec ( " << ScatHypMom.x()
	    <<  ", " << ScatHypMom.y()
	    <<  ", " << ScatHypMom.z()  << " )"  << std::endl;
	  */

          double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
          double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;

          double cdistDecayProton1;
          ThreeVector vertexDecayProton1
            = VertexPoint3D( VertScat1, ProtonVec2.pos0, 
			     ScatHypMom, ProtonVec2.mom, 
			     cdistDecayProton1 );

          vertexDecayProton1_2p[0] = vertexDecayProton1.x();
          vertexDecayProton1_2p[1] = vertexDecayProton1.y();
          vertexDecayProton1_2p[2] = vertexDecayProton1.z();
          cdistDecayProton1_2p = cdistDecayProton1;

          /*  
	      std::cout << "Decay Point of Scattered Lambda (1) ( " << vertexDecayProton1.x()
	      << ", " << vertexDecayProton1.y()
	      << ", " << vertexDecayProton1.z()
	      << " ), cdist = " << cdistDecayProton1 << std::endl;
	  */
          LorentzVector LvProton(ProtonVec2.mom , 
				 sqrt(ProtonMass*ProtonMass+ProtonVec2.mom.mag2()) );
          LorentzVector LvBeam(ScatHypMom, 
			       sqrt(LambdaMass*LambdaMass+ScatHypMom.mag2()) );
          LorentzVector LvPi = LvBeam-LvProton;
          double MissMassPi0 = LvPi.mag2();

          MissMassPi_1_2p = MissMassPi0;

          //std::cout << "Missing mass of Pi- (1) = "  << MissMassPi0 << std::endl; 
          double costDecayProton=ScatHypMom*ProtonVec2.mom/
            (ScatHypMom.mag()*ProtonVec2.mom.mag());
          double thetaDecayProton=acos(costDecayProton)*Rad2Deg;

          double M;
          double p1, E1;
          double phi = thetaDecayProton;
          //double p4, E4, Ekin4; 

          if (thetaDecayProton>0. && thetaDecayProton<90.) {

            double p4_1, p4_2, E4, Ekin4, E4_2, Ekin4_2;
            double thetaCM;
            bool flagLambdaDecay =
              calc2BodyInelastic( LambdaMass,
				  ScatHypMom.mag(),
				  0.,
				  ProtonMass,
				  PionMass,
				  thetaDecayProton,
				  &p4_1,
				  &p4_2,
				  &thetaCM);

	    if (flagLambdaDecay) {
              E4 = sqrt(ProtonMass*ProtonMass+p4_1*p4_1);
              Ekin4 = E4 - ProtonMass;

	      E4_2 = sqrt(ProtonMass*ProtonMass+p4_2*p4_2);
              Ekin4_2 = E4_2 - ProtonMass;
            }
          }

        }
      }

      /* Analysis of 2nd proton*/

      double costScat2=DecayLambdaVec.mom*ProtonVec2.mom/(DecayLambdaVec.mom.mag()*ProtonVec2.mom.mag());
      double thetaScat2=acos(costScat2)*Rad2Deg;

      //std::cout << "Lambda p scat (2) : theta = " << thetaScat2 << std::endl;

      double vertexScat2[3];
      vertexScat2[0] = VertScat2.x();
      vertexScat2[1] = VertScat2.y();
      vertexScat2[2] = VertScat2.z();

      double scatMomCal2 = -999.;
      double scatEkinCal2 = -999.;
      double thetaScatCM2 = -999.;

      double vertexDecayProton2_2p[3] = {-999., -999., -999.};
      double cdistDecayProton2_2p = -999.;

      double MissMassPi_2_2p = -999.;

      if (thetaScat2>0. && thetaScat2<90.) {

        bool flagLambdaPScat =
          calc2BodyKinema(LambdaMass, DecayLambdaVec.mom.mag(),
			  ProtonMass, thetaScat2,
			  &scatMomCal2, &scatEkinCal2, &thetaScatCM2);

	if (flagLambdaPScat) {
          scatEkinCal2 *= 1000;
          /*
            event.scatMomCal = scatMomCal;
            event.scatEkinCal = scatEkinCal*1000;
            event.thetaScatCM = thetaScatCM;
	  */

          double dE_LambdaPScat = EkinCorP2 - scatEkinCal2;
          //std::cout << "Lambda p scat (2) : dE = " << dE_LambdaPScat  << std::endl;
	  
          ThreeVector ScatProtonMom = ProtonVec2.mom;
          ThreeVector ScatHypMom = DecayLambdaVec.mom-ScatProtonMom;
          /*
	    std::cout << "Proton Mom Vec ( " << ScatProtonMom.x()
	    <<  ", " << ScatProtonMom.y()
	    <<  ", " << ScatProtonMom.z()  << " )"  << std::endl;
	    
	    std::cout << "Recoil Sigma Mom Vec ( " << ScatHypMom.x()
	    <<  ", " << ScatHypMom.y()
	    <<  ", " << ScatHypMom.z()  << " )"  << std::endl;
	  */
          double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
          double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;
          //std::cout << "Opening angle (2) = " << thetaScatOpenAngle << std::endl;
	  
          ThreeVector recoilPoint = VertScat2 + ScatHypMom*(100/ScatHypMom.mag());
	  
          /*
	    std::cout << "Proton1 Mom Vec ( " << ProtonVec1.mom.x() 
	    <<  ", " << ProtonVec1.mom.y()
	    <<  ", " << ProtonVec1.mom.z()  << " )"  << std::endl;
	    
	    std::cout << "Proton1 Mom Pos0 ( " << ProtonVec1.pos0.x() 
	    <<  ", " << ProtonVec1.pos0.y()
	    <<  ", " << ProtonVec1.pos0.z()  << " )"  << std::endl;
	  */

          double cdistDecayProton2;
          ThreeVector vertexDecayProton2
            = VertexPoint3D( VertScat2, ProtonVec1.pos0, ScatHypMom, ProtonVec1.mom, cdistDecayProton2 );

          vertexDecayProton2_2p[0] = vertexDecayProton2.x();
          vertexDecayProton2_2p[1] = vertexDecayProton2.y();
          vertexDecayProton2_2p[2] = vertexDecayProton2.z();
          cdistDecayProton2_2p = cdistDecayProton2;

          //std::cout << "Decay Point of Scattered Lambda (2) ( " << vertexDecayProton2.x()
          //<< ", " << vertexDecayProton2.y()
          //<< ", " << vertexDecayProton2.z()
          //<< " ), cdist = " << cdistDecayProton2 << std::endl; 
	  
          LorentzVector LvProton(ProtonVec1.mom , 
				 sqrt(ProtonMass*ProtonMass+ProtonVec1.mom.mag2()) );
          LorentzVector LvBeam(ScatHypMom, 
			       sqrt(LambdaMass*LambdaMass+ScatHypMom.mag2()) );
          LorentzVector LvPi = LvBeam-LvProton;
          double MissMassPi0 = LvPi.mag2();

          MissMassPi_2_2p = MissMassPi0;

          //std::cout << "Missing mass of Pi- (2) = "  << MissMassPi0 << std::endl;

          double costDecayProton=ScatHypMom*ProtonVec1.mom/
            (ScatHypMom.mag()*ProtonVec1.mom.mag());
          double thetaDecayProton=acos(costDecayProton)*Rad2Deg;

          double M;
          double p1, E1;
          double phi = thetaDecayProton;
          //double p4, E4, Ekin4;                               

	  if (thetaDecayProton>0. && thetaDecayProton<90.) {

            double p4_1, p4_2, E4, Ekin4, E4_2, Ekin4_2;
            double thetaCM;
            bool flagLambdaDecay =
              calc2BodyInelastic( LambdaMass,
				  ScatHypMom.mag(),
				  0.,
				  ProtonMass,
				  PionMass,
				  thetaDecayProton,
				  &p4_1,
				  &p4_2,
				  &thetaCM);

            if (flagLambdaDecay) {
              E4 = sqrt(ProtonMass*ProtonMass+p4_1*p4_1);
              Ekin4 = E4 - ProtonMass;

              E4_2 = sqrt(ProtonMass*ProtonMass+p4_2*p4_2);
              Ekin4_2 = E4_2 - ProtonMass;
            }
          }
        }
      }

      double dE_LambdaPScat1 = -999.;
      if (EkinCorP1>0 && scatEkinCal1>0)
        dE_LambdaPScat1 = EkinCorP1 - scatEkinCal1;

      double dE_LambdaPScat2 = -999.;
      if (EkinCorP2>0 && scatEkinCal2>0)
        dE_LambdaPScat2 = EkinCorP2 - scatEkinCal2;


      if (fabs(dE_LambdaPScat1) <= fabs(dE_LambdaPScat2)) {
        //std::cout << "dE (1st) = " << dE_LambdaPScat1
        //<< ", dE (2nd) = " << dE_LambdaPScat2 
        //<< std::endl;

        event.EkinP_Lp = ProtonCFT1->GetTotalE();
        event.vertexLambdaPScat_p[0] = VertScat1.x();
        event.vertexLambdaPScat_p[1] = VertScat1.y();
        event.vertexLambdaPScat_p[2] = VertScat1.z();
        event.cdistLambdaPScat_p = cdist1;
        event.thetaLambdaPScat_Lambda = thetaScat1;
        event.scatMomCalLambdaPScat_p = scatMomCal1;
        event.scatEkinCalLambdaPScat_p = scatEkinCal1;
        event.thetaScatCMLambdaPScat_p = thetaScatCM1;
        event.EkinCorP_Lp = EkinCorP1;
        event.vertexDecayProton1_Lp_2p[0] = vertexDecayProton1_2p[0];
        event.vertexDecayProton1_Lp_2p[1] = vertexDecayProton1_2p[1];
        event.vertexDecayProton1_Lp_2p[2] = vertexDecayProton1_2p[2];
        event.cdistDecayProton1_Lp_2p     = cdistDecayProton1_2p;
        event.MissMassPi_1_Lp_2p         = MissMassPi_1_2p;

      } else {
        //std::cout << "dE (2nd) = " << dE_LambdaPScat2
        //<< ", dE (1st) = " << dE_LambdaPScat1
        //<< std::endl;
        event.EkinP_Lp = ProtonCFT2->GetTotalE();
        event.vertexLambdaPScat_p[0] = VertScat2.x();
        event.vertexLambdaPScat_p[1] = VertScat2.y();
        event.vertexLambdaPScat_p[2] = VertScat2.z();
        event.cdistLambdaPScat_p = cdist2;
        event.thetaLambdaPScat_Lambda  = thetaScat2;
        event.scatMomCalLambdaPScat_p  = scatMomCal2;
        event.scatEkinCalLambdaPScat_p = scatEkinCal2;
        event.thetaScatCMLambdaPScat_p = thetaScatCM2;
        event.EkinCorP_Lp = EkinCorP2;

        event.vertexDecayProton1_Lp_2p[0] = vertexDecayProton2_2p[0];
        event.vertexDecayProton1_Lp_2p[1] = vertexDecayProton2_2p[1];
        event.vertexDecayProton1_Lp_2p[2] = vertexDecayProton2_2p[2];
        event.cdistDecayProton1_Lp_2p     = cdistDecayProton2_2p;
        event.MissMassPi_1_Lp_2p         = MissMassPi_2_2p;

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
  double E1 = sqrt(KaonMass*KaonMass + bmom*bmom);
  double A = KaonMass*KaonMass + ProtonMass*ProtonMass
    + KaonMass*KaonMass - XiMinusMass*XiMinusMass
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
  event.nXi=0;
  for (int i=0; i<MaxHits; i++) {
    event.p[i]=-999.0;
    event.pK[i]=-999.0;
    event.pKCal[i]=-999.0;
    event.theta[i]=-999.0;
    event.MissMass[i]=-999.0;
    event.pXiCal[i]=-999.0;
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
  event.EkinP_2nd = -999.0;
  event.EkinCorP_2nd = -999.0;

  event.vertexDecayPi[0] = -999.9;
  event.vertexDecayPi[1] = -999.9;
  event.vertexDecayPi[2] = -999.9;
  event.cdistDecayPi = -999.9;
  event.thetaDecayPi = -999.9;
  event.decayPiMomCal = -999.9;
  event.decayPiMomVecCal[0] = -999.9;
  event.decayPiMomVecCal[1] = -999.9;
  event.decayPiMomVecCal[2] = -999.9;
  event.decayLambdaMomCal = -999.9;
  event.decayLambdaMomVecCal[0] = -999.9;
  event.decayLambdaMomVecCal[1] = -999.9;
  event.decayLambdaMomVecCal[2] = -999.9;

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

  event.thetaDecayPi_Xi1 = -999.9;
  event.decayPiMomCal_Xi1 = -999.9;
  event.vertexDecayPi_Xi1[0] = -999.9;
  event.vertexDecayPi_Xi1[1] = -999.9;
  event.vertexDecayPi_Xi1[2] = -999.9;
  event.cdistDecayPi_Xi1  = -999.;

  event.thetaDecayPi_Xi2 = -999.9;
  event.decayPiMomCal_Xi2 = -999.9;
  event.vertexDecayPi_Xi2[0] = -999.9;
  event.vertexDecayPi_Xi2[1] = -999.9;
  event.vertexDecayPi_Xi2[2] = -999.9;
  event.cdistDecayPi_Xi2  = -999.;

  event.vertexScat_2nd[0] = -999.9;
  event.vertexScat_2nd[1] = -999.9;
  event.vertexScat_2nd[2] = -999.9;

  event.hypBeamVec_2nd[0] = -999.9;
  event.hypBeamVec_2nd[1] = -999.9;
  event.hypBeamVec_2nd[2] = -999.9;
  
  event.cdistScat_2nd = -999.9;
  event.thetaScat_2nd = -999.9;
  event.scatMomCal_2nd = -999.9;
  event.scatEkinCal_2nd = -999.9;
  event.thetaScatCM_2nd = -999.9;

  event.thetaDecayPi_2nd_Xi1 = -999.9;
  event.decayPiMomCal_2nd_Xi1 = -999.9;
  event.vertexDecayPi_2nd_Xi1[0] = -999.9;
  event.vertexDecayPi_2nd_Xi1[1] = -999.9;
  event.vertexDecayPi_2nd_Xi1[2] = -999.9;
  event.cdistDecayPi_2nd_Xi1  = -999.;

  event.thetaDecayPi_2nd_Xi2 = -999.9;
  event.decayPiMomCal_2nd_Xi2 = -999.9;
  event.vertexDecayPi_2nd_Xi2[0] = -999.9;
  event.vertexDecayPi_2nd_Xi2[1] = -999.9;
  event.vertexDecayPi_2nd_Xi2[2] = -999.9;
  event.cdistDecayPi_2nd_Xi2  = -999.;


  event.thetaDecayPi2 = -999.9;
  event.decayPiMomCal2 = -999.9;
  event.vertexDecayPi2[0] = -999.9;
  event.vertexDecayPi2[1] = -999.9;
  event.vertexDecayPi2[2] = -999.9;
  event.cdistDecayPi2  = -999.;



  event.vertexNpScat[0] = -999.9;
  event.vertexNpScat[1] = -999.9;
  event.vertexNpScat[2] = -999.9;

  event.cdistNpScat = -999.9;
  event.thetaNpScat = -999.9;
  event.scatNpMomCal = -999.9;
  event.scatNpEkinCal = -999.9;
  event.thetaNpScatCM = -999.9;

  event.nLambda = 0;
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

  event.cdistLambdaLambdaConv = -999.9;
  event.thetaLambdaLambdaConv = -999.9;
  event.vertexLambdaLambdaConv[0] = -999.9;
  event.vertexLambdaLambdaConv[1] = -999.9;
  event.vertexLambdaLambdaConv[2] = -999.9;

  event.momCalLambda = -999.9;
  event.momCalLambda2 = -999.9;
  event.thetaCMLambdaLambdaConv = -999.9;
  event.invMass_PiP = -999.9;

  event.cdistLambdaDecay_2nd = -999.9;
  event.thetaLambdaDecay_2nd = -999.9;
  event.vertexLambdaDecay_2nd[0] = -999.9;
  event.vertexLambdaDecay_2nd[1] = -999.9;
  event.vertexLambdaDecay_2nd[2] = -999.9;

  event.momPiFromLambda_2nd = -999.9;
  event.momProtonFromLambda_2nd = -999.9;
  event.momLambda_2nd = -999.9;
  event.momVecPiFromLambda_2nd[0] = -999.9;
  event.momVecPiFromLambda_2nd[1] = -999.9;
  event.momVecPiFromLambda_2nd[2] = -999.9;
  event.momVecProtonFromLambda_2nd[0] = -999.9;
  event.momVecProtonFromLambda_2nd[1] = -999.9;
  event.momVecProtonFromLambda_2nd[2] = -999.9;
  event.momVecLambda_2nd[0] = -999.9;
  event.momVecLambda_2nd[1] = -999.9;
  event.momVecLambda_2nd[2] = -999.9;

  event.cdistLambdaLambdaConv_2nd = -999.9;
  event.thetaLambdaLambdaConv_2nd = -999.9;
  event.vertexLambdaLambdaConv_2nd[0] = -999.9;
  event.vertexLambdaLambdaConv_2nd[1] = -999.9;
  event.vertexLambdaLambdaConv_2nd[2] = -999.9;

  event.momCalLambda_2nd = -999.9;
  event.momCalLambda2_2nd = -999.9;
  event.thetaCMLambdaLambdaConv_2nd = -999.9;
  event.invMass_PiP_2nd = -999.9;


  event.momPiFromXi       = -999.9;
  event.momXi             = -999.9;
  event.momVecPiFromXi[0] = -999.9;
  event.momVecPiFromXi[1] = -999.9;
  event.momVecPiFromXi[2] = -999.9;
  event.momVecXi[0]       = -999.9;
  event.momVecXi[1]       = -999.9;
  event.momVecXi[2]       = -999.9;
  event.cdistLambdaPi     = -999.9;
  event.thetaLambdaPi     = -999.9;
  event.vertexLambdaPi[0] = -999.9;
  event.vertexLambdaPi[1] = -999.9;
  event.vertexLambdaPi[2] = -999.9;

  event.cdistXiPScat_Xi = -999.;
  event.vertexXiPScat_Xi[0] = -999.;
  event.vertexXiPScat_Xi[1] = -999.;
  event.vertexXiPScat_Xi[2] = -999.;
  event.costXiPScat_Xi = -999.;
  event.thetaXiPScat_Xi = -999.;
  event.thetaCMXiPScat_Xi = -999.;
  event.momCalXi1 = -999.;
  event.momCalXi2 = -999.;
  //event.momCalRecoilProton = -999.;
  //event.momRecoilProton = -999.;
  //event.deltaThetaRecoilProton = -999.;

  event.EkinP_Lp = -999.;
  event.EkinCorP_Lp = -999.;
  event.EkinP_2nd_Lp = -999.;
  event.EkinCorP_2nd_Lp = -999.;
  event.invMass_PiP_Lp = -999.;
  event.momLambda1_Lp = -999.;
  event.momLambda2_Lp = -999.;
  event.vertexLambdaDecay_Lp[0] = -999.;
  event.vertexLambdaDecay_Lp[1] = -999.;
  event.vertexLambdaDecay_Lp[2] = -999.;
  event.cdistLambdaDecay_Lp = -999.;
  event.cdistLambdaPScat_Lambda = -999.;
  event.vertexLambdaPScat_Lambda[0] = -999.;
  event.vertexLambdaPScat_Lambda[1] = -999.;
  event.vertexLambdaPScat_Lambda[2] = -999.;
  event.costLambdaPScat_Lambda = -999.;
  event.thetaLambdaPScat_Lambda = -999.;
  event.thetaCMLambdaPScat_Lambda = -999.;
  event.momCalLambda1_Lp = -999.;
  event.momCalLambda2_Lp = -999.;
  event.momCalRecoilProton_Lp = -999.;
  event.momRecoilProton_Lp = -999.;
  event.deltaThetaRecoilProton_Lp = -999.;
  event.vertexLambdaPScat_p[0] = -999.;
  event.vertexLambdaPScat_p[1] = -999.;
  event.vertexLambdaPScat_p[2] = -999.;
  event.cdistLambdaPScat_p = -999.;
  event.thetaScatLambdaPScat_p = -999.;
  event.scatMomCalLambdaPScat_p = -999.;
  event.scatEkinCalLambdaPScat_p = -999.;
  event.thetaScatCMLambdaPScat_p = -999.;
  event.costScatOpenAngle_Lp = -999.;
  event.thetaScatOpenAngle_Lp = -999.;
  event.cdistDecayPi2_Lp = -999.;
  event.vertexDecayPi2_Lp[0] = -999.;
  event.vertexDecayPi2_Lp[1] = -999.;
  event.vertexDecayPi2_Lp[2] = -999.;
  event.cdistDecayProton1_Lp_2p = -999.;
  event.vertexDecayProton1_Lp_2p[0] = -999.;
  event.vertexDecayProton1_Lp_2p[1] = -999.;
  event.vertexDecayProton1_Lp_2p[2] = -999.;
  event.MissMassPi_1_Lp_2p = -999.;

  event.momLambda1_Lp_1st = -999.;
  event.momLambda2_Lp_1st = -999.;
  event.vertexLambdaDecay_Lp_1st[0] = -999.;
  event.vertexLambdaDecay_Lp_1st[1] = -999.;
  event.vertexLambdaDecay_Lp_1st[2] = -999.;
  event.cdistLambdaDecay_Lp_1st = -999.;
  event.cdistLambdaPScat_Lambda_1st = -999.;
  event.vertexLambdaPScat_Lambda_1st[0] = -999.;
  event.vertexLambdaPScat_Lambda_1st[1] = -999.;
  event.vertexLambdaPScat_Lambda_1st[2] = -999.;
  event.costLambdaPScat_Lambda_1st = -999.;
  event.thetaLambdaPScat_Lambda_1st = -999.;
  event.thetaCMLambdaPScat_Lambda_1st = -999.;
  event.momCalLambda1_Lp_1st = -999.;
  event.momCalLambda2_Lp_1st = -999.;
  event.momCalRecoilProton_Lp_1st = -999.;
  event.momRecoilProton_Lp_1st = -999.;
  event.deltaThetaRecoilProton_Lp_1st = -999.;

  event.momLambda1_Lp_2nd = -999.;
  event.momLambda2_Lp_2nd = -999.;
  event.vertexLambdaDecay_Lp_2nd[0] = -999.;
  event.vertexLambdaDecay_Lp_2nd[1] = -999.;
  event.vertexLambdaDecay_Lp_2nd[2] = -999.;
  event.cdistLambdaDecay_Lp_2nd = -999.;
  event.cdistLambdaPScat_Lambda_2nd = -999.;
  event.vertexLambdaPScat_Lambda_2nd[0] = -999.;
  event.vertexLambdaPScat_Lambda_2nd[1] = -999.;
  event.vertexLambdaPScat_Lambda_2nd[2] = -999.;
  event.costLambdaPScat_Lambda_2nd = -999.;
  event.thetaLambdaPScat_Lambda_2nd = -999.;
  event.thetaCMLambdaPScat_Lambda_2nd = -999.;
  event.momCalLambda1_Lp_2nd = -999.;
  event.momCalLambda2_Lp_2nd = -999.;
  event.momCalRecoilProton_Lp_2nd = -999.;
  event.momRecoilProton_Lp_2nd = -999.;
  event.deltaThetaRecoilProton_Lp_2nd = -999.;


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

  tree->Branch("nXi",   &event.nXi,  "nXi/I");
  tree->Branch("pKCal",   event.pKCal,  "pKCal[nXi]/D");
  tree->Branch("theta",   event.theta,  "theta[nXi]/D");
  tree->Branch("pXiCal",   event.pXiCal,  "pXiCal[nXi]/D");
  tree->Branch("Vertex_x",   event.Vertex_x,  "Vertex_x[nXi]/D");
  tree->Branch("Vertex_y",   event.Vertex_y,  "Vertex_y[nXi]/D");
  tree->Branch("Vertex_z",   event.Vertex_z,  "Vertex_z[nXi]/D");

  tree->Branch("nPi_CFT",   &event.nPi_CFT,  "nPi_CFT/I");
  tree->Branch("nK_CFT",   &event.nK_CFT,  "nK_CFT/I");
  tree->Branch("nP_CFT",   &event.nP_CFT,  "nP_CFT/I");
  tree->Branch("EkinP",   &event.EkinP,  "EkinP/D");
  tree->Branch("EkinCorP",   &event.EkinCorP,  "EkinCorP/D");
  tree->Branch("EkinP_2nd",   &event.EkinP_2nd,  "EkinP_2nd/D");
  tree->Branch("EkinCorP_2nd",   &event.EkinCorP_2nd,  "EkinCorP_2nd/D");

  tree->Branch("vertexDecayPi",   event.vertexDecayPi,  "vertexDecayPi[3]/D");
  tree->Branch("cdistDecayPi",   &event.cdistDecayPi,  "cdistDecayPi/D");
  tree->Branch("thetaDecayPi",   &event.thetaDecayPi,  "thetaDecayPi/D");
  tree->Branch("decayPiMomCal",   &event.decayPiMomCal,  "decayPiMomCal/D");
  tree->Branch("decayPiMomVecCal",   event.decayPiMomVecCal,  "decayPiMomVecCal/D");
  tree->Branch("decayLambdaMomCal",   &event.decayLambdaMomCal,  "decayLambdaMomCal/D");
  tree->Branch("decayLambdaMomVecCal",   event.decayLambdaMomVecCal,  "decayLambdaMomVecCal[3]/D");

  tree->Branch("vertexScat",   event.vertexScat,  "vertexScat[3]/D");
  tree->Branch("hypBeamVec",   event.hypBeamVec,  "hypBeamVec[3]/D");
  tree->Branch("cdistScat",   &event.cdistScat,  "cdistScat/D");
  tree->Branch("thetaScat",   &event.thetaScat,  "thetaScat/D");
  tree->Branch("scatMomCal",   &event.scatMomCal,  "scatMomCal/D");
  tree->Branch("scatEkinCal",   &event.scatEkinCal,  "scatEkinCal/D");
  tree->Branch("thetaScatCM",   &event.thetaScatCM,  "thetaScatCM/D");

  tree->Branch("thetaDecayPi_Xi1",   &event.thetaDecayPi_Xi1,  "thetaDecayPi_Xi1/D");
  tree->Branch("decayPiMomCal_Xi1",   &event.decayPiMomCal_Xi1,  "decayPiMomCal_Xi1/D");
  tree->Branch("vertexDecayPi_Xi1",   event.vertexDecayPi_Xi1,  "vertexDecayPi_Xi1[3]/D");
  tree->Branch("cdistDecayPi_Xi1",   &event.cdistDecayPi_Xi1,  "cdistDecayPi_Xi1/D");
  tree->Branch("thetaDecayPi_Xi2",   &event.thetaDecayPi_Xi2,  "thetaDecayPi_Xi2/D");
  tree->Branch("decayPiMomCal_Xi2",   &event.decayPiMomCal_Xi2,  "decayPiMomCal_Xi2/D");
  tree->Branch("vertexDecayPi_Xi2",   event.vertexDecayPi_Xi2,  "vertexDecayPi_Xi2[3]/D");
  tree->Branch("cdistDecayPi_Xi2",   &event.cdistDecayPi_Xi2,  "cdistDecayPi_Xi2/D");


  tree->Branch("vertexScat_2nd",   event.vertexScat_2nd,  "vertexScat_2nd[3]/D");
  tree->Branch("hypBeamVec_2nd",   event.hypBeamVec_2nd,  "hypBeamVec_2nd[3]/D");
  tree->Branch("cdistScat_2nd",   &event.cdistScat_2nd,  "cdistScat_2nd/D");
  tree->Branch("thetaScat_2nd",   &event.thetaScat_2nd,  "thetaScat_2nd/D");
  tree->Branch("scatMomCal_2nd",   &event.scatMomCal_2nd,  "scatMomCal_2nd/D");
  tree->Branch("scatEkinCal_2nd",   &event.scatEkinCal_2nd,  "scatEkinCal_2nd/D");
  tree->Branch("thetaScatCM_2nd",   &event.thetaScatCM_2nd,  "thetaScatCM_2nd/D");

  tree->Branch("thetaDecayPi_2nd_Xi1",   &event.thetaDecayPi_2nd_Xi1,  "thetaDecayPi_2nd_Xi1/D");
  tree->Branch("decayPiMomCal_2nd_Xi1",   &event.decayPiMomCal_2nd_Xi1,  "decayPiMomCal_2nd_Xi1/D");
  tree->Branch("vertexDecayPi_2nd_Xi1",   event.vertexDecayPi_2nd_Xi1,  "vertexDecayPi_2nd_Xi1[3]/D");
  tree->Branch("cdistDecayPi_2nd_Xi1",   &event.cdistDecayPi_2nd_Xi1,  "cdistDecayPi_2nd_Xi1/D");
  tree->Branch("thetaDecayPi_2nd_Xi2",   &event.thetaDecayPi_2nd_Xi2,  "thetaDecayPi_2nd_Xi2/D");
  tree->Branch("decayPiMomCal_2nd_Xi2",   &event.decayPiMomCal_2nd_Xi2,  "decayPiMomCal_2nd_Xi2/D");
  tree->Branch("vertexDecayPi_2nd_Xi2",   event.vertexDecayPi_2nd_Xi2,  "vertexDecayPi_2nd_Xi2[3]/D");
  tree->Branch("cdistDecayPi_2nd_Xi2",   &event.cdistDecayPi_2nd_Xi2,  "cdistDecayPi_2nd_Xi2/D");


  tree->Branch("thetaDecayPi2",   &event.thetaDecayPi2,  "thetaDecayPi2/D");
  tree->Branch("decayPiMomCal2",   &event.decayPiMomCal2,  "decayPiMomCal2/D");
  tree->Branch("vertexDecayPi2",   event.vertexDecayPi2,  "vertexDecayPi2[3]/D");
  tree->Branch("cdistDecayPi2",   &event.cdistDecayPi2,  "cdistDecayPi2/D");


  tree->Branch("vertexNpScat",   event.vertexNpScat,  "vertexNpScat[3]/D");
  tree->Branch("cdistNpScat",   &event.cdistNpScat,  "cdistNpScat/D");
  tree->Branch("thetaNpScat",   &event.thetaNpScat,  "thetaNpScat/D");
  tree->Branch("scatNpMomCal",   &event.scatNpMomCal,  "scatNpMomCal/D");
  tree->Branch("scatNpEkinCal",   &event.scatNpEkinCal,  "scatNpEkinCal/D");
  tree->Branch("thetaNpScatCM",   &event.thetaNpScatCM,  "thetaNpScatCM/D");


  tree->Branch("nLambda",   &event.nLambda,  "nLambda/I");
  tree->Branch("cdistLambdaDecay",   &event.cdistLambdaDecay,  "cdistLambdaDecay/D");
  tree->Branch("thetaLambdaDecay",   &event.thetaLambdaDecay,  "thetaLambdaDecay/D");
  tree->Branch("vertexLambdaDecay",   event.vertexLambdaDecay,  "vertexLambdaDecay[3]/D");
  tree->Branch("momPiFromLambda",   &event.momPiFromLambda,  "momPiFromLambda/D");
  tree->Branch("momProtonFromLambda",   &event.momProtonFromLambda,  "momProtonFromLambda/D");
  tree->Branch("momLambda",   &event.momLambda,  "momLambda/D");
  tree->Branch("momVecPiFromLambda",   event.momVecPiFromLambda,  "momVecPiFromLambda/D");
  tree->Branch("momVecProtonFromLambda",   event.momVecProtonFromLambda,  "momVecProtonFromLambda/D");
  tree->Branch("momVecLambda",   event.momVecLambda,  "momVecLambda/D");
  tree->Branch("cdistLambdaLambdaConv",   &event.cdistLambdaLambdaConv,  "cdistLambdaLambdaConv/D");
  tree->Branch("thetaLambdaLambdaConv",   &event.thetaLambdaLambdaConv,  "thetaLambdaLambdaConv/D");
  tree->Branch("vertexLambdaLambdaConv",   event.vertexLambdaLambdaConv,  "vertexLambdaLambdaConv[3]/D");
  tree->Branch("momCalLambda",   &event.momCalLambda,  "momCalLambda/D");
  tree->Branch("momCalLambda2",   &event.momCalLambda2,  "momCalLambda2/D");
  tree->Branch("thetaCMLambdaLambdaConv",   &event.thetaCMLambdaLambdaConv,  "thetaCMLambdaLambdaConv/D");
  tree->Branch("invMass_PiP",   &event.invMass_PiP,  "invMass_PiP/D");

  tree->Branch("cdistLambdaDecay_2nd",   &event.cdistLambdaDecay_2nd,  "cdistLambdaDecay_2nd/D");
  tree->Branch("thetaLambdaDecay_2nd",   &event.thetaLambdaDecay_2nd,  "thetaLambdaDecay_2nd/D");
  tree->Branch("vertexLambdaDecay_2nd",   event.vertexLambdaDecay_2nd,  "vertexLambdaDecay_2nd[3]/D");
  tree->Branch("momPiFromLambda_2nd",   &event.momPiFromLambda_2nd,  "momPiFromLambda_2nd/D");
  tree->Branch("momProtonFromLambda_2nd",   &event.momProtonFromLambda_2nd,  "momProtonFromLambda_2nd/D");
  tree->Branch("momLambda_2nd",   &event.momLambda_2nd,  "momLambda_2nd/D");
  tree->Branch("momVecPiFromLambda_2nd",   event.momVecPiFromLambda_2nd,  "momVecPiFromLambda_2nd/D");
  tree->Branch("momVecProtonFromLambda_2nd",   event.momVecProtonFromLambda_2nd,  "momVecProtonFromLambda_2nd/D");
  tree->Branch("momVecLambda_2nd",   event.momVecLambda_2nd,  "momVecLambda_2nd/D");
  tree->Branch("cdistLambdaLambdaConv_2nd",   &event.cdistLambdaLambdaConv_2nd,  "cdistLambdaLambdaConv_2nd/D");
  tree->Branch("thetaLambdaLambdaConv_2nd",   &event.thetaLambdaLambdaConv_2nd,  "thetaLambdaLambdaConv_2nd/D");
  tree->Branch("vertexLambdaLambdaConv_2nd",   event.vertexLambdaLambdaConv_2nd,  "vertexLambdaLambdaConv_2nd[3]/D");
  tree->Branch("momCalLambda_2nd",   &event.momCalLambda_2nd,  "momCalLambda_2nd/D");
  tree->Branch("momCalLambda2_2nd",   &event.momCalLambda2_2nd,  "momCalLambda2_2nd/D");
  tree->Branch("thetaCMLambdaLambdaConv_2nd",   &event.thetaCMLambdaLambdaConv_2nd,  "thetaCMLambdaLambdaConv_2nd/D");
  tree->Branch("invMass_PiP_2nd",   &event.invMass_PiP_2nd,  "invMass_PiP_2nd/D");

  tree->Branch("momPiFromXi",   &event.momPiFromXi,  "momPiFromXi/D");
  tree->Branch("momXi",   &event.momXi,  "momXi/D");
  tree->Branch("momVecPiFromXi",   event.momVecPiFromXi,  "momVecPiFromXi[3]/D");
  tree->Branch("momVecXi",   event.momVecXi,  "momVecXi[3]/D");
  tree->Branch("cdistLambdaPi",   &event.cdistLambdaPi,  "cdistLambdaPi/D");
  tree->Branch("thetaLambdaPi",   &event.thetaLambdaPi,  "thetaLambdaPi/D");
  tree->Branch("vertexLambdaPi",   event.vertexLambdaPi,  "vertexLambdaPi[3]/D");

  tree->Branch("cdistXiPScat_Xi",   &event.cdistXiPScat_Xi, "cdistXiPScat_Xi/D");
  tree->Branch("vertexXiPScat_Xi",   event.vertexXiPScat_Xi, "vertexXiPScat_Xi[3]/D");
  tree->Branch("costXiPScat_Xi",   &event.costXiPScat_Xi, "costXiPScat_Xi/D");
  tree->Branch("thetaXiPScat_Xi",   &event.thetaXiPScat_Xi, "thetaXiPScat_Xi/D");
  tree->Branch("thetaCMXiPScat_Xi",   &event.thetaCMXiPScat_Xi, "thetaCMXiPScat_Xi/D");
  tree->Branch("momCalXi1",   &event.momCalXi1, "momCalXi1/D");
  tree->Branch("momCalXi2",   &event.momCalXi2, "momCalXi2/D");
  //tree->Branch("momCalRecoilProton",   &event.momCalRecoilProton, "momCalRecoilProton/D");
  //tree->Branch("momRecoilProton",   &event.momRecoilProton, "momRecoilProton/D");
  //tree->Branch("deltaThetaRecoilProton",   &event.deltaThetaRecoilProton, "deltaThetaRecoilProton/D");

  tree->Branch("EkinP_Lp",                   &event.EkinP_Lp,			"EkinP_Lp/D");			
  tree->Branch("EkinCorP_Lp", 		     &event.EkinCorP_Lp,		"EkinCorP_Lp/D");				
  tree->Branch("EkinP_2nd_Lp", 		     &event.EkinP_2nd_Lp,		"EkinP_2nd_Lp/D");				
  tree->Branch("EkinCorP_2nd_Lp", 	     &event.EkinCorP_2nd_Lp,		"EkinCorP_2nd_Lp/D");		
  tree->Branch("invMass_PiP_Lp", 	     &event.invMass_PiP_Lp,		"invMass_PiP_Lp/D");		
  tree->Branch("momLambda1_Lp", 	     &event.momLambda1_Lp,		"momLambda1_Lp/D");				
  tree->Branch("momLambda2_Lp", 	     &event.momLambda2_Lp,		"momLambda2_Lp/D");				
  tree->Branch("vertexLambdaDecay_Lp", 	     event.vertexLambdaDecay_Lp,	"vertexLambdaDecay_Lp[3]/D");	
  tree->Branch("cdistLambdaDecay_Lp", 	     &event.cdistLambdaDecay_Lp,	"cdistLambdaDecay_Lp/D");			
  tree->Branch("cdistLambdaPScat_Lambda",    &event.cdistLambdaPScat_Lambda,	"cdistLambdaPScat_Lambda/D");	
  tree->Branch("vertexLambdaPScat_Lambda",   event.vertexLambdaPScat_Lambda,	"vertexLambdaPScat_Lambda[3]/D");	
  tree->Branch("costLambdaPScat_Lambda",     &event.costLambdaPScat_Lambda,	"costLambdaPScat_Lambda/D");	
  tree->Branch("thetaLambdaPScat_Lambda",    &event.thetaLambdaPScat_Lambda,	"thetaLambdaPScat_Lambda/D");	
  tree->Branch("thetaCMLambdaPScat_Lambda",  &event.thetaCMLambdaPScat_Lambda,"thetaCMLambdaPScat_Lambda/D");		
  tree->Branch("momCalLambda1_Lp", 	     &event.momCalLambda1_Lp,		"momCalLambda1_Lp/D");		
  tree->Branch("momCalLambda2_Lp", 	     &event.momCalLambda2_Lp,		"momCalLambda2_Lp/D");		
  tree->Branch("momCalRecoilProton_Lp",      &event.momCalRecoilProton_Lp,	"momCalRecoilProton_Lp/D");			
  tree->Branch("momRecoilProton_Lp", 	     &event.momRecoilProton_Lp,	"momRecoilProton_Lp/D");			
  tree->Branch("deltaThetaRecoilProton_Lp",  &event.deltaThetaRecoilProton_Lp,"deltaThetaRecoilProton_Lp/D");		
  tree->Branch("vertexLambdaPScat_p", 	     event.vertexLambdaPScat_p,	"vertexLambdaPScat_p[3]/D");	
  tree->Branch("cdistLambdaPScat_p", 	     &event.cdistLambdaPScat_p,	"cdistLambdaPScat_p/D");			
  tree->Branch("thetaScatLambdaPScat_p",     &event.thetaScatLambdaPScat_p,	"thetaScatLambdaPScat_p/D");	
  tree->Branch("scatMomCalLambdaPScat_p",    &event.scatMomCalLambdaPScat_p,	"scatMomCalLambdaPScat_p/D");	
  tree->Branch("scatEkinCalLambdaPScat_p",   &event.scatEkinCalLambdaPScat_p,	"scatEkinCalLambdaPScat_p/D");	
  tree->Branch("thetaScatCMLambdaPScat_p",   &event.thetaScatCMLambdaPScat_p,	"thetaScatCMLambdaPScat_p/D");	
  tree->Branch("costScatOpenAngle_Lp", 	     &event.costScatOpenAngle_Lp,	"costScatOpenAngle_Lp/D");			
  tree->Branch("thetaScatOpenAngle_Lp",      &event.thetaScatOpenAngle_Lp,	"thetaScatOpenAngle_Lp/D");			
  tree->Branch("cdistDecayPi2_Lp", 	     &event.cdistDecayPi2_Lp,		"cdistDecayPi2_Lp/D");		
  tree->Branch("vertexDecayPi2_Lp", 	     event.vertexDecayPi2_Lp,		"vertexDecayPi2_Lp[3]/D");		
  tree->Branch("cdistDecayProton1_Lp_2p",    &event.cdistDecayProton1_Lp_2p,	"cdistDecayProton1_Lp_2p/D");	
  tree->Branch("vertexDecayProton1_Lp_2p",   event.vertexDecayProton1_Lp_2p,	"vertexDecayProton1_Lp_2p[3]/D");	
  tree->Branch("MissMassPi_1_Lp_2p", 	     &event.MissMassPi_1_Lp_2p,       "MissMassPi_1_Lp_2p/D");                 

  tree->Branch("momLambda1_Lp_1st",              &event.momLambda1_Lp_1st,		     "momLambda1_Lp_1st/D");		     
  tree->Branch("momLambda2_Lp_1st",		 &event.momLambda2_Lp_1st,		     "momLambda2_Lp_1st/D");		     
  tree->Branch("vertexLambdaDecay_Lp_1st",	 event.vertexLambdaDecay_Lp_1st,	     "vertexLambdaDecay_Lp_1st[3]/D");	     
  tree->Branch("cdistLambdaDecay_Lp_1st",	 &event.cdistLambdaDecay_Lp_1st,	     "cdistLambdaDecay_Lp_1st/D");	     
  tree->Branch("cdistLambdaPScat_Lambda_1st",	 &event.cdistLambdaPScat_Lambda_1st,	     "cdistLambdaPScat_Lambda_1st/D");	     
  tree->Branch("vertexLambdaPScat_Lambda_1st",	 event.vertexLambdaPScat_Lambda_1st,       "vertexLambdaPScat_Lambda_1st[3]/D");    
  tree->Branch("costLambdaPScat_Lambda_1st",	 &event.costLambdaPScat_Lambda_1st,	     "costLambdaPScat_Lambda_1st/D");	     
  tree->Branch("thetaLambdaPScat_Lambda_1st",	 &event.thetaLambdaPScat_Lambda_1st,	     "thetaLambdaPScat_Lambda_1st/D");	     
  tree->Branch("thetaCMLambdaPScat_Lambda_1st",	 &event.thetaCMLambdaPScat_Lambda_1st,     "thetaCMLambdaPScat_Lambda_1st/D");	     	     
  tree->Branch("momCalLambda1_Lp_1st",		 &event.momCalLambda1_Lp_1st,		     "momCalLambda1_Lp_1st/D");		     
  tree->Branch("momCalLambda2_Lp_1st",		 &event.momCalLambda2_Lp_1st,		     "momCalLambda2_Lp_1st/D");		     
  tree->Branch("momCalRecoilProton_Lp_1st",	 &event.momCalRecoilProton_Lp_1st,	     "momCalRecoilProton_Lp_1st/D");	     
  tree->Branch("momRecoilProton_Lp_1st",	 &event.momRecoilProton_Lp_1st,	     "momRecoilProton_Lp_1st/D");	     
  tree->Branch("deltaThetaRecoilProton_Lp_1st",	 &event.deltaThetaRecoilProton_Lp_1st,     "deltaThetaRecoilProton_Lp_1st/D");	     	     
						 				     					     	     
  tree->Branch("momLambda1_Lp_2nd",		 &event.momLambda1_Lp_2nd,		     "momLambda1_Lp_2nd/D");		     
  tree->Branch("momLambda2_Lp_2nd",		 &event.momLambda2_Lp_2nd,		     "momLambda2_Lp_2nd/D");		     
  tree->Branch("vertexLambdaDecay_Lp_2nd",	 event.vertexLambdaDecay_Lp_2nd,	     "vertexLambdaDecay_Lp_2nd[3]/D");	     
  tree->Branch("cdistLambdaDecay_Lp_2nd",	 &event.cdistLambdaDecay_Lp_2nd,	     "cdistLambdaDecay_Lp_2nd/D");	     
  tree->Branch("cdistLambdaPScat_Lambda_2nd",	 &event.cdistLambdaPScat_Lambda_2nd,	     "cdistLambdaPScat_Lambda_2nd/D");	     
  tree->Branch("vertexLambdaPScat_Lambda_2nd",	 event.vertexLambdaPScat_Lambda_2nd,       "vertexLambdaPScat_Lambda_2nd[3]/D");    
  tree->Branch("costLambdaPScat_Lambda_2nd",	 &event.costLambdaPScat_Lambda_2nd,	     "costLambdaPScat_Lambda_2nd/D");	     
  tree->Branch("thetaLambdaPScat_Lambda_2nd",	 &event.thetaLambdaPScat_Lambda_2nd,	     "thetaLambdaPScat_Lambda_2nd/D");	     
  tree->Branch("thetaCMLambdaPScat_Lambda_2nd",	 &event.thetaCMLambdaPScat_Lambda_2nd,     "thetaCMLambdaPScat_Lambda_2nd/D");	     	     
  tree->Branch("momCalLambda1_Lp_2nd",		 &event.momCalLambda1_Lp_2nd,		     "momCalLambda1_Lp_2nd/D");		     
  tree->Branch("momCalLambda2_Lp_2nd",		 &event.momCalLambda2_Lp_2nd,		     "momCalLambda2_Lp_2nd/D");		     
  tree->Branch("momCalRecoilProton_Lp_2nd",	 &event.momCalRecoilProton_Lp_2nd,	     "momCalRecoilProton_Lp_2nd/D");	     
  tree->Branch("momRecoilProton_Lp_2nd",	 &event.momRecoilProton_Lp_2nd,	     "momRecoilProton_Lp_2nd/D");	     
  tree->Branch("deltaThetaRecoilProton_Lp_2nd",	 &event.deltaThetaRecoilProton_Lp_2nd,     "deltaThetaRecoilProton_Lp_2nd/D");       

									
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
