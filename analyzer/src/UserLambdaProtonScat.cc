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
#include "CFTTrackingEffMan.hh"

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
const double GammaMass       = 0.0;
const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double Kaon0Mass       = 0.497614;
const double ProtonMass      = 0.93827200;
const double NeutronMass     = 0.93956563;
const double LambdaMass      = 1.115648;
const double SigmaMinusMass  = 1.197449;
const double SigmaZeroMass   = 1.192642;
////////////////////////////////////////////////////

const double MinOfMassSquareK = 0.1;
const double MaxOfMassSquareK = 0.4;
const double MinOfMassSquarePi = -0.1;
const double MaxOfMassSquarePi = 0.1;
const double MinOfMassSquareP = 0.6;
const double MaxOfMassSquareP = 1.2;

const double MinOfMassSigma = 1.15;
const double MaxOfMassSigma = 1.25;

const double MinOfMassLambda = 1.;
const double MaxOfMassLambda = 1.25;

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
double calcKMomFromThetaLambda(double bmom, double cost, double p0);
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
  double egamma;
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
  double costDecayProtonCM0;
  double momVectorPiMinusFromK0[3];
  double momPiMinusFromK0;

  // Tracking
  int ntSks;
  double p[MaxHits];
  double q[MaxHits];

  int ntSksPart;
  double pSksPart[MaxHits];
  double yTofHyps[MaxHits];
  double xTofHyps[MaxHits];
  double m2[MaxHits];

  int nK;
  double pK[MaxHits];
  double u0SdcIn;
  double v0SdcIn;

  double resSftU;
  double resSftV;

  double MissMass[MaxHits];
  double MissMass_others1[MaxHits];
  double MissMass_others2[MaxHits];
  double HypsInvMass_PiP;
  double HypsInvMass_Phi;
  // double HypsInvMass_PiP;

  int nLambda;
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
  double EkinP_2nd;
  double EkinCorP_2nd;

  double vertexDecayPi[3];
  double cdistDecayPi;
  double thetaDecayPi;
  double decayPiMomCal;
  double decayPiMomVecCal[3];
  double decayPMomCal;
  double decayPMomVecCal[3];
  double costDecayPMomCal;
  double costDecayPMomCalCM;

  double vertexScat[3];
  double hypBeamVec[3];
  double cdistScat;
  double thetaScat;
  double phiScat;
  double scatMomCal;
  double scatEkinCal;
  double thetaScatCM;
  double costScatOpenAngle;
  double thetaScatOpenAngle;

  double cdistDecayPi2;
  double vertexDecayPi2[3];
  double costDecayProtonLambdaPScat_cal;
  double costDecayProtonLambdaPScatCM_cal;

  double vertexScat_2nd[3];
  double hypBeamVec_2nd[3];
  double cdistScat_2nd;
  double thetaScat_2nd;
  double phiScat_2nd;
  double scatMomCal_2nd;
  double scatEkinCal_2nd;
  double thetaScatCM_2nd;

  double vertexDecayProton1_2p[3];
  double cdistDecayProton1_2p;
  double MissMassPi_1_2p;
  double costDecayProtonLambdaPScat1_2p;
  double costDecayProtonLambdaPScatCM1_2p;

  double vertexDecayProton2_2p[3];
  double cdistDecayProton2_2p;
  double MissMassPi_2_2p;
  double costSigmaDecayProton_2p;
  double costDecayProtonLambdaPScat2_2p;
  double costDecayProtonLambdaPScatCM2_2p;

  double scatMomCal_S0pConv1;
  double scatMomCal_S0pConv2;
  double scatEkinCal_S0pConv1;
  double scatEkinCal_S0pConv2;
  double thetaS0pConvCM;
  double costScatOpenAngle_S0pConv;
  double thetaScatOpenAngle_S0pConv;
  double cdistDecayPi2_S0pConv;
  double vertexDecayPi2_S0pConv[3];

  double vertexDecayProton_S0pConv1_2p[3];
  double cdistDecayProton_S0pConv1_2p;
  double MissMassPiGamma_1_2p;

  double scatMomCal_S0pConv1_2nd;
  double scatMomCal_S0pConv2_2nd;
  double scatEkinCal_S0pConv1_2nd;
  double scatEkinCal_S0pConv2_2nd;
  double thetaS0pConvCM_2nd;
  double costScatOpenAngle_S0pConv_2nd;
  double thetaScatOpenAngle_S0pConv_2nd;
  double cdistDecayPi2_S0pConv_2nd;
  double vertexDecayPi2_S0pConv_2nd[3];

  double vertexDecayProton_S0pConv1_2p_2nd[3];
  double cdistDecayProton_S0pConv1_2p_2nd;
  double MissMassPiGamma_1_2p_2nd;



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

  // Missing Mass
  int nPiPi;
  double vtx_pipi_x[MaxHits];
  double vtx_pipi_y[MaxHits];
  double vtx_pipi_z[MaxHits];
  double invmass[MaxHits];
  double theta_pipi[MaxHits];
  double cdist_pipi[MaxHits];
  double cdistBeam_pipi[MaxHits];
  double missmass[MaxHits];
  double missmom[MaxHits];
  double momK0[MaxHits];
  double vtx_k0_x[MaxHits];
  double vtx_k0_y[MaxHits];
  double vtx_k0_z[MaxHits];
  double theta_k0[MaxHits];
  double cdist_k0[MaxHits];
  double momCalLambda[MaxHits];
  double momCalLambda_x[MaxHits];
  double momCalLambda_y[MaxHits];
  double momCalLambda_z[MaxHits];
  double missmass0;
  double invmass0;

  double momLambda1;
  double momLambda2;
  double vertexLambdaDecay[3];
  double cdistLambdaDecay;
  double cdistLambdaPScat_Lambda;
  double vertexLambdaPScat_Lambda[3];
  double costLambdaPScat_Lambda;
  double thetaLambdaPScat_Lambda;
  double thetaCMLambdaPScat_Lambda;
  double phiLambdaPScat_Lambda;
  double momCalLambda1;
  double momCalLambda2;
  double momCalRecoilProton;
  double momRecoilProton;
  double deltaThetaRecoilProton;
  double costDecayProtonLambdaPScat;
  double costDecayProtonLambdaBeam;
  double costDecayProtonLambdaPScatCM;
  double costDecayProtonLambdaBeamCM;

  double momLambda1_1st;
  double momLambda2_1st;
  double vertexLambdaDecay_1st[3];
  double cdistLambdaDecay_1st;
  double cdistLambdaPScat_Lambda_1st;
  double vertexLambdaPScat_Lambda_1st[3];
  double costLambdaPScat_Lambda_1st;
  double thetaLambdaPScat_Lambda_1st;
  double phiLambdaPScat_Lambda_1st;
  double thetaCMLambdaPScat_Lambda_1st;
  double momCalLambda1_1st;
  double momCalLambda2_1st;
  double momCalRecoilProton_1st;
  double momRecoilProton_1st;
  double deltaThetaRecoilProton_1st;
  double costDecayProtonLambdaPScat_1st;
  double costDecayProtonLambdaPScatCM_1st;

  double momLambda1_2nd;
  double momLambda2_2nd;
  double vertexLambdaDecay_2nd[3];
  double cdistLambdaDecay_2nd;
  double cdistLambdaPScat_Lambda_2nd;
  double vertexLambdaPScat_Lambda_2nd[3];
  double costLambdaPScat_Lambda_2nd;
  double thetaLambdaPScat_Lambda_2nd;
  double phiLambdaPScat_Lambda_2nd;
  double thetaCMLambdaPScat_Lambda_2nd;
  double momCalLambda1_2nd;
  double momCalLambda2_2nd;
  double momCalRecoilProton_2nd;
  double momRecoilProton_2nd;
  double deltaThetaRecoilProton_2nd;
  double costDecayProtonLambdaPScat_2nd;
  double costDecayProtonLambdaPScatCM_2nd;

  double vertexK0Production[3];
  double cdistK0Production;
  double thetaK0Production;
  double momK0Production;
  double vertexK0decay[3];
  double cdistK0decay;
  double openingAngleK0decay;
  double momCalLambdaProduction;
  double momCalLambdaProductionVec[3];

  double energyBGO_p_lambda;
  double energyBGO_pi_lambda;
  int    segPiID_pi_lambda;
  double invMass_PiP;
  double invMass_PiP2;
  double theta_PiP;

  double vertexPPScat[3];
  double cdistTwoProton;
  double thetaTwoProton;
  double momLambdaProtonDecay;
  double momLambdaProtonDecayVec[3];
  double vertexLambdaProtonDecay[3];
  double cdistLambdaProtonDecay;
  double thetaLambdaProtonDecay;
  double momCalLambdaProtonDecay1;
  double momCalLambdaProtonDecay2;

  // from decay pi
  double vertexPpScat[3];
  double cdistPpScat;
  double thetaPpScat;
  double scatPpMomCal;
  double scatPpEkinCal;
  double thetaPpScatCM;


  int nPiPi_CFT;
  double vtx_pipi_cft_x[MaxHits];
  double vtx_pipi_cft_y[MaxHits];
  double vtx_pipi_cft_z[MaxHits];
  double theta_cft_pipi[MaxHits];
  double cdist_cft_pipi[MaxHits];
  double cdistBeam_cft_pipi[MaxHits];

  int nPiP_CFT;
  double vtx_pip_cft_x[MaxHits];
  double vtx_pip_cft_y[MaxHits];
  double vtx_pip_cft_z[MaxHits];
  double theta_cft_pip[MaxHits];
  double cdist_cft_pip[MaxHits];
  double cdistBeam_cft_pip[MaxHits];
  double invmass_cft_pip[MaxHits];

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

  // newly added for pi- momentum from K0
  /*
  In.getline(line, sizeof(line));
  sscanf(line,"momVectorPiMinusFromK0 %lf %lf %lf",
	 &event.momVectorPiMinusFromK0[0],
	 &event.momVectorPiMinusFromK0[1],
	 &event.momVectorPiMinusFromK0[2]);
  */

  event.momPiMinusFromK0 =
    sqrt(event.momVectorPiMinusFromK0[0]*event.momVectorPiMinusFromK0[0] +
	 event.momVectorPiMinusFromK0[1]*event.momVectorPiMinusFromK0[1] +
	 event.momVectorPiMinusFromK0[2]*event.momVectorPiMinusFromK0[2]);
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

  // primary analysis for polarization measurement
  {
    ThreeVector NormX, NormY, NormZ;
    ThreeVector LambdaMom0(event.momVectorHypBeam[0],
			   event.momVectorHypBeam[1],
			   event.momVectorHypBeam[2]);
    ThreeVector K0Mom0(event.momVectorScatMeson[0],
		       event.momVectorScatMeson[1],
		       event.momVectorScatMeson[2]);
    ThreeVector ProtonMom0(event.momVectorDecayNucleon[0],
			   event.momVectorDecayNucleon[1],
			   event.momVectorDecayNucleon[2]);
    NormZ = LambdaMom0/(LambdaMom0.mag());
    NormY = gaiseki(LambdaMom0, K0Mom0);
    NormY *= (1./NormY.mag());
    NormX = gaiseki(NormY, NormZ);
    NormX *= (1./NormX.mag());
    double E = std::sqrt( LambdaMass*LambdaMass+LambdaMom0.mag2());
    ThreeVector BetaLambdaBeam = (LambdaMom0/LambdaMom0.mag())*(LambdaMom0.mag()/E);
    LorentzVector LvDecayProtonVec( ProtonMom0, std::sqrt( ProtonMass*ProtonMass+ProtonMom0.mag2()));
    LorentzVector LvDecayProtonVecCM = LvDecayProtonVec.boost(-BetaLambdaBeam);
    event.costDecayProtonCM0 =  LvDecayProtonVecCM.vect()*NormY/(LvDecayProtonVecCM.vect().mag()*NormY.mag());
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
    if (layer >= PlMinSdcIn && layer <= PlMaxSdcIn){
      rawData->AddSdcInRawHit(layer, segment, dtime);
      // std::cout << "layer:" << layer << " segment:" << segment << " dtime:" << dtime << std::endl;
      HF1( 600+layer, segment );
    } else if (layer >= PlMinSdcOut && layer <= PlMaxSdcOut){
      rawData->AddSdcOutRawHit(layer, segment, dtime);
      // std::cout << "layer:" << layer << " segment:" << segment << " dtime:" << dtime << std::endl;
      HF1( 700+layer, segment );
      // std::cout << "layer:" << layer << " segment:" << segment << " dtime:" << dtime << std::endl;
    } else if (layer >= PlOffsBc+PlMinBcOut && layer <= PlOffsBc+PlMaxBcOut)
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
    } // else if  (layer == 10) {
    //   rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 0, 0, edep);
    //   rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 1, 0, time);
    //   rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 0, 1, edep);
    //   rawData->AddCHHodoRawHit(DetIdCH, layer, segment, 1, 1, time);
    // }
    else if (layer == 70) {
      int Layer = segment/NumOfSegPiV;
      int Segment = segment%NumOfSegPiV;

      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 0, 0, edep);
      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 1, 0, time);
      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 0, 1, edep);
      rawData->AddPiVHodoRawHit(DetIdPiV, Layer, Segment, 1, 1, time);

      if (FlagEvDisp) {
	if (std::abs(time)<30 && edep>0.1)
	  std::cout << "PiV Seg : " << Segment << ", Layer = " << Layer << ", time : " << time << ", dE : "
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
      double time = hit->MeanTime();

      if (std::abs(time)<30 && dE>0.1) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.ShowHitPiV(layer, segment, dE);
      }
    }
  }


  int ncTof=hodoAna->GetNClustersTOF();

  DCAnalyzer *DCAna = new DCAnalyzer;
  DCAna->DecodeSdcInHits(rawData);
  DCAna->DecodeSdcOutHits(rawData);
  DCAna->DecodeBcOutHits(rawData);

  DCAna->TrackSearchSdcInHyps();
  // DCAna->TrackSearchSdcInFiber();
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

    for( int ih=0; ih<nh; ++ih ){
      DCLTrackHit *hit=tp->GetHit(ih);
      int layerId=hit->GetLayer();

      HF1( 13, layerId );
      double wire=hit->GetWire();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1( 140+layerId, res );
      HF1( 160+layerId, pos );
    }

  }

  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof );

  DCAna->TrackSearchSdcOutHyps();
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
  // DCAna->TrackSearchSksTmp(event.momScatMeson);
  DCAna->TrackSearchHypsTmp(event.momScatMeson);
  int ntSks = DCAna->GetNTracksSks();
  event.ntSks = ntSks;

  std::vector <SksParticle *> SksPartCont;
  std::vector<ThreeVector> HypsMom;
  std::vector<Double_t> HypsM2;
  std::vector<Double_t> HypsPol;
  std::vector<int> good;
  bool okk1;
  bool okk2;

  HF1( 50, double(ntSks) );

  // if(ntSks==1) okk1 = true;
  // if(ntSks==2) okk2 = true;


  for( int i=0; i<ntSks; ++i ){
    SksTrack *tp=DCAna->GetSksTrack(i);
    if(!tp) continue;
    bool ok = false;
    int nh=tp->GetNHits();
    double chisqr=tp->chisqr();
    ThreeVector Ppos=tp->PrimaryPosition();
    ThreeVector Pmom=(-1)*tp->PrimaryMomentum();
    ThreeVector bPos(0, 0, 0);
    ThreeVector bMom(0, 0, 1);
    ThreeVector vtx = VertexPoint(bPos, Ppos, bMom, Pmom);

    double polarity=tp->polarity();
    ThreeVector Mom = (-1.)*Pmom;
    // Double_t ut_ = Mom.x()/Mom.z(), vt_ = Mom.y()/Mom.z();
    // Double_t pt = polarity/std::sqrt(1+ut_*ut_+vt_*vt_);
    // ThreeVector Hyps_mom(pt*ut_, pt*vt_, pt);
    // std::cout << "ntSks = " << ntSks << ", i = "  << i << ", Hyps_mom = " << Hyps_mom << std::endl;
    // HypsMom.push_back(Mom);
    // HypsPol.push_back(polarity);

    double pathL=tp->PathLengthToTOF() + (0-vtx.z())*Pmom.mag()/Pmom.z();
    tp->SetPathLengthToTOF(pathL);
    //double xt=Ppos.x(), yt=Ppos.y();
    double xt=vtx.x(), yt=vtx.y();
    double p=Pmom.mag();
    double ut=Pmom.x()/p, vt=Pmom.y()/p;
    double cost = 1./sqrt(1.+ut*ut+vt*vt);
    double theta = acos(cost)*Rad2Deg;
    double phi   = atan2( ut, vt );

    event.p[i] = p;
    event.q[i] = -polarity;

    HF1( 51, double(nh) );
    HF1( 52, chisqr );
    HF1( 54, xt ); HF1( 55, yt ); HF1( 56, ut ); HF1( 57,vt );
    HF2( 58, xt, ut ); HF2( 59, yt, vt ); HF2( 60, xt, yt );
    HF1( 61, p ); HF1( 62, pathL );
    HF1(63, theta);

    for( int j=0; j<nh; ++j ){
      TrackHit *hit=tp->GetHit(j);
      if(!hit) continue;
      int layerId = hit->GetLayer();
      HF1( 53, hit->GetLayer() );
      double res  = hit->GetResidual();
      HF1( 100+layerId, res );
      double pos  = hit->GetLocalHitPos();
      HF1( 150+layerId, pos );
    }
    if (chisqr>MaxChiSqrSks) continue;

    DCLocalTrack *tout = tp->GetLocalTrackOut();
    double xtof = tout->GetFTOFX();
    std::cout << "ncTof" << ncTof << std::endl;
    for (int j=0; j<ncTof; ++j) {
      HodoCluster *clTof=hodoAna->GetClusterTOF(j);
      if (!clTof || !clTof->GoodForAnalysis() ) continue;

      double meanSeg = clTof->MeanSeg();
      double tofSegPos = DCGeomMan::GetInstance().calcWirePosition(IdTof, meanSeg);

      double diffX = xtof-tofSegPos;
      HF1(64, diffX);
      std::cout << "i = " << i << ", diffX = " << diffX <<std::endl;
      // if (!(fabs(diffX)<6000.)) ok = false;
      // if (fabs(diffX)<60.) {
      if (fabs(diffX)<700.) {
      // if (fabs(diffX)<6000.) {
	SksParticle *SksPart = new SksParticle(tp, clTof, tout);
	SksPartCont.push_back(SksPart);

	double sigma_time = gRandom->Gaus(0, 0.15);
	double m2 = MassSquare(p, pathL, clTof->CMeanTime() + sigma_time);
	HF1(65, m2);
	HypsM2.push_back(m2);
	HypsMom.push_back(Mom);
	HypsPol.push_back(polarity);
	std::cout << "i = " << i << ", m2 = " << m2 << ", Mom = " << Mom << ", polaroty" << polarity << std::endl;
	ok =true;
	//if (m2>0)
	if (((int)meanSeg)%2==1)
	  HF1(66, m2);
	  //HF1(66, sqrt(m2));
      }
      if(ok) good.push_back(i);
      ok =false;
    }
  }

  int ntSksPart = SksPartCont.size();
  event.ntSksPart = ntSksPart;
  for (int i=0; i<ntSksPart; i++) {
    event.pSksPart[i] = SksPartCont[i]->Momentum().mag();
    event.m2[i] = SksPartCont[i]->MassSquare();
    event.yTofHyps[i] = SksPartCont[i]->yTof();
    event.xTofHyps[i] = SksPartCont[i]->xTof();
  }

  std::cout << "good size = " << good.size() << std::endl;
  std::cout << "HypsM2.size = " << HypsM2.size() << std::endl;
  if(good.size()==2){
    // SksParticle *t_a = SksPartCont[0];
    // double m2_a = SksPart_a->MassSquare();
    // double Pol_a = SksPart_a->Polarity();
    // ThreeVector mom_a = (-1.)*SksPartCont[0]->GetTrack()->PrimaryMomentum();
    // SksParticle *SksPart_b = SksPartCont[1];
    // double m2_b = SksPart_b->MassSquare();
    // double Pol_b = SksPart_b->Polarity();
    // ThreeVector mom_b = (-1.)*SksPartCont[1]->GetTrack()->PrimaryMomentum();
    ThreeVector mom_a = HypsMom[0];
    ThreeVector mom_b = HypsMom[1];
    Double_t m2_a = HypsM2[0];
    Double_t m2_b = HypsM2[1];
    Double_t Pol_a = HypsPol[0];
    Double_t Pol_b = HypsPol[1];
    ThreeVector Mom1(-999.9, -999.9, -999.9);
    ThreeVector Mom2(-999.9, -999.9, -999.9);
    bool Flag = false;
    bool Flag_PiP = false;
    bool Flag_Phi = false;
    std::cout << "mom_a = " << mom_a << ", mom_b = " << mom_b << std::endl;
    std::cout << "m2_a = " << m2_a << ", m2_b = " << m2_b << std::endl;
    // std::cout << "Pol_a = " << Pol_a << ", Pol_b = " << Pol_b << std::endl;
    if(Pol_a*Pol_b==-1) Flag = true;
    if(Flag){
      // if(true){
      std::cout << "OK!OK!" << "mom_a = " << mom_a << ", mom_b = " << mom_b << std::endl;
      if(std::isnan(m2_a)){
	Flag = false;
	// invMass_PiP = qnan;
	// invMass_PiP_select = qnan;
	// invMass_Phi = qnan;
	// invMass_Phi_select = qnan;
      }else if(std::isnan(m2_b)){
	Flag = false;
	// invMass_PiP = qnan;
	// invMass_PiP_select = qnan;
	// invMass_Phi = qnan;
	// invMass_Phi_select = qnan;
      }
      // std::cout << "Flag is OK!" << std::endl;
      if(Flag){
	if(m2_b>MinOfMassSquarePi && m2_b<MaxOfMassSquarePi && Pol_b==1){
	// if(m2_b>MinOfMassSquarePi && m2_b<MaxOfMassSquarePi && m2_a>MinOfMassSquareP && m2_a<MaxOfMassSquareP){
	  // if(m2_a>MinOfMassSquareP && m2_a<MaxOfMassSquareP){
	  // if(m2_a < m2_b){
	  // Mom2 = Pol_a*mom_a;
	  // Mom1 = Pol_b*mom_b;
	  Mom2 = (-1.)*mom_a;
	  Mom1 = (-1.)*mom_b;
	  Flag_PiP = true;
	  // std::cout << "first Mom1 = " << Mom1 << ", Mom2 = " << Mom2 << std::endl;
	  // }

	}else if(m2_a>MinOfMassSquarePi && m2_a<MaxOfMassSquarePi && Pol_a==1){
	// }else if(m2_a>MinOfMassSquarePi && m2_a<MaxOfMassSquarePi && m2_b>MinOfMassSquareP && m2_b<MaxOfMassSquareP){
	  // if( m2_b>MinOfMassSquareP && m2_b<MaxOfMassSquareP){
	  // }else if(m2_a > m2_b){
	  // Mom2 = Pol_b*mom_b;
	  // Mom1 = Pol_a*mom_a;
	  Mom2 = (-1.)*mom_b;
	  Mom1 = (-1.)*mom_a;
	  // Mom2 = (-1.)*SksPart_b->Momentum();
	  // Mom1 = (1.)*SksPart_a->Momentum();
	  Flag_PiP = true;
	  // std::cout << "second Mom1 = " << Mom1 << ", Mom2 = " << Mom2 << std::endl;
	  // }
	}else if(m2_a>MinOfMassSquareK && m2_a<MaxOfMassSquareK && m2_b>MinOfMassSquareK && m2_b<MaxOfMassSquareK){
	  // if(m2_b>MinOfMassSquareK && m2_b<MaxOfMassSquareK){
	  Mom2 = (-1.)*mom_b;
	  Mom1 = (-1.)*mom_a;
	  // Mom1 = mom_b;
	  // Mom2 = mom_a;
	  Flag_Phi = true;
	  // std::cout << "Phi Mom1 = " << Mom1 << ", Mom2 = " << Mom2 << std::endl;
	  // }
	}
      }
    // std::cout << "ntSksPart = " << ntSksPart << std::endl;
    // std::cout << "final Mom1 = " << Mom1 << ", Mom2 = " << Mom2 << std::endl;
      if(Flag_PiP){
      LorentzVector LvPi( Mom1, std::sqrt( PionMass*PionMass+Mom1.mag2()) );
      LorentzVector LvP( Mom2, std::sqrt( ProtonMass*ProtonMass+Mom2.mag2()) );
      LorentzVector LvPiP = LvPi + LvP;
      double InvMass = LvPiP.mag();
      double PiMass = LvP.mag();
      double ProtonMass = LvP.mag();
      event.HypsInvMass_PiP = InvMass;
      // event.HypsInvMass_PiP = Mom2.mag2();
      HF1(69, InvMass);
      std::cout << "InvMass_PiP = " << InvMass << std::endl;
      }
      if(Flag_Phi){
      LorentzVector LvPi( Mom1, std::sqrt( KaonMass*KaonMass+Mom1.mag2()) );
      LorentzVector LvP( Mom2, std::sqrt( KaonMass*KaonMass+Mom2.mag2()) );
      LorentzVector LvPiP = LvPi + LvP;
      double InvMass = LvPiP.mag();
      event.HypsInvMass_Phi = InvMass;
      HF1(70, InvMass);
      }
    }
  }


  // DCAna->TrackSearchBcOut();
  // int ntBcOut = DCAna->GetNtracksBcOut();

  std::vector <ThreeVector> BmomCont, BposCont;

  double beam_energy_spread = 0.0085; // GeV
  double egamma = event.beammom+gRandom->Uniform(-beam_energy_spread, beam_energy_spread);
  ThreeVector bmom(0, 0, egamma);
  double theta = gRandom->Uniform(0, 0.00006);
  double phi = gRandom->Uniform(0, 360.*Deg2Rad);
  bmom.setTheta(theta);
  bmom.setPhi(phi);
  event.egamma = egamma;

  double r_sigma = 3.0; // mm
  double r_beam = gRandom->Gaus(0, r_sigma);
  if (r_beam<0)
    r_beam *= -1.0;
  double theta_beam = gRandom->Uniform(0, 360.*Deg2Rad);

  double x_beam = r_beam*cos(theta_beam);
  double y_beam = r_beam*sin(theta_beam);

  ThreeVector bpos(x_beam, y_beam, 0);


  BmomCont.push_back(bmom);
  BposCont.push_back(bpos);

  int nBeam=BmomCont.size();

  if (nBeam != 1 || ntSksPart != 1) {
    tree->Fill();

    for_each(SksPartCont.begin(), SksPartCont.end(), DeleteObject());

    delete hodoAna;
    delete DCAna;
    delete rawData;

    if (FlagEvDisp) {
      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      /*
      if (1) {
	//if (event.scatFlag==1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1&&nPi_CFT==1) {
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

  // Kaon means Track in Kurama
  SksParticle *KaonPart = SksPartCont[0];
  event.u0SdcIn = KaonPart->GetTrack()->GetLocalTrackIn()->GetU0();
  event.v0SdcIn = KaonPart->GetTrack()->GetLocalTrackIn()->GetV0();


  if (FlagEvDisp) {

    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    //evDisp.DrawVertex(event.Vertex_x[0], event.Vertex_y[0], event.Vertex_z[0]);

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
  HF1(401, ntCFT);

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

      HF1(437, dist);
      HF1(438, layer, dist);

      double u=Dir.x(), v=Dir.y();
      double x0=Pos0.x(), y0=Pos0.y();
      double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);

      if (t>=0) {
	//if (fabs(dist)<30 && time>0 && time<5) {
	if (fabs(dist)<80 && time>0 && time<5) {
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
      event.zDirCFT[it] = Dir.z();
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
    else if (CFTPart->GetMass() > 0.9) {
      bool flagTracking = true;
      ThreeVector momScatP(event.momVectorProtonScat[0], event.momVectorProtonScat[1], event.momVectorProtonScat[2]);
      ThreeVector momDecayP(event.momVectorDecayNucleon[0], event.momVectorDecayNucleon[1], event.momVectorDecayNucleon[2]);
      double costScat = -999.;
      if (momScatP.x()>-990.) {
	costScat = (Dir*momScatP)/(Dir.mag()*momScatP.mag());
      }
      double costDecay = -999.;
      if (momDecayP.x()>-990.) {
	costDecay = (Dir*momDecayP)/(Dir.mag()*momDecayP.mag());
      }
      if (costScat>0.995) {
	double theta = momScatP.theta()*Rad2Deg;
	double p=momScatP.mag();
	double vtz = event.scatPos0[2];
	flagTracking = CFTTrackingEffMan::GetInstance().CheckTracking_P(vtz, theta, p);
      } else if (costDecay>0.995) {
	double theta = momDecayP.theta()*Rad2Deg;
	double p=momDecayP.mag();
	double vtz = event.scatPos0[2];
	flagTracking = CFTTrackingEffMan::GetInstance().CheckTracking_P(vtz, theta, p);
      }
      if (flagTracking)
	CFTProtonCont.push_back(CFTPart);
    } else if (CFTPart->GetMass() > 0.0 && CFTPart->GetMass() < 0.2) {
      bool flagTracking = true;
      ThreeVector momDecayPi(event.momVectorDecayPi[0], event.momVectorDecayPi[1], event.momVectorDecayPi[2]);
      double costDecay = -999.;
      if (momDecayPi.x()>-990.) {
	costDecay = (Dir*momDecayPi)/(Dir.mag()*momDecayPi.mag());
      }
      if (costDecay>0.995) {
	double theta = momDecayPi.theta()*Rad2Deg;
	double p=momDecayPi.mag();
	double vtz = event.scatPos0[2];
	flagTracking = CFTTrackingEffMan::GetInstance().CheckTracking_Pi(vtz, theta, p);
      }

      if (flagTracking)
	CFTPionCont.push_back(CFTPart);
    }

  }

  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    int nCFTVtx = CFTVtxCont.size();
    for (int iv=0; iv<nCFTVtx; iv++) {
      ThreeVector CFTVtx = CFTVtxCont[iv];
      evDisp.DrawCFTVertex(CFTVtx.x(), CFTVtx.y(), CFTVtx.z());
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
      event.zDirCFT_P[it] = Dir.z();
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
      event.zDirCFT_Pi[it] = Dir.z();
      event.xPos0CFT_Pi[it] = Pos0.x();
      event.yPos0CFT_Pi[it] = Pos0.y();
      event.zPos0CFT_Pi[it] = Pos0.z();

      event.zBGOCFT_Pi[it] = zBGO;

      event.CFT_TotalEdep_Pi[it] = CFTPionCont[it]->GetFiberTotal_E();
      event.CFT_NormTotalEdep_Pi[it] = CFTPionCont[it]->GetNormFiberTotal_E();
      event.BGO_Edep_Pi[it] = CFTPionCont[it]->GetBGO_E();
      event.PiV_Edep_Pi[it] = CFTPionCont[it]->GetPiV_E();
      event.TotalEdep_Pi[it] = CFTPionCont[it]->GetTotalE();
      if (FlagEvDisp) {
	std::cout << "Emeasure = "
		  << event.TotalEdep_Pi[it] + event.PiV_Edep_Pi[it]
		  << std::endl;
      }

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

  std::vector <PartVec> LambdaVecCont;
  std::vector <SksParticle *> KaonPartCont;

  int nK=0;

  // Missing Mass Analysis
  int nSksPart = SksPartCont.size();
  int nLambda=0;
  for (int i=0; i<nSksPart; i++) {
    SksParticle *SksPart = SksPartCont[i];
    // double m2 = SksPart->MassSquare();
    double m2 = HypsM2[i];
    if (m2>MinOfMassSquareK && m2<MaxOfMassSquareP) {
    // if (m2>MinOfMassSquareK && m2<MaxOfMassSquareK) {
      // ThreeVector priMom = (-1.)*SksPart->Momentum();
      ThreeVector priMom = (-1.)*HypsMom[i];

      LorentzVector LvBeam(  bMom, bMom.mag() );
      LorentzVector LvTgt( 0., 0., 0., ProtonMass );
      LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
      LorentzVector LvRc = LvBeam+LvTgt-LvScat;
      double mismass = LvRc.mag();
      HF1(67, mismass);
      event.MissMass[nK] = mismass;
      event.pK[nK] = priMom.mag();

      nK++;

      if (mismass>MinOfMassLambda && mismass<MaxOfMassLambda) {
	KaonPartCont.push_back(SksPart);

	ThreeVector pPos=SksPart->Position();
	ThreeVector vtx = VertexPoint(bPos, pPos, bMom, priMom);
	event.Vertex_x[nLambda] = vtx.x();
	event.Vertex_y[nLambda] = vtx.y();
	event.Vertex_z[nLambda] = vtx.z();

	double cost = bMom*priMom/(bMom.mag()*priMom.mag());
	event.theta[nLambda] = acos(cost)*Rad2Deg;


	double p0 = priMom.mag();
	double p_Cor, E_Cor;
	CorrElossOut(&p_Cor, &E_Cor, p0, KAON, priMom/p0, vtx);
	ThreeVector KaonMomCor = priMom* p_Cor/p0;
	ThreeVector pSigmaCor = bMom-KaonMomCor;

	double calMom = calcKMomFromThetaLambda(event.beammom, cost, priMom.mag());
	event.pKCal[nLambda] = calMom;
	ThreeVector KaonMomCal = priMom*(calMom/priMom.mag());
	ThreeVector pSigmaCal = bMom-KaonMomCal;
	event.pSigmaCal[nLambda] = pSigmaCal.mag();
	event.momCalLambda[nLambda] = pSigmaCal.mag();
	event.momCalLambda_x[nLambda] = pSigmaCal.x();
	event.momCalLambda_y[nLambda] = pSigmaCal.y();
	event.momCalLambda_z[nLambda] = pSigmaCal.z();
	event.missmass[nLambda] = mismass;

	PartVec lambdaVec;
	lambdaVec.pos0 = vtx;
	lambdaVec.mom  = pSigmaCal; // original
	//sigmaVec.mom  = pSigmaCal*event.momHypBeam/pSigmaCal.mag();
	//sigmaVec.mom  = ThreeVector(event.momVectorHypBeam[0],event.momVectorHypBeam[1],event.momVectorHypBeam[2]);
	//sigmaVec.mom  = pSigmaCor;
	LambdaVecCont.push_back(lambdaVec);

	nLambda++;
      }
    }
  }
  event.nK = nK;
  event.nLambda = nLambda;

  // std::cout << "HypsM2[0] = " << HypsM2[0] << ", HypsM2[1] = " << HypsM2[1] << std::endl;
  // std::cout << "HypsPol[0] = " << HypsPol[0] << ", HypsPol[1] = " << HypsPol[1] << std::endl;
  // if(HypsM2.size()==2){
  // if(nSksPart==2){
  // if(good.size()==2){

  // if(good.size()==1){
  // if(okk1){
  std::cout << "evenrt.ntSksPart2 = " << event.ntSksPart << std::endl;
  if(event.ntSksPart==1){
    Double_t Pol = HypsPol[0];
    Double_t m2 = Pol*HypsM2[0];
    ThreeVector priMom = HypsMom[0];
    LorentzVector LvBeam(  bMom, bMom.mag() );
    LorentzVector LvTgt( 0., 0., 0., ProtonMass );
    LorentzVector LvScat( priMom, sqrt(KaonMass*KaonMass+priMom.mag2()) );
    LorentzVector LvRc = LvBeam+LvTgt-LvScat;
    double mismass = LvRc.mag();
    std::cout << "missmass = " << mismass << std::endl;
    if (m2>MinOfMassSquarePi && m2<MaxOfMassSquarePi) {
      // ThreeVector priMom = (1.)*SksPart->Momentum();
      event.MissMass_others1[nK] = mismass;
    }else if (m2>MinOfMassSquareP && m2<MaxOfMassSquareP) {
      // ThreeVector priMom = (-1.)*SksPart->Momentum();
      event.MissMass_others2[nK] = mismass;
      // event.pK[nK] = priMom.mag();
    }
  }

  ThreeVector LambdaMissMom(0, 0, 0);
  ThreeVector MomK(0, 0, 0);
  ThreeVector VertexK(0, 0, 0);
  if (LambdaVecCont.size() == 1) {
    LambdaMissMom = LambdaVecCont[0].mom;
    MomK = (-1.)*KaonPartCont[0]->Momentum();
    VertexK = ThreeVector(event.Vertex_x[0], event.Vertex_y[0], event.Vertex_z[0]);
    event.momCalLambdaProduction = event.momCalLambda[0];
    event.momCalLambdaProductionVec[0] = event.momCalLambda_x[0];
    event.momCalLambdaProductionVec[1] = event.momCalLambda_y[0];
    event.momCalLambdaProductionVec[2] = event.momCalLambda_z[0];
  }

  bool flagLambdaPScat = false;

  ThreeVector NormX, NormY, NormZ;
  // NormZ : direction of Lambda beam
  // NormY : normal vector of(pi, K0) plane
  // NormX : vector product of Y x Z
  if (LambdaVecCont.size() == 1 ) {

    NormZ = LambdaMissMom/(LambdaMissMom.mag());
    NormY = gaiseki(LambdaMissMom, MomK);
    NormY *= (1./NormY.mag());
    NormX = gaiseki(NormY, NormZ);
    NormX *= (1./NormX.mag());

  }
  ThreeVector BetaLambdaBeam;
  if (LambdaVecCont.size() == 1 ) {
    double E = std::sqrt( LambdaMass*LambdaMass+LambdaMissMom.mag2());
    BetaLambdaBeam = (LambdaMissMom/LambdaMissMom.mag())*(LambdaMissMom.mag()/E);
  }

  bool    flagDecayCal = false;
  PartVec PiDecayVec;
  PartVec DecayProtonVec;

  // Lambda decay assumption
  if (LambdaVecCont.size() == 1 ) {
    if (nPi_CFT ==1) {
      for (int i1 = 0; i1 < nPi_CFT; i1++ ) {
        CFTParticle *PiCFT = CFTPionCont[i1];

	PiDecayVec.pos0 = PiCFT->GetPos0();
	PiDecayVec.mom = PiCFT->GetDir();

	double cdistDecayPi;
	ThreeVector VertLambdaDecay = VertexPoint3D(VertexK, PiDecayVec.pos0,
						    LambdaMissMom, PiDecayVec.mom,
						    cdistDecayPi);
	double costDecayPi=LambdaMissMom*PiDecayVec.mom/(LambdaMissMom.mag()*PiDecayVec.mom.mag());
	double thetaDecayPi=acos(costDecayPi)*Rad2Deg;

	event.vertexDecayPi[0] = VertLambdaDecay.x();
	event.vertexDecayPi[1] = VertLambdaDecay.y();
	event.vertexDecayPi[2] = VertLambdaDecay.z();

	event.cdistDecayPi = cdistDecayPi;
	event.thetaDecayPi = thetaDecayPi;

	double decayPiMomCal;

	flagDecayCal = calcDecayPiMom( LambdaMissMom.mag(),
				       LambdaMass,
				       ProtonMass,
				       PionMass,
				       costDecayPi,
				       &decayPiMomCal);

	if (flagDecayCal) {
	  event.decayPiMomCal = decayPiMomCal;
	  event.decayPiMomVecCal[0] = decayPiMomCal*PiDecayVec.mom.x()/PiDecayVec.mom.mag();
	  event.decayPiMomVecCal[1] = decayPiMomCal*PiDecayVec.mom.y()/PiDecayVec.mom.mag();
	  event.decayPiMomVecCal[2] = decayPiMomCal*PiDecayVec.mom.z()/PiDecayVec.mom.mag();

	  PiDecayVec.mom *= (decayPiMomCal/PiDecayVec.mom.mag());

	  DecayProtonVec.mom = LambdaMissMom-PiDecayVec.mom;
	  DecayProtonVec.pos0 = VertLambdaDecay;

	  event.decayPMomCal = DecayProtonVec.mom.mag();
	  event.decayPMomVecCal[0] = DecayProtonVec.mom.x();
	  event.decayPMomVecCal[1] = DecayProtonVec.mom.y();
	  event.decayPMomVecCal[2] = DecayProtonVec.mom.z();

	  event.costDecayPMomCal =  DecayProtonVec.mom*NormY/(DecayProtonVec.mom.mag()*NormY.mag());
	  LorentzVector LvDecayProtonVec( DecayProtonVec.mom, std::sqrt( ProtonMass*ProtonMass+DecayProtonVec.mom.mag2()));
	  //std::cout << "1 ProtonVec(lab) = " << DecayProtonVec.mom << std::endl;
	  //std::cout << "1 betaVec = " << BetaLambdaBeam << std::endl;
	  LorentzVector LvDecayProtonVecCM = LvDecayProtonVec.boost(-BetaLambdaBeam);
	  //std::cout << "1 ProtonVec(CM) = " << LvDecayProtonVecCM.vect() << std::endl;
	  event.costDecayPMomCalCM =  LvDecayProtonVecCM.vect()*NormY/(LvDecayProtonVecCM.vect().mag()*NormY.mag());
	}
      }
    }
  }

  if (LambdaVecCont.size() == 1 ) {
    if (nPi_CFT ==1 && nP_CFT ==1) {
      for (int i1 = 0; i1 < nPi_CFT; i1++ ) {

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

	  double cdistBeanAndVertex = -999.;
	  ThreeVector ppi = bMom;
	  ThreeVector xpi = bPos;

	  cdistBeanAndVertex =
	    closeDistLineAndPoint(xpi, ppi, vertex);

          event.energyBGO_p_lambda = ProtonCFT2->GetBGO_E();
          event.energyBGO_pi_lambda = PiCFT1->GetBGO_E();
          //event.segPiID_pi_lambda = PiCFT1->GetTrackPiIDSeg();

          double cost = Dir1*Dir2/(Dir1.mag()*Dir2.mag());
          double theta = std::acos(cost)*Rad2Deg;
	  event.theta_PiP = theta;

          double totE_pi = PiCFT1->GetTotalE() + PiCFT1->GetPiV_E();
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
          event.invMass_PiP = InvMass;

	  if (FlagEvDisp) {
	    std::cout << "invmass_pip = " << InvMass << std::endl;
	  }

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
          event.momLambda1 = (Mom1+Mom2).mag();
          event.momLambda2 = momVecLambda.mag();
          event.vertexLambdaDecay[0] = vertex.x();
          event.vertexLambdaDecay[1] = vertex.y();
          event.vertexLambdaDecay[2] = vertex.z();
          event.cdistLambdaDecay     = cdist;

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
              VertexPoint3D( VertexK, vertex,
			     LambdaMissMom, momVecLambda, cdistLambdaPScat_Lambda);

            double costLambdaPScat_Lambda=LambdaMissMom*momVecLambda/
              (LambdaMissMom.mag()*momVecLambda.mag());
            double thetaLambdaPScat_Lambda=acos(costLambdaPScat_Lambda)*Rad2Deg;

	    // phi angle for polarization
	    double p_NormX1 = momVecLambda*NormX;
	    double p_NormY1 = momVecLambda*NormY;
	    double p_NormZ1 = momVecLambda*NormZ;
	    ThreeVector pvec_Norm1(p_NormX1, p_NormY1, p_NormZ1);
	    double theta_norm1, phi_norm1;
	    calcThetaPhi(pvec_Norm1, &theta_norm1, &phi_norm1);
	    if (phi_norm1>270)
	      phi_norm1 -= 360.;

	    // Proton direction in L production plane
	    double costDecayProtonLambdaBeam =
	      Mom2*NormY/(Mom2.mag()*NormY.mag());

	    LorentzVector LvDecayProtonVec( Mom2,
					    std::sqrt( ProtonMass*ProtonMass+Mom2.mag2()));

	    /*
	    ThreeVector Mom0(event.momVectorDecayNucleon[0],
			     event.momVectorDecayNucleon[1],
			     event.momVectorDecayNucleon[2]);
	    LorentzVector LvDecayProtonVec( Mom0, std::sqrt( ProtonMass*ProtonMass+Mom0.mag2()));
	    */
	    /*
	    std::cout << "2 ProtonVec(lab) = " << Mom2 << std::endl;
	    std::cout << "2 betaVec = " << BetaLambdaBeam << std::endl;
	    */
	    LorentzVector LvDecayProtonVecCM = LvDecayProtonVec.boost(-BetaLambdaBeam);
	    //std::cout << "2 ProtonVec(CM) = " << LvDecayProtonVecCM.vect() << std::endl;
	    double costDecayProtonLambdaBeamCM =
	      LvDecayProtonVecCM.vect()*NormY/(LvDecayProtonVecCM.vect().mag()*NormY.mag());

	    /*
	    // test
	    double theta1, phi1;
	    calcThetaPhi(momVecLambda, &theta1, &phi1);
	    ThreeVector Xdir(1, 0, 0);
	    ThreeVector Ydir(0, 1, 0);
	    ThreeVector Zdir(0, 0, 1);
	    Xdir.rotateY(theta1*Deg2Rad);
	    Xdir.rotateZ(phi1*Deg2Rad);
	    Ydir.rotateY(theta1*Deg2Rad);
	    Ydir.rotateZ(phi1*Deg2Rad);
	    Zdir.rotateY(theta1*Deg2Rad);
	    Zdir.rotateZ(phi1*Deg2Rad);
	    std::cout << "theta1 = " << theta1
		      << "phi1 = " << phi1 << std::endl;
	    std::cout << "Xdir = " << Xdir << std::endl;
	    std::cout << "Ydir = " << Ydir << std::endl;
	    std::cout << "Zdir = " << Zdir << std::endl;
	    */
	    // Proton direction in Lp scattering plane
	    ThreeVector NormY_LpScatPlane = gaiseki(NormZ, momVecLambda);
	    NormY_LpScatPlane *= (1./NormY_LpScatPlane.mag());
	    /*
	    std::cout << "momVecLambda(norm) = " << momVecLambda*(1/momVecLambda.mag()) << std::endl;
	    std::cout << "NormY_LpScatPlane = " << NormY_LpScatPlane << std::endl;
	    double cost1 = Ydir*NormY_LpScatPlane/(Ydir.mag()*NormY_LpScatPlane.mag());
	    double theta2 = acos(cost1)*Rad2Deg;
	    std::cout << "theta2 = " << theta2 << std::endl;
	    Xdir.rotate(-theta2*Deg2Rad, Zdir);
	    Ydir.rotate(-theta2*Deg2Rad, Zdir);
	    std::cout << "Xdir = " << Xdir << std::endl;
	    std::cout << "Ydir = " << Ydir << std::endl;
	    std::cout << "Zdir = " << Zdir << std::endl;
	    std::cout << "cos(Ydir*NormY) = " << Ydir*NormY_LpScatPlane << std::endl;
	    if (std::abs(Ydir*NormY_LpScatPlane - 1)>0.02) {
	      ThreeVector Xdir2(1, 0, 0);
	      ThreeVector Ydir2(0, 1, 0);
	      Xdir2.rotateY(theta1*Deg2Rad);
	      Xdir2.rotateZ(phi1*Deg2Rad);
	      Ydir2.rotateY(theta1*Deg2Rad);
	      Ydir2.rotateZ(phi1*Deg2Rad);
	      Xdir2.rotate(theta2*Deg2Rad, Zdir);
	      Ydir2.rotate(theta2*Deg2Rad, Zdir);
	      std::cout << "Xdir2 = " << Xdir2 << std::endl;
	      std::cout << "Ydir2 = " << Ydir2 << std::endl;
	      std::cout << "cos(Ydir2*NormY) = " << Ydir2*NormY_LpScatPlane << std::endl;
	    }
	    */

	    double costDecayProtonLambdaPScat =
	      Mom2*NormY_LpScatPlane/(Mom2.mag()*NormY_LpScatPlane.mag());

	    double E = std::sqrt( LambdaMass*LambdaMass+momVecLambda.mag2());
	    ThreeVector BetaScatteredLambda = (momVecLambda/momVecLambda.mag())*(momVecLambda.mag()/E);
	    LorentzVector LvDecayProtonVec2( Mom2, std::sqrt( ProtonMass*ProtonMass+Mom2.mag2()));
	    LorentzVector LvDecayProtonVecCM2 = LvDecayProtonVec2.boost(-BetaScatteredLambda);
	    double costDecayProtonLambdaPScatCM  =
	      LvDecayProtonVecCM2.vect()*NormY_LpScatPlane/(LvDecayProtonVecCM2.vect().mag()*NormY_LpScatPlane.mag());

            double momCalLambda1 =-999.;
            double momCalLambda2 =-999.;
            double thetaCMLambdaPScat_Lambda =-999.;

            flagLambdaPScat_Lambda =
              calc2BodyInelastic(LambdaMass, LambdaMissMom.mag(),
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
            event.phiLambdaPScat_Lambda = phi_norm1;
            event.momCalLambda1 = momCalLambda1;
            event.momCalLambda2 = momCalLambda2;

	    event.costDecayProtonLambdaPScat = costDecayProtonLambdaPScat;
	    event.costDecayProtonLambdaBeam = costDecayProtonLambdaBeam;
	    event.costDecayProtonLambdaPScatCM = costDecayProtonLambdaPScatCM;
	    event.costDecayProtonLambdaBeamCM = costDecayProtonLambdaBeamCM;

            if (flagLambdaPScat_Lambda) {
              double dp = momCalLambda1 - momVecLambda.mag();
              //std::cout << "Lambda p scat : dp = " << dp << std::endl;
            }
          }
	  // recoil proton assumption
	  {
	    double cdistLambdaPScat_p;
	    ThreeVector VertLambdaPScat_p =
	      VertexPoint3D( VertexK, Pos2,
			     LambdaMissMom, Dir2, cdistLambdaPScat_p);

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
	    event.EkinP    = totE_proton;
	    event.EkinCorP = EkinCorP1;

	    double costScat1=LambdaMissMom*Mom2/(LambdaMissMom.mag()*Mom2.mag());
	    double thetaScat1=acos(costScat1)*Rad2Deg;

	    //std::cout << "Lambda p scat (1) : theta = " << thetaScat1 << std::endl;

	    event.vertexScat[0] = VertLambdaPScat_p.x();
	    event.vertexScat[1] = VertLambdaPScat_p.y();
	    event.vertexScat[2] = VertLambdaPScat_p.z();
	    event.cdistScat     = cdistLambdaPScat_p;

	    double scatMomCal1 = -999.;
	    double scatEkinCal1 = -999.;
	    double thetaScatCM1 = -999.;

	    if (thetaScat1>0. && thetaScat1<90.) {

	      bool flagLambdaPScat_p =
		calc2BodyKinema(LambdaMass, LambdaMissMom.mag(),
				ProtonMass, thetaScat1,
				&scatMomCal1, &scatEkinCal1, &thetaScatCM1);


	      if (flagLambdaPScat_p) {
		scatEkinCal1 *= 1000;
		event.scatMomCal  = scatMomCal1;
		event.scatEkinCal = scatEkinCal1;
		event.thetaScat   = thetaScat1;
		event.thetaScatCM = thetaScatCM1;

		double dE_LambdaPScat = EkinCorP1 - scatEkinCal1;
		//std::cout << "Lambda p scat (1) : dE = " << dE_LambdaPScat  << std::endl;

		ThreeVector ScatProtonMom = Mom2;
		ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

		/*
		std::cout << "Lambda Mom Vec ( " << LambdaMissMom.x()
			  <<  ", " << LambdaMissMom.y()
			  <<  ", " << LambdaMissMom.z()  << " )"  << std::endl;

		std::cout << "Proton Mom Vec ( " << ScatProtonMom.x()
			  <<  ", " << ScatProtonMom.y()
			  <<  ", " << ScatProtonMom.z()  << " )"  << std::endl;

		std::cout << "Recoil Sigma Mom Vec ( " << ScatHypMom.x()
			  <<  ", " << ScatHypMom.y()
			  <<  ", " << ScatHypMom.z()  << " )"  << std::endl;
		*/

		double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
		double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;

		event.costScatOpenAngle  = costScatOpenAngle;
		event.thetaScatOpenAngle = thetaScatOpenAngle;

		double cdistDecayPi;
		ThreeVector vertexDecayPi
		  = VertexPoint3D( VertLambdaPScat_p, Pos1,
				   ScatHypMom, Dir1,
				   cdistDecayPi );

		event.vertexDecayPi2[0] = vertexDecayPi.x();
		event.vertexDecayPi2[1] = vertexDecayPi.y();
		event.vertexDecayPi2[2] = vertexDecayPi.z();
		event.cdistDecayPi2     = cdistDecayPi;

		double cost = ScatHypMom*Dir1/(ScatHypMom.mag()*Dir1.mag());
		double decayPiMomCal;

		bool flagDecayCal = calcDecayPiMom( ScatHypMom.mag(),
						    LambdaMass,
						    ProtonMass,
						    PionMass,
						    cost,
						    &decayPiMomCal);

		if (flagDecayCal) {
		  ThreeVector PiDecayVecMom = Dir1 * (decayPiMomCal/Dir1.mag());
		  ThreeVector DecayProtonVecMom = ScatHypMom-PiDecayVecMom;

		  // Proton direction in Lp scattering plane
		  ThreeVector NormY_LpScatPlane = gaiseki(NormZ, ScatHypMom);
		  NormY_LpScatPlane *= (1./NormY_LpScatPlane.mag());
		  double costDecayProtonLambdaPScat_cal =
		    DecayProtonVecMom*NormY_LpScatPlane/(DecayProtonVecMom.mag()*NormY_LpScatPlane.mag());
		  event.costDecayProtonLambdaPScat_cal = costDecayProtonLambdaPScat_cal;

		  double E = std::sqrt( LambdaMass*LambdaMass+ScatHypMom.mag2());
		  ThreeVector BetaScatteredLambda = ScatHypMom*(1./E);
		  LorentzVector LvDecayProtonVec( DecayProtonVecMom, std::sqrt( ProtonMass*ProtonMass+DecayProtonVecMom.mag2()));
		  LorentzVector LvDecayProtonVecCM2 = LvDecayProtonVec.boost(-BetaScatteredLambda);
		  event.costDecayProtonLambdaPScatCM_cal  =
		    LvDecayProtonVecCM2.vect()*NormY_LpScatPlane/(LvDecayProtonVecCM2.vect().mag()*NormY_LpScatPlane.mag());
		}

		/*
		std::cout << "Decay Point of Scattered Lambda (1) ( " << vertexDecayPi.x()
			  << ", " << vertexDecayPi.y()
			  << ", " << vertexDecayPi.z()
			  << " ), cdist = " << cdistDecayPi << std::endl;
		*/
	      }


	      double momCal1 = -999.;
	      double momCal2 = -999.;
	      double thetaS0pConvCM1 = -999.;

	      bool flagSigma0p_p =
		calc2BodyInelastic(LambdaMass, LambdaMissMom.mag(),ProtonMass,
				   ProtonMass, SigmaZeroMass,
				   thetaScat1,
				   &momCal1, &momCal2, &thetaS0pConvCM1);

	      if (flagSigma0p_p) {
		thetaS0pConvCM1 = 180. - thetaS0pConvCM1;
		double scatEkinCal_S0pConv1 = (sqrt(momCal1*momCal1+ProtonMass*ProtonMass) - ProtonMass)*1000;
		double scatEkinCal_S0pConv2 = (sqrt(momCal2*momCal2+ProtonMass*ProtonMass) - ProtonMass)*1000;

		event.scatMomCal_S0pConv1 = momCal1;
		event.scatMomCal_S0pConv2 = momCal2;
		event.scatEkinCal_S0pConv1 = scatEkinCal_S0pConv1;
		event.scatEkinCal_S0pConv2 = scatEkinCal_S0pConv2;
		event.thetaS0pConvCM = thetaS0pConvCM1;

		double dE_Sigma0p_conv = EkinCorP1 - scatEkinCal_S0pConv1;
		//std::cout << "Lambda p scat (1) : dE = " << dE_LambdaPScat  << std::endl;

		ThreeVector ScatProtonMom = Mom2;
		ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

		/*
		std::cout << "Lambda Mom Vec ( " << LambdaMissMom.x()
			  <<  ", " << LambdaMissMom.y()
			  <<  ", " << LambdaMissMom.z()  << " )"  << std::endl;

		std::cout << "Proton Mom Vec ( " << ScatProtonMom.x()
			  <<  ", " << ScatProtonMom.y()
			  <<  ", " << ScatProtonMom.z()  << " )"  << std::endl;

		std::cout << "Recoil Sigma Mom Vec ( " << ScatHypMom.x()
			  <<  ", " << ScatHypMom.y()
			  <<  ", " << ScatHypMom.z()  << " )"  << std::endl;
		*/

		double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
		double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;

		event.costScatOpenAngle_S0pConv = costScatOpenAngle;
		event.thetaScatOpenAngle_S0pConv = thetaScatOpenAngle;

		double cdistDecayPi;
		ThreeVector vertexDecayPi
		  = VertexPoint3D( VertLambdaPScat_p, Pos1,
				   ScatHypMom, Dir1,
				   cdistDecayPi );

		event.vertexDecayPi2_S0pConv[0] = vertexDecayPi.x();
		event.vertexDecayPi2_S0pConv[1] = vertexDecayPi.y();
		event.vertexDecayPi2_S0pConv[2] = vertexDecayPi.z();
		event.cdistDecayPi2_S0pConv     = cdistDecayPi;

	      }

	    }
	  }

        }
      }
    } else if (nPi_CFT ==1 && nP_CFT ==2) {
      for (int i1 = 0; i1 < nPi_CFT; i1++ ) {

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

	  double cdistBeanAndVertex = -999.;
	  ThreeVector ppi = bMom;
	  ThreeVector xpi = bPos;

	  cdistBeanAndVertex =
	    closeDistLineAndPoint(xpi, ppi, vertex);

          event.energyBGO_p_lambda = ProtonCFT2->GetBGO_E();
          event.energyBGO_pi_lambda = PiCFT1->GetBGO_E();
          //event.segPiID_pi_lambda = PiCFT1->GetTrackPiIDSeg();

          double cost = Dir1*Dir2/(Dir1.mag()*Dir2.mag());
          double theta = std::acos(cost)*Rad2Deg;
	  event.theta_PiP = theta;

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
	  if (i2==0)
	    event.invMass_PiP = InvMass;
	  else if (i2==1)
	    event.invMass_PiP2 = InvMass;

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
          event.momLambda1 = (Mom1+Mom2).mag();
          event.momLambda2 = momVecLambda.mag();
          event.vertexLambdaDecay[0] = vertex.x();
          event.vertexLambdaDecay[1] = vertex.y();
          event.vertexLambdaDecay[2] = vertex.z();
          event.cdistLambdaDecay     = cdist;

	  if (i2==0) {
	    event.momLambda1_1st = (Mom1+Mom2).mag();
	    event.momLambda2_1st = momVecLambda.mag();
	    event.vertexLambdaDecay_1st[0] = vertex.x();
	    event.vertexLambdaDecay_1st[1] = vertex.y();
	    event.vertexLambdaDecay_1st[2] = vertex.z();
	    event.cdistLambdaDecay_1st     = cdist;
	  } else if (i2==1) {
	    event.momLambda1_2nd = (Mom1+Mom2).mag();
	    event.momLambda2_2nd = momVecLambda.mag();
	    event.vertexLambdaDecay_2nd[0] = vertex.x();
	    event.vertexLambdaDecay_2nd[1] = vertex.y();
	    event.vertexLambdaDecay_2nd[2] = vertex.z();
	    event.cdistLambdaDecay_2nd     = cdist;
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
              VertexPoint3D( VertexK, vertex,
			     LambdaMissMom, momVecLambda, cdistLambdaPScat_Lambda);

            double costLambdaPScat_Lambda=LambdaMissMom*momVecLambda/
              (LambdaMissMom.mag()*momVecLambda.mag());
            double thetaLambdaPScat_Lambda=acos(costLambdaPScat_Lambda)*Rad2Deg;

	    // phi angle for polarization
	    double p_NormX1 = momVecLambda*NormX;
	    double p_NormY1 = momVecLambda*NormY;
	    double p_NormZ1 = momVecLambda*NormZ;
	    ThreeVector pvec_Norm1(p_NormX1, p_NormY1, p_NormZ1);
	    double theta_norm1, phi_norm1;
	    calcThetaPhi(pvec_Norm1, &theta_norm1, &phi_norm1);
	    if (phi_norm1>270)
	      phi_norm1 -= 360.;

	    // Proton direction in Lp scattering plane
	    ThreeVector NormY_LpScatPlane = gaiseki(NormZ, momVecLambda);
	    NormY_LpScatPlane *= (1./NormY_LpScatPlane.mag());
	    double costDecayProtonLambdaPScat =
	      Mom2*NormY_LpScatPlane/(Mom2.mag()*NormY_LpScatPlane.mag());

	    double E = std::sqrt( LambdaMass*LambdaMass+momVecLambda.mag2());
	    ThreeVector BetaScatteredLambda = (momVecLambda/momVecLambda.mag())*(momVecLambda.mag()/E);
	    LorentzVector LvDecayProtonVec( Mom2, std::sqrt( ProtonMass*ProtonMass+Mom2.mag2()));
	    LorentzVector LvDecayProtonVecCM2 = LvDecayProtonVec.boost(-BetaScatteredLambda);
	    double costDecayProtonLambdaPScatCM  =
	      LvDecayProtonVecCM2.vect()*NormY_LpScatPlane/(LvDecayProtonVecCM2.vect().mag()*NormY_LpScatPlane.mag());

            double momCalLambda1 =-999.;
            double momCalLambda2 =-999.;
            double thetaCMLambdaPScat_Lambda =-999.;

            flagLambdaPScat_Lambda =
              calc2BodyInelastic(LambdaMass, LambdaMissMom.mag(),
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
            event.phiLambdaPScat_Lambda = phi_norm1;
            event.thetaCMLambdaPScat_Lambda = thetaCMLambdaPScat_Lambda;
            event.momCalLambda1 = momCalLambda1;
            event.momCalLambda2 = momCalLambda2;

	    event.costDecayProtonLambdaPScat = costDecayProtonLambdaPScat;
	    event.costDecayProtonLambdaPScatCM = costDecayProtonLambdaPScatCM;

	    if (i2==0) {
	      event.cdistLambdaPScat_Lambda_1st = cdistLambdaPScat_Lambda;
	      event.vertexLambdaPScat_Lambda_1st[0] = VertLambdaPScat_Lambda.x();
	      event.vertexLambdaPScat_Lambda_1st[1] = VertLambdaPScat_Lambda.y();
	      event.vertexLambdaPScat_Lambda_1st[2] = VertLambdaPScat_Lambda.z();

	      event.costLambdaPScat_Lambda_1st = costLambdaPScat_Lambda;
	      event.thetaLambdaPScat_Lambda_1st = thetaLambdaPScat_Lambda;
	      event.phiLambdaPScat_Lambda_1st = phi_norm1;
	      event.thetaCMLambdaPScat_Lambda_1st = thetaCMLambdaPScat_Lambda;
	      event.momCalLambda1_1st = momCalLambda1;
	      event.momCalLambda2_1st = momCalLambda2;

	      event.costDecayProtonLambdaPScat_1st = costDecayProtonLambdaPScat;
	      event.costDecayProtonLambdaPScatCM_1st = costDecayProtonLambdaPScatCM;
	    } else if (i2==1) {
	      event.cdistLambdaPScat_Lambda_2nd = cdistLambdaPScat_Lambda;
	      event.vertexLambdaPScat_Lambda_2nd[0] = VertLambdaPScat_Lambda.x();
	      event.vertexLambdaPScat_Lambda_2nd[1] = VertLambdaPScat_Lambda.y();
	      event.vertexLambdaPScat_Lambda_2nd[2] = VertLambdaPScat_Lambda.z();

	      event.costLambdaPScat_Lambda_2nd = costLambdaPScat_Lambda;
	      event.thetaLambdaPScat_Lambda_2nd = thetaLambdaPScat_Lambda;
	      event.phiLambdaPScat_Lambda_2nd = phi_norm1;
	      event.thetaCMLambdaPScat_Lambda_2nd = thetaCMLambdaPScat_Lambda;
	      event.momCalLambda1_2nd = momCalLambda1;
	      event.momCalLambda2_2nd = momCalLambda2;

	      event.costDecayProtonLambdaPScat_2nd = costDecayProtonLambdaPScat;
	      event.costDecayProtonLambdaPScatCM_2nd = costDecayProtonLambdaPScatCM;
	    }

            if (flagLambdaPScat_Lambda) {
              double dp = momCalLambda1 - momVecLambda.mag();
              //std::cout << "Lambda p scat : dp = " << dp << std::endl;
            }


	    // recoil proton assumption
	    ThreeVector momVecCalRecoilProton = LambdaMissMom-momVecLambda;

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
	      event.momCalRecoilProton_1st = momVecCalRecoilProton.mag();
	      event.momRecoilProton_1st = MomRP.mag();
	      event.deltaThetaRecoilProton_1st = theta_diff;
	    } else if (i2==1) {
	      event.momCalRecoilProton_2nd = momVecCalRecoilProton.mag();
	      event.momRecoilProton_2nd = MomRP.mag();
	      event.deltaThetaRecoilProton_2nd = theta_diff;
	    }
          }
	}
      }


      if (event.momCalLambda1_1st>=0 && event.momCalLambda1_2nd<0) {
	event.momLambda1           = event.momLambda1_1st;
	event.momLambda2           = event.momLambda2_1st;
	event.vertexLambdaDecay[0] = event.vertexLambdaDecay_1st[0];
	event.vertexLambdaDecay[1] = event.vertexLambdaDecay_1st[1];
	event.vertexLambdaDecay[2] = event.vertexLambdaDecay_1st[2];
	event.cdistLambdaDecay     = event.cdistLambdaDecay_1st;

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_1st;
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_1st[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_1st[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_1st[2];

	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_1st;
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_1st;
	event.phiLambdaPScat_Lambda       = event.phiLambdaPScat_Lambda_1st;
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_1st;
	event.momCalLambda1               = event.momCalLambda1_1st;
	event.momCalLambda2               = event.momCalLambda2_1st;

	event.momCalRecoilProton     = event.momCalRecoilProton_1st;
	event.momRecoilProton        = event.momRecoilProton_1st;
	event.deltaThetaRecoilProton = event.deltaThetaRecoilProton_1st;

	event.costDecayProtonLambdaPScat = event.costDecayProtonLambdaPScat_1st;
	event.costDecayProtonLambdaPScatCM = event.costDecayProtonLambdaPScatCM_1st;
      } else if (event.momCalLambda1_1st<0 && event.momCalLambda1_2nd>=0) {
	event.momLambda1           = event.momLambda1_2nd;
	event.momLambda2           = event.momLambda2_2nd;
	event.vertexLambdaDecay[0] = event.vertexLambdaDecay_2nd[0];
	event.vertexLambdaDecay[1] = event.vertexLambdaDecay_2nd[1];
	event.vertexLambdaDecay[2] = event.vertexLambdaDecay_2nd[2];
	event.cdistLambdaDecay     = event.cdistLambdaDecay_2nd;

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_2nd;
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_2nd[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_2nd[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_2nd[2];

	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_2nd;
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_2nd;
	event.phiLambdaPScat_Lambda       = event.phiLambdaPScat_Lambda_2nd;
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_2nd;
	event.momCalLambda1               = event.momCalLambda1_2nd;
	event.momCalLambda2               = event.momCalLambda2_2nd;

	event.momCalRecoilProton     = event.momCalRecoilProton_2nd;
	event.momRecoilProton        = event.momRecoilProton_2nd;
	event.deltaThetaRecoilProton = event.deltaThetaRecoilProton_2nd;

	event.costDecayProtonLambdaPScat = event.costDecayProtonLambdaPScat_2nd;
	event.costDecayProtonLambdaPScatCM = event.costDecayProtonLambdaPScatCM_2nd;

	double val1 = event.invMass_PiP;
	double val2 = event.invMass_PiP2;
	event.invMass_PiP = val2;
	event.invMass_PiP2 = val1;
      } else if (event.momCalLambda1_1st>=0 && event.momCalLambda1_2nd>=0
		 && std::abs(event.momCalLambda1_1st-event.momLambda2_1st) <=
		 std::abs(event.momCalLambda1_2nd-event.momLambda2_2nd)) {
	event.momLambda1           = event.momLambda1_1st;
	event.momLambda2           = event.momLambda2_1st;
	event.vertexLambdaDecay[0] = event.vertexLambdaDecay_1st[0];
	event.vertexLambdaDecay[1] = event.vertexLambdaDecay_1st[1];
	event.vertexLambdaDecay[2] = event.vertexLambdaDecay_1st[2];
	event.cdistLambdaDecay     = event.cdistLambdaDecay_1st;

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_1st;
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_1st[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_1st[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_1st[2];

	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_1st;
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_1st;
	event.phiLambdaPScat_Lambda       = event.phiLambdaPScat_Lambda_1st;
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_1st;
	event.momCalLambda1               = event.momCalLambda1_1st;
	event.momCalLambda2               = event.momCalLambda2_1st;

	event.momCalRecoilProton     = event.momCalRecoilProton_1st;
	event.momRecoilProton        = event.momRecoilProton_1st;
	event.deltaThetaRecoilProton = event.deltaThetaRecoilProton_1st;

	event.costDecayProtonLambdaPScat = event.costDecayProtonLambdaPScat_1st;
	event.costDecayProtonLambdaPScatCM = event.costDecayProtonLambdaPScatCM_1st;
      }  else if (event.momCalLambda1_1st>=0 && event.momCalLambda1_2nd>=0
		 && std::abs(event.momCalLambda1_1st-event.momLambda2_1st) >=
		 std::abs(event.momCalLambda1_2nd-event.momLambda2_2nd)) {
	event.momLambda1           = event.momLambda1_2nd;
	event.momLambda2           = event.momLambda2_2nd;
	event.vertexLambdaDecay[0] = event.vertexLambdaDecay_2nd[0];
	event.vertexLambdaDecay[1] = event.vertexLambdaDecay_2nd[1];
	event.vertexLambdaDecay[2] = event.vertexLambdaDecay_2nd[2];
	event.cdistLambdaDecay     = event.cdistLambdaDecay_2nd;

	event.cdistLambdaPScat_Lambda     = event.cdistLambdaPScat_Lambda_2nd;
	event.vertexLambdaPScat_Lambda[0] = event.vertexLambdaPScat_Lambda_2nd[0];
	event.vertexLambdaPScat_Lambda[1] = event.vertexLambdaPScat_Lambda_2nd[1];
	event.vertexLambdaPScat_Lambda[2] = event.vertexLambdaPScat_Lambda_2nd[2];

	event.costLambdaPScat_Lambda      = event.costLambdaPScat_Lambda_2nd;
	event.thetaLambdaPScat_Lambda     = event.thetaLambdaPScat_Lambda_2nd;
	event.phiLambdaPScat_Lambda       = event.phiLambdaPScat_Lambda_2nd;
	event.thetaCMLambdaPScat_Lambda   = event.thetaCMLambdaPScat_Lambda_2nd;
	event.momCalLambda1               = event.momCalLambda1_2nd;
	event.momCalLambda2               = event.momCalLambda2_2nd;

	event.momCalRecoilProton     = event.momCalRecoilProton_2nd;
	event.momRecoilProton        = event.momRecoilProton_2nd;
	event.deltaThetaRecoilProton = event.deltaThetaRecoilProton_2nd;

	event.costDecayProtonLambdaPScat = event.costDecayProtonLambdaPScat_2nd;
	event.costDecayProtonLambdaPScatCM = event.costDecayProtonLambdaPScatCM_2nd;

	double val1 = event.invMass_PiP;
	double val2 = event.invMass_PiP2;
	event.invMass_PiP = val2;
	event.invMass_PiP2 = val1;
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

      ThreeVector vtx = VertexK;

      ThreeVector VertScat1 =
        VertexPoint3D( vtx, ProtonVec1.pos0,
		       LambdaMissMom, ProtonVec1.mom, cdist1);

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
		       LambdaMissMom, ProtonVec2.mom, cdist2);

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
      double costScat1=LambdaMissMom*ProtonVec1.mom/(LambdaMissMom.mag()*ProtonVec1.mom.mag());
      double thetaScat1=acos(costScat1)*Rad2Deg;

      // phi angle for polarization
      double p_NormX1 = ProtonVec1.mom*NormX;
      double p_NormY1 = ProtonVec1.mom*NormY;
      double p_NormZ1 = ProtonVec1.mom*NormZ;
      ThreeVector pvec_Norm1(p_NormX1, p_NormY1, p_NormZ1);
      double theta_norm1, phi_norm1;
      calcThetaPhi(pvec_Norm1, &theta_norm1, &phi_norm1);
      if (phi_norm1>270)
	phi_norm1 -= 360.;

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
      double costDecayProtonLambdaPScat1 = -999.;
      double costDecayProtonLambdaPScatCM1 = -999.;


      double scatMomCal_S0pConv1_1;
      double scatMomCal_S0pConv2_1;
      double scatEkinCal_S0pConv1_1;
      double scatEkinCal_S0pConv2_1;
      double thetaS0pConvCM_1;

      double vertexDecayProton_S0pConv1_2p[3] = {-999., -999., -999.};
      double cdistDecayProton_S0pConv1_2p = -999.;
      double MissMassPiGamma_1_2p = -999.;

      if (thetaScat1>0. && thetaScat1<90.) {

        flagLambdaPScat =
          calc2BodyKinema(LambdaMass, LambdaMissMom.mag(),
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
          ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

          /*
	    std::cout << "Lambda Mom Vec ( " << LambdaMissMom.x()
	    <<  ", " << LambdaMissMom.y()
	    <<  ", " << LambdaMissMom.z()  << " )"  << std::endl;

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

	  // Proton direction in Lp scattering plane
	  ThreeVector NormY_LpScatPlane = gaiseki(NormZ, ScatHypMom);
	  NormY_LpScatPlane *= (1./NormY_LpScatPlane.mag());
	  costDecayProtonLambdaPScat1 =
	    ProtonVec2.mom*NormY_LpScatPlane/(ProtonVec2.mom.mag()*NormY_LpScatPlane.mag());

	  double E = std::sqrt( LambdaMass*LambdaMass+ScatHypMom.mag2());
	  ThreeVector BetaScatteredLambda = ScatHypMom*(1./E);
	  LorentzVector LvDecayProtonVec( ProtonVec2.mom, std::sqrt( ProtonMass*ProtonMass+ProtonVec2.mom.mag2()));
	  LorentzVector LvDecayProtonVecCM2 = LvDecayProtonVec.boost(-BetaScatteredLambda);
	  costDecayProtonLambdaPScatCM1 =
	    LvDecayProtonVecCM2.vect()*NormY_LpScatPlane/(LvDecayProtonVecCM2.vect().mag()*NormY_LpScatPlane.mag());


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

	double momCal1 = -999.;
	double momCal2 = -999.;
	double thetaS0pConvCM1 = -999.;

	bool flagSigma0p_p =
	  calc2BodyInelastic(LambdaMass, LambdaMissMom.mag(),ProtonMass,
			     ProtonMass, SigmaZeroMass,
			     thetaScat1,
			     &momCal1, &momCal2, &thetaS0pConvCM1);

        if (flagSigma0p_p) {
	  thetaS0pConvCM1 = 180. - thetaS0pConvCM1;
	  double scatEkinCal_S0pConv1 = (sqrt(momCal1*momCal1+ProtonMass*ProtonMass) - ProtonMass)*1000;
	  double scatEkinCal_S0pConv2 = (sqrt(momCal2*momCal2+ProtonMass*ProtonMass) - ProtonMass)*1000;

	  scatMomCal_S0pConv1_1 = momCal1;
	  scatMomCal_S0pConv2_1  = momCal2;
	  scatEkinCal_S0pConv1_1  = scatEkinCal_S0pConv1;
	  scatEkinCal_S0pConv2_1  = scatEkinCal_S0pConv2;
	  thetaS0pConvCM_1  = thetaS0pConvCM1;

          ThreeVector ScatProtonMom = ProtonVec1.mom;
          ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

          double cdistDecayProton1;
          ThreeVector vertexDecayProton1
            = VertexPoint3D( VertScat1, ProtonVec2.pos0,
			     ScatHypMom, ProtonVec2.mom,
			     cdistDecayProton1 );

          vertexDecayProton_S0pConv1_2p[0] = vertexDecayProton1.x();
          vertexDecayProton_S0pConv1_2p[1] = vertexDecayProton1.y();
          vertexDecayProton_S0pConv1_2p[2] = vertexDecayProton1.z();
          cdistDecayProton_S0pConv1_2p = cdistDecayProton1;

          LorentzVector LvProton(ProtonVec2.mom ,
				 sqrt(ProtonMass*ProtonMass+ProtonVec2.mom.mag2()) );
          LorentzVector LvBeam(ScatHypMom,
			       sqrt(SigmaZeroMass*SigmaZeroMass+ScatHypMom.mag2()) );
          LorentzVector LvPi = LvBeam-LvProton;
          double MissMassPi0 = LvPi.mag2();

          MissMassPiGamma_1_2p = MissMassPi0;

        }

      }

      /* Analysis of 2nd proton*/

      double costScat2=LambdaMissMom*ProtonVec2.mom/(LambdaMissMom.mag()*ProtonVec2.mom.mag());
      double thetaScat2=acos(costScat2)*Rad2Deg;

      // phi angle for polarization
      double p_NormX2 = ProtonVec2.mom*NormX;
      double p_NormY2 = ProtonVec2.mom*NormY;
      double p_NormZ2 = ProtonVec2.mom*NormZ;
      ThreeVector pvec_Norm2(p_NormX2, p_NormY2, p_NormZ2);
      double theta_norm2, phi_norm2;
      calcThetaPhi(pvec_Norm2, &theta_norm2, &phi_norm2);
      if (phi_norm2>270)
	phi_norm2 -= 360.;

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

      double costDecayProtonLambdaPScat2 = -999.;
      double costDecayProtonLambdaPScatCM2 = -999.;

      double scatMomCal_S0pConv1_2;
      double scatMomCal_S0pConv2_2;
      double scatEkinCal_S0pConv1_2;
      double scatEkinCal_S0pConv2_2;
      double thetaS0pConvCM_2;

      double vertexDecayProton_S0pConv2_2p[3] = {-999., -999., -999.};
      double cdistDecayProton_S0pConv2_2p = -999.;
      double MissMassPiGamma_2_2p = -999.;

      if (thetaScat2>0. && thetaScat2<90.) {

        flagLambdaPScat =
          calc2BodyKinema(LambdaMass, LambdaMissMom.mag(),
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
          ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;
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


	  // Proton direction in Lp scattering plane
	  ThreeVector NormY_LpScatPlane = gaiseki(NormZ, ScatHypMom);
	  NormY_LpScatPlane *= (1./NormY_LpScatPlane.mag());
	  costDecayProtonLambdaPScat2 =
	    ProtonVec1.mom*NormY_LpScatPlane/(ProtonVec1.mom.mag()*NormY_LpScatPlane.mag());

	  double E = std::sqrt( LambdaMass*LambdaMass+ScatHypMom.mag2());
	  ThreeVector BetaScatteredLambda = ScatHypMom*(1./E);
	  LorentzVector LvDecayProtonVec( ProtonVec1.mom, std::sqrt( ProtonMass*ProtonMass+ProtonVec1.mom.mag2()));
	  LorentzVector LvDecayProtonVecCM2 = LvDecayProtonVec.boost(-BetaScatteredLambda);
	  costDecayProtonLambdaPScatCM2 =
	    LvDecayProtonVecCM2.vect()*NormY_LpScatPlane/(LvDecayProtonVecCM2.vect().mag()*NormY_LpScatPlane.mag());


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

	double momCal1 = -999.;
	double momCal2 = -999.;
	double thetaS0pConvCM1 = -999.;

	bool flagSigma0p_p =
	  calc2BodyInelastic(LambdaMass, LambdaMissMom.mag(),ProtonMass,
			     ProtonMass, SigmaZeroMass,
			     thetaScat2,
			     &momCal1, &momCal2, &thetaS0pConvCM1);

        if (flagSigma0p_p) {
	  thetaS0pConvCM1 = 180. - thetaS0pConvCM1;
	  double scatEkinCal_S0pConv1 = (sqrt(momCal1*momCal1+ProtonMass*ProtonMass) - ProtonMass)*1000;
	  double scatEkinCal_S0pConv2 = (sqrt(momCal2*momCal2+ProtonMass*ProtonMass) - ProtonMass)*1000;

	  scatMomCal_S0pConv1_2 = momCal1;
	  scatMomCal_S0pConv2_2  = momCal2;
	  scatEkinCal_S0pConv1_2  = scatEkinCal_S0pConv1;
	  scatEkinCal_S0pConv2_2  = scatEkinCal_S0pConv2;
	  thetaS0pConvCM_2  = thetaS0pConvCM1;

          ThreeVector ScatProtonMom = ProtonVec2.mom;
          ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

          double cdistDecayProton1;
          ThreeVector vertexDecayProton1
            = VertexPoint3D( VertScat2, ProtonVec1.pos0,
			     ScatHypMom, ProtonVec1.mom,
			     cdistDecayProton1 );

          vertexDecayProton_S0pConv2_2p[0] = vertexDecayProton1.x();
          vertexDecayProton_S0pConv2_2p[1] = vertexDecayProton1.y();
          vertexDecayProton_S0pConv2_2p[2] = vertexDecayProton1.z();
          cdistDecayProton_S0pConv2_2p = cdistDecayProton1;

          LorentzVector LvProton(ProtonVec1.mom ,
				 sqrt(ProtonMass*ProtonMass+ProtonVec1.mom.mag2()) );
          LorentzVector LvBeam(ScatHypMom,
			       sqrt(SigmaZeroMass*SigmaZeroMass+ScatHypMom.mag2()) );
          LorentzVector LvPi = LvBeam-LvProton;
          double MissMassPi0 = LvPi.mag2();

          MissMassPiGamma_2_2p = MissMassPi0;

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

        event.EkinP = ProtonCFT1->GetTotalE();
        event.vertexScat[0] = VertScat1.x();
        event.vertexScat[1] = VertScat1.y();
        event.vertexScat[2] = VertScat1.z();
        event.cdistScat = cdist1;
        event.thetaScat = thetaScat1;
        event.phiScat   = phi_norm1;
        event.scatMomCal = scatMomCal1;
        event.scatEkinCal = scatEkinCal1;
        event.thetaScatCM = thetaScatCM1;
        event.EkinCorP = EkinCorP1;
        event.vertexDecayProton1_2p[0] = vertexDecayProton1_2p[0];
        event.vertexDecayProton1_2p[1] = vertexDecayProton1_2p[1];
        event.vertexDecayProton1_2p[2] = vertexDecayProton1_2p[2];
        event.cdistDecayProton1_2p     = cdistDecayProton1_2p;
        event.MissMassPi_1_2p         = MissMassPi_1_2p;
        event.costDecayProtonLambdaPScat1_2p = costDecayProtonLambdaPScat1;
        event.costDecayProtonLambdaPScatCM1_2p = costDecayProtonLambdaPScatCM1;

        event.EkinP_2nd = ProtonCFT2->GetTotalE();
        event.vertexScat_2nd[0] = VertScat2.x();
        event.vertexScat_2nd[1] = VertScat2.y();
        event.vertexScat_2nd[2] = VertScat2.z();
        event.cdistScat_2nd = cdist2;
        event.thetaScat_2nd = thetaScat2;
        event.phiScat_2nd   = phi_norm2;
        event.scatMomCal_2nd = scatMomCal2;
        event.scatEkinCal_2nd = scatEkinCal2;
        event.thetaScatCM_2nd = thetaScatCM2;
        event.EkinCorP_2nd = EkinCorP2;
        event.vertexDecayProton2_2p[0] = vertexDecayProton2_2p[0];
        event.vertexDecayProton2_2p[1] = vertexDecayProton2_2p[1];
        event.vertexDecayProton2_2p[2] = vertexDecayProton2_2p[2];
        event.cdistDecayProton2_2p     = cdistDecayProton2_2p;
        event.MissMassPi_2_2p         = MissMassPi_2_2p;
        event.costDecayProtonLambdaPScat2_2p = costDecayProtonLambdaPScat2;
        event.costDecayProtonLambdaPScatCM2_2p = costDecayProtonLambdaPScatCM2;
      } else {
        //std::cout << "dE (2nd) = " << dE_LambdaPScat2
        //<< ", dE (1st) = " << dE_LambdaPScat1
        //<< std::endl;
        event.EkinP = ProtonCFT2->GetTotalE();
        event.vertexScat[0] = VertScat2.x();
        event.vertexScat[1] = VertScat2.y();
        event.vertexScat[2] = VertScat2.z();
        event.cdistScat = cdist2;
        event.thetaScat = thetaScat2;
        event.phiScat   = phi_norm2;
        event.scatMomCal = scatMomCal2;
        event.scatEkinCal = scatEkinCal2;
        event.thetaScatCM = thetaScatCM2;
        event.EkinCorP = EkinCorP2;

        event.vertexDecayProton1_2p[0] = vertexDecayProton2_2p[0];
        event.vertexDecayProton1_2p[1] = vertexDecayProton2_2p[1];
        event.vertexDecayProton1_2p[2] = vertexDecayProton2_2p[2];
        event.cdistDecayProton1_2p     = cdistDecayProton2_2p;
        event.MissMassPi_1_2p         = MissMassPi_2_2p;
        event.costDecayProtonLambdaPScat1_2p = costDecayProtonLambdaPScat2;
        event.costDecayProtonLambdaPScatCM1_2p = costDecayProtonLambdaPScatCM2;

        event.EkinP_2nd = ProtonCFT1->GetTotalE();
        event.vertexScat_2nd[0] = VertScat1.x();
        event.vertexScat_2nd[1] = VertScat1.y();
        event.vertexScat_2nd[2] = VertScat1.z();
        event.cdistScat_2nd = cdist1;
        event.thetaScat_2nd = thetaScat1;
        event.phiScat_2nd   = phi_norm1;
        event.scatMomCal_2nd = scatMomCal1;
        event.scatEkinCal_2nd = scatEkinCal1;
        event.thetaScatCM_2nd = thetaScatCM1;
        event.EkinCorP_2nd = EkinCorP1;

        event.vertexDecayProton2_2p[0] = vertexDecayProton1_2p[0];
        event.vertexDecayProton2_2p[1] = vertexDecayProton1_2p[1];
        event.vertexDecayProton2_2p[2] = vertexDecayProton1_2p[2];
        event.cdistDecayProton2_2p     = cdistDecayProton1_2p;
        event.MissMassPi_2_2p         = MissMassPi_1_2p;
        event.costDecayProtonLambdaPScat2_2p = costDecayProtonLambdaPScat1;
        event.costDecayProtonLambdaPScatCM2_2p = costDecayProtonLambdaPScatCM1;
      }

      double dE_S0pConv1 = -999.;
      if (EkinCorP1>0 && scatEkinCal_S0pConv1_1>0)
        dE_S0pConv1 = EkinCorP1 - scatEkinCal_S0pConv1_1;

      double dE_S0pConv2 = -999.;
      if (EkinCorP2>0 && scatEkinCal_S0pConv1_2>0)
        dE_S0pConv2 = EkinCorP2 - scatEkinCal_S0pConv1_2;

      if (fabs(dE_S0pConv1) <= fabs(dE_S0pConv2)) {
	event.scatMomCal_S0pConv1 = scatMomCal_S0pConv1_1;
	event.scatMomCal_S0pConv2 = scatMomCal_S0pConv2_1;
	event.scatEkinCal_S0pConv1 = scatEkinCal_S0pConv1_1;
	event.scatEkinCal_S0pConv2 = scatEkinCal_S0pConv2_1;
	event.thetaS0pConvCM = thetaS0pConvCM_1;

	event.vertexDecayProton_S0pConv1_2p[0] = vertexDecayProton_S0pConv1_2p[0];
	event.vertexDecayProton_S0pConv1_2p[1] = vertexDecayProton_S0pConv1_2p[1];
	event.vertexDecayProton_S0pConv1_2p[2] = vertexDecayProton_S0pConv1_2p[2];
	event.cdistDecayProton_S0pConv1_2p = vertexDecayProton_S0pConv1_2p[2];
	event.MissMassPiGamma_1_2p = cdistDecayProton_S0pConv1_2p;

	event.scatMomCal_S0pConv1_2nd = scatMomCal_S0pConv1_2;
	event.scatMomCal_S0pConv2_2nd = scatMomCal_S0pConv2_2;
	event.scatEkinCal_S0pConv1_2nd = scatEkinCal_S0pConv1_2;
	event.scatEkinCal_S0pConv2_2nd = scatEkinCal_S0pConv2_2;
	event.thetaS0pConvCM_2nd = thetaS0pConvCM_2;

	event.vertexDecayProton_S0pConv1_2p_2nd[0] = vertexDecayProton_S0pConv2_2p[0];
	event.vertexDecayProton_S0pConv1_2p_2nd[1] = vertexDecayProton_S0pConv2_2p[1];
	event.vertexDecayProton_S0pConv1_2p_2nd[2] = vertexDecayProton_S0pConv2_2p[2];
	event.cdistDecayProton_S0pConv1_2p_2nd = vertexDecayProton_S0pConv2_2p[2];
	event.MissMassPiGamma_1_2p_2nd = cdistDecayProton_S0pConv2_2p;
      } else {
	event.scatMomCal_S0pConv1 = scatMomCal_S0pConv1_2;
	event.scatMomCal_S0pConv2 = scatMomCal_S0pConv2_2;
	event.scatEkinCal_S0pConv1 = scatEkinCal_S0pConv1_2;
	event.scatEkinCal_S0pConv2 = scatEkinCal_S0pConv2_2;
	event.thetaS0pConvCM = thetaS0pConvCM_2;

	event.vertexDecayProton_S0pConv1_2p[0] = vertexDecayProton_S0pConv2_2p[0];
	event.vertexDecayProton_S0pConv1_2p[1] = vertexDecayProton_S0pConv2_2p[1];
	event.vertexDecayProton_S0pConv1_2p[2] = vertexDecayProton_S0pConv2_2p[2];
	event.cdistDecayProton_S0pConv1_2p = vertexDecayProton_S0pConv2_2p[2];
	event.MissMassPiGamma_1_2p = cdistDecayProton_S0pConv2_2p;

	event.scatMomCal_S0pConv1_2nd = scatMomCal_S0pConv1_1;
	event.scatMomCal_S0pConv2_2nd = scatMomCal_S0pConv2_1;
	event.scatEkinCal_S0pConv1_2nd = scatEkinCal_S0pConv1_1;
	event.scatEkinCal_S0pConv2_2nd = scatEkinCal_S0pConv2_1;
	event.thetaS0pConvCM_2nd = thetaS0pConvCM_1;

	event.vertexDecayProton_S0pConv1_2p_2nd[0] = vertexDecayProton_S0pConv1_2p[0];
	event.vertexDecayProton_S0pConv1_2p_2nd[1] = vertexDecayProton_S0pConv1_2p[1];
	event.vertexDecayProton_S0pConv1_2p_2nd[2] = vertexDecayProton_S0pConv1_2p[2];
	event.cdistDecayProton_S0pConv1_2p_2nd = vertexDecayProton_S0pConv1_2p[2];
	event.MissMassPiGamma_1_2p_2nd = cdistDecayProton_S0pConv1_2p;

      }


      // pp scattering assumption
      double cdistTwoProton;
      ThreeVector VertPPScat = VertexPoint3D( ProtonVec1.pos0, ProtonVec2.pos0,
					      ProtonVec1.mom, ProtonVec2.mom, cdistTwoProton);

      double costTwoProton = ProtonVec1.mom*ProtonVec2.mom/(ProtonVec1.mom.mag()*ProtonVec2.mom.mag());
      double thetaTwoProton = acos(costTwoProton)*Rad2Deg;

      double p_cor1, e_cor1;
      double totE_p1 = ProtonCFT1->GetTotalE();
      double p_p1 = sqrt(totE_p1*totE_p1+2.*totE_p1*ProtonMass*1000.);
      p_p1 /= 1000.; // GeV/c

      CorrElossOutWithCFRP(&p_cor1, &e_cor1, p_p1, PROTON,
			   ProtonVec1.mom/ProtonVec1.mom.mag(), VertPPScat,
			   ProtonVec1.pos0);

      double p_cor2, e_cor2;
      double totE_p2 = ProtonCFT2->GetTotalE();
      double p_p2 = sqrt(totE_p2*totE_p2+2.*totE_p2*ProtonMass*1000.);
      p_p2 /= 1000.; // GeV/c

      CorrElossOutWithCFRP(&p_cor2, &e_cor2, p_p2, PROTON,
			   ProtonVec2.mom/ProtonVec2.mom.mag(), VertPPScat,
			   ProtonVec2.pos0);

      ProtonVec1.mom *= p_cor1/ProtonVec1.mom.mag();
      ProtonVec2.mom *= p_cor2/ProtonVec2.mom.mag();

      // primary proton
      ThreeVector momPrimP = ProtonVec1.mom + ProtonVec2.mom;
      double momp0 = momPrimP.mag();

      double cdistLambdaProtonDecay;
      ThreeVector VertLambdaProtonDecay = VertexPoint3D( VertexK, VertPPScat,
							 LambdaMissMom, momPrimP, cdistLambdaProtonDecay);
      double costLambdaProtonDecay=LambdaMissMom*momPrimP/(LambdaMissMom.mag()*momPrimP.mag());
      double thetaLambdaProtonDecay=acos(costLambdaProtonDecay)*Rad2Deg;
      double momCalLambdaProtonDecay1 = -999.9;
      double momCalLambdaProtonDecay2 = -999.9;

      {
	double p1 = LambdaMissMom.mag();
	double m1 = LambdaMass;
	double E1 = sqrt(p1*p1+m1*m1);
	double m2 = PionMass;
	double m3 = ProtonMass;
	double A  = (m1*m1+m3*m3-m2*m2)/2.;
	double hanbetu = (A*p1*costLambdaProtonDecay)*(A*p1*costLambdaProtonDecay)-(E1*E1-p1*p1*costLambdaProtonDecay*costLambdaProtonDecay)*(E1*E1*m3*m3-A*A);
	if (hanbetu >= 0) {
	  double ans1 = (A*p1*costLambdaProtonDecay+sqrt(hanbetu))/(E1*E1-p1*p1*costLambdaProtonDecay*costLambdaProtonDecay);
	  double ans2 = (A*p1*costLambdaProtonDecay-sqrt(hanbetu))/(E1*E1-p1*p1*costLambdaProtonDecay*costLambdaProtonDecay);
	  if (ans1>=0 || ans2>=0) {
	    if (fabs(momp0-ans1) <= fabs(momp0-ans2)) {
	      momCalLambdaProtonDecay1 = ans1;
	      momCalLambdaProtonDecay2 = ans2;
	    } else {
	      momCalLambdaProtonDecay2 = ans1;
	      momCalLambdaProtonDecay1 = ans2;
	    }
	  } else if (ans1<0 && ans2<0) {
	    std::cout << "reconstructed decayProtonMomCal two negative answers: "
		      << ans1 << ", " << ans2
		      << std::endl;
	  }
	}
      }

      event.vertexPPScat[0] = VertPPScat.x();
      event.vertexPPScat[1] = VertPPScat.y();
      event.vertexPPScat[2] = VertPPScat.z();
      event.cdistTwoProton = cdistTwoProton;
      event.thetaTwoProton = thetaTwoProton;
      event.momLambdaProtonDecay = momp0;
      event.momLambdaProtonDecayVec[0] = momPrimP.x();
      event.momLambdaProtonDecayVec[1] = momPrimP.y();
      event.momLambdaProtonDecayVec[2] = momPrimP.z();
      event.vertexLambdaProtonDecay[0] = VertLambdaProtonDecay.x();
      event.vertexLambdaProtonDecay[1] = VertLambdaProtonDecay.y();
      event.vertexLambdaProtonDecay[2] = VertLambdaProtonDecay.z();
      event.cdistLambdaProtonDecay = cdistLambdaProtonDecay;
      event.thetaLambdaProtonDecay = thetaLambdaProtonDecay;
      event.momCalLambdaProtonDecay1 = momCalLambdaProtonDecay1;
      event.momCalLambdaProtonDecay2 = momCalLambdaProtonDecay2;


    } else if ( nP_CFT ==1) {
      CFTParticle *ProtonCFT1 = CFTProtonCont[0];
      double EkinP1 = ProtonCFT1->GetTotalE();
      PartVec ProtonVec1;
      ProtonVec1.pos0 = ProtonCFT1->GetPos0();
      ProtonVec1.mom = ProtonCFT1->GetDir();

      /* Lambda p scattering assumption */
      double cdist1;

      ThreeVector vtx = VertexK;

      ThreeVector VertScat1 =
        VertexPoint3D( vtx, ProtonVec1.pos0,
		       LambdaMissMom, ProtonVec1.mom, cdist1);

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

      double costScat1=LambdaMissMom*ProtonVec1.mom/(LambdaMissMom.mag()*ProtonVec1.mom.mag());
      double thetaScat1=acos(costScat1)*Rad2Deg;

      // phi angle for polarization
      double p_NormX1 = ProtonVec1.mom*NormX;
      double p_NormY1 = ProtonVec1.mom*NormY;
      double p_NormZ1 = ProtonVec1.mom*NormZ;
      ThreeVector pvec_Norm1(p_NormX1, p_NormY1, p_NormZ1);
      double theta_norm1, phi_norm1;
      calcThetaPhi(pvec_Norm1, &theta_norm1, &phi_norm1);
      if (phi_norm1>270)
	phi_norm1 -= 360.;

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

        flagLambdaPScat =
          calc2BodyKinema(LambdaMass, LambdaMissMom.mag(),
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
          ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

          /*
	    std::cout << "Lambda Mom Vec ( " << LambdaMissMom.x()
	    <<  ", " << LambdaMissMom.y()
	    <<  ", " << LambdaMissMom.z()  << " )"  << std::endl;

	    std::cout << "Proton Mom Vec ( " << ScatProtonMom.x()
	    <<  ", " << ScatProtonMom.y()
	    <<  ", " << ScatProtonMom.z()  << " )"  << std::endl;

	    std::cout << "Recoil Sigma Mom Vec ( " << ScatHypMom.x()
	    <<  ", " << ScatHypMom.y()
	    <<  ", " << ScatHypMom.z()  << " )"  << std::endl;
	  */

          double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
          double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;
	}

	double momCal1 = -999.;
	double momCal2 = -999.;
	double thetaS0pConvCM1 = -999.;

	bool flagSigma0p_p =
	  calc2BodyInelastic(LambdaMass, LambdaMissMom.mag(),ProtonMass,
			     ProtonMass, SigmaZeroMass,
			     thetaScat1,
			     &momCal1, &momCal2, &thetaS0pConvCM1);

	if (flagSigma0p_p) {
	  thetaS0pConvCM1 = 180. - thetaS0pConvCM1;
	  double scatEkinCal_S0pConv1 = (sqrt(momCal1*momCal1+ProtonMass*ProtonMass) - ProtonMass)*1000;
	  double scatEkinCal_S0pConv2 = (sqrt(momCal2*momCal2+ProtonMass*ProtonMass) - ProtonMass)*1000;

	  event.scatMomCal_S0pConv1 = momCal1;
	  event.scatMomCal_S0pConv2 = momCal2;
	  event.scatEkinCal_S0pConv1 = scatEkinCal_S0pConv1;
	  event.scatEkinCal_S0pConv2 = scatEkinCal_S0pConv2;
	  event.thetaS0pConvCM = thetaS0pConvCM1;

	  ThreeVector ScatProtonMom = ProtonVec1.mom;
	  ThreeVector ScatHypMom = LambdaMissMom-ScatProtonMom;

	  double costScatOpenAngle=ScatHypMom*ScatProtonMom/(ScatHypMom.mag()*ScatProtonMom.mag());
	  double thetaScatOpenAngle=acos(costScatOpenAngle)*Rad2Deg;

	  event.costScatOpenAngle_S0pConv = costScatOpenAngle;
	  event.thetaScatOpenAngle_S0pConv = thetaScatOpenAngle;

	}
      }

      event.EkinP = ProtonCFT1->GetTotalE();
      event.vertexScat[0] = VertScat1.x();
      event.vertexScat[1] = VertScat1.y();
      event.vertexScat[2] = VertScat1.z();
      event.cdistScat = cdist1;
      event.thetaScat = thetaScat1;
      event.phiScat   = phi_norm1;
      event.scatMomCal = scatMomCal1;
      event.scatEkinCal = scatEkinCal1;
      event.thetaScatCM = thetaScatCM1;
      event.EkinCorP = EkinCorP1;

      /* pp scattering assumption */
      if (flagDecayCal) {

	double cdistPpScat;

	ThreeVector VertPpScat = VertexPoint3D( DecayProtonVec.pos0, ProtonVec1.pos0,
						DecayProtonVec.mom, ProtonVec1.mom,
						cdistPpScat);
	double costPpScat=DecayProtonVec.mom*ProtonVec1.mom/
	  (DecayProtonVec.mom.mag()*ProtonVec1.mom.mag());
	double thetaPpScat=acos(costPpScat)*Rad2Deg;

	event.vertexPpScat[0] = VertPpScat.x();
	event.vertexPpScat[1] = VertPpScat.y();
	event.vertexPpScat[2] = VertPpScat.z();
	event.cdistPpScat = cdistPpScat;
	event.thetaPpScat = thetaPpScat;

	if (thetaPpScat>0.&&thetaPpScat<90.) {
	  double scatPpMomCal;
	  double scatPpEkinCal;
	  double thetaPpScatCM;

	  bool flagPpScat = calc2BodyKinema(ProtonMass, DecayProtonVec.mom.mag(),
					    ProtonMass, thetaPpScat,
					    &scatPpMomCal, &scatPpEkinCal, &thetaPpScatCM);

	  if (flagPpScat) {
	    event.scatPpMomCal  = scatPpMomCal;
	    event.scatPpEkinCal = scatPpEkinCal*1000;
	    event.thetaPpScatCM = thetaPpScatCM;
	  }
	}
      }
      /* End of pp scattering assumption */


    }
  }


  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();

    double mom_pi = event.momPiMinusFromK0*1000;
    double mass_pi = 139.0;
    double Ekin_pi = sqrt(mass_pi*mass_pi+mom_pi*mom_pi)-mass_pi;
    std::cout << "Ekin_pi(from K0) = " << Ekin_pi << ", p = " << mom_pi << std::endl;

    if (1) {
    //if (index_pi_from_K0!=-1) {
    //if (index_pi_from_K0!=-1 && event.invMass_PiP>1) {
      //if (ntCFT_2nd>0) {
    //if (nPi_CFT>=2 && nP_CFT==1) {
    //if (event.scatFlag==-1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1&&nPi_CFT==1) {
    //if (event.scatFlag==-1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1) {
    //if (nP_CFT==1&&nPi_CFT==1) {
    //if (event.scatFlag==1) {
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

double calcKMomFromThetaLambda(double bmom, double cost, double p0)
{
  double E1 = sqrt(GammaMass*GammaMass + bmom*bmom);
  double A = GammaMass*GammaMass + ProtonMass*ProtonMass
    + KaonMass*KaonMass - LambdaMass*LambdaMass
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
  event.costDecayProtonCM0 = -999.9;

  event.egamma=0;
  event.ntSks=0;
  event.ntSksPart=0;
  event.nK=0;
  event.nLambda=0;
  event.HypsInvMass_PiP=-999.0;
  event.HypsInvMass_Phi=-999.0;
  for (int i=0; i<MaxHits; i++) {
    event.p[i]=-999.0;
    event.pK[i]=-999.0;
    event.pKCal[i]=-999.0;
    event.theta[i]=-999.0;
    event.MissMass[i]=-999.0;
    event.MissMass_others1[i]=-999.0;
    event.MissMass_others2[i]=-999.0;
    event.pSigmaCal[i]=-999.0;
    event.Vertex_x[i]=-999.0;
    event.Vertex_y[i]=-999.0;
    event.Vertex_z[i]=-999.0;

    event.pSksPart[i]=-999.0;
    event.yTofHyps[i]=-999.0;
    event.xTofHyps[i]=-999.0;
    event.m2[i]=-999.0;
  }

  event.u0SdcIn = -999.;
  event.v0SdcIn = -999.;

  event.resSftU = -999.;
  event.resSftV = -999.;

  event.nPiPi      =  0;
  event.nPiPi_CFT  =  0;
  event.nPiP_CFT   =  0;

  for(int i = 0; i<MaxHits; ++i){
    event.theta[i]  =  -999.;
    event.vtx_pipi_x[i]  =  -999.;
    event.vtx_pipi_y[i]  =  -999.;
    event.vtx_pipi_z[i]  =  -999.;
    event.cdist_pipi[i]  =  -999.;
    event.cdistBeam_pipi[i]  =  -999.;
    event.invmass[i]  =  -999.;
    event.theta_pipi[i]  =  -999.;
    event.missmass[i]  =  -999.;
    event.missmom[i]  =  -999.;
    event.momK0[i]  =  -999.;
    event.vtx_k0_x[i]  =  -999.;
    event.vtx_k0_y[i]  =  -999.;
    event.vtx_k0_z[i]  =  -999.;
    event.cdist_k0[i]  =  -999.;
    event.theta_k0[i]  =  -999.;
    event.momCalLambda[i]  =  -999.;
    event.momCalLambda_x[i]  =  -999.;
    event.momCalLambda_y[i]  =  -999.;
    event.momCalLambda_z[i]  =  -999.;
    event.missmass0  =  -999.;
    event.invmass0  =  -999.;

    event.vtx_pipi_cft_x[i]  =  -999.;
    event.vtx_pipi_cft_y[i]  =  -999.;
    event.vtx_pipi_cft_z[i]  =  -999.;
    event.cdist_cft_pipi[i]  =  -999.;
    event.cdistBeam_cft_pipi[i]  =  -999.;
    event.theta_cft_pipi[i]  =  -999.;

    event.vtx_pip_cft_x[i]  =  -999.;
    event.vtx_pip_cft_y[i]  =  -999.;
    event.vtx_pip_cft_z[i]  =  -999.;
    event.cdist_cft_pip[i]  =  -999.;
    event.cdistBeam_cft_pip[i]  =  -999.;
    event.theta_cft_pip[i]  =  -999.;
    event.invmass_cft_pip[i]  =  -999.;
  }


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
  event.decayPMomCal = -999.9;
  event.decayPMomVecCal[0] = -999.9;
  event.decayPMomVecCal[1] = -999.9;
  event.decayPMomVecCal[2] = -999.9;
  event.costDecayPMomCal = -999.9;
  event.costDecayPMomCalCM = -999.9;

  event.vertexScat[0] = -999.9;
  event.vertexScat[1] = -999.9;
  event.vertexScat[2] = -999.9;

  event.hypBeamVec[0] = -999.9;
  event.hypBeamVec[1] = -999.9;
  event.hypBeamVec[2] = -999.9;

  event.cdistScat = -999.9;
  event.thetaScat = -999.9;
  event.phiScat = -999.9;
  event.scatMomCal = -999.9;
  event.scatEkinCal = -999.9;
  event.thetaScatCM = -999.9;
  event.costScatOpenAngle = -999.9;
  event.thetaScatOpenAngle = -999.9;

  event.cdistDecayPi2 = -999.9;
  event.vertexDecayPi2[0] = -999.9;
  event.vertexDecayPi2[1] = -999.9;
  event.vertexDecayPi2[2] = -999.9;
  event.costDecayProtonLambdaPScat_cal = -999.9;
  event.costDecayProtonLambdaPScatCM_cal = -999.9;

  event.vertexScat_2nd[0] = -999.9;
  event.vertexScat_2nd[1] = -999.9;
  event.vertexScat_2nd[2] = -999.9;

  event.hypBeamVec_2nd[0] = -999.9;
  event.hypBeamVec_2nd[1] = -999.9;
  event.hypBeamVec_2nd[2] = -999.9;

  event.cdistScat_2nd = -999.9;
  event.thetaScat_2nd = -999.9;
  event.phiScat_2nd = -999.9;
  event.scatMomCal_2nd = -999.9;
  event.scatEkinCal_2nd = -999.9;
  event.thetaScatCM_2nd = -999.9;

  event.vertexDecayProton1_2p[0] = -999.9;
  event.vertexDecayProton1_2p[1] = -999.9;
  event.vertexDecayProton1_2p[2] = -999.9;
  event.cdistDecayProton1_2p = -999.9;
  event.MissMassPi_1_2p = -999.9;
  event.costDecayProtonLambdaPScat1_2p = -999.9;
  event.costDecayProtonLambdaPScatCM1_2p = -999.9;
  event.vertexDecayProton2_2p[0] = -999.9;
  event.vertexDecayProton2_2p[1] = -999.9;
  event.vertexDecayProton2_2p[2] = -999.9;
  event.cdistDecayProton2_2p = -999.9;
  event.MissMassPi_2_2p = -999.9;
  event.costDecayProtonLambdaPScat2_2p = -999.9;
  event.costDecayProtonLambdaPScatCM2_2p = -999.9;
  event.costSigmaDecayProton_2p = -999.9;

  event.momLambda1 = -999.;
  event.momLambda2 = -999.;
  event.vertexLambdaDecay[0] = -999.;
  event.vertexLambdaDecay[1] = -999.;
  event.vertexLambdaDecay[2] = -999.;
  event.cdistLambdaDecay     = -999.;
  event.cdistLambdaPScat_Lambda = -999.;
  event.vertexLambdaPScat_Lambda[0] = -999.;
  event.vertexLambdaPScat_Lambda[1] = -999.;
  event.vertexLambdaPScat_Lambda[2] = -999.;
  event.costLambdaPScat_Lambda = -999.;
  event.thetaLambdaPScat_Lambda = -999.;
  event.phiLambdaPScat_Lambda = -999.;
  event.thetaCMLambdaPScat_Lambda = -999.;
  event.momCalLambda1 = -999.;
  event.momCalLambda2 = -999.;
  event.momCalRecoilProton = -999.;
  event.momRecoilProton = -999.;
  event.deltaThetaRecoilProton = -999.;
  event.costDecayProtonLambdaPScat = -999.;
  event.costDecayProtonLambdaBeam = -999.;
  event.costDecayProtonLambdaPScatCM = -999.;
  event.costDecayProtonLambdaBeamCM = -999.;

  event.momLambda1_1st = -999.;
  event.momLambda2_1st = -999.;
  event.vertexLambdaDecay_1st[0] = -999.;
  event.vertexLambdaDecay_1st[1] = -999.;
  event.vertexLambdaDecay_1st[2] = -999.;
  event.cdistLambdaDecay_1st     = -999.;
  event.cdistLambdaPScat_Lambda_1st = -999.;
  event.vertexLambdaPScat_Lambda_1st[0] = -999.;
  event.vertexLambdaPScat_Lambda_1st[1] = -999.;
  event.vertexLambdaPScat_Lambda_1st[2] = -999.;
  event.costLambdaPScat_Lambda_1st = -999.;
  event.thetaLambdaPScat_Lambda_1st = -999.;
  event.phiLambdaPScat_Lambda_1st = -999.;
  event.thetaCMLambdaPScat_Lambda_1st = -999.;
  event.momCalLambda1_1st = -999.;
  event.momCalLambda2_1st = -999.;
  event.momCalRecoilProton_1st = -999.;
  event.momRecoilProton_1st = -999.;
  event.deltaThetaRecoilProton_1st = -999.;
  event.costDecayProtonLambdaPScat_1st = -999.;
  event.costDecayProtonLambdaPScatCM_1st = -999.;

  event.momLambda1_2nd = -999.;
  event.momLambda2_2nd = -999.;
  event.vertexLambdaDecay_2nd[0] = -999.;
  event.vertexLambdaDecay_2nd[1] = -999.;
  event.vertexLambdaDecay_2nd[2] = -999.;
  event.cdistLambdaDecay_2nd     = -999.;
  event.cdistLambdaPScat_Lambda_2nd = -999.;
  event.vertexLambdaPScat_Lambda_2nd[0] = -999.;
  event.vertexLambdaPScat_Lambda_2nd[1] = -999.;
  event.vertexLambdaPScat_Lambda_2nd[2] = -999.;
  event.costLambdaPScat_Lambda_2nd = -999.;
  event.thetaLambdaPScat_Lambda_2nd = -999.;
  event.phiLambdaPScat_Lambda_2nd = -999.;
  event.thetaCMLambdaPScat_Lambda_2nd = -999.;
  event.momCalLambda1_2nd = -999.;
  event.momCalLambda2_2nd = -999.;
  event.momCalRecoilProton_2nd = -999.;
  event.momRecoilProton_2nd = -999.;
  event.deltaThetaRecoilProton_2nd = -999.;
  event.costDecayProtonLambdaPScat_2nd = -999.;
  event.costDecayProtonLambdaPScatCM_2nd = -999.;

  event.scatMomCal_S0pConv1 = -999.;
  event.scatMomCal_S0pConv2 = -999.;
  event.scatEkinCal_S0pConv1 = -999.;
  event.scatEkinCal_S0pConv2 = -999.;
  event.thetaS0pConvCM = -999.;
  event.costScatOpenAngle_S0pConv = -999.;
  event.thetaScatOpenAngle_S0pConv = -999.;
  event.cdistDecayPi2_S0pConv = -999.;
  event.vertexDecayPi2_S0pConv[0] = -999.;
  event.vertexDecayPi2_S0pConv[1] = -999.;
  event.vertexDecayPi2_S0pConv[2] = -999.;

  event.vertexDecayProton_S0pConv1_2p[0] = -999.;
  event.vertexDecayProton_S0pConv1_2p[1] = -999.;
  event.vertexDecayProton_S0pConv1_2p[2] = -999.;
  event.cdistDecayProton_S0pConv1_2p = -999.;
  event.MissMassPiGamma_1_2p = -999.;

  event.scatMomCal_S0pConv1_2nd = -999.;
  event.scatMomCal_S0pConv2_2nd = -999.;
  event.scatEkinCal_S0pConv1_2nd = -999.;
  event.scatEkinCal_S0pConv2_2nd = -999.;
  event.thetaS0pConvCM_2nd = -999.;
  event.costScatOpenAngle_S0pConv_2nd = -999.;
  event.thetaScatOpenAngle_S0pConv_2nd = -999.;
  event.cdistDecayPi2_S0pConv_2nd = -999.;
  event.vertexDecayPi2_S0pConv_2nd[0] = -999.;
  event.vertexDecayPi2_S0pConv_2nd[1] = -999.;
  event.vertexDecayPi2_S0pConv_2nd[2] = -999.;

  event.vertexDecayProton_S0pConv1_2p_2nd[0] = -999.;
  event.vertexDecayProton_S0pConv1_2p_2nd[1] = -999.;
  event.vertexDecayProton_S0pConv1_2p_2nd[2] = -999.;
  event.cdistDecayProton_S0pConv1_2p_2nd = -999.;
  event.MissMassPiGamma_1_2p_2nd = -999.;


  event.energyBGO_p_lambda = -999.;
  event.energyBGO_pi_lambda = -999.;
  event.segPiID_pi_lambda = -1;
  event.invMass_PiP = -999.;
  event.invMass_PiP2 = -999.;
  event.theta_PiP = -999.;

  event.vertexK0Production[0] = -999.;
  event.vertexK0Production[1] = -999.;
  event.vertexK0Production[2] = -999.;
  event.cdistK0Production = -999.;
  event.thetaK0Production = -999.;
  event.momK0Production = -999.;
  event.vertexK0decay[0] = -999.;
  event.vertexK0decay[1] = -999.;
  event.vertexK0decay[2] = -999.;
  event.cdistK0decay = -999.;
  event.openingAngleK0decay = -999.;
  event.momCalLambdaProduction = -999.;
  event.momCalLambdaProductionVec[0] = -999.;
  event.momCalLambdaProductionVec[1] = -999.;
  event.momCalLambdaProductionVec[2] = -999.;

  event.vertexPPScat[0] =  -999.9;
  event.vertexPPScat[1] =  -999.9;
  event.vertexPPScat[2] =  -999.9;
  event.cdistTwoProton =  -999.9;
  event.thetaTwoProton =  -999.9;
  event.momLambdaProtonDecay =  -999.9;
  event.momLambdaProtonDecayVec[0] =  -999.9;
  event.momLambdaProtonDecayVec[1] =  -999.9;
  event.momLambdaProtonDecayVec[2] =  -999.9;
  event.vertexLambdaProtonDecay[0] =  -999.9;
  event.vertexLambdaProtonDecay[1] =  -999.9;
  event.vertexLambdaProtonDecay[2] =  -999.9;
  event.cdistLambdaProtonDecay =  -999.9;
  event.thetaLambdaProtonDecay =  -999.9;
  event.momCalLambdaProtonDecay1 =  -999.9;
  event.momCalLambdaProtonDecay2 =  -999.9;

  event.vertexPpScat[0] = -999.9;
  event.vertexPpScat[1] = -999.9;
  event.vertexPpScat[2] = -999.9;

  event.cdistPpScat = -999.9;
  event.thetaPpScat = -999.9;
  event.scatPpMomCal = -999.9;
  event.scatPpEkinCal = -999.9;
  event.thetaPpScatCM = -999.9;

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
  HB1(11, "N Hit SdcIn tracking", 15, 0, 15);
  HB1(12, "Chisquare SdcIn tracking", 100, 0, 50);
  HB1(13, "Hit Layer  SdcIn tracking", 15, 0, 15);
  HB1(14, "X0  SdcIn tracking", 200, -50, 50);
  HB1(15, "Y0  SdcIn tracking", 200, -50, 50);
  HB1(16, "U0  SdcIn tracking", 200, -0.5, 0.5);
  HB1(17, "V0  SdcIn tracking", 200, -0.5, 0.5);
  HB2(18, "X0%U0  SdcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(19, "Y0%V0  SdcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(20, "X0%Y0  SdcIn tracking", 200, -50, 50, 200, -50, 50);

  for (int l=1; l<NumOfLayersSdcIn+1; l++) {
    sprintf(buf, "Residual layer = %d", l);
    HB1(140+l, buf, 100, -1, 1);
    sprintf(buf, "Hit position layer = %d", l);
    HB1(160+l, buf, 100, -250, 250);
    sprintf(buf, "Hitpat layer = %d", l);
    HB1(600+l, buf, 130, 0, 130);
  }


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
  HB1(51, "N Hit Sks tracking", 20, 0, 20);
  HB1(52, "Chisquare Sks tracking", 400, 0, 100);
  HB1(53, "Hit Layer  Sks tracking", 40, 0, 40);
  HB1(54, "X0 (Target)  Sks tracking", 200, -50, 50);
  HB1(55, "Y0 (Target)  Sks tracking", 200, -50, 50);
  HB1(56, "U0  Sks tracking", 200, -0.5, 0.5);
  HB1(57, "V0  Sks tracking", 200, -0.5, 0.5);
  HB2(58, "X0%U0  Sks tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(59, "Y0%V0  Sks tracking", 200, -50, 50, 200, -0.5, 0.5);
  HB2(60, "X0%Y0  Sks tracking", 200, -50, 50, 200, -50, 50);
  HB1(61, "Momentum", 200, 0.4, 2.0);
  HB1(62, "Path Length", 1000, 0., 5000.);
  HB1(63, "Theta", 300, 0., 30.);
  HB1(64, "X Diff at TOF", 200, -1000., 1000.);
  HB1(65, "Mass square", 400, -0.5, 1.5);
  HB1(66, "Mass", 400, -0.5, 1.5);
  HB1(67, "MissMass", 400, 0.5, 1.5);
  HB1(68, "MissMass (Eloss Corr)", 400, 0.5, 1.5);
  HB1(69, "HypsInvMass_PiP", 400, 0.5, 2.0);
  HB1(70, "HypsInvMass_Phi", 400, 0.5, 1.5);

  for (int l=PlMinSdcOut; l<PlMaxSdcOut+1; l++) {
    sprintf(buf, "Residual layer = %d", l);
    HB1(200+l, buf, 100, -1, 1);
    sprintf(buf, "Hitpat layer = %d", l);
    HB1(700+l, buf, 150, 0, 150);
  }

  // HB1(70, "N track BcOut", 10, 0, 10);
  // HB1(71, "N Hit BcOut tracking", 20, 0, 20);
  // HB1(72, "Chisquare BcOut tracking", 100, 0, 50);
  // HB1(73, "Hit Layer  BcOut tracking", 20, 0, 20);
  // HB1(74, "X0  BcOut tracking", 200, -50, 50);
  // HB1(75, "Y0  BcOut tracking", 200, -50, 50);
  // HB1(76, "U0  BcOut tracking", 200, -0.5, 0.5);
  // HB1(77, "V0  BcOut tracking", 200, -0.5, 0.5);
  // HB2(78, "X0%U0  BcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  // HB2(79, "Y0%V0  BcIn tracking", 200, -50, 50, 200, -0.5, 0.5);
  // HB2(80, "X0%Y0  BcIn tracking", 200, -50, 50, 200, -50, 50);
  // HB1(81, "Residual SFT-U  BcOut tracking", 400, -50, 50);
  // HB1(82, "Residual SFT-V  BcOut tracking", 400, -50, 50);

  // for (int l=1; l<NumOfLayersBcOut+1; l++) {
  //   sprintf(buf, "Residual layer = %d", l);
  //   HB1(100+l, buf, 100, -2, 2);
  //   sprintf(buf, "Hit position layer = %d", l);
  //   HB1(120+l, buf, 100, -120, 120);
  // }

  for (int l=1; l<NumOfLayersSdcIn+1; l++) {
    sprintf(buf, "Residual layer = %d", l);
    HB1(100+l, buf, 100, -2, 2);
    sprintf(buf, "Hit position layer = %d", l);
    HB1(150+l, buf, 100, -120, 120);
  }

  for (int l=PlMinSdcOut; l<PlMaxSdcOut+1; l++) {
    sprintf(buf, "Residual layer = %d", l);
    HB1(100+l, buf, 100, -1, 1);
    sprintf(buf, "Hit position layer = %d", l);
    HB1(150+l, buf, 100, -120, 120);
  }

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
  HB1(437, "Delta Pos at PiV", 400, -100, 100);
  HB2(438, "Delta Pos at PiV", 100, 0, 100, 200, -200, 200);

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
  tree->Branch("costDecayProtonCM0",   &event.costDecayProtonCM0,  "costDecayProtonCM0/D");
  tree->Branch("momVectorPiMinusFromK0",   event.momVectorPiMinusFromK0,  "momVectorPiMinusFromK0[3]/D");
  tree->Branch("momPiMinusFromK0",   &event.momPiMinusFromK0,  "momPiMinusFromK0/D");

  tree->Branch("egamma", &event.egamma, "egamma/D");
  tree->Branch("ntSks",   &event.ntSks,  "ntSks/I");
  tree->Branch("p",   &event.p,  "p[ntSks]/D");
  tree->Branch("q",   &event.q,  "q[ntSks]/D");

  tree->Branch("ntSksPart",   &event.ntSksPart,  "ntSksPart/I");
  tree->Branch("pSksPart",   event.pSksPart,  "pSksPart[ntSksPart]/D");
  tree->Branch("yTofHyps",   event.yTofHyps,  "yTofHyps[ntSksPart]/D");
  tree->Branch("xTofHyps",   event.xTofHyps,  "xTofHyps[ntSksPart]/D");
  tree->Branch("m2",   event.m2,  "m2[ntSksPart]/D");

  tree->Branch("nK",   &event.nK,  "nK/I");
  tree->Branch("pK",   event.pK,  "pK[nK]/D");
  tree->Branch("MissMass",   event.MissMass,  "MissMass[nK]/D");
  tree->Branch("MissMass_others1",   event.MissMass_others1,  "MissMass_others1[nK]/D");
  tree->Branch("MissMass_others2",   event.MissMass_others2,  "MissMass_others2[nK]/D");
  tree->Branch("HypsInvMass_PiP",   &event.HypsInvMass_PiP,  "HypsInvMass_PiP/D");
  tree->Branch("HypsInvMass_Phi",   &event.HypsInvMass_Phi,  "HypsInvMass_Phi/D");
  tree->Branch("u0SdcIn",   &event.u0SdcIn,  "u0SdcIn/D");
  tree->Branch("v0SdcIn",   &event.v0SdcIn,  "v0SdcIn/D");

  tree->Branch("resSftU",   &event.resSftU,  "resSftU/D");
  tree->Branch("resSftV",   &event.resSftV,  "resSftV/D");

  tree->Branch("nLambda",  &event.nLambda,  "nLambda/I");
  tree->Branch("pKCal",   event.pKCal,  "pKCal[nLambda]/D");
  tree->Branch("theta",   event.theta,  "theta[nLambda]/D");
  tree->Branch("pSigmaCal",   event.pSigmaCal,  "pSigmaCal[nLambda]/D");
  tree->Branch("Vertex_x",   event.Vertex_x,  "Vertex_x[nLambda]/D");
  tree->Branch("Vertex_y",   event.Vertex_y,  "Vertex_y[nLambda]/D");
  tree->Branch("Vertex_z",   event.Vertex_z,  "Vertex_z[nLambda]/D");
  tree->Branch("missmass",   event.missmass, "missmass[nLambda]/D");
  tree->Branch("missmom",   event.missmom, "missmom[nLambda]/D");
  tree->Branch("momCalLambda",   event.momCalLambda, "momCalLambda[nLambda]/D");
  tree->Branch("momCalLambda_x",   event.momCalLambda_x, "momCalLambda_x[nLambda]/D");
  tree->Branch("momCalLambda_y",   event.momCalLambda_y, "momCalLambda_y[nLambda]/D");
  tree->Branch("momCalLambda_z",   event.momCalLambda_z, "momCalLambda_z[nLambda]/D");


  tree->Branch("nPiPi",      &event.nPiPi,     "nPiPi/I");
  tree->Branch("vtx_pipi_x",   event.vtx_pipi_x, "vtx_pipi_x[nPiPi]/D");
  tree->Branch("vtx_pipi_y",   event.vtx_pipi_y, "vtx_pipi_y[nPiPi]/D");
  tree->Branch("vtx_pipi_z",   event.vtx_pipi_z, "vtx_pipi_z[nPiPi]/D");
  tree->Branch("invmass",   event.invmass, "invmass[nPiPi]/D");
  tree->Branch("theta_pipi",   event.theta_pipi, "theta_pipi[nPiPi]/D");
  tree->Branch("cdist_pipi",   event.cdist_pipi, "cdist_pipi[nPiPi]/D");
  tree->Branch("cdistBeam_pipi",   event.cdistBeam_pipi, "cdistBeam_pipi[nPiPi]/D");
  tree->Branch("momK0",   event.momK0, "momK0[nPiPi]/D");
  tree->Branch("vtx_k0_x",   event.vtx_k0_x, "vtx_k0_x[nPiPi]/D");
  tree->Branch("vtx_k0_y",   event.vtx_k0_y, "vtx_k0_y[nPiPi]/D");
  tree->Branch("vtx_k0_z",   event.vtx_k0_z, "vtx_k0_z[nPiPi]/D");
  tree->Branch("theta_k0",   event.theta_k0, "theta_k0[nPiPi]/D");
  tree->Branch("cdist_k0",   event.cdist_k0, "cdist_k0[nPiPi]/D");
  tree->Branch("missmass0",   &event.missmass0, "missmass0/D");
  tree->Branch("invmass0",    &event.invmass0,  "invmass0/D");

  tree->Branch("momLambda1",   &event.momLambda1, "momLambda1/D");
  tree->Branch("momLambda2",   &event.momLambda2, "momLambda2/D");
  tree->Branch("vertexLambdaDecay",   event.vertexLambdaDecay, "vertexLambdaDecay[3]/D");
  tree->Branch("cdistLambdaDecay",   &event.cdistLambdaDecay, "cdistLambdaDecay/D");
  tree->Branch("cdistLambdaPScat_Lambda",   &event.cdistLambdaPScat_Lambda, "cdistLambdaPScat_Lambda/D");
  tree->Branch("vertexLambdaPScat_Lambda",   event.vertexLambdaPScat_Lambda, "vertexLambdaPScat_Lambda[3]/D");
  tree->Branch("costLambdaPScat_Lambda",   &event.costLambdaPScat_Lambda, "costLambdaPScat_Lambda/D");
  tree->Branch("thetaLambdaPScat_Lambda",   &event.thetaLambdaPScat_Lambda, "thetaLambdaPScat_Lambda/D");
  tree->Branch("phiLambdaPScat_Lambda",     &event.phiLambdaPScat_Lambda,   "phiLambdaPScat_Lambda/D");
  tree->Branch("thetaCMLambdaPScat_Lambda",   &event.thetaCMLambdaPScat_Lambda, "thetaCMLambdaPScat_Lambda/D");
  tree->Branch("momCalLambda1",   &event.momCalLambda1, "momCalLambda1/D");
  tree->Branch("momCalLambda2",   &event.momCalLambda2, "momCalLambda2/D");
  tree->Branch("momCalRecoilProton",   &event.momCalRecoilProton, "momCalRecoilProton/D");
  tree->Branch("momRecoilProton",   &event.momRecoilProton, "momRecoilProton/D");
  tree->Branch("deltaThetaRecoilProton",   &event.deltaThetaRecoilProton, "deltaThetaRecoilProton/D");
  tree->Branch("costDecayProtonLambdaPScat",   &event.costDecayProtonLambdaPScat, "costDecayProtonLambdaPScat/D");
  tree->Branch("costDecayProtonLambdaBeam",   &event.costDecayProtonLambdaBeam, "costDecayProtonLambdaBeam/D");
  tree->Branch("costDecayProtonLambdaPScatCM",   &event.costDecayProtonLambdaPScatCM, "costDecayProtonLambdaPScatCM/D");
  tree->Branch("costDecayProtonLambdaBeamCM",   &event.costDecayProtonLambdaBeamCM, "costDecayProtonLambdaBeamCM/D");

  tree->Branch("momLambda1_1st",   &event.momLambda1_1st, "momLambda1_1st/D");
  tree->Branch("momLambda2_1st",   &event.momLambda2_1st, "momLambda2_1st/D");
  tree->Branch("vertexLambdaDecay_1st",   event.vertexLambdaDecay_1st, "vertexLambdaDecay_1st[3]/D");
  tree->Branch("cdistLambdaDecay_1st",   &event.cdistLambdaDecay_1st, "cdistLambdaDecay_1st/D");
  tree->Branch("cdistLambdaPScat_Lambda_1st",   &event.cdistLambdaPScat_Lambda_1st, "cdistLambdaPScat_Lambda_1st/D");
  tree->Branch("vertexLambdaPScat_Lambda_1st",   event.vertexLambdaPScat_Lambda_1st, "vertexLambdaPScat_Lambda_1st[3]/D");
  tree->Branch("costLambdaPScat_Lambda_1st",   &event.costLambdaPScat_Lambda_1st, "costLambdaPScat_Lambda_1st/D");
  tree->Branch("thetaLambdaPScat_Lambda_1st",   &event.thetaLambdaPScat_Lambda_1st, "thetaLambdaPScat_Lambda_1st/D");
  tree->Branch("phiLambdaPScat_Lambda_1st",     &event.phiLambdaPScat_Lambda_1st,   "phiLambdaPScat_Lambda_1st/D");
  tree->Branch("thetaCMLambdaPScat_Lambda_1st",   &event.thetaCMLambdaPScat_Lambda_1st, "thetaCMLambdaPScat_Lambda_1st/D");
  tree->Branch("momCalLambda1_1st",   &event.momCalLambda1_1st, "momCalLambda1_1st/D");
  tree->Branch("momCalLambda2_1st",   &event.momCalLambda2_1st, "momCalLambda2_1st/D");
  tree->Branch("momCalRecoilProton_1st",   &event.momCalRecoilProton_1st, "momCalRecoilProton_1st/D");
  tree->Branch("momRecoilProton_1st",   &event.momRecoilProton_1st, "momRecoilProton_1st/D");
  tree->Branch("deltaThetaRecoilProton_1st",   &event.deltaThetaRecoilProton_1st, "deltaThetaRecoilProton_1st/D");
  tree->Branch("costDecayProtonLambdaPScat_1st",   &event.costDecayProtonLambdaPScat_1st, "costDecayProtonLambdaPScat_1st/D");
  tree->Branch("costDecayProtonLambdaPScatCM_1st",   &event.costDecayProtonLambdaPScatCM_1st, "costDecayProtonLambdaPScatCM_1st/D");

  tree->Branch("momLambda1_2nd",   &event.momLambda1_2nd, "momLambda1_2nd/D");
  tree->Branch("momLambda2_2nd",   &event.momLambda2_2nd, "momLambda2_2nd/D");
  tree->Branch("vertexLambdaDecay_2nd",   event.vertexLambdaDecay_2nd, "vertexLambdaDecay_2nd[3]/D");
  tree->Branch("cdistLambdaDecay_2nd",   &event.cdistLambdaDecay_2nd, "cdistLambdaDecay_2nd/D");
  tree->Branch("cdistLambdaPScat_Lambda_2nd",   &event.cdistLambdaPScat_Lambda_2nd, "cdistLambdaPScat_Lambda_2nd/D");
  tree->Branch("vertexLambdaPScat_Lambda_2nd",   event.vertexLambdaPScat_Lambda_2nd, "vertexLambdaPScat_Lambda_2nd[3]/D");
  tree->Branch("costLambdaPScat_Lambda_2nd",   &event.costLambdaPScat_Lambda_2nd, "costLambdaPScat_Lambda_2nd/D");
  tree->Branch("thetaLambdaPScat_Lambda_2nd",   &event.thetaLambdaPScat_Lambda_2nd, "thetaLambdaPScat_Lambda_2nd/D");
  tree->Branch("phiLambdaPScat_Lambda_2nd",     &event.phiLambdaPScat_Lambda_2nd,   "phiLambdaPScat_Lambda_2nd/D");
  tree->Branch("thetaCMLambdaPScat_Lambda_2nd",   &event.thetaCMLambdaPScat_Lambda_2nd, "thetaCMLambdaPScat_Lambda_2nd/D");
  tree->Branch("momCalLambda1_2nd",   &event.momCalLambda1_2nd, "momCalLambda1_2nd/D");
  tree->Branch("momCalLambda2_2nd",   &event.momCalLambda2_2nd, "momCalLambda2_2nd/D");
  tree->Branch("momCalRecoilProton_2nd",   &event.momCalRecoilProton_2nd, "momCalRecoilProton_2nd/D");
  tree->Branch("momRecoilProton_2nd",   &event.momRecoilProton_2nd, "momRecoilProton_2nd/D");
  tree->Branch("deltaThetaRecoilProton_2nd",   &event.deltaThetaRecoilProton_2nd, "deltaThetaRecoilProton_2nd/D");
  tree->Branch("costDecayProtonLambdaPScat_2nd",   &event.costDecayProtonLambdaPScat_2nd, "costDecayProtonLambdaPScat_2nd/D");
  tree->Branch("costDecayProtonLambdaPScatCM_2nd",   &event.costDecayProtonLambdaPScatCM_2nd, "costDecayProtonLambdaPScatCM_2nd/D");

  tree->Branch("energyBGO_p_lambda",   &event.energyBGO_p_lambda, "energyBGO_p_lambda/D");
  tree->Branch("energyBGO_pi_lambda",   &event.energyBGO_pi_lambda, "energyBGO_pi_lambda/D");
  tree->Branch("segPiID_pi_lambda",   &event.segPiID_pi_lambda, "segPiID_pi_lambda/I");
  tree->Branch("invMass_PiP",   &event.invMass_PiP, "invMass_PiP/D");
  tree->Branch("invMass_PiP2",   &event.invMass_PiP2, "invMass_PiP2/D");
  tree->Branch("theta_PiP",     &event.theta_PiP,   "theta_PiP/D");

  tree->Branch("vertexK0Production",   event.vertexK0Production, "vertexK0Production[3]/D");
  tree->Branch("cdistK0Production",   &event.cdistK0Production, "cdistK0Production/D");
  tree->Branch("thetaK0Production",   &event.thetaK0Production, "thetaK0Production/D");
  tree->Branch("momK0Production",   &event.momK0Production, "momK0Production/D");
  tree->Branch("vertexK0decay",   event.vertexK0decay, "vertexK0decay[3]/D");
  tree->Branch("cdistK0decay",   &event.cdistK0decay, "cdistK0decay/D");
  tree->Branch("openingAngleK0decay",   &event.openingAngleK0decay, "openingAngleK0decay/D");
  tree->Branch("momCalLambdaProduction",   &event.momCalLambdaProduction, "momCalLambdaProduction/D");
  tree->Branch("momCalLambdaProductionVec",   event.momCalLambdaProductionVec, "momCalLambdaProductionVec[3]/D");

  tree->Branch("vertexPPScat",   event.vertexPPScat,  "vertexPPScat[3]/D");
  tree->Branch("cdistTwoProton",   &event.cdistTwoProton,  "cdistTwoProton/D");
  tree->Branch("thetaTwoProton",   &event.thetaTwoProton,  "thetaTwoProton/D");
  tree->Branch("momLambdaProtonDecay",   &event.momLambdaProtonDecay,  "momLambdaProtonDecay/D");
  tree->Branch("momLambdaProtonDecayVec",   event.momLambdaProtonDecayVec,  "momLambdaProtonDecayVec[3]/D");
  tree->Branch("vertexLambdaProtonDecay",   event.vertexLambdaProtonDecay,  "vertexLambdaProtonDecay[3]/D");
  tree->Branch("cdistLambdaProtonDecay",   &event.cdistLambdaProtonDecay,  "cdistLambdaProtonDecay/D");
  tree->Branch("thetaLambdaProtonDecay",   &event.thetaLambdaProtonDecay,  "thetaLambdaProtonDecay/D");
  tree->Branch("momCalLambdaProtonDecay1",   &event.momCalLambdaProtonDecay1,  "momCalLambdaProtonDecay1/D");
  tree->Branch("momCalLambdaProtonDecay2",   &event.momCalLambdaProtonDecay2,  "momCalLambdaProtonDecay2/D");

  tree->Branch("vertexPpScat",   event.vertexPpScat,  "vertexPpScat[3]/D");
  tree->Branch("cdistPpScat",   &event.cdistPpScat,  "cdistPpScat/D");
  tree->Branch("thetaPpScat",   &event.thetaPpScat,  "thetaPpScat/D");
  tree->Branch("scatPpMomCal",   &event.scatPpMomCal,  "scatPpMomCal/D");
  tree->Branch("scatPpEkinCal",   &event.scatPpEkinCal,  "scatPpEkinCal/D");
  tree->Branch("thetaPpScatCM",   &event.thetaPpScatCM,  "thetaPpScatCM/D");

  tree->Branch("nPiPi_CFT",      &event.nPiPi_CFT,     "nPiPi_CFT/I");
  tree->Branch("vtx_pipi_cft_x",   event.vtx_pipi_cft_x, "vtx_pipi_cft_x[nPiPi_CFT]/D");
  tree->Branch("vtx_pipi_cft_y",   event.vtx_pipi_cft_y, "vtx_pipi_cft_y[nPiPi_CFT]/D");
  tree->Branch("vtx_pipi_cft_z",   event.vtx_pipi_cft_z, "vtx_pipi_cft_z[nPiPi_CFT]/D");
  tree->Branch("theta_cft_pipi",   event.theta_cft_pipi, "theta_cft_pipi[nPiPi_CFT]/D");
  tree->Branch("cdist_cft_pipi",   event.cdist_cft_pipi, "cdist_cft_pipi[nPiPi_CFT]/D");
  tree->Branch("cdistBeam_cft_pipi",   event.cdistBeam_cft_pipi, "cdistBeam_cft_pipi[nPiPi_CFT]/D");

  tree->Branch("nPiP_CFT",      &event.nPiP_CFT,     "nPiP_CFT/I");
  tree->Branch("vtx_pip_cft_x",   event.vtx_pip_cft_x, "vtx_pip_cft_x[nPiP_CFT]/D");
  tree->Branch("vtx_pip_cft_y",   event.vtx_pip_cft_y, "vtx_pip_cft_y[nPiP_CFT]/D");
  tree->Branch("vtx_pip_cft_z",   event.vtx_pip_cft_z, "vtx_pip_cft_z[nPiP_CFT]/D");
  tree->Branch("theta_cft_pip",   event.theta_cft_pip, "theta_cft_pip[nPiP_CFT]/D");
  tree->Branch("cdist_cft_pip",   event.cdist_cft_pip, "cdist_cft_pip[nPiP_CFT]/D");
  tree->Branch("cdistBeam_cft_pip",   event.cdistBeam_cft_pip, "cdistBeam_cft_pip[nPiP_CFT]/D");
  tree->Branch("invmass_cft_pip",   event.invmass_cft_pip, "invmass_cft_pip[nPiP_CFT]/D");


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
  tree->Branch("decayPMomCal",   &event.decayPMomCal,  "decayPMomCal/D");
  tree->Branch("decayPMomVecCal",   event.decayPMomVecCal,  "decayPMomVecCal[3]/D");
  tree->Branch("costDecayPMomCal",  &event.costDecayPMomCal,  "costDecayPMomCal/D");
  tree->Branch("costDecayPMomCalCM",  &event.costDecayPMomCalCM,  "costDecayPMomCalCM/D");

  tree->Branch("vertexScat",   event.vertexScat,  "vertexScat[3]/D");
  tree->Branch("hypBeamVec",   event.hypBeamVec,  "hypBeamVec[3]/D");
  tree->Branch("cdistScat",   &event.cdistScat,  "cdistScat/D");
  tree->Branch("thetaScat",   &event.thetaScat,  "thetaScat/D");
  tree->Branch("phiScat",     &event.phiScat,    "phiScat/D");
  tree->Branch("scatMomCal",   &event.scatMomCal,  "scatMomCal/D");
  tree->Branch("scatEkinCal",   &event.scatEkinCal,  "scatEkinCal/D");
  tree->Branch("thetaScatCM",   &event.thetaScatCM,  "thetaScatCM/D");
  tree->Branch("costScatOpenAngle", &event.costScatOpenAngle,  "costScatOpenAngle/D");
  tree->Branch("thetaScatOpenAngle", &event.thetaScatOpenAngle,  "thetaScatOpenAngle/D");

  tree->Branch("cdistDecayPi2",    &event.cdistDecayPi2,  "cdistDecayPi2/D");
  tree->Branch("vertexDecayPi2",   event.vertexDecayPi2,  "vertexDecayPi2[3]/D");
  tree->Branch("costDecayProtonLambdaPScat_cal",    &event.costDecayProtonLambdaPScat_cal,  "costDecayProtonLambdaPScat_cal/D");
  tree->Branch("costDecayProtonLambdaPScatCM_cal",    &event.costDecayProtonLambdaPScatCM_cal,  "costDecayProtonLambdaPScatCM_cal/D");


  tree->Branch("vertexScat_2nd",   event.vertexScat_2nd,  "vertexScat_2nd[3]/D");
  tree->Branch("hypBeamVec_2nd",   event.hypBeamVec_2nd,  "hypBeamVec_2nd[3]/D");
  tree->Branch("cdistScat_2nd",   &event.cdistScat_2nd,  "cdistScat_2nd/D");
  tree->Branch("thetaScat_2nd",   &event.thetaScat_2nd,  "thetaScat_2nd/D");
  tree->Branch("phiScat_2nd",     &event.phiScat_2nd,    "phiScat_2nd/D");
  tree->Branch("scatMomCal_2nd",   &event.scatMomCal_2nd,  "scatMomCal_2nd/D");
  tree->Branch("scatEkinCal_2nd",   &event.scatEkinCal_2nd,  "scatEkinCal_2nd/D");
  tree->Branch("thetaScatCM_2nd",   &event.thetaScatCM_2nd,  "thetaScatCM_2nd/D");

  tree->Branch("vertexDecayProton1_2p",   event.vertexDecayProton1_2p,  "vertexDecayProton1_2p[3]/D");
  tree->Branch("cdistDecayProton1_2p",   &event.cdistDecayProton1_2p,  "cdistDecayProton1_2p/D");
  tree->Branch("MissMassPi_1_2p",   &event.MissMassPi_1_2p,  "MissMassPi_1_2p/D");
  tree->Branch("costDecayProtonLambdaPScat1_2p",   &event.costDecayProtonLambdaPScat1_2p,  "costDecayProtonLambdaPScat1_2p/D");
  tree->Branch("costDecayProtonLambdaPScatCM1_2p",   &event.costDecayProtonLambdaPScatCM1_2p,  "costDecayProtonLambdaPScatCM1_2p/D");
  tree->Branch("vertexDecayProton2_2p",   event.vertexDecayProton2_2p,  "vertexDecayProton2_2p[3]/D");
  tree->Branch("cdistDecayProton2_2p",   &event.cdistDecayProton2_2p,  "cdistDecayProton2_2p/D");
  tree->Branch("MissMassPi_2_2p",   &event.MissMassPi_2_2p,  "MissMassPi_2_2p/D");
  tree->Branch("costDecayProtonLambdaPScat2_2p",   &event.costDecayProtonLambdaPScat2_2p,  "costDecayProtonLambdaPScat2_2p/D");
  tree->Branch("costDecayProtonLambdaPScatCM2_2p",   &event.costDecayProtonLambdaPScatCM2_2p,  "costDecayProtonLambdaPScatCM2_2p/D");
  tree->Branch("costSigmaDecayProton_2p",   &event.costSigmaDecayProton_2p,  "costSigmaDecayProton_2p/D");


  tree->Branch("scatMomCal_S0pConv1", 	     &event.scatMomCal_S0pConv1, 	   "scatMomCal_S0pConv1/D");
  tree->Branch("scatMomCal_S0pConv2", 	     &event.scatMomCal_S0pConv2, 	   "scatMomCal_S0pConv2/D");
  tree->Branch("scatEkinCal_S0pConv1", 	     &event.scatEkinCal_S0pConv1, 	   "scatEkinCal_S0pConv1/D");
  tree->Branch("scatEkinCal_S0pConv2", 	     &event.scatEkinCal_S0pConv2, 	   "scatEkinCal_S0pConv2/D");
  tree->Branch("thetaS0pConvCM", 	     &event.thetaS0pConvCM, 	     	   "thetaS0pConvCM/D");
  tree->Branch("costScatOpenAngle_S0pConv",  &event.costScatOpenAngle_S0pConv,  "costScatOpenAngle_S0pConv/D");
  tree->Branch("thetaScatOpenAngle_S0pConv", &event.thetaScatOpenAngle_S0pConv, "thetaScatOpenAngle_S0pConv/D");
  tree->Branch("cdistDecayPi2_S0pConv",      &event.cdistDecayPi2_S0pConv,      "cdistDecayPi2_S0pConv/D");
  tree->Branch("vertexDecayPi2_S0pConv",     event.vertexDecayPi2_S0pConv,      "vertexDecayPi2_S0pConv[3]/D");

  tree->Branch("vertexDecayProton_S0pConv1_2p", event.vertexDecayProton_S0pConv1_2p, "vertexDecayProton_S0pConv1_2p[3]/D");
  tree->Branch("cdistDecayProton_S0pConv1_2p",	&event.cdistDecayProton_S0pConv1_2p,     "cdistDecayProton_S0pConv1_2p/D");
  tree->Branch("MissMassPiGamma_1_2p",          &event.MissMassPiGamma_1_2p,             "MissMassPiGamma_1_2p/D");

  tree->Branch("scatMomCal_S0pConv1_2nd",	&event.scatMomCal_S0pConv1_2nd,	  "scatMomCal_S0pConv1_2nd/D");
  tree->Branch("scatMomCal_S0pConv2_2nd",	&event.scatMomCal_S0pConv2_2nd,	  "scatMomCal_S0pConv2_2nd/D");
  tree->Branch("scatEkinCal_S0pConv1_2nd",	&event.scatEkinCal_S0pConv1_2nd,	  "scatEkinCal_S0pConv1_2nd/D");
  tree->Branch("scatEkinCal_S0pConv2_2nd",	&event.scatEkinCal_S0pConv2_2nd,	  "scatEkinCal_S0pConv2_2nd/D");
  tree->Branch("thetaS0pConvCM_2nd",		&event.thetaS0pConvCM_2nd,		  "thetaS0pConvCM_2nd/D");
  tree->Branch("costScatOpenAngle_S0pConv_2nd",	&event.costScatOpenAngle_S0pConv_2nd,  "costScatOpenAngle_S0pConv_2nd/D");
  tree->Branch("thetaScatOpenAngle_S0pConv_2nd",&event.thetaScatOpenAngle_S0pConv_2nd, "thetaScatOpenAngle_S0pConv_2nd/D");
  tree->Branch("cdistDecayPi2_S0pConv_2nd",	&event.cdistDecayPi2_S0pConv_2nd,	  "cdistDecayPi2_S0pConv_2nd/D");
  tree->Branch("vertexDecayPi2_S0pConv_2nd",    event.vertexDecayPi2_S0pConv_2nd,  "vertexDecayPi2_S0pConv_2nd[3]/D");

  tree->Branch("vertexDecayProton_S0pConv1_2p_2nd", event.vertexDecayProton_S0pConv1_2p_2nd, "vertexDecayProton_S0pConv1_2p_2nd[3]/D");
  tree->Branch("cdistDecayProton_S0pConv1_2p_2nd",  &event.cdistDecayProton_S0pConv1_2p_2nd,     "cdistDecayProton_S0pConv1_2p_2nd/D");
  tree->Branch("MissMassPiGamma_1_2p_2nd",          &event.MissMassPiGamma_1_2p_2nd,             "MissMassPiGamma_1_2p_2nd/D");


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
  std::cout << "DCGeomFileName_ = "<< DCGeomFileName_ << std::endl;
  std::cout << "CftEffFileName_ = "<< CftEffFileName_ << std::endl;

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

  // if( mtTrigFileName_!="" )
  //   MatrixTrigMan::GetInstance().Initialize( mtTrigFileName_ );

  // if( sftTrigFileName_!="" )
  //   MatrixTrigMan::GetInstance().InitializeSftTrig( sftTrigFileName_ );

  // if( massTrigFileName_!="" )
  //   MatrixTrigMan::GetInstance().InitializeMassTrig( massTrigFileName_ );


  if( CftEffFileName_!="" )
    CFTTrackingEffMan::GetInstance().Initialize( CftEffFileName_ );


  return true;
}
