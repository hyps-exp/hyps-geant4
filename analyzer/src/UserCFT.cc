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
#include "CFTLocalTrack.hh"
#include "CFTParticle.hh"
#include "TemplateLib.hh"

#include "HistHelper.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "EvDispCFT.hh"

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


bool ProcessOneEvent( std::ifstream & );
void DefineHistograms( const char * );
void InitializeEvent( void );

const double MaxChiSqr = 100.;
const double PIni = 1.0;

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);


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

  int ntCFT;
  double BGO_Edep[MaxHits2];
  double TotalEdep[MaxHits2];
  double CFT_TotalEdep[MaxHits2];
  double CFT_NormTotalEdep[MaxHits2];  
  double CFT_MaxFiberEdep [NumOfPlaneCFT][MaxHits2];
  double PiV_Edep[MaxHits2];
  double theta[MaxHits2];

  int nP_CFT;
  int nPi_CFT;

  int BGO_ClustNum_Pi[MaxHits2];
  int  BGO_ClustNumNoTrack_Pi[MaxHits2];

  int    nhBGO_NoTrack;
  double MaxBGOEdep_NoTrack;
  double TotalBGOEdep_NoTrack;

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

  const int GeneratedNumInEffStudy = 1000;
  int EffStudyEvent = events/GeneratedNumInEffStudy;

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
  hodoAna->DecodeBGOHits(rawData);
  hodoAna->DecodePiVHits(rawData);

  for(int l = 0; l<NumOfLayersCFT; ++l){
    int nhit = hodoAna->GetNHitsCFT(l);

    int hid = 100+10*l;
    HF1(hid+1, nhit);

    for(int i = 0; i<nhit; ++i){
      const CFTFiberHit* hit = hodoAna->GetHitCFT(l, i);
      int mhit = hit->GetNumOfHit();
      int seg = hit->SegmentId();
      int pair_id = hit->PairId();

      bool flag_ctime=false;
      for(int m = 0; m<mhit; ++m){
	double ctime = hit->GetCTime(m);

	HF1 (hid+2, pair_id);
	HF1 (hid+3, ctime);
      }

      double pe_hi = hit->GetAHiGain();
      HF1 (hid+4, pe_hi);

      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	evDisp.ShowHitFiber(l, seg, pe_hi*100);
      }

    }
  }


  for(int l = 0; l<NumOfLayersCFT; ++l){
    int nhit = hodoAna->GetNClustersCFT(l);
    event.FiberHits[l] = nhit;

    int hid = 200+10*l;
    HF1(hid+1, nhit);

    for(int n = 0; n<nhit; ++n){
      CFTFiberCluster* cl = hodoAna->GetClusterCFT(l, n);
      int sizeCl    = cl->ClusterSize();
      double cmt    = cl->CMeanTime();
      //double cpos   = cl->MeanPosition();

      double max_pe   = cl->MaxPhotonNum();
      double total_pe   = cl->TotalPhotonNum();

      event.FiberSeg[l][n] = cl->MeanPairId();
      event.FiberTime[l][n] = cmt;
      event.FiberEdep[l][n] = total_pe;

      HF1 (hid+3, sizeCl);
      HF1 (hid+4, cmt);
      HF1 (hid+5, total_pe);
      HF1 (hid+6, max_pe);

      double mean_seg = cl->MeanPairId();
      HF1 (hid+2, mean_seg);
      
    }
  }

  int nhitBGO = hodoAna->GetNHitsBGO();
  HF1(301, nhitBGO);
  for (int n=0; n<nhitBGO; n++) {
    const Hodo2Hit* hit = hodoAna->GetHitBGO(n);
    int segment = hit->SegmentId();
    double time = hit->CMeanTime();
    double dE   = hit->DeltaE();

    HF1(302, segment);
    HF1(303, time);
    HF1(304, dE);
    HF2(305, event.thetaMeson, dE); // thetaMeson is proton angle here

    if (FlagEvDisp) {
      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      evDisp.ShowHitBGO(segment, dE);
    }

  }


  int nclBGO = hodoAna->GetNClustersBGO();
  HF1(311, nclBGO);
  event.CrystalHits = nclBGO;
  for (int n=0; n<nclBGO; n++) {
    const HodoCluster* cl = hodoAna->GetClusterBGO(n);
    int csize   = cl->ClusterSize();
    double segment = cl->MeanSeg();
    double time = cl->CMeanTime();
    double dE   = cl->DeltaE();

    event.CrystalSeg[n] = segment;
    event.CrystalTime[n] = time;
    event.CrystalEdep[n] = dE;

    HF1(312, segment);
    HF1(313, csize);
    HF1(314, time);
    HF1(315, dE);
    HF2(316, event.thetaMeson, dE); // thetaMeson is proton angle here
  }

  int nhitPiV = hodoAna->GetNHitsPiV();
  HF1(321, nhitPiV);
  for (int n=0; n<nhitPiV; n++) {
    const Hodo2Hit* hit = hodoAna->GetHitPiV(n);
    int segment = hit->SegmentId();
    int layer   = hit->PlaneId();
    double time = hit->CMeanTime();
    double dE   = hit->DeltaE();

    HF1(322, segment);
    HF1(323, time);
    HF1(324, dE);
    HF2(325, event.thetaMeson, dE); // thetaMeson is proton angle here


    if (FlagEvDisp) {
      const EvDispCFT & evDisp = EvDispCFT::GetInstance();
      evDisp.ShowHitPiV(layer, segment, dE);
    }

  }


  DCAnalyzer *DCAna = new DCAnalyzer;
  DCAna->DecodeCFTHits(rawData);
  DCAna->TrackSearchCFT();

  std::vector <CFTParticle *> CFTPartCont;
  std::vector <CFTParticle *> CFTProtonCont;
  std::vector <CFTParticle *> CFTPionCont;

  int ntCFT = DCAna->GetNtracksCFT();
  HF1(401, ntCFT);
  event.ntCFT = ntCFT;
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

    double cost = Dir.z()/Dir.mag();
    double theta = acos(cost)*Rad2Deg;


    double R=160.;
    double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());
    double B=(Dir.x()*Pos0.x()+Dir.y()*Pos0.y());
    double C=Pos0.x()*Pos0.x()+Pos0.y()*Pos0.y()-R*R;
    double t1, t2;
    double trackPosBGO[3]={-999., -999., -999};
    if (B*B-A*C>=0) {
      t1 = (-B+sqrt(B*B-A*C))/A;
      t2 = (-B-sqrt(B*B-A*C))/A;

      if (t1>=0 && t2<=0) {
	trackPosBGO[0] = Dir.x()*t1 + Pos0.x();
	trackPosBGO[1] = Dir.y()*t1 + Pos0.y();
	trackPosBGO[2] = Dir.z()*t1 + Pos0.z();
      } else if (t2>=0 && t1<=0) {
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
	//evDisp.DrawTrackInXYPlane(x1, y1, x2, y2);
	evDisp.DrawTrackInXYPlane(Pos0.x(), Pos0.y(), trackPosBGO[0], trackPosBGO[1]);
      }
    } else if (xyFitFlag==1) {
      double y1=-100., y2=100.;
      double x1 = Axy*y1+Bxy;
      double x2 = Axy*y2+Bxy;
      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	//evDisp.DrawTrackInXYPlane(x1, y1, x2, y2);
	evDisp.DrawTrackInXYPlane(Pos0.x(), Pos0.y(), trackPosBGO[0], trackPosBGO[1]);
      }
    }

    CFTParticle * CFTPart = new CFTParticle(tp);

    double dE_BGO=0.;

    for (int n=0; n<nhitBGO; n++) {
      //const Hodo2Hit* hit = hodoAna->GetHitBGO(n);
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

      double u=Dir.x(), v=Dir.y();
      double x0=Pos0.x(), y0=Pos0.y();
      double t = (u*(x-x0)+v*(y-y0))/(u*u+v*v);

      if (t>=0) {
	HF1(430, dist);
	
	if (fabs(dist)<30 && time>0 && time<5) {
	  dE_PiV += dE;
	  CFTPart->AddPiVHit(hit);
	}
      }
    }    
    if (nhitPiV>0)
      HF1(431, dE_PiV);

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
      event.CFT_TotalEdep[it] = totalDE;
      event.CFT_NormTotalEdep[it] = normalizedTotalDE;
      event.BGO_Edep[it] = dE_BGO;
      event.PiV_Edep[it] = dE_PiV;
      event.TotalEdep[it] = dE_BGO+totalDE;

      event.theta[it] = theta;

      int nh = tp->GetNHit();
      for (int j=0; j<nh; j++) {
	CFTFiberCluster *fcl = tp->GetHit(j);
	int layer = fcl->GetTrackingLayer() - 61;

	event.CFT_MaxFiberEdep[layer][it] = fcl->MaxDEHiGain();
      }

      int nhU = tp->GetNHitU();
      for (int j=0; j<nhU; j++) {
	CFTFiberCluster *fcl = tp->GetHitU(j);
	int layer = fcl->GetTrackingLayer() - 61;

	event.CFT_MaxFiberEdep[layer][it] = fcl->MaxDEHiGain();
      }


    }

    if (dE_PiV>1.) {
        HF2(432, dE_BGO, totalDE);
	HF2(433, dE_BGO, maxDE);
	HF2(434, dE_BGO, normalizedTotalDE);
	HF2(435, dE_BGO, normalizedMaxDE);
    } else {
        HF2(436, dE_BGO, totalDE);
	HF2(437, dE_BGO, maxDE);
	HF2(438, dE_BGO, normalizedTotalDE);
	HF2(439, dE_BGO, normalizedMaxDE);
    }

    if (tp->GetZTrackFlag()==0 || tp->GetZTrackFlag()==1) {
      double x1 = 100., x2 = -100;
      double z1 = (x1-x0)/u0, z2 = (x2-x0)/u0;

      double y1 = v0*z1+y0, y2 = v0*z2+y0;

      if (FlagEvDisp) {
	const EvDispCFT & evDisp = EvDispCFT::GetInstance();
	//evDisp.DrawTrackInZXPlane(z1, x1, z2, x2);
	//evDisp.DrawTrackInZYPlane(z1, y1, z2, y2);
	evDisp.DrawTrackInZXPlane(Pos0.z(), Pos0.x(), trackPosBGO[2], trackPosBGO[0]);
	evDisp.DrawTrackInZYPlane(Pos0.z(), Pos0.y(), trackPosBGO[2], trackPosBGO[1]);
      }
    }


    if (CFTPart->GetMass() > 0.9)
      CFTProtonCont.push_back(CFTPart);
    else if (CFTPart->GetMass() > 0.0 && CFTPart->GetMass() < 0.2)
      CFTPionCont.push_back(CFTPart);

  }

  std::vector <ThreeVector> CFTVtxCont;
  int ntCFT_2nd = 0;

  if (CFTProtonCont.size() == 0 || CFTPionCont.size() == 0 ) {
    CFTVtxCont.push_back(ThreeVector(event.scatPos0[0], event.scatPos0[1], event.scatPos0[2]));
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
	/*	
	HF1(502, nhXY);
	HF1(503, nhZ);
	HF1(504, chisqrXY);
	HF1(505, chisqrZ);
	HF1(506, x0);    HF1(507, y0);
	HF1(508, u0);    HF1(509, v0);
	HF2(510, x0, y0);HF1(511, u0, v0);
	*/	
	double totalDE = tp->TotalDEHiGain();
	double maxDE = tp->MaxDEHiGain();
	double pathlength = tp->GetTotalPathLength();
	double normalizedTotalDE = tp->NormalizedTotalDEHiGain();
	double normalizedMaxDE = tp->NormalizedMaxDEHiGain();
	/*	
	HF1(512, totalDE);
	HF1(513, maxDE);
	HF1(514, pathlength);
	HF1(515, normalizedTotalDE);
	HF1(516, normalizedMaxDE);
	HF2(517, pathlength, totalDE);
	HF2(518, pathlength, maxDE);
	HF2(519, pathlength, normalizedTotalDE);
	HF2(520, pathlength, normalizedMaxDE);
	*/	
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
	
	
	if (CFTPart->GetMass() > 0.9 && CFTProtonCont.size() == 0)
	  CFTProtonCont.push_back(CFTPart);
	else if (CFTPart->GetMass() > 0.0 && CFTPart->GetMass() < 0.2 && CFTPionCont.size() == 0)
	  CFTPionCont.push_back(CFTPart);
      }
    }
  }



  int nP_CFT = CFTProtonCont.size();
  int nPi_CFT = CFTPionCont.size();
  
  event.nP_CFT = nP_CFT;
  event.nPi_CFT = nPi_CFT;

  for (int it=0; it<nP_CFT; it++) {
    int nhBGO = CFTProtonCont[it]->NHitBGO();
    for (int j=0; j<nhBGO; j++) {
      CFTProtonCont[it]->GetBGOHit(j)->SetBGO_Pid(10);
    }
  }


  for (int it=0; it<nPi_CFT; it++) {
    int nh = CFTPionCont[it]->NHitBGO();
    int minSeg, maxSeg;
    for (int j=0; j<nh; j++) {
      CFTPionCont[it]->GetBGOHit(j)->SetBGO_Pid(1);
      if (j==0) {
	minSeg = CFTPionCont[it]->GetBGOHit(j)->SegmentId();
	maxSeg = CFTPionCont[it]->GetBGOHit(j)->SegmentId();
      } else {
	int seg=CFTPionCont[it]->GetBGOHit(j)->SegmentId();
	if (minSeg == 0 && seg == 23)
	  minSeg = seg;
	else if (seg < minSeg)
	  minSeg = seg;

	if (maxSeg == 23 && seg ==0)
	  maxSeg = seg;
	else if (seg > maxSeg)
	  maxSeg = seg;
	  
      }
    }

    int clustNum=0;
    for (int j=1; j<=23; j++) {
      int TargetSeg = minSeg - j;
      if (TargetSeg<0)
	TargetSeg += 24;

      bool hitFlag = false;
      for (int n=0; n<nhitBGO; n++) {
	Hodo2Hit* hit = hodoAna->GetHitBGO(n);
	int seg = hit->SegmentId();
	double time = hit->CMeanTime();
	double dE   = hit->DeltaE();
	if ((time>0 && time<5 && dE > 0.2) && (seg == TargetSeg))
	  hitFlag = true;
      }
      if (hitFlag)
	clustNum++;
      else
	break;
    }

    for (int j=1; j<=23; j++) {
      int TargetSeg = maxSeg + j;
      if (TargetSeg>=24)
	TargetSeg -= 24;

      bool hitFlag = false;
      for (int n=0; n<nhitBGO; n++) {
	Hodo2Hit* hit = hodoAna->GetHitBGO(n);
	int seg = hit->SegmentId();
	double time = hit->CMeanTime();
	double dE   = hit->DeltaE();
	if ((time>0 && time<5 && dE > 0.2) && (seg == TargetSeg))
	  hitFlag = true;
      }
      if (hitFlag)
	clustNum++;
      else
	break;
    }

    if (it<MaxHits2) {
      event.BGO_ClustNum_Pi[it] = CFTPionCont[it]->NHitBGO();
      event.BGO_ClustNumNoTrack_Pi[it] = clustNum;
    }
  }

  event.nhBGO_NoTrack=0;
  event.MaxBGOEdep_NoTrack=0.;
  event.TotalBGOEdep_NoTrack=0.;
  for (int i=0; i<nhitBGO; i++) {
    Hodo2Hit* hit = hodoAna->GetHitBGO(i);
    int pid = hit->GetBGO_Pid();
    double time = hit->CMeanTime();
    double dE   = hit->DeltaE();
    //std::cout << "pid : " << pid << ", time : " << time << ", dE : " << dE << std::endl;
    if (pid == 0 && time>0 && time <5 && dE >= 0.2) {
      event.nhBGO_NoTrack++;
      if (dE >= event.MaxBGOEdep_NoTrack)
	event.MaxBGOEdep_NoTrack = dE;
      event.TotalBGOEdep_NoTrack += dE;
    }
  }



  HF1(1000, EffStudyEvent);
  if (nP_CFT == 1)
    HF1(1001, EffStudyEvent);
  if (nPi_CFT == 1)
    HF1(1002, EffStudyEvent);
  if (nP_CFT == 1 && nPi_CFT == 1) {
    HF1(1003, EffStudyEvent);
    if (event.BGO_ClustNumNoTrack_Pi[0]==0 && event.MaxBGOEdep_NoTrack<10)
      HF1(1004, EffStudyEvent);
  }

  if (nP_CFT == 2) {
    HF1(1005, EffStudyEvent);
    if (nPi_CFT == 0)
      HF1(1006, EffStudyEvent);
  } 


  for (int nt=0; nt<nPi_CFT; nt++) {
    CFTParticle *PiCFT = CFTPionCont[nt];
    CFTLocalTrack *tp = PiCFT->GetTrack();
    int NHit = tp->GetNHit();
    for (int i=0; i<NHit; i++) {
      CFTFiberCluster *fcl = tp->GetHit(i);
      int clSize = fcl->ClusterSize();
      for (int j=0; j<clSize; j++) {
	CFTFLHit *flhit = fcl->GetHit(j);
	int seg = flhit->PairId();
	int l = flhit->GetTrackingLayer() -61;
	//std::cout << l << ", " << seg << std::endl;

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
      for (int j=0; j<clSize; j++) {
	CFTFLHit *flhit = fcl->GetHit(j);
	int seg = flhit->PairId();
	int l = flhit->GetTrackingLayer()-61;
	//std::cout << l << ", " << seg << std::endl;

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

  for (int nt=0; nt<nP_CFT; nt++) {
    CFTParticle *ProtonCFT = CFTProtonCont[nt];
    CFTLocalTrack *tp = ProtonCFT->GetTrack();
    int NHit = tp->GetNHit();
    for (int i=0; i<NHit; i++) {
      CFTFiberCluster *fcl = tp->GetHit(i);
      int clSize = fcl->ClusterSize();
      for (int j=0; j<clSize; j++) {
	CFTFLHit *flhit = fcl->GetHit(j);
	int seg = flhit->PairId();
	int l = flhit->GetTrackingLayer() -61;
	//std::cout << l << ", " << seg << std::endl;

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
      for (int j=0; j<clSize; j++) {
	CFTFLHit *flhit = fcl->GetHit(j);
	int seg = flhit->PairId();
	int l = flhit->GetTrackingLayer()-61;
	//std::cout << l << ", " << seg << std::endl;

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
  
  if (FlagEvDisp) {
    const EvDispCFT & evDisp = EvDispCFT::GetInstance();
    //if (1) {
    if (ntCFT_2nd==1) {
    //if (event.scatFlag==-1&&event.NNscatFlag==-1&&event.PiNscatFlag==-1&&nP_CFT==1&&nPi_CFT==1) {
    //if (nP_CFT==1&&nPi_CFT==1) {
      evDisp.UpdateCanvas();
      evDisp.get_command();
    }
    evDisp.EndOfEvent();
  }


  tree->Fill();


  for_each(CFTPartCont.begin(), CFTPartCont.end(), DeleteObject());

  delete DCAna;
  delete hodoAna;
  delete rawData;

  events++;

  return true;
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

  event.ntCFT = 0;
  for (int i=0; i<MaxHits2; i++) {
    event.BGO_Edep[i]=-999.;
    event.TotalEdep[i]=-999.;
    event.CFT_TotalEdep[i]=-999.;
    event.CFT_NormTotalEdep[i]=-999.;
    event.PiV_Edep[i]=-999.;
    event.theta[i]=-999.;

    event.BGO_ClustNum_Pi[i]=0;
    event.BGO_ClustNumNoTrack_Pi[i]=0;
  }

  for (int i=0; i<NumOfPlaneCFT; i++) {
    for (int j=0; j<MaxHits2; j++) {
      event.CFT_MaxFiberEdep [i][j] = -999.;
    }
  }

  event.nP_CFT = 0;
  event.nPi_CFT = 0;

  event.nhBGO_NoTrack = 0;
  event.MaxBGOEdep_NoTrack = -999.;
  event.TotalBGOEdep_NoTrack = -999.;
}


void DefineHistograms( const char *filename )
{
  new TFile( filename, "recreate" );

  char buf[100];
  for (int l=0; l<NumOfLayersCFT; l++) {
    int hid = 100+10*l;
    sprintf(buf, "N hits layer %d", l);
    HB1( hid+1, buf,  20, 0, 20 );
    sprintf(buf, "Hit pattern layer %d", l);
    HB1( hid+2, buf, NumOfSegCFT[l], 0, NumOfSegCFT[l] );
    sprintf(buf, "Time layer %d", l);
    HB1( hid+3, buf, 1000, 0, 10 );
    sprintf(buf, "dE layer %d", l);
    HB1( hid+4, buf, 1000, 0, 10 );
  }

  for (int l=0; l<NumOfLayersCFT; l++) {
    int hid = 200+10*l;
    sprintf(buf, "N Cluster layer %d", l);
    HB1( hid+1, buf,  20, 0, 20 );
    sprintf(buf, "Cluster Hit pattern layer %d", l);
    HB1( hid+2, buf, NumOfSegCFT[l], 0, NumOfSegCFT[l] );
    sprintf(buf, "Cluster size layer %d", l);
    HB1( hid+3, buf,  20, 0, 20 );
    sprintf(buf, "Time (Cluster) layer %d", l);
    HB1( hid+4, buf, 1000, 0, 10 );
    sprintf(buf, "Total dE (Cluster) layer %d", l);
    HB1( hid+5, buf, 1000, 0, 10 );
    sprintf(buf, "Max dE (Cluster) layer %d", l);
    HB1( hid+6, buf, 1000, 0, 10 );
  }

  HB1(301, "N Hit BGO", 10, 0, 10);
  HB1(302, "Hit pattern BGO", NumOfSegBGO, 0, NumOfSegBGO);
  HB1(303, "Time BGO", 1000, 0, 10);
  HB1(304, "dE BGO", 1000, 0, 200);
  HB2(305, "dE % theta", 180, 0, 90, 200, 0, 200);

  HB1(311, "N Cluster BGO", 10, 0, 10);
  HB1(312, "Hit pattern BGO", NumOfSegBGO, 0, NumOfSegBGO);
  HB1(313, "Cluster size BGO", 10, 0, 10);
  HB1(314, "Time BGO", 1000, 0, 10);
  HB1(315, "dE BGO", 1000, 0, 200);
  HB2(316, "dE % theta", 180, 0, 90, 200, 0, 200);

  HB1(321, "N Hit PiV", 10, 0, 10);
  HB1(322, "Hit pattern PiV", NumOfSegPiV, 0, NumOfSegPiV);
  HB1(323, "Time PiV", 1000, 0, 10);
  HB1(324, "dE PiV", 1000, 0, 200);
  HB2(325, "dE % theta", 180, 0, 90, 200, 0, 200);

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
  HB1(421, "Delta Pos at BGO", 200, -50, 50);
  HB2(422, "Total DE%dE BGO",150, 0, 150,  100, 0, 50);
  HB2(423, "Total Maximum DE%dE BGO", 150, 0, 150, 100, 0, 50);
  HB2(424, "Normalized Total DE%dE BGO",150, 0, 150,  100, 0, 10);
  HB2(425, "Normalized Total Maximum DE%dE BGO", 150, 0, 150, 100, 0, 10);
  HB2(426, "Total DE%dE BGO+Fiber(total)",150, 0, 150,  100, 0, 50);
  HB2(427, "Total Maximum DE%dE BGO+Fiber(max)", 150, 0, 150, 100, 0, 50);
  HB2(428, "Normalized Total DE%dE BGO+Fiber(total)",150, 0, 150,  100, 0, 5);
  HB2(429, "Normalized Total Maximum DE%dE BGO+Fiber(max)", 150, 0, 150, 100, 0, 5);
  HB1(430, "Delta Pos at PiV", 200, -50, 50);
  HB1(431, "Total DE (PiV)", 100, 0, 20);
  HB2(432, "Total DE%dE BGO (PiV Hit)",150, 0, 150,  100, 0, 50);
  HB2(433, "Total Maximum DE%dE BGO (PiV Hit)", 150, 0, 150, 100, 0, 50);
  HB2(434, "Normalized Total DE%dE BGO (PiV Hit)",150, 0, 150,  100, 0, 10);
  HB2(435, "Normalized Total Maximum DE%dE BGO (PiV Hit)", 150, 0, 150, 100, 0, 10);
  HB2(436, "Total DE%dE BGO (No PiV Hit)",150, 0, 150,  100, 0, 50);
  HB2(437, "Total Maximum DE%dE BGO (No PiV Hit)", 150, 0, 150, 100, 0, 50);
  HB2(438, "Normalized Total DE%dE BGO (No PiV Hit)",150, 0, 150,  100, 0, 10);
  HB2(439, "Normalized Total Maximum DE%dE BGO (No PiV Hit)", 150, 0, 150, 100, 0, 10);

  HB1(1000, "Generated Number", 500, 0, 500);
  HB1(1001, "Proton Detection Number", 500, 0, 500);
  HB1(1002, "Pion Detection Number", 500, 0, 500);
  HB1(1003, "Proton & Pion Detection Number", 500, 0, 500);
  HB1(1004, "Proton & Pion Detection Number (w/ Decay Cut)", 500, 0, 500);
  HB1(1005, "2 Proton Detection Number", 500, 0, 500);
  HB1(1006, "2 Proton and 0 Pion Detection Number", 500, 0, 500);


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

  tree->Branch("ntCFT",   &event.ntCFT,  "ntCFT/I");
  tree->Branch("BGO_Edep",   event.BGO_Edep,  "BGO_Edep[ntCFT]/D");
  tree->Branch("TotalEdep",   event.TotalEdep,  "TotalEdep[ntCFT]/D");
  tree->Branch("CFT_TotalEdep",   event.CFT_TotalEdep,  "CFT_TotalEdep[ntCFT]/D");
  tree->Branch("CFT_NormTotalEdep",   event.CFT_NormTotalEdep,  "CFT_NormTotalEdep[ntCFT]/D");
  tree->Branch("PiV_Edep",   event.PiV_Edep,  "PiV_Edep[ntCFT]/D");
  tree->Branch("theta",   event.theta,  "theta[ntCFT]/D");

  for (int i=0; i<NumOfPlaneCFT; i++) {
    char buf1[100], buf2[100];

    sprintf(buf1, "CFT_MaxFiberEdep%d", i);
    sprintf(buf2, "CFT_MaxFiberEdep%d[ntCFT]/D", i);
    tree->Branch(buf1,   event.CFT_MaxFiberEdep[i],  buf2);
  }


  tree->Branch("nP_CFT",   &event.nP_CFT,  "nP_CFT/I");
  tree->Branch("nPi_CFT",   &event.nPi_CFT,  "nPi_CFT/I");
}

bool ConfMan::InitializeParameterFiles(void)
{

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( DCGeomFileName_!="" )
    DCGeomManager_->Initialize(DCGeomFileName_);
  else
    DCGeomManager_->Initialize();


  return true;
}
