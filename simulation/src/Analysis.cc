/*
  Analysis.cc
  2007/4  K.Shirotori
*/

#include "Analysis.hh"
//#include "PrimaryInfo.hh"
#include "CFTFiberHit.hh"
#include "CrystalHit.hh"
#include "SKSChamberHit.hh"
#include "SKSCounterHit.hh"
#include "SKSVirtualPlaneHit.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "Randomize.hh"
#include "RootHelper.hh"

#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>

using namespace CLHEP;

FILE *fpEff;

Analysis::Analysis( G4String file )
  : filename_(file), fActive_(true), fTriggered(false), fPrimary_(false),
    DataFile_()
{
  ConfMan *confMan = ConfMan::GetConfManager();
  ReactionMode_ = confMan->ReactionMode();

  DefineHistograms();
}

Analysis::~Analysis()
{
  SaveFile();
}
//static Event event;

void Analysis::BeginOfRun( const G4Run *aRun )
{
  G4cout<<"*****************************"<< G4endl;  
  G4cout << "Analysis::BeginOfRun::RunNo = " << aRun->GetRunID() << G4endl;
  
  trigNum=0;
}

void Analysis::EndOfRun( const G4Run *aRun )
{
  G4cout<<"*****************************"<< G4endl;  
  G4cout << "Analysis::EndOfRun::RunNo = " << aRun->GetRunID() << ", Total # of SksMinus triggers =" 
	 << std::setw(8) << trigNum << G4endl;
}

void Analysis::BeginOfPrimaryAction()
{
  event.thetaMeson          = -999.9;
  event.phiMeson            = -999.9;
  event.thetaMesonCM        = -999.9;
  event.phiMesonCM          = -999.9;
  event.thetaScatHypCM   = -999.9;
  event.phiScatHypCM     = -999.9;
  event.thetaScatPLab    = -999.9;
  event.thetaScatHypLab  = -999.9;

  event.beammom     = -999.9;

  for (int i=0; i<3; i++) {
    event.momVectorScatMeson[i]      = -999.9;
    event.momVectorHypBeam[i]        = -999.9;
    event.momVectorHypScat[i]        = -999.9;
    event.momVectorProtonScat[i]     = -999.9;
    event.momVectorDecayPi[i]        = -999.9;
    event.momVectorDecayNucleon[i]   = -999.9;
    event.primaryVertex[i]           = -999.9;
    event.scatPos[i]                 = -999.9;
    event.NNscatPos[i]                 = -999.9;
    event.PiNscatPos[i]                 = -999.9;
    event.decayPos[i]                = -999.9;
    event.momVectorPiMinusFromK0[i]  = -999.9;
    //event.fermiMom[i]                = -999.9;
    //event.protonStopPos[i]           = -999.9;
  }
  /*
  event.scatPartTpcHitNum = -1;
  event.scatPartCrystalHitNum = -1;
  event.scatPartSi1HitNum = -1;
  event.scatPartSi2HitNum = -1;
  event.scatPartSi3HitNum = -1;
  */


  event.decayEventFlag = -1;
  event.decayFlag = -1;
  event.scatFlag = -1;
  event.scatTarget = -1;
  event.NNscatFlag = -1;
  event.NNscatTarget = -1;
  event.PiNscatFlag = -1;
  event.PiNscatTarget = -1;

  event.fLengthInH = -999.9;
  event.fLengthInD = -999.9;
  event.fLengthTotal = -999.9;

  event.nDaughter = 0;
  for (int i=0; i<MaxHit; i++) {
    event.daughterPID[i] = -1;
    event.daughterMom[i] = G4ThreeVector(-999, -999, -999);
  }

}

void Analysis::BeginOfEvent( )
{
  fTriggered = false;
}

void Analysis::SetScatMesonMomentum(G4ThreeVector mesonMom)
{
  event.momVectorScatMeson[0] = mesonMom.x();
  event.momVectorScatMeson[1] = mesonMom.y();
  event.momVectorScatMeson[2] = mesonMom.z();

  //event.momMeson = mesonMom.mag();

  return;
}

void Analysis::GetScatMesonMomentum(G4double *mesonMom)
{
  mesonMom[0] = event.momVectorScatMeson[0];
  mesonMom[1] = event.momVectorScatMeson[1];
  mesonMom[2] = event.momVectorScatMeson[2];

  return;
}

void Analysis::SetHypBeamMomentum(G4ThreeVector hypBeamMom)
{
  event.momVectorHypBeam[0] = hypBeamMom.x();
  event.momVectorHypBeam[1] = hypBeamMom.y();
  event.momVectorHypBeam[2] = hypBeamMom.z();

  return;
}

void Analysis::GetHypBeamMomentum(G4double *hypBeamMom)
{
  hypBeamMom[0] = event.momVectorHypBeam[0];
  hypBeamMom[1] = event.momVectorHypBeam[1];
  hypBeamMom[2] = event.momVectorHypBeam[2];

  return;
}

void Analysis::SetScatHypMomentum(G4ThreeVector hypScatMom)
{
  event.momVectorHypScat[0] = hypScatMom.x();
  event.momVectorHypScat[1] = hypScatMom.y();
  event.momVectorHypScat[2] = hypScatMom.z();

  return;
}

void Analysis::GetScatHypMomentum(G4double *hypScatMom)
{
  hypScatMom[0] = event.momVectorHypScat[0];
  hypScatMom[1] = event.momVectorHypScat[1];
  hypScatMom[2] = event.momVectorHypScat[2];

  return;
}

void Analysis::SetScatProtonMomentum(G4ThreeVector protonScatMom)
{
  event.momVectorProtonScat[0] = protonScatMom.x();
  event.momVectorProtonScat[1] = protonScatMom.y();
  event.momVectorProtonScat[2] = protonScatMom.z();

  return;
}

void Analysis::GetScatProtonMomentum(G4double *protonScatMom)
{
  protonScatMom[0] = event.momVectorProtonScat[0];
  protonScatMom[1] = event.momVectorProtonScat[1];
  protonScatMom[2] = event.momVectorProtonScat[2];

  return;
}

void Analysis::SetDecayPiMomentum(G4ThreeVector piMom)
{
  event.momVectorDecayPi[0] = piMom.x();
  event.momVectorDecayPi[1] = piMom.y();
  event.momVectorDecayPi[2] = piMom.z();

  return;
}

void Analysis::GetDecayPiMomentum(G4double *piMom)
{
  piMom[0] = event.momVectorDecayPi[0];
  piMom[1] = event.momVectorDecayPi[1];
  piMom[2] = event.momVectorDecayPi[2];

  return;
}

void Analysis::SetDecayNucleonMomentum(G4ThreeVector nucleonMom)
{
  event.momVectorDecayNucleon[0] = nucleonMom.x();
  event.momVectorDecayNucleon[1] = nucleonMom.y();
  event.momVectorDecayNucleon[2] = nucleonMom.z();

  return;
}

void Analysis::GetDecayNucleonMomentum(G4double *nucleonMom)
{
  nucleonMom[0] = event.momVectorDecayNucleon[0];
  nucleonMom[1] = event.momVectorDecayNucleon[1];
  nucleonMom[2] = event.momVectorDecayNucleon[2];

  return;
}

void Analysis::SetDaughterMomentum(G4int type, G4ThreeVector Mom)
{
  event.daughterPID[event.nDaughter] = type;
  event.daughterMom[event.nDaughter] = Mom;

  event.nDaughter++;

  return;
}

int Analysis::GetDaughterPID(G4int i)
{
  if (i>=0 && i<event.nDaughter)
    return event.daughterPID[i];

  return -1;
}

G4ThreeVector Analysis::GetDaughterMomentum(G4int i)
{
  if (i>=0 && i<event.nDaughter)
    return event.daughterMom[i];

  return G4ThreeVector(-999, -999, -999);
}


void Analysis::SetThetaMeson(double theta)
{
  event.thetaMeson = theta;

  return;
}

G4double Analysis::GetThetaMeson()
{
  return event.thetaMeson;
}

void Analysis::SetPhiMeson(double phi)
{
  event.phiMeson = phi;

  return;
}

void Analysis::SetThetaMesonCM(double theta)
{
  event.thetaMesonCM = theta;

  return;
}

void Analysis::SetPhiMesonCM(double phi)
{
  event.phiMesonCM = phi;

  return;
}

void Analysis::SetThetaScatHypCM(double theta)
{
  event.thetaScatHypCM = theta;

  return;
}

double Analysis::GetThetaScatHypCM()
{
  return event.thetaScatHypCM;
}

void Analysis::SetPhiScatHypCM(double phi)
{
  event.phiScatHypCM = phi;

  return;
}

void Analysis::SetThetaScatPLab(double theta)
{
  event.thetaScatPLab = theta;

  return;
}

double Analysis::GetThetaScatPLab()
{
  return event.thetaScatPLab;

}

void Analysis::SetThetaScatHypLab(double theta)
{
  event.thetaScatHypLab = theta;

  return;
}

double Analysis::GetThetaScatHypLab()
{
  return event.thetaScatHypLab;

}

void Analysis::SetBeamMomentum(G4double mom)
{
  event.beammom = mom;
}

double Analysis::GetBeamMomentum()
{
  return event.beammom;
}

void Analysis::SetPrimaryVertex(G4ThreeVector pos)
{
  event.primaryVertex[0] = pos.x();
  event.primaryVertex[1] = pos.y();
  event.primaryVertex[2] = pos.z();

  return;
}

G4ThreeVector Analysis::GetPrimaryVertex()
{
  return G4ThreeVector(event.primaryVertex[0], event.primaryVertex[1], event.primaryVertex[2]);
}

void Analysis::SetScatPos(G4ThreeVector pos)
{
  event.scatPos[0] = pos.x();
  event.scatPos[1] = pos.y();
  event.scatPos[2] = pos.z();

  return;
}

G4ThreeVector Analysis::GetScatPos()
{
  return G4ThreeVector(event.scatPos[0], event.scatPos[1], event.scatPos[2]);
}

void Analysis::SetNNScatPos(G4ThreeVector pos)
{
  event.NNscatPos[0] = pos.x();
  event.NNscatPos[1] = pos.y();
  event.NNscatPos[2] = pos.z();

  return;
}

void Analysis::SetPiNScatPos(G4ThreeVector pos)
{
  event.PiNscatPos[0] = pos.x();
  event.PiNscatPos[1] = pos.y();
  event.PiNscatPos[2] = pos.z();

  return;
}

void Analysis::SetDecayPos(G4ThreeVector pos)
{
  event.decayPos[0] = pos.x();
  event.decayPos[1] = pos.y();
  event.decayPos[2] = pos.z();

  return;
}

G4ThreeVector Analysis::GetDecayPos()
{
  return G4ThreeVector(event.decayPos[0], event.decayPos[1], event.decayPos[2]);
}

void Analysis::SetDecayEventFlag()
{
  event.decayEventFlag = 1;
  return;
}

void Analysis::SetDecayFlag()
{
  event.decayFlag = 1;
  return;
}

void Analysis::SetDecayFlag(int decayMode)
{
  event.decayFlag = decayMode;
  return;
}

void Analysis::SetScatFlag()
{
  event.scatFlag = 1;
  return;
}

void Analysis::SetScatFlag(int reactMode)
{
  event.scatFlag = reactMode;
  return;
}

void Analysis::SetScatTarget(int tgt)
{
  event.scatTarget = tgt;
  return;
}


void Analysis::SetNNScatFlag()
{
  event.NNscatFlag = 1;
  return;
}

void Analysis::SetNNScatTarget(int tgt)
{
  event.NNscatTarget = tgt;
  return;
}

void Analysis::SetPiNScatFlag()
{
  event.PiNscatFlag = 1;
  return;
}

void Analysis::SetPiNScatTarget(int tgt)
{
  event.PiNscatTarget = tgt;
  return;
}

void Analysis::SetFlightLengthInTarget(double dxInH, double dxInD)
{
  event.fLengthInH = dxInH;
  event.fLengthInD = dxInD;
  event.fLengthTotal = event.fLengthInH+event.fLengthInD;

  return;
}

void Analysis::SetPiMinusMomentumFromK0(G4ThreeVector piMom)
{
  event.momVectorPiMinusFromK0[0] = piMom.x();
  event.momVectorPiMinusFromK0[1] = piMom.y();
  event.momVectorPiMinusFromK0[2] = piMom.z();

  return;
}


void Analysis::EndOfEvent( const G4Event *anEvent )
{
  //Tree
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  //Analysis Part for each Event
#if 0
  PrintHitsInformation( anEvent, G4cout );
#endif

  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhFiber=0, nhCrystal=0,  nhCounter=0, nhVirtualPlane=0;

  CFTFiberHitsCollection   *FiberHC;
  CrystalHitsCollection    *CrystalHC;
  SKSChamberHitsCollection *ChamberHC;
  SKSCounterHitsCollection *CounterHC;
  SKSVirtualPlaneHitsCollection *VirtualPlaneHC;


  G4int colIdFiber   = SDMan->GetCollectionID( "CFTFiberCollection" );
  G4int colIdCrystal   = SDMan->GetCollectionID( "CrystalCollection" );
  G4int colIdChamber   = SDMan->GetCollectionID( "chamberSD" );
  G4int colIdCounter   = SDMan->GetCollectionID( "counterSD" );
  G4int colIdVirtualPlane   = SDMan->GetCollectionID( "virtualPlaneSD" );

  FiberHC   = dynamic_cast<CFTFiberHitsCollection *>(HCE->GetHC(colIdFiber));
  CrystalHC = dynamic_cast<CrystalHitsCollection *>(HCE->GetHC(colIdCrystal));
  ChamberHC = dynamic_cast<SKSChamberHitsCollection *>(HCE->GetHC(colIdChamber));
  CounterHC = dynamic_cast<SKSCounterHitsCollection *>(HCE->GetHC(colIdCounter));
  VirtualPlaneHC = dynamic_cast<SKSVirtualPlaneHitsCollection *>(HCE->GetHC(colIdVirtualPlane));


  if( FiberHC )  nhFiber  = FiberHC->entries();
  if( CrystalHC )  nhCrystal  = CrystalHC->entries();
  if( CounterHC )  nhCounter  = CounterHC->entries();
  if( VirtualPlaneHC )  nhVirtualPlane  = VirtualPlaneHC->entries();

  int nhFiberForScat[NumOfPlaneCFT];
  for (int i=0; i<NumOfPlaneCFT; i++)
    nhFiberForScat[i] = 0;

  int nhFiberForDecayProton[NumOfPlaneCFT];
  for (int i=0; i<NumOfPlaneCFT; i++)
    nhFiberForDecayProton[i] = 0;

  if (FiberHC) {
    for (G4int i=0; i<nhFiber; i++) {
      int planeId = (*FiberHC)[i]->GetLayerID();
      int segmentId = (*FiberHC)[i]->GetSegmentID();
      double de = (*FiberHC)[i]->GetEdepAbs();
      double time = (*FiberHC)[i]->GetTime();
      G4String name = (*FiberHC)[i]->GetParticleName();

      /*
      G4cout << "planeId : " << planeId << ", segmentId : " << segmentId 
	     << ", de : " << de << ", time : " << time << ", name : "
	     << name << G4endl;
      */
      int particleId = -1;
      if (name == "pi+" || name == "spi+")
	particleId = 1;
      else if (name == "pi-" || name == "spi-")
	particleId = 2;
      else if (name == "kaon+" || name == "skaon+")
	particleId = 3;
      else if (name == "kaon-")
	particleId = 4;
      else if (name == "e-")
	particleId = 5;
      else if (name == "e+")
	particleId = 6;
      else if (name == "mu-")
	particleId = 7;
      else if (name == "mu+")
	particleId = 8;
      else if (name == "gamma")
	particleId = 9;
      else if (name == "proton")
	particleId = 10;
      else if (name == "deuteron")
	particleId = 11;
      else if (name == "sigma+")
	particleId = 12;
      else if (name == "sigma-")
	particleId = 13;
      else if (name == "proton_d") 
	particleId = 14;
      else {
	;//std::cout << "Particle Name : " << name << std::endl;
      }

      if (time/ns < 200) {
	if (event.FiberHits[planeId]<MaxHit) {
	  event.FiberSeg[planeId][event.FiberHits[planeId]] = segmentId;
	  event.FiberTime[planeId][event.FiberHits[planeId]] = time;
	  event.FiberEdep[planeId][event.FiberHits[planeId]] = de;
	  event.FiberPID[planeId][event.FiberHits[planeId]]  = particleId;
	  event.FiberHits[planeId]++;
	}
	if (particleId == 10 || particleId == 11)
	  nhFiberForScat[planeId]++;
	else if (particleId == 14)
	  nhFiberForDecayProton[planeId]++;
      }
    }
  }

  int nhCrystalForScat=0;
  int nhCrystalForDecayProton=0;
  if (CrystalHC) {
    for (G4int i=0; i<nhCrystal; i++) {
      //int planeId = (*CrystalHC)[i]->GetLayerID();
      int segmentId = (*CrystalHC)[i]->GetChannelID();
      double de = (*CrystalHC)[i]->GetEdepAbs();
      double time = (*CrystalHC)[i]->GetTime();
      G4String name = (*CrystalHC)[i]->GetParticleName();
      /*
      G4cout << "BGO segmentId : " << segmentId 
	     << ", de : " << de << ", time : " << time << ", name : "
	     << name << G4endl;
      */

      int particleId = -1;
      if (name == "pi+")
	particleId = 1;
      else if (name == "pi-")
	particleId = 2;
      else if (name == "kaon+" || name == "skaon+")
	particleId = 3;
      else if (name == "kaon-")
	particleId = 4;
      else if (name == "e-")
	particleId = 5;
      else if (name == "e+")
	particleId = 6;
      else if (name == "mu-")
	particleId = 7;
      else if (name == "mu+")
	particleId = 8;
      else if (name == "gamma")
	particleId = 9;
      else if (name == "proton")
	particleId = 10;
      else if (name == "deuteron")
	particleId = 11;
      else if (name == "sigma+")
	particleId = 12;
      else if (name == "sigma-")
	particleId = 13;
      else if (name == "proton_d")
	particleId = 14;
      else {
	;//std::cout << "Particle Name : " << name << std::endl;
      }
      if (time/ns < 200) {
	if (event.CrystalHits<MaxHit) {
	  event.CrystalSeg[event.CrystalHits] = segmentId;
	  event.CrystalTime[event.CrystalHits] = time;
	  event.CrystalEdep[event.CrystalHits] = de;
	  event.CrystalPID[event.CrystalHits] = particleId;
	  event.CrystalHits++;
	}
	if (particleId == 10 || particleId == 11)
	  nhCrystalForScat++;
	if (particleId == 14)
	  nhCrystalForDecayProton++;
      }
    }
  }

  // Chamber Hit
  if (HCE) ChamberHC = (SKSChamberHitsCollection*)(HCE->GetHC(colIdChamber));
  if (ChamberHC) {
    int n_hit = ChamberHC->entries();
    // G4cout << "Chamber hit : " << n_hit << G4endl;
    for (G4int i=0;i<n_hit;i++){
#if 1
      int planeID = (*ChamberHC)[i]->GetDetectorID();
      int hitWire = (*ChamberHC)[i]->GetHitWire();
      double dlength = (*ChamberHC)[i]->GetDRlength();
      bool sensitiveHit = (*ChamberHC)[i]->GetSensitiveHit();
      //DcPassMulti[planeID]++; /* not take into account the efficiency of chamber*/

      // G4cout << "planeId : " << planeID << " hitWire : "
      // << hitWire << " dlength : " << dlength << G4endl;

      // if (planeID == 1) event.HitPat1 = hitWire;
      // else if (planeID == 2) event.HitPat2 = hitWire;
      // else if (planeID == 3) event.HitPat3 = hitWire;
      // else if (planeID == 4) event.HitPat4 = hitWire;
      // else if (planeID == 5) event.HitPat5 = hitWire;
      // else if (planeID == 6) event.HitPat6 = hitWire;
      // else if (planeID == 7) event.HitPat7 = hitWire;
      // else if (planeID == 8) event.HitPat8 = hitWire;
      // else if (planeID == 9) event.HitPat9 = hitWire;
      // else if (planeID == 10) event.HitPat10 = hitWire;
      // else if (planeID == 31) event.HitPat31 = hitWire;
      // else if (planeID == 32) event.HitPat32 = hitWire;
      // else if (planeID == 33) event.HitPat33 = hitWire;
      // else if (planeID == 34) event.HitPat34 = hitWire;
      // else if (planeID == 35) event.HitPat35 = hitWire;
      // else if (planeID == 36) event.HitPat36 = hitWire;
      // else if (planeID == 37) event.HitPat37 = hitWire;
      // else if (planeID == 38) event.HitPat38 = hitWire;
      // else if (planeID == 39) event.HitPat39 = hitWire;
      // else if (planeID == 40) event.HitPat40 = hitWire;

      // if (DcMulti[planeID] < DCMultiMax && sensitiveHit) {
      // 	/* take into account the efficiency of chamber*/
      // 	Dc1ev[planeID][DcMulti[planeID]].hitWire = hitWire;
      // 	Dc1ev[planeID][DcMulti[planeID]].dlength = dlength;
      // 	DcMulti[planeID]++; 
      // } else if (DcMulti[planeID] >= DCMultiMax){
      // 	G4cerr << "SKSEventAction::EndOfEventAction too many mulitiplicity SDC" << planeID << " Multi : "  << DcMulti[planeID] << G4endl;
      // }	else if (!sensitiveHit){
      // 	G4cerr << "SKSEventAction::EndOfEventAction inefficiency of SDC" << planeID  << G4endl;
      // }
#endif
    }
  }

  int MultiFTOF=0, MultiCH=0, MultiAC=0, MultiT0=0;
  // Counter Hit

  if (CounterHC) {
    for (G4int i=0;i<nhCounter;i++){
      int planeID = (*CounterHC)[i]->GetDetectorID();

      const DCGeomMan &geomMan = DCGeomMan::GetInstance();
      const int FTofId = geomMan.GetDetectorId("FTOF");
      const int T0Id = geomMan.GetDetectorId("T0");
      const int ACId = geomMan.GetDetectorId("AC");
      //const int CHId = geomMan.GetDetectorId("CH");
      if (planeID == FTofId) {
	/*
	G4cout << "FTOF : segment " << (*CounterHC)[i]->GetChannelID() 
	       << " time : " <<  (*CounterHC)[i]->GetTime()
	       << " de : "   <<  (*CounterHC)[i]->GetEdepAbs() << G4endl;
	*/
	event.HitSegFTOF[MultiFTOF] = (*CounterHC)[i]->GetChannelID() ;
	event.lcCherenkov[MultiFTOF] = (*CounterHC)[i]->GetlcCherenkov();
	MultiFTOF++;
      } else if (planeID == T0Id) {
	MultiT0++;
      } else if (planeID == ACId) {
	MultiAC++;
      }  else if (planeID == 71) { // BGO VP
	double pos[3];
	(*CounterHC)[i]->GetPos(pos);
	const DCGeomMan &geomMan = DCGeomMan::GetInstance();
	const int lnum = geomMan.GetDetectorId("CFT_U1");
	int seg =  (*CounterHC)[i]->GetChannelID() ;
	G4ThreeVector gloPos = geomMan.GetGlobalPosition(lnum);
	pos[0] = pos[0] - gloPos.x();
	pos[1] = pos[1] - gloPos.y();
	pos[2] = pos[2] - gloPos.z();

	event.xPosBGOVP[seg] = pos[0];
	event.yPosBGOVP[seg] = pos[1];
	event.zPosBGOVP[seg] = pos[2];

	/*
	G4cout << "BGO_VP : segment " << (*CounterHC)[i]->GetChannelID() 
	       << " time : " <<  (*CounterHC)[i]->GetTime()
	       << " (x, y, z) = ( "   <<  pos[0] << ", " 
	       << pos[1] << ", " << pos[2] << ")" << G4endl;
	*/	
      }
    }
  }


  event.nHitCH = MultiCH;
  event.nHitFTOF = MultiFTOF;
  event.nHitT0 = MultiT0;
  event.nHitAC = MultiAC;

  if (VirtualPlaneHC) {
    for (G4int i=0;i<nhVirtualPlane;i++){
      int planeID = (*VirtualPlaneHC)[i]->GetDetectorID();

      const DCGeomMan &geomMan = DCGeomMan::GetInstance();
      const int VP1Id = geomMan.GetDetectorId("VP1");
      const int VP2Id = geomMan.GetDetectorId("VP2");
      const int VP3Id = geomMan.GetDetectorId("VP3");
      const int VP4Id = geomMan.GetDetectorId("VP4");
      const int VP5Id = geomMan.GetDetectorId("VP5");
      if (planeID == VP1Id) {
	double x, y;
	(*VirtualPlaneHC)[i]->GetHitPosition(&x, &y);
	if (event.nHitVP1<MaxHit) {
	  event.xPosVP1[event.nHitVP1] = x;
	  event.yPosVP1[event.nHitVP1] = y;
	  event.nHitVP1++;
	} 
      } else if (planeID == VP2Id) {
	double x, y;
	(*VirtualPlaneHC)[i]->GetHitPosition(&x, &y);
	if (event.nHitVP2<MaxHit) {
	  event.xPosVP2[event.nHitVP2] = x;
	  event.yPosVP2[event.nHitVP2] = y;
	  event.nHitVP2++;
	} 
      } else if (planeID == VP3Id) {
	double x, y;
	(*VirtualPlaneHC)[i]->GetHitPosition(&x, &y);
	if (event.nHitVP3<MaxHit) {
	  event.xPosVP3[event.nHitVP3] = x;
	  event.yPosVP3[event.nHitVP3] = y;
	  event.nHitVP3++;
	} 
      }else if (planeID == VP4Id) {
	double x, y;
	(*VirtualPlaneHC)[i]->GetHitPosition(&x, &y);
	if (event.nHitVP4<MaxHit) {
	  event.xPosVP4[event.nHitVP4] = x;
	  event.yPosVP4[event.nHitVP4] = y;
	  event.nHitVP4++;
	} 
      }else if (planeID == VP5Id) {
	double x, y;
	(*VirtualPlaneHC)[i]->GetHitPosition(&x, &y);
	if (event.nHitVP5<MaxHit) {
	  event.xPosVP5[event.nHitVP5] = x;
	  event.yPosVP5[event.nHitVP5] = y;
	  event.nHitVP5++;
	} 
      }
    }
  }


  if (ReactionMode_ == 4 || ReactionMode_ == 5 || ReactionMode_ == 6  /*||  ReactionMode_ == 7  || ReactionMode_ == 8 */ || ReactionMode_ == 9 || ReactionMode_ == 10) {
    if (MultiFTOF>0 && MultiCH>0)
      fTriggered = true;
  } else if (ReactionMode_ == 14 || ReactionMode_ == 15 || ReactionMode_ == 16  ||  ReactionMode_ == 17  
	     || ReactionMode_ == 18  || ReactionMode_ == 19 || ReactionMode_ == 20
	     || ReactionMode_ == 21  || ReactionMode_ == 22 || ReactionMode_ == 32
	     || ReactionMode_ == 33  || ReactionMode_ == 34 || ReactionMode_ == 35
	     || ReactionMode_ == 36  || ReactionMode_ == 38 || ReactionMode_ == 39
	     ) {
    if (MultiFTOF>0)
      fTriggered = true;
  } else if ( ReactionMode_ == 37
	     ) {
    if (!fPrimary_ && MultiFTOF>0)
      fTriggered = true;
  } else if (ReactionMode_ == 95 || ReactionMode_ == 99 || ReactionMode_ == 100 || ReactionMode_ == 101) {
    if (event.scatFlag==1)
      fTriggered = true;
    /*
      static int num=0;
      static int nDetect=0;
      static int nDetectDecayProton=0;
      const int Nstudy = 1000;
      
      if (nhFiberForScat[1]>0 && nhFiberForScat[2]>0 && nhFiberForScat[3]>0 &&
      nhFiberForScat[4]>0 && nhCrystalForScat>0) {
      nDetect++;
      }
      if (nhFiberForDecayProton[1]>0 && nhFiberForDecayProton[2]>0 && nhFiberForDecayProton[3]>0 &&
      nhFiberForDecayProton[4]>0 && nhCrystalForDecayProton>0) {
      nDetectDecayProton++;
      }
      num++;
      if (num>=Nstudy) {
      std::cout << nDetect << "/" << num << ", " << nDetectDecayProton << "/" << num << std::endl;
      fprintf(fpEff, "%f %d %d\n", (double)nDetect/(double)num, nDetect, num);
      num=0;
      nDetect=0;
      nDetectDecayProton=0;
      }
    */
  } else if (ReactionMode_ == 96 || ReactionMode_ == 102) {
    fTriggered = false;
  }
  else {
    fTriggered = true;
  }
  

  if( DataFile_.is_open() && fTriggered){
    PrintHitsInformation( anEvent, DataFile_ );
  }

  tree->Fill();
}

void Analysis::SaveFile( void ) const
{
  if( fActive_ )
    gFile->Write();
}

void Analysis::Terminate( void ) const
{
  G4cout << "[Analysis] Terminate()" << G4endl;   
  if( fActive_ ){
    gFile->Write();
    gFile->Close();
  }
}

void Analysis::ShowStatus( void ) const
{
  G4cout << "Analyzer Status\n"
	 << "  File : " << filename_ << "\n"
	 << "  Analysis : " << fActive_ << "\n"
	 << "  Trigger : " << fTriggered << "\n"
	 << "  # of triggered events : " << trigNum 
	 << G4endl;
}

void Analysis::InitializeEvent( void )
{
  for (int i=0; i<NumOfPlaneCFT; i++) {
    event.FiberHits[i]=0;
    for (int j=0; j<MaxHit; j++) {
      event.FiberSeg[i][j] = -1;
      event.FiberTime[i][j] = -999.;
      event.FiberEdep[i][j] = -999.;
      event.FiberPID[i][j] = -1;
    }
  }

  event.CrystalHits=0;
  for (int i=0; i<MaxHit; i++) {
    event.CrystalSeg[i] = -1;
    event.CrystalTime[i] = -999.;
    event.CrystalEdep[i] = -999.;
    event.CrystalPID[i] = -1;
  }

  event.nHitVP1 = 0;
  event.nHitVP2 = 0;
  event.nHitVP3 = 0;
  event.nHitVP4 = 0;
  event.nHitVP5 = 0;
  for (int i=0; i<MaxHit; i++) {
    event.xPosVP1[i] = -999.;
    event.yPosVP1[i] = -999.;
    event.xPosVP2[i] = -999.;
    event.yPosVP2[i] = -999.;
    event.xPosVP3[i] = -999.;
    event.yPosVP3[i] = -999.;
    event.xPosVP4[i] = -999.;
    event.yPosVP4[i] = -999.;
    event.xPosVP5[i] = -999.;
    event.yPosVP5[i] = -999.;
  }

  // event.HitPat1 = -999;
  // event.HitPat2 = -999;
  // event.HitPat3 = -999;
  // event.HitPat4 = -999;
  // event.HitPat5 = -999;
  // event.HitPat6 = -999;
  // event.HitPat7 = -999;
  // event.HitPat8 = -999;
  // event.HitPat9 = -999;
  // event.HitPat10 = -999;
  // event.HitPat31 = -999;
  // event.HitPat32 = -999;
  // event.HitPat33 = -999;
  // event.HitPat34 = -999;
  // event.HitPat35 = -999;
  // event.HitPat36 = -999;
  // event.HitPat37 = -999;
  // event.HitPat38 = -999;
  // event.HitPat39 = -999;
  // event.HitPat40 = -999;
  
  event.nHitCH = 0;
  event.nHitFTOF = 0;
  event.nHitT0 = 0;
  event.nHitAC = 0;

  for (int i=0; i<MaxHit; i++) {
    event.HitSegCH[i] = -1;
    event.HitSegFTOF[i] = -1;
    event.lcCherenkov[i] = -1;
  }

  for (int i=0; i<5; i++) {
    event.xPosBGOVP[i] = -999.;
    event.yPosBGOVP[i] = -999.;
    event.zPosBGOVP[i] = -999.;
  }
}

void Analysis::DefineHistograms()
{
  new TFile( filename_, "recreate" );
  fActive_=true;

  // Primaray Generation


  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Yd scat primary
  /*
  tree->Branch("momMeson",&event.momMeson, "momMeson/D");
  tree->Branch("momSigma",&event.momSigma, "momSigma/D");
  tree->Branch("momScatSigma",&event.momScatSigma, "momScatSigma/D");
  tree->Branch("momScatProton",&event.momScatProton, "momScatProton/D");
  tree->Branch("momDecayNucleon",&event.momDecayNucleon, "momDecayNucleon/D");
  tree->Branch("momDecayPi",&event.momDecayPi, "momDecayPi/D");
  */
  tree->Branch("thetaMeson",&event.thetaMeson, "thetaMeson/D");
  tree->Branch("phiMeson",&event.phiMeson, "phiMeson/D");
  tree->Branch("thetaMesonCM",&event.thetaMesonCM, "thetaMesonCM/D");
  tree->Branch("phiMesonCM",&event.phiMesonCM, "phiMesonCM/D");
  tree->Branch("thetaScatHypCM",&event.thetaScatHypCM, "thetaScatHypCM/D");
  tree->Branch("phiScatHypCM",&event.phiScatHypCM, "phiScatHypCM/D");
  tree->Branch("thetaScatPLab",&event.thetaScatPLab, "thetaScatPLab/D");
  tree->Branch("thetaScatHypLab",&event.thetaScatHypLab, "thetaScatHypLab/D");
  tree->Branch("momVectorScatMeson",event.momVectorScatMeson, "momVectorScatMeson[3]/D");
  tree->Branch("momVectorHypBeam",event.momVectorHypBeam, "momVectorHypBeam[3]/D");
  tree->Branch("momVectorHypScat",event.momVectorHypScat, "momVectorHypScat[3]/D");
  tree->Branch("momVectorProtonScat",event.momVectorProtonScat, "momVectorProtonScat[3]/D");
  tree->Branch("momVectorDecayPi",event.momVectorDecayPi, "momVectorDecayPi[3]/D");
  tree->Branch("momVectorDecayNucleon",event.momVectorDecayNucleon, "momVectorDecayNucleon[3]/D");
  /*
  tree->Branch("fermiMom",event.fermiMom, "fermiMom[3]/D");
  tree->Branch("protonStopPos",event.protonStopPos, "protonStopPos[3]/D");
  */
  tree->Branch("primaryVertex",event.primaryVertex, "primaryVertex[3]/D");
  tree->Branch("scatPos",event.scatPos, "scatPos[3]/D");
  tree->Branch("NNscatPos",event.NNscatPos, "NNscatPos[3]/D");
  tree->Branch("PiNscatPos",event.PiNscatPos, "PiNscatPos[3]/D");
  tree->Branch("decayPos",event.decayPos, "decayPos[3]/D");

  tree->Branch("decayEventFlag",&event.decayEventFlag, "decayEventFlag/I");
  tree->Branch("decayFlag",&event.decayFlag, "decayFlag/I");
  tree->Branch("scatFlag",&event.scatFlag, "scatFlag/I");
  tree->Branch("scatTarget",&event.scatTarget, "scatTarget/I");
  tree->Branch("NNscatFlag",&event.NNscatFlag, "NNscatFlag/I");
  tree->Branch("NNscatTarget",&event.NNscatTarget, "NNscatTarget/I");
  tree->Branch("PiNscatFlag",&event.PiNscatFlag, "PiNscatFlag/I");
  tree->Branch("PiNscatTarget",&event.PiNscatTarget, "PiNscatTarget/I");

  tree->Branch("fLengthInH",&event.fLengthInH, "fLengthInH/D");
  tree->Branch("fLengthInD",&event.fLengthInD, "fLengthInD/D");
  tree->Branch("fLengthTotal",&event.fLengthTotal, "fLengthTotal/D");

  char buf1[100], buf2[100];
  for (int i=0; i<NumOfPlaneCFT; i++) {
    sprintf(buf1, "Fiber%dHits", i);
    sprintf(buf2, "Fiber%dHits/I", i);
    tree->Branch(buf1, &event.FiberHits[i],  buf2);
  }

  for (int i=0; i<NumOfPlaneCFT; i++) {
    sprintf(buf1, "Fiber%dSeg", i);
    sprintf(buf2, "Fiber%dSeg[Fiber%dHits]/I", i, i);
    tree->Branch(buf1,  event.FiberSeg[i],  buf2);
    sprintf(buf1, "Fiber%dTime", i);
    sprintf(buf2, "Fiber%dTime[Fiber%dHits]/D", i, i);
    tree->Branch(buf1, event.FiberTime[i], buf2);
    sprintf(buf1, "Fiber%dEdep", i);
    sprintf(buf2, "Fiber%dEdep[Fiber%dHits]/D", i, i);
    tree->Branch(buf1, event.FiberEdep[i], buf2);
    sprintf(buf1, "Fiber%dPID", i);
    sprintf(buf2, "Fiber%dPID[Fiber%dHits]/I", i, i);
    tree->Branch(buf1, event.FiberPID[i], buf2);
  }

  tree->Branch("CrystalHits", &event.CrystalHits,  "CrystalHits/I");

  sprintf(buf1, "CrystalSeg");
  sprintf(buf2, "CrystalSeg[CrystalHits]/I");
  tree->Branch(buf1,  event.CrystalSeg,  buf2);
  sprintf(buf1, "CrystalTime");
  sprintf(buf2, "CrystalTime[CrystalHits]/D");
  tree->Branch(buf1, event.CrystalTime, buf2);
  sprintf(buf1, "CrystalEdep");
  sprintf(buf2, "CrystalEdep[CrystalHits]/D");
  tree->Branch(buf1, event.CrystalEdep, buf2);
  sprintf(buf1, "CrystalPID");
  sprintf(buf2, "CrystalPID[CrystalHits]/I");
  tree->Branch(buf1, event.CrystalPID, buf2);

  tree->Branch("nHitVP1", &event.nHitVP1,  "nHitVP1/I");
  tree->Branch("xPosVP1", event.xPosVP1, "xPosVP1[nHitVP1]/D");
  tree->Branch("yPosVP1", event.yPosVP1, "yPosVP1[nHitVP1]/D");
  tree->Branch("nHitVP2", &event.nHitVP2,  "nHitVP2/I");
  tree->Branch("xPosVP2", event.xPosVP2, "xPosVP2[nHitVP2]/D");
  tree->Branch("yPosVP2", event.yPosVP2, "yPosVP2[nHitVP2]/D");
  tree->Branch("nHitVP3", &event.nHitVP3,  "nHitVP3/I");
  tree->Branch("xPosVP3", event.xPosVP3, "xPosVP3[nHitVP3]/D");
  tree->Branch("yPosVP3", event.yPosVP3, "yPosVP3[nHitVP3]/D");
  tree->Branch("nHitVP4", &event.nHitVP4,  "nHitVP4/I");
  tree->Branch("xPosVP4", event.xPosVP4, "xPosVP4[nHitVP4]/D");
  tree->Branch("yPosVP4", event.yPosVP4, "yPosVP4[nHitVP4]/D");
  tree->Branch("nHitVP5", &event.nHitVP5,  "nHitVP5/I");
  tree->Branch("xPosVP5", event.xPosVP5, "xPosVP5[nHitVP5]/D");
  tree->Branch("yPosVP5", event.yPosVP5, "yPosVP5[nHitVP5]/D");

  // tree->Branch("hitPat1", &event.HitPat1,  "HitPat1/I");
  // tree->Branch("hitPat2", &event.HitPat2,  "HitPat2/I");
  // tree->Branch("hitPat3", &event.HitPat3,  "HitPat3/I");
  // tree->Branch("hitPat4", &event.HitPat4,  "HitPat4/I");
  // tree->Branch("hitPat5", &event.HitPat5,  "HitPat5/I");
  // tree->Branch("hitPat6", &event.HitPat6,  "HitPat6/I");
  // tree->Branch("hitPat7", &event.HitPat7,  "HitPat7/I");
  // tree->Branch("hitPat8", &event.HitPat8,  "HitPat8/I");
  // tree->Branch("hitPat9", &event.HitPat9,  "HitPat9/I");
  // tree->Branch("hitPat10", &event.HitPat10,  "HitPat10/I");
  // tree->Branch("hitPat31", &event.HitPat31,  "HitPat31/I");
  // tree->Branch("hitPat32", &event.HitPat32,  "HitPat32/I");
  // tree->Branch("hitPat33", &event.HitPat33,  "HitPat33/I");
  // tree->Branch("hitPat34", &event.HitPat34,  "HitPat34/I");
  // tree->Branch("hitPat35", &event.HitPat35,  "HitPat35/I");
  // tree->Branch("hitPat36", &event.HitPat36,  "HitPat36/I");
  // tree->Branch("hitPat37", &event.HitPat37,  "HitPat37/I");
  // tree->Branch("hitPat38", &event.HitPat38,  "HitPat38/I");
  // tree->Branch("hitPat39", &event.HitPat39,  "HitPat39/I");
  // tree->Branch("hitPat40", &event.HitPat40,  "HitPat40/I");


  tree->Branch("nHitCH", &event.nHitCH,  "nHitCH/I");
  tree->Branch("HitSegCH", event.HitSegCH, "HitSegCH[nHitCH]/I");
  tree->Branch("nHitFTOF", &event.nHitFTOF,  "nHitFTOF/I");
  tree->Branch("HitSegFTOF", event.HitSegFTOF, "HitSegCH[nHitFTOF]/I");
  tree->Branch("lcCherenkov", event.lcCherenkov, "lcCherenkov[nHitFTOF]/I");
  tree->Branch("nHitT0", &event.nHitT0,  "nHitT0/I");
  tree->Branch("nHitAC", &event.nHitAC,  "nHitAC/I");

  tree->Branch("xPosBGOVP", event.xPosBGOVP, "xPosBGOVP[5]/D");
  tree->Branch("yPosBGOVP", event.yPosBGOVP, "yPosBGOVP[5]/D");
  tree->Branch("zPosBGOVP", event.zPosBGOVP, "zPosBGOVP[5]/D");
}

void Analysis::
PrintHitsInformation( const G4Event *anEvent, 
		      std::ostream &ost ) const
{
  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  std::ios::fmtflags oldFlags = ost.flags();
  std::size_t preSiz = ost.precision();
  ost.setf( std::ios::fixed );
  ost << "=== Event#" << std::setw(7) << anEvent->GetEventID() 
      << " ===" << std::endl;

  ost << " --- PrimaryInformation(Yd scattering) ---" << std::endl; 
  //thetaMeson  phiMeson thetaMesonCM  phiMesonCM thetaScatHypCM phiScatHypCM
  ost.precision(2);
  /*
  ost << event.thetaMeson << std::setw(8) 
      << event.phiMeson   << std::setw(8)
      << event.thetaMesonCM   << std::setw(8) 
      << event.phiMesonCM   << std::setw(8) 
      << event.thetaScatHypCM   << std::setw(8) 
      << event.phiScatHypCM     << std::setw(8) 
      << event.thetaScatPLab   << std::setw(8) 
      << event.thetaScatHypLab   << std::setw(8) 
      << std::endl; 
  */
  ost << event.thetaMeson << "  " 
      << event.phiMeson   << "  "
      << event.thetaMesonCM   << "  " 
      << event.phiMesonCM   << "  " 
      << event.thetaScatHypCM   << "  " 
      << event.phiScatHypCM     << "  " 
      << event.thetaScatPLab   << "  " 
      << event.thetaScatHypLab   << "  " 
      << event.beammom   << "  " 
      << std::endl; 

  ost.precision(7);
  ost << "momVectorScatMeson" << std::setw(15) 
      << event.momVectorScatMeson[0]   << std::setw(15) 
      << event.momVectorScatMeson[1]   << std::setw(15) 
      << event.momVectorScatMeson[2]   << std::endl;
  ost << "momVectorHypBeam" << std::setw(15) 
      << event.momVectorHypBeam[0]   << std::setw(15) 
      << event.momVectorHypBeam[1]   << std::setw(15) 
      << event.momVectorHypBeam[2]   << std::endl;
  ost << "momVectorHypScat" << std::setw(15) 
      << event.momVectorHypScat[0]   << std::setw(15) 
      << event.momVectorHypScat[1]   << std::setw(15) 
      << event.momVectorHypScat[2]   << std::endl;
  ost << "momVectorProtonScat" << std::setw(15) 
      << event.momVectorProtonScat[0]   << std::setw(15) 
      << event.momVectorProtonScat[1]   << std::setw(15) 
      << event.momVectorProtonScat[2]   << std::endl;
  ost << "momVectorDecayPi" << std::setw(15) 
      << event.momVectorDecayPi[0]   << std::setw(15) 
      << event.momVectorDecayPi[1]   << std::setw(15) 
      << event.momVectorDecayPi[2]   << std::endl;
  ost << "momVectorDecayNucleon" << std::setw(15) 
      << event.momVectorDecayNucleon[0]   << std::setw(15) 
      << event.momVectorDecayNucleon[1]   << std::setw(15) 
      << event.momVectorDecayNucleon[2]   << std::endl;
  ost << "primaryVertex" << std::setw(15) 
      << event.primaryVertex[0]   << std::setw(15) 
      << event.primaryVertex[1]   << std::setw(15) 
      << event.primaryVertex[2]   << std::endl;
  ost << "scatPos" << std::setw(15) 
      << event.scatPos[0]   << std::setw(15) 
      << event.scatPos[1]   << std::setw(15) 
      << event.scatPos[2]   << std::endl;
  ost << "NNscatPos" << std::setw(15) 
      << event.NNscatPos[0]   << std::setw(15) 
      << event.NNscatPos[1]   << std::setw(15) 
      << event.NNscatPos[2]   << std::endl;
  ost << "PiNscatPos" << std::setw(15) 
      << event.PiNscatPos[0]   << std::setw(15) 
      << event.PiNscatPos[1]   << std::setw(15) 
      << event.PiNscatPos[2]   << std::endl;
  ost << "decayPos" << std::setw(15) 
      << event.decayPos[0]   << std::setw(15) 
      << event.decayPos[1]   << std::setw(15) 
      << event.decayPos[2]   << std::endl;
  ost << "Flags" << std::setw(15) 
      << event.decayFlag   << std::setw(15) 
      << event.scatFlag   << std::setw(15) 
      << event.scatTarget   << std::setw(15) 
      << event.NNscatFlag   << std::setw(15) 
      << event.NNscatTarget   << std::setw(15) 
      << event.PiNscatFlag   << std::setw(15) 
      << event.PiNscatTarget   << std::endl;

  /*
  //Generated point infomation
  //Particle Mass, Gene Pos, Gene Mom, Gene Theta, Gene Phi, Beam Mom, Beam u&v
  ost << " --- PrimaryInformation ---" << std::endl; 
  ost.precision(4);
  ost << ( (pInfo_->charge>=0.) ? " " : " ")
      << std::setw(7) << pInfo_->mass/GeV << " [GeV] ";
  ost.precision(2);
  ost << std::setw(9) << pInfo_->x/mm 
      << std::setw(9) << pInfo_->y/mm
      << std::setw(9) << pInfo_->z/mm
      << " [mm] ";
  ost.precision(5);
  ost << std::setw(7) << pInfo_->p/GeV << " [GeV/c] ";
  ost.precision(2);
  ost << std::setw(6) << pInfo_->theta
      << std::setw(8) << pInfo_->phi << " [degree] ";
  ost.precision(5);
  ost << std::setw(7) << pInfo_->pbeam/GeV << " [GeV/c] "
      << std::setw(8) << pInfo_->ubeam << " "
      << std::setw(8) << pInfo_->vbeam << std::endl;
  ost.precision(3);
  */

  int nhCFT=0;

  for (int i=0; i<NumOfPlaneCFT;i++) {
    nhCFT += event.FiberHits[i];
  }

  ost << " --- CFT Fiber --- Total # of Hits : "
      << std::setw(3) << nhCFT << std::endl;

  for (int i=0; i<NumOfPlaneCFT;i++) {
    for (int j=0; j<event.FiberHits[i]; j++) {
      ost << std::setw(3) << 61+i
	  << std::setw(5) << event.FiberSeg[i][j]
	  << std::setw(15) << event.FiberTime[i][j]/ns << " [ns]"
	  << std::setw(15) << event.FiberEdep[i][j]/MeV << " [MeV]"
	  << std::setw(7) << event.FiberPID[i][j] 
	   << std::endl;

    }
  }


  ost << " --- Crystal" << " --- Total # of Hits : "
      << std::setw(3) << event.CrystalHits << std::endl;
  for (int j=0; j<event.CrystalHits; j++) {
    ost << std::setw(3) << 81
	<< std::setw(5) << event.CrystalSeg[j]
	<< std::setw(15) << event.CrystalTime[j]/ns << " [ns]"
	<< std::setw(15) << event.CrystalEdep[j]/MeV << " [MeV]"
	<< std::setw(7) << event.CrystalPID[j] 
	<< std::endl;
  }

  G4int colIdChamber   = SDMan->GetCollectionID( "chamberSD" );

  SKSChamberHitsCollection *ChamberHC;
  // Chamber Hit
  if (HCE) ChamberHC = (SKSChamberHitsCollection*)(HCE->GetHC(colIdChamber));
  if (ChamberHC) {
    int n_hit = ChamberHC->entries();
    ost << " --- DC" << " --- Total # of Hits : "
	<< std::setw(3) << n_hit << std::endl;

    for (G4int i=0;i<n_hit;i++){
      int planeID = (*ChamberHC)[i]->GetDetectorID();
      int hitWire = (*ChamberHC)[i]->GetHitWire();
      double dlength = (*ChamberHC)[i]->GetDRlength();
      //bool sensitiveHit = (*ChamberHC)[i]->GetSensitiveHit();

      ost << planeID << " " << hitWire << " " << dlength <<std::endl;
    }
  }

  G4int colIdCounter   = SDMan->GetCollectionID( "counterSD" );

  SKSCounterHitsCollection *CounterHC;
  CounterHC = dynamic_cast<SKSCounterHitsCollection *>(HCE->GetHC(colIdCounter));
  int nhPiV = 0;
  if (CounterHC) {
    nhPiV = CounterHC->entries();
  }
  

  if (CounterHC) {
    int n_hit = CounterHC->entries();
    ost << " --- Counter" << " --- Total # of Hits : "
	<< std::setw(3) << n_hit << std::endl;

    for (G4int i=0;i<n_hit;i++){
      int planeID = (*CounterHC)[i]->GetDetectorID();
      int channelID = (*CounterHC)[i]->GetChannelID();
      double time= (*CounterHC)[i]->GetTime();
      double de = (*CounterHC)[i]->GetEdepAbs();

      ost << planeID << " " << channelID << " " << time << " " << de << std::endl;
    }
  }


  ost.flags( oldFlags );
  ost.precision( preSiz );
}

void Analysis::SetDataFile( const char *datafile )
{
  DataFile_.open( datafile );
}

void Analysis::OpenEffStudyFile( const char *datafile )
{
  ConfMan *confMan = ConfMan::GetConfManager();
  ReactionMode_ = confMan->ReactionMode();

  if (ReactionMode_ == 95 || ReactionMode_ == 99 || ReactionMode_ == 100 || ReactionMode_ == 101) {
    fpEff = fopen(datafile,"w");
    if (!fpEff) {
      std::cerr << "cannot open OutputEffStudy.txt\n" << std::endl;
      exit(-1);
    }
  }

}


