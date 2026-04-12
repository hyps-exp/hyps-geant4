/*
  Analysis.hh
  2007/4  K.Shirotori
*/

#ifndef Analysis_h
#define Analysis_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "RootHelper.hh"
#include "DetectorID.hh"

#include <fstream>

class G4Run;
class G4Event;
//class PrimaryInfo;

#define MaxHit 30

struct Event{
  //Primary for Yd scattering
  G4double thetaMeson;
  G4double phiMeson;
  G4double thetaMesonCM;
  G4double phiMesonCM;
  G4double thetaScatHypCM;
  G4double phiScatHypCM;
  G4double thetaScatPLab;
  G4double thetaScatHypLab;
  G4double momVectorScatMeson[3];
  G4double momVectorHypBeam[3];
  G4double momVectorHypScat[3];
  G4double momVectorProtonScat[3];
  G4double momVectorDecayPi[3];
  G4double momVectorDecayNucleon[3];
  G4double primaryVertex[3];
  G4double scatPos[3];
  G4double NNscatPos[3];
  G4double PiNscatPos[3];
  G4double decayPos[3];
  G4double beammom;
  /*
  G4double fermiMom[3];
  G4double protonStopPos[3];
  G4int    scatPartTpcHitNum;
  G4int    scatPartCrystalHitNum;
  G4int    scatPartSi1HitNum;
  G4int    scatPartSi2HitNum;
  G4int    scatPartSi3HitNum;
  */
  G4int    decayEventFlag;
  G4int    decayFlag;
  G4int    scatFlag;
  G4int    scatTarget;
  G4int    NNscatFlag;
  G4int    NNscatTarget;
  G4int    PiNscatFlag;
  G4int    PiNscatTarget;

  G4double momVectorPiMinusFromK0[3];

  G4double fLengthInH;
  G4double fLengthInD;
  G4double fLengthTotal;

  //Fiber
  int      FiberHits[NumOfPlaneCFT];
  int      FiberSeg[NumOfPlaneCFT][MaxHit];
  double   FiberTime[NumOfPlaneCFT][MaxHit];
  double   FiberEdep[NumOfPlaneCFT][MaxHit];
  int      FiberPID[NumOfPlaneCFT][MaxHit];


  //Crystal
  int      CrystalHits;
  int      CrystalSeg[MaxHit];
  double   CrystalTime[MaxHit];
  double   CrystalEdep[MaxHit];
  int      CrystalPID[MaxHit];

  // Chamber HitPat
  int HitPat1;
  int HitPat2;
  int HitPat3;
  int HitPat4;
  int HitPat5;
  int HitPat6;
  int HitPat7;
  int HitPat8;
  int HitPat9;
  int HitPat10;
  int HitPat31;
  int HitPat32;
  int HitPat33;
  int HitPat34;
  int HitPat35;
  int HitPat36;
  int HitPat37;
  int HitPat38;
  int HitPat39;
  int HitPat40;

  // Virtual Plane
  int      nHitVP1;
  double   xPosVP1[MaxHit];
  double   yPosVP1[MaxHit];
  int      nHitVP2;
  double   xPosVP2[MaxHit];
  double   yPosVP2[MaxHit];
  int      nHitVP3;
  double   xPosVP3[MaxHit];
  double   yPosVP3[MaxHit];
  int      nHitVP4;
  double   xPosVP4[MaxHit];
  double   yPosVP4[MaxHit];
  int      nHitVP5;
  double   xPosVP5[MaxHit];
  double   yPosVP5[MaxHit];


  // CH, FTOF 
  int nHitCH;
  int   HitSegCH[MaxHit];
  int nHitFTOF;
  int   HitSegFTOF[MaxHit];
  int   lcCherenkov[MaxHit];
  int nHitT0;
  int nHitAC;

  double xPosBGOVP[5];
  double yPosBGOVP[5];
  double zPosBGOVP[5];

  //for virtual particle daughters
  int nDaughter;
  int daughterPID[MaxHit];
  G4ThreeVector daughterMom[MaxHit];

};

class Analysis
{
public:
  Analysis( G4String histname );
  virtual ~Analysis();

public:
  void BeginOfRun( const G4Run *aRun );
  void EndOfRun( const G4Run *aRun );
  void BeginOfPrimaryAction();
  void BeginOfEvent( );
  void EndOfEvent( const G4Event *anEvent );
  //void PrimaryGeneration( const PrimaryInfo *info );
  
  void SetScatMesonMomentum(G4ThreeVector mesonMom);
  void GetScatMesonMomentum(G4double  *mesonMom);
  void SetHypBeamMomentum(G4ThreeVector hypBeamMom);
  void GetHypBeamMomentum(G4double *hypBeamMom);
  void SetScatHypMomentum(G4ThreeVector hypScatMom);
  void GetScatHypMomentum(G4double *hypScatMom);
  void SetScatProtonMomentum(G4ThreeVector protonScatMom);
  void GetScatProtonMomentum(G4double *protonScatMom);
  void SetDecayPiMomentum(G4ThreeVector piMom);
  void GetDecayPiMomentum(G4double *piMom);
  void SetDecayNucleonMomentum(G4ThreeVector nucleonMom);
  void GetDecayNucleonMomentum(G4double *nucleonMom);
  void SetDaughterMomentum(G4int type, G4ThreeVector Mom);
  int  GetDaughterPID(G4int i);
  G4ThreeVector GetDaughterMomentum(G4int i);
  int  GetNDaughter() { return event.nDaughter;}
  void SetThetaMeson(double theta);
  void SetPhiMeson(double phi);
  void SetThetaMesonCM(double theta);
  G4double GetThetaMeson();
  void SetPhiMesonCM(double phi);
  void SetThetaScatHypCM(double theta);
  double GetThetaScatHypCM();
  void SetPhiScatHypCM(double phi);
  void SetThetaScatPLab(double theta);
  double GetThetaScatPLab();
  void SetThetaScatHypLab(double theta);
  double GetThetaScatHypLab();
  void SetBeamMomentum(G4double mom);
  double GetBeamMomentum();
  void SetPrimaryVertex(G4ThreeVector pos);
  G4ThreeVector GetPrimaryVertex();
  void SetScatPos(G4ThreeVector pos);
  G4ThreeVector GetScatPos();
  void SetNNScatPos(G4ThreeVector pos);
  void SetPiNScatPos(G4ThreeVector pos);
  void SetDecayPos(G4ThreeVector pos);
  G4ThreeVector GetDecayPos();
  void SetDecayEventFlag();
  void SetDecayFlag();
  void SetDecayFlag(int decayMode);
  void SetScatFlag();
  void SetScatFlag(int reactMode);
  void SetScatTarget(int tgt);
  void SetNNScatFlag();
  G4int  GetNNScatFlag() { return event.NNscatFlag;}
  void SetNNScatTarget(int tgt);
  void SetPiNScatFlag();
  void SetPiNScatTarget(int tgt);
  void SetFlightLengthInTarget(double dxInH, double dxInD);
  void SetPiMinusMomentumFromK0(G4ThreeVector piMom);

  void SetFileName( const G4String &filename ) { filename_=filename; }
  void DefineHistograms( void );
  G4bool GetTriggerStatus( void ) const { return fTriggered; }
  void SaveFile( void ) const;
  void Terminate( void ) const;
  const G4String &GetFileName( void ) const { return filename_; }
  void SetActive( void ) { fActive_=true; }
  void SetInActive( void ) { fActive_=false; }
  void ShowStatus( void ) const;

  void SetDataFile( const char *datafile );
  //   void SetDataFile( const G4String &datafile ) { datafile_=datafile; }
  //   const G4String &GetDataFile( void ) const { return datafile_; }
  void OpenEffStudyFile( const char *datafile );
  void SetPrimaryFlag( G4bool flag) { fPrimary_ = flag; }
  G4bool GetPrimaryFlag( void ) { return fPrimary_; }

private:
  Event event;
  G4String filename_;
  G4bool fActive_;
  G4bool fTriggered;
  G4int ReactionMode_;
  G4bool fPrimary_;

  G4int trigNum;
  //PrimaryInfo *pInfo_;
  std::ofstream DataFile_;
  //   G4String datafile_;

private:
  void PrintHitsInformation( const G4Event *anEvent, 
			    std::ostream &ost ) const;
public:
  void InitializeEvent(void);

};

#endif
