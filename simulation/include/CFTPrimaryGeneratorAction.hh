//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN03PrimaryGeneratorAction.hh,v 1.6 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CFTPrimaryGeneratorAction_h
#define CFTPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class G4Event;

class Analysis;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CFTPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
  CFTPrimaryGeneratorAction(Analysis *anaMan);
   ~CFTPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void OpenEffStudyFile( const char *datafile );

  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4  class
    G4ParticleDefinition* particle;

  G4ParticleGun *gunKM_, *gunKP_;
  G4ParticleGun *gunPiM_, *gunPiP_, *gunPiZ_;
  G4ParticleGun  *gunMu_, *gunG_; 
  G4ParticleGun *gunLambda_;
  G4ParticleGun *gunProton_, *gunNeutron_;

  int ReactionMode_;
  double PrimaryE_;
  double GenerateThetaCM_;
  Analysis *anaMan_;

  void MakeGuns();
  void GenerateProton(G4Event*);
  void GeneratePiMinus(G4Event* anEvent);
  void GenerateScatProtonCheck(G4Event* anEvent);
  void GenerateScatProtonCheck_SigmaPlusP(G4Event* anEvent);
  void GeneratePiKSigma(G4Event* anEvent);
  void GeneratePiKSigmaPlus(G4Event* anEvent);
  void GeneratePiKSigmaStar(G4Event* anEvent);
  void GeneratePiKSigmaPlusStar(G4Event* anEvent);
  void GeneratePiMinusP_Elastic(G4Event* anEvent);
  void GeneratePiMinusP_ElasticChargeExchange(G4Event* anEvent);
  void GeneratePiPlusP_Elastic(G4Event* anEvent);
  void GeneratePiMinusP_InElastic1(G4Event* anEvent);
  void GeneratePiMinusP_InElastic2(G4Event* anEvent);
  void GeneratePiMinusP_InElastic3(G4Event* anEvent);
  void GeneratePiMinusP_InElastic4(G4Event* anEvent);
  void GeneratePiMinusP_InElastic5(G4Event* anEvent);
  void GeneratePiMinusP_InElastic6(G4Event* anEvent);
  void GeneratePPScat(G4Event* anEvent);
  void GeneratePiMinusP_KPlusPiLambda(G4Event* anEvent);
  void GeneratePiMinusP_K0Lambda(G4Event* anEvent);
  void GeneratePiMinusP_K0Lambda_Scat(G4Event* anEvent);
  void GenerateVPFourBodyDecay(G4Event* anEvent);
  void GenerateVPFourBodyDecay_Scat(G4Event* anEvent);
  bool NeutronP_Scattering(G4Event* anEvent, G4ThreeVector momentumDecayN, G4ThreeVector Pos0);
  bool ProtonP_Scattering(G4Event* anEvent, G4ThreeVector momentumDecayN, G4ThreeVector Pos0);
  bool PiMinusP_Scattering(G4Event* anEvent, G4ThreeVector momentumDecayPi, G4ThreeVector Pos0);
  bool LambdaP_Scattering(G4Event* anEvent, G4ThreeVector momentumLambda, G4ThreeVector Pos0);
  void GenerateKMinusP_KPlusXi(G4Event* anEvent);
  void GenerateKMinusP_KPlusXi_Scat(G4Event* anEvent);
  void GenerateGammaP_KPlusLambda(G4Event* anEvent);
  void GenerateGammaP_KPlusLambda_Scat(G4Event* anEvent);
  void GenerateGammaP_KPlusSigma(G4Event* anEvent);  
  void GenerateGammaP_PhiP(G4Event* anEvent);
  void GenerateUniform_PiPlus(G4Event* anEvent);
  void GenerateUniform_Proton(G4Event* anEvent);
  void GeneratePiKBG(G4Event* anEvent);
  void GenerateUniform(G4Event* anEvent);
  void GeneratePiKSigmaScat2(G4Event* anEvent);
  void GeneratePiKSigmaPlusScat2(G4Event* anEvent);
  void GeneratePiMinusBeam(G4Event* anEvent);
  void GeneratePiPlusBeam(G4Event* anEvent);
  void EfficiencyStudyForScatPart(G4Event* anEvent, int reactionMode);
  void EfficiencyStudyForLambda(G4Event* anEvent);
  void GenerateSigmaBeam();
  void calcThetaPhi(G4ThreeVector vec, double *theta, double *phi);
  int  getTargetFlag(G4ThreeVector pos);
  G4bool decayCheck(G4double ctau, G4double momentum, G4double mass, G4double dx);
  G4bool scatteringCheck(G4double rate, G4double dx);
  G4double calcEnergyDeposit(G4double momentum, G4double mass, G4double dx, G4ThreeVector pos, G4int flagTgtType);
  G4double calc_dE_dx(double beta, int flagTgtType);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


