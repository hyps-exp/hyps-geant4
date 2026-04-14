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
// $Id: ExN03PrimaryGeneratorAction.cc,v 1.7 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CFTPrimaryGeneratorAction.hh"

#include "CFTDetectorConstruction.hh"

#include "RadDeg.hh"
#include "DCGeomMan.hh"
#include "Analysis.hh"
#include "ConfMan.hh"

#include "Kinema2Body.hh"
#include "Kinema3Body.hh"
#include "Kinema3Resonance.hh"
#include "DecayWithPolarization.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

#include "Kinematics.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FILE *fpIn;
FILE *fpSigBeam;

const double AsymPara_Lambda = 0.75;

using namespace CLHEP;

CFTPrimaryGeneratorAction::CFTPrimaryGeneratorAction(Analysis *analysisManager)
  : anaMan_(analysisManager)
{
  MakeGuns();

  ConfMan *confMan = ConfMan::GetConfManager();
  ReactionMode_ = confMan->ReactionMode();
  PrimaryE_     = confMan->GetPrimaryE();
  GenerateThetaCM_ = confMan->GetGenerateThetaCM();
  
  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particle = particleTable->FindParticle(particleName="proton");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(70.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.*cm,0.*cm));

    if (ReactionMode_ == 96 || ReactionMode_ == 102) {
    fpSigBeam = fopen("SigmaBeamStudy.txt","r");
    if (!fpSigBeam) {
      std::cerr << "Cannot open SigmaBeamStudy.txt" << std::endl;
      exit(-1);
    }
  } 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTPrimaryGeneratorAction::~CFTPrimaryGeneratorAction()
{
  delete particleGun;
  //delete particle;
}

void CFTPrimaryGeneratorAction::OpenEffStudyFile( const char *datafile )
{
  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  if (ReactionMode == 95 || ReactionMode == 99 || ReactionMode == 100 || ReactionMode == 101) {
    fpIn = fopen(datafile,"r");
    if (!fpIn) {
      std::cerr << "Cannot open InputEffStudy.txt" << std::endl;
      exit(-1);
    }
  }
}


void CFTPrimaryGeneratorAction::MakeGuns()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  gunKM_  = new G4ParticleGun( G4KaonMinus::KaonMinusDefinition() );
  gunKP_  = new G4ParticleGun( G4KaonPlus::KaonPlusDefinition() );
  gunPiM_ = new G4ParticleGun( G4PionMinus::PionMinusDefinition() );
  gunPiP_ = new G4ParticleGun( G4PionPlus::PionPlusDefinition() );
  gunPiZ_ = new G4ParticleGun( G4PionZero::PionZeroDefinition() );
  gunMu_  = new G4ParticleGun( G4MuonMinus::MuonMinusDefinition() );
  gunG_   = new G4ParticleGun( G4Gamma::GammaDefinition() );
  gunLambda_ = new G4ParticleGun( G4Lambda::LambdaDefinition() );
  gunProton_ = new G4ParticleGun( G4Proton::ProtonDefinition() );
  gunNeutron_ = new G4ParticleGun( G4Neutron::NeutronDefinition() );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if (ReactionMode_ != 37)
    anaMan_->BeginOfPrimaryAction();

  switch(ReactionMode_) {
  case 1 :
    GenerateProton(anEvent);
    break;
  case 2 :
    GenerateScatProtonCheck(anEvent);
    break;
  case 3 :
    GenerateScatProtonCheck_SigmaPlusP( anEvent);
    break;
  case 4 :
    GeneratePiKSigma(anEvent);
    break;
  case 5 :
    GenerateUniform(anEvent);
    break;
  case 6 :
    GeneratePiKSigmaScat2(anEvent);
    break;
  case 7 :
    GeneratePiMinusBeam(anEvent);
    break;
  case 8 :
    GeneratePiPlusBeam(anEvent);
    break;
  case 9 :
    GeneratePiKSigmaPlus(anEvent);
    break;
  case 10 :
    GeneratePiKSigmaPlusScat2(anEvent);
    break;
  case 11 :
    GeneratePiKBG(anEvent);
    break;
  case 12 :
    GeneratePiKSigmaStar(anEvent);
    break;
  case 13 :
    GeneratePiKSigmaPlusStar(anEvent);
    break;
  case 14 :
    GeneratePiMinusP_Elastic(anEvent);
    break;
  case 15 :
    GeneratePiPlusP_Elastic(anEvent);
    break;
  case 16 :
    GeneratePiMinusP_InElastic1(anEvent);
    break;
  case 17 :
    GeneratePiMinusP_InElastic2(anEvent);
    break;
  case 18 :
    GeneratePiMinusP_InElastic3(anEvent);
    break;
  case 19 :
    GeneratePiMinusP_InElastic4(anEvent);
    break;
  case 20 :
    GeneratePiMinusP_InElastic5(anEvent);
    break;
  case 21 :
    GeneratePiMinusP_ElasticChargeExchange(anEvent);
    break;
  case 22 :
    GeneratePiMinusP_InElastic6(anEvent);
    break;
  case 30 :
    GeneratePPScat(anEvent);
    break;
  case 31 :
    GeneratePiMinusP_KPlusPiLambda(anEvent);
    break;
  case 32 :
    GeneratePiMinusP_K0Lambda(anEvent);
    break;
  case 33 :
    GeneratePiMinusP_K0Lambda_Scat(anEvent);
    break;
  case 34 :
    GenerateKMinusP_KPlusXi(anEvent);
    break;
  case 35 :
    GenerateKMinusP_KPlusXi_Scat(anEvent);
    break;
  case 36 :
    GenerateVPFourBodyDecay(anEvent);
    break;
  case 37 :
    GenerateVPFourBodyDecay_Scat(anEvent);
    break;
  case 38 :
    GenerateGammaP_KPlusLambda(anEvent);
    break;
  case 39 :
    GenerateGammaP_KPlusLambda_Scat(anEvent);
    break;
  case 40 :
   GenerateGammaP_PhiP(anEvent);
    break;
  case 41 :
    GenerateUniform_PiPlus(anEvent);
    break;
  case 42 :
    GenerateUniform_Proton(anEvent);
    break;
  case 43 :
    GenerateGammaP_KPlusSigma(anEvent);
    break;
  case 103 :
    GeneratePiMinus(anEvent);
    break;
  case 95 :
  case 99 :   // Sigma+->pi+n
  case 100 :  // Sigma+->pi0p
    EfficiencyStudyForScatPart(anEvent, ReactionMode_);
    break;
  case 101 :  // Lambda
    EfficiencyStudyForLambda(anEvent);
    break;
  case 96 :
  case 102 :
    GenerateSigmaBeam();
    break;
  default:
    break;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTPrimaryGeneratorAction::GenerateProton(G4Event* anEvent)
{


  G4double phi = (G4double)RandFlat::shoot(0., 360.);
  //G4double theta = (G4double)RandFlat::shoot(10., 90.);
  G4double theta = 45;

  G4ThreeVector localDir(sin(theta*Deg2Rad)*cos(phi*Deg2Rad),
			 sin(theta*Deg2Rad)*sin(phi*Deg2Rad), 
			 cos(theta*Deg2Rad));
  G4ThreeVector localPos(0, 0, 0.*mm);

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();  
  int lnum = 0; // target
  G4ThreeVector globalDir = geomMan.Local2GlobalDir(lnum, localDir);
  G4ThreeVector globalPos = geomMan.Local2GlobalPos(lnum, localPos);
  
  //gunProton_->SetParticleEnergy(PrimaryE_*MeV);
  G4double Ekin = (G4double)RandFlat::shoot(10., 150.);
  /*
  gunProton_->SetParticleEnergy(Ekin*MeV);
  gunProton_->SetParticleMomentumDirection(globalDir);
  gunProton_->SetParticlePosition(globalPos);

  gunProton_->GeneratePrimaryVertex(anEvent);
  */
  gunPiP_->SetParticleEnergy(Ekin*MeV);
  gunPiP_->SetParticleMomentumDirection(globalDir);
  gunPiP_->SetParticlePosition(globalPos);

  gunPiP_->GeneratePrimaryVertex(anEvent);

  // Tempolary
  anaMan_->SetThetaMeson(theta);
  anaMan_->SetPhiMeson(phi);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTPrimaryGeneratorAction::GeneratePiMinus(G4Event* anEvent)
{


  //G4double phi = (G4double)RandFlat::shoot(0., 360.);
  G4double phi = 0.;
  //G4double theta = (G4double)RandFlat::shoot(10., 90.);
  G4double theta = 90;

  G4ThreeVector localDir(sin(theta*Deg2Rad)*cos(phi*Deg2Rad),
			 sin(theta*Deg2Rad)*sin(phi*Deg2Rad), 
			 cos(theta*Deg2Rad));
  G4ThreeVector localPos(0, 0, 0.*mm);

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();  
  int lnum = 0; // target
  G4ThreeVector globalDir = geomMan.Local2GlobalDir(lnum, localDir);
  G4ThreeVector globalPos = geomMan.Local2GlobalPos(lnum, localPos);
  
  //gunProton_->SetParticleEnergy(PrimaryE_*MeV);
  //G4double Ekin = (G4double)RandFlat::shoot(10., 150.);
  G4double Ekin = 150.;
  /*
  gunProton_->SetParticleEnergy(Ekin*MeV);
  gunProton_->SetParticleMomentumDirection(globalDir);
  gunProton_->SetParticlePosition(globalPos);

  gunProton_->GeneratePrimaryVertex(anEvent);
  */
  gunPiM_->SetParticleEnergy(Ekin*MeV);
  gunPiM_->SetParticleMomentumDirection(globalDir);
  gunPiM_->SetParticlePosition(globalPos);

  gunPiM_->GeneratePrimaryVertex(anEvent);

  // Tempolary
  anaMan_->SetThetaMeson(theta);
  anaMan_->SetPhiMeson(phi);


}


void CFTPrimaryGeneratorAction::GenerateScatProtonCheck(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  //G4ParticleDefinition* spiMinus = particleTable->FindParticle("spi-");
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  //G4ParticleDefinition* kPlus = particleTable->FindParticle("kaon+");
  G4ParticleDefinition* kPlus = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* sigma = particleTable->FindParticle("sigma-");
  //G4ParticleDefinition* ssigma = particleTable->FindParticle("ssigma+");
  //G4ParticleDefinition* usigma = particleTable->FindParticle("usigma-");
  //G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  //G4double m_neutron = neutron->GetPDGMass()/GeV;

  double mass_sigma = 1.197;
  double beammom = 1.3;

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  primary_vertex_x = 0.*mm;
  primary_vertex_y = 0.*mm;
  primary_vertex_z = -150.*mm;

  Kinema2Body Sigma(piMinus->GetPDGMass()/GeV, 
		    proton->GetPDGMass()/GeV,
		    kPlus->GetPDGMass()/GeV,
		    mass_sigma);

  Sigma.SetMomentum(1, beammom);

  double thetaCM = GenerateThetaCM_;
  Sigma.SetThetaCM(thetaCM);
  Sigma.calc_kinema();

  double phi = (G4double)RandFlat::shoot(0., 180.);

  /* kaon+ */
  //double Energy_k = Sigma.GetEnergyLab(3);
  double momentum_k = Sigma.GetMomentumLab(3);
  double ThetaK = Sigma.GetThetaLab();
  G4ThreeVector momentumKPlus(momentum_k*sin(ThetaK*Deg2Rad)*cos(phi*Deg2Rad),
			      momentum_k*sin(ThetaK*Deg2Rad)*sin(phi*Deg2Rad),
			      momentum_k*cos(ThetaK*Deg2Rad));
  double PhiK = phi;
  double ThetaKCM = Sigma.GetThetaCM();

  //G4cout << "pi- (" << momentumPiMinus.x() << ", "
  //<< momentumPiMinus.y() << ", " << momentumPiMinus.z() << ") "
  //<< G4endl;

  /* sigma */
  //double Energy_sig = Sigma.GetEnergyLab(4);
  double momentum_sig = Sigma.GetMomentumLab(4);
  double ThetaSig = Sigma.GetPhiLab();

  G4ThreeVector momentumSigma(momentum_sig*sin(ThetaSig*Deg2Rad)*cos((phi+180)*Deg2Rad),
			      momentum_sig*sin(ThetaSig*Deg2Rad)*sin((phi+180)*Deg2Rad),
			      momentum_sig*cos(ThetaSig*Deg2Rad));

  double PhiSig = phi+180.;
  /*
  G4cout << "Sigma (" << momentumSigma.x() << ", "
	 << momentumSigma.y() << ", " << momentumSigma.z() << ") "
	 << G4endl;
  */

  /* pi- beam */
  //double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam( 0., 0., -beammom);

  /*
  PrimaryInfo pInfo;
  pInfo.p     = momentum_k; 
  pInfo.theta = ThetaK;
  pInfo.phi   = PhiK;
  pInfo.x     = primary_vertex_x;
  pInfo.y     = primary_vertex_y;
  pInfo.z     = primary_vertex_z;
  pInfo.pbeam = momentumBeam.mag();
  pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
  pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
  pInfo.mass = kPlus->GetPDGMass()/GeV;
  pInfo.charge = kPlus->GetPDGMass()/GeV;
  anaMan_->PrimaryGeneration( &pInfo );
  */

  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  /* Sigma-p scatt */
  //G4double flength=20.*mm;
  //G4double ctau=44.34; /*mm*/
  G4double p_sigma = momentum_sig;
  G4double m_sigma = sigma->GetPDGMass()/GeV;
  //G4double E_sigma;
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_p= 0.1274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localSigmaPos = localVertexPos;
  G4bool flagScattering=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);
    // Multiple scattering
    /*
    if (flagTgtType == 0 || flagTgtType == 1) {
      momentumSigma = calcMultipleScattering(p_sigma, m_sigma, dx/mm, momentumSigma, flagTgtType);

    }
    */
    totalx += dx;
    localSigmaPos += momentumSigma*dx/momentumSigma.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomSigma ( " << momentumSigma.x() << ", "  
	      << momentumSigma.y() << ", " << momentumSigma.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localSigmaPos.x() << ", "  
	      << localSigmaPos.y() << ", " << localSigmaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumSigma, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagScattering=false;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);

    if (flagTgtType == 0) {
      react_rate = react_rate_p;
      dxInH += dx;
    }

    flagScattering  = scatteringCheck(react_rate, dx/mm);

    p_sigma = calcEnergyDeposit(p_sigma, m_sigma, dx/mm, localSigmaPos, flagTgtType);
    //E_sigma = sqrt(p_sigma*p_sigma + m_sigma*m_sigma);

    /*
    G4cout << "Length : " << totalx << ", p = " << p_sigma*1000. 
	   << " MeV/c, Ekin = "	   << (E_sigma-m_sigma)*1000. 
	   << " MeV" << G4endl;
    */

    if ( flagScattering )
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }

  momentumSigma *= p_sigma/momentumSigma.mag();

  if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localSigmaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance SigmaScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return;
    }

    SigmaScat = Kinema3Resonance(mass_sigma, 
				 scatParticle->GetPDGMass()/GeV,
				 0.0 , mass_sigma,
				 scatParticle->GetPDGMass()/GeV,
				 mass_sigma, 0.0, p_sigma, 0.0,
				 scatDistFlag);

    
    /* scattered Sigma */
    //G4double Energy_scatSigma = SigmaScat.GetEnergy(4);
    //G4double momentum_scatSigma = SigmaScat.GetMomentum(4);
    double ThetaScatHyp = SigmaScat.GetTheta(4);
    SigmaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatSigma( mom[1], mom[2], mom[0]);
    momentumScatSigma.rotateY(ThetaSig*deg);
    momentumScatSigma.rotateZ(PhiSig*deg);
    
    double ThetaScatSigCM, PhiScatSigCM;
    ThetaScatSigCM = SigmaScat.GetThetaCM(1);
    PhiScatSigCM = SigmaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = SigmaScat.GetEnergy(5);
    //G4double momentum_scatPart = SigmaScat.GetMomentum(5);
    SigmaScat.GetMomentum(5,mom);
    double ThetaScatP = SigmaScat.GetTheta(5);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaSig*deg);
    momentumScatPart.rotateZ(PhiSig*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /* 
       G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatPart.x() << ", "
      << momentumScatPart.y() << ", " << momentumScatPart.z() << ") "
      << G4endl;

      G4cout << "Momentum conservation check (x, y, z) = ( " 
	     << momentumSigma.x() - (momentumScatSigma.x()+momentumScatPart.x())
	     << ", "
 	     << momentumSigma.y() - (momentumScatSigma.y()+momentumScatPart.y())
	     << ", "
 	     << momentumSigma.z() - (momentumScatSigma.z()+momentumScatPart.z())
	     << ")" << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatHypMomentum(momentumScatSigma);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    //anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatSigCM);
    anaMan_->SetPhiScatHypCM(PhiScatSigCM);
    anaMan_->SetThetaScatPLab(ThetaScatP);
    anaMan_->SetThetaScatHypLab(ThetaScatHyp);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    /*
    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /* beam */
    /*
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /*sigma*/    
    /*
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /*scatt sigma*/    
    /*
    particleGun->SetParticleDefinition(sigma);
    G4ThreeVector gloMomScatSigma = geomMan.Local2GlobalDir(TgtId, momentumScatSigma);
    particleGun->SetParticleMomentumDirection(gloMomScatSigma);
    particleGun->SetParticleEnergy((Energy_scatSigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  }
}



void CFTPrimaryGeneratorAction::GenerateScatProtonCheck_SigmaPlusP(G4Event* anEvent)
{

  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  //G4ParticleDefinition* spiMinus = particleTable->FindParticle("spi-");
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  //G4ParticleDefinition* kPlus = particleTable->FindParticle("kaon+");
  G4ParticleDefinition* kPlus = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* sigma = particleTable->FindParticle("sigma+");
  //G4ParticleDefinition* ssigma = particleTable->FindParticle("ssigma+");
  //G4ParticleDefinition* usigma = particleTable->FindParticle("usigma+");
  //G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  //G4double m_neutron = neutron->GetPDGMass()/GeV;

  double mass_sigma = 1.189;
  //beammom = 1.4;
  double beammom = 1.3;

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  primary_vertex_x = 0.*mm;
  primary_vertex_y = 0.*mm;
  primary_vertex_z = -150.*mm;

  Kinema2Body Sigma(piMinus->GetPDGMass()/GeV, 
		    proton->GetPDGMass()/GeV,
		    kPlus->GetPDGMass()/GeV,
		    mass_sigma);

  Sigma.SetMomentum(1, beammom);

  double thetaCM = GenerateThetaCM_;
  Sigma.SetThetaCM(thetaCM);
  Sigma.calc_kinema();

  double phi = (G4double)RandFlat::shoot(0., 180.);

  /* kaon+ */
  //double Energy_k = Sigma.GetEnergyLab(3);
  double momentum_k = Sigma.GetMomentumLab(3);
  double ThetaK = Sigma.GetThetaLab();
  G4ThreeVector momentumKPlus(momentum_k*sin(ThetaK*Deg2Rad)*cos(phi*Deg2Rad),
			      momentum_k*sin(ThetaK*Deg2Rad)*sin(phi*Deg2Rad),
			      momentum_k*cos(ThetaK*Deg2Rad));
  double PhiK = phi;
  double ThetaKCM = Sigma.GetThetaCM();

  //G4cout << "pi- (" << momentumPiMinus.x() << ", "
  //<< momentumPiMinus.y() << ", " << momentumPiMinus.z() << ") "
  //<< G4endl;

  /* sigma */
  //double Energy_sig = Sigma.GetEnergyLab(4);
  double momentum_sig = Sigma.GetMomentumLab(4);
  double ThetaSig = Sigma.GetPhiLab();

  G4ThreeVector momentumSigma(momentum_sig*sin(ThetaSig*Deg2Rad)*cos((phi+180)*Deg2Rad),
			      momentum_sig*sin(ThetaSig*Deg2Rad)*sin((phi+180)*Deg2Rad),
			      momentum_sig*cos(ThetaSig*Deg2Rad));

  double PhiSig = phi+180.;
  /*
  G4cout << "Sigma (" << momentumSigma.x() << ", "
	 << momentumSigma.y() << ", " << momentumSigma.z() << ") "
	 << G4endl;
  */

  /* pi- beam */
  //double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam( 0., 0., -beammom);

  /*
  PrimaryInfo pInfo;
  pInfo.p     = momentum_k; 
  pInfo.theta = ThetaK;
  pInfo.phi   = PhiK;
  pInfo.x     = primary_vertex_x;
  pInfo.y     = primary_vertex_y;
  pInfo.z     = primary_vertex_z;
  pInfo.pbeam = momentumBeam.mag();
  pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
  pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
  pInfo.mass = kPlus->GetPDGMass()/GeV;
  pInfo.charge = kPlus->GetPDGMass()/GeV;
  anaMan_->PrimaryGeneration( &pInfo );
  */

  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  /* Sigma-p scatt */
  //G4double flength=20.*mm;
  //G4double ctau=24.04; /*mm*/
  G4double p_sigma = momentum_sig;
  G4double m_sigma = sigma->GetPDGMass()/GeV;
  //G4double E_sigma;
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_p= 0.1274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localSigmaPos = localVertexPos;
  G4bool flagScattering=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);
    // Multiple scattering
    /*
    if (flagTgtType == 0 || flagTgtType == 1) {
      momentumSigma = calcMultipleScattering(p_sigma, m_sigma, dx/mm, momentumSigma, flagTgtType);

    }
    */
    totalx += dx;
    localSigmaPos += momentumSigma*dx/momentumSigma.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomSigma ( " << momentumSigma.x() << ", "  
	      << momentumSigma.y() << ", " << momentumSigma.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localSigmaPos.x() << ", "  
	      << localSigmaPos.y() << ", " << localSigmaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumSigma, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagScattering=false;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);

    if (flagTgtType == 0) {
      react_rate = react_rate_p;
      dxInH += dx;
    }

    flagScattering  = scatteringCheck(react_rate, dx/mm);

    p_sigma = calcEnergyDeposit(p_sigma, m_sigma, dx/mm, localSigmaPos, flagTgtType);
    //E_sigma = sqrt(p_sigma*p_sigma + m_sigma*m_sigma);

    /*
    G4cout << "Length : " << totalx << ", p = " << p_sigma*1000. 
	   << " MeV/c, Ekin = "	   << (E_sigma-m_sigma)*1000. 
	   << " MeV" << G4endl;
    */

    if ( flagScattering )
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }

  momentumSigma *= p_sigma/momentumSigma.mag();

  if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localSigmaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance SigmaScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return;
    }

    SigmaScat = Kinema3Resonance(mass_sigma, 
				 scatParticle->GetPDGMass()/GeV,
				 0.0 , mass_sigma,
				 scatParticle->GetPDGMass()/GeV,
				 mass_sigma, 0.0, p_sigma, 0.0,
				 scatDistFlag);

    
    /* scattered Sigma */
    //G4double Energy_scatSigma = SigmaScat.GetEnergy(4);
    //G4double momentum_scatSigma = SigmaScat.GetMomentum(4);
    double ThetaScatHyp = SigmaScat.GetTheta(4);
    SigmaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatSigma( mom[1], mom[2], mom[0]);
    momentumScatSigma.rotateY(ThetaSig*deg);
    momentumScatSigma.rotateZ(PhiSig*deg);
    
    double ThetaScatSigCM = SigmaScat.GetThetaCM(1);
    double PhiScatSigCM = SigmaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = SigmaScat.GetEnergy(5);
    //G4double momentum_scatPart = SigmaScat.GetMomentum(5);
    SigmaScat.GetMomentum(5,mom);
    double ThetaScatP = SigmaScat.GetTheta(5);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaSig*deg);
    momentumScatPart.rotateZ(PhiSig*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /* 
       G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatPart.x() << ", "
      << momentumScatPart.y() << ", " << momentumScatPart.z() << ") "
      << G4endl;

      G4cout << "Momentum conservation check (x, y, z) = ( " 
	     << momentumSigma.x() - (momentumScatSigma.x()+momentumScatPart.x())
	     << ", "
 	     << momentumSigma.y() - (momentumScatSigma.y()+momentumScatPart.y())
	     << ", "
 	     << momentumSigma.z() - (momentumScatSigma.z()+momentumScatPart.z())
	     << ")" << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatHypMomentum(momentumScatSigma);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    //anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatSigCM);
    anaMan_->SetPhiScatHypCM(PhiScatSigCM);
    anaMan_->SetThetaScatPLab(ThetaScatP);
    anaMan_->SetThetaScatHypLab(ThetaScatHyp);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    /*
    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /* beam */
    /*
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /*sigma*/    
    /*
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /*scatt sigma*/    
    /*
    particleGun->SetParticleDefinition(sigma);
    G4ThreeVector gloMomScatSigma = geomMan.Local2GlobalDir(TgtId, momentumScatSigma);
    particleGun->SetParticleMomentumDirection(gloMomScatSigma);
    particleGun->SetParticleEnergy((Energy_scatSigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  }
}

void CFTPrimaryGeneratorAction::GeneratePiKSigma(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  //G4ParticleDefinition* sigma = particleTable->FindParticle("sigma-");

  double width = 0.0;
  double mass_sigma = 1.197;
  double beammom = 1.32;

  Kinema3Resonance Sigma;
  int DistFlag=3; // 1.3GeV/c (pi-, K+)

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_sigma,
			   kaonPlus->GetPDGMass()/GeV,
			   mass_sigma, width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumKaonPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  /* sigma */
  //double Energy_sig = Sigma.GetEnergy(4);
  //double momentum_sig = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma(mom[1], mom[2], mom[0]);
  //double ThetaSig = Sigma.GetTheta(4);
  //double PhiSig = Sigma.GetPhi(4);

  /* pi+ beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePiKSigmaPlus(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piPlus = particleTable->FindParticle("pi+");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  //G4ParticleDefinition* sigma = particleTable->FindParticle("sigma+");

  double width = 0.0;
  double mass_sigma = 1.189;
  double beammom = 1.42;

  Kinema3Resonance Sigma;
  int DistFlag=2; // (pi+, K+)

  Sigma = Kinema3Resonance(piPlus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_sigma,
			   kaonPlus->GetPDGMass()/GeV,
			   mass_sigma, width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumKaonPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  /* sigma */
  //double Energy_sig = Sigma.GetEnergy(4);
  //double momentum_sig = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma(mom[1], mom[2], mom[0]);
  //double ThetaSig = Sigma.GetTheta(4);
  //double PhiSig = Sigma.GetPhi(4);

  /* pi+ beam */
  double Energy_beam = sqrt(pow(piPlus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<50.)) 
    //if (!(ThetaK>0.&&ThetaK<180.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiKSigmaStar(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  //G4ParticleDefinition* sigma = particleTable->FindParticle("sigma-");

  double width = 39.4/1000.;

  double mass_sigma_star = 1.387;
  double beammom = 1.3;

  Kinema3Resonance Sigma;
  int DistFlag=3; // 1.3GeV/c (pi-, K+)

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   1.116/*Lambda*/ , piMinus->GetPDGMass()/GeV, 
			   kaonPlus->GetPDGMass()/GeV,
			   mass_sigma_star, width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumKaonPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  /* sigma */
  //double Energy_sig = Sigma.GetEnergy(4);
  //double momentum_sig = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma(mom[1], mom[2], mom[0]);
  //double ThetaSig = Sigma.GetTheta(4);
  //double PhiSig = Sigma.GetPhi(4);

  /* pi+ beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePiKSigmaPlusStar(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piPlus = particleTable->FindParticle("pi+");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  //G4ParticleDefinition* sigma = particleTable->FindParticle("sigma+");

  double width = 36.0/1000.;
  double mass_sigma_star = 1.382;
  double beammom = 1.42;
  //beammom = 1.80;

  Kinema3Resonance Sigma;
  int DistFlag=2; // (pi+, K+)

  Sigma = Kinema3Resonance(piPlus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   1.116 /*Lambda*/ , piPlus->GetPDGMass()/GeV, 
			   kaonPlus->GetPDGMass()/GeV,
			   mass_sigma_star, width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumKaonPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  /* sigma */
  //double Energy_sig = Sigma.GetEnergy(4);
  //double momentum_sig = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma(mom[1], mom[2], mom[0]);
  //double ThetaSig = Sigma.GetTheta(4);
  //double PhiSig = Sigma.GetPhi(4);

  /* pi+ beam */
  double Energy_beam = sqrt(pow(piPlus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<50.)) 
    //if (!(ThetaK>0.&&ThetaK<180.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePiMinusP_Elastic(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double width = 0.0;
  double mass_p = 0.938;
  //beammom = 1.3;
  double beammom = (G4double)RandGauss::shoot(1.334, 0.0115);

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  ThreeVector beamMomVec(0, 0, beammom);
  ThreeVector beamPosVec(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  ThreeVector Vertex(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  ThreeVector beamMomVecCor = CorrElossIn(beamMomVec, beamPosVec, Vertex, piMinus->GetPDGMass()/GeV );
  //G4cout << "Vertex ( " << primary_vertex_x << ", " << primary_vertex_y << ", "<< primary_vertex_z << "), Pini = " << beamMomVec.mag() << ", Pcor = " << beamMomVecCor.mag() << G4endl;

  beammom = beamMomVecCor.mag();

  Kinema3Resonance Sigma;
  int DistFlag=5; // 1.32GeV/c pi-p elastic

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_p,
			   piMinus->GetPDGMass()/GeV,
			   mass_p, width, beammom, 0.0, DistFlag);

  /* pi- */
  double Energy_pi = Sigma.GetEnergy(5);
  //double momentum_pi = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiMinus(mom[1], mom[2], mom[0]);
  double ThetaPi = Sigma.GetTheta(5);
  double PhiPi = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* proton */
  double Energy_p = Sigma.GetEnergy(4);
  //double momentum_p = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumProton(mom[1], mom[2], mom[0]);
  //double ThetaP = Sigma.GetTheta(4);
  //double PhiP = Sigma.GetPhi(4);

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  /*
  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;
  */

  if (1) {
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomPiMinus = geomMan.Local2GlobalDir(TgtId, momentumPiMinus);
    particleGun->SetParticleMomentumDirection(gloMomPiMinus);
    particleGun->SetParticleEnergy((Energy_pi - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomProton = geomMan.Local2GlobalDir(TgtId, momentumProton);
    particleGun->SetParticleMomentumDirection(gloMomProton);
    particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi);
    anaMan_->SetPhiMeson(PhiPi);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiMinusP_ElasticChargeExchange(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* piZero = particleTable->FindParticle("pi0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  double width = 0.0;
  //double mass_p = 0.938;
  double beammom = 1.3;

  Kinema3Resonance Sigma;
  int DistFlag=7; // 1.32GeV/c pi-p -> pi0n

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , neutron->GetPDGMass()/GeV,
			   piZero->GetPDGMass()/GeV,
			   neutron->GetPDGMass()/GeV, width, beammom, 0.0, DistFlag);

  /* pi0 */
  double Energy_pi = Sigma.GetEnergy(5);
  //double momentum_pi = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiZero(mom[1], mom[2], mom[0]);
  double ThetaPi = Sigma.GetTheta(5);
  double PhiPi = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* neutron */
  double Energy_n = Sigma.GetEnergy(4);
  //double momentum_n = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumNeutron(mom[1], mom[2], mom[0]);
  //double ThetaN = Sigma.GetTheta(4);
  //double PhiN = Sigma.GetPhi(4);

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  /*
  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;
  */

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(piZero);
    G4ThreeVector gloMomPiZero = geomMan.Local2GlobalDir(TgtId, momentumPiZero);
    particleGun->SetParticleMomentumDirection(gloMomPiZero);
    particleGun->SetParticleEnergy((Energy_pi - piZero->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(neutron);
    G4ThreeVector gloMomNeutron = geomMan.Local2GlobalDir(TgtId, momentumNeutron);
    particleGun->SetParticleMomentumDirection(gloMomNeutron);
    particleGun->SetParticleEnergy((Energy_n - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumNeutron);
    anaMan_->SetScatMesonMomentum(momentumPiZero);
    anaMan_->SetThetaMeson(ThetaPi);
    anaMan_->SetPhiMeson(PhiPi);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiPlusP_Elastic(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piPlus = particleTable->FindParticle("pi+");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double width = 0.0;
  double mass_p = 0.938;
  double beammom = 1.4;

  Kinema3Resonance Sigma;
  int DistFlag=6; // 1.45GeV/c pi+p elastic

  Sigma = Kinema3Resonance(piPlus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_p,
			   piPlus->GetPDGMass()/GeV,
			   mass_p, width, beammom, 0.0, DistFlag);

  /* pi+ */
  double Energy_pi = Sigma.GetEnergy(5);
  //double momentum_pi = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiPlus(mom[1], mom[2], mom[0]);
  double ThetaPi = Sigma.GetTheta(5);
  double PhiPi = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* proton */
  double Energy_p = Sigma.GetEnergy(4);
  //double momentum_p = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumProton(mom[1], mom[2], mom[0]);
  //double ThetaP = Sigma.GetTheta(4);
  //double PhiP = Sigma.GetPhi(4);

  /* pi+ beam */
  double Energy_beam = sqrt(pow(piPlus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  /*
  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;
  */

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomPiPlus = geomMan.Local2GlobalDir(TgtId, momentumPiPlus);
    particleGun->SetParticleMomentumDirection(gloMomPiPlus);
    particleGun->SetParticleEnergy((Energy_pi - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomProton = geomMan.Local2GlobalDir(TgtId, momentumProton);
    particleGun->SetParticleMomentumDirection(gloMomProton);
    particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    /* beam */
    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiPlus);
    anaMan_->SetThetaMeson(ThetaPi);
    anaMan_->SetPhiMeson(PhiPi);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}



void CFTPrimaryGeneratorAction::GeneratePiMinusP_InElastic1(G4Event* anEvent)
/* pi-p --> pi+Delta- --> pi+pi-n*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* piPlus  = particleTable->FindParticle("pi+");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double width = 0.130;
  double mass_delta = 1.236;
  double beammom = 1.05;

  Kinema3Resonance Sigma;
  int DistFlag=0; // flat 

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   piMinus->GetPDGMass()/GeV , neutron->GetPDGMass()/GeV,
			   piPlus->GetPDGMass()/GeV,
			   mass_delta, width, beammom, 0.0, DistFlag);

  /* pi+ */
  double Energy_pi1 = Sigma.GetEnergy(5);
  //double momentum_pi1 = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiPlus(mom[1], mom[2], mom[0]);
  double ThetaPi1 = Sigma.GetTheta(5);
  double PhiPi1 = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* pi- */
  double Energy_pi2 = Sigma.GetEnergy(3);
  //double momentum_pi2 = Sigma.GetMomentum(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumPiMinus(mom[1], mom[2], mom[0]);
  //double ThetaPi2 = Sigma.GetTheta(3);
  //double PhiPi2 = Sigma.GetPhi(3);
  
  /* neutron */
  double Energy_n = Sigma.GetEnergy(4);
  //double momentum_n = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumNeutron(mom[1], mom[2], mom[0]);
  //double ThetaN = Sigma.GetTheta(4);
  //double PhiN = Sigma.GetPhi(4);

  /*
  std::cout << "Momentum x : " << momentumPiPlus.x()+momentumPiMinus.x()+momentumNeutron.x() << std::endl;
  std::cout << "Momentum y : " << momentumPiPlus.y()+momentumPiMinus.y()+momentumNeutron.y() << std::endl;
  std::cout << "Momentum z : " << momentumPiPlus.z()+momentumPiMinus.z()+momentumNeutron.z() << std::endl;
  */

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaPi1>0.&&ThetaPi1<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomPiPlus = geomMan.Local2GlobalDir(TgtId, momentumPiPlus);
    particleGun->SetParticleMomentumDirection(gloMomPiPlus);
    particleGun->SetParticleEnergy((Energy_pi1 - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    bool flagPiPScattering = PiMinusP_Scattering(anEvent, momentumPiMinus, localVertexPos);

    if (!flagPiPScattering) {
      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomPiMinus = geomMan.Local2GlobalDir(TgtId, momentumPiMinus);
      particleGun->SetParticleMomentumDirection(gloMomPiMinus);
      particleGun->SetParticleEnergy((Energy_pi2 - piMinus->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
      anaMan_->SetPiMinusMomentumFromK0(momentumPiMinus);
    }

    bool flagNPScattering = NeutronP_Scattering(anEvent, momentumNeutron, localVertexPos);
    if (!flagNPScattering) {
      particleGun->SetParticleDefinition(neutron);
      G4ThreeVector gloMomNeutron = geomMan.Local2GlobalDir(TgtId, momentumNeutron);
      particleGun->SetParticleMomentumDirection(gloMomNeutron);
      particleGun->SetParticleEnergy((Energy_n - neutron->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    }
    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumNeutron);
    anaMan_->SetScatMesonMomentum(momentumPiPlus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetDecayPiMomentum(momentumPiMinus);
    anaMan_->SetDecayNucleonMomentum(momentumNeutron);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiMinusP_InElastic2(G4Event* anEvent)
/* pi-p --> pi-Delta+ --> pi-pi+n*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* piPlus  = particleTable->FindParticle("pi+");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double width = 0.130;
  double mass_delta = 1.236;
  double beammom = 1.05;

  Kinema3Resonance Sigma;
  int DistFlag=0; // flat 

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   piPlus->GetPDGMass()/GeV , neutron->GetPDGMass()/GeV,
			   piMinus->GetPDGMass()/GeV,
			   mass_delta, width, beammom, 0.0, DistFlag);

  /* pi- */
  double Energy_pi1 = Sigma.GetEnergy(5);
  //double momentum_pi1 = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiMinus(mom[1], mom[2], mom[0]);
  double ThetaPi1 = Sigma.GetTheta(5);
  double PhiPi1 = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* pi+ */
  double Energy_pi2 = Sigma.GetEnergy(3);
  //double momentum_pi2 = Sigma.GetMomentum(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumPiPlus(mom[1], mom[2], mom[0]);
  //double ThetaPi2 = Sigma.GetTheta(3);
  //double PhiPi2 = Sigma.GetPhi(3);

  /* neutron */
  double Energy_n = Sigma.GetEnergy(4);
  //double momentum_n = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumNeutron(mom[1], mom[2], mom[0]);
  //double ThetaN = Sigma.GetTheta(4);
  //double PhiN = Sigma.GetPhi(4);

  /*
  std::cout << "Momentum x : " << momentumPiPlus.x()+momentumPiMinus.x()+momentumNeutron.x() << std::endl;
  std::cout << "Momentum y : " << momentumPiPlus.y()+momentumPiMinus.y()+momentumNeutron.y() << std::endl;
  std::cout << "Momentum z : " << momentumPiPlus.z()+momentumPiMinus.z()+momentumNeutron.z() << std::endl;
  */

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaPi1>0.&&ThetaPi1<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    bool flagPiPScattering = PiMinusP_Scattering(anEvent, momentumPiMinus, localVertexPos);
    if (!flagPiPScattering) {
      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomPiMinus = geomMan.Local2GlobalDir(TgtId, momentumPiMinus);
      particleGun->SetParticleMomentumDirection(gloMomPiMinus);
      particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
      anaMan_->SetPiMinusMomentumFromK0(momentumPiMinus);
    }

    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomPiPlus = geomMan.Local2GlobalDir(TgtId, momentumPiPlus);
    particleGun->SetParticleMomentumDirection(gloMomPiPlus);
    particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    bool flagNPScattering = NeutronP_Scattering(anEvent, momentumNeutron, localVertexPos);

    if (!flagNPScattering) {
      particleGun->SetParticleDefinition(neutron);
      G4ThreeVector gloMomNeutron = geomMan.Local2GlobalDir(TgtId, momentumNeutron);
      particleGun->SetParticleMomentumDirection(gloMomNeutron);
      particleGun->SetParticleEnergy((Energy_n - neutron->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    }

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumNeutron);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetDecayPiMomentum(momentumPiPlus);
    anaMan_->SetDecayNucleonMomentum(momentumNeutron);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiMinusP_InElastic3(G4Event* anEvent)
/* pi-p --> pi-Delta+ --> pi-pi0p*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  //G4ParticleDefinition* piPlus  = particleTable->FindParticle("pi+");
  G4ParticleDefinition* piZero  = particleTable->FindParticle("pi0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double width = 0.130;
  double mass_delta = 1.236;
  double beammom = 1.05;

  Kinema3Resonance Sigma;
  int DistFlag=0; // flat 

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   piZero->GetPDGMass()/GeV , proton->GetPDGMass()/GeV,
			   piMinus->GetPDGMass()/GeV,
			   mass_delta, width, beammom, 0.0, DistFlag);

  /* pi- */
  double Energy_pi1 = Sigma.GetEnergy(5);
  //double momentum_pi1 = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiMinus(mom[1], mom[2], mom[0]);
  double ThetaPi1 = Sigma.GetTheta(5);
  double PhiPi1 = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* pi0 */
  double Energy_pi2 = Sigma.GetEnergy(3);
  //double momentum_pi2 = Sigma.GetMomentum(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumPiZero(mom[1], mom[2], mom[0]);
  //double ThetaPi2 = Sigma.GetTheta(3);
  //double PhiPi2 = Sigma.GetPhi(3);

  /* proton */
  double Energy_p = Sigma.GetEnergy(4);
  //double momentum_p = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumProton(mom[1], mom[2], mom[0]);
  double ThetaP = Sigma.GetTheta(4);
  //double PhiP = Sigma.GetPhi(4);

  /*
  std::cout << "Momentum x : " << momentumPiZero.x()+momentumPiMinus.x()+momentumProton.x() << std::endl;
  std::cout << "Momentum y : " << momentumPiZero.y()+momentumPiMinus.y()+momentumProton.y() << std::endl;
  std::cout << "Momentum z : " << momentumPiZero.z()+momentumPiMinus.z()+momentumProton.z() << std::endl;
  */

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaP>0.&&ThetaP<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomPiMinus = geomMan.Local2GlobalDir(TgtId, momentumPiMinus);
    particleGun->SetParticleMomentumDirection(gloMomPiMinus);
    particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(piZero);
    G4ThreeVector gloMomPiZero = geomMan.Local2GlobalDir(TgtId, momentumPiZero);
    particleGun->SetParticleMomentumDirection(gloMomPiZero);
    particleGun->SetParticleEnergy((Energy_pi2 - piZero->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomProton = geomMan.Local2GlobalDir(TgtId, momentumProton);
    particleGun->SetParticleMomentumDirection(gloMomProton);
    particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetDecayPiMomentum(momentumPiZero);
    anaMan_->SetDecayNucleonMomentum(momentumProton);
    anaMan_->SetPiMinusMomentumFromK0(momentumPiMinus);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiMinusP_InElastic4(G4Event* anEvent)
/* pi-p -->  pi-pi0p (phase space)*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  //G4ParticleDefinition* piPlus  = particleTable->FindParticle("pi+");
  G4ParticleDefinition* piZero  = particleTable->FindParticle("pi0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double beammom = 1.05;

  Kinema3Body Sigma = Kinema3Body(piMinus->GetPDGMass()/GeV, 
		      proton->GetPDGMass()/GeV,
		      piZero->GetPDGMass()/GeV , proton->GetPDGMass()/GeV,
		      piMinus->GetPDGMass()/GeV,
		      beammom, 0.0);

  /* pi- */
  double Energy_pi1 = Sigma.GetEnergy(5);
  //double momentum_pi1 = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiMinus(mom[1], mom[2], mom[0]);
  double ThetaPi1 = Sigma.GetTheta(5);
  double PhiPi1 = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* pi0 */
  double Energy_pi2 = Sigma.GetEnergy(3);
  //double momentum_pi2 = Sigma.GetMomentum(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumPiZero(mom[1], mom[2], mom[0]);
  //double ThetaPi2 = Sigma.GetTheta(3);
  //double PhiPi2 = Sigma.GetPhi(3);

  /* proton */
  double Energy_p = Sigma.GetEnergy(4);
  //double momentum_p = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumProton(mom[1], mom[2], mom[0]);
  double ThetaP = Sigma.GetTheta(4);
  //double PhiP = Sigma.GetPhi(4);

  /*
  std::cout << "Momentum x : " << momentumPiZero.x()+momentumPiMinus.x()+momentumProton.x() << std::endl;
  std::cout << "Momentum y : " << momentumPiZero.y()+momentumPiMinus.y()+momentumProton.y() << std::endl;
  std::cout << "Momentum z : " << momentumPiZero.z()+momentumPiMinus.z()+momentumProton.z() << std::endl;
  */

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaP>0.&&ThetaP<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomPiMinus = geomMan.Local2GlobalDir(TgtId, momentumPiMinus);
    particleGun->SetParticleMomentumDirection(gloMomPiMinus);
    particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(piZero);
    G4ThreeVector gloMomPiZero = geomMan.Local2GlobalDir(TgtId, momentumPiZero);
    particleGun->SetParticleMomentumDirection(gloMomPiZero);
    particleGun->SetParticleEnergy((Energy_pi2 - piZero->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomProton = geomMan.Local2GlobalDir(TgtId, momentumProton);
    particleGun->SetParticleMomentumDirection(gloMomProton);
    particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetDecayPiMomentum(momentumPiZero);
    anaMan_->SetDecayNucleonMomentum(momentumProton);
    anaMan_->SetPiMinusMomentumFromK0(momentumPiMinus);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePiMinusP_InElastic5(G4Event* anEvent)
/* pi-p -->  pi-pi+n (phase space)*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* piPlus  = particleTable->FindParticle("pi+");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  double beammom = 1.05;

  Kinema3Body Sigma = Kinema3Body(piMinus->GetPDGMass()/GeV, 
		      proton->GetPDGMass()/GeV,
		      piPlus->GetPDGMass()/GeV , neutron->GetPDGMass()/GeV,
		      piMinus->GetPDGMass()/GeV,
		      beammom, 0.0);

  /* pi- */
  double Energy_pi1 = Sigma.GetEnergy(5);
  //double momentum_pi1 = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiMinus(mom[1], mom[2], mom[0]);
  double ThetaPi1 = Sigma.GetTheta(5);
  double PhiPi1 = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* pi+ */
  double Energy_pi2 = Sigma.GetEnergy(3);
  //double momentum_pi2 = Sigma.GetMomentum(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumPiPlus(mom[1], mom[2], mom[0]);
  double ThetaPi2 = Sigma.GetTheta(3);
  //double PhiPi2 = Sigma.GetPhi(3);

  /* neutron */
  double Energy_n = Sigma.GetEnergy(4);
  //double momentum_n = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumNeutron(mom[1], mom[2], mom[0]);
  //double ThetaN = Sigma.GetTheta(4);
  //double PhiN = Sigma.GetPhi(4);

  /*
  std::cout << "Momentum x : " << momentumPiPlus.x()+momentumPiMinus.x()+momentumNeutron.x() << std::endl;
  std::cout << "Momentum y : " << momentumPiPlus.y()+momentumPiMinus.y()+momentumNeutron.y() << std::endl;
  std::cout << "Momentum z : " << momentumPiPlus.z()+momentumPiMinus.z()+momentumNeutron.z() << std::endl;
  */

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaPi2>0.&&ThetaPi2<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    G4bool flagPiPScattering = PiMinusP_Scattering(anEvent, momentumPiMinus, localVertexPos);
    if (!flagPiPScattering) {
      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomPiMinus = geomMan.Local2GlobalDir(TgtId, momentumPiMinus);
      particleGun->SetParticleMomentumDirection(gloMomPiMinus);
      particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
      anaMan_->SetPiMinusMomentumFromK0(momentumPiMinus);
    }

    particleGun->SetParticleDefinition(piPlus);
    G4ThreeVector gloMomPiPlus = geomMan.Local2GlobalDir(TgtId, momentumPiPlus);
    particleGun->SetParticleMomentumDirection(gloMomPiPlus);
    particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    G4bool flagNPScattering = NeutronP_Scattering(anEvent, momentumNeutron, localVertexPos);
    if (!flagNPScattering) {
      particleGun->SetParticleDefinition(neutron);
      G4ThreeVector gloMomNeutron = geomMan.Local2GlobalDir(TgtId, momentumNeutron);
      particleGun->SetParticleMomentumDirection(gloMomNeutron);
      particleGun->SetParticleEnergy((Energy_n - neutron->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    }

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumNeutron);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetDecayPiMomentum(momentumPiPlus);
    anaMan_->SetDecayNucleonMomentum(momentumNeutron);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePiMinusP_InElastic6(G4Event* anEvent)
/* pi-p -->  pi0pi0n (phase space)*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* piZero  = particleTable->FindParticle("pi0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  double beammom = 1.05;

  Kinema3Body Sigma = Kinema3Body(piMinus->GetPDGMass()/GeV, 
		      proton->GetPDGMass()/GeV,
		      piZero->GetPDGMass()/GeV , neutron->GetPDGMass()/GeV,
		      piZero->GetPDGMass()/GeV,
		      beammom, 0.0);

  /* pi0 */
  double Energy_pi1 = Sigma.GetEnergy(5);
  //double momentum_pi1 = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumPiZero1(mom[1], mom[2], mom[0]);
  double ThetaPi1 = Sigma.GetTheta(5);
  double PhiPi1 = Sigma.GetPhi(5);
  double ThetaPiCM = 180.-Sigma.GetThetaCM(1);
  double PhiPiCM = Sigma.GetPhiCM(1);

  /* pi0 */
  double Energy_pi2 = Sigma.GetEnergy(3);
  //double momentum_pi2 = Sigma.GetMomentum(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumPiZero2(mom[1], mom[2], mom[0]);
  //double ThetaPi2 = Sigma.GetTheta(3);
  //double PhiPi2 = Sigma.GetPhi(3);

  /* neutron */
  double Energy_n = Sigma.GetEnergy(4);
  //double momentum_n = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumNeutron(mom[1], mom[2], mom[0]);
  //double ThetaN = Sigma.GetTheta(4);
  //double PhiN = Sigma.GetPhi(4);

  /*
  std::cout << "Momentum x : " << momentumPiZero1.x()+momentumPiZero2.x()+momentumNeutron.x() << std::endl;
  std::cout << "Momentum y : " << momentumPiZero1.y()+momentumPiZero2.y()+momentumNeutron.y() << std::endl;
  std::cout << "Momentum z : " << momentumPiZero1.z()+momentumPiZero2.z()+momentumNeutron.z() << std::endl;
  */

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  /*
  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;
  */

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(piZero);
    G4ThreeVector gloMomPiZero1 = geomMan.Local2GlobalDir(TgtId, momentumPiZero1);
    particleGun->SetParticleMomentumDirection(gloMomPiZero1);
    particleGun->SetParticleEnergy((Energy_pi1 - piZero->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(piZero);
    G4ThreeVector gloMomPiZero2 = geomMan.Local2GlobalDir(TgtId, momentumPiZero2);
    particleGun->SetParticleMomentumDirection(gloMomPiZero2);
    particleGun->SetParticleEnergy((Energy_pi2 - piZero->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(neutron);
    G4ThreeVector gloMomNeutron = geomMan.Local2GlobalDir(TgtId, momentumNeutron);
    particleGun->SetParticleMomentumDirection(gloMomNeutron);
    particleGun->SetParticleEnergy((Energy_n - neutron->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumNeutron);
    anaMan_->SetScatMesonMomentum(momentumPiZero1);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetDecayPiMomentum(momentumPiZero2);
    anaMan_->SetDecayNucleonMomentum(momentumNeutron);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePPScat(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double width = 0.0;
  double mass_p = 0.938;
  //double beammom = 0.6;
  //double beammom = (G4double)RandGauss::shoot(0.465, 0.0049);
  //double beammom = (G4double)RandGauss::shoot(0.50, 0.0062);
  //double beammom = (G4double)RandGauss::shoot(0.55, 0.0062);
  //double beammom = (G4double)RandGauss::shoot(0.598, 0.0062);
  //double beammom = (G4double)RandGauss::shoot(0.652, 0.0062);
  //double beammom = (G4double)RandGauss::shoot(0.70, 0.0062);
  double beammom = (G4double)RandGauss::shoot(0.75, 0.0062);


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);

  ThreeVector beamMomVec(0, 0, beammom);
  ThreeVector beamPosVec(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  ThreeVector Vertex(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  ThreeVector beamMomVecCor = CorrElossIn(beamMomVec, beamPosVec, Vertex, mass_p);
  //G4cout << "Vertex ( " << primary_vertex_x << ", " << primary_vertex_y << ", "<< primary_vertex_z << "), Pini = " << beamMomVec.mag() << ", Pcor = " << beamMomVecCor.mag() << G4endl;

  beammom = beamMomVecCor.mag();

  Kinema3Resonance PPScat;
  int DistFlag=0; // flat distribution

  PPScat = Kinema3Resonance(proton->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_p,
			   proton->GetPDGMass()/GeV,
			   mass_p, width, beammom, 0.0, DistFlag);

  /* p1 */
  double Energy_p1 = PPScat.GetEnergy(5);
  //double momentum_p1 = PPScat.GetMomentum(5);
  PPScat.GetMomentum(5,mom);
  G4ThreeVector momentumP1(mom[1], mom[2], mom[0]);
  double ThetaP1 = PPScat.GetTheta(5);
  double PhiP1   = PPScat.GetPhi(5);
  double ThetaP1CM = 180.-PPScat.GetThetaCM(1);
  double PhiP1CM = PPScat.GetPhiCM(1);

  /* recoil proton */
  double Energy_p2 = PPScat.GetEnergy(4);
  //double momentum_p2 = PPScat.GetMomentum(4);
  PPScat.GetMomentum(4,mom);
  G4ThreeVector momentumP2(mom[1], mom[2], mom[0]);
  //double ThetaP2 = PPScat.GetTheta(4);
  //double PhiP2 = PPScat.GetPhi(4);

  /* p beam */
  double Energy_beam = sqrt(pow(proton->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  /*
  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;
  */

  
  if (1) {
    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomP1 = geomMan.Local2GlobalDir(TgtId, momentumP1);
    particleGun->SetParticleMomentumDirection(gloMomP1);
    particleGun->SetParticleEnergy((Energy_p1 - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomP2 = geomMan.Local2GlobalDir(TgtId, momentumP2);
    particleGun->SetParticleMomentumDirection(gloMomP2);
    particleGun->SetParticleEnergy((Energy_p2 - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    /* beam */
    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumP1);
    anaMan_->SetScatMesonMomentum(momentumP1);
    anaMan_->SetThetaMeson(ThetaP1);
    anaMan_->SetPhiMeson(PhiP1);
    anaMan_->SetThetaMesonCM(ThetaP1CM);
    anaMan_->SetPhiMesonCM(PhiP1CM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GeneratePiMinusP_KPlusPiLambda(G4Event* anEvent)
/* pi-p -->  K+pi-Lambda (phase space)*/
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* kaonPlus  = particleTable->FindParticle("kaon+");
  G4ParticleDefinition* lambda  = particleTable->FindParticle("lambda");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double beammom = 1.32;

  Kinema3Body Sigma = Kinema3Body(piMinus->GetPDGMass()/GeV, 
		      proton->GetPDGMass()/GeV,
		      piMinus->GetPDGMass()/GeV , lambda->GetPDGMass()/GeV,
		      kaonPlus->GetPDGMass()/GeV,
		      beammom, 0.0);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumK(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  //double PhiK = Sigma.GetPhi(5);
  //double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  //double PhiKCM = Sigma.GetPhiCM(1);

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    /*
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    */
    anaMan_->SetBeamMomentum(beammom);
    //anaMan_->SetDecayPiMomentum(momentumPiZero);
    //anaMan_->SetDecayNucleonMomentum(momentumProton);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

void CFTPrimaryGeneratorAction::GeneratePiMinusP_K0Lambda(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  //G4ParticleDefinition* kaon0  = particleTable->FindParticle("kaon0S");
  G4ParticleDefinition* kaon0  = particleTable->FindParticle("kaon0");
  G4ParticleDefinition* lambda  = particleTable->FindParticle("lambda");
  //G4ParticleDefinition* lambda  = particleTable->FindParticle("sigma0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double beammom = 1.05;
  double width = 0.0;

  Kinema3Resonance Sigma;
  int DistFlag=8; // 1.05 K0Lambda

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , 
			   lambda->GetPDGMass()/GeV,
			   kaon0->GetPDGMass()/GeV,
			   lambda->GetPDGMass()/GeV,
			   width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumK(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  // Lambda
  double E_lambda = Sigma.GetEnergy(4);
  //double momentum_lambda = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumLambda( mom[1], mom[2], mom[0]);
  //double ThetaLambda = Sigma.GetTheta(4);
  //double PhiLambda = Sigma.GetPhi(4);

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;

  //std::cout << "momentumK      = " << momentumK << std::endl;
  //std::cout << "momentumLambda = " << momentumLambda << std::endl;

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaon0);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    G4bool flagPolarization = false;

    if (!flagPolarization) {
      particleGun->SetParticleDefinition(lambda);
      G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
      particleGun->SetParticleMomentumDirection(gloMomLambda);
      particleGun->SetParticleEnergy((E_lambda - lambda->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    } else {
      G4ThreeVector NormZ = momentumLambda*(1/momentumLambda.mag());
      G4ThreeVector NormY = NormZ.cross(momentumK);
      NormY *= (1./NormY.mag());

      double ctau = 78.9; /*mm*/
      DecayWithPolarization LambdaDecay(lambda->GetPDGMass()/GeV,
					proton->GetPDGMass()/GeV,
					piMinus->GetPDGMass()/GeV,
					//1.0, // Polarization
					0.0, // Polarization
					AsymPara_Lambda, // Asymmetry parameter
					NormY, 
					momentumLambda,
					localVertexPos,
					false, // phi flag, lambda is opposite for x axis
					ctau,
					false // charge flag
					);

      //std::cout << "ThetaLambda = " << ThetaLambda 
      //<< "PhiLambda" << PhiLambda << std::endl;
      G4ThreeVector localDecayPos  = LambdaDecay.GetDecayPos();
      G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

      G4ThreeVector momDecayProton = LambdaDecay.GetMomentum2();

      //std::cout << "momDecayProton = " <<  momDecayProton << std::endl;

      particleGun->SetParticleDefinition(proton);
      G4ThreeVector gloMomDecayProton = geomMan.Local2GlobalDir(TgtId, momDecayProton);
      particleGun->SetParticleMomentumDirection(gloMomDecayProton);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin2()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      G4ThreeVector momDecayPi     = LambdaDecay.GetMomentum3();
      //std::cout << "momDecayPi = " <<  momDecayPi << std::endl;
      //std::cout << "momentumLambda(gene) = " << momDecayProton+momDecayPi << std::endl;
      //std::cout << "momentumLambda(org)  = " << momentumLambda << std::endl;

      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momDecayPi);
      particleGun->SetParticleMomentumDirection(gloMomDecayPi);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin3()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      anaMan_->SetDecayPiMomentum(momDecayPi);
      anaMan_->SetDecayNucleonMomentum(momDecayProton);
      anaMan_->SetDecayPos(localDecayPos);
    }



    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    */
    anaMan_->SetBeamMomentum(beammom);
    //anaMan_->SetDecayPiMomentum(momentumPiZero);
    //anaMan_->SetDecayNucleonMomentum(momentumProton);

    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}

#if 1
void CFTPrimaryGeneratorAction::GeneratePiMinusP_K0Lambda_Scat(G4Event* anEvent)
{
  G4bool flagPolarization = true;
  //G4bool flagPolarization = false;

  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* spiMinus = particleTable->FindParticle("spi-");
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* kaon0  = particleTable->FindParticle("kaon0");
  G4ParticleDefinition* lambda  = particleTable->FindParticle("lambda");  
  G4ParticleDefinition* ulambda  = particleTable->FindParticle("ulambda");
  //G4ParticleDefinition* lambda  = particleTable->FindParticle("sigma0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  G4double m_proton  = proton->GetPDGMass()/GeV;
  G4double m_neutron = neutron->GetPDGMass()/GeV;

  double beammom = 1.05;
  double width = 0.0;

  Kinema3Resonance Sigma;
  int DistFlag=8; // 1.05 K0Lambda

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , 
			   lambda->GetPDGMass()/GeV,
			   kaon0->GetPDGMass()/GeV,
			   lambda->GetPDGMass()/GeV,
			   width, beammom, 0.0, DistFlag);

  /* K0 */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumK(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  // Lambda
  double momentum_lambda = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumLambda( mom[1], mom[2], mom[0]);
  double ThetaLambda = Sigma.GetTheta(4);
  double PhiLambda = Sigma.GetPhi(4);

  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  /*
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  */

  primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
  //primary_vertex_x = (G4double)RandGauss::shoot(0., 8.)*mm;
  //primary_vertex_y = (G4double)RandGauss::shoot(0., 8.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
  primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
					       primary_vertex_y, 
					       primary_vertex_z));
  if (primaryTgtType != 0)
    return;


  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  

  /* Sigma-p scatt */
  //G4double flength=20.*mm;
  G4double ctau=78.9; /*mm*/
  G4double p_lambda = momentum_lambda;
  G4double m_lambda = lambda->GetPDGMass()/GeV;
  G4double E_lambda = sqrt(p_lambda*p_lambda+m_lambda*m_lambda);
  //G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  G4double react_rate_p= 0.01274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;

  //G4double react_rate_sigmaPn = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_sigma0p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  G4double react_rate_sigmaPn = 0.01274; // 1/mm /* 30mb and LH2 target*/
  G4double react_rate_sigma0p = 0.01274; // 1/mm /* 30mb and LH2 target*/

  G4double pth_sigmaPn = 0.6345; 
  G4double pth_sigma0p = 0.6435; 
  

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localLambdaPos = localVertexPos;
  G4bool flagDecay=false;
  G4bool flagScattering=false;
  G4bool flagSigmaPlusN=false;
  G4bool flagSigma0P=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localLambdaPos);
    totalx += dx;
    localLambdaPos += momentumLambda*dx/momentumLambda.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomLambda ( " << momentumLambda.x() << ", "  
	      << momentumLambda.y() << ", " << momentumLambda.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localLambdaPos.x() << ", "  
	      << localLambdaPos.y() << ", " << localLambdaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumLambda, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;
    flagScattering=false;
    flagSigmaPlusN=false;
    flagSigma0P=false;
    
    flagDecay = decayCheck(ctau, p_lambda, m_lambda, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localLambdaPos);

    if (flagTgtType == 0 ) {
      react_rate = react_rate_p;
      dxInH += dx;

      flagScattering  = scatteringCheck(react_rate, dx/mm);
      if (p_lambda > pth_sigmaPn)
	flagSigmaPlusN  = scatteringCheck(react_rate_sigmaPn, dx/mm);
      if (p_lambda > pth_sigma0p)
	flagSigma0P  = scatteringCheck(react_rate_sigma0p, dx/mm);

      //p_sigma = calcEnergyDeposit(p_sigma, m_sigma, dx/mm, localSigmaPos, flagTgtType);
      //E_sigma = sqrt(p_sigma*p_sigma + m_sigma*m_sigma);
    }
    /*
    G4cout << "Length : " << totalx << ", p = " << p_lambda*1000. 
	   << " MeV/c, Ekin = "	   << (E_lambda-m_lambda)*1000. 
	   << " MeV" << G4endl;
    */
    if (fabs(p_lambda)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }
    if (flagDecay || flagScattering || flagSigmaPlusN || flagSigma0P)
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }

  momentumLambda *= p_lambda/momentumLambda.mag();

  if (flagDecay) {
    //std::cout << "decay" << std::endl;
    G4ThreeVector localDecayPos=localLambdaPos;
    G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

    if (1) {
      static int Ndecay=0;
      int decayMode=-1;
      double randval = RandFlat::shoot(0., 1.);
      if (randval>0.639)
	decayMode = 0; // pi0n
      else
	decayMode = 1; // pi-p

      Ndecay++;

      G4ParticleDefinition* decayNucl;
      G4ParticleDefinition* decayPi;
      if (decayMode==1) {
	decayNucl = particleTable->FindParticle("proton_d");
	decayPi = particleTable->FindParticle("pi-");
      } else if (decayMode==0) {
	decayNucl = particleTable->FindParticle("neutron");
	decayPi = particleTable->FindParticle("pi0");
      }

      //neutron scattering
      Kinema3Resonance LambdaDecay;
      LambdaDecay = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				    0.0,
				    0.0 , decayNucl->GetPDGMass()/GeV,
				    decayPi->GetPDGMass()/GeV,
				    decayNucl->GetPDGMass()/GeV,
				    0.0, p_lambda, 0.0);


      /* Neucleon */
      double Energy_decayN = LambdaDecay.GetEnergy(4);
      double momentum_decayN = LambdaDecay.GetMomentum(4);
      LambdaDecay.GetMomentum(4,mom);
      G4ThreeVector momentumDecayN(mom[1], mom[2], mom[0]);
      //double ThetaDecayN = LambdaDecay.GetTheta(4);
      //double PhiDecayN = LambdaDecay.GetPhi(4);
      
      /* pi */
      double Energy_decayPi = LambdaDecay.GetEnergy(5);
      double momentum_decayPi = LambdaDecay.GetMomentum(5);
      LambdaDecay.GetMomentum(5,mom);
      G4ThreeVector momentumDecayPi( mom[1], mom[2], mom[0]);
      //double ThetaDecayPi = LambdaDecay.GetTheta(5);
      //double PhiDecayPi = LambdaDecay.GetPhi(5);
      
      // Neucleon momentum at beam frame
      momentumDecayN.rotateY(ThetaLambda*deg);
      momentumDecayN.rotateZ(PhiLambda*deg);

      // pi momentum at beam frame
      momentumDecayPi.rotateY(ThetaLambda*deg);
      momentumDecayPi.rotateZ(PhiLambda*deg);

      /*
      std::cout << "LambdaMom = ( " << momentumLambda.x() << ", " 
		<< momentumLambda.y() << ", " << momentumLambda.z()
		<< ")" << std::endl;
      std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
		<< momentumDecayN.y() << ", " << momentumDecayN.z()
		<< ")" << std::endl;
      std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
		<< momentumDecayPi.y() << ", " << momentumDecayPi.z()
		<< ")" << std::endl;

      std::cout << "Delta = ( " 
		<< momentumLambda.x()-momentumDecayN.x()-momentumDecayPi.x()
		<< ", " 
		<< momentumLambda.y()-momentumDecayN.y()-momentumDecayPi.y() 
		<< ", " 
		<< momentumLambda.z()-momentumDecayN.z()-momentumDecayPi.z()
		<< ")" << std::endl;
      */
      double ThetaDecayNatBeamFrame, PhiDecayNatBeamFrame; 
      calcThetaPhi(momentumDecayN, 
		   &ThetaDecayNatBeamFrame, 
		   &PhiDecayNatBeamFrame);

      double ThetaDecayPiatBeamFrame, PhiDecayPiatBeamFrame; 
      calcThetaPhi(momentumDecayPi, 
		   &ThetaDecayPiatBeamFrame, 
		   &PhiDecayPiatBeamFrame);
      /*
      std::cout << "ThetaDecayNatBeamFrame: " << ThetaDecayNatBeamFrame
		<< " PhiDecayNatBeamFrame: "  << PhiDecayNatBeamFrame
		<< std::endl;

      std::cout << "ThetaDecayPiatBeamFrame: " << ThetaDecayPiatBeamFrame
		<< " PhiDecayPiatBeamFrame: "  << PhiDecayPiatBeamFrame
		<< std::endl;
      */
      if (PhiDecayNatBeamFrame<0)
	return;
      if (PhiDecayPiatBeamFrame<0)
	return;


      bool flagProtonScattering=false;
      bool flagNeutronScattering=false;
      bool flagPionScattering=false;


      if (decayMode==1) {
	/* pp scatt */
	G4double proton_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
	
	double Ekin_p[14] = {18.2, 19.8, 25.63, 30.14, 39.4, 68.3, 95., 98., 118., 142.,
			     147., 172., 250., 312.};
	
	double pp_cs_table[14] = { 351.8, 314.1, 238.7, 188.4, 138.2, 81.6,
				   56.5, 56.5, 52.7, 52.7, 51.5, 50.2, 50.2, 46.4};
	
	double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_proton*m_proton)
		       -m_proton)*1000. ;// MeV
	
	int index_Ekin=-1;
	double pp_cs;
	for (int i=0; i<14; i++) {
	  if (i==0) {
	    if (Ekin < Ekin_p[i]) {
	      index_Ekin = i;
	      pp_cs = pp_cs_table[i];
	      break;
	    } else if (Ekin >= Ekin_p[i] && Ekin < Ekin_p[i+1]) {
	      index_Ekin = i;
	      pp_cs = (pp_cs_table[i+1]*(Ekin-Ekin_p[i]) + pp_cs_table[i]*(Ekin_p[i+1]-Ekin))/(Ekin_p[i+1]-Ekin_p[i]);
	      break;
	    }
	  } else if (i>=1 && i<=12) {
	    if (Ekin >= Ekin_p[i] && Ekin < Ekin_p[i+1]) {
	      index_Ekin = i;
	      pp_cs = (pp_cs_table[i+1]*(Ekin-Ekin_p[i]) + pp_cs_table[i]*(Ekin_p[i+1]-Ekin))/(Ekin_p[i+1]-Ekin_p[i]);
	      break;
	    }
	  } else if (i==13) {
	    if (Ekin >= Ekin_p[i]) {
	      index_Ekin = i;
	      pp_cs = pp_cs_table[i];
	      break;
	    }
	  }
	}
	
	if (index_Ekin<0 || index_Ekin>=14) {
	  fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
	  exit(-1);
	}
	proton_react_rate *= pp_cs/30.;
	
	G4double dy = 0.5*mm; //mm
	G4double totaly=0.0;  //mm
	
	G4ThreeVector localProtonPos = localDecayPos;
	int flagTgtType2 = -1;
	
	int nIte=0;
	while (1) {
	  flagTgtType2 = -1;
	  totaly += dy;
	  localProtonPos += momentumDecayN*dy/momentumDecayN.mag();
	  
	  flagProtonScattering=false;
	  
	  flagTgtType2 = getTargetFlag(localProtonPos);
	  
	  if (flagTgtType2 == 0 || flagTgtType2 == 1) {
	    flagProtonScattering  = scatteringCheck(proton_react_rate, dy/mm);
	  }
	  
	  if (flagProtonScattering)
	    break;
	  
	  nIte++;
	  if (nIte>1000)
	    break;
	}

	/* pi-p scatt */
	G4double pion_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
	//G4double pion_react_rate = 0.0002548; // 1/mm /* 60mb and LH2 target*/
	
	double p_pi_table[8] = {0.1479, 0.1738, 0.1883, 0.2123, 0.2379, 0.2712, 0.2983, 0.3228};

	double pip_cs_table[8] = { 3.4, 5.12, 6.64, 9.8, 15.2, 19.4, 18.3, 14.8};
	double p = momentum_decayPi;
	int index_p=-1;
	double pip_cs;
	for (int i=0; i<8; i++) {
	  if (i==0) {
	    if (p < p_pi_table[i]) {
	      index_p = i;
	      pip_cs = pip_cs_table[i];
	      break;
	    } else if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	      index_p = i;
	      pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	      break;
	    }
	  } else if (i>=1 && i<=6) {
	    if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	      index_p = i;
	      pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	      break;
	    }
	  } else if (i==7) {
	    if (p >= p_pi_table[i]) {
	      index_p = i;
	      pip_cs = pip_cs_table[i];
	      break;
	    }
	  }
	}

	if (index_p<0 || index_p>=8) {
	  fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_p : %d", index_p);
	  exit(-1);
	}
	pion_react_rate *= pip_cs/30.;
	
	G4double dz = 0.5*mm; //mm
	G4double totalz=0.0;  //mm

	G4double ctau_pi=7804.5; /*mm*/
	G4double p_pi = momentumDecayPi.mag();
	G4double m_pi = piMinus->GetPDGMass()/GeV;

	G4ThreeVector localPionPos = localDecayPos;
	int flagTgtType3 = -1;

	flagPionScattering=false;
	G4bool flagPionDecay=false;

	nIte=0;
	while (1) {
	  flagTgtType3 = -1;
	  totalz += dz;
	  localPionPos += momentumDecayPi*dz/momentumDecayPi.mag();
	  
	  flagPionScattering=false;
	  flagPionDecay=false;

	  flagTgtType3 = getTargetFlag(localPionPos);

	  if (flagTgtType3 == 0) {
	    flagPionScattering  = scatteringCheck(pion_react_rate, dz/mm);
	  }
	  
	  if (flagPionScattering)
	    break;

	  flagPionDecay = decayCheck(ctau_pi, p_pi, m_pi, dz/mm );
	  if (flagPionDecay)
	    break;

	  nIte++;
	  if (nIte>1000)
	    break;
	}


	if (flagProtonScattering) {
	  //std::cout << "*********************** pp scat ***********************" << std::endl;
	  G4ThreeVector localNNScatPos=localProtonPos;
	  G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
	  int scatDistFlag = 0;
	  G4ParticleDefinition* scatParticle;
	  if ( flagTgtType2 == 0 ) {
	    // inside the LH2 target;
	    scatParticle = particleTable->FindParticle("proton");
	    scatDistFlag = 1;  // flat distribution
	  } else {
	    return;
	  }
	  //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	  Kinema3Resonance NNScat;
	  NNScat = Kinema3Resonance(decayNucl->GetPDGMass()/GeV, 
				    scatParticle->GetPDGMass()/GeV,
				    0.0 , decayNucl->GetPDGMass()/GeV,
				    scatParticle->GetPDGMass()/GeV,
				    decayNucl->GetPDGMass()/GeV,
				    0.0, momentum_decayN, 0.0, scatDistFlag);
	  
	  /* scattered Necleon */
	  G4double Energy_scatN = NNScat.GetEnergy(4);
	  //G4double momentum_scatN = NNScat.GetMomentum(4);
	  NNScat.GetMomentum(4,mom);
	  G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);
	  
	  momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /* scattered Proton or Deuteron*/
	  G4double Energy_scatPart = NNScat.GetEnergy(5);
	  //G4double momentum_scatPart = NNScat.GetMomentum(5);
	  NNScat.GetMomentum(5,mom);
	  G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	  momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /*
	    std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
	    << momentumDecayN.y() << ", " << momentumDecayN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
	    << momentumScatN.y() << ", " << momentumScatN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
	    << momentumScatPart.y() << ", " << momentumScatPart.z()
	    << ")" << std::endl;
	    
	    std::cout << "Delta = ( " 
	    << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
	    << ", " 
	    << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
	    << ", " 
	    << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
	    << ")" << std::endl;
	  */
	  
	  particleGun->SetParticleDefinition(kaon0);
	  G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
	  particleGun->SetParticleMomentumDirection(gloMomK);
	  particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* beam */
	  particleGun->SetParticleDefinition(spiMinus);
	  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	  particleGun->SetParticleMomentumDirection(gloMomBeam);
	  particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*sigma*/    
	  /*
	  particleGun->SetParticleDefinition(ssigma);
	  G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	  particleGun->SetParticleMomentumDirection(-gloMomSigma);
	  // Energy_sig is original Energy without energy deposit
	  particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  */	  

	  /* decay pi */
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	  particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	  particleGun->SetParticleEnergy((Energy_decayPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scat proton*/    
	  //G4cout << "pp scat" << G4endl;

	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
	  particleGun->SetParticleMomentumDirection(gloMomentumScatN);
	  particleGun->SetParticleEnergy((Energy_scatN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scatt proton*/    
	  particleGun->SetParticleDefinition(scatParticle);
	  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	  particleGun->SetParticleMomentumDirection(gloMomScatPart);
	  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  
	  anaMan_->SetPrimaryVertex(localVertexPos);
	  anaMan_->SetHypBeamMomentum(momentumLambda);
	  anaMan_->SetScatMesonMomentum(momentumK);
	  anaMan_->SetThetaMeson(ThetaK);
	  anaMan_->SetPhiMeson(PhiK);
	  anaMan_->SetThetaMesonCM(ThetaKCM);
	  anaMan_->SetPhiMesonCM(PhiKCM);
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetNNScatPos(localNNScatPos);
	  anaMan_->SetDecayFlag(decayMode);
	  anaMan_->SetNNScatFlag();
	  anaMan_->SetNNScatTarget(flagTgtType2);	
	  anaMan_->SetScatProtonMomentum(momentumScatPart);
	  anaMan_->SetFlightLengthInTarget(dxInH, 0);
	  anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  anaMan_->SetBeamMomentum(beammom);
	  
	  return;
	}

	if (flagPionScattering) {
	  //std::cout << "piNscat" << std::endl;
	  G4ThreeVector localPiNScatPos=localPionPos;
	  G4ThreeVector globalPiNScatPos = geomMan.Local2GlobalPos(TgtId, localPiNScatPos);
	  int scatDistFlag = 0;
	  G4ParticleDefinition* scatParticle;
	  if ( flagTgtType3 == 0 ) {
	    // inside the LH2 target;
	    scatParticle = particleTable->FindParticle("proton");
	    scatDistFlag = 0;
	  } else {
	    return;
	  }
	  //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	  Kinema3Resonance PiNScat;
	  PiNScat = Kinema3Resonance(decayPi->GetPDGMass()/GeV, 
				     scatParticle->GetPDGMass()/GeV,
				     0.0 , decayPi->GetPDGMass()/GeV,
				     scatParticle->GetPDGMass()/GeV,
				     decayPi->GetPDGMass()/GeV,
				     0.0, momentum_decayPi, 0.0, scatDistFlag);
	  
	  /* scattered pion */
	  G4double Energy_scatPi = PiNScat.GetEnergy(4);
	  //G4double momentum_scatPi = PiNScat.GetMomentum(4);
	  PiNScat.GetMomentum(4,mom);
	  G4ThreeVector momentumScatPi( mom[1], mom[2], mom[0]);
	  
	  momentumScatPi.rotateY(ThetaDecayPiatBeamFrame*deg);
	  momentumScatPi.rotateZ(PhiDecayPiatBeamFrame*deg);
	  
	  /* scattered Proton or Deuteron*/
	  G4double Energy_scatPart = PiNScat.GetEnergy(5);
	  //G4double momentum_scatPart = PiNScat.GetMomentum(5);
	  PiNScat.GetMomentum(5,mom);
	  G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	  momentumScatPart.rotateY(ThetaDecayPiatBeamFrame*deg);
	  momentumScatPart.rotateZ(PhiDecayPiatBeamFrame*deg);
	  
	  /*
	    std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
	    << momentumDecayPi.y() << ", " << momentumDecayPi.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPi = ( " << momentumScatPi.x() << ", " 
	    << momentumScatPi.y() << ", " << momentumScatPi.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
	    << momentumScatPart.y() << ", " << momentumScatPart.z()
	    << ")" << std::endl;
	    
	    std::cout << "Delta = ( " 
	    << momentumDecayPi.x()-momentumScatPi.x()-momentumScatPart.x()
	    << ", " 
	    << momentumDecayPi.y()-momentumScatPi.y()-momentumScatPart.y()
	    << ", " 
	    << momentumDecayPi.z()-momentumScatPi.z()-momentumScatPart.z()
	    << ")" << std::endl;
	  */
	  
	  particleGun->SetParticleDefinition(kaon0);
	  G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
	  particleGun->SetParticleMomentumDirection(gloMomK);
	  particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* beam */
	  particleGun->SetParticleDefinition(spiMinus);
	  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	  particleGun->SetParticleMomentumDirection(gloMomBeam);
	  particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*sigma*/    
	  /*
	  particleGun->SetParticleDefinition(ssigma);
	  G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	  particleGun->SetParticleMomentumDirection(-gloMomSigma);
	  // Energy_sig is original Energy without energy deposit
	  particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  */
	  /* decay neucleon */
	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomDecayN = geomMan.Local2GlobalDir(TgtId, momentumDecayN);
	  particleGun->SetParticleMomentumDirection(gloMomDecayN);
	  particleGun->SetParticleEnergy((Energy_decayN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scat pi*/    
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomentumScatPi = geomMan.Local2GlobalDir(TgtId, momentumScatPi);
	  particleGun->SetParticleMomentumDirection(gloMomentumScatPi);
	  particleGun->SetParticleEnergy((Energy_scatPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalPiNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scatt proton*/    
	  particleGun->SetParticleDefinition(scatParticle);
	  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	  particleGun->SetParticleMomentumDirection(gloMomScatPart);
	  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalPiNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  
	  anaMan_->SetPrimaryVertex(localVertexPos);
	  anaMan_->SetHypBeamMomentum(momentumLambda);
	  anaMan_->SetScatMesonMomentum(momentumK);
	  anaMan_->SetThetaMeson(ThetaK);
	  anaMan_->SetPhiMeson(PhiK);
	  anaMan_->SetThetaMesonCM(ThetaKCM);
	  anaMan_->SetPhiMesonCM(PhiKCM);
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetPiNScatPos(localPiNScatPos);
	  anaMan_->SetDecayFlag(decayMode);
	  anaMan_->SetPiNScatFlag();
	  anaMan_->SetPiNScatTarget(flagTgtType3);	
	  anaMan_->SetScatProtonMomentum(momentumScatPart);
	  anaMan_->SetFlightLengthInTarget(dxInH, 0);
	  //anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayPiMomentum(momentumScatPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  anaMan_->SetBeamMomentum(beammom);
	  
	  return;
	}
      } else if (decayMode==0) {
	/* n-p scatt */
	G4double neutron_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
	
	double Ekin_n[30] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
			     110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
		   210., 220., 230., 240., 250., 260., 270., 280., 290., 300.};
	
	double np_cs_table[30] = { 951.827, 488.477, 312.839, 223.243, 170.507,
				   136.696, 113.728,  97.449,  85.517,  76.527,
				   69.593,  64.133,  59.756,  56.191,  53.244,
				   50.779,  48.691,  46.905,  45.364,  44.022,
				   42.844,  41.803,  40.876,  40.046,  39.298,
				   38.620,  38.002,  37.436,  36.915,  36.432};
	
	double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_neutron*m_neutron)
		       -m_neutron)*1000. ;// MeV
	
	
	int index_Ekin=-1;
	double np_cs;
	for (int i=0; i<30; i++) {
	  if (i==0) {
	    if (Ekin < Ekin_n[i]) {
	      index_Ekin = i;
	      np_cs = np_cs_table[i];
	      break;
	    } else if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	      index_Ekin = i;
	      np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	      break;
	    }
	  } else if (i>=1 && i<=28) {
	    if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	      index_Ekin = i;
	      np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	      break;
	    }
	  } else if (i==29) {
	    if (Ekin >= Ekin_n[i]) {
	      index_Ekin = i;
	      np_cs = np_cs_table[i];
	      break;
	    }
	  }
	}
	
	if (index_Ekin<0 || index_Ekin>=30) {
	  fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
	  exit(-1);
	}
	neutron_react_rate *= np_cs/30.;
	
	G4double dy = 0.5*mm; //mm
	G4double totaly=0.0;  //mm
	
	G4ThreeVector localNeutronPos = localDecayPos;
	int flagTgtType2 = -1;
	
	int nIte=0;
	while (1) {
	  flagTgtType2 = -1;
	  totaly += dy;
	  localNeutronPos += momentumDecayN*dy/momentumDecayN.mag();
	  
	  flagNeutronScattering=false;
	  
	  flagTgtType2 = getTargetFlag(localNeutronPos);
	  
	  if (flagTgtType2 == 0 || flagTgtType2 == 1) {
	    flagNeutronScattering  = scatteringCheck(neutron_react_rate, dy/mm);
	  }
	  
	  if (flagNeutronScattering)
	    break;
	  
	  nIte++;
	  if (nIte>1000)
	    break;
	}

	if (flagNeutronScattering) {
	  //std::cout << "NNscat" << std::endl;
	  G4ThreeVector localNNScatPos=localNeutronPos;
	  G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
	  int scatDistFlag = 0;
	  G4ParticleDefinition* scatParticle;
	  if ( flagTgtType2 == 0 ) {
	    // inside the LH2 target;
	    scatParticle = particleTable->FindParticle("proton");
	    scatDistFlag = 100+index_Ekin;
	  } else {
	    return;
	  }
	  //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	  Kinema3Resonance NNScat;
	  NNScat = Kinema3Resonance(decayNucl->GetPDGMass()/GeV, 
				    scatParticle->GetPDGMass()/GeV,
				    0.0 , decayNucl->GetPDGMass()/GeV,
				    scatParticle->GetPDGMass()/GeV,
				    decayNucl->GetPDGMass()/GeV,
				    0.0, momentum_decayN, 0.0, scatDistFlag);
	  
	  /* scattered Necleon */
	  G4double Energy_scatN = NNScat.GetEnergy(4);
	  //G4double momentum_scatN = NNScat.GetMomentum(4);
	  NNScat.GetMomentum(4,mom);
	  G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);
	  
	  momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /* scattered Proton or Deuteron*/
	  G4double Energy_scatPart = NNScat.GetEnergy(5);
	  //G4double momentum_scatPart = NNScat.GetMomentum(5);
	  NNScat.GetMomentum(5,mom);
	  G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	  momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /*
	    std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
	    << momentumDecayN.y() << ", " << momentumDecayN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
	    << momentumScatN.y() << ", " << momentumScatN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
	    << momentumScatPart.y() << ", " << momentumScatPart.z()
	    << ")" << std::endl;
	    
	    std::cout << "Delta = ( " 
	    << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
	    << ", " 
	    << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
	    << ", " 
	    << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
	    << ")" << std::endl;
	  */
	  
	  particleGun->SetParticleDefinition(kaon0);
	  G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
	  particleGun->SetParticleMomentumDirection(gloMomK);
	  particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* beam */
	  particleGun->SetParticleDefinition(spiMinus);
	  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	  particleGun->SetParticleMomentumDirection(gloMomBeam);
	  particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*sigma*/    
	  /*
	  particleGun->SetParticleDefinition(ssigma);
	  G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	  particleGun->SetParticleMomentumDirection(-gloMomSigma);
	  // Energy_sig is original Energy without energy deposit
	  particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  */
	  
	  /* decay pi */
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	  particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	  particleGun->SetParticleEnergy((Energy_decayPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scat proton*/    
	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
	  particleGun->SetParticleMomentumDirection(gloMomentumScatN);
	  particleGun->SetParticleEnergy((Energy_scatN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scatt proton*/    
	  particleGun->SetParticleDefinition(scatParticle);
	  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	  particleGun->SetParticleMomentumDirection(gloMomScatPart);
	  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  
	  anaMan_->SetPrimaryVertex(localVertexPos);
	  anaMan_->SetHypBeamMomentum(momentumLambda);
	  anaMan_->SetScatMesonMomentum(momentumK);
	  anaMan_->SetThetaMeson(ThetaK);
	  anaMan_->SetPhiMeson(PhiK);
	  anaMan_->SetThetaMesonCM(ThetaKCM);
	  anaMan_->SetPhiMesonCM(PhiKCM);
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetNNScatPos(localNNScatPos);
	  anaMan_->SetDecayFlag(decayMode);
	  anaMan_->SetNNScatFlag();
	  anaMan_->SetNNScatTarget(flagTgtType2);	
	  anaMan_->SetScatProtonMomentum(momentumScatPart);
	  anaMan_->SetFlightLengthInTarget(dxInH, 0);
	  anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  //anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  anaMan_->SetDecayNucleonMomentum(momentumScatN);
	  anaMan_->SetBeamMomentum(beammom);
	  
	  return;
	}

      }
    }

    static int nDecay = 0;
    const int prescale = 100;
    //const int prescale = 1;
    nDecay++;

    if (nDecay%prescale != 0)
      return;


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetDecayPos(localLambdaPos);
    anaMan_->SetDecayEventFlag();
    anaMan_->SetDecayFlag();
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kaon0);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*decaying lambda*/    
    if (!flagPolarization) {
      particleGun->SetParticleDefinition(ulambda);
      
      G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
      particleGun->SetParticleMomentumDirection(gloMomLambda);
      particleGun->SetParticleEnergy((E_lambda - lambda->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    }  else {
      G4ThreeVector NormZ = momentumLambda*(1/momentumLambda.mag());
      G4ThreeVector NormY = NormZ.cross(momentumK);
      NormY *= (1./NormY.mag());

      double ctau_ulambda = 0.; /*mm*/
      DecayWithPolarization LambdaDecay(lambda->GetPDGMass()/GeV,
					proton->GetPDGMass()/GeV,
					piMinus->GetPDGMass()/GeV,
					1.0, // Polarization
					AsymPara_Lambda, // Asymmetry parameter
					NormY, 
					momentumLambda,
					localDecayPos,
					false, // phi flag, lambda is opposite for x axis
					ctau_ulambda,
					false // charge flag
					);

      //std::cout << "ThetaLambda = " << ThetaLambda 
      //<< "PhiLambda" << PhiLambda << std::endl;
      G4ThreeVector momDecayProton = LambdaDecay.GetMomentum2();

      //std::cout << "momDecayProton = " <<  momDecayProton << std::endl;

      particleGun->SetParticleDefinition(proton);
      G4ThreeVector gloMomDecayProton = geomMan.Local2GlobalDir(TgtId, momDecayProton);
      particleGun->SetParticleMomentumDirection(gloMomDecayProton);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin2()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      G4ThreeVector momDecayPi     = LambdaDecay.GetMomentum3();
      //std::cout << "momDecayPi = " <<  momDecayPi << std::endl;
      //std::cout << "momentumLambda(gene) = " << momDecayProton+momDecayPi << std::endl;
      //std::cout << "momentumLambda(org)  = " << momentumLambda << std::endl;

      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momDecayPi);
      particleGun->SetParticleMomentumDirection(gloMomDecayPi);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin3()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      anaMan_->SetDecayPiMomentum(momDecayPi);
      anaMan_->SetDecayNucleonMomentum(momDecayProton);
      anaMan_->SetDecayPos(localDecayPos);
    }

  } else if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localLambdaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance LambdaScat;
    G4ParticleDefinition* scatParticle;
    //int scatDistFlag = 10; // 10 " : 90 < phi < 270, theta : flat
    int scatDistFlag = 0; // 10 " : 0 < phi < 360, theta : flat
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
    } else {
      return;
    }

    LambdaScat = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , 
				  lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  lambda->GetPDGMass()/GeV,
				  0.0, p_lambda, 0.0,
				  scatDistFlag);
    
    
    /* scattered lambda */
    G4double Energy_scatLambda = LambdaScat.GetEnergy(4);
    //G4double momentum_scatLambda = LambdaScat.GetMomentum(4);
    LambdaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatLambda( mom[1], mom[2], mom[0]);
    momentumScatLambda.rotateY(ThetaLambda*deg);
    momentumScatLambda.rotateZ(PhiLambda*deg);
    
    double ThetaScatLambdaCM = LambdaScat.GetThetaCM(1);
    double PhiScatLambdaCM   = LambdaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = LambdaScat.GetEnergy(5);
    //G4double momentum_scatPart = LambdaScat.GetMomentum(5);
    LambdaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaLambda*deg);
    momentumScatPart.rotateZ(PhiLambda*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Lambda (" << momentumLambda.x() << ", "
      << momentumLambda.y() << ", " << momentumLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatHypMomentum(momentumScatLambda);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatLambdaCM);
    anaMan_->SetPhiScatHypCM(PhiScatLambdaCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kaon0);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt lambda*/    
    if (!flagPolarization) {
      G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
      particleGun->SetParticleDefinition(lambda);
      G4ThreeVector gloMomScatLambda = geomMan.Local2GlobalDir(TgtId, momentumScatLambda);
      particleGun->SetParticleMomentumDirection(gloMomScatLambda);
      particleGun->SetParticleEnergy((Energy_scatLambda - lambda->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalScatPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    } else {
      G4ThreeVector NormZ = momentumLambda*(1/momentumLambda.mag());
      G4ThreeVector NormY = NormZ.cross(momentumScatLambda);
      NormY *= (1./NormY.mag());


      double ctau_lambda = 78.9; /*mm*/
      DecayWithPolarization LambdaDecay(lambda->GetPDGMass()/GeV,
					proton->GetPDGMass()/GeV,
					piMinus->GetPDGMass()/GeV,
					1.0, // Polarization
					AsymPara_Lambda, // Asymmetry parameter
					NormY,
					momentumScatLambda,
					localScatPos,
					true, // phi flag, lambda is same with x axis
					ctau_lambda,
					false // charge flag
					);

      //std::cout << "ThetaLambda = " << ThetaLambda 
      //<< "PhiLambda" << PhiLambda << std::endl;
      G4ThreeVector localDecayPos  = LambdaDecay.GetDecayPos();
      G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

      G4ThreeVector momDecayProton = LambdaDecay.GetMomentum2();

      //std::cout << "momDecayProton = " <<  momDecayProton << std::endl;

      particleGun->SetParticleDefinition(proton);
      G4ThreeVector gloMomDecayProton = geomMan.Local2GlobalDir(TgtId, momDecayProton);
      particleGun->SetParticleMomentumDirection(gloMomDecayProton);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin2()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      G4ThreeVector momDecayPi     = LambdaDecay.GetMomentum3();
      //std::cout << "momDecayPi = " <<  momDecayPi << std::endl;
      //std::cout << "momentumLambda(gene) = " << momDecayProton+momDecayPi << std::endl;
      //std::cout << "momentumLambda(org)  = " << momentumLambda << std::endl;

      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momDecayPi);
      particleGun->SetParticleMomentumDirection(gloMomDecayPi);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin3()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      anaMan_->SetDecayPiMomentum(momDecayPi);
      anaMan_->SetDecayNucleonMomentum(momDecayProton);
      anaMan_->SetDecayPos(localDecayPos);
    }
    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  } else if (flagSigmaPlusN || flagSigma0P) {
    //std::cout << "conversion scat" << std::endl;
    G4ThreeVector localScatPos=localLambdaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance LambdaScat;
    G4ParticleDefinition* scatParticle;
    //int scatDistFlag = 10; // 10 " : 90 < phi < 270, theta : flat
    int scatDistFlag = 0; // 10 " : 0 < phi < 360, theta : flat
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
    } else {
      return;
    }
    G4ParticleDefinition* scatHyperon = particleTable->FindParticle("sigma+");
    G4ParticleDefinition* recoilNucleon = particleTable->FindParticle("neutron");      
    int reactMode=2;

    if (flagSigma0P) {
      scatHyperon = particleTable->FindParticle("sigma0");
      recoilNucleon = particleTable->FindParticle("proton");      
      reactMode=3;
    }


    LambdaScat = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , 
				  scatHyperon->GetPDGMass()/GeV,
				  recoilNucleon->GetPDGMass()/GeV,
				  scatHyperon->GetPDGMass()/GeV,
				  0.0, p_lambda, 0.0,
				  scatDistFlag);
    
    
    /* scattered Hyperon */
    G4double Energy_scatHyperon = LambdaScat.GetEnergy(4);
    //G4double momentum_scatHyperon = LambdaScat.GetMomentum(4);
    LambdaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatHyperon( mom[1], mom[2], mom[0]);
    momentumScatHyperon.rotateY(ThetaLambda*deg);
    momentumScatHyperon.rotateZ(PhiLambda*deg);
    
    double ThetaScatHyperonCM = LambdaScat.GetThetaCM(1);
    double PhiScatHyperonCM   = LambdaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = LambdaScat.GetEnergy(5);
    //G4double momentum_scatPart = LambdaScat.GetMomentum(5);
    LambdaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaLambda*deg);
    momentumScatPart.rotateZ(PhiLambda*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Lambda (" << momentumLambda.x() << ", "
      << momentumLambda.y() << ", " << momentumLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    //std::cout << "reactMomde = "  << reactMode << std::endl;

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatHypMomentum(momentumScatHyperon);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatHyperonCM);
    anaMan_->SetPhiScatHypCM(PhiScatHyperonCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag(reactMode);
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kaon0);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaon0->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt lambda*/    
    G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
    particleGun->SetParticleDefinition(scatHyperon);
    G4ThreeVector gloMomScatHyperon = geomMan.Local2GlobalDir(TgtId, momentumScatHyperon);
    particleGun->SetParticleMomentumDirection(gloMomScatHyperon);
    particleGun->SetParticleEnergy((Energy_scatHyperon - scatHyperon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(recoilNucleon);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - recoilNucleon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  }


}

#endif


void CFTPrimaryGeneratorAction::GenerateVPFourBodyDecay(G4Event* anEvent)
{
  // Basic Settings //

  //G4cout << "*-----START CFTPrimaryGeneratorAction::GenerateVPFourBodyDecay-----*" << G4endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition*  VP_PiMinusP = particleTable->FindParticle("VP_PiMinusP_1.05");

  // Mass
  double VPMass = VP_PiMinusP->GetPDGMass()/GeV;

  double beammom = 1.05;
  // VP_PiMinusP Gun Settings
  double Energy_VP = sqrt(pow(VPMass, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumVP(0., 0., beammom);

  G4ParticleDefinition* piMinus = particleTable->FindParticle("spi-");
  /* pi- beam */
  G4double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(VP_PiMinusP);
    G4ThreeVector localMomVP(0.,0.,1.);
    G4ThreeVector gloMomVP = geomMan.Local2GlobalDir(TgtId, localMomVP);
    particleGun->SetParticleMomentumDirection(gloMomVP);
    particleGun->SetParticleEnergy((Energy_VP - VPMass)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetBeamMomentum(beammom);
    /*
    anaMan_->SetHypBeamMomentum(momentumProton);
    anaMan_->SetScatMesonMomentum(momentumPiMinus);
    anaMan_->SetThetaMeson(ThetaPi1);
    anaMan_->SetPhiMeson(PhiPi1);
    anaMan_->SetThetaMesonCM(ThetaPiCM);
    anaMan_->SetPhiMesonCM(PhiPiCM);
    */
    anaMan_->SetBeamMomentum(beammom);
    //anaMan_->SetDecayPiMomentum(momentumPiZero);
    //anaMan_->SetDecayNucleonMomentum(momentumProton);

  }


}

void CFTPrimaryGeneratorAction::GenerateVPFourBodyDecay_Scat(G4Event* anEvent)
{
  // Basic Settings //
  //G4cout << "*-----START CFTPrimaryGeneratorAction::GenerateVPFourBodyDecay-----*" << G4endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  static int GenerateNum=0;

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  
  if (GenerateNum%2 == 0) {
    anaMan_->BeginOfPrimaryAction();
    anaMan_->SetPrimaryFlag(true);

    G4ParticleDefinition* VP_PiMinusP = particleTable->FindParticle("VP_PiMinusP_1.05");

    // Mass
    double VPMass = VP_PiMinusP->GetPDGMass()/GeV;
    double beammom = 1.05;
    // VP_PiMinusP Gun Settings
    double Energy_VP = sqrt(pow(VPMass, 2.0) + pow(beammom, 2.0));
    G4ThreeVector momentumVP(0., 0., beammom);
    
    //G4ParticleDefinition* piMinus = particleTable->FindParticle("spi-");

    /* pi- beam */
    //G4double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
    G4ThreeVector momentumBeam(0., 0., -beammom);
    
    //G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
    G4ThreeVector localVertexPos(0., 0., 0.);
    //Rotate to global coordinate
    G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);

    if (1) {
      particleGun->SetParticleDefinition(VP_PiMinusP);
      G4ThreeVector localMomVP(0.,0.,1.);
      G4ThreeVector gloMomVP = geomMan.Local2GlobalDir(TgtId, localMomVP);
      particleGun->SetParticleMomentumDirection(gloMomVP);
      particleGun->SetParticleEnergy((Energy_VP - VPMass)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
      
      anaMan_->SetPrimaryVertex(localVertexPos);
      anaMan_->SetBeamMomentum(beammom);
      /*
	anaMan_->SetHypBeamMomentum(momentumProton);
	anaMan_->SetScatMesonMomentum(momentumPiMinus);
	anaMan_->SetThetaMeson(ThetaPi1);
	anaMan_->SetPhiMeson(PhiPi1);
	anaMan_->SetThetaMesonCM(ThetaPiCM);
	anaMan_->SetPhiMesonCM(PhiPiCM);
      */
      anaMan_->SetBeamMomentum(beammom);
      //anaMan_->SetDecayPiMomentum(momentumPiZero);
      //anaMan_->SetDecayNucleonMomentum(momentumProton);
      
    }
  } else {
    anaMan_->SetPrimaryFlag(false);

    G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
    G4int primaryTgtType = -1;
    // generate primary vertex position
    // This is continued untill the vertex position is LH2
    do {
      primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
      primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
      primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
      primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						   primary_vertex_y, 
						   primary_vertex_z));
    } while (primaryTgtType != 0 && primaryTgtType != 1);

    G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
    //Rotate to global coordinate
    G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);


    int nDaughter = anaMan_->GetNDaughter();
    for (int i=0; i<nDaughter; i++) {
      int pid = anaMan_->GetDaughterPID(i);
      G4ThreeVector globalMom = anaMan_->GetDaughterMomentum(i);
      G4ThreeVector localMom = geomMan.Global2LocalDir(TgtId, globalMom);

      G4ParticleDefinition* Particle;

      if (pid == 1)
	Particle = particleTable->FindParticle("pi+");
      else if (pid == 2)
	Particle = particleTable->FindParticle("pi-");
      else if (pid == 3)
	Particle = particleTable->FindParticle("proton");
      else if (pid == 4)
	Particle = particleTable->FindParticle("neutron");
      else {
	std::cerr << "CFTPrimaryGeneratorAction::GenerateVPFourBodyDecay_Scat No such particle type " 
		  << pid << std::endl;
	return;
      }
      /*
      std::cout << "pid " << pid << ", Mom = ( "
		<< localMom.x() << ", "
		<< localMom.y() << ", "
		<< localMom.z() << ") " << std::endl;
      */
      bool flagScattering = false;

      if (pid == 2)
	flagScattering = PiMinusP_Scattering(anEvent, localMom, localVertexPos);
      else if (pid == 3)
	flagScattering = ProtonP_Scattering(anEvent, localMom, localVertexPos);
      else if (pid == 4)
	flagScattering = NeutronP_Scattering(anEvent, localMom, localVertexPos);

      if (!flagScattering) {
	particleGun->SetParticleDefinition(Particle);
	G4double Mass    = Particle->GetPDGMass()/GeV;
	G4double Energy = sqrt(Mass*Mass + localMom.mag2());
	particleGun->SetParticleMomentumDirection(globalMom);
	particleGun->SetParticleEnergy((Energy - Mass)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	/*
	std::cout << "pid = " << pid
		  << "Mass = " << Mass
		  << "Ekin = " << (Energy - Mass)
		  << ", Mom = ( "
		  << localMom.x() << ", "
		  << localMom.y() << ", "
		  << localMom.z() << ") " << std::endl;
	*/
      }
    }
    double beammom = 1.05;
    G4ParticleDefinition* piMinus = particleTable->FindParticle("spi-");
    /* pi- beam */
    G4double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
    G4ThreeVector momentumBeam(0., 0., -beammom);

    /* beam */
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetBeamMomentum(beammom);
    /*
      anaMan_->SetHypBeamMomentum(momentumProton);
      anaMan_->SetScatMesonMomentum(momentumPiMinus);
      anaMan_->SetThetaMeson(ThetaPi1);
      anaMan_->SetPhiMeson(PhiPi1);
      anaMan_->SetThetaMesonCM(ThetaPiCM);
      anaMan_->SetPhiMesonCM(PhiPiCM);
    */
    anaMan_->SetBeamMomentum(beammom);
    //anaMan_->SetDecayPiMomentum(momentumPiZero);
    //anaMan_->SetDecayNucleonMomentum(momentumProton);

  }

  GenerateNum++;
}




void CFTPrimaryGeneratorAction::GenerateKMinusP_KPlusXi(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* kaonMinus = particleTable->FindParticle("skaon-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  //G4ParticleDefinition* xi = particleTable->FindParticle("xi-");

  double width = 0.0;
  double mass_xi = 1.321;
  double beammom = 1.8;

  Kinema3Resonance Xi;
  int DistFlag=3; // 1.3GeV/c (pi-, K+)

  Xi = Kinema3Resonance(kaonMinus->GetPDGMass()/GeV, 
			proton->GetPDGMass()/GeV,
			0.0 , mass_xi,
			kaonPlus->GetPDGMass()/GeV,
			mass_xi, width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Xi.GetEnergy(5);
  //double momentum_k = Xi.GetMomentum(5);
  Xi.GetMomentum(5,mom);
  G4ThreeVector momentumKaonPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Xi.GetTheta(5);
  double PhiK = Xi.GetPhi(5);
  double ThetaKCM = 180.-Xi.GetThetaCM(1);
  double PhiKCM = Xi.GetPhiCM(1);

  /* xi */
  //double Energy_xi = Xi.GetEnergy(4);
  //double momentum_xi = Xi.GetMomentum(4);
  Xi.GetMomentum(4,mom);
  G4ThreeVector momentumXi(mom[1], mom[2], mom[0]);
  //double ThetaXi = Xi.GetTheta(4);
  //double PhiXi = Xi.GetPhi(4);

  /* K- beam */
  double Energy_beam = sqrt(pow(kaonMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(kaonMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumXi);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
    /*
    PrimaryInfo pInfo;
    pInfo.p     = momentum_k; 
    pInfo.theta = ThetaK;
    pInfo.phi   = PhiK;
    pInfo.x     = primary_vertex_x;
    pInfo.y     = primary_vertex_y;
    pInfo.z     = primary_vertex_z;
    pInfo.pbeam = momentumBeam.mag();
    pInfo.ubeam = momentumBeam.x()/(-momentumBeam.z());
    pInfo.vbeam = momentumBeam.y()/(-momentumBeam.z());
    pInfo.mass = kaonPlus->GetPDGMass()/GeV;
    pInfo.charge = kaonPlus->GetPDGMass()/GeV;
    anaMan_->PrimaryGeneration( &pInfo );
    */
  }

}


void CFTPrimaryGeneratorAction::GenerateKMinusP_KPlusXi_Scat(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* kaonMinus = particleTable->FindParticle("skaon-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* xi = particleTable->FindParticle("xi-");
  G4ParticleDefinition* uxi = particleTable->FindParticle("uxi-");
  //G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");
  G4ParticleDefinition* lambda = particleTable->FindParticle("lambda");
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");

  double width = 0.0;
  double mass_xi = 1.321;
  double beammom = 1.48;

  Kinema3Resonance Xi;
  int DistFlag=3; // 1.3GeV/c (pi-, K+)

  Xi = Kinema3Resonance(kaonMinus->GetPDGMass()/GeV, 
			proton->GetPDGMass()/GeV,
			0.0 , mass_xi,
			kaonPlus->GetPDGMass()/GeV,
			mass_xi, width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Xi.GetEnergy(5);
  //double momentum_k = Xi.GetMomentum(5);
  Xi.GetMomentum(5,mom);
  G4ThreeVector momentumKaonPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Xi.GetTheta(5);
  double PhiK = Xi.GetPhi(5);
  double ThetaKCM = 180.-Xi.GetThetaCM(1);
  double PhiKCM = Xi.GetPhiCM(1);

  /* xi */
  //double Energy_xi = Xi.GetEnergy(4);
  double momentum_xi = Xi.GetMomentum(4);
  Xi.GetMomentum(4,mom);
  G4ThreeVector momentumXi(mom[1], mom[2], mom[0]);
  double ThetaXi = Xi.GetTheta(4);
  double PhiXi = Xi.GetPhi(4);

  /* K- beam */
  double Energy_beam = sqrt(pow(kaonMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  

  /* xi-p scatt */
  //G4double flength=20.*mm;
  G4double ctau=49.17; /*mm*/
  G4double p_xi = momentum_xi;
  G4double m_xi = xi->GetPDGMass()/GeV;
  G4double E_xi = 0;
  // Xi-p --> Xi-p
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_p= 0.01274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;
  // Xi-p --> LL
  G4double react_rate_ll_conv = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_lp_conv = 0.01274; // 1/mm /* temporary */
  // Xi-p --> Xi0n
  G4double react_rate_xi0n_conv = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_sp_conv = 0.01274; // 1/mm /* temporary */

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localXiPos = localVertexPos;
  G4bool flagDecay=false;
  G4bool flagScattering=false;
  G4bool flagLambdaLambdaConv=false;
  G4bool flagXi0NConv=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localXiPos);
    // Multiple scattering
    /*
    if (flagTgtType == 0 || flagTgtType == 1) {
      momentumSigma = calcMultipleScattering(p_sigma, m_sigma, dx/mm, momentumSigma, flagTgtType);

    }
    */
    totalx += dx;
    localXiPos += momentumXi*dx/momentumXi.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomXi ( " << momentumXi.x() << ", "  
	      << momentumXi.y() << ", " << momentumXi.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localXiPos.x() << ", "  
	      << localXiPos.y() << ", " << localXiPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumXi, &theta1, &phi1);
      std::cout << "Theta " << ThetaXi << "--> "  << theta1
		<< ", Phi " << PhiXi << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;
    flagScattering=false;

    flagDecay = decayCheck(ctau, p_xi, m_xi, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localXiPos);

    if (flagTgtType == 0 ) {
      react_rate = react_rate_p;
      dxInH += dx;

      flagScattering  = scatteringCheck(react_rate, dx/mm);
      flagLambdaLambdaConv  = scatteringCheck(react_rate_ll_conv, dx/mm);
      flagXi0NConv  = scatteringCheck(react_rate_xi0n_conv, dx/mm);

      p_xi = calcEnergyDeposit(p_xi, m_xi, dx/mm, localXiPos, flagTgtType);
      E_xi = sqrt(p_xi*p_xi + m_xi*m_xi);
    }
    /*
    G4cout << "Length : " << totalx << ", p = " << p_xi*1000. 
	   << " MeV/c, Ekin = "	   << (E_xi-m_xi)*1000. 
	   << " MeV" << G4endl;
    */
    if (fabs(p_xi)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }
    if (flagDecay || flagScattering || flagLambdaLambdaConv || flagXi0NConv)
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }

  momentumXi *= p_xi/momentumXi.mag();

  if (flagDecay) {
    //std::cout << "decay" << std::endl;
    G4ThreeVector localDecayPos=localXiPos;
    G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

    if (1) {
      Kinema3Resonance XiDecay;
      XiDecay = Kinema3Resonance(xi->GetPDGMass()/GeV, 
				 0.0,
				 0.0 , 
				 lambda->GetPDGMass()/GeV,
				 piMinus->GetPDGMass()/GeV,
				 lambda->GetPDGMass()/GeV,
				 0.0, p_xi, 0.0);


      /* Lambda */
      double Energy_decayL = XiDecay.GetEnergy(4);
      //double momentum_decayL = XiDecay.GetMomentum(4);
      XiDecay.GetMomentum(4,mom);
      G4ThreeVector momentumDecayL(mom[1], mom[2], mom[0]);
      //double ThetaDecayL = XiDecay.GetTheta(4);
      //double PhiDecayL = XiDecay.GetPhi(4);
      
      /* pi- */
      double Energy_decayPi = XiDecay.GetEnergy(5);
      //double momentum_decayPi = XiDecay.GetMomentum(5);
      XiDecay.GetMomentum(5,mom);
      G4ThreeVector momentumDecayPi( mom[1], mom[2], mom[0]);
      //double ThetaDecayPi = XiDecay.GetTheta(5);
      //double PhiDecayPi = XiDecay.GetPhi(5);
      
      // Lambda momentum at beam frame
      momentumDecayL.rotateY(ThetaXi*deg);
      momentumDecayL.rotateZ(PhiXi*deg);

      // pi- momentum at beam frame
      momentumDecayPi.rotateY(ThetaXi*deg);
      momentumDecayPi.rotateZ(PhiXi*deg);

      /*
      std::cout << "XiMom = ( " << momentumXi.x() << ", " 
		<< momentumXi.y() << ", " << momentumXi.z()
		<< ")" << std::endl;
      std::cout << "DecayL = ( " << momentumDecayL.x() << ", " 
		<< momentumDecayL.y() << ", " << momentumDecayL.z()
		<< ")" << std::endl;
      std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
		<< momentumDecayPi.y() << ", " << momentumDecayPi.z()
		<< ")" << std::endl;

      std::cout << "Delta = ( " 
		<< momentumXi.x()-momentumDecayL.x()-momentumDecayPi.x()
		<< ", " 
		<< momentumXi.y()-momentumDecayL.y()-momentumDecayPi.y() 
		<< ", " 
		<< momentumXi.z()-momentumDecayL.z()-momentumDecayPi.z()
		<< ")" << std::endl;
      */
      double ThetaDecayLatBeamFrame, PhiDecayLatBeamFrame; 
      calcThetaPhi(momentumDecayL, 
		   &ThetaDecayLatBeamFrame, 
		   &PhiDecayLatBeamFrame);

      double ThetaDecayPiatBeamFrame, PhiDecayPiatBeamFrame; 
      calcThetaPhi(momentumDecayPi, 
		   &ThetaDecayPiatBeamFrame, 
		   &PhiDecayPiatBeamFrame);
      /*
      std::cout << "ThetaDecayLatBeamFrame: " << ThetaDecayLatBeamFrame
		<< " PhiDecayLatBeamFrame: "  << PhiDecayLatBeamFrame
		<< std::endl;

      std::cout << "ThetaDecayPiatBeamFrame: " << ThetaDecayPiatBeamFrame
		<< " PhiDecayPiatBeamFrame: "  << PhiDecayPiatBeamFrame
		<< std::endl;
      */
      if (PhiDecayLatBeamFrame<0)
	return;
      if (PhiDecayPiatBeamFrame<0)
	return;


      /* lambda p scat */ 

      bool flagLambdaPScat = LambdaP_Scattering(anEvent, momentumDecayL, localDecayPos);
      if (flagLambdaPScat) {
	// generate K-K+ event
	particleGun->SetParticleDefinition(kaonPlus);
	G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
	particleGun->SetParticleMomentumDirection(gloMomKPlus);
	particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	/* beam */
	particleGun->SetParticleDefinition(kaonMinus);
	G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	particleGun->SetParticleMomentumDirection(gloMomBeam);
	particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	// generate pi- from X- decay
	/* decay pi- */
	particleGun->SetParticleDefinition(piMinus);
	G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	particleGun->SetParticleEnergy((Energy_decayPi - piMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);


	anaMan_->SetPrimaryVertex(localVertexPos);
	anaMan_->SetHypBeamMomentum(momentumXi);
	anaMan_->SetScatMesonMomentum(momentumKaonPlus);
	anaMan_->SetThetaMeson(ThetaK);
	anaMan_->SetPhiMeson(PhiK);
	anaMan_->SetThetaMesonCM(ThetaKCM);
	anaMan_->SetPhiMesonCM(PhiKCM);
	anaMan_->SetDecayFlag();
	anaMan_->SetBeamMomentum(beammom);
	anaMan_->SetFlightLengthInTarget(dxInH, 0);

	return;
      }

      bool flagPiPScat = PiMinusP_Scattering(anEvent, momentumDecayPi, localDecayPos);
      if (flagPiPScat) {
	// generate K-K+ event
	particleGun->SetParticleDefinition(kaonPlus);
	G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
	particleGun->SetParticleMomentumDirection(gloMomKPlus);
	particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	/* beam */
	particleGun->SetParticleDefinition(kaonMinus);
	G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	particleGun->SetParticleMomentumDirection(gloMomBeam);
	particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	// generate Lambda- from X- decay
	particleGun->SetParticleDefinition(lambda);
	G4ThreeVector gloMomDecayL = geomMan.Local2GlobalDir(TgtId, momentumDecayL);
	particleGun->SetParticleMomentumDirection(gloMomDecayL);
	particleGun->SetParticleEnergy((Energy_decayL - lambda->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	anaMan_->SetPrimaryVertex(localVertexPos);
	anaMan_->SetHypBeamMomentum(momentumXi);
	anaMan_->SetScatMesonMomentum(momentumKaonPlus);
	anaMan_->SetThetaMeson(ThetaK);
	anaMan_->SetPhiMeson(PhiK);
	anaMan_->SetThetaMesonCM(ThetaKCM);
	anaMan_->SetPhiMesonCM(PhiKCM);
	anaMan_->SetDecayFlag();
	anaMan_->SetBeamMomentum(beammom);
	anaMan_->SetFlightLengthInTarget(dxInH, 0);

	return;
      }
      

    }


    static int nDecay = 0;
    const int prescale = 100;
    //const int prescale = 1;
    nDecay++;

    if (nDecay%prescale != 0)
      return;


    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(kaonMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*decaying xi*/    
    particleGun->SetParticleDefinition(uxi);
    G4ThreeVector gloMomXi = geomMan.Local2GlobalDir(TgtId, momentumXi);
    particleGun->SetParticleMomentumDirection(gloMomXi);
    // E_sigma is energy after energy deposit
    particleGun->SetParticleEnergy((E_xi - xi->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumXi);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetDecayPos(localXiPos);
    anaMan_->SetDecayEventFlag();
    anaMan_->SetDecayFlag();
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

  } else if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localXiPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance XiPScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return;
    }

    XiPScat = Kinema3Resonance(xi->GetPDGMass()/GeV,
			       scatParticle->GetPDGMass()/GeV,
			       0.0 , 
			       xi->GetPDGMass()/GeV,
			       scatParticle->GetPDGMass()/GeV,
			       xi->GetPDGMass()/GeV, 0.0, 
			       p_xi, 0.0,
			       scatDistFlag);
    
    
    /* scattered Xi */
    G4double Energy_scatXi = XiPScat.GetEnergy(4);
    //G4double momentum_scatXi = XiPScat.GetMomentum(4);
    XiPScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatXi( mom[1], mom[2], mom[0]);
    momentumScatXi.rotateY(ThetaXi*deg);
    momentumScatXi.rotateZ(PhiXi*deg);
    
    double ThetaScatXiCM = XiPScat.GetThetaCM(1);
    double PhiScatXiCM = XiPScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatXi (" << momentumScatXi.x() << ", "
      << momentumScatXi.y() << ", " << momentumScatXi.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = XiPScat.GetEnergy(5);
    //G4double momentum_scatPart = XiPScat.GetMomentum(5);
    XiPScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaXi*deg);
    momentumScatPart.rotateZ(PhiXi*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumXi);
    anaMan_->SetScatHypMomentum(momentumScatXi);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatXiCM);
    anaMan_->SetPhiScatHypCM(PhiScatXiCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);


    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(kaonMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt xi*/    
    particleGun->SetParticleDefinition(xi);
    G4ThreeVector gloMomScatXi = geomMan.Local2GlobalDir(TgtId, momentumScatXi);
    particleGun->SetParticleMomentumDirection(gloMomScatXi);
    particleGun->SetParticleEnergy((Energy_scatXi - xi->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  }  else if (flagLambdaLambdaConv || flagXi0NConv) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localXiPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);

    Kinema3Resonance XiPScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    G4ParticleDefinition* hyperon;
    int reactMode=-100;
    if (flagLambdaLambdaConv) {
      hyperon = particleTable->FindParticle("lambda");
      scatParticle = particleTable->FindParticle("lambda");
      reactMode=2;
    } else { 
      hyperon = particleTable->FindParticle("xi0");
      scatParticle = particleTable->FindParticle("neutron");
      reactMode=3;
    }
    double mass_proton = 0.93827200;

    XiPScat = Kinema3Resonance(xi->GetPDGMass()/GeV,
				 mass_proton,
				 0.0 , 
				 hyperon->GetPDGMass()/GeV,
				 scatParticle->GetPDGMass()/GeV,
				 hyperon->GetPDGMass()/GeV, 
				 0.0, 
				 p_xi, 0.0,
				 scatDistFlag);

    
    /* scattered Hyperon */
    G4double Energy_scatXi = XiPScat.GetEnergy(4);
    //G4double momentum_scatXi = XiPScat.GetMomentum(4);
    XiPScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatXi( mom[1], mom[2], mom[0]);
    momentumScatXi.rotateY(ThetaXi*deg);
    momentumScatXi.rotateZ(PhiXi*deg);
    
    double ThetaScatXiCM = XiPScat.GetThetaCM(1);
    double PhiScatXiCM = XiPScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatXi (" << momentumScatXi.x() << ", "
      << momentumScatXi.y() << ", " << momentumScatXi.z() << ") "
      << G4endl;
    */

    /* scattered Particle*/
    G4double Energy_scatPart = XiPScat.GetEnergy(5);
    //G4double momentum_scatPart = XiPScat.GetMomentum(5);
    XiPScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaXi*deg);
    momentumScatPart.rotateZ(PhiXi*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Xi (" << momentumXi.x() << ", "
      << momentumXi.y() << ", " << momentumXi.z() << ") "
      << G4endl;
      
      G4cout << "ScatXi (" << momentumScatXi.x() << ", "
      << momentumScatXi.y() << ", " << momentumScatXi.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatPart.x() << ", "
      << momentumScatPart.y() << ", " << momentumScatPart.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumXi);
    anaMan_->SetScatHypMomentum(momentumScatXi);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKaonPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatXiCM);
    anaMan_->SetPhiScatHypCM(PhiScatXiCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag(reactMode);
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(kaonMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt hyperon*/    
    particleGun->SetParticleDefinition(hyperon);
    G4ThreeVector gloMomScatXi = geomMan.Local2GlobalDir(TgtId, momentumScatXi);
    particleGun->SetParticleMomentumDirection(gloMomScatXi);
    particleGun->SetParticleEnergy((Energy_scatXi - hyperon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
  }

}

void CFTPrimaryGeneratorAction::GenerateGammaP_KPlusLambda(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* gamma = particleTable->FindParticle("gamma");
  G4ParticleDefinition* kaon  = particleTable->FindParticle("kaon+");
  // G4ParticleDefinition* kaon  = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* lambda  = particleTable->FindParticle("lambda");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double beammom = RandFlat::shoot(1.5, 2.4);
  //double beammom = RandFlat::shoot(2.25, 2.4);
  double width = 0.0;

  Kinema3Resonance Sigma;
  int DistFlag=9; // E dep gamma p --> K+ Lambda
  // int DistFlag=0; // Random K+ direction

  Sigma = Kinema3Resonance(gamma->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , 
			   lambda->GetPDGMass()/GeV,
			   kaon->GetPDGMass()/GeV,
			   lambda->GetPDGMass()/GeV,
			   width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumK(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  // Lambda
  double E_lambda = Sigma.GetEnergy(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumLambda( mom[1], mom[2], mom[0]);

  /* gamma beam */
  double Energy_beam = sqrt(pow(gamma->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  // if (!(ThetaK>0.&&ThetaK<25.)) 
  //   return;

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    //primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_x = 0.*mm;    
    //primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = 0.*mm;        
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaon);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    G4bool flagPolarization = false; 

    if (!flagPolarization) {
      particleGun->SetParticleDefinition(lambda);
      G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
      particleGun->SetParticleMomentumDirection(gloMomLambda);
      particleGun->SetParticleEnergy((E_lambda - lambda->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalVertexPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    } else {
      G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
      
      G4ThreeVector NormZ = momentumLambda*(1/momentumLambda.mag());
      G4ThreeVector NormY = NormZ.cross(momentumK);
      NormY *= (1./NormY.mag());

      double ctau = 78.9; /*mm*/
      DecayWithPolarization LambdaDecay(lambda->GetPDGMass()/GeV,
					proton->GetPDGMass()/GeV,
					piMinus->GetPDGMass()/GeV,
					//1.0, // Polarization
					0.0, // Polarization
					AsymPara_Lambda, // Asymmetry parameter
					NormY, 
					momentumLambda,
					localVertexPos,
					false, // phi flag, lambda is opposite for x axis
					ctau,
					false // charge flag
					);

      //std::cout << "ThetaLambda = " << ThetaLambda 
      //<< "PhiLambda" << PhiLambda << std::endl;
      G4ThreeVector localDecayPos  = LambdaDecay.GetDecayPos();
      G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

      G4ThreeVector momDecayProton = LambdaDecay.GetMomentum2();

      //std::cout << "momDecayProton = " <<  momDecayProton << std::endl;

      particleGun->SetParticleDefinition(proton);
      G4ThreeVector gloMomDecayProton = geomMan.Local2GlobalDir(TgtId, momDecayProton);
      particleGun->SetParticleMomentumDirection(gloMomDecayProton);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin2()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      G4ThreeVector momDecayPi     = LambdaDecay.GetMomentum3();
      //std::cout << "momDecayPi = " <<  momDecayPi << std::endl;
      //std::cout << "momentumLambda(gene) = " << momDecayProton+momDecayPi << std::endl;
      //std::cout << "momentumLambda(org)  = " << momentumLambda << std::endl;

      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momDecayPi);
      particleGun->SetParticleMomentumDirection(gloMomDecayPi);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin3()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      anaMan_->SetDecayPiMomentum(momDecayPi);
      anaMan_->SetDecayNucleonMomentum(momDecayProton);
      anaMan_->SetDecayPos(localDecayPos);
    }

    /* beam */
    particleGun->SetParticleDefinition(gamma);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
  }

}

void CFTPrimaryGeneratorAction::GenerateGammaP_KPlusLambda_Scat(G4Event* anEvent)
{
  double mom[3];

  //G4bool flagPolarization = true;
  G4bool flagPolarization = false;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* gamma = particleTable->FindParticle("gamma");;
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* kaonPlus = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* lambda = particleTable->FindParticle("lambda");

  double beammom = RandFlat::shoot(1.5, 2.4);
  double width = 0.0;

  int DistFlag=9; // E dep gamma p --> K+ Lambda

  Kinema3Resonance Sigma =
    Kinema3Resonance(gamma->GetPDGMass()/GeV, 
		     proton->GetPDGMass()/GeV,
		     0.0 , 
		     lambda->GetPDGMass()/GeV,
		     kaonPlus->GetPDGMass()/GeV,
		     lambda->GetPDGMass()/GeV,
		     width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumK(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  // Lambda
  double E_lambda = Sigma.GetEnergy(4);
  double momentum_lambda = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumLambda( mom[1], mom[2], mom[0]);
  double ThetaLambda = Sigma.GetTheta(4);
  double PhiLambda = Sigma.GetPhi(4);

  /* gamma beam */
  double Energy_beam = sqrt(pow(gamma->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);

  if (!(ThetaK>0.&&ThetaK<25.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = 0.*mm;
    primary_vertex_y = 0.*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);

  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  

  /* Lambda p scatt */
  //G4double flength=20.*mm;
  G4double ctau=78.9; /*mm*/
  G4double p_lambda = momentum_lambda;
  G4double m_lambda = lambda->GetPDGMass()/GeV;
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_p= 0.01274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;

  G4double react_rate_sigmaPn = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  G4double react_rate_sigma0p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_sigmaPn = 0.01274; // 1/mm /* temporary */
  //G4double react_rate_sigma0p = 0.01274; // 1/mm /* temporary */

  G4double pth_sigmaPn = 0.6345; 
  G4double pth_sigma0p = 0.6435; 

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localLambdaPos = localVertexPos;
  G4bool flagDecay=false;
  G4bool flagScattering=false;
  G4bool flagSigmaPlusN=false;
  G4bool flagSigma0P=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localLambdaPos);
    totalx += dx;
    localLambdaPos += momentumLambda*dx/momentumLambda.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomLambda ( " << momentumLambda.x() << ", "  
	      << momentumLambda.y() << ", " << momentumLambda.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localLambdaPos.x() << ", "  
	      << localLambdaPos.y() << ", " << localLambdaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumLambda, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;
    flagScattering=false;
    flagSigmaPlusN=false;
    flagSigma0P=false;

    flagDecay = decayCheck(ctau, p_lambda, m_lambda, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localLambdaPos);

    if (flagTgtType == 0 ) {
      react_rate = react_rate_p;
      dxInH += dx;

      flagScattering  = scatteringCheck(react_rate, dx/mm);
      if (p_lambda > pth_sigmaPn)
	flagSigmaPlusN  = scatteringCheck(react_rate_sigmaPn, dx/mm);
      if (p_lambda > pth_sigma0p)
	flagSigma0P  = scatteringCheck(react_rate_sigma0p, dx/mm);
    }
    /*
    G4cout << "Length : " << totalx << ", p = " << p_xi*1000. 
	   << " MeV/c, Ekin = "	   << (E_xi-m_xi)*1000. 
	   << " MeV" << G4endl;
    */
    if (fabs(p_lambda)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }
    if (flagDecay || flagScattering || flagSigmaPlusN || flagSigma0P)
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }

  momentumLambda *= p_lambda/momentumLambda.mag();

  if (flagDecay) {
    //std::cout << "decay" << std::endl;
    G4ThreeVector localDecayPos=localLambdaPos;
    G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

    if (1) {
      static int Ndecay=0;
      int decayMode=-1;
      double randval = RandFlat::shoot(0., 1.);
      if (randval>0.639)
	decayMode = 0; // pi0n
      else
	decayMode = 1; // pi-p

      Ndecay++;

      G4ParticleDefinition* decayNucl;
      G4ParticleDefinition* decayPi;
      if (decayMode==1) {
	decayNucl = particleTable->FindParticle("proton_d");
	decayPi = particleTable->FindParticle("pi-");
      } else if (decayMode==0) {
	decayNucl = particleTable->FindParticle("neutron");
	decayPi = particleTable->FindParticle("pi0");
      }

      Kinema3Resonance LambdaDecay =
	Kinema3Resonance(lambda->GetPDGMass()/GeV,
			 0.0,
			 0.0 , decayNucl->GetPDGMass()/GeV,
			 decayPi->GetPDGMass()/GeV,
			 decayNucl->GetPDGMass()/GeV,
			 0.0, p_lambda, 0.0);

      /* Neucleon */
      double Energy_decayN = LambdaDecay.GetEnergy(4);
      //double momentum_decayN = LambdaDecay.GetMomentum(4);
      LambdaDecay.GetMomentum(4,mom);
      G4ThreeVector momentumDecayN(mom[1], mom[2], mom[0]);
      //double ThetaDecayN = LambdaDecay.GetTheta(4);
      //double PhiDecayN = LambdaDecay.GetPhi(4);
    
      /* pi */
      double Energy_decayPi = LambdaDecay.GetEnergy(5);
      //double momentum_decayPi = LambdaDecay.GetMomentum(5);
      LambdaDecay.GetMomentum(5,mom);
      G4ThreeVector momentumDecayPi( mom[1], mom[2], mom[0]);
      //double ThetaDecayPi = LambdaDecay.GetTheta(5);
      //double PhiDecayPi = LambdaDecay.GetPhi(5);

      //double ThetaLambda=0, PhiLambda=0; 
      //calcThetaPhi(momentumLambda, 
      //&ThetaLambda, 
      //&PhiLambda);
      
      // Neucleon momentum at beam frame
      momentumDecayN.rotateY(ThetaLambda*deg);
      momentumDecayN.rotateZ(PhiLambda*deg);
    
      // pi momentum at beam frame
      momentumDecayPi.rotateY(ThetaLambda*deg);
      momentumDecayPi.rotateZ(PhiLambda*deg);

      bool flagNpScattering = false;
      bool flagPpScattering = false;
      bool flagPipScattering = false;

      if (decayMode==0) {
	// pi0 n
	flagNpScattering = NeutronP_Scattering(anEvent, momentumDecayN, localDecayPos);
      
	if (flagNpScattering) {
	  // generate decay pi0
	
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	  particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	  particleGun->SetParticleEnergy((Energy_decayPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	
	  //return;
	}
      } else if (decayMode==1) {
      
	flagPpScattering = ProtonP_Scattering(anEvent, momentumDecayN, localDecayPos);
	flagPipScattering = PiMinusP_Scattering(anEvent, momentumDecayPi, localDecayPos);
	if (flagPpScattering) {
	  // generate decay pi-
	
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	  particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	  particleGun->SetParticleEnergy((Energy_decayPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  
	  //return;
	} else if (flagPipScattering) {
	  // generate decay n
	
	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomDecayN = geomMan.Local2GlobalDir(TgtId, momentumDecayN);
	  particleGun->SetParticleMomentumDirection(gloMomDecayN);
	  particleGun->SetParticleEnergy((Energy_decayN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  
	  //return;
	
	}
      }

      if (flagNpScattering || flagPpScattering || flagPipScattering) {
	particleGun->SetParticleDefinition(kaonPlus);
	G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumK);
	particleGun->SetParticleMomentumDirection(gloMomKPlus);
	particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	/* beam */
	particleGun->SetParticleDefinition(gamma);
	G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	particleGun->SetParticleMomentumDirection(gloMomBeam);
	particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	anaMan_->SetPrimaryVertex(localVertexPos);
	anaMan_->SetHypBeamMomentum(momentumLambda);
	anaMan_->SetScatMesonMomentum(momentumK);
	anaMan_->SetThetaMeson(ThetaK);
	anaMan_->SetPhiMeson(PhiK);
	anaMan_->SetThetaMesonCM(ThetaKCM);
	anaMan_->SetPhiMesonCM(PhiKCM);
	anaMan_->SetDecayPos(localLambdaPos);
	anaMan_->SetDecayEventFlag();
	anaMan_->SetDecayFlag();
	anaMan_->SetFlightLengthInTarget(dxInH, 0);
	anaMan_->SetBeamMomentum(beammom);
	
	return;
      }

      
    }

    static int nDecay = 0;
    const int prescale = 100;
    //const int prescale = 1;
    nDecay++;

    if (nDecay%prescale != 0)
      return;


    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(gamma);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetDecayPos(localLambdaPos);
    anaMan_->SetDecayEventFlag();
    anaMan_->SetDecayFlag();
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);
    
    /*decaying lambda*/    
    if (!flagPolarization) {
      G4ParticleDefinition*  ulambda  = particleTable->FindParticle("ulambda");
      particleGun->SetParticleDefinition(ulambda);
      
      G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
      particleGun->SetParticleMomentumDirection(gloMomLambda);
      particleGun->SetParticleEnergy((E_lambda - lambda->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    }  else {
      G4ThreeVector NormZ = momentumLambda*(1/momentumLambda.mag());
      G4ThreeVector NormY = NormZ.cross(momentumK);
      NormY *= (1./NormY.mag());

      G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
      double ctau_ulambda = 0.; /*mm*/
      DecayWithPolarization LambdaDecay(lambda->GetPDGMass()/GeV,
					proton->GetPDGMass()/GeV,
					piMinus->GetPDGMass()/GeV,
					1.0, // Polarization
					AsymPara_Lambda, // Asymmetry parameter
					NormY, 
					momentumLambda,
					localDecayPos,
					false, // phi flag, lambda is opposite for x axis
					ctau_ulambda,
					false // charge flag
					);

      //std::cout << "ThetaLambda = " << ThetaLambda 
      //<< "PhiLambda" << PhiLambda << std::endl;
      G4ThreeVector momDecayProton = LambdaDecay.GetMomentum2();

      //std::cout << "momDecayProton = " <<  momDecayProton << std::endl;

      particleGun->SetParticleDefinition(proton);
      G4ThreeVector gloMomDecayProton = geomMan.Local2GlobalDir(TgtId, momDecayProton);
      particleGun->SetParticleMomentumDirection(gloMomDecayProton);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin2()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      G4ThreeVector momDecayPi     = LambdaDecay.GetMomentum3();
      //std::cout << "momDecayPi = " <<  momDecayPi << std::endl;
      //std::cout << "momentumLambda(gene) = " << momDecayProton+momDecayPi << std::endl;
      //std::cout << "momentumLambda(org)  = " << momentumLambda << std::endl;

      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momDecayPi);
      particleGun->SetParticleMomentumDirection(gloMomDecayPi);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin3()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      anaMan_->SetDecayPiMomentum(momDecayPi);
      anaMan_->SetDecayNucleonMomentum(momDecayProton);
      anaMan_->SetDecayPos(localDecayPos);
    }

  } else if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localLambdaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance LambdaScat;
    G4ParticleDefinition* scatParticle;
    //int scatDistFlag = 10; // 10 " : 90 < phi < 270, theta : flat
    int scatDistFlag = 0; // 10 " : 0 < phi < 360, theta : flat
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
    } else {
      return;
    }

    LambdaScat = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , 
				  lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  lambda->GetPDGMass()/GeV,
				  0.0, p_lambda, 0.0,
				  scatDistFlag);
    
    
    /* scattered lambda */
    G4double Energy_scatLambda = LambdaScat.GetEnergy(4);
    //G4double momentum_scatLambda = LambdaScat.GetMomentum(4);
    LambdaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatLambda( mom[1], mom[2], mom[0]);
    momentumScatLambda.rotateY(ThetaLambda*deg);
    momentumScatLambda.rotateZ(PhiLambda*deg);
    
    G4double ThetaScatLambdaCM = LambdaScat.GetThetaCM(1);
    G4double PhiScatLambdaCM   = LambdaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = LambdaScat.GetEnergy(5);
    //G4double momentum_scatPart = LambdaScat.GetMomentum(5);
    LambdaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaLambda*deg);
    momentumScatPart.rotateZ(PhiLambda*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Lambda (" << momentumLambda.x() << ", "
      << momentumLambda.y() << ", " << momentumLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatHypMomentum(momentumScatLambda);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatLambdaCM);
    anaMan_->SetPhiScatHypCM(PhiScatLambdaCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(gamma);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt lambda*/    
    if (!flagPolarization) {
      G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
      particleGun->SetParticleDefinition(lambda);
      G4ThreeVector gloMomScatLambda = geomMan.Local2GlobalDir(TgtId, momentumScatLambda);
      particleGun->SetParticleMomentumDirection(gloMomScatLambda);
      particleGun->SetParticleEnergy((Energy_scatLambda - lambda->GetPDGMass()/GeV)*GeV);
      particleGun->SetParticlePosition(globalScatPos);
      particleGun->GeneratePrimaryVertex(anEvent);
    } else {
      G4ThreeVector NormZ = momentumLambda*(1/momentumLambda.mag());
      G4ThreeVector NormY = NormZ.cross(momentumScatLambda);
      NormY *= (1./NormY.mag());

      G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
      double ctau_lambda = 78.9; /*mm*/
      DecayWithPolarization LambdaDecay(lambda->GetPDGMass()/GeV,
					proton->GetPDGMass()/GeV,
					piMinus->GetPDGMass()/GeV,
					1.0, // Polarization
					AsymPara_Lambda, // Asymmetry parameter
					NormY,
					momentumScatLambda,
					localScatPos,
					true, // phi flag, lambda is same with x axis
					ctau_lambda,
					false // charge flag
					);

      //std::cout << "ThetaLambda = " << ThetaLambda 
      //<< "PhiLambda" << PhiLambda << std::endl;
      G4ThreeVector localDecayPos  = LambdaDecay.GetDecayPos();
      G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);

      G4ThreeVector momDecayProton = LambdaDecay.GetMomentum2();

      //std::cout << "momDecayProton = " <<  momDecayProton << std::endl;

      particleGun->SetParticleDefinition(proton);
      G4ThreeVector gloMomDecayProton = geomMan.Local2GlobalDir(TgtId, momDecayProton);
      particleGun->SetParticleMomentumDirection(gloMomDecayProton);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin2()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      G4ThreeVector momDecayPi     = LambdaDecay.GetMomentum3();
      //std::cout << "momDecayPi = " <<  momDecayPi << std::endl;
      //std::cout << "momentumLambda(gene) = " << momDecayProton+momDecayPi << std::endl;
      //std::cout << "momentumLambda(org)  = " << momentumLambda << std::endl;

      particleGun->SetParticleDefinition(piMinus);
      G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momDecayPi);
      particleGun->SetParticleMomentumDirection(gloMomDecayPi);
      particleGun->SetParticleEnergy(LambdaDecay.GetEkin3()*GeV);
      particleGun->SetParticlePosition(globalDecayPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      anaMan_->SetDecayPiMomentum(momDecayPi);
      anaMan_->SetDecayNucleonMomentum(momDecayProton);
      anaMan_->SetDecayPos(localDecayPos);
    }
    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    
  }  else if (flagSigmaPlusN || flagSigma0P) {
    //std::cout << "conversion scat" << std::endl;
    G4ThreeVector localScatPos=localLambdaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance LambdaScat;
    G4ParticleDefinition* scatParticle;
    //int scatDistFlag = 10; // 10 " : 90 < phi < 270, theta : flat
    int scatDistFlag = 0; // 10 " : 0 < phi < 360, theta : flat
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
    } else {
      return;
    }
    G4ParticleDefinition* scatHyperon = particleTable->FindParticle("sigma+");
    G4ParticleDefinition* recoilNucleon = particleTable->FindParticle("neutron");      
    int reactMode=2;

    if (flagSigma0P) {
      scatHyperon = particleTable->FindParticle("sigma0");
      recoilNucleon = particleTable->FindParticle("proton");      
      reactMode=3;
    }


    LambdaScat = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , 
				  scatHyperon->GetPDGMass()/GeV,
				  recoilNucleon->GetPDGMass()/GeV,
				  scatHyperon->GetPDGMass()/GeV,
				  0.0, p_lambda, 0.0,
				  scatDistFlag);
    
    
    /* scattered Hyperon */
    G4double Energy_scatHyperon = LambdaScat.GetEnergy(4);
    //G4double momentum_scatHyperon = LambdaScat.GetMomentum(4);
    LambdaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatHyperon( mom[1], mom[2], mom[0]);
    momentumScatHyperon.rotateY(ThetaLambda*deg);
    momentumScatHyperon.rotateZ(PhiLambda*deg);
    
    G4double ThetaScatHyperonCM = LambdaScat.GetThetaCM(1);
    G4double PhiScatHyperonCM   = LambdaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = LambdaScat.GetEnergy(5);
    //G4double momentum_scatPart = LambdaScat.GetMomentum(5);
    LambdaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaLambda*deg);
    momentumScatPart.rotateZ(PhiLambda*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Lambda (" << momentumLambda.x() << ", "
      << momentumLambda.y() << ", " << momentumLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    //std::cout << "reactMomde = "  << reactMode << std::endl;

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatHypMomentum(momentumScatHyperon);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatHyperonCM);
    anaMan_->SetPhiScatHypCM(PhiScatHyperonCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag(reactMode);
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kaonPlus);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(gamma);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt lambda*/    
    G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
    particleGun->SetParticleDefinition(scatHyperon);
    G4ThreeVector gloMomScatHyperon = geomMan.Local2GlobalDir(TgtId, momentumScatHyperon);
    particleGun->SetParticleMomentumDirection(gloMomScatHyperon);
    particleGun->SetParticleEnergy((Energy_scatHyperon - scatHyperon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(recoilNucleon);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - recoilNucleon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);


  }

}

void CFTPrimaryGeneratorAction::GenerateGammaP_KPlusSigma(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* gamma = particleTable->FindParticle("gamma");
  G4ParticleDefinition* kaon  = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* sigma  = particleTable->FindParticle("sigma0");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double beammom = RandFlat::shoot(1.5, 2.4);
  double width = 0.0;

  Kinema3Resonance Sigma;
  int DistFlag=9; // E dep gamma p --> K+ Lambda

  Sigma = Kinema3Resonance(gamma->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , 
			   sigma->GetPDGMass()/GeV,
			   kaon->GetPDGMass()/GeV,
			   sigma->GetPDGMass()/GeV,
			   width, beammom, 0.0, DistFlag);

  /* K+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumK(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  // Sigma
  double E_sigma = Sigma.GetEnergy(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma( mom[1], mom[2], mom[0]);

  /* gamma beam */
  double Energy_beam = sqrt(pow(gamma->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  if (!(ThetaK>0.&&ThetaK<25.)) 
    return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    //primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_x = 0.*mm;    
    //primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = 0.*mm;        
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaon);
    G4ThreeVector gloMomK = geomMan.Local2GlobalDir(TgtId, momentumK);
    particleGun->SetParticleMomentumDirection(gloMomK);
    particleGun->SetParticleEnergy((Energy_k - kaon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);


    particleGun->SetParticleDefinition(sigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(gloMomSigma);
    particleGun->SetParticleEnergy((E_sigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(gamma);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumK);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetBeamMomentum(beammom);
  }

}



void CFTPrimaryGeneratorAction::GenerateGammaP_PhiP(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* gamma = particleTable->FindParticle("gamma");
  G4ParticleDefinition* kaonP = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* kaonM = particleTable->FindParticle("skaon-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");

  double beammom = RandFlat::shoot(1.6, 2.4);
  double mass = 1019.456/1000.;
  double width = 4.26/1000.;

  Kinema3Resonance Sigma;
  int DistFlag=0; // 

  Sigma = Kinema3Resonance(gamma->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   kaonP->GetPDGMass()/GeV, 
			   kaonM->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   mass,
			   width, beammom, 0.0, DistFlag);

  /* proton */
  double Energy_p = Sigma.GetEnergy(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumP(mom[1], mom[2], mom[0]);
  double ThetaP = Sigma.GetTheta(5);
  double PhiP = Sigma.GetPhi(5);

  /* K+ */
  double Energy_KP = Sigma.GetEnergy(3);
  Sigma.GetMomentum(3,mom);
  G4ThreeVector momentumKP(mom[1], mom[2], mom[0]);
  double ThetaKP = Sigma.GetTheta(3);
  double PhiKP = Sigma.GetPhi(3);

  /* K- */
  double Energy_KM = Sigma.GetEnergy(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumKM(mom[1], mom[2], mom[0]);
  double ThetaKM = Sigma.GetTheta(4);
  double PhiKM = Sigma.GetPhi(4);

  /* gamma beam */
  double Energy_beam = sqrt(pow(gamma->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam(0., 0., -beammom);


  //if (!(ThetaK>0.&&ThetaK<25.)) 
  //return;


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    //primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_x = 0.*mm;    
    //primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = 0.*mm;        
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);
  /*
  primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 100.)*mm;
  */
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  if (1) {
    particleGun->SetParticleDefinition(kaonP);
    G4ThreeVector gloMomKP = geomMan.Local2GlobalDir(TgtId, momentumKP);
    particleGun->SetParticleMomentumDirection(gloMomKP);
    particleGun->SetParticleEnergy((Energy_KP - kaonP->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(kaonM);
    G4ThreeVector gloMomKM = geomMan.Local2GlobalDir(TgtId, momentumKM);
    particleGun->SetParticleMomentumDirection(gloMomKM);
    particleGun->SetParticleEnergy((Energy_KM - kaonM->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleDefinition(proton);
    G4ThreeVector gloMomP = geomMan.Local2GlobalDir(TgtId, momentumP);
    particleGun->SetParticleMomentumDirection(gloMomP);
    particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(gamma);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - gamma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetBeamMomentum(beammom);
  }

}

void CFTPrimaryGeneratorAction::GenerateUniform_PiPlus(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* piPlus;

  piPlus = particleTable->FindParticle("pi+");

  //double mom = 0.8;
  double mom = (G4double)RandFlat::shoot(0.9, 2.5);
  double Energy_pi = sqrt(mom*mom+
			 (piPlus->GetPDGMass()/GeV)*(piPlus->GetPDGMass()/GeV));

  double cost = (G4double)RandFlat::shoot(sqrt(3.)/2., 1.);
  double theta = acos(cost)*Rad2Deg;
  //double theta = 0.;
  double phi =  (G4double)RandFlat::shoot(0., 360.);
  G4ThreeVector momentumPiPlus(sin(theta*Deg2Rad)*cos(phi*Deg2Rad),
				 sin(theta*Deg2Rad)*sin(phi*Deg2Rad),
				 cos(theta*Deg2Rad));


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = 0.*mm;
    primary_vertex_y = 0.*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);


  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  particleGun->SetParticleDefinition(piPlus);
  G4ThreeVector gloMomPiPlus = geomMan.Local2GlobalDir(TgtId, momentumPiPlus);
  particleGun->SetParticleMomentumDirection(gloMomPiPlus);
  particleGun->SetParticleEnergy((Energy_pi - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalVertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetScatMesonMomentum(momentumPiPlus);
  anaMan_->SetThetaMeson(theta);
  anaMan_->SetPhiMeson(phi);


}

void CFTPrimaryGeneratorAction::GenerateUniform_Proton(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton;

  proton = particleTable->FindParticle("proton");

  //double mom = 0.8;
  double mom = (G4double)RandFlat::shoot(0.9, 2.5);
  double Energy_p = sqrt(mom*mom+
			 (proton->GetPDGMass()/GeV)*(proton->GetPDGMass()/GeV));

  double cost = (G4double)RandFlat::shoot(sqrt(3.)/2., 1.);
  double theta = acos(cost)*Rad2Deg;
  //double theta = 0.;
  double phi =  (G4double)RandFlat::shoot(0., 360.);
  G4ThreeVector momentumProton(sin(theta*Deg2Rad)*cos(phi*Deg2Rad),
				 sin(theta*Deg2Rad)*sin(phi*Deg2Rad),
				 cos(theta*Deg2Rad));


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = 0.*mm;
    primary_vertex_y = 0.*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);


  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  particleGun->SetParticleDefinition(proton);
  G4ThreeVector gloMomProton = geomMan.Local2GlobalDir(TgtId, momentumProton);
  particleGun->SetParticleMomentumDirection(gloMomProton);
  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalVertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetScatMesonMomentum(momentumProton);
  anaMan_->SetThetaMeson(theta);
  anaMan_->SetPhiMeson(phi);


}


bool CFTPrimaryGeneratorAction::NeutronP_Scattering(G4Event* anEvent, G4ThreeVector momentumDecayN, G4ThreeVector Pos0)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");

  G4ParticleDefinition* neutron;
  neutron = particleTable->FindParticle("neutron");
  G4double m_neutron = neutron->GetPDGMass()/GeV;


  double ThetaDecayNatBeamFrame, PhiDecayNatBeamFrame; 
  calcThetaPhi(momentumDecayN, 
	       &ThetaDecayNatBeamFrame, 
	       &PhiDecayNatBeamFrame);

  if (PhiDecayNatBeamFrame<0)
    return false;

  /* n-p scatt */
  G4double neutron_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/

  double Ekin_n[30] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
		       110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
		       210., 220., 230., 240., 250., 260., 270., 280., 290., 300.};

  double np_cs_table[30] = { 951.827, 488.477, 312.839, 223.243, 170.507,
			     136.696, 113.728,  97.449,  85.517,  76.527,
			     69.593,  64.133,  59.756,  56.191,  53.244,
			     50.779,  48.691,  46.905,  45.364,  44.022,
			     42.844,  41.803,  40.876,  40.046,  39.298,
			     38.620,  38.002,  37.436,  36.915,  36.432};
  
  double momentum_decayN = momentumDecayN.mag();

  double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_neutron*m_neutron)
		     -m_neutron)*1000. ;// MeV


  int index_Ekin=-1;
  double np_cs;
  for (int i=0; i<30; i++) {
    if (i==0) {
      if (Ekin < Ekin_n[i]) {
	index_Ekin = i;
	np_cs = np_cs_table[i];
	break;
      } else if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	index_Ekin = i;
	np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	break;
      }
    } else if (i>=1 && i<=28) {
      if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	index_Ekin = i;
	np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	break;
      }
    } else if (i==29) {
      if (Ekin >= Ekin_n[i]) {
	index_Ekin = i;
	np_cs = np_cs_table[i];
	break;
      }
    }
  }
  
  if (index_Ekin<0 || index_Ekin>=30) {
    fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
    exit(-1);
  }
  neutron_react_rate *= np_cs/30.;
  
  G4double dy = 0.5*mm; //mm
  G4double totaly=0.0;  //mm
  
  G4ThreeVector localNeutronPos = Pos0;
  int flagTgtType2 = -1;
  bool flagNeutronScattering=false;
  int nIte=0;
  while (1) {
    flagTgtType2 = -1;
    totaly += dy;
    localNeutronPos += momentumDecayN*dy/momentumDecayN.mag();
    
    flagNeutronScattering=false;
    
    flagTgtType2 = getTargetFlag(localNeutronPos);
    
    if (flagTgtType2 == 0 ) {
      flagNeutronScattering  = scatteringCheck(neutron_react_rate, dy/mm);
    }
    
    if (flagNeutronScattering)
      break;
    
    nIte++;
    if (nIte>1000)
      break;
  }

  if (flagNeutronScattering) {
    //std::cout << "NNscat" << std::endl;
    G4ThreeVector localNNScatPos=localNeutronPos;
    G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
    int scatDistFlag = 0;
    G4ParticleDefinition* scatParticle;
    if ( flagTgtType2 == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 100+index_Ekin;
    } else {
      return false;
    }
    //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
    Kinema3Resonance NNScat;
    NNScat = Kinema3Resonance(neutron->GetPDGMass()/GeV, 
			      scatParticle->GetPDGMass()/GeV,
			      0.0 , neutron->GetPDGMass()/GeV,
			      scatParticle->GetPDGMass()/GeV,
			      neutron->GetPDGMass()/GeV,
			      0.0, momentum_decayN, 0.0, scatDistFlag);

    double ThetaScatNeutronCM, PhiScatNeutronCM;
    ThetaScatNeutronCM = NNScat.GetThetaCM(1);
    PhiScatNeutronCM = NNScat.GetPhiCM(1);

    double mom[3];
    
    /* scattered Neutron */
    G4double Energy_scatN = NNScat.GetEnergy(4);
    //G4double momentum_scatN = NNScat.GetMomentum(4);
    NNScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);
    
    momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
    momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
    
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = NNScat.GetEnergy(5);
    //G4double momentum_scatPart = NNScat.GetMomentum(5);
    NNScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
    momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);
    
    /*
      std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
      << momentumDecayN.y() << ", " << momentumDecayN.z()
      << ")" << std::endl;
      
      std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
      << momentumScatN.y() << ", " << momentumScatN.z()
      << ")" << std::endl;
      
      std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
      << momentumScatPart.y() << ", " << momentumScatPart.z()
      << ")" << std::endl;
      
      std::cout << "Delta = ( " 
      << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
      << ", " 
      << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
      << ", " 
      << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
      << ")" << std::endl;
    */


    /* decay pi- */
    /*
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
    particleGun->SetParticleMomentumDirection(gloMomDecayPi);
    particleGun->SetParticleEnergy((Energy_decayPi - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */
    
    /*scat neutron*/    
    particleGun->SetParticleDefinition(neutron);
    G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
    particleGun->SetParticleMomentumDirection(gloMomentumScatN);
    particleGun->SetParticleEnergy((Energy_scatN - neutron->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalNNScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalNNScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    //anaMan_->SetPrimaryVertex(localVertexPos);
    //anaMan_->SetHypBeamMomentum(momentumSigma);
    //anaMan_->SetScatMesonMomentum(momentumKPlus);
    //anaMan_->SetThetaMeson(ThetaK);
    //anaMan_->SetPhiMeson(PhiK);
    //anaMan_->SetThetaMesonCM(ThetaKCM);
    //anaMan_->SetPhiMesonCM(PhiKCM);
    //anaMan_->SetDecayPos(localDecayPos);
    anaMan_->SetNNScatPos(localNNScatPos);
    anaMan_->SetDecayFlag();
    anaMan_->SetNNScatFlag();
    anaMan_->SetNNScatTarget(flagTgtType2);	
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    //anaMan_->SetFlightLengthInTarget(dxInH, 0);
    //anaMan_->SetDecayPiMomentum(momentumDecayPi);
    //anaMan_->SetDecayNucleonMomentum(momentumDecayN);
    //anaMan_->SetBeamMomentum(beammom);
    anaMan_->SetThetaScatHypCM(ThetaScatNeutronCM);
    anaMan_->SetPhiScatHypCM(PhiScatNeutronCM);

    return true;
  }
  
  return false;

}


bool CFTPrimaryGeneratorAction::ProtonP_Scattering(G4Event* anEvent, G4ThreeVector momentumDecayN, G4ThreeVector Pos0)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");

  G4ParticleDefinition* proton;
  G4ParticleDefinition* decayNucl;

  proton = particleTable->FindParticle("proton");
  decayNucl = particleTable->FindParticle("proton_d");

  G4double m_proton = proton->GetPDGMass()/GeV;


  double ThetaDecayNatBeamFrame, PhiDecayNatBeamFrame; 
  calcThetaPhi(momentumDecayN, 
	       &ThetaDecayNatBeamFrame, 
	       &PhiDecayNatBeamFrame);

  if (PhiDecayNatBeamFrame<0)
    return false;

  bool flagProtonScattering=false;

  double momentum_decayN = momentumDecayN.mag();

  G4double proton_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
	
  double Ekin_p[14] = {18.2, 19.8, 25.63, 30.14, 39.4, 68.3, 95., 98., 118., 142.,
		       147., 172., 250., 312.};
  
  double pp_cs_table[14] = { 351.8, 314.1, 238.7, 188.4, 138.2, 81.6,
			     56.5, 56.5, 52.7, 52.7, 51.5, 50.2, 50.2, 46.4};
  
  double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_proton*m_proton)
		 -m_proton)*1000. ;// MeV
	
  int index_Ekin=-1;
  double pp_cs;
  for (int i=0; i<14; i++) {
    if (i==0) {
      if (Ekin < Ekin_p[i]) {
	index_Ekin = i;
	pp_cs = pp_cs_table[i];
	break;
      } else if (Ekin >= Ekin_p[i] && Ekin < Ekin_p[i+1]) {
	index_Ekin = i;
	pp_cs = (pp_cs_table[i+1]*(Ekin-Ekin_p[i]) + pp_cs_table[i]*(Ekin_p[i+1]-Ekin))/(Ekin_p[i+1]-Ekin_p[i]);
	break;
      }
    } else if (i>=1 && i<=12) {
      if (Ekin >= Ekin_p[i] && Ekin < Ekin_p[i+1]) {
	index_Ekin = i;
	pp_cs = (pp_cs_table[i+1]*(Ekin-Ekin_p[i]) + pp_cs_table[i]*(Ekin_p[i+1]-Ekin))/(Ekin_p[i+1]-Ekin_p[i]);
	break;
      }
    } else if (i==13) {
      if (Ekin >= Ekin_p[i]) {
	index_Ekin = i;
	pp_cs = pp_cs_table[i];
	break;
      }
    }
  }
  
  if (index_Ekin<0 || index_Ekin>=14) {
    fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
    exit(-1);
  }
  proton_react_rate *= pp_cs/30.;
  
  G4double dy = 0.5*mm; //mm
  G4double totaly=0.0;  //mm
  
  G4ThreeVector localProtonPos = Pos0;
  int flagTgtType2 = -1;
	
  int nIte=0;
  while (1) {
    flagTgtType2 = -1;
    totaly += dy;
    localProtonPos += momentumDecayN*dy/momentumDecayN.mag();
    
    flagProtonScattering=false;
    
    flagTgtType2 = getTargetFlag(localProtonPos);
    
    if (flagTgtType2 == 0 || flagTgtType2 == 1) {
      flagProtonScattering  = scatteringCheck(proton_react_rate, dy/mm);
    }
    
    if (flagProtonScattering)
      break;
    
    nIte++;
    if (nIte>1000)
      break;
  }
  
  if (flagProtonScattering) {
    //std::cout << "*********************** pp scat ***********************" << std::endl;
    G4ThreeVector localNNScatPos=localProtonPos;
    G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
    int scatDistFlag = 0;
    G4ParticleDefinition* scatParticle;
    if ( flagTgtType2 == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 1;  // flat distribution
    } else {
      return false;
    }
    //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
    Kinema3Resonance NNScat;
    NNScat = Kinema3Resonance(decayNucl->GetPDGMass()/GeV, 
			      scatParticle->GetPDGMass()/GeV,
			      0.0 , decayNucl->GetPDGMass()/GeV,
			      scatParticle->GetPDGMass()/GeV,
			      decayNucl->GetPDGMass()/GeV,
			      0.0, momentum_decayN, 0.0, scatDistFlag);
	  
    double mom[3];

    /* scattered Necleon */
    G4double Energy_scatN = NNScat.GetEnergy(4);
    //G4double momentum_scatN = NNScat.GetMomentum(4);
    NNScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);
	  
    momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
    momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
	  
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = NNScat.GetEnergy(5);
    //G4double momentum_scatPart = NNScat.GetMomentum(5);
    NNScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
    momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);
    
    /*
      std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
      << momentumDecayN.y() << ", " << momentumDecayN.z()
      << ")" << std::endl;
      
      std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
      << momentumScatN.y() << ", " << momentumScatN.z()
      << ")" << std::endl;
      
      std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
      << momentumScatPart.y() << ", " << momentumScatPart.z()
      << ")" << std::endl;
      
      std::cout << "Delta = ( " 
      << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
      << ", " 
      << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
      << ", " 
      << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
      << ")" << std::endl;
    */
	  
    /*scat proton*/    
    //G4cout << "pp scat" << G4endl;

    particleGun->SetParticleDefinition(decayNucl);
    G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
    particleGun->SetParticleMomentumDirection(gloMomentumScatN);
    particleGun->SetParticleEnergy((Energy_scatN - decayNucl->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalNNScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalNNScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
	  
    //anaMan_->SetPrimaryVertex(localVertexPos);
    //anaMan_->SetHypBeamMomentum(momentumSigma);
    //anaMan_->SetScatMesonMomentum(momentumKPlus);
    //anaMan_->SetThetaMeson(ThetaK);
    //anaMan_->SetPhiMeson(PhiK);
    //anaMan_->SetThetaMesonCM(ThetaKCM);
    //anaMan_->SetPhiMesonCM(PhiKCM);
    //anaMan_->SetDecayPos(localDecayPos);
    anaMan_->SetNNScatPos(localNNScatPos);
    //anaMan_->SetDecayFlag(decayMode);
    anaMan_->SetNNScatFlag();
    anaMan_->SetNNScatTarget(flagTgtType2);	
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    //anaMan_->SetFlightLengthInTarget(dxInH, 0);
    //anaMan_->SetDecayPiMomentum(momentumDecayPi);
    anaMan_->SetDecayNucleonMomentum(momentumDecayN);
    //anaMan_->SetBeamMomentum(beammom);
    
    return true;
  }

  return false;
}


bool CFTPrimaryGeneratorAction::PiMinusP_Scattering(G4Event* anEvent, G4ThreeVector momentumDecayPi, G4ThreeVector Pos0)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");


  G4ParticleDefinition* piMinus;
  piMinus = particleTable->FindParticle("pi-");

  double ThetaDecayPiatBeamFrame, PhiDecayPiatBeamFrame; 
  calcThetaPhi(momentumDecayPi, 
	       &ThetaDecayPiatBeamFrame, 
	       &PhiDecayPiatBeamFrame);

  if (PhiDecayPiatBeamFrame<0)
    return false;

  /* pi-p scatt */
  G4double pion_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double pion_react_rate = 0.0002548; // 1/mm /* 60mb and LH2 target*/
  
  double p_pi_table[8] = {0.1479, 0.1738, 0.1883, 0.2123, 0.2379, 0.2712, 0.2983, 0.3228};
  
  double pip_cs_table[8] = { 3.4, 5.12, 6.64, 9.8, 15.2, 19.4, 18.3, 14.8};
  double p = momentumDecayPi.mag();

  int index_p=-1;
  double pip_cs;
  for (int i=0; i<8; i++) {
    if (i==0) {
      if (p < p_pi_table[i]) {
	index_p = i;
	pip_cs = pip_cs_table[i];
	break;
      } else if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	index_p = i;
	pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	break;
      }
    } else if (i>=1 && i<=6) {
      if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	index_p = i;
	pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	break;
      }
    } else if (i==7) {
      if (p >= p_pi_table[i]) {
	index_p = i;
	pip_cs = pip_cs_table[i];
	break;
      }
    }
  }
  
  if (index_p<0 || index_p>=8) {
    fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_p : %d", index_p);
    exit(-1);
  }
  pion_react_rate *= pip_cs/30.;
  
  G4double dz = 0.5*mm; //mm
  G4double totalz=0.0;  //mm
  
  G4double ctau_pi=7804.5; /*mm*/
  G4double p_pi = momentumDecayPi.mag();
  G4double m_pi = piMinus->GetPDGMass()/GeV;
  
  G4ThreeVector localPionPos = Pos0;
  int flagTgtType3 = -1;
  
  G4bool flagPionScattering=false;
  G4bool flagPionDecay=false;
  
  int nIte=0;
  while (1) {
    flagTgtType3 = -1;
    totalz += dz;
    localPionPos += momentumDecayPi*dz/momentumDecayPi.mag();
    
    flagPionScattering=false;
    flagPionDecay=false;
    
    flagTgtType3 = getTargetFlag(localPionPos);
    
    if (flagTgtType3 == 0) {
      flagPionScattering  = scatteringCheck(pion_react_rate, dz/mm);
    }
    
    if (flagPionScattering)
      break;
    
    flagPionDecay = decayCheck(ctau_pi, p_pi, m_pi, dz/mm );
    if (flagPionDecay)
      break;
    
    nIte++;
    if (nIte>1000)
      break;
  }

  if (flagPionScattering) {
    //std::cout << "piNscat" << std::endl;
    G4ThreeVector localPiNScatPos=localPionPos;
    G4ThreeVector globalPiNScatPos = geomMan.Local2GlobalPos(TgtId, localPiNScatPos);
    int scatDistFlag = 0;
    G4ParticleDefinition* scatParticle;
    if ( flagTgtType3 == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return false;
    }


    //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
    Kinema3Resonance PiNScat;
    PiNScat = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			       scatParticle->GetPDGMass()/GeV,
			       0.0 , piMinus->GetPDGMass()/GeV,
			       scatParticle->GetPDGMass()/GeV,
			       piMinus->GetPDGMass()/GeV,
			       0.0, momentumDecayPi.mag(), 0.0, scatDistFlag);
    
    double mom[3];

    /* scattered pion */
    G4double Energy_scatPi = PiNScat.GetEnergy(4);
    //G4double momentum_scatPi = PiNScat.GetMomentum(4);
    PiNScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatPi( mom[1], mom[2], mom[0]);
    
    momentumScatPi.rotateY(ThetaDecayPiatBeamFrame*deg);
    momentumScatPi.rotateZ(PhiDecayPiatBeamFrame*deg);
    
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = PiNScat.GetEnergy(5);
    //G4double momentum_scatPart = PiNScat.GetMomentum(5);
    PiNScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaDecayPiatBeamFrame*deg);
    momentumScatPart.rotateZ(PhiDecayPiatBeamFrame*deg);
    
    /*
      std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
      << momentumDecayPi.y() << ", " << momentumDecayPi.z()
      << ")" << std::endl;
      
      std::cout << "ScatPi = ( " << momentumScatPi.x() << ", " 
      << momentumScatPi.y() << ", " << momentumScatPi.z()
      << ")" << std::endl;
      
      std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
      << momentumScatPart.y() << ", " << momentumScatPart.z()
      << ")" << std::endl;
      
      std::cout << "Delta = ( " 
      << momentumDecayPi.x()-momentumScatPi.x()-momentumScatPart.x()
      << ", " 
      << momentumDecayPi.y()-momentumScatPi.y()-momentumScatPart.y()
      << ", " 
      << momentumDecayPi.z()-momentumScatPi.z()-momentumScatPart.z()
      << ")" << std::endl;
    */
    
    /* decay neutron */
    /*
    particleGun->SetParticleDefinition(neutron);
    G4ThreeVector gloMomDecayN = geomMan.Local2GlobalDir(TgtId, momentumDecayN);
    particleGun->SetParticleMomentumDirection(gloMomDecayN);
    particleGun->SetParticleEnergy((Energy_decayN - neutron->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    */    
    
    /*scat pi-*/    
    particleGun->SetParticleDefinition(piMinus);
    G4ThreeVector gloMomentumScatPi = geomMan.Local2GlobalDir(TgtId, momentumScatPi);
    particleGun->SetParticleMomentumDirection(gloMomentumScatPi);
    particleGun->SetParticleEnergy((Energy_scatPi - piMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalPiNScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalPiNScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    //anaMan_->SetPrimaryVertex(localVertexPos);
    //anaMan_->SetHypBeamMomentum(momentumSigma);
    //anaMan_->SetScatMesonMomentum(momentumKPlus);
    //anaMan_->SetThetaMeson(ThetaK);
    //anaMan_->SetPhiMeson(PhiK);
    //anaMan_->SetThetaMesonCM(ThetaKCM);
    //anaMan_->SetPhiMesonCM(PhiKCM);
    //anaMan_->SetDecayPos(localDecayPos);
    anaMan_->SetPiNScatPos(localPiNScatPos);
    anaMan_->SetDecayFlag();
    anaMan_->SetPiNScatFlag();
    anaMan_->SetPiNScatTarget(flagTgtType3);	
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    //anaMan_->SetFlightLengthInTarget(dxInH, 0);
    //anaMan_->SetDecayPiMomentum(momentumDecayPi);
    anaMan_->SetDecayPiMomentum(momentumScatPi); // momentum after scattering
    //anaMan_->SetDecayNucleonMomentum(momentumDecayN);
    //anaMan_->SetBeamMomentum(beammom);

    return true;
  }

  return false;

}


bool CFTPrimaryGeneratorAction::LambdaP_Scattering(G4Event* anEvent, G4ThreeVector momentumLambda, G4ThreeVector Pos0)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");

  G4ParticleDefinition* lambda  = particleTable->FindParticle("lambda");  
  //G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");

  double ThetaLambda, PhiLambda; 
  calcThetaPhi(momentumLambda, 
	       &ThetaLambda, 
	       &PhiLambda);



  /* Lambda p scatt */
  G4double ctau=78.9; /*mm*/
  G4double p_lambda = momentumLambda.mag();
  G4double m_lambda = lambda->GetPDGMass()/GeV;
  //G4double E_lambda = sqrt(p_lambda*p_lambda+m_lambda*m_lambda);
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  G4double react_rate = react_rate_p;

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localLambdaPos = Pos0;
  G4bool flagDecay=false;
  G4bool flagScattering=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localLambdaPos);
    totalx += dx;
    localLambdaPos += momentumLambda*dx/momentumLambda.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomLambda ( " << momentumLambda.x() << ", "  
	      << momentumLambda.y() << ", " << momentumLambda.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localLambdaPos.x() << ", "  
	      << localLambdaPos.y() << ", " << localLambdaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumLambda, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;
    flagScattering=false;

    flagDecay = decayCheck(ctau, p_lambda, m_lambda, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localLambdaPos);

    if (flagTgtType == 0 ) {
      react_rate = react_rate_p;
      dxInH += dx;

      flagScattering  = scatteringCheck(react_rate, dx/mm);

      //p_sigma = calcEnergyDeposit(p_sigma, m_sigma, dx/mm, localSigmaPos, flagTgtType);
      //E_sigma = sqrt(p_sigma*p_sigma + m_sigma*m_sigma);
    }
    /*
    G4cout << "Length : " << totalx << ", p = " << p_lambda*1000. 
	   << " MeV/c, Ekin = "	   << (E_lambda-m_lambda)*1000. 
	   << " MeV" << G4endl;
    */
    if (fabs(p_lambda)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }
    if (flagDecay || flagScattering)
      break;

    nIteration++;
    if (nIteration>5000)
      return false;
  }

  momentumLambda *= p_lambda/momentumLambda.mag();


  if (flagDecay) {
    //std::cout << "decay" << std::endl;
    G4ThreeVector localDecayPos=localLambdaPos;
    G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);


    //static int Ndecay=0;
    int decayMode=-1;
    double randval = RandFlat::shoot(0., 1.);
    if (randval>0.639)
      decayMode = 0; // pi0n
    else
      decayMode = 1; // pi-p

    G4ParticleDefinition* decayNucl;
    G4ParticleDefinition* decayPi;
    if (decayMode==1) {
      decayNucl = particleTable->FindParticle("proton_d");
      decayPi = particleTable->FindParticle("pi-");
    } else if (decayMode==0) {
      decayNucl = particleTable->FindParticle("neutron");
      decayPi = particleTable->FindParticle("pi0");
    }

    Kinema3Resonance LambdaDecay;
    LambdaDecay = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				   0.0,
				   0.0 , decayNucl->GetPDGMass()/GeV,
				   decayPi->GetPDGMass()/GeV,
				   decayNucl->GetPDGMass()/GeV,
				   0.0, p_lambda, 0.0);
    
    double mom[3];

    /* Neucleon */
    double Energy_decayN = LambdaDecay.GetEnergy(4);
    //double momentum_decayN = LambdaDecay.GetMomentum(4);
    LambdaDecay.GetMomentum(4,mom);
    G4ThreeVector momentumDecayN(mom[1], mom[2], mom[0]);
    //double ThetaDecayN = LambdaDecay.GetTheta(4);
    //double PhiDecayN = LambdaDecay.GetPhi(4);
    
    /* pi */
    double Energy_decayPi = LambdaDecay.GetEnergy(5);
    //double momentum_decayPi = LambdaDecay.GetMomentum(5);
    LambdaDecay.GetMomentum(5,mom);
    G4ThreeVector momentumDecayPi( mom[1], mom[2], mom[0]);
    //double ThetaDecayPi = LambdaDecay.GetTheta(5);
    //double PhiDecayPi = LambdaDecay.GetPhi(5);
    
    // Neucleon momentum at beam frame
    momentumDecayN.rotateY(ThetaLambda*deg);
    momentumDecayN.rotateZ(PhiLambda*deg);
    
    // pi momentum at beam frame
    momentumDecayPi.rotateY(ThetaLambda*deg);
    momentumDecayPi.rotateZ(PhiLambda*deg);
    
    
    
    if (decayMode==0) {
      // pi0 n
      bool flagNpScattering = NeutronP_Scattering(anEvent, momentumDecayN, localDecayPos);
      
      if (flagNpScattering) {
	// generate decay pi0
	
	particleGun->SetParticleDefinition(decayPi);
	G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	particleGun->SetParticleEnergy((Energy_decayPi - piMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	anaMan_->SetDecayPos(localDecayPos);
	anaMan_->SetDecayPiMomentum(momentumDecayPi);
	anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	
	return true;
      }
    } else if (decayMode==1) {
      
      bool flagPpScattering = ProtonP_Scattering(anEvent, momentumDecayN, localDecayPos);
      if (flagPpScattering) {
	// generate decay pi-
	
	particleGun->SetParticleDefinition(decayPi);
	G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	particleGun->SetParticleEnergy((Energy_decayPi - piMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	anaMan_->SetDecayPos(localDecayPos);
	anaMan_->SetDecayPiMomentum(momentumDecayPi);
	anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	
	return true;
      }


      bool flagPipScattering = PiMinusP_Scattering(anEvent, momentumDecayPi, localDecayPos);
      
      if (flagPipScattering) {
	// generate decay n
	
	particleGun->SetParticleDefinition(decayNucl);
	G4ThreeVector gloMomDecayN = geomMan.Local2GlobalDir(TgtId, momentumDecayN);
	particleGun->SetParticleMomentumDirection(gloMomDecayN);
	particleGun->SetParticleEnergy((Energy_decayN - neutron->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	anaMan_->SetDecayPos(localDecayPos);
	anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	
	return true;
	
      }
    }

    return false;
    
  } else if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localLambdaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance LambdaScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return false;
    }

    LambdaScat = Kinema3Resonance(lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , 
				  lambda->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  lambda->GetPDGMass()/GeV,
				  0.0, p_lambda, 0.0,
				  scatDistFlag);
    
    double mom[3];

    /* scattered lambda */
    G4double Energy_scatLambda = LambdaScat.GetEnergy(4);
    //G4double momentum_scatLambda = LambdaScat.GetMomentum(4);
    LambdaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatLambda( mom[1], mom[2], mom[0]);
    momentumScatLambda.rotateY(ThetaLambda*deg);
    momentumScatLambda.rotateZ(PhiLambda*deg);
    
    double ThetaScatLambdaCM, PhiScatLambdaCM;
    ThetaScatLambdaCM = LambdaScat.GetThetaCM(1);
    PhiScatLambdaCM   = LambdaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = LambdaScat.GetEnergy(5);
    //G4double momentum_scatPart = LambdaScat.GetMomentum(5);
    LambdaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaLambda*deg);
    momentumScatPart.rotateZ(PhiLambda*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Lambda (" << momentumLambda.x() << ", "
      << momentumLambda.y() << ", " << momentumLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatLambda (" << momentumScatLambda.x() << ", "
      << momentumScatLambda.y() << ", " << momentumScatLambda.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    //anaMan_->SetPrimaryVertex(localVertexPos);
    //anaMan_->SetHypBeamMomentum(momentumLambda);
    anaMan_->SetScatHypMomentum(momentumScatLambda);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    //anaMan_->SetScatMesonMomentum(momentumK);
    //anaMan_->SetThetaMeson(ThetaK);
    //anaMan_->SetPhiMeson(PhiK);
    //anaMan_->SetThetaMesonCM(ThetaKCM);
    //anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatLambdaCM);
    anaMan_->SetPhiScatHypCM(PhiScatLambdaCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag(11);
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    //anaMan_->SetBeamMomentum(beammom);


    /*scatt lambda*/    
    G4ThreeVector gloMomLambda = geomMan.Local2GlobalDir(TgtId, momentumLambda);
    particleGun->SetParticleDefinition(lambda);
    G4ThreeVector gloMomScatLambda = geomMan.Local2GlobalDir(TgtId, momentumScatLambda);
    particleGun->SetParticleMomentumDirection(gloMomScatLambda);
    particleGun->SetParticleEnergy((Energy_scatLambda - lambda->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    return true;

  }

  return false;

}


#if 1
void CFTPrimaryGeneratorAction::GeneratePiKSigmaPlusScat2(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* spiPlus = particleTable->FindParticle("spi+");
  G4ParticleDefinition* piPlus = particleTable->FindParticle("pi+");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  //G4ParticleDefinition* kPlus = particleTable->FindParticle("kaon+");
  G4ParticleDefinition* kPlus = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* sigma = particleTable->FindParticle("sigma+");
  G4ParticleDefinition* ssigma = particleTable->FindParticle("ssigma+");
  G4ParticleDefinition* usigma = particleTable->FindParticle("usigma+");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");
  //G4ParticleDefinition* proton_d = particleTable->FindParticle("proton_d");

  G4double m_neutron = neutron->GetPDGMass()/GeV;
  G4double m_proton = proton->GetPDGMass()/GeV;

  double width = 0.0;
  double mass_sigma = 1.189;
  double beammom = 1.42;

  Kinema3Resonance Sigma;
  int DistFlag=2; // (pi+, K+)

  Sigma = Kinema3Resonance(piPlus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_sigma,
			   kPlus->GetPDGMass()/GeV,
			   mass_sigma, width, beammom, 0.0, DistFlag);

  /* kaon+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumKPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  //G4cout << "pi- (" << momentumPiMinus.x() << ", "
  //<< momentumPiMinus.y() << ", " << momentumPiMinus.z() << ") "
  //<< G4endl;

  /* sigma */
  double Energy_sig = Sigma.GetEnergy(4);
  double momentum_sig = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma( mom[1], mom[2], mom[0]);
  double ThetaSig = Sigma.GetTheta(4);
  double PhiSig = Sigma.GetPhi(4);
  /*
  G4cout << "Sigma (" << momentumSigma.x() << ", "
	 << momentumSigma.y() << ", " << momentumSigma.z() << ") "
	 << G4endl;
  */

  /* pi+ beam */
  double Energy_beam = sqrt(pow(piPlus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam( 0., 0., -beammom);


  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;



  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;

  primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
  primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
					       primary_vertex_y, 
					       primary_vertex_z));
  if (primaryTgtType != 0)
    return;

  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  /* Sigma-p scatt */
  //G4double flength=20.*mm;
  G4double ctau=24.04; /*mm*/
  G4double p_sigma = momentum_sig;
  G4double m_sigma = sigma->GetPDGMass()/GeV;
  G4double E_sigma = 0; // initialization
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_p= 0.01274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;
  double dxInD=0.;

  G4ThreeVector localSigmaPos = localVertexPos;
  G4bool flagDecay=false;
  G4bool flagScattering=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);
    // Multiple scattering
    /*
    if (flagTgtType == 0 || flagTgtType == 1) {
      momentumSigma = calcMultipleScattering(p_sigma, m_sigma, dx/mm, momentumSigma, flagTgtType);

    }
    */
    totalx += dx;
    localSigmaPos += momentumSigma*dx/momentumSigma.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomSigma ( " << momentumSigma.x() << ", "  
	      << momentumSigma.y() << ", " << momentumSigma.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localSigmaPos.x() << ", "  
	      << localSigmaPos.y() << ", " << localSigmaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumSigma, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;
    flagScattering=false;

    flagDecay = decayCheck(ctau, p_sigma, m_sigma, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);

    if (flagTgtType == 0) {
      react_rate = react_rate_p;
      dxInH += dx;

      flagScattering  = scatteringCheck(react_rate, dx/mm);
      p_sigma = calcEnergyDeposit(p_sigma, m_sigma, dx/mm, localSigmaPos, flagTgtType);
      E_sigma = sqrt(p_sigma*p_sigma + m_sigma*m_sigma);
    }
    /*
    G4cout << "Length : " << totalx << ", p = " << p_sigma*1000. 
	   << " MeV/c, Ekin = "	   << (E_sigma-m_sigma)*1000. 
	   << " MeV" << G4endl;
    */
    if (fabs(p_sigma)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }
    if (flagDecay || flagScattering)
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }


  if (flagDecay) {
    //std::cout << "decay" << std::endl;
    G4ThreeVector localDecayPos=localSigmaPos;
    G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    momentumSigma *= p_sigma/momentumSigma.mag();

    if (1) {
      static int Ndecay=0;
      int decayMode=-1;
      decayMode = Ndecay%2;
      Ndecay++;

      G4ParticleDefinition* decayNucl = 0;
      G4ParticleDefinition* decayPi = 0;
      if (decayMode==0) {
	decayNucl = particleTable->FindParticle("neutron");
	decayPi = particleTable->FindParticle("pi+");
      } else if (decayMode==1) {
	decayNucl = particleTable->FindParticle("proton_d");
	decayPi = particleTable->FindParticle("pi0");
      }

      //neutron scattering
      Kinema3Resonance SigmaDecay;
      SigmaDecay = Kinema3Resonance(mass_sigma, 
				    0.0,
				    0.0 , decayNucl->GetPDGMass()/GeV,
				    decayPi->GetPDGMass()/GeV,
				    decayNucl->GetPDGMass()/GeV,
				    0.0, p_sigma, 0.0);


      /* Neucleon */
      double Energy_decayN = SigmaDecay.GetEnergy(4);
      double momentum_decayN = SigmaDecay.GetMomentum(4);
      SigmaDecay.GetMomentum(4,mom);
      G4ThreeVector momentumDecayN(mom[1], mom[2], mom[0]);
      //double ThetaDecayN = SigmaDecay.GetTheta(4);
      //double PhiDecayN = SigmaDecay.GetPhi(4);
      
      /* pi */
      double Energy_decayPi = SigmaDecay.GetEnergy(5);
      double momentum_decayPi = SigmaDecay.GetMomentum(5);
      SigmaDecay.GetMomentum(5,mom);
      G4ThreeVector momentumDecayPi( mom[1], mom[2], mom[0]);
      //double ThetaDecayPi = SigmaDecay.GetTheta(5);
      //double PhiDecayPi = SigmaDecay.GetPhi(5);
      
      // Neucleon momentum at beam frame
      momentumDecayN.rotateY(ThetaSig*deg);
      momentumDecayN.rotateZ(PhiSig*deg);

      // pi momentum at beam frame
      momentumDecayPi.rotateY(ThetaSig*deg);
      momentumDecayPi.rotateZ(PhiSig*deg);

      /*
      std::cout << "SigmaMom = ( " << momentumSigma.x() << ", " 
		<< momentumSigma.y() << ", " << momentumSigma.z()
		<< ")" << std::endl;
      std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
		<< momentumDecayN.y() << ", " << momentumDecayN.z()
		<< ")" << std::endl;
      std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
		<< momentumDecayPi.y() << ", " << momentumDecayPi.z()
		<< ")" << std::endl;

      std::cout << "Delta = ( " 
		<< momentumSigma.x()-momentumDecayN.x()-momentumDecayPi.x()
		<< ", " 
		<< momentumSigma.y()-momentumDecayN.y()-momentumDecayPi.y() 
		<< ", " 
		<< momentumSigma.z()-momentumDecayN.z()-momentumDecayPi.z()
		<< ")" << std::endl;
      */
      double ThetaDecayNatBeamFrame, PhiDecayNatBeamFrame; 
      calcThetaPhi(momentumDecayN, 
		   &ThetaDecayNatBeamFrame, 
		   &PhiDecayNatBeamFrame);

      double ThetaDecayPiatBeamFrame, PhiDecayPiatBeamFrame; 
      calcThetaPhi(momentumDecayPi, 
		   &ThetaDecayPiatBeamFrame, 
		   &PhiDecayPiatBeamFrame);
      /*
      std::cout << "ThetaDecayNatBeamFrame: " << ThetaDecayNatBeamFrame
		<< " PhiDecayNatBeamFrame: "  << PhiDecayNatBeamFrame
		<< std::endl;

      std::cout << "ThetaDecayPiatBeamFrame: " << ThetaDecayPiatBeamFrame
		<< " PhiDecayPiatBeamFrame: "  << PhiDecayPiatBeamFrame
		<< std::endl;
      */
      if (PhiDecayNatBeamFrame<0)
	return;
      if (PhiDecayPiatBeamFrame<0)
	return;


      bool flagProtonScattering=false;
      bool flagNeutronScattering=false;
      bool flagPionScattering=false;


      if (decayMode==0) {
	/* n-p scatt */
	G4double neutron_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
	
	double Ekin_n[30] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
			     110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
		   210., 220., 230., 240., 250., 260., 270., 280., 290., 300.};
	
	double np_cs_table[30] = { 951.827, 488.477, 312.839, 223.243, 170.507,
				   136.696, 113.728,  97.449,  85.517,  76.527,
				   69.593,  64.133,  59.756,  56.191,  53.244,
				   50.779,  48.691,  46.905,  45.364,  44.022,
				   42.844,  41.803,  40.876,  40.046,  39.298,
				   38.620,  38.002,  37.436,  36.915,  36.432};
	
	double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_neutron*m_neutron)
		       -m_neutron)*1000. ;// MeV
	
	
	int index_Ekin=-1;
	double np_cs;
	for (int i=0; i<30; i++) {
	  if (i==0) {
	    if (Ekin < Ekin_n[i]) {
	      index_Ekin = i;
	      np_cs = np_cs_table[i];
	      break;
	    } else if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	      index_Ekin = i;
	      np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	      break;
	    }
	  } else if (i>=1 && i<=28) {
	    if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	      index_Ekin = i;
	      np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	      break;
	    }
	  } else if (i==29) {
	    if (Ekin >= Ekin_n[i]) {
	      index_Ekin = i;
	      np_cs = np_cs_table[i];
	      break;
	    }
	  }
	}
	
	if (index_Ekin<0 || index_Ekin>=30) {
	  fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
	  exit(-1);
	}
	neutron_react_rate *= np_cs/30.;
	
	G4double dy = 0.5*mm; //mm
	G4double totaly=0.0;  //mm
	
	G4ThreeVector localNeutronPos = localDecayPos;
	int flagTgtType2 = -1;
	
	int nIte=0;
	while (1) {
	  flagTgtType2 = -1;
	  totaly += dy;
	  localNeutronPos += momentumDecayN*dy/momentumDecayN.mag();
	  
	  flagNeutronScattering=false;
	  
	  flagTgtType2 = getTargetFlag(localNeutronPos);
	  
	  if (flagTgtType2 == 0 || flagTgtType2 == 1) {
	    flagNeutronScattering  = scatteringCheck(neutron_react_rate, dy/mm);
	  }
	  
	  if (flagNeutronScattering)
	    break;
	  
	  nIte++;
	  if (nIte>1000)
	    break;
	}


	/* pi+p scatt */
	G4double pion_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/

	double p_pi_table[20] = {0.137, 0.158, 0.176, 0.189, 0.204,
				0.216, 0.228, 0.238, 0.251, 0.260,
				0.2718, 0.2903, 0.298, 0.3110, 0.327,
				0.3406, 0.3529, 0.3710, 0.3850, 0.401};

	double pip_cs_table[20] = { 17.6, 26.09, 37.3, 50.3, 71.4, 
				   91.0, 122., 139.6, 167.0, 189.0,
				   187.0, 202.9, 194.2, 178.3, 156.8,
				   134.81, 125.6, 96.2, 85.2, 73.8};
	double p = momentum_decayPi;
	int index_p=-1;
	double pip_cs;
	for (int i=0; i<20; i++) {
	  if (i==0) {
	    if (p < p_pi_table[i]) {
	      index_p = i;
	      pip_cs = pip_cs_table[i];
	      break;
	    } else if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	      index_p = i;
	      pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	      break;
	    }
	  } else if (i>=1 && i<=18) {
	    if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	      index_p = i;
	      pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	      break;
	    }
	  } else if (i==19) {
	    if (p >= p_pi_table[i]) {
	      index_p = i;
	      pip_cs = pip_cs_table[i];
	      break;
	    }
	  }
	}
	
	if (index_p<0 || index_p>=20) {
	  fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_p : %d", index_p);
	  exit(-1);
	}
	pion_react_rate *= pip_cs/30.;

	G4double dz = 0.5*mm; //mm
	G4double totalz=0.0;  //mm
	
	G4double ctau_pi=7804.5; /*mm*/
	G4double p_pi = momentumDecayPi.mag();
	G4double m_pi = piPlus->GetPDGMass()/GeV;
	
	G4ThreeVector localPionPos = localDecayPos;
	int flagTgtType3 = -1;
	
	G4bool flagPionDecay=false;
	
	nIte=0;
	while (1) {
	  flagTgtType3 = -1;
	  totalz += dz;
	  localPionPos += momentumDecayPi*dz/momentumDecayPi.mag();
	  
	  flagPionScattering=false;
	  flagPionDecay=false;
	  
	  flagTgtType3 = getTargetFlag(localPionPos);
	  
	  if (flagTgtType3 == 0 || flagTgtType3 == 1) {
	    flagPionScattering  = scatteringCheck(pion_react_rate, dz/mm);
	  }
	  
	  if (flagPionScattering)
	    break;
	  
	  flagPionDecay = decayCheck(ctau_pi, p_pi, m_pi, dz/mm );
	  if (flagPionDecay)
	    break;
	  
	  nIte++;
	  if (nIte>1000)
	    break;
	}

	if (flagNeutronScattering) {
	  //std::cout << "NNscat" << std::endl;
	  G4ThreeVector localNNScatPos=localNeutronPos;
	  G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
	  int scatDistFlag = 0;
	  G4ParticleDefinition* scatParticle;
	  if ( flagTgtType2 == 0 ) {
	    // inside the LH2 target;
	    scatParticle = particleTable->FindParticle("proton");
	    scatDistFlag = 100+index_Ekin;
	  } else {
	    return;
	  }
	  //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	  Kinema3Resonance NNScat;
	  NNScat = Kinema3Resonance(decayNucl->GetPDGMass()/GeV, 
				    scatParticle->GetPDGMass()/GeV,
				    0.0 , decayNucl->GetPDGMass()/GeV,
				    scatParticle->GetPDGMass()/GeV,
				    decayNucl->GetPDGMass()/GeV,
				    0.0, momentum_decayN, 0.0, scatDistFlag);
	  
	  /* scattered Necleon */
	  G4double Energy_scatN = NNScat.GetEnergy(4);
	  //G4double momentum_scatN = NNScat.GetMomentum(4);
	  NNScat.GetMomentum(4,mom);
	  G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);
	  
	  momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /* scattered Proton or Deuteron*/
	  G4double Energy_scatPart = NNScat.GetEnergy(5);
	  //G4double momentum_scatPart = NNScat.GetMomentum(5);
	  NNScat.GetMomentum(5,mom);
	  G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	  momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /*
	    std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
	    << momentumDecayN.y() << ", " << momentumDecayN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
	    << momentumScatN.y() << ", " << momentumScatN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
	    << momentumScatPart.y() << ", " << momentumScatPart.z()
	    << ")" << std::endl;
	    
	    std::cout << "Delta = ( " 
	    << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
	    << ", " 
	    << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
	    << ", " 
	    << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
	    << ")" << std::endl;
	  */
	  
	  particleGun->SetParticleDefinition(kPlus);
	  G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
	  particleGun->SetParticleMomentumDirection(gloMomKPlus);
	  particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* beam */
	  particleGun->SetParticleDefinition(spiPlus);
	  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	  particleGun->SetParticleMomentumDirection(gloMomBeam);
	  particleGun->SetParticleEnergy((Energy_beam - spiPlus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*sigma*/    
	  particleGun->SetParticleDefinition(ssigma);
	  G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	  particleGun->SetParticleMomentumDirection(-gloMomSigma);
	  // Energy_sig is original Energy without energy deposit
	  particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* decay pi */
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	  particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	  particleGun->SetParticleEnergy((Energy_decayPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scat proton*/    
	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
	  particleGun->SetParticleMomentumDirection(gloMomentumScatN);
	  particleGun->SetParticleEnergy((Energy_scatN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scatt proton*/    
	  particleGun->SetParticleDefinition(scatParticle);
	  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	  particleGun->SetParticleMomentumDirection(gloMomScatPart);
	  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  
	  anaMan_->SetPrimaryVertex(localVertexPos);
	  anaMan_->SetHypBeamMomentum(momentumSigma);
	  anaMan_->SetScatMesonMomentum(momentumKPlus);
	  anaMan_->SetThetaMeson(ThetaK);
	  anaMan_->SetPhiMeson(PhiK);
	  anaMan_->SetThetaMesonCM(ThetaKCM);
	  anaMan_->SetPhiMesonCM(PhiKCM);
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetNNScatPos(localNNScatPos);
	  anaMan_->SetDecayFlag(decayMode);
	  anaMan_->SetNNScatFlag();
	  anaMan_->SetNNScatTarget(flagTgtType2);	
	  anaMan_->SetScatProtonMomentum(momentumScatPart);
	  anaMan_->SetFlightLengthInTarget(dxInH, dxInD);
	  anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  anaMan_->SetBeamMomentum(beammom);
	  
	  return;
	}

	if (flagPionScattering) {
	  //std::cout << "piNscat" << std::endl;
	  G4ThreeVector localPiNScatPos=localPionPos;
	  G4ThreeVector globalPiNScatPos = geomMan.Local2GlobalPos(TgtId, localPiNScatPos);
	  int scatDistFlag = 0;
	  G4ParticleDefinition* scatParticle;
	  if ( flagTgtType3 == 0 ) {
	    // inside the LH2 target;
	    scatParticle = particleTable->FindParticle("proton");
	    scatDistFlag = 0;
	  } else {
	    return;
	  }
	  //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	  Kinema3Resonance PiNScat;
	  PiNScat = Kinema3Resonance(decayPi->GetPDGMass()/GeV, 
				     scatParticle->GetPDGMass()/GeV,
				     0.0 , decayPi->GetPDGMass()/GeV,
				     scatParticle->GetPDGMass()/GeV,
				     decayPi->GetPDGMass()/GeV,
				     0.0, momentum_decayPi, 0.0, scatDistFlag);
	  
	  /* scattered pion */
	  G4double Energy_scatPi = PiNScat.GetEnergy(4);
	  //G4double momentum_scatPi = PiNScat.GetMomentum(4);
	  PiNScat.GetMomentum(4,mom);
	  G4ThreeVector momentumScatPi( mom[1], mom[2], mom[0]);
	  
	  momentumScatPi.rotateY(ThetaDecayPiatBeamFrame*deg);
	  momentumScatPi.rotateZ(PhiDecayPiatBeamFrame*deg);
	  
	  /* scattered Proton or Deuteron*/
	  G4double Energy_scatPart = PiNScat.GetEnergy(5);
	  //G4double momentum_scatPart = PiNScat.GetMomentum(5);
	  PiNScat.GetMomentum(5,mom);
	  G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	  momentumScatPart.rotateY(ThetaDecayPiatBeamFrame*deg);
	  momentumScatPart.rotateZ(PhiDecayPiatBeamFrame*deg);
	  
	  /*
	    std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
	    << momentumDecayPi.y() << ", " << momentumDecayPi.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPi = ( " << momentumScatPi.x() << ", " 
	    << momentumScatPi.y() << ", " << momentumScatPi.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
	    << momentumScatPart.y() << ", " << momentumScatPart.z()
	    << ")" << std::endl;
	    
	    std::cout << "Delta = ( " 
	    << momentumDecayPi.x()-momentumScatPi.x()-momentumScatPart.x()
	    << ", " 
	    << momentumDecayPi.y()-momentumScatPi.y()-momentumScatPart.y()
	    << ", " 
	    << momentumDecayPi.z()-momentumScatPi.z()-momentumScatPart.z()
	    << ")" << std::endl;
	  */
	  
	  particleGun->SetParticleDefinition(kPlus);
	  G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
	  particleGun->SetParticleMomentumDirection(gloMomKPlus);
	  particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* beam */
	  particleGun->SetParticleDefinition(spiPlus);
	  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	  particleGun->SetParticleMomentumDirection(gloMomBeam);
	  particleGun->SetParticleEnergy((Energy_beam - spiPlus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*sigma*/    
	  particleGun->SetParticleDefinition(ssigma);
	  G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	  particleGun->SetParticleMomentumDirection(-gloMomSigma);
	  // Energy_sig is original Energy without energy deposit
	  particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* decay neucleon */
	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomDecayN = geomMan.Local2GlobalDir(TgtId, momentumDecayN);
	  particleGun->SetParticleMomentumDirection(gloMomDecayN);
	  particleGun->SetParticleEnergy((Energy_decayN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scat pi*/    
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomentumScatPi = geomMan.Local2GlobalDir(TgtId, momentumScatPi);
	  particleGun->SetParticleMomentumDirection(gloMomentumScatPi);
	  particleGun->SetParticleEnergy((Energy_scatPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalPiNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scatt proton*/    
	  particleGun->SetParticleDefinition(scatParticle);
	  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	  particleGun->SetParticleMomentumDirection(gloMomScatPart);
	  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalPiNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  
	  anaMan_->SetPrimaryVertex(localVertexPos);
	  anaMan_->SetHypBeamMomentum(momentumSigma);
	  anaMan_->SetScatMesonMomentum(momentumKPlus);
	  anaMan_->SetThetaMeson(ThetaK);
	  anaMan_->SetPhiMeson(PhiK);
	  anaMan_->SetThetaMesonCM(ThetaKCM);
	  anaMan_->SetPhiMesonCM(PhiKCM);
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetPiNScatPos(localPiNScatPos);
	  anaMan_->SetDecayFlag(decayMode);
	  anaMan_->SetPiNScatFlag();
	  anaMan_->SetPiNScatTarget(flagTgtType3);	
	  anaMan_->SetScatProtonMomentum(momentumScatPart);
	  anaMan_->SetFlightLengthInTarget(dxInH, dxInD);
	  //anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayPiMomentum(momentumScatPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  anaMan_->SetBeamMomentum(beammom);
	  
	  return;
	}
      } else if (decayMode==1) {

	G4double proton_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
	
	double Ekin_p[14] = {18.2, 19.8, 25.63, 30.14, 39.4, 68.3, 95., 98., 118., 142.,
			     147., 172., 250., 312.};
	
	double pp_cs_table[14] = { 351.8, 314.1, 238.7, 188.4, 138.2, 81.6,
				   56.5, 56.5, 52.7, 52.7, 51.5, 50.2, 50.2, 46.4};
	
	double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_proton*m_proton)
		       -m_proton)*1000. ;// MeV
	
	int index_Ekin=-1;
	double pp_cs;
	for (int i=0; i<14; i++) {
	  if (i==0) {
	    if (Ekin < Ekin_p[i]) {
	      index_Ekin = i;
	      pp_cs = pp_cs_table[i];
	      break;
	    } else if (Ekin >= Ekin_p[i] && Ekin < Ekin_p[i+1]) {
	      index_Ekin = i;
	      pp_cs = (pp_cs_table[i+1]*(Ekin-Ekin_p[i]) + pp_cs_table[i]*(Ekin_p[i+1]-Ekin))/(Ekin_p[i+1]-Ekin_p[i]);
	      break;
	    }
	  } else if (i>=1 && i<=12) {
	    if (Ekin >= Ekin_p[i] && Ekin < Ekin_p[i+1]) {
	      index_Ekin = i;
	      pp_cs = (pp_cs_table[i+1]*(Ekin-Ekin_p[i]) + pp_cs_table[i]*(Ekin_p[i+1]-Ekin))/(Ekin_p[i+1]-Ekin_p[i]);
	      break;
	    }
	  } else if (i==13) {
	    if (Ekin >= Ekin_p[i]) {
	      index_Ekin = i;
	      pp_cs = pp_cs_table[i];
	      break;
	    }
	  }
	}
	
	if (index_Ekin<0 || index_Ekin>=14) {
	  fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
	  exit(-1);
	}
	proton_react_rate *= pp_cs/30.;
	
	G4double dy = 0.5*mm; //mm
	G4double totaly=0.0;  //mm
	
	G4ThreeVector localProtonPos = localDecayPos;
	int flagTgtType2 = -1;
	
	int nIte=0;
	while (1) {
	  flagTgtType2 = -1;
	  totaly += dy;
	  localProtonPos += momentumDecayN*dy/momentumDecayN.mag();
	  
	  flagProtonScattering=false;
	  
	  flagTgtType2 = getTargetFlag(localProtonPos);
	  
	  if (flagTgtType2 == 0 || flagTgtType2 == 1) {
	    flagProtonScattering  = scatteringCheck(proton_react_rate, dy/mm);
	  }
	  
	  if (flagProtonScattering)
	    break;
	  
	  nIte++;
	  if (nIte>1000)
	    break;
	}

	if (flagProtonScattering) {
	  //std::cout << "*********************** pp scat ***********************" << std::endl;
	  G4ThreeVector localNNScatPos=localProtonPos;
	  G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
	  int scatDistFlag = 0;
	  G4ParticleDefinition* scatParticle;
	  if ( flagTgtType2 == 0 ) {
	    // inside the LH2 target;
	    scatParticle = particleTable->FindParticle("proton");
	    scatDistFlag = 1;  // flat distribution
	  } else {
	    return;
	  }
	  //std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	  Kinema3Resonance NNScat;
	  NNScat = Kinema3Resonance(decayNucl->GetPDGMass()/GeV, 
				    scatParticle->GetPDGMass()/GeV,
				    0.0 , decayNucl->GetPDGMass()/GeV,
				    scatParticle->GetPDGMass()/GeV,
				    decayNucl->GetPDGMass()/GeV,
				    0.0, momentum_decayN, 0.0, scatDistFlag);
	  
	  /* scattered Necleon */
	  G4double Energy_scatN = NNScat.GetEnergy(4);
	  //G4double momentum_scatN = NNScat.GetMomentum(4);
	  NNScat.GetMomentum(4,mom);
	  G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);
	  
	  momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /* scattered Proton or Deuteron*/
	  G4double Energy_scatPart = NNScat.GetEnergy(5);
	  //G4double momentum_scatPart = NNScat.GetMomentum(5);
	  NNScat.GetMomentum(5,mom);
	  G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	  momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
	  momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);
	  
	  /*
	    std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
	    << momentumDecayN.y() << ", " << momentumDecayN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
	    << momentumScatN.y() << ", " << momentumScatN.z()
	    << ")" << std::endl;
	    
	    std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
	    << momentumScatPart.y() << ", " << momentumScatPart.z()
	    << ")" << std::endl;
	    
	    std::cout << "Delta = ( " 
	    << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
	    << ", " 
	    << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
	    << ", " 
	    << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
	    << ")" << std::endl;
	  */
	  
	  particleGun->SetParticleDefinition(kPlus);
	  G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
	  particleGun->SetParticleMomentumDirection(gloMomKPlus);
	  particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* beam */
	  particleGun->SetParticleDefinition(spiPlus);
	  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	  particleGun->SetParticleMomentumDirection(gloMomBeam);
	  particleGun->SetParticleEnergy((Energy_beam - spiPlus->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalVertexPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*sigma*/    
	  particleGun->SetParticleDefinition(ssigma);
	  G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	  particleGun->SetParticleMomentumDirection(-gloMomSigma);
	  // Energy_sig is original Energy without energy deposit
	  particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /* decay pi */
	  particleGun->SetParticleDefinition(decayPi);
	  G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	  particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	  particleGun->SetParticleEnergy((Energy_decayPi - decayPi->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalDecayPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scat proton*/    
	  //G4cout << "pp scat" << G4endl;

	  particleGun->SetParticleDefinition(decayNucl);
	  G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
	  particleGun->SetParticleMomentumDirection(gloMomentumScatN);
	  particleGun->SetParticleEnergy((Energy_scatN - decayNucl->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  /*scatt proton*/    
	  particleGun->SetParticleDefinition(scatParticle);
	  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	  particleGun->SetParticleMomentumDirection(gloMomScatPart);
	  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	  particleGun->SetParticlePosition(globalNNScatPos);
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  
	  anaMan_->SetPrimaryVertex(localVertexPos);
	  anaMan_->SetHypBeamMomentum(momentumSigma);
	  anaMan_->SetScatMesonMomentum(momentumKPlus);
	  anaMan_->SetThetaMeson(ThetaK);
	  anaMan_->SetPhiMeson(PhiK);
	  anaMan_->SetThetaMesonCM(ThetaKCM);
	  anaMan_->SetPhiMesonCM(PhiKCM);
	  anaMan_->SetDecayPos(localDecayPos);
	  anaMan_->SetNNScatPos(localNNScatPos);
	  anaMan_->SetDecayFlag(decayMode);
	  anaMan_->SetNNScatFlag();
	  anaMan_->SetNNScatTarget(flagTgtType2);	
	  anaMan_->SetScatProtonMomentum(momentumScatPart);
	  anaMan_->SetFlightLengthInTarget(dxInH, dxInD);
	  anaMan_->SetDecayPiMomentum(momentumDecayPi);
	  anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	  anaMan_->SetBeamMomentum(beammom);
	  
	  return;
	}

      }
    }

    static int nDecay = 0;
    const int prescale = 100;
    nDecay++;


    if (nDecay%prescale != 0)
      return;

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetDecayPos(localSigmaPos);
    anaMan_->SetDecayEventFlag();
    anaMan_->SetDecayFlag();
    anaMan_->SetFlightLengthInTarget(dxInH, dxInD);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /* beam */
    particleGun->SetParticleDefinition(spiPlus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*sigma*/    
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*decaying sigma*/    
    particleGun->SetParticleDefinition(usigma);
    particleGun->SetParticleMomentumDirection(gloMomSigma);
    // E_sigma is energy after energy deposit
    particleGun->SetParticleEnergy((E_sigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);
  } else if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localSigmaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance SigmaScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return;
    }

    
    SigmaScat = Kinema3Resonance(mass_sigma, 
				 scatParticle->GetPDGMass()/GeV,
				 0.0 , mass_sigma,
				 scatParticle->GetPDGMass()/GeV,
				 mass_sigma, 0.0, p_sigma, 0.0,
				 scatDistFlag);

    
    /* scattered Sigma */
    G4double Energy_scatSigma = SigmaScat.GetEnergy(4);
    //G4double momentum_scatSigma = SigmaScat.GetMomentum(4);
    SigmaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatSigma( mom[1], mom[2], mom[0]);
    momentumScatSigma.rotateY(ThetaSig*deg);
    momentumScatSigma.rotateZ(PhiSig*deg);
    
    double ThetaScatSigCM, PhiScatSigCM;
    ThetaScatSigCM = SigmaScat.GetThetaCM(1);
    PhiScatSigCM = SigmaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = SigmaScat.GetEnergy(5);
    //G4double momentum_scatPart = SigmaScat.GetMomentum(5);
    SigmaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaSig*deg);
    momentumScatPart.rotateZ(PhiSig*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatHypMomentum(momentumScatSigma);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatSigCM);
    anaMan_->SetPhiScatHypCM(PhiScatSigCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, dxInD);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(spiPlus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*sigma*/    
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt sigma*/    
    particleGun->SetParticleDefinition(sigma);
    G4ThreeVector gloMomScatSigma = geomMan.Local2GlobalDir(TgtId, momentumScatSigma);
    particleGun->SetParticleMomentumDirection(gloMomScatSigma);
    particleGun->SetParticleEnergy((Energy_scatSigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    

    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  }
}



#endif

void CFTPrimaryGeneratorAction::GeneratePiKBG(G4Event* anEvent)
{

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* part;

  static int event_num = 0;
  double mom = 0;

  if (event_num%2 == 0) {
    std::cout << "pi+" << std::endl;
    part = particleTable->FindParticle("pi+");
    mom = (G4double)RandFlat::shoot(0.3, 1.4);
  } else {
    std::cout << "proton" << std::endl;
    part = particleTable->FindParticle("proton");
    mom = (G4double)RandFlat::shoot(0.7, 1.4);
  }

  event_num++;

  double Energy = sqrt(mom*mom+
		       (part->GetPDGMass()/GeV)*(part->GetPDGMass()/GeV));

  double cost = (G4double)RandFlat::shoot(0.64, 1.); // 0 - 50 degree
  double theta = acos(cost)*Rad2Deg;
  double phi =  (G4double)RandFlat::shoot(0., 360.);
  G4ThreeVector momentum(sin(theta*Deg2Rad)*cos(phi*Deg2Rad),
			 sin(theta*Deg2Rad)*sin(phi*Deg2Rad),
			 cos(theta*Deg2Rad));


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);


  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  particleGun->SetParticleDefinition(part);
  G4ThreeVector gloMom = geomMan.Local2GlobalDir(TgtId, momentum);
  particleGun->SetParticleMomentumDirection(gloMom);
  particleGun->SetParticleEnergy((Energy - part->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalVertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetScatMesonMomentum(momentum);
  anaMan_->SetThetaMeson(theta);
  anaMan_->SetPhiMeson(phi);

}

void CFTPrimaryGeneratorAction::GenerateUniform(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonPlus;

  kaonPlus = particleTable->FindParticle("skaon+");

  //double mom = 0.8;
  double mom = 1.3;
  double Energy_k = sqrt(mom*mom+
			 (kaonPlus->GetPDGMass()/GeV)*(kaonPlus->GetPDGMass()/GeV));

  double cost = (G4double)RandFlat::shoot(sqrt(3.)/2., 1.);
  double theta = acos(cost)*Rad2Deg;
  //double theta = 0.;
  double phi =  (G4double)RandFlat::shoot(0., 360.);
  G4ThreeVector momentumKaonPlus(sin(theta*Deg2Rad)*cos(phi*Deg2Rad),
				 sin(theta*Deg2Rad)*sin(phi*Deg2Rad),
				 cos(theta*Deg2Rad));


  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;
  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
  } while (primaryTgtType != 0 && primaryTgtType != 1);


  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  particleGun->SetParticleDefinition(kaonPlus);
  G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKaonPlus);
  particleGun->SetParticleMomentumDirection(gloMomKPlus);
  particleGun->SetParticleEnergy((Energy_k - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalVertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetScatMesonMomentum(momentumKaonPlus);
  anaMan_->SetThetaMeson(theta);
  anaMan_->SetPhiMeson(phi);


}


void CFTPrimaryGeneratorAction::GeneratePiKSigmaScat2(G4Event* anEvent)
{
  double mom[3];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* spiMinus = particleTable->FindParticle("spi-");
  G4ParticleDefinition* piMinus = particleTable->FindParticle("pi-");
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  //G4ParticleDefinition* kPlus = particleTable->FindParticle("kaon+");
  G4ParticleDefinition* kPlus = particleTable->FindParticle("skaon+");
  G4ParticleDefinition* sigma = particleTable->FindParticle("sigma-");
  G4ParticleDefinition* ssigma = particleTable->FindParticle("ssigma+");
  G4ParticleDefinition* usigma = particleTable->FindParticle("usigma-");
  G4ParticleDefinition* neutron = particleTable->FindParticle("neutron");

  G4double m_neutron = neutron->GetPDGMass()/GeV;

  double width = 0.0;
  double mass_sigma = 1.197;
  double beammom = 1.32;

  Kinema3Resonance Sigma;
  int DistFlag=3; // 1.3GeV/c (pi-, K+)

  Sigma = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
			   proton->GetPDGMass()/GeV,
			   0.0 , mass_sigma,
			   kPlus->GetPDGMass()/GeV,
			   mass_sigma, width, beammom, 0.0, DistFlag);

  /* kaon+ */
  double Energy_k = Sigma.GetEnergy(5);
  //double momentum_k = Sigma.GetMomentum(5);
  Sigma.GetMomentum(5,mom);
  G4ThreeVector momentumKPlus(mom[1], mom[2], mom[0]);
  double ThetaK = Sigma.GetTheta(5);
  double PhiK = Sigma.GetPhi(5);
  double ThetaKCM = 180.-Sigma.GetThetaCM(1);
  double PhiKCM = Sigma.GetPhiCM(1);

  //G4cout << "pi- (" << momentumPiMinus.x() << ", "
  //<< momentumPiMinus.y() << ", " << momentumPiMinus.z() << ") "
  //<< G4endl;

  /* sigma */
  double Energy_sig = Sigma.GetEnergy(4);
  double momentum_sig = Sigma.GetMomentum(4);
  Sigma.GetMomentum(4,mom);
  G4ThreeVector momentumSigma( mom[1], mom[2], mom[0]);
  double ThetaSig = Sigma.GetTheta(4);
  double PhiSig = Sigma.GetPhi(4);
  /*
  G4cout << "Sigma (" << momentumSigma.x() << ", "
	 << momentumSigma.y() << ", " << momentumSigma.z() << ") "
	 << G4endl;
  */
  /* pi- beam */
  double Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  G4ThreeVector momentumBeam( 0., 0., -beammom);



  if (!(ThetaK>0.&&ThetaK<50.)) 
    return;

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  G4int primaryTgtType = -1;

  // generate primary vertex position
  // This is continued untill the vertex position is LH2
  /*
  do {
    primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
    primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
    primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
						 primary_vertex_y, 
						 primary_vertex_z));
    //} while (primaryTgtType != 0 && primaryTgtType != 1);
  } while (primaryTgtType != 0);
  */

  primary_vertex_x = (G4double)RandGauss::shoot(0., 5.)*mm;
  primary_vertex_y = (G4double)RandGauss::shoot(0., 5.)*mm;
  //primary_vertex_x = (G4double)RandGauss::shoot(0., 8.)*mm;
  //primary_vertex_y = (G4double)RandGauss::shoot(0., 8.)*mm;
  primary_vertex_z = (G4double)RandFlat::shoot(-150., 150.)*mm;
  primaryTgtType = getTargetFlag(G4ThreeVector(primary_vertex_x, 
					       primary_vertex_y, 
					       primary_vertex_z));
  if (primaryTgtType != 0)
    return;
  
  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  /* Sigma-p scatt */
  //G4double flength=20.*mm;
  G4double ctau=44.34; /*mm*/
  G4double p_sigma = momentum_sig;
  G4double m_sigma = sigma->GetPDGMass()/GeV;
  G4double E_sigma = 0; // initialization
  G4double react_rate_p = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_p= 0.01274; // 1/mm /* temporary */
  G4double react_rate = react_rate_p;

  G4double react_rate_lp_conv = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_lp_conv = 0.01274; // 1/mm /* temporary */
  G4double react_rate_sp_conv = 0.0001274; // 1/mm /* 30mb and LH2 target*/
  //G4double react_rate_sp_conv = 0.01274; // 1/mm /* temporary */

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;

  G4ThreeVector localSigmaPos = localVertexPos;
  G4bool flagDecay=false;
  G4bool flagScattering=false;
  G4bool flagLambdaPConv=false;
  G4bool flagSigmaPConv=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);
    // Multiple scattering
    /*
    if (flagTgtType == 0 || flagTgtType == 1) {
      momentumSigma = calcMultipleScattering(p_sigma, m_sigma, dx/mm, momentumSigma, flagTgtType);

    }
    */
    totalx += dx;
    localSigmaPos += momentumSigma*dx/momentumSigma.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomSigma ( " << momentumSigma.x() << ", "  
	      << momentumSigma.y() << ", " << momentumSigma.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localSigmaPos.x() << ", "  
	      << localSigmaPos.y() << ", " << localSigmaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumSigma, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;
    flagScattering=false;

    flagDecay = decayCheck(ctau, p_sigma, m_sigma, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);

    if (flagTgtType == 0 ) {
      react_rate = react_rate_p;
      dxInH += dx;

      flagScattering  = scatteringCheck(react_rate, dx/mm);
      flagLambdaPConv  = scatteringCheck(react_rate_lp_conv, dx/mm);
      flagSigmaPConv  = scatteringCheck(react_rate_sp_conv, dx/mm);

      p_sigma = calcEnergyDeposit(p_sigma, m_sigma, dx/mm, localSigmaPos, flagTgtType);
      E_sigma = sqrt(p_sigma*p_sigma + m_sigma*m_sigma);
    }
    /*
    G4cout << "Length : " << totalx << ", p = " << p_sigma*1000. 
	   << " MeV/c, Ekin = "	   << (E_sigma-m_sigma)*1000. 
	   << " MeV" << G4endl;
    */
    if (fabs(p_sigma)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }
    if (flagDecay || flagScattering || flagLambdaPConv || flagSigmaPConv)
      break;

    nIteration++;
    if (nIteration>5000)
      return;
  }

  momentumSigma *= p_sigma/momentumSigma.mag();

  if (flagDecay) {
    //std::cout << "decay" << std::endl;
    G4ThreeVector localDecayPos=localSigmaPos;
    G4ThreeVector globalDecayPos = geomMan.Local2GlobalPos(TgtId, localDecayPos);
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    if (1) {
      //neutron scattering
      Kinema3Resonance SigmaDecay;
      SigmaDecay = Kinema3Resonance(mass_sigma, 
				    0.0,
				    0.0 , neutron->GetPDGMass()/GeV,
				    piMinus->GetPDGMass()/GeV,
				    neutron->GetPDGMass()/GeV,
				    0.0, p_sigma, 0.0);


      /* Neutron */
      double Energy_decayN = SigmaDecay.GetEnergy(4);
      double momentum_decayN = SigmaDecay.GetMomentum(4);
      SigmaDecay.GetMomentum(4,mom);
      G4ThreeVector momentumDecayN(mom[1], mom[2], mom[0]);
      //double ThetaDecayN = SigmaDecay.GetTheta(4);
      //double PhiDecayN = SigmaDecay.GetPhi(4);
      
      /* pi- */
      double Energy_decayPi = SigmaDecay.GetEnergy(5);
      double momentum_decayPi = SigmaDecay.GetMomentum(5);
      SigmaDecay.GetMomentum(5,mom);
      G4ThreeVector momentumDecayPi( mom[1], mom[2], mom[0]);
      //double ThetaDecayPi = SigmaDecay.GetTheta(5);
      //double PhiDecayPi = SigmaDecay.GetPhi(5);
      
      // Neutron momentum at beam frame
      momentumDecayN.rotateY(ThetaSig*deg);
      momentumDecayN.rotateZ(PhiSig*deg);

      // pi- momentum at beam frame
      momentumDecayPi.rotateY(ThetaSig*deg);
      momentumDecayPi.rotateZ(PhiSig*deg);

      /*
      std::cout << "SigmaMom = ( " << momentumSigma.x() << ", " 
		<< momentumSigma.y() << ", " << momentumSigma.z()
		<< ")" << std::endl;
      std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
		<< momentumDecayN.y() << ", " << momentumDecayN.z()
		<< ")" << std::endl;
      std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
		<< momentumDecayPi.y() << ", " << momentumDecayPi.z()
		<< ")" << std::endl;

      std::cout << "Delta = ( " 
		<< momentumSigma.x()-momentumDecayN.x()-momentumDecayPi.x()
		<< ", " 
		<< momentumSigma.y()-momentumDecayN.y()-momentumDecayPi.y() 
		<< ", " 
		<< momentumSigma.z()-momentumDecayN.z()-momentumDecayPi.z()
		<< ")" << std::endl;
      */
      double ThetaDecayNatBeamFrame, PhiDecayNatBeamFrame; 
      calcThetaPhi(momentumDecayN, 
		   &ThetaDecayNatBeamFrame, 
		   &PhiDecayNatBeamFrame);

      double ThetaDecayPiatBeamFrame, PhiDecayPiatBeamFrame; 
      calcThetaPhi(momentumDecayPi, 
		   &ThetaDecayPiatBeamFrame, 
		   &PhiDecayPiatBeamFrame);
      /*
      std::cout << "ThetaDecayNatBeamFrame: " << ThetaDecayNatBeamFrame
		<< " PhiDecayNatBeamFrame: "  << PhiDecayNatBeamFrame
		<< std::endl;

      std::cout << "ThetaDecayPiatBeamFrame: " << ThetaDecayPiatBeamFrame
		<< " PhiDecayPiatBeamFrame: "  << PhiDecayPiatBeamFrame
		<< std::endl;
      */
      if (PhiDecayNatBeamFrame<0)
	return;
      if (PhiDecayPiatBeamFrame<0)
	return;


      /* n-p scatt */
      G4double neutron_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/

      double Ekin_n[30] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
		   110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
		   210., 220., 230., 240., 250., 260., 270., 280., 290., 300.};

      double np_cs_table[30] = { 951.827, 488.477, 312.839, 223.243, 170.507,
			   136.696, 113.728,  97.449,  85.517,  76.527,
			   69.593,  64.133,  59.756,  56.191,  53.244,
			   50.779,  48.691,  46.905,  45.364,  44.022,
			   42.844,  41.803,  40.876,  40.046,  39.298,
			   38.620,  38.002,  37.436,  36.915,  36.432};

      double Ekin = (sqrt(momentum_decayN*momentum_decayN+m_neutron*m_neutron)
		     -m_neutron)*1000. ;// MeV


      int index_Ekin=-1;
      double np_cs;
      for (int i=0; i<30; i++) {
	if (i==0) {
	  if (Ekin < Ekin_n[i]) {
	    index_Ekin = i;
	    np_cs = np_cs_table[i];
	    break;
	  } else if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	    index_Ekin = i;
	    np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	    break;
	  }
	} else if (i>=1 && i<=28) {
	  if (Ekin >= Ekin_n[i] && Ekin < Ekin_n[i+1]) {
	    index_Ekin = i;
	    np_cs = (np_cs_table[i+1]*(Ekin-Ekin_n[i]) + np_cs_table[i]*(Ekin_n[i+1]-Ekin))/(Ekin_n[i+1]-Ekin_n[i]);
	    break;
	  }
	} else if (i==29) {
	  if (Ekin >= Ekin_n[i]) {
	    index_Ekin = i;
	    np_cs = np_cs_table[i];
	    break;
	  }
	}
      }

      if (index_Ekin<0 || index_Ekin>=30) {
	fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_Ekin : %d", index_Ekin);
	exit(-1);
      }
      neutron_react_rate *= np_cs/30.;

      G4double dy = 0.5*mm; //mm
      G4double totaly=0.0;  //mm

      G4ThreeVector localNeutronPos = localDecayPos;
      int flagTgtType2 = -1;
      bool flagNeutronScattering=false;
      int nIte=0;
      while (1) {
	flagTgtType2 = -1;
	totaly += dy;
	localNeutronPos += momentumDecayN*dy/momentumDecayN.mag();

	flagNeutronScattering=false;

	flagTgtType2 = getTargetFlag(localNeutronPos);

	if (flagTgtType2 == 0 ) {
	  flagNeutronScattering  = scatteringCheck(neutron_react_rate, dy/mm);
	}
	
	if (flagNeutronScattering)
	  break;

	nIte++;
	if (nIte>1000)
	  break;
      }


      /* pi-p scatt */
      G4double pion_react_rate = 0.0001274; // 1/mm /* 30mb and LH2 target*/
      //G4double pion_react_rate = 0.0002548; // 1/mm /* 60mb and LH2 target*/

      double p_pi_table[8] = {0.1479, 0.1738, 0.1883, 0.2123, 0.2379, 0.2712, 0.2983, 0.3228};

      double pip_cs_table[8] = { 3.4, 5.12, 6.64, 9.8, 15.2, 19.4, 18.3, 14.8};
      double p = momentum_decayPi;
      int index_p=-1;
      double pip_cs;
      for (int i=0; i<8; i++) {
	if (i==0) {
	  if (p < p_pi_table[i]) {
	    index_p = i;
	    pip_cs = pip_cs_table[i];
	    break;
	  } else if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	    index_p = i;
	    pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	    break;
	  }
	} else if (i>=1 && i<=6) {
	  if (p >= p_pi_table[i] && p < p_pi_table[i+1]) {
	    index_p = i;
	    pip_cs = (pip_cs_table[i+1]*(p-p_pi_table[i]) + pip_cs_table[i]*(p_pi_table[i+1]-p))/(p_pi_table[i+1]-p_pi_table[i]);
	    break;
	  }
	} else if (i==7) {
	  if (p >= p_pi_table[i]) {
	    index_p = i;
	    pip_cs = pip_cs_table[i];
	    break;
	  }
	}
      }

      if (index_p<0 || index_p>=8) {
	fprintf(stderr, "PrimaryGeneratorAction::GeneratePiKSigmaScat2 invalid index_p : %d", index_p);
	exit(-1);
      }
      pion_react_rate *= pip_cs/30.;

      G4double dz = 0.5*mm; //mm
      G4double totalz=0.0;  //mm

      G4double ctau_pi=7804.5; /*mm*/
      G4double p_pi = momentumDecayPi.mag();
      G4double m_pi = piMinus->GetPDGMass()/GeV;

      G4ThreeVector localPionPos = localDecayPos;
      int flagTgtType3 = -1;

      G4bool flagPionScattering=false;
      G4bool flagPionDecay=false;

      nIte=0;
      while (1) {
	flagTgtType3 = -1;
	totalz += dz;
	localPionPos += momentumDecayPi*dz/momentumDecayPi.mag();

	flagPionScattering=false;
	flagPionDecay=false;

	flagTgtType3 = getTargetFlag(localPionPos);

	if (flagTgtType3 == 0) {
	  flagPionScattering  = scatteringCheck(pion_react_rate, dz/mm);
	}
	
	if (flagPionScattering)
	  break;

	flagPionDecay = decayCheck(ctau_pi, p_pi, m_pi, dz/mm );
	if (flagPionDecay)
	  break;

	nIte++;
	if (nIte>1000)
	  break;
      }

      if (flagNeutronScattering) {
	//std::cout << "NNscat" << std::endl;
	G4ThreeVector localNNScatPos=localNeutronPos;
	G4ThreeVector globalNNScatPos = geomMan.Local2GlobalPos(TgtId, localNNScatPos);
	int scatDistFlag = 0;
	G4ParticleDefinition* scatParticle;
	if ( flagTgtType2 == 0 ) {
	  // inside the LH2 target;
	  scatParticle = particleTable->FindParticle("proton");
	  scatDistFlag = 100+index_Ekin;
	} else {
	  return;
	}
	//std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	Kinema3Resonance NNScat;
	NNScat = Kinema3Resonance(neutron->GetPDGMass()/GeV, 
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , neutron->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  neutron->GetPDGMass()/GeV,
				  0.0, momentum_decayN, 0.0, scatDistFlag);

	double ThetaScatNeutronCM = NNScat.GetThetaCM(1);
	double PhiScatNeutronCM = NNScat.GetPhiCM(1);

	/* scattered Neutron */
	G4double Energy_scatN = NNScat.GetEnergy(4);
	//G4double momentum_scatN = NNScat.GetMomentum(4);
	NNScat.GetMomentum(4,mom);
	G4ThreeVector momentumScatN( mom[1], mom[2], mom[0]);

	momentumScatN.rotateY(ThetaDecayNatBeamFrame*deg);
	momentumScatN.rotateZ(PhiDecayNatBeamFrame*deg);
	
	/* scattered Proton or Deuteron*/
	G4double Energy_scatPart = NNScat.GetEnergy(5);
	//G4double momentum_scatPart = NNScat.GetMomentum(5);
	NNScat.GetMomentum(5,mom);
	G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	momentumScatPart.rotateY(ThetaDecayNatBeamFrame*deg);
	momentumScatPart.rotateZ(PhiDecayNatBeamFrame*deg);

	/*
	std::cout << "DecayN = ( " << momentumDecayN.x() << ", " 
		  << momentumDecayN.y() << ", " << momentumDecayN.z()
		  << ")" << std::endl;

	std::cout << "ScatN = ( " << momentumScatN.x() << ", " 
		  << momentumScatN.y() << ", " << momentumScatN.z()
		  << ")" << std::endl;

	std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
		  << momentumScatPart.y() << ", " << momentumScatPart.z()
		  << ")" << std::endl;

	std::cout << "Delta = ( " 
		  << momentumDecayN.x()-momentumScatN.x()-momentumScatPart.x()
		  << ", " 
		  << momentumDecayN.y()-momentumScatN.y()-momentumScatPart.y()
		  << ", " 
		  << momentumDecayN.z()-momentumScatN.z()-momentumScatPart.z()
		  << ")" << std::endl;
	*/

	particleGun->SetParticleDefinition(kPlus);
	G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
	particleGun->SetParticleMomentumDirection(gloMomKPlus);
	particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/* beam */
	particleGun->SetParticleDefinition(spiMinus);
	G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	particleGun->SetParticleMomentumDirection(gloMomBeam);
	particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/*sigma*/    
	particleGun->SetParticleDefinition(ssigma);
	G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	particleGun->SetParticleMomentumDirection(-gloMomSigma);
	// Energy_sig is original Energy without energy deposit
	particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/* decay pi- */
	particleGun->SetParticleDefinition(piMinus);
	G4ThreeVector gloMomDecayPi = geomMan.Local2GlobalDir(TgtId, momentumDecayPi);
	particleGun->SetParticleMomentumDirection(gloMomDecayPi);
	particleGun->SetParticleEnergy((Energy_decayPi - piMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/*scat neutron*/    
	particleGun->SetParticleDefinition(neutron);
	G4ThreeVector gloMomentumScatN = geomMan.Local2GlobalDir(TgtId, momentumScatN);
	particleGun->SetParticleMomentumDirection(gloMomentumScatN);
	particleGun->SetParticleEnergy((Energy_scatN - neutron->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalNNScatPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/*scatt proton*/    
	particleGun->SetParticleDefinition(scatParticle);
	G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	particleGun->SetParticleMomentumDirection(gloMomScatPart);
	particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalNNScatPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	anaMan_->SetPrimaryVertex(localVertexPos);
	anaMan_->SetHypBeamMomentum(momentumSigma);
	anaMan_->SetScatMesonMomentum(momentumKPlus);
	anaMan_->SetThetaMeson(ThetaK);
	anaMan_->SetPhiMeson(PhiK);
	anaMan_->SetThetaMesonCM(ThetaKCM);
	anaMan_->SetPhiMesonCM(PhiKCM);
	anaMan_->SetDecayPos(localDecayPos);
	anaMan_->SetNNScatPos(localNNScatPos);
	anaMan_->SetDecayFlag();
	anaMan_->SetNNScatFlag();
	anaMan_->SetNNScatTarget(flagTgtType2);	
	anaMan_->SetScatProtonMomentum(momentumScatPart);
	anaMan_->SetFlightLengthInTarget(dxInH, 0);
	anaMan_->SetDecayPiMomentum(momentumDecayPi);
	anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	anaMan_->SetBeamMomentum(beammom);
	anaMan_->SetThetaScatHypCM(ThetaScatNeutronCM);
	anaMan_->SetPhiScatHypCM(PhiScatNeutronCM);

	return;
      }

      if (flagPionScattering) {
	//std::cout << "piNscat" << std::endl;
	G4ThreeVector localPiNScatPos=localPionPos;
	G4ThreeVector globalPiNScatPos = geomMan.Local2GlobalPos(TgtId, localPiNScatPos);
	int scatDistFlag = 0;
	G4ParticleDefinition* scatParticle;
	if ( flagTgtType3 == 0 ) {
	  // inside the LH2 target;
	  scatParticle = particleTable->FindParticle("proton");
	  scatDistFlag = 0;
	} else {
	  return;
	}
	//std::cout << scatParticle->GetPDGMass()/GeV << std::endl; 
	Kinema3Resonance PiNScat;
	PiNScat = Kinema3Resonance(piMinus->GetPDGMass()/GeV, 
				  scatParticle->GetPDGMass()/GeV,
				  0.0 , piMinus->GetPDGMass()/GeV,
				  scatParticle->GetPDGMass()/GeV,
				  piMinus->GetPDGMass()/GeV,
				  0.0, momentum_decayPi, 0.0, scatDistFlag);

	/* scattered pion */
	G4double Energy_scatPi = PiNScat.GetEnergy(4);
	//G4double momentum_scatPi = PiNScat.GetMomentum(4);
	PiNScat.GetMomentum(4,mom);
	G4ThreeVector momentumScatPi( mom[1], mom[2], mom[0]);

	momentumScatPi.rotateY(ThetaDecayPiatBeamFrame*deg);
	momentumScatPi.rotateZ(PhiDecayPiatBeamFrame*deg);
	
	/* scattered Proton or Deuteron*/
	G4double Energy_scatPart = PiNScat.GetEnergy(5);
	//G4double momentum_scatPart = PiNScat.GetMomentum(5);
	PiNScat.GetMomentum(5,mom);
	G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
	momentumScatPart.rotateY(ThetaDecayPiatBeamFrame*deg);
	momentumScatPart.rotateZ(PhiDecayPiatBeamFrame*deg);

	/*
	std::cout << "DecayPi = ( " << momentumDecayPi.x() << ", " 
		  << momentumDecayPi.y() << ", " << momentumDecayPi.z()
		  << ")" << std::endl;

	std::cout << "ScatPi = ( " << momentumScatPi.x() << ", " 
		  << momentumScatPi.y() << ", " << momentumScatPi.z()
		  << ")" << std::endl;

	std::cout << "ScatPart = ( " << momentumScatPart.x() << ", " 
		  << momentumScatPart.y() << ", " << momentumScatPart.z()
		  << ")" << std::endl;

	std::cout << "Delta = ( " 
		  << momentumDecayPi.x()-momentumScatPi.x()-momentumScatPart.x()
		  << ", " 
		  << momentumDecayPi.y()-momentumScatPi.y()-momentumScatPart.y()
		  << ", " 
		  << momentumDecayPi.z()-momentumScatPi.z()-momentumScatPart.z()
		  << ")" << std::endl;
	*/

	particleGun->SetParticleDefinition(kPlus);
	G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
	particleGun->SetParticleMomentumDirection(gloMomKPlus);
	particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/* beam */
	particleGun->SetParticleDefinition(spiMinus);
	G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
	particleGun->SetParticleMomentumDirection(gloMomBeam);
	particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalVertexPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/*sigma*/    
	particleGun->SetParticleDefinition(ssigma);
	G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
	particleGun->SetParticleMomentumDirection(-gloMomSigma);
	// Energy_sig is original Energy without energy deposit
	particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);
    
	/* decay neutron */
	particleGun->SetParticleDefinition(neutron);
	G4ThreeVector gloMomDecayN = geomMan.Local2GlobalDir(TgtId, momentumDecayN);
	particleGun->SetParticleMomentumDirection(gloMomDecayN);
	particleGun->SetParticleEnergy((Energy_decayN - neutron->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalDecayPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/*scat pi-*/    
	particleGun->SetParticleDefinition(piMinus);
	G4ThreeVector gloMomentumScatPi = geomMan.Local2GlobalDir(TgtId, momentumScatPi);
	particleGun->SetParticleMomentumDirection(gloMomentumScatPi);
	particleGun->SetParticleEnergy((Energy_scatPi - piMinus->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalPiNScatPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	/*scatt proton*/    
	particleGun->SetParticleDefinition(scatParticle);
	G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
	particleGun->SetParticleMomentumDirection(gloMomScatPart);
	particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
	particleGun->SetParticlePosition(globalPiNScatPos);
	particleGun->GeneratePrimaryVertex(anEvent);

	anaMan_->SetPrimaryVertex(localVertexPos);
	anaMan_->SetHypBeamMomentum(momentumSigma);
	anaMan_->SetScatMesonMomentum(momentumKPlus);
	anaMan_->SetThetaMeson(ThetaK);
	anaMan_->SetPhiMeson(PhiK);
	anaMan_->SetThetaMesonCM(ThetaKCM);
	anaMan_->SetPhiMesonCM(PhiKCM);
	anaMan_->SetDecayPos(localDecayPos);
	anaMan_->SetPiNScatPos(localPiNScatPos);
	anaMan_->SetDecayFlag();
	anaMan_->SetPiNScatFlag();
	anaMan_->SetPiNScatTarget(flagTgtType3);	
	anaMan_->SetScatProtonMomentum(momentumScatPart);
	anaMan_->SetFlightLengthInTarget(dxInH, 0);
	//anaMan_->SetDecayPiMomentum(momentumDecayPi);
	anaMan_->SetDecayPiMomentum(momentumScatPi); // momentum after scattering
	anaMan_->SetDecayNucleonMomentum(momentumDecayN);
	anaMan_->SetBeamMomentum(beammom);

	return;
      }
    }

    static int nDecay = 0;
    const int prescale = 100;
    //const int prescale = 1;
    nDecay++;

    if (nDecay%prescale != 0)
      return;

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetDecayPos(localSigmaPos);
    anaMan_->SetDecayEventFlag();
    anaMan_->SetDecayFlag();
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /* beam */
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*sigma*/    
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);
    
    /*decaying sigma*/    
    particleGun->SetParticleDefinition(usigma);
    particleGun->SetParticleMomentumDirection(gloMomSigma);
    // E_sigma is energy after energy deposit
    particleGun->SetParticleEnergy((E_sigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalDecayPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  } else if (flagScattering) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localSigmaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);
    Kinema3Resonance SigmaScat;
    G4ParticleDefinition* scatParticle;
    int scatDistFlag = 0;
    if ( flagTgtType == 0 ) {
      // inside the LH2 target;
      scatParticle = particleTable->FindParticle("proton");
      scatDistFlag = 0;
    } else {
      return;
    }

    SigmaScat = Kinema3Resonance(mass_sigma, 
				 scatParticle->GetPDGMass()/GeV,
				 0.0 , mass_sigma,
				 scatParticle->GetPDGMass()/GeV,
				 mass_sigma, 0.0, p_sigma, 0.0,
				 scatDistFlag);
    
    
    /* scattered Sigma */
    G4double Energy_scatSigma = SigmaScat.GetEnergy(4);
    //G4double momentum_scatSigma = SigmaScat.GetMomentum(4);
    SigmaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatSigma( mom[1], mom[2], mom[0]);
    momentumScatSigma.rotateY(ThetaSig*deg);
    momentumScatSigma.rotateZ(PhiSig*deg);
    
    double ThetaScatSigCM = SigmaScat.GetThetaCM(1);
    double PhiScatSigCM = SigmaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = SigmaScat.GetEnergy(5);
    //G4double momentum_scatPart = SigmaScat.GetMomentum(5);
    SigmaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaSig*deg);
    momentumScatPart.rotateZ(PhiSig*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatHypMomentum(momentumScatSigma);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatSigCM);
    anaMan_->SetPhiScatHypCM(PhiScatSigCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag();
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*sigma*/    
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt sigma*/    
    particleGun->SetParticleDefinition(sigma);
    G4ThreeVector gloMomScatSigma = geomMan.Local2GlobalDir(TgtId, momentumScatSigma);
    particleGun->SetParticleMomentumDirection(gloMomScatSigma);
    particleGun->SetParticleEnergy((Energy_scatSigma - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

  }  else if (flagLambdaPConv || flagSigmaPConv) {
    //std::cout << "scat" << std::endl;
    G4ThreeVector localScatPos=localSigmaPos;
    G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);

    Kinema3Resonance SigmaScat;
    G4ParticleDefinition* scatParticle= particleTable->FindParticle("neutron");
    int scatDistFlag = 0;
    G4ParticleDefinition* hyperon;
    int reactMode=-100;
    if (flagLambdaPConv) {
      hyperon = particleTable->FindParticle("lambda");
      reactMode=2;
    } else { 
      hyperon = particleTable->FindParticle("sigma0");
      reactMode=3;
    }
    double mass_proton = 0.93827200;

    SigmaScat = Kinema3Resonance(mass_sigma, 
				 mass_proton,
				 0.0 , hyperon->GetPDGMass()/GeV,
				 scatParticle->GetPDGMass()/GeV,
				 hyperon->GetPDGMass()/GeV, 0.0, p_sigma, 0.0,
				 scatDistFlag);

    
    /* scattered Hyperon */
    G4double Energy_scatSigma = SigmaScat.GetEnergy(4);
    //G4double momentum_scatSigma = SigmaScat.GetMomentum(4);
    SigmaScat.GetMomentum(4,mom);
    G4ThreeVector momentumScatSigma( mom[1], mom[2], mom[0]);
    momentumScatSigma.rotateY(ThetaSig*deg);
    momentumScatSigma.rotateZ(PhiSig*deg);
    
    double ThetaScatSigCM = SigmaScat.GetThetaCM(1);
    double PhiScatSigCM = SigmaScat.GetPhiCM(1);
    
    /*
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
    */
    /* scattered Proton or Deuteron*/
    G4double Energy_scatPart = SigmaScat.GetEnergy(5);
    //G4double momentum_scatPart = SigmaScat.GetMomentum(5);
    SigmaScat.GetMomentum(5,mom);
    G4ThreeVector momentumScatPart(mom[1], mom[2], mom[0]);
    momentumScatPart.rotateY(ThetaSig*deg);
    momentumScatPart.rotateZ(PhiSig*deg);
    /*
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */
    //G4cout << "ThetaPi : " << ThetaPi  << ", "
    //<< "PhiPi : " << PhiPi  << G4endl;
    
    
    //if (1) {
    /*
      G4cout << "Flight length = " << flength/cm << G4endl;
      G4cout << "Sigma (" << momentumSigma.x() << ", "
      << momentumSigma.y() << ", " << momentumSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatSigma (" << momentumScatSigma.x() << ", "
      << momentumScatSigma.y() << ", " << momentumScatSigma.z() << ") "
      << G4endl;
      
      G4cout << "ScatProton (" << momentumScatProton.x() << ", "
      << momentumScatProton.y() << ", " << momentumScatProton.z() << ") "
      << G4endl;
    */

    anaMan_->SetPrimaryVertex(localVertexPos);
    anaMan_->SetHypBeamMomentum(momentumSigma);
    anaMan_->SetScatHypMomentum(momentumScatSigma);
    anaMan_->SetScatProtonMomentum(momentumScatPart);
    anaMan_->SetScatMesonMomentum(momentumKPlus);
    anaMan_->SetThetaMeson(ThetaK);
    anaMan_->SetPhiMeson(PhiK);
    anaMan_->SetThetaMesonCM(ThetaKCM);
    anaMan_->SetPhiMesonCM(PhiKCM);
    anaMan_->SetThetaScatHypCM(ThetaScatSigCM);
    anaMan_->SetPhiScatHypCM(PhiScatSigCM);
    anaMan_->SetScatPos(localScatPos);
    anaMan_->SetScatFlag(reactMode);
    anaMan_->SetScatTarget(flagTgtType);
    anaMan_->SetFlightLengthInTarget(dxInH, 0);
    anaMan_->SetBeamMomentum(beammom);

    particleGun->SetParticleDefinition(kPlus);
    G4ThreeVector gloMomKPlus = geomMan.Local2GlobalDir(TgtId, momentumKPlus);
    particleGun->SetParticleMomentumDirection(gloMomKPlus);
    particleGun->SetParticleEnergy((Energy_k - kPlus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /* beam */
    particleGun->SetParticleDefinition(spiMinus);
    G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
    particleGun->SetParticleMomentumDirection(gloMomBeam);
    particleGun->SetParticleEnergy((Energy_beam - spiMinus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalVertexPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*sigma*/    
    particleGun->SetParticleDefinition(ssigma);
    G4ThreeVector gloMomSigma = geomMan.Local2GlobalDir(TgtId, momentumSigma);
    particleGun->SetParticleMomentumDirection(-gloMomSigma);
    // Energy_sig is original Energy without energy deposit
    particleGun->SetParticleEnergy((Energy_sig - sigma->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt hyperon*/    
    particleGun->SetParticleDefinition(hyperon);
    G4ThreeVector gloMomScatSigma = geomMan.Local2GlobalDir(TgtId, momentumScatSigma);
    particleGun->SetParticleMomentumDirection(gloMomScatSigma);
    particleGun->SetParticleEnergy((Energy_scatSigma - hyperon->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);

    /*scatt proton*/    
    particleGun->SetParticleDefinition(scatParticle);
    G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
    particleGun->SetParticleMomentumDirection(gloMomScatPart);
    particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(globalScatPos);
    particleGun->GeneratePrimaryVertex(anEvent);
  }
}


void CFTPrimaryGeneratorAction::GeneratePiMinusBeam(G4Event* anEvent)
{
  double beammom;
  double Energy_beam;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* piMinus;

  piMinus = particleTable->FindParticle("pi-");
  beammom = 1.32;

  static FILE *fp = NULL;
  if (!fp) {
    fp = fopen("param/beamPosAndDir.txt","r");
    if (!fp) {
      fprintf(stderr, "cannot open param/beamPosAndDir.txt\n");
    }
  }

  char buf[100];
  if (fgets(buf, sizeof(buf), fp) == NULL) {
    fp = NULL;
    return;
  }
  double xtgt, ytgt, u0, v0;
  sscanf(buf, "%lf %lf %lf %lf", &xtgt, &ytgt, &u0, &v0);

  /* pi- beam */
  Energy_beam = sqrt(pow(piMinus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  //G4ThreeVector momentumBeam( 0., 0., beammom);
  G4ThreeVector momentumBeam( u0, v0, 1.0);

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;

  //primary_vertex_z = -1200.*mm;
  //primary_vertex_z = -150.*mm;
  primary_vertex_z = -210.*mm; // original
  primary_vertex_x = xtgt + u0*primary_vertex_z;
  primary_vertex_y = ytgt + v0*primary_vertex_z;
  //primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  //primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;


  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  /* beam */
  particleGun->SetParticleDefinition(piMinus);
  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
  particleGun->SetParticleMomentumDirection(gloMomBeam);
  particleGun->SetParticleEnergy((Energy_beam - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalVertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetBeamMomentum(beammom);

}

void CFTPrimaryGeneratorAction::GeneratePiPlusBeam(G4Event* anEvent)
{
  double beammom;
  double Energy_beam;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* piPlus;

  piPlus = particleTable->FindParticle("pi+");
  beammom = 1.40;


  static FILE *fp = NULL;
  if (!fp) {
    fp = fopen("param/beamPosAndDir.txt","r");
    if (!fp) {
      fprintf(stderr, "cannot open param/beamPosAndDir.txt\n");
    }
  }

  char buf[100];
  if (fgets(buf, sizeof(buf), fp) == NULL) {
    fp = NULL;
    return;
  }
  double xtgt, ytgt, u0, v0;
  sscanf(buf, "%lf %lf %lf %lf", &xtgt, &ytgt, &u0, &v0);

  /* pi+ beam */
  Energy_beam = sqrt(pow(piPlus->GetPDGMass()/GeV, 2.0) + pow(beammom, 2.0));
  //G4ThreeVector momentumBeam( 0., 0., beammom);
  G4ThreeVector momentumBeam( u0, v0, 1.0);

  G4double primary_vertex_x, primary_vertex_y, primary_vertex_z;

  //primary_vertex_z = -1000.*mm;
  //primary_vertex_z = -150.*mm;
  primary_vertex_z = -210.*mm;
  primary_vertex_x = xtgt + u0*primary_vertex_z;
  primary_vertex_y = ytgt + v0*primary_vertex_z;
  //primary_vertex_x = (G4double)RandGauss::shoot(-5., 5.)*mm;
  //primary_vertex_y = (G4double)RandGauss::shoot(5., 5.)*mm;

  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  //Rotate to global coordinate
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalVertexPos = geomMan.Local2GlobalPos(TgtId, localVertexPos);
  
  /* beam */
  particleGun->SetParticleDefinition(piPlus);
  G4ThreeVector gloMomBeam = geomMan.Local2GlobalDir(TgtId, momentumBeam);
  particleGun->SetParticleMomentumDirection(gloMomBeam);
  particleGun->SetParticleEnergy((Energy_beam - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalVertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetBeamMomentum(beammom);

}


void CFTPrimaryGeneratorAction::EfficiencyStudyForScatPart(G4Event* anEvent, int reactionMode)
{

  static int num=0;
  const int Nstudy = 1000;

  char buf[300];
  static double scatPosX, scatPosY, scatPosZ;
  static double hypBeamMomX, hypBeamMomY, hypBeamMomZ;
  static double thetaScatPart, momScatPart;
  static int flagTgtType;

  static bool flagFileEnd = false;


  if (flagFileEnd)
    return;

  if (num==0) {
    if (fgets(buf, sizeof(buf), fpIn) == NULL) {
      flagFileEnd = true;
      return;
    }
    sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf %d ",
	   &scatPosX, &scatPosY, &scatPosZ,
	   &hypBeamMomX, &hypBeamMomY, &hypBeamMomZ,
	   &thetaScatPart, &momScatPart, &flagTgtType);
    printf("scatPos = (%lf, %lf, %lf),  hypBeamMom=(%lf, %lf, %lf), theta=%lf, momScat=%lf,  flagTgt=%d\n",
	   scatPosX, scatPosY, scatPosZ,
	   hypBeamMomX, hypBeamMomY, hypBeamMomZ,
	   thetaScatPart, momScatPart, flagTgtType);
  }


  num++;
  if (num>=Nstudy)
    num=0;

  if (momScatPart<0)
    momScatPart = 0.;

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ThreeVector hypBeamVec(hypBeamMomX, hypBeamMomY, hypBeamMomZ);
  double ThetaSig, PhiSig;
  calcThetaPhi(hypBeamVec, &ThetaSig, &PhiSig);

  G4ThreeVector localScatPos(scatPosX, scatPosY, scatPosZ);
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);

  double phiScatPart = RandFlat::shoot(0., 360.);

  G4ParticleDefinition* scatParticle;
  if ( flagTgtType == 0 || flagTgtType == -1) {
    // inside the LH2 target;
    scatParticle = particleTable->FindParticle("proton");
    
  } else if (flagTgtType == 1) {
    // inside the LD2 target;
    scatParticle = particleTable->FindParticle("deuteron");      
  } else {
    return;
  }

  G4ThreeVector 
    momentumScatPart( momScatPart*sin(thetaScatPart*Deg2Rad)*cos(phiScatPart*Deg2Rad), 
		      momScatPart*sin(thetaScatPart*Deg2Rad)*sin(phiScatPart*Deg2Rad), 
		      momScatPart*cos(thetaScatPart*Deg2Rad)); 
  momentumScatPart.rotateY(ThetaSig*deg);
  momentumScatPart.rotateZ(PhiSig*deg);
  double Energy_scatPart =sqrt(momScatPart*momScatPart+
			       (scatParticle->GetPDGMass()/GeV)*(scatParticle->GetPDGMass()/GeV));

  anaMan_->SetScatProtonMomentum(momentumScatPart);
  anaMan_->SetScatPos(localScatPos);
  anaMan_->SetScatFlag();
  anaMan_->SetScatTarget(flagTgtType);
  
  /*scatt proton*/    
  
  particleGun->SetParticleDefinition(scatParticle);
  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
  particleGun->SetParticleMomentumDirection(gloMomScatPart);
  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalScatPos);
  particleGun->GeneratePrimaryVertex(anEvent);


  /* scattered hyperon */
  G4ParticleDefinition* hyperon = 0;
  if ( reactionMode == 95) {
    hyperon = particleTable->FindParticle("sigma-");
  } else if ( reactionMode == 99) {
    hyperon = particleTable->FindParticle("sigma1+");
  } else if ( reactionMode == 100) {
    hyperon = particleTable->FindParticle("sigma2+");
  }

  G4ThreeVector momentumHyperon = hypBeamVec - momentumScatPart;
  double momHyperon = momentumHyperon.mag();
  double massHyperon = hyperon->GetPDGMass()/GeV;
  double Energy_hyperon =sqrt(momHyperon*momHyperon+massHyperon*massHyperon);

  
  
  particleGun->SetParticleDefinition(hyperon);
  G4ThreeVector gloMomHyperon = geomMan.Local2GlobalDir(TgtId, momentumHyperon);
  particleGun->SetParticleMomentumDirection(gloMomHyperon);
  particleGun->SetParticleEnergy((Energy_hyperon - massHyperon)*GeV);
  particleGun->SetParticlePosition(globalScatPos);
  particleGun->GeneratePrimaryVertex(anEvent);


}

void CFTPrimaryGeneratorAction::EfficiencyStudyForLambda(G4Event* anEvent)
{

  static int num=0;
  const int Nstudy = 1000;

  char buf[300];
  static double scatPosX, scatPosY, scatPosZ;
  static double hypBeamMomX, hypBeamMomY, hypBeamMomZ;
  static double thetaScatPart, momScatPart, momCalScatPart;
  static int flagTgtType;

  if (num==0) {
    if (fgets(buf, sizeof(buf), fpIn) == NULL)
      return;

    sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d ",
	   &scatPosX, &scatPosY, &scatPosZ,
	   &hypBeamMomX, &hypBeamMomY, &hypBeamMomZ,
	   &thetaScatPart, &momScatPart, &momCalScatPart, &flagTgtType);
    printf("scatPos = (%lf, %lf, %lf),  hypBeamMom=(%lf, %lf, %lf), theta=%lf, momScat=%lf, momCalScat=%lf, flagTgt=%d\n",
	   scatPosX, scatPosY, scatPosZ,
	   hypBeamMomX, hypBeamMomY, hypBeamMomZ,
	   thetaScatPart, momScatPart, momCalScatPart, flagTgtType);
  }
  num++;
  if (num>=Nstudy)
    num=0;


  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ThreeVector hypBeamVec(hypBeamMomX, hypBeamMomY, hypBeamMomZ);
  double ThetaSig, PhiSig;
  calcThetaPhi(hypBeamVec, &ThetaSig, &PhiSig);

  G4ThreeVector localScatPos(scatPosX, scatPosY, scatPosZ);
  G4int TgtId = geomMan.GetDetectorId("Target");
  G4ThreeVector globalScatPos = geomMan.Local2GlobalPos(TgtId, localScatPos);

  double phiScatPart = RandFlat::shoot(0., 360.);

  G4ParticleDefinition* scatParticle;
  scatParticle = particleTable->FindParticle("lambda");

  G4ThreeVector 
    momentumScatPart( momScatPart*sin(thetaScatPart*Deg2Rad)*cos(phiScatPart*Deg2Rad), 
		      momScatPart*sin(thetaScatPart*Deg2Rad)*sin(phiScatPart*Deg2Rad), 
		      momScatPart*cos(thetaScatPart*Deg2Rad)); 
  momentumScatPart.rotateY(ThetaSig*deg);
  momentumScatPart.rotateZ(PhiSig*deg);
  double Energy_scatPart =sqrt(momScatPart*momScatPart+
			       (scatParticle->GetPDGMass()/GeV)*(scatParticle->GetPDGMass()/GeV));

  anaMan_->SetScatProtonMomentum(momentumScatPart);
  anaMan_->SetScatPos(localScatPos);
  anaMan_->SetScatFlag();
  anaMan_->SetScatTarget(flagTgtType);
  
  /*scatt lambda*/    
  
  particleGun->SetParticleDefinition(scatParticle);
  G4ThreeVector gloMomScatPart = geomMan.Local2GlobalDir(TgtId, momentumScatPart);
  particleGun->SetParticleMomentumDirection(gloMomScatPart);
  particleGun->SetParticleEnergy((Energy_scatPart - scatParticle->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(globalScatPos);
  particleGun->GeneratePrimaryVertex(anEvent);


}

void CFTPrimaryGeneratorAction::GenerateSigmaBeam()
{

  char buf[300];
  if (fgets(buf, sizeof(buf), fpSigBeam) == NULL)
    return;
  
  double primary_vertex_x, primary_vertex_y, primary_vertex_z;
  double mom_x, mom_y, mom_z;

  sscanf(buf, "%lf %lf %lf %lf %lf %lf",
	 &primary_vertex_x, &primary_vertex_y, &primary_vertex_z,
	 &mom_x, &mom_y, &mom_z);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* sigma;
  if (ReactionMode_ == 96)
    sigma = particleTable->FindParticle("sigma-");
  else if (ReactionMode_ == 102)
    sigma = particleTable->FindParticle("sigma+");
  else {
    fprintf(stderr, "CFTPrimaryGeneratorAction::GenerateSigmaBeam::Invalid reaction mode %d\n", ReactionMode_);
    exit(-1);
  }

  G4ThreeVector localVertexPos(primary_vertex_x, primary_vertex_y, primary_vertex_z);
  G4ThreeVector momentumSigma( mom_x, mom_y, mom_z);
  
  /* Sigma-p scatt */
  //G4double flength=20.*mm;
  G4double ctau=0;
  if (ReactionMode_ == 96)
    ctau=44.34; /*mm*/
  else if (ReactionMode_ == 102)
    ctau=24.04; /*mm*/

  G4double p_sigma = momentumSigma.mag();
  G4double m_sigma = sigma->GetPDGMass()/GeV;

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  double dxInH=0.;
  double dxInD=0.;

  G4ThreeVector localSigmaPos = localVertexPos;
  G4bool flagDecay=false;
  int flagTgtType = -1;

  int nIteration=0;
  while (1) {
    flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);

    totalx += dx;
    localSigmaPos += momentumSigma*dx/momentumSigma.mag();

#if 0
    std::cout << "------" << std::endl;
    std::cout << "MomSigma ( " << momentumSigma.x() << ", "  
	      << momentumSigma.y() << ", " << momentumSigma.z()
	      << ")" << std::endl;
    std::cout << "LocalPos ( " << localSigmaPos.x() << ", "  
	      << localSigmaPos.y() << ", " << localSigmaPos.z()
	      << ")" << std::endl;
    {
      double theta1, phi1;
      calcThetaPhi(momentumSigma, &theta1, &phi1);
      std::cout << "Theta " << ThetaSig << "--> "  << theta1
		<< ", Phi " << PhiSig << "--> "  << phi1
		<<   std::endl;
    }
#endif
    flagDecay=false;

    // tmp
    flagDecay = decayCheck(ctau, p_sigma, m_sigma, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localSigmaPos);


    if (flagTgtType == 0)
      dxInH += dx;
    else if (flagTgtType == 1)
      dxInD += dx;

    if (fabs(p_sigma)<0.000001) {
      flagDecay = true;
    }
    if (flagDecay)
      break;

    nIteration++;
    if (nIteration>5000) {
      G4cout << "Momentum = " << p_sigma << G4endl;
      break;
    }
  }
  anaMan_->SetPrimaryVertex(localVertexPos);
  anaMan_->SetFlightLengthInTarget(dxInH, dxInD);
  anaMan_->SetHypBeamMomentum(momentumSigma);
}


void CFTPrimaryGeneratorAction::calcThetaPhi(G4ThreeVector vec, double *theta, double *phi)
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


int CFTPrimaryGeneratorAction::getTargetFlag(G4ThreeVector pos)
{
  int flagTgtType = -1;

  double r = sqrt((pos.x()/mm)*(pos.x()/mm)+(pos.y()/mm)*(pos.y()/mm));
  if ( r <= LH2TgtR/2. 
       && fabs(pos.z()/mm) <= Target_Length/2. ) {
    // inside the LH2 target;
    flagTgtType = 0;
    return flagTgtType;
  }

  return flagTgtType;	

}

G4bool CFTPrimaryGeneratorAction::decayCheck(G4double ctau, G4double momentum, G4double mass, G4double dx)
{
  G4double lambda=ctau*momentum/mass;
  G4double value;

  value = exp(-dx/lambda);

  G4double x;
  x=(G4double)RandFlat::shoot();
  if (x > value)
    return true;
  else 
    return false;  

}

G4bool CFTPrimaryGeneratorAction::scatteringCheck(G4double rate, G4double dx)
{
  G4double value;

  value = rate*dx;

  G4double x;
  x=(G4double)RandFlat::shoot();
  if (x < value)
    return true;
  else 
    return false;  

}


G4double CFTPrimaryGeneratorAction::calcEnergyDeposit(G4double momentum, G4double mass, G4double dx, G4ThreeVector pos, G4int flagTgtType)
{

  G4double E0 = sqrt(momentum*momentum+mass*mass) * 1000.; //GeV --> MeV
  G4double beta = momentum/sqrt(momentum*momentum+mass*mass);
  
  G4double dEdx = calc_dE_dx(beta, flagTgtType);    // MeV/cm
  if (dEdx < 0.) {
    return -1.;
  }
  G4double dE = dEdx * dx / 10.; // dx (mm), MeV
  G4double E = E0 - dE;          // MeV
  if (E-mass*1000.<0.005) {
    printf("#particle stops in material at (%f, %f, %f)\n",pos.x(), pos.y(), pos.z());
    return 0.0;
  }

  G4double p = sqrt((E/1000.)*(E/1000.) -mass*mass); //GeV/c
  
  return p;
}

G4double CFTPrimaryGeneratorAction::calc_dE_dx(double beta, int flagTgtType)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;

  int z = 1;

  double rho;   /*g/cm^3*/
  double I;     /*eV*/
  double Z_A;

  if ( flagTgtType == 0 ) {
    // inside the LH2 target;
    rho = 0.0708; /* g/cm^3 */
    I   = 21.8;   /* eV */
    Z_A = 0.99212;
  } else if (flagTgtType == 1) {
    // inside the LD2 target;
    rho = 0.169; /* g/cm^3 */
    I   = 21.8;  /* eV */
    Z_A = 0.49652;

  } else {
    //fprintf(stderr, "Out of the target : (%f, %f, %f)\n", pos.x(), pos.y(), pos.z());
    return 0.0;
  }

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}


