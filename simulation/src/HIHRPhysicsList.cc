//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file HIHRPhysicsList.cc
/// \brief Implementation of the HIHRPhysicsList class

#include "HIHRPhysicsList.hh"
#include "Transportation.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsXS.hh"

#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsShielding.hh"

#include "G4IonElasticPhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonPhysicsPHP.hh"
#include "G4IonINCLXXPhysics.hh"

//#include "GammaNuclearPhysics.hh"


#include "G4ProcessManager.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HIHRPhysicsList::HIHRPhysicsList() 
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // EM physics
  RegisterPhysics(new G4EmStandardPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  G4int verb = 0;  

  // Hadron Elastic scattering
  //
  // This
  RegisterPhysics( new G4HadronElasticPhysicsHP(verb));
  ///RegisterPhysics( new G4HadronElasticPhysicsXS(verb));  

  // Hadron Inelastic physics
  // This
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));

  //RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC_AllHP(verb));
  ////RegisterPhysics( new G4HadronPhysicsQGSP_BIC(verb));  
  ////RegisterPhysics( new G4HadronInelasticQBBC(verb));
  ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));
  ////RegisterPhysics( new G4HadronPhysicsShielding(verb));
    
  // Ion Elastic scattering
  //
  RegisterPhysics( new G4IonElasticPhysics(verb));
  
  // Ion Inelastic physics
  // This
  RegisterPhysics( new G4IonPhysicsXS(verb));
  ////RegisterPhysics( new G4IonPhysicsPHP(verb));
  ////RegisterPhysics( new G4IonQMDPhysics(verb));
  ////RegisterPhysics( new G4IonINCLXXPhysics(verb));

  // Gamma physics
  //
  //RegisterPhysics( new GammaNuclearPhysics("gamma"));

  //AddTransportationSks();
  ConstructSpecialParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HIHRPhysicsList::~HIHRPhysicsList()
{ 
}

///////Transportation
void HIHRPhysicsList::AddTransportationSks()
{
  Transportation* theTransportationProcess= new Transportation();

  // loop over all particles in G4ParticleTable
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if ( !particle->IsShortLived() ) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        G4Exception("K18PhysicsList::AddTransportation : no process manager!","", FatalException, "");
      } 
      else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager->AddProcess(theTransportationProcess);
        pmanager->SetProcessOrderingToFirst(theTransportationProcess,
                                            idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess,
                                             idxPostStep);
      }
    }
    else {
      // shortlived particle case
    }
  }
}

void HIHRPhysicsList::ConstructSpecialParticle()
{
  G4DecayTable* decayTable;
  G4VDecayChannel* mode;
  G4ParticleDefinition* particle;


  // skaon- non-decay kaon-
  particle = new G4ParticleDefinition(
				      "skaon-",    0.493677*GeV,   5.315e-14*MeV,    -1.*eplus,
				      0,              -1,             0,
				      1,              -1,             0,
				      "meson",               0,             0,        -321,
				      true,       0,          NULL);

  // skaon+ non-decay kaon+
  particle = new G4ParticleDefinition(
				      "skaon+",    0.493677*GeV,  5.315e-14*MeV,    +1.*eplus,
				      0,              -1,             0,
				      1,              +1,             0,
				      "meson",               0,             0,         321,
				      true,       0,          NULL);


  // spi+ non-decay kaon-
  particle = new G4ParticleDefinition(
				      "spi+",    0.1395700*GeV,  2.5284e-14*MeV,    +1.*eplus,
				      0,              -1,             0,
				      2,              +2,            -1,
				      "meson",               0,             0,         211,
				      true,       0,          NULL);

  // spi- non-decay kaon-
  particle = new G4ParticleDefinition(
				      "spi-",    0.1395700*GeV, 2.5284e-14*MeV,    -1.*eplus,
				      0,              -1,             0,
				      2,              -2,            -1,
				      "meson",               0,             0,        -211,
				      true,               0,          NULL);


  // usigma- life time is 0
  particle = new G4ParticleDefinition(
				      "usigma-",    1.19744*GeV,  4.45e-12*MeV,    -1*eplus,
				      1,              +1,             0,
				      2,              -2,             0,
				      "baryon",               0,            +1,        3112,
				      false,       0.0*ns,          NULL);

  //create Decay Table 
  decayTable = new G4DecayTable();
  // create decay channels
  // sigma- -> neutron + pi-
  mode = new G4PhaseSpaceDecayChannel("usigma-",1.00,2,"neutron","pi-");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);



  // usigma+ life time is 0
  particle = new G4ParticleDefinition(
           "usigma+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0,          NULL);

 //create Decay Table
  decayTable = new G4DecayTable();

  // create decay channels
  G4VDecayChannel** modetable = new G4VDecayChannel*[2];
  // sigma+ -> proton + pi0
  modetable[0] = new G4PhaseSpaceDecayChannel("usigma+",0.516,2,"proton","pi0");
  // sigma+ -> neutron + pi+
  modetable[1] = new G4PhaseSpaceDecayChannel("usigma+",0.483,2,"neutron","pi+");

  for (G4int index=0; index <2; index++ ) decayTable->Insert(modetable[index]);
  delete [] modetable;
  particle->SetDecayTable(decayTable);
  
  // usigma1+  decay only to sigma+ -> pi+ neutron channel
  particle = new G4ParticleDefinition(
           "usigma1+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // usigma+ -> neutron + pi+
  mode = new G4PhaseSpaceDecayChannel("usigma1+",1.0,2,"neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // usigma2+  decay only to sigma+ -> pi0 proton channel
  particle = new G4ParticleDefinition(
           "usigma2+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> proton + pi0
  mode = new G4PhaseSpaceDecayChannel("usigma2+", 1.0,2,"proton","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  // decaying lambda
  particle = new G4ParticleDefinition(
				      "ulambda",    1.115684*GeV,  2.501e-12*MeV,         0.0,
				      1,              +1,             0,
				      0,               0,             0,
				      "baryon",               0,            +1,        3122,
				      false,       0.0*ns,          NULL);

  //create Decay Table
  decayTable = new G4DecayTable();

  // create decay channels
  G4VDecayChannel** mode2 = new G4VDecayChannel*[2];
  // lambda -> proton + pi-
  mode2[0] = new G4PhaseSpaceDecayChannel("ulambda",0.639,2,"proton","pi-");
  // lambda -> neutron + pi0
  mode2[1] = new G4PhaseSpaceDecayChannel("ulambda",0.358,2,"neutron","pi0");

  for (G4int index=0; index <2; index++ ) decayTable->Insert(mode2[index]);
  delete [] mode2;

  particle->SetDecayTable(decayTable);


  particle = new G4ParticleDefinition(
				      "uxi-",    1.32132*GeV,  4.02e-12*MeV,    -1*eplus,
				      1,              +1,             0,
				      1,              -1,             0,
				      "baryon",               0,            +1,        3312,
				      false,       0.0*ns,          NULL);

  //create Decay Table
  decayTable =  new G4DecayTable();
  // Xi- -> lambda + pi-
  mode = new G4PhaseSpaceDecayChannel("uxi-",1.000,2,"lambda","pi-");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);



  //ConstructStableHyperons()
    // ssigma+ non-decay sigma+
  particle = new G4ParticleDefinition(
           "ssigma+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   true,                0,          NULL);


  // sigma1+  decay only to sigma+ -> pi+ neutron channel
  particle = new G4ParticleDefinition(
           "sigma1+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0799*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> neutron + pi+
  mode = new G4PhaseSpaceDecayChannel("sigma1+",1.0,2,"neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // sigma2+  decay only to sigma+ -> pi0 proton channel
  particle = new G4ParticleDefinition(
           "sigma2+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0799*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> proton + pi0
  mode = new G4PhaseSpaceDecayChannel("sigma2+", 1.0,2,"proton","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  // Construct tagged proton
  G4Ions* particle_p;
  particle_p = new G4Ions(
		"proton_d",   0.9382723*GeV,       0.0*MeV,       eplus, 
		    1,              +1,             0,          
		    1,              +1,             0,             
	     "baryon",               0,            +1,        2212,
		 true,            -1.0,          NULL,
		false,       "nucleon",             0,
                  0.0 
		);


  
  //difine two pi- for VP
  particle = new G4ParticleDefinition(
				      "pi-1",     0.13957*GeV,  2.5284e-14*MeV,      -1.*eplus,
				      0,              -1,              0,
				      2,              -2,             -1,
				      "meson",               0,              0,            -211,
				      false,       26.030*ns,           NULL,
				      false,            "pi");

  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("pi-1",1.00,2,"mu-","anti_nu_mu");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  particle = new G4ParticleDefinition(
				      "pi-2",     0.13957*GeV,  2.5284e-14*MeV,      -1.*eplus,
				      0,              -1,              0,
				      2,              -2,             -1,
				      "meson",               0,              0,            -211,
				      false,       26.030*ns,           NULL,
				      false,            "pi");

  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("pi-2",1.00,2,"mu-","anti_nu_mu");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);




  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  G4double mp  = 0.9382723;
  G4double mpi = 0.13957;
  G4double p = 1.05;
  G4double sqrt_s = sqrt(mpi*mpi + mp*mp + 2*mp*sqrt(mpi*mpi+p*p));

  particle = new G4ParticleDefinition(
				      "VP_PiMinusP_1.05",     sqrt_s*GeV,  2.501e-12*MeV,      0.0,
				      1,              +1,              0,
				      0,               0,              0,
				      "baryon",               0,             +1,     3122,
				      false,       0.*ns,           NULL,
				      false,  "VP_PiMinusP_1.05");

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("VP_PiMinusP_1.05",1.0,4,"pi-1","pi-2","pi+","proton");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HIHRPhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
}  
