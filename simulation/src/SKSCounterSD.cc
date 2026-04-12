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
// $Id: ExN03CalorimeterSD.cc,v 1.5 2001/10/10 14:58:12 maire Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SKSCounterSD.hh"

#include "SKSCounterHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

#include "DCGeomMan.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSCounterSD::SKSCounterSD(G4String name)
:G4VSensitiveDetector(name)
{
  detectorname = name;
  collectionName.insert(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSCounterSD::~SKSCounterSD()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterSD::Initialize(G4HCofThisEvent*HCE)
{
  //G4cout << "####SKSCalorimeterSD::Initialize " << detectorname << G4endl;

  CalCollection = new SKSCounterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for (G4int i=0;i<NumOfCounter;i++)
    for (G4int j=0; j<CounterSegmentMax; j++)
      {HitID[i][j] = -1;};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SKSCounterSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  const DCGeomMan &geomMan = DCGeomMan::GetInstance();

  //G4cout << "####SKSCalorimeterSD::ProcessHits " << detectorname  << G4endl;
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double time = aStep->GetTrack()->GetGlobalTime();

  G4double stepl = 0.;
  G4String particleName;
  particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();
      
  if ((edep==0.)&&(stepl==0.)) return false;      
  //if (particleName == "e-") return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4ThreeVector position = physVol->GetObjectTranslation();

  G4double mom  = aStep->GetTrack()->GetDynamicParticle()->GetTotalMomentum();
  G4double mass = aStep->GetTrack()->GetDefinition()->GetPDGMass();
  G4double beta = -1.0;
  if (mass>0)
    beta = fabs(mom/mass)/sqrt(1.0+(mom/mass)*(mom/mass));

  G4int detectorID = 0;
  G4int lcCherenkov = 0;
  G4int channelID = -1;
  G4String name = physVol->GetName();
  char candidate[30];
  G4int i;
  G4int TofId = geomMan.GetDetectorId("FTOF");
  G4int flagCherenkov=0;
  G4int detectorIndex = 0;


  if (name.find("sFTOF_phys")==0) {
    int ch, sh;
    sscanf(name.c_str(),"sFTOF_phys%d", &ch);
    detectorID = geomMan.GetDetectorId("FTOF");
    channelID = ch-1;
    detectorIndex = 0;
    G4cout << particleName << G4endl;
    if (particleName == "pi-")
      lcCherenkov = 0;
    else if (particleName == "kaon+")
      lcCherenkov = 1;
    else if (particleName == "proton")
      lcCherenkov = 2;
    else {
      lcCherenkov = 3;
    }
    // SKSCounterHit::SetLcCherenkov(lcCherenkov);
  } else if (name.find("sCH_phys")==0) {
    int ch;
    if (sscanf(name.c_str(),"sCH_phys%d", &ch) == 1) {
      detectorID = geomMan.GetDetectorId("CH");
      channelID = ch-1;
      detectorIndex = 1;
    } else if (sscanf(name.c_str(),"sCH_physD%d", &ch) == 1) {
      detectorID = geomMan.GetDetectorId("CH");
      channelID = ch-1;
      detectorIndex = 1;
    } 
  } else if (name.find("physPiV")==0) {
    int ch;
    sscanf(name.c_str(),"physPiV%d", &ch);
    detectorID = geomMan.GetDetectorId("PIV");
    channelID = ch;
    detectorIndex = 2;
  } else if (name.find("sAC_phys")==0) {
    double index = 1.06;
    if (beta > 1./index) {
      int ch;
      sscanf(name.c_str(),"sAC_phys%d", &ch);
      detectorID = geomMan.GetDetectorId("AC");
      channelID = ch;
      detectorIndex = 3;
      flagCherenkov = 1;
    } else {
      return false;
    }
  } else if (name.find("physBGO_VP")==0) {
    int ch;
    detectorID = 71;
    if (particleName == "proton")
      channelID = 0;
    else if (particleName == "pi+")
      channelID = 1;
    else if (particleName == "pi-")
      channelID = 2;
    else if (particleName == "kaon+" || particleName == "skaon+")
      channelID = 3;
    else {
      channelID = 4;
      //G4cout << particleName << G4endl;
    }
    detectorIndex = 4;
  } else {
    G4cout <<  "SKSCounterSD::ProcessHits No such detector " << name << G4endl;
    return false;
  } 
    

 if (detectorIndex >= CounterSegmentMax) {
    fprintf(stderr, "ChannelID  %d is greater than Max Segment at %d Counter\n", channelID, detectorID);
    return false;
  }

  G4ThreeVector trackpos = aStep->GetTrack()->GetPosition();

  if (HitID[detectorIndex][channelID]==-1)
    { 
      SKSCounterHit* calHit = new SKSCounterHit();
      calHit->SetDetectorID(detectorID);
      calHit->SetChannelID(channelID);
      calHit->AddAbs(edep, stepl);
      calHit->SetTime(time);
      calHit->SetPos(trackpos.x(), trackpos.y(), trackpos.z());
      calHit->SetLcCherenkov(flagCherenkov);
      HitID[detectorIndex][channelID] = CalCollection->insert(calHit) - 1;
    }
  else
    { 
      (*CalCollection)[HitID[detectorIndex][channelID]]->AddAbs(edep,stepl);
      //(*CalCollection)[HitID[0]]->SetTime(time);
    }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //G4cout << "####SKSCalorimeterSD::EndOfEvent " << detectorname  << G4endl;
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  //G4cout << detectorname << " HCID = " << HCID  << G4endl;
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

