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

#include "CrystalSD.hh"

#include "CrystalHit.hh"
#include "DetectorID.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

#include "DCGeomMan.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalSD::CrystalSD(G4String name)
:G4VSensitiveDetector(name)
{
  detectorname = name;
  collectionName.insert( G4String( "CrystalCollection") );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalSD::~CrystalSD()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalSD::Initialize(G4HCofThisEvent*HCE)
{
  //G4cout << "####SKSCalorimeterSD::Initialize " << detectorname << G4endl;

  CalCollection = new CrystalHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  HitID.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CrystalSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
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


  if (particleName == "proton") {
    G4int parentId = aStep->GetTrack()->GetParentID();
    if (parentId>0) {
      particleName += "_d";
    }
  }

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4int replicaNum = theTouchable->GetReplicaNumber();

  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4ThreeVector position = physVol->GetObjectTranslation();
  
  G4int layerID = -1, channelID = -1;
  G4String name = physVol->GetName();
  char candidate[30];
  G4int i;

  int layer;
  int ch;

  if (sscanf(name.c_str(),"physBGO%d", &ch) == 1) {
    channelID = ch;
  } else {
    G4cout <<  "CrystalSD::ProcessHits No such detector " <<  name << G4endl;
    return false;
  }

  if (channelID >= NumOfBGOUnit*(NumOfBGOInOneUnit+NumOfBGOInOneUnit2) || channelID <0) {
    //if (channelID >= 25 || channelID <0) {
    fprintf(stderr, "CrystalSD::ProcessHits, ChannelID  %d is greater than Max Segment\n", channelID);
    return false;
  }

  int lnum = geomMan.GetDetectorId("Target");
  G4ThreeVector trackpos = aStep->GetTrack()->GetPosition();
  G4ThreeVector localTrackPos = geomMan.Global2LocalPos(lnum, trackpos);

  std::map<int , int>::iterator it = HitID.begin();
  while (it != HitID.end()) {
    int seg = (*it).first;
    int index = (*it).second;
    if (seg == channelID) {
      (*CalCollection)[index]->AddAbs(edep,stepl);
      return true;
    }
    ++it;
  }

  CrystalHit* calHit = new CrystalHit();
  //calHit->SetLayerID(layerID);
  calHit->SetChannelID(channelID);
  calHit->AddAbs(edep, stepl);
  calHit->SetTime(time);
  calHit->SetParticleName(particleName);
  calHit->SetPos(localTrackPos.x(), localTrackPos.y(), localTrackPos.z());
  
  int index =  CalCollection->insert(calHit) - 1;
  HitID.insert( std::pair<int, int>(channelID, index));

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //G4cout << "####SKSCalorimeterSD::EndOfEvent " << detectorname  << G4endl;
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  //G4cout << detectorname << " HCID = " << HCID  << G4endl;
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

