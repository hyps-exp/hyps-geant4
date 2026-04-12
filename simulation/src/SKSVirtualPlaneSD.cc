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

#include "SKSVirtualPlaneSD.hh"

#include "SKSVirtualPlaneHit.hh"
#include "CFTDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"

#include "G4ios.hh"

#include "DCGeomMan.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace CLHEP;

SKSVirtualPlaneSD::SKSVirtualPlaneSD(G4String name)
  :G4VSensitiveDetector(name)
{
  detectorname = name;
  collectionName.insert(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSVirtualPlaneSD::~SKSVirtualPlaneSD()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSVirtualPlaneSD::Initialize(G4HCofThisEvent* HCE)
{
  //G4cout << "####SKSVirtualPlaneSD::Initialize " << detectorname << G4endl;

  CalCollection = new SKSVirtualPlaneHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 

  for (G4int i=0;i<DCPlaneMax;i++)
    for (G4int j=0; j<DCMultiMax; j++) {
      HitID[i][j] = -1;
    }
  
  for (G4int i=0;i<DCPlaneMax;i++)
    Multi[i]=0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SKSVirtualPlaneSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  const DCGeomMan &geomMan = DCGeomMan::GetInstance();

  G4double dlength;
  G4bool   MultiFlag = false;
  G4int    hitWire;

  G4double edep = aStep->GetTotalEnergyDeposit();

  G4double stepl = 0.;
  G4String particleName;

  particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();

  if ((edep==0.)&&(stepl==0.)) return false;      
  if (particleName == "e-") return false;      

  G4ThreeVector position = aStep->GetTrack()->GetPosition();

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4String name = physVol->GetName();
  //G4cout << name << ", " << particleName << G4endl;

  int lnum = geomMan.GetDetectorId(name);
  G4ThreeVector localPos = geomMan.Global2LocalPos(lnum, position);

  if (Multi[lnum] >= DCMultiMax) {
    fprintf(stderr, "Too much multiplicity %d at plane%d\n", Multi[lnum], lnum);
    return false;
  }

  if (HitID[lnum][Multi[lnum]]==-1)
    { 
      SKSVirtualPlaneHit* calHit = new SKSVirtualPlaneHit();
      calHit->SetDetectorID(lnum);
      calHit->SetHitPosition(localPos.x(), localPos.y());
      HitID[lnum][Multi[lnum]] = CalCollection->insert(calHit) - 1;
      Multi[lnum]++;
    }
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSVirtualPlaneSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //G4cout << "####SKSVirtualPlaneSD::EndOfEvent " << detectorname  << G4endl;
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  //G4cout << detectorname << " HCID = " << HCID  << G4endl;
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSVirtualPlaneSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSVirtualPlaneSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSVirtualPlaneSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
