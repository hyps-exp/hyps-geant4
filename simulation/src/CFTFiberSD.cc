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

#include "CFTFiberSD.hh"
#include "CFTFiberHit.hh"
#include "DetectorID.hh"
#include "DCGeomMan.hh"
#include "RadDeg.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTFiberSD::CFTFiberSD(G4String name)
:G4VSensitiveDetector(name)
{
  detectorname = name;
  collectionName.insert(G4String( "CFTFiberCollection" ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTFiberSD::~CFTFiberSD()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberSD::Initialize(G4HCofThisEvent*HCE)
{
  //G4cout << "####CFTFiberSD::Initialize " << detectorname << G4endl;

  FiberCollection = new CFTFiberHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for (G4int layer=0;layer<NumOfPlaneCFT;layer++) {
    HitID[layer].clear();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CFTFiberSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  //G4cout << "####CFTFiberSD::ProcessHits " << detectorname  << G4endl;
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double time = aStep->GetTrack()->GetGlobalTime();

  G4double stepl = 0.;
  G4String particleName;
  particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();

   
  if ((edep==0.)&&(stepl==0.)) return false;      
  //if (particleName == "e-") return false;      

  G4int particleId = aStep->GetTrack()->GetTrackID();
  G4int parentId = aStep->GetTrack()->GetParentID();

  if (particleName == "proton") {
    if (parentId>0) {
      particleName += "_d";
    }
  }


  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 

  //theTouchable->MoveUpHistory();
  G4int layerID = -1;
  G4int segmentID = -1;
  G4String name = physVol->GetName();

  char candidate[100];

  bool flagFind=false;

  int layer, segment;
  if (sscanf(name.c_str(),"physCFTFiberPHI%d_%d", &layer, &segment) == 2) {
    /*
    const DCGeomMan & geomMan=DCGeomMan::GetInstance();
    int lnum = 0; // target
    G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localPos = geomMan.Global2LocalPos(lnum, hitpos);
    double phi = calcPhi(localPos);
    */

    layerID = (layer-1)*2+1;
    segmentID = segment;
    //G4cout << layerID << " : pii " << phi << ", segment " << segmentID << std::endl;
  } else if (sscanf(name.c_str(),"physCFTFiberUV%d", &layer) == 1) {
    layerID = (layer-1)*2;

    const DCGeomMan & geomMan=DCGeomMan::GetInstance();
    int lnum = 0; // target
    G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localPos = geomMan.Global2LocalPos(lnum, hitpos);

    double phi = calcPhi(localPos);
    double z = localPos.z() + Target_Length/2.;

    //G4cout << layerID << " : pii-" << phi << ", z-" << z << std::endl;

    if (layerID == CFT_U1 || layerID == CFT_U3) {
      double z0 = z/mm - CFT_Z/360.*phi;
      if (z0 <= -CFT_Z/(double)NumOfSegCFT[layerID]/2)
	z0 += CFT_Z;

      segmentID = 
	(int)((z0+(CFT_Z/(double)NumOfSegCFT[layerID]/2.))/(CFT_Z/(double)NumOfSegCFT[layerID]));

    } else if (layerID == CFT_V2 || layerID == CFT_V4) {

      double z0 = z/mm - CFT_Z/360.*(360.-phi);

      if (z0 <= -CFT_Z/(double)NumOfSegCFT[layerID]/2)
	z0 += CFT_Z;

      segmentID = 
	(int)((z0+(CFT_Z/(double)NumOfSegCFT[layerID]/2.))/(CFT_Z/(double)NumOfSegCFT[layerID]));

    } else {
      fprintf(stderr, "CFTFiberSD::ProcessHits Invalid layerId \n");
      fprintf(stderr, "layerID = %d\n", layerID);
      return false;
    }
  }

  //G4cout << "layerID = " << layerID << ", segmentID = " << segmentID << G4endl;
  if (!((layerID >=0 && layerID<=NumOfPlaneCFT) && (segmentID >=0 && segmentID <NumOfSegCFT[layerID]) )) {
    fprintf(stderr, "CFTFiberSD::ProcessHits Invalid layerId or Invalid segmentID \n");
    fprintf(stderr, "layerID = %d, segmentID = %d \n", layerID, segmentID);
    return false;
  }

  std::map<int , int>::iterator it = HitID[layerID].begin();
  while (it != HitID[layerID].end()) {
    int seg = (*it).first;
    int index = (*it).second;
    //G4cout << "layerID : " << layerID << " segmentID : " << segmentID 
    //<< " seg : " << seg << ", index : " << index << G4endl;
    if (seg == segmentID) {
      (*FiberCollection)[index]->AddAbs(edep,stepl);
      return true;
    }
    ++it;
  }

  CFTFiberHit* fiberHit = new CFTFiberHit();
  fiberHit->SetLayerID(layerID);
  fiberHit->SetSegmentID(segmentID);
  fiberHit->AddAbs(edep, stepl);
  fiberHit->SetTime(time);
  fiberHit->SetParticleID(particleId);
  fiberHit->SetParentID(parentId);
  fiberHit->SetParticleName(particleName);
  int index = FiberCollection->insert(fiberHit) - 1;

  HitID[layerID].insert( std::pair<int, int>(segmentID, index));
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //G4cout << "####CFTFiberSD::EndOfEvent " << detectorname  << G4endl;
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  //G4cout << detectorname << " HCID = " << HCID  << G4endl;
  HCE->AddHitsCollection(HCID,FiberCollection);
}

G4double CFTFiberSD::calcPhi(G4ThreeVector pos)
{
  if (pos.x()>=0. && pos.y()>=0)
    return atan(pos.y()/pos.x())*Rad2Deg;
  else if (pos.x()<0. && pos.y()>=0. )
    return 180. + atan(pos.y()/pos.x())*Rad2Deg;
  else if (pos.x()<0. && pos.y()<0. )
    return 180. + atan(pos.y()/pos.x())*Rad2Deg;
  else if (pos.x()>=0. && pos.y()<0. )
    return 360. + atan(pos.y()/pos.x())*Rad2Deg;
}

G4double CFTFiberSD::FiberLineFunction(double lx, double theta, double x)
{
  double x0 = lx*cos(theta*Deg2Rad);
  double z0 = lx*sin(theta*Deg2Rad);

  return -1./tan(theta*Deg2Rad)*(x-x0)+z0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

