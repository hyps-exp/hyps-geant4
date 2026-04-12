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

#include "SKSChamberSD.hh"

#include "SKSChamberHit.hh"
#include "CFTDetectorConstruction.hh"
#include "ConfMan.hh"

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

SKSChamberSD::SKSChamberSD(G4String name)
  :G4VSensitiveDetector(name), ReactionMode(-1), lnum_SFTX(-1), lnum_SFTU(-1), lnum_SFTV(-1)
{
  detectorname = name;
  collectionName.insert(name);

  for (int i=0; i<DCPlaneMax; i++) {
    cham_eff[i]=1.0;
  }

  //lnum_SFTX= DCGeomMan::GetInstance().GetDetectorId("SFT-x-1");
  //lnum_SFTU= DCGeomMan::GetInstance().GetDetectorId("SFT-u-1");
  //lnum_SFTV= DCGeomMan::GetInstance().GetDetectorId("SFT-v-1");

  ConfMan *confMan = ConfMan::GetConfManager();
  ReactionMode = confMan->ReactionMode();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSChamberSD::~SKSChamberSD()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSChamberSD::Initialize(G4HCofThisEvent* HCE)
{
  //G4cout << "####SKSChamberSD::Initialize " << detectorname << G4endl;

  CalCollection = new SKSChamberHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 

  for (G4int i=0;i<DCPlaneMax;i++)
    for (G4int j=0; j<DCMultiMax; j++) {
      HitID[i][j] = -1;
      HitWire[i][j] = -1;
      dlength[i][j] = -1.0;
    }
  
  for (G4int i=0;i<DCPlaneMax;i++)
    Multi[i]=0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SKSChamberSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  const DCGeomMan &geomMan = DCGeomMan::GetInstance();

  G4double dlength;
  G4bool   MultiFlag = false;
  G4int    hitWire = -1;

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
  //G4cout << name << G4endl;

  if (name.find("sSH1")==0) {
    int ch, layer;
    sscanf(name.c_str(),"sSH1_%d_phys%d", &layer, &ch);
    hitWire = ch;
    dlength = 0.;

    if (layer==1)
      name = "SH1-0";
    else if (layer==2)
      name = "SH1-1";
    else if (layer==3)
      name = "SH1-2";
    else if (layer==4)
      name = "SH1-3";
  } else if (name.find("sSH2")==0) {
    int ch, layer;
    sscanf(name.c_str(),"sSH2_%d_phys%d", &layer, &ch);
    hitWire = ch;
    dlength = 0.;

    if (layer==1)
      name = "SH2-0";
    else if (layer==2)
      name = "SH2-1";
    else if (layer==3)
      name = "SH2-2";
    else if (layer==4)
      name = "SH2-3";
  }

  int lnum = geomMan.GetDetectorId(name);
  // G4cerr << "Det.ID = " << lnum <<G4endl;

 /* efficiency check */
  bool sensitiveHit=false;
  if (CheckEfficiency(lnum)==0)
    sensitiveHit=false; /* do not hit due to low efficiency */
  else
    sensitiveHit=true;

  if (Multi[lnum] >= DCMultiMax) {
    fprintf(stderr, "Too much multiplicity %d at plane%d\n", Multi[lnum], lnum);
    return false;
  }

  if (hitWire<0)
  // if (hitWire>0)
    hitWire = getHitWire(lnum, position, &dlength);

  if (hitWire < 0) {
    G4cerr <<  "Invalid hitWire" << hitWire << ", Det.ID = " << lnum <<G4endl;
    return false;
  }

  /*
  // pi- 
  if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7 || ReactionMode == 14 || (ReactionMode >= 16 && ReactionMode <= 22)) {
    if (name == "SDC3-x-1" || name == "SDC3-x-2") {
      if (hitWire>=40 && hitWire<= 46)
	return false;
    } else if (name == "SDC4-x-1" || name == "SDC4-x-2") {
      if (hitWire>=59 && hitWire<= 65)
	return false;
    } else if (name == "SDC3-u-1" || name == "SDC3-u-2" || name == "SDC4-u-1" || name == "SDC4-u-2") {
      if (hitWire>=57 && hitWire<= 65)
	return false;
    }
  } else { // pi+
    if (name == "SDC3-x-1" || name == "SDC3-x-2") {
      if (hitWire>=26 && hitWire<= 32)
	return false;
    } else if (name == "SDC4-x-1" || name == "SDC4-x-2") {
      if (hitWire>=79 && hitWire<= 84)
	return false;
    } else if (name == "SDC3-u-1" || name == "SDC3-u-2" || name == "SDC4-u-1" || name == "SDC4-u-2") {
      if (hitWire>=57 && hitWire<= 65)
	return false;
    }
  }
  */
  for (int i=0; i<Multi[lnum]; i++) {
    /*
    G4cout << "hitWire = " << hitWire << ", HitWire[" << lnum << "][" << i
	   << "] = " << HitWire[lnum][i] << G4endl;
    */
    if (hitWire == HitWire[lnum][i]) {
      MultiFlag = true;
      break;
    }
  }
  //G4cout << "MultiFlag =" << MultiFlag << G4endl;

  if (HitID[lnum][Multi[lnum]]==-1 && MultiFlag != true )
    { 
      SKSChamberHit* calHit = new SKSChamberHit();
      calHit->SetDetectorID(lnum);
      calHit->SetHitWire(hitWire);
      calHit->SetDRlength(dlength);
      calHit->SetSensitiveHit(sensitiveHit);
      HitID[lnum][Multi[lnum]] = CalCollection->insert(calHit) - 1;
      HitWire[lnum][Multi[lnum]] = hitWire;
      // here
      // printf("%d-th cham, multi=%d, hitWire=%d dlength=%f\n", 
      // lnum, Multi[lnum], HitWire[lnum][Multi[lnum]], dlength);
      // printf("lnum=%d, pos=(%f, %f, %f)\n", lnum, position.x(),position.y(),position.z());
      // here
      Multi[lnum]++;
    }
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSChamberSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //G4cout << "####SKSChamberSD::EndOfEvent " << detectorname  << G4endl;
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  //G4cout << detectorname << " HCID = " << HCID  << G4endl;
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSChamberSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSChamberSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSChamberSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SKSChamberSD::getHitWire(G4int lnum, G4ThreeVector position, 
			       G4double *dlength)
{
  const DCGeomMan &geomMan = DCGeomMan::GetInstance();

  G4ThreeVector localPos = geomMan.Global2LocalPos(lnum, position);
  G4int hitWire = geomMan.calcWireNumber(lnum, localPos.x());

  //if (lnum>=lnum_SFTX && lnum<=lnum_SFTV)
  if (lnum>=lnum_SFTU && lnum<=lnum_SFTX)
    *dlength = 0;
  else {
    *dlength = fabs(geomMan.calcWirePosition(lnum, hitWire) - localPos.x());
    *dlength += RandGauss::shoot(0.0,geomMan.GetResolution(lnum));
  }

  //if (lnum>=1&&lnum<=10) {
  if (0) {
    //G4cout << "lnum=" << lnum << " ,wire=" << hitWire << ", dlength=" 
      //<< geomMan.calcWirePosition(lnum, hitWire) - localPos.x() 
    //<< *dlength
    //<< G4endl;
    G4cout << "lnum=" << lnum << ", wirePos = "<< geomMan.calcWirePosition(lnum, hitWire)
    <<", localPos=(" << localPos.x() << ", " 
    << localPos.y() << ", " << localPos.z() << ")" << G4endl;

  }



  return hitWire;
}

G4int SKSChamberSD::CheckEfficiency(G4int planeID)
{
  double x;

  x = RandFlat::shoot();
  //printf("plane=%d x=%f eff=%f\n", planeID, x, cham_eff[planeID]);
  if (x>cham_eff[planeID])
    return 0;
  else
    return 1;
}
