/*
  SteppingAction.cc
  2007/4  K.Shirotori
*/

#include "SteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"

#include "Analysis.hh"
#include "DCGeomMan.hh"

#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction( Analysis* ana)
  :  anaMan(ana)
{ 
  primaryVertex = G4ThreeVector(-9999., -9999., -9999.); 
}

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* theStep )
{ 
  
  int i;
  G4Track * theTrack = theStep->GetTrack();
  
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(theStep->GetPreStepPoint()->GetTouchable());
  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4String name = physVol->GetName();
  //G4cout<< "Name=" << name << G4endl;
  G4String particle = theTrack->GetDefinition()->GetParticleName();
  G4ThreeVector pos = theTrack->GetPosition();
  
  G4StepPoint *preStepPoint = theStep->GetPreStepPoint();
  G4TouchableHandle Touchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector posl = Touchable->GetHistory()->
    GetTopTransform().TransformPoint( pos );
  
  G4ThreeVector position = theStep->GetPreStepPoint()->GetPosition();
  int  pid = theTrack->GetTrackID();  
  int  parentid = theTrack->GetParentID();  
  double charge = theStep->GetTrack()->GetDefinition()->GetPDGCharge();
  
  
#if 1
  if (theTrack->GetParentID()==0&&particle=="ssigma+"){
    primaryVertex = anaMan->GetPrimaryVertex();
    const DCGeomMan & geomMan=DCGeomMan::GetInstance();
    int lnum = geomMan.GetDetectorId("Target");      
    primaryVertex = geomMan.Local2GlobalPos(lnum, primaryVertex);
    
    if (position.y() > primaryVertex.y() || position.x() < primaryVertex.x()) {
      theTrack->SetTrackStatus(fStopAndKill);      
      primaryVertex = G4ThreeVector(-9999., -9999., -9999.);
    }
    return;
  }
#endif
  

  //All particle which hit KURAMA are killed.
  if (    name == "upGuard_U_phys" || name == "upGuard_D_phys"
	  || name == "upGuard_L_phys" || name == "upGuard_R_phys" 
	  || name == "Yoke_U_phys" || name == "Yoke_D_phys" 
	  || name == "Yoke_L_phys" || name == "Yoke_R_phys" 
	  || name == "Yoke_L_GapSpace_phys" || name == "Yoke_R_GapSpace_phys" 
	  || name == "downGuard_U_phys" || name == "downGuard_D_phys"
	  || name == "downGuard_L_phys" || name == "downGuard_R_phys" 
	  ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  
  /* 10mよりも長く走る粒子はストップさせる */
  if (theTrack->GetTrackLength() >= 10.0*m) {
    /*    
    std::cout << "SteppingAction::UserSteppingAction  Too long track length :"
	      << theTrack->GetTrackLength() << std::endl;
    std::cout << "This Track is FORCED KILLED." << std::endl;
    */
    theTrack->SetTrackStatus(fStopAndKill);      
    return;
  }
}

