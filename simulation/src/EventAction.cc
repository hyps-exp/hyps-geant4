/*
  EventAction.cc
  2007/4  K.Shirotori
*/

#include "EventAction.hh"
#include "Analysis.hh"
#include "CFTFiberHit.hh"
#include "CrystalHit.hh"
#include "SKSChamberHit.hh"
#include "SKSCounterHit.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4VTrajectory.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VTrajectoryPoint.hh"

#include "ConfMan.hh"

#include <iomanip>

EventAction::EventAction( Analysis *analysisManager )
  : G4UserEventAction(), anaMan(analysisManager), fAlwaysDraw(false),
    CFTFibercolID(-1), CrystalcolID(-1), ChambercolID(-1), CountercolID(-1), VirtualPlanecolID(-1)
{}


EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction( const G4Event *anEvent )
{
  G4int eventID = anEvent->GetEventID();
  if( eventID%1000==0 )
    G4cout << "Event: " << std::setw(8) << eventID << G4endl;
  
  G4SDManager *SDManager = G4SDManager::GetSDMpointer();
  if( CFTFibercolID<0 || CrystalcolID<0 || ChambercolID || VirtualPlanecolID){

    CFTFibercolID = SDManager->GetCollectionID( G4String( "CFTFiberCollection" ) );
    CrystalcolID = SDManager->GetCollectionID( G4String( "CrystalCollection" ) ); 
    ChambercolID = SDManager->GetCollectionID( G4String( "chamberSD" ) ); 
    CountercolID = SDManager->GetCollectionID( G4String( "counterSD" ) ); 
    VirtualPlanecolID = SDManager->GetCollectionID( G4String( "virtualPlaneSD" ) ); 
  }
#if 0
  G4cout << "SksEventAction::BeginOfEventAction "
	 << "Collection IDs :" 
	 << CFTFibercolID << " " << CrystalcolID << " "
	 << G4endl;
#endif
  
  if( anaMan ) anaMan->BeginOfEvent( );

}

void EventAction::EndOfEventAction( const G4Event *anEvent )
{

  G4VVisManager *visManager = G4VVisManager::GetConcreteInstance();

  //G4cout << "EventAction::EndOfEventAction" << G4endl;

  anaMan->EndOfEvent( anEvent );

  G4bool fTriggered = anaMan->GetTriggerStatus();


  if( visManager && ( fAlwaysDraw || fTriggered ) ){
    G4UImanager *UI = G4UImanager::GetUIpointer();
    // Trajectories
    G4TrajectoryContainer *trajectoryContainer =
      anEvent->GetTrajectoryContainer();
    if( trajectoryContainer ){
      G4int nTrajectories = trajectoryContainer->entries();
      for( G4int i=0; i<nTrajectories; i++ ){
        G4VTrajectory *trj = (*trajectoryContainer)[i];
	//if (trj->GetParticleName() != "skaon+")
	//trj->DrawTrajectory(10);
	trj->DrawTrajectory();
#if 0
	G4int nPoint=trj->GetPointEntries();
	for( G4int j=0; j<nPoint; ++j ){
	  G4VTrajectoryPoint *tp=trj->GetPoint(j);
	  G4ThreeVector position=tp->GetPosition();
	  G4cout << position << G4endl;
	}
#endif

      }
    }
  }

}
