/*
  EventAction.hh
  2007/4  K.Shirotori
*/

#ifndef EventAction_h 
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class Analysis;

class EventAction : public G4UserEventAction
{
public:
  EventAction( Analysis *analysisManager=0 );
  ~EventAction();

public:
  void BeginOfEventAction( const G4Event *anEvent );
  void EndOfEventAction( const G4Event *anEvent );
  
private:
  Analysis *anaMan;
  G4bool fAlwaysDraw;
  
  G4int CFTFibercolID;
  G4int CrystalcolID;
  G4int ChambercolID;
  G4int CountercolID;
  G4int VirtualPlanecolID;
};

#endif

  
