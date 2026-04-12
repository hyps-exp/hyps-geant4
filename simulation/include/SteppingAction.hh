/*
  SteppingAction.hh
  2007/4  K.Shirotori
*/

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"

class Analysis;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction( Analysis*);
  ~SteppingAction();

  void UserSteppingAction(const G4Step*);

private:
  Analysis*             anaMan;  
  G4ThreeVector primaryVertex;
};

#endif
