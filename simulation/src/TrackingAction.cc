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

#include "TrackingAction.hh"
#include "Analysis.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4Trajectory.hh"

#include "G4SystemOfUnits.hh"

TrackingAction::TrackingAction(Analysis* ana)
  : anaMan(ana)
{
  // messenger defaults
  StoreTrajectoriesOfSecondariesFlag = true;
}

TrackingAction::~TrackingAction()
{
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // if !StoreTrajectoriesOfSecondariesFlag, create trajectory only for primaries
  if ( (StoreTrajectoriesOfSecondariesFlag) || (aTrack->GetParentID()==0) )
    { fpTrackingManager->SetStoreTrajectory(true); }
  else
    { fpTrackingManager->SetStoreTrajectory(false); }
  
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //G4cout << "TrackingAction::PostUserTrackingAction" << G4endl;
  G4String particleName = aTrack->GetDefinition()->GetParticleName();

  if (particleName == "sigma+" || particleName == "sigma1+" || particleName == "sigma2+" || particleName == "sigma-" || particleName == "lambda") {
    sigmaID = aTrack->GetTrackID();
    G4ThreeVector decayPos = aTrack->GetPosition();
    anaMan->SetDecayPos(decayPos);    
    return;
  } 

  if (particleName == "kaon0" || particleName == "kaon0S") {
    k0ID = aTrack->GetTrackID();
    G4ThreeVector decayPos = aTrack->GetPosition();
    anaMan->SetDecayPos(decayPos);    
    return;
  } 


  if (aTrack->GetParentID() == sigmaID) {
    //G4cout << particleName << G4endl;
    if (particleName == "pi+" || particleName == "pi0" || particleName == "pi-"){
      G4ThreeVector piMomentum=aTrack->GetVertexMomentumDirection();
      G4double Ekin = aTrack->GetVertexKineticEnergy()/GeV;
      G4double M = aTrack->GetDefinition()->GetPDGMass()/GeV;
      G4double mom = sqrt((Ekin+M)*(Ekin+M)-M*M);

      anaMan->SetDecayPiMomentum(mom*piMomentum);
      //G4cout << particleName << " : ( " << piMomentum.x() << ", "
      //<< piMomentum.y() << ", " << piMomentum.z() << ") " << G4endl;
    } else if (particleName == "proton" || particleName == "neutron"){
      G4ThreeVector nucleonMomentum=aTrack->GetVertexMomentumDirection();
      G4double Ekin = aTrack->GetVertexKineticEnergy()/GeV;
      G4double M = aTrack->GetDefinition()->GetPDGMass()/GeV;
      G4double mom = sqrt((Ekin+M)*(Ekin+M)-M*M);
      anaMan->SetDecayNucleonMomentum(mom*nucleonMomentum);

      if (particleName == "neutron")
	anaMan->SetDecayFlag(0);  
      else if (particleName == "proton")
	anaMan->SetDecayFlag(1);

      //G4cout << particleName << " : ( " 
      //<< nucleonMomentum.x() << ", "
      //<< nucleonMomentum.y() << ", " 
      //<< nucleonMomentum.z() << ") " << G4endl;
    }
  } else if (aTrack->GetParentID() == k0ID) {
    if (particleName == "pi+"){
      G4ThreeVector piMomentum=aTrack->GetVertexMomentumDirection();
      G4double Ekin = aTrack->GetVertexKineticEnergy()/GeV;
      G4double M = aTrack->GetDefinition()->GetPDGMass()/GeV;
      G4double mom = sqrt((Ekin+M)*(Ekin+M)-M*M);

      anaMan->SetDecayPiMomentum(mom*piMomentum);
      //G4cout << particleName << " : ( " << piMomentum.x() << ", "
      //<< piMomentum.y() << ", " << piMomentum.z() << ") " << G4endl;
    } else if (particleName == "pi-"){
      G4ThreeVector nucleonMomentum=aTrack->GetVertexMomentumDirection();
      G4double Ekin = aTrack->GetVertexKineticEnergy()/GeV;
      G4double M = aTrack->GetDefinition()->GetPDGMass()/GeV;
      G4double mom = sqrt((Ekin+M)*(Ekin+M)-M*M);
      anaMan->SetDecayNucleonMomentum(mom*nucleonMomentum);

    }


  }

  /*
  if (particleName == "proton" && aTrack->GetTrackID() == 1) {
    G4ThreeVector stopPos = aTrack->GetPosition();
    anaRoot->SetProtonStopPos(stopPos);    
    return;
  }
  */
}


