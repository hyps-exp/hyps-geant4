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
// $Id: ExN03CalorimeterSD.hh,v 1.3 2001/07/11 09:58:21 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SKSChamberSD_h
#define SKSChamberSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//class CFTDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "SKSChamberHit.hh"

#define DCPlaneMax 130
#define DCMultiMax 20

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SKSChamberSD : public G4VSensitiveDetector
{
  public:
  
     SKSChamberSD(G4String);
     ~SKSChamberSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();
      G4int getHitWire(G4int lnum, G4ThreeVector position, 
		       G4double *dlength);
      G4int CheckEfficiency(G4int planeID);
  private:
      SKSChamberHitsCollection*  CalCollection;      
      G4int                    HitID[DCPlaneMax][DCMultiMax];
      G4int                    HitWire[DCPlaneMax][DCMultiMax];
      G4double                 dlength[DCPlaneMax][DCMultiMax];
      G4int                    Multi[DCPlaneMax];
      G4String                 detectorname;
      G4double                 cham_eff[DCPlaneMax];
      G4int                    lnum_SFTX;
      G4int                    lnum_SFTU;
      G4int                    lnum_SFTV;
      G4int                    ReactionMode;
};

#endif

