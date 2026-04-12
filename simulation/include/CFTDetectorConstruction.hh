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
// $Id: ExN03DetectorConstruction.hh,v 1.6 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CFTDetectorConstruction_h
#define CFTDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DetectorID.hh"
#include "G4ThreeVector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;

class G4Material;
class SKSChamberSD;
class SKSCounterSD;
class SKSVirtualPlaneSD;

struct MaterialList;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CFTDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    CFTDetectorConstruction();
   ~CFTDetectorConstruction();

  public:
     
     G4VPhysicalVolume* Construct();

  private:
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physWorld;    //pointer to the physical World

  SKSChamberSD *chamberSD;
  SKSCounterSD *counterSD;
  SKSVirtualPlaneSD *virtualPlaneSD;

  MaterialList *mList_;
  MaterialList *DefineMaterials();
  void  ConstructCFT(G4VPhysicalVolume *pMother);
  void  ConstructTarget(G4VPhysicalVolume *pMother);
  void  ConstructKURAMA(G4VPhysicalVolume *pMother);
  void  ConstructSdcIn(G4VPhysicalVolume *pMother);
  void  ConstructSdcOut(G4VPhysicalVolume* pMother);
  void  ConstructBcOut(G4VPhysicalVolume *pMother);
  void  ConstructVP(G4VPhysicalVolume *pMother);
  void  ConstructCH(G4VPhysicalVolume* pMother);
  void  ConstructCHwithBeamHole(G4VPhysicalVolume* pMother);
  void  ConstructSH(G4VPhysicalVolume* pMother);
  void  ConstructT0(G4VPhysicalVolume* pMother);
  void  ConstructAC(G4VPhysicalVolume* pMother);
  void  ConstructeVeto(G4VPhysicalVolume* pMother);
  void  ConstructTOF(G4VPhysicalVolume* pMother);
  void  ConstructBeamDump(G4VPhysicalVolume* pMother);
  void  ConstructS2S(G4VPhysicalVolume* pMother);
  G4double GetDistance(G4ThreeVector pos1, G4ThreeVector pos2);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

