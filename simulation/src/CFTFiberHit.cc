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
// $Id: ExN03CalorHit.cc,v 1.4 2001/10/10 14:58:12 maire Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CFTFiberHit.hh"

G4Allocator<CFTFiberHit> CFTFiberHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTFiberHit::CFTFiberHit()
{
  layerID = -1;
  segmentID = -1;
  EdepAbs = 0.;
  TrackLengthAbs = 0.;
  Time = 0.;
  particleId = -1;
  parentId = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTFiberHit::~CFTFiberHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTFiberHit::CFTFiberHit(const CFTFiberHit& right)
  : G4VHit()
{
  layerID = right.layerID;
  segmentID = right.segmentID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  particleId = right.particleId;
  parentId = right.parentId;
  particleName = right.particleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const CFTFiberHit& CFTFiberHit::operator=(const CFTFiberHit& right)
{
  layerID = right.layerID;
  segmentID = right.segmentID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  particleId = right.particleId;
  parentId = right.parentId;
  particleName = right.particleName;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int CFTFiberHit::operator==(const CFTFiberHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberHit::Draw()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CFTFiberHit::Print()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

