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

#include "SKSCounterHit.hh"

G4Allocator<SKSCounterHit> SKSCounterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSCounterHit::SKSCounterHit()
{
  detectorID = -1;
  channelID = -1;
  EdepAbs = 0.;
  TrackLengthAbs = 0.;
  Time = 0.;
  Position[XCOORD]=-9999.9;
  Position[YCOORD]=-9999.9;
  Position[ZCOORD]=-9999.9;
  lcCherenkov = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSCounterHit::~SKSCounterHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SKSCounterHit::SKSCounterHit(const SKSCounterHit& right)
  : G4VHit()
{
  detectorID = right.detectorID;
  channelID = right.channelID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  Position[XCOORD]=right.Position[XCOORD];
  Position[YCOORD]=right.Position[YCOORD];
  Position[ZCOORD]=right.Position[ZCOORD];
  lcCherenkov = right.lcCherenkov;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SKSCounterHit& SKSCounterHit::operator=(const SKSCounterHit& right)
{
  detectorID = right.detectorID;
  channelID = right.channelID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  Position[XCOORD]=right.Position[XCOORD];
  Position[YCOORD]=right.Position[YCOORD];
  Position[ZCOORD]=right.Position[ZCOORD];
  lcCherenkov = right.lcCherenkov;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SKSCounterHit::operator==(const SKSCounterHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterHit::Draw()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SKSCounterHit::Print()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

