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

#include "CrystalHit.hh"

G4Allocator<CrystalHit> CrystalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalHit::CrystalHit()
{
  layerID = -1;
  channelID = -1;
  EdepAbs = 0.;
  TrackLengthAbs = 0.;
  Time = 0.;
  Position[0]=-9999.9;
  Position[1]=-9999.9;
  Position[2]=-9999.9;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalHit::~CrystalHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalHit::CrystalHit(const CrystalHit& right)
  : G4VHit()
{
  layerID = right.layerID;
  channelID = right.channelID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  Name = right.Name;
  Position[0]=right.Position[0];
  Position[1]=right.Position[1];
  Position[2]=right.Position[2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const CrystalHit& CrystalHit::operator=(const CrystalHit& right)
{
  layerID = right.layerID;
  channelID = right.channelID;
  EdepAbs = right.EdepAbs;
  TrackLengthAbs = right.TrackLengthAbs;
  Time = right.Time;
  Name = right.Name;
  Position[0]=right.Position[0];
  Position[1]=right.Position[1];
  Position[2]=right.Position[2];

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int CrystalHit::operator==(const CrystalHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalHit::Draw()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalHit::Print()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

