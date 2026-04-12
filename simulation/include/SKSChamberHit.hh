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
// $Id: ExN03CalorHit.hh,v 1.4 2001/10/10 14:58:11 maire Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SKSChamberHit_h
#define SKSChamberHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "common.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SKSChamberHit : public G4VHit
{
 public:

   SKSChamberHit();
  ~SKSChamberHit();
   SKSChamberHit(const SKSChamberHit&);
   const SKSChamberHit& operator=(const SKSChamberHit&);
   G4int operator==(const SKSChamberHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   void Draw();
   void Print();
      
 public:
   void SetDetectorID(G4int i) {detectorID = i;};
   void SetHitWire(G4int wire) {hitWire = wire;};
   void SetDRlength(G4double dl) {dlength = dl;};
   void SetSensitiveHit(bool sensitive) {sensitiveHit = sensitive;};
   G4int GetDetectorID() { return detectorID;};
   G4int GetHitWire()    { return hitWire;};
   G4double GetDRlength() {return dlength;};
   bool GetSensitiveHit() { return sensitiveHit;};
 private:
   G4int detectorID;
   G4int hitWire;
   G4double dlength;
  bool     sensitiveHit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<SKSChamberHit> SKSChamberHitsCollection;

extern G4Allocator<SKSChamberHit> SKSChamberHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* SKSChamberHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) SKSChamberHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void SKSChamberHit::operator delete(void* aHit)
{
  SKSChamberHitAllocator.FreeSingle((SKSChamberHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
