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

#ifndef SKSVirtualPlaneHit_h
#define SKSVirtualPlaneHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "common.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SKSVirtualPlaneHit : public G4VHit
{
 public:

   SKSVirtualPlaneHit();
  ~SKSVirtualPlaneHit();
   SKSVirtualPlaneHit(const SKSVirtualPlaneHit&);
   const SKSVirtualPlaneHit& operator=(const SKSVirtualPlaneHit&);
   G4int operator==(const SKSVirtualPlaneHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   void Draw();
   void Print();
      
 public:
   void SetDetectorID(G4int i) {detectorID = i;};
   G4int GetDetectorID() { return detectorID;};
   void SetHitPosition(double x, double y) {xpos = x; ypos = y;};
   void GetHitPosition(double *x, double *y) {*x = xpos; *y = ypos;};
 private:
   G4int detectorID;
   G4double xpos, ypos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<SKSVirtualPlaneHit> SKSVirtualPlaneHitsCollection;

extern G4Allocator<SKSVirtualPlaneHit> SKSVirtualPlaneHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* SKSVirtualPlaneHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) SKSVirtualPlaneHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void SKSVirtualPlaneHit::operator delete(void* aHit)
{
  SKSVirtualPlaneHitAllocator.FreeSingle((SKSVirtualPlaneHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
