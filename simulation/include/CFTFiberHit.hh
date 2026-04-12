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

#ifndef CFTFiberHit_h
#define CFTFiberHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CFTFiberHit : public G4VHit
{
 public:

   CFTFiberHit();
  ~CFTFiberHit();
   CFTFiberHit(const CFTFiberHit&);
   const CFTFiberHit& operator=(const CFTFiberHit&);
   G4int operator==(const CFTFiberHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   void Draw();
   void Print();
      
 public:
   void SetLayerID(G4int i) {layerID = i;};
   void SetSegmentID(G4int i) {segmentID = i;};
   void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
   void SetTime(G4double time) {Time =time;};
   void SetParticleID(G4int i) {particleId = i;};
   void SetParentID(G4int i) {parentId = i;};
   void SetParticleName(G4String name) {particleName = name;};

   G4int GetLayerID() { return layerID;};
   G4int GetSegmentID() { return segmentID;};
   G4double GetEdepAbs()     { return EdepAbs; };
   G4double GetTrakAbs()     { return TrackLengthAbs; };
   G4double GetTime()     { return Time; };
   G4int GetParticleID() {return particleId;};
   G4int GetParentID()   {return parentId;};
   G4String GetParticleName() {return particleName;};

 private:
   G4int layerID;
   G4int segmentID;
   G4double EdepAbs, TrackLengthAbs;
   G4double Time;
   G4String particleName;
   G4int    particleId;
   G4int    parentId;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<CFTFiberHit> CFTFiberHitsCollection;

extern G4Allocator<CFTFiberHit> CFTFiberHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CFTFiberHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) CFTFiberHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void CFTFiberHit::operator delete(void* aHit)
{
  CFTFiberHitAllocator.FreeSingle((CFTFiberHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


