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

#ifndef SKSCounterHit_h
#define SKSCounterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "common.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SKSCounterHit : public G4VHit
{
 public:

   SKSCounterHit();
  ~SKSCounterHit();
   SKSCounterHit(const SKSCounterHit&);
   const SKSCounterHit& operator=(const SKSCounterHit&);
   G4int operator==(const SKSCounterHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   void Draw();
   void Print();
      
 public:
  
   void SetDetectorID(G4int i) {detectorID = i;};
   void SetChannelID(G4int i) {channelID = i;};
   void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
   void SetTime(G4double time) {Time =time;};
   void SetPos(G4double x, G4double y, G4double z) {Position[XCOORD]=x;Position[YCOORD]=y;Position[ZCOORD]=z;};
   void SetLcCherenkov(G4int hit) {lcCherenkov =hit;};
   G4int GetDetectorID() { return detectorID;};
   G4int GetlcCherenkov() { return lcCherenkov;};
   G4int GetChannelID() { return channelID;};
   G4double GetEdepAbs()     { return EdepAbs; };
   G4double GetTrakAbs()     { return TrackLengthAbs; };
   G4double GetTime()     { return Time; };
   void GetPos(G4double *pos) 
  { pos[XCOORD]=Position[XCOORD]; pos[YCOORD]=Position[YCOORD]; pos[ZCOORD]=Position[ZCOORD];};
   G4int GetLcCherenkov() {return lcCherenkov;};
 private:
   G4int detectorID;
   G4int channelID;
   G4double EdepAbs, TrackLengthAbs;
   G4double Time;
   G4double Position[XYZ];
   G4int lcCherenkov;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<SKSCounterHit> SKSCounterHitsCollection;

extern G4Allocator<SKSCounterHit> SKSCounterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* SKSCounterHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) SKSCounterHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void SKSCounterHit::operator delete(void* aHit)
{
  SKSCounterHitAllocator.FreeSingle((SKSCounterHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


