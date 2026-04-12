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

#ifndef CrystalHit_h
#define CrystalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CrystalHit : public G4VHit
{
 public:

   CrystalHit();
  ~CrystalHit();
   CrystalHit(const CrystalHit&);
   const CrystalHit& operator=(const CrystalHit&);
   G4int operator==(const CrystalHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   void Draw();
   void Print();
      
 public:
  
   void SetLayerID(G4int i) {layerID = i;};
   void SetChannelID(G4int i) {channelID = i;};
   void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
   void SetTime(G4double time) {Time =time;};
   void SetParticleName(G4String name) {Name = name;};
   void SetPos(G4double x, G4double y, G4double z) {Position[0]=x;Position[1]=y;Position[2]=z;};
   G4int GetLayerID() { return layerID;};
   G4int GetChannelID() { return channelID;};
   G4double GetEdepAbs()     { return EdepAbs; };
   G4double GetTrakAbs()     { return TrackLengthAbs; };
   G4double GetTime()     { return Time; };
   G4String GetParticleName() { return Name;};
   void GetPos(G4double *pos) 
  { pos[0]=Position[0]; pos[1]=Position[1]; pos[2]=Position[2];};

 private:
   G4int layerID;
   G4int channelID;
   G4double EdepAbs, TrackLengthAbs;
   G4double Time;
   G4double Position[3];
   G4String Name;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<CrystalHit> CrystalHitsCollection;

extern G4Allocator<CrystalHit> CrystalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CrystalHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) CrystalHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void CrystalHit::operator delete(void* aHit)
{
  CrystalHitAllocator.FreeSingle((CrystalHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


