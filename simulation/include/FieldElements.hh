/*
  FieldElements.hh

  2004/7/5  T.Takahashi

*/

#ifndef FieldElements_h 

#define FieldElements_h

#include "DCGeomRecord.hh"

enum FldElemReg { FERSurface=0, FERInside, FEROutside };

class FieldElements 
{
public:
  FieldElements( const char *name, const G4ThreeVector &pos,
                 double ta, double ra1, double ra2 );
  virtual ~FieldElements() {}

private:
  DCGeomRecord record_;

public:
  G4ThreeVector Local2GlobalPos( const G4ThreeVector &in ) const;
  G4ThreeVector Local2GlobalDir( const G4ThreeVector &in ) const;
  G4ThreeVector Global2LocalPos( const G4ThreeVector &in ) const;
  G4ThreeVector Global2LocalDir( const G4ThreeVector &in ) const;

  virtual G4ThreeVector GetField( const G4ThreeVector &gPos ) const = 0;
  virtual bool ExistField( const G4ThreeVector &gPos ) const = 0;

  virtual FldElemReg checkRegion( const G4ThreeVector &gPos, 
                                  double Tolerance ) const = 0;
};



#endif 
