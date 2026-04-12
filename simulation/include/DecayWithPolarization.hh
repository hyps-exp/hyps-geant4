#ifndef __DECAY_WITH_POLARIZATION_HH__
#define __DECAY_WITH_POLARIZATION_HH__

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

class DecayWithPolarization {
private:
  double M1;     // Mother particle

// morther particle direction z axis
// morther particle polarization y axis
  double Pol;    // Polarization
  double Asym; // Asymmetry parameter
  G4ThreeVector PolVec;  // Polarization vector (y direction) 
  double CTau; // ctau (mm)
  bool   ChargeFlag; 

  double M2, M3; // Daughter particle
  double PCM;
  double ThetaCM;
  double PhiCM;
  G4ThreeVector P2vec_cm;
  G4ThreeVector P3vec_cm;
  G4LorentzVector Lv2;
  G4LorentzVector Lv3;

  G4ThreeVector P1Vec;  
  G4ThreeVector PosVec;  // production position
  G4ThreeVector DecayPosVec;  // production position

  G4ThreeVector P2vec_lab;
  G4ThreeVector P3vec_lab;
  double Ekin2;
  double Ekin3;
public:
  DecayWithPolarization(double m1, double m2, double m3, 
			double pol, double asym, G4ThreeVector polVec,
			G4ThreeVector p1vec, G4ThreeVector posvec, 
			bool phiFlag, double ctau, bool chargeFlag);
  double RandAngularDistribution();
  void calcThetaPhi(G4ThreeVector vec, double *theta, double *phi);
  void DecayCheck();
  G4bool decayCheck(G4double ctau, G4double momentum, G4double mass, G4double dx);
  int  getTargetFlag(G4ThreeVector pos);
  G4double calcEnergyDeposit(G4double momentum, G4double mass, G4double dx, G4ThreeVector pos, G4int flagTgtType);
  G4double calc_dE_dx(double beta, int flagTgtType);
  G4ThreeVector GetMomentum2() {return P2vec_lab;}
  G4ThreeVector GetMomentum3() {return P3vec_lab;}
  double        GetEkin2()     {return Ekin2;}
  double        GetEkin3()     {return Ekin3;}
  G4ThreeVector GetDecayPos()  {return DecayPosVec;}

};
#endif
