/*
  PhysicsList.hh
  2007/4  K.Shirotori
*/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList
{
public:
  PhysicsList();
  ~PhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructHeavyIon();
  void ConstructStableMesons();
  void ConstructStableHyperons();
  void ConstructDecayingHyperons();
  void ConstructVirtualResonance();  
private:
  // these methods Construct physics processes and register them
  void AddTransportationSks();
  void ConstructDecay();
  void ConstructEM();
  //void ConstructHadronic();
};

#endif

 
