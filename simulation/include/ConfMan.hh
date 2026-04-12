/*
  ConfMan.hh
  2007/4  K.Shirotori
*/

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>

class DCGeomMan;
//class EvDisp;

class ConfMan
{
public:
  explicit ConfMan( const std::string & filename );
  ~ConfMan();
private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );

private:
  std::string ConfFileName_;
  static ConfMan *confManager_;
  std::string ConfPath_;

private:
  std::string DCGeomFileName_;
  int SksJMode_;
  DCGeomMan *DCGeomManager_;

  std::string FieldMapFileName_;
  double K18Momentum_;
  int GeomFlag_;
  int FieldMapFlag_;
  double MapScale_;
  double Charge_;
  int TgtFlag_;
  double TgtLength_;
  int SKSMaterialFlag_;
  int PhysFlag_;
  int StepFlag_;
  int BeamFlag_;
  int ReactionMode_;
  double PrimaryE_;
  double GenerateThetaCM_;

public:
  bool Initialize( void );
  static ConfMan *GetConfManager( void ) { return confManager_; }
  DCGeomMan *GetDCGeomManager( void ) { return DCGeomManager_; }
  std::string FieldMapName(void) { return FieldMapFileName_; }
  double K18Momentum( void ) const { return K18Momentum_; }
  int GeomFlag( void ) const { return GeomFlag_; }
  int FieldMapFlag( void ) const { return FieldMapFlag_; }
  double MapScale( void ) const { return MapScale_; }
  double Charge( void ) const { return Charge_; }
  double TgtLength( void ) const { return TgtLength_; }
  int SKSMaterialFlag( void ) const { return SKSMaterialFlag_; }
  int Target( void ) const { return TgtFlag_; }
  int TgtFlag( void ) const { return TgtFlag_; }
  int PhysFlag( void ) const { return PhysFlag_; }
  int StepFlag( void ) const { return StepFlag_; }
  int BeamFlag( void ) const { return BeamFlag_; }
  int ReactionMode( void ) const { return ReactionMode_; }
  double GetPrimaryE( void ) const { return PrimaryE_; }
  double GetGenerateThetaCM( void ) const { return GenerateThetaCM_; }
  bool RegisterParamFile(const char* target, std::string& dest);
private:
  bool InitializeParameterFiles( void );
  bool EndAnalysis( void );
  bool ResolveConfDirPath( void );

};

#endif
