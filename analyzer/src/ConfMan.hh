/*
  ConfMan.hh

  2012/1/24
*/

#ifndef Confman_h
#define Confman_h 1

#include <string>

class HodoParamMan;
class HodoPHCMan;

class ScalerAnalyzer_;
class ScalerAna;

class DCGeomMan;
class DCTdcCalibMan;
class DCDriftParamMan;
class K18TransMatrix;
class FieldMan;
class SimuData;
class MatrixTrigMan;

class EvDispCFT;

class ConfMan
{

private:
  //Conf
  std::string   m_ConfFileName;
  static ConfMan *confManager_;
  std::string m_ConfPath;

  //Hodo
  std::string   HodoParamFileName_;
  HodoParamMan *HodoParamManager_;
  std::string   HodoPHCFileName_;
  HodoPHCMan   *HodoPHCManager_;

  //DC
  std::string      DCGeomFileName_;
  DCGeomMan       *DCGeomManager_;
  std::string      DCTdcCalibFileName_;
  DCTdcCalibMan   *DCTdcCalibManager_;
  std::string      DCDriftParamFileName_;
  DCDriftParamMan *DCDriftParamManager_;
  std::string      K18MatrixFileName_;
  K18TransMatrix  *K18Matrix_;
  std::string      FieldMapFileName_;
  double           MapScale_;

  std::string      bh1FilterFileName_;
  std::string      bh2FilterFileName_;

  // Matrix trig
  std::string      mtTrigFileName_;
  std::string      sftTrigFileName_;
  std::string      massTrigFileName_;
  MatrixTrigMan   *mtTrigMan_;

  // CFTTrackingEff
  std::string      CftEffFileName_;

  // h_Utility
  std::string      FileManagerName_;

  double K18Momentum_;
  double SKSFieldNMR_;
  double SKSFieldCalc_;

  //Scaler
  std::string ScalerDefinitionFileName_;
  ScalerAna  *ScalerAnalyzer_;

  //MHTDC event suppression
  bool FlagMHTDC_;

  EvDispCFT *evDispCFT_;
  mutable bool FlagEvDispCFT_;

public:
  ~ConfMan();
  ConfMan(const std::string &confFile);
  //  static ConfMan& GetInstance();
  bool            Initialize();
  // bool RegisterParamFile(const char* target, std::string& dest);
  bool RegisterParamFile(const std::string& target, std::string& dest);

  //bool InitializeEvDisp( int RunNum );
  bool InitializeEvDispCFT( int RunNum );
  static ConfMan *GetConfManager( void ) { return confManager_; }

  //Hodo
  HodoParamMan *GetHodoParamManager( void ) { return HodoParamManager_;}
  HodoPHCMan   *GetHodoPHCManager( void ) { return HodoPHCManager_;}

  //DC
  DCGeomMan       *GetDCGeomManager( void ) { return DCGeomManager_; }
  DCTdcCalibMan   *GetDCTdcCalibManager( void ) { return DCTdcCalibManager_; }
  DCDriftParamMan *GetDCDriftParamManager( void ) { return DCDriftParamManager_; }
  K18TransMatrix  *GetK18Matrix( void ) { return K18Matrix_; }

  double MapScale( void ) const { return MapScale_; }
  double K18Momentum( void ) const { return K18Momentum_; }
  double SKSFieldNMR( void ) const { return SKSFieldNMR_; }
  double SKSFieldCalc( void ) const { return SKSFieldCalc_; }

  bool GetMHTDCFlag( void ) const { return FlagMHTDC_; }

  // Matrix Trigger
  MatrixTrigMan       *GetMatrixTrigManager( void ) { return mtTrigMan_; }


  bool GetEvDispCftFlag( void ) const { return FlagEvDispCFT_; }
  bool SetEvDispCftFlag( bool flag ) const { FlagEvDispCFT_=flag; return true;}
  EvDispCFT *GetEvDispCFT( void ) { return evDispCFT_; }

private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator=(const ConfMan&);
  bool ResolveConfDirPath( void );

public:
  //bool InitializeHistograms();
  bool InitializeParameterFiles();
};

#endif
