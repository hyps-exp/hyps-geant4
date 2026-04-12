/*
  MatrixTrigMan.hh

  2012/1/24
*/

#ifndef MatrixTrigMan_h
#define MatrixTrigMan_h 1

#include "ThreeVector.hh"
#include "DetectorID.hh"
#include <string>
#include <vector>
#include <map>

class DCGeomRecord;

class MatrixTrigMan 
{
private:
  MatrixTrigMan();
public:
  ~MatrixTrigMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; bool ret=Initialize(); return ret;}
  bool Initialize( const std::string &filename )
  { filename_=filename; bool ret=Initialize(); return ret;}

  bool InitializeSftTrig( void );
  bool InitializeSftTrig( const char *filename )
  { filenameSftTrig_=filename; bool ret=InitializeSftTrig(); return ret;}
  bool InitializeSftTrig( const std::string &filename )
  { filenameSftTrig_=filename; bool ret=InitializeSftTrig(); return ret;}

  bool InitializeMassTrig( void );
  bool InitializeMassTrig( const char *filename )
  { filenameMassTrig_=filename; bool ret=InitializeMassTrig(); return ret;}
  bool InitializeMassTrig( const std::string &filename )
  { filenameMassTrig_=filename; bool ret=InitializeMassTrig(); return ret;}

  static MatrixTrigMan & GetInstance( void );
  bool MatrixTrigger( int ch, int tof );
  bool SftTrigger( int ch, int tof, int sft );  
  bool MassTrigger( int ch, int tof, double time );
private:
  static MatrixTrigMan *mtMan_;
  std::string filename_;
  std::string filenameSftTrig_;
  std::string filenameMassTrig_;

  int TOF_Table_min_[NumOfSegCH];
  int TOF_Table_max_[NumOfSegCH];

  int SFT_Table_min_[NumOfSegCH][NumOfSegTOF];
  int SFT_Table_max_[NumOfSegCH][NumOfSegTOF];

  double MsT_Table_min_[NumOfSegCH][NumOfSegTOF];
  double MsT_Table_max_[NumOfSegCH][NumOfSegTOF];
private:
};


#endif
