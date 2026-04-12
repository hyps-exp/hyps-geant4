/*
  MatrixTrigMan.cc

  2012/1/24
*/

#include "MatrixTrigMan.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>

const int MaxChar = 200;

MatrixTrigMan *MatrixTrigMan::mtMan_=0;

MatrixTrigMan::MatrixTrigMan()
{
  for (int i=0; i<NumOfSegCH; i++) {
    TOF_Table_min_[i] = -1;
    TOF_Table_max_[i] = -1;
    for (int j=0; j<NumOfSegTOF; j++) {
      SFT_Table_min_[i][j] = -1;
      SFT_Table_max_[i][j] = -1;

      MsT_Table_min_[i][j] = -1.;
      MsT_Table_max_[i][j] = -1.;
    }
  }

}

MatrixTrigMan::~MatrixTrigMan()
{}


MatrixTrigMan & MatrixTrigMan::GetInstance( void )
{
  if( !mtMan_ ){
    mtMan_ = new MatrixTrigMan();
  }
  return *mtMan_;
}


bool MatrixTrigMan::Initialize( void )
{
  static const std::string funcname = "[MatrixTrigMan::Initialize]";
  char str[MaxChar];
  char cname[MaxChar];
  int ch_id, tof_min, tof_max;

  FILE *fp;

  if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %d %d",
		  &ch_id, &tof_min, &tof_max) == 3 ){
	if (ch_id>=0 && ch_id<NumOfSegCH) {
	  TOF_Table_min_[ch_id] = tof_min;
	  TOF_Table_max_[ch_id] = tof_max;
	} else {
	  std::cerr << funcname << ": Invalid CH ID " << ch_id << std::endl;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << funcname << " Initialization finished." << std::endl;

  return true;
}

bool MatrixTrigMan::InitializeSftTrig( void )
{
  static const std::string funcname = "[MatrixTrigMan::InitializeSftTrig]";
  char str[MaxChar];
  char cname[MaxChar];
  int ch_id, tof_id, sft_min, sft_max;

  FILE *fp;

  if( ( fp = fopen( filenameSftTrig_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %d %d %d",
		  &ch_id, &tof_id, &sft_min, &sft_max) == 4 ){
	if (ch_id>=0 && ch_id<NumOfSegCH && tof_id>=0 && tof_id <NumOfSegTOF) {
	  SFT_Table_min_[ch_id][tof_id] = sft_min;
	  SFT_Table_max_[ch_id][tof_id] = sft_max;
	} else {
	  std::cerr << funcname << ": Invalid CH ID " << ch_id  << ", or TOF ID " << tof_id << std::endl;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << funcname << " Initialization finished." << std::endl;

  return true;
}

bool MatrixTrigMan::InitializeMassTrig( void )
{
  static const std::string funcname = "[MatrixTrigMan::InitializeMassTrig]";
  char str[MaxChar];
  char cname[MaxChar];
  int ch_id, tof_id;
  double time_min, time_max;

  FILE *fp;

  if( ( fp = fopen( filenameMassTrig_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %d %lf %lf",
		  &ch_id, &tof_id, &time_min, &time_max) == 4 ){
	if (ch_id>=0 && ch_id<NumOfSegCH && tof_id>=0 && tof_id <NumOfSegTOF) {
	  MsT_Table_min_[ch_id][tof_id] = time_min;
	  MsT_Table_max_[ch_id][tof_id] = time_max;
	} else {
	  std::cerr << funcname << ": Invalid CH ID " << ch_id  << ", or TOF ID " << tof_id << std::endl;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << funcname << " Initialization finished." << std::endl;

  return true;
}


bool MatrixTrigMan::MatrixTrigger( int ch, int tof )
{
  static const std::string funcname = "[MatrixTrigMan::MatrixTrigger]";

  if ( !(ch>=0 && ch<NumOfSegCH) ) {
    std::cerr << funcname << ": Invalid CH ID " << ch << std::endl;
    return false;
  }

  if (tof>=TOF_Table_min_[ch] && tof<=TOF_Table_max_[ch])
    return true;
  else
    return false;

}

bool MatrixTrigMan::SftTrigger( int ch, int tof, int sft )
{
  static const std::string funcname = "[MatrixTrigMan::SftTrigger]";

  if ( !(ch>=0 && ch<NumOfSegCH) ) {
    std::cerr << funcname << ": Invalid CH ID " << ch << std::endl;
    return false;
  }
  if ( !(tof>=0 && tof<NumOfSegTOF) ) {
    std::cerr << funcname << ": Invalid TOF ID " << ch << std::endl;
    return false;
  }

  if (sft>=SFT_Table_min_[ch][tof] && sft<=SFT_Table_max_[ch][tof])
    return true;
  else
    return false;

}

bool MatrixTrigMan::MassTrigger( int ch, int tof, double time )
{
  static const std::string funcname = "[MatrixTrigMan::MassTrigger]";

  if ( !(ch>=0 && ch<NumOfSegCH) ) {
    std::cerr << funcname << ": Invalid CH ID " << ch << std::endl;
    return false;
  }
  if ( !(tof>=0 && tof<NumOfSegTOF) ) {
    std::cerr << funcname << ": Invalid TOF ID " << ch << std::endl;
    return false;
  }

  if (time>=MsT_Table_min_[ch][tof] && time<=MsT_Table_max_[ch][tof])
    return true;
  else
    return false;

}
