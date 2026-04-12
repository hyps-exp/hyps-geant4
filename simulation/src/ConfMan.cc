
/*
  ConfMan.cc
  2007/4  K.Shirotori
*/

#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <filesystem>

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),DCGeomManager_(0), K18Momentum_(1.5),
    GeomFlag_(0), FieldMapFlag_(0), MapScale_(0), Charge_(1),
    SKSMaterialFlag_(0), TgtFlag_(0), TgtLength_(0.0),
    PhysFlag_(0), StepFlag_(0), ReactionMode_(0), BeamFlag_(0),
    PrimaryE_(0.), GenerateThetaCM_(0.)
{
  static const std::string funcname = "[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constructing twice" << std::endl;
    exit(-1);
  }
  confManager_ = this;
}

ConfMan::~ConfMan()
{
  EndAnalysis();
  confManager_=0;
}

bool ConfMan::EndAnalysis( void )
{

  if(DCGeomManager_){
    delete DCGeomManager_; DCGeomManager_=0;
  }

  return true;
}

//______________________________________________________________________________
bool ConfMan::ResolveConfDirPath( void )
{
  try {
    std::filesystem::path conf_path(ConfFileName_);
    std::filesystem::path conf_dir = conf_path.parent_path();

    if (std::filesystem::exists(conf_dir)) {
      ConfPath_ = std::filesystem::canonical(conf_dir).string();
    } else {
      ConfPath_ = std::filesystem::absolute(conf_dir).string();
    }
    return true;

  } catch (const std::filesystem::filesystem_error& e) {
    std::cerr << "#E: [ConfMan::" << __func__ << "] " << e.what() << std::endl;
    return false;
  }
}

//______________________________________________________________________________
bool ConfMan::RegisterParamFile( const char* target, std::string& dest )
{
  std::string buf_str(target);

  if (!buf_str.empty()
      && (buf_str.back() == '\n' || buf_str.back() == '\r'))
    {
      buf_str.erase(buf_str.size() - 1);
    }

  if (!target || target[0] == '\0') {
    return false;
  }

  try {
    std::string tgt_str(target);
    tgt_str.erase(tgt_str.find_last_not_of("\n\r") + 1);

    std::filesystem::path base(ConfPath_);
    std::filesystem::path combined = base / tgt_str;

    dest = std::filesystem::weakly_canonical(combined).string();

    return true;
  } catch (const std::exception& e) {
    std::cerr << "Path error: " << e.what() << std::endl;
    return false;
  }
}


const int BufSize = 144;

bool ConfMan::Initialize( void )
{
  static const std::string funcname = "[ConfMan::Initialize]";

  FILE *fp;
  char buf[BufSize], buf1[BufSize];
  int    intval;
  double val;

  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  if (!ResolveConfDirPath()) return false;

  // read configuration file
  while( fgets(buf,BufSize,fp)!=0 ){
    if( buf[0]!='#' ){
      //Geometry
      if( sscanf(buf,"DCGEO: %s",buf1)==1 )
	// DCGeomFileName_=buf1;
	RegisterParamFile(buf1, DCGeomFileName_);
      if( sscanf(buf,"FLDMAP: %s",buf1)==1 )
	// FieldMapFileName_=buf1;
	RegisterParamFile(buf1, FieldMapFileName_);
      else if( sscanf(buf,"CHECKMAP: %lf", &val )==1 )
	FieldMapFlag_ = (int)val;
      else if( sscanf(buf,"MAP_SCALE: %lf", &val )==1 )
	MapScale_=val;
      else if( sscanf(buf,"CHARGE: %lf", &val )==1 )
	Charge_=val;
      //Target
      else if( sscanf(buf,"TGTMATER: %lf", &val )==1 )
	TgtFlag_ = (int)val;
      else if( sscanf(buf,"TGTLENGTH: %lf", &val )==1 )
	TgtLength_=val;
      //Material
      else if( sscanf(buf,"MFLAGSKS: %lf", &val )==1 )
	SKSMaterialFlag_ = (int)val;
      //Phys process
      else if( sscanf(buf,"EM: %lf",  &val )==1 )
	PhysFlag_ = ( PhysFlag_ | ( (int)val<<0 ));
      else if( sscanf(buf,"DECAY: %lf",   &val )==1 )
	PhysFlag_ = ( PhysFlag_ | ( (int)val<<1 ));
      else if( sscanf(buf,"HADRON: %lf", &val )==1 )
	PhysFlag_ = ( PhysFlag_ | ( (int)val<<2 ));
      //Stopping action
      else if( sscanf(buf,"SKSSTOP: %lf",  &val )==1 )
	StepFlag_ = ( StepFlag_ | ( (int)val<<0 ));
      else if( sscanf(buf,"SKSGSTOP: %lf",   &val )==1 )
	StepFlag_ = ( StepFlag_ | ( (int)val<<1 ));
      else if( sscanf(buf,"NUSTOP: %lf", &val )==1 )
	StepFlag_ = ( StepFlag_ | ( (int)val<<2 ));
      else if( sscanf(buf,"GSTOP: %lf", &val )==1 )
	StepFlag_ = ( StepFlag_ | ( (int)val<<3 ));
      else if( sscanf(buf,"ESTOP: %lf", &val )==1 )
	StepFlag_ = ( StepFlag_ | ( (int)val<<4 ));
      //Primary generator
      else if( sscanf(buf,"PK18: %lf", &val )==1 )
	K18Momentum_=val;
      else if( sscanf(buf,"PK18WIDE: %lf",  &val )==1 )
	BeamFlag_ = ( BeamFlag_ | ( (int)val<<0 ));
      else if( sscanf(buf,"STBEAM: %lf",  &val )==1 )
	BeamFlag_ = ( BeamFlag_ | ( (int)val<<1 ));
      else if( sscanf(buf, "REACTION: %lf", &val )==1)
	ReactionMode_=(int)val;
      else if( sscanf(buf, "PRIME: %lf", &val )==1)
	PrimaryE_=val;
      else if( sscanf(buf, "THETACM: %lf", &val )==1)
	GenerateThetaCM_=val;
    } /* if( buf[0]... ) */
  } /* while(...) */

  fclose(fp);

  std::cout << "------------------"  << ConfFileName_ << "------" << std::endl;
  std::cout << "**********Geometry**********" << std::endl;
  std::cout << "DC Geom. Param.:  "  << DCGeomFileName_           << std::endl;
  std::cout << "SKS Material:     "      << SKSMaterialFlag_      << std::endl;
  std::cout << "**********Field Map**********" << std::endl;
  std::cout << "FieldMap:         "  << FieldMapFileName_         << std::endl;
  std::cout << "Map Scale:        "  << MapScale_                 << std::endl;
  std::cout << "Charge:           "  << Charge_                   << std::endl;
  std::cout << "Check Map:        "  << FieldMapFlag_             << std::endl;
  std::cout << "**********Target**********" << std::endl;
  std::cout << "Tgt Material:     "      << TgtFlag_      << std::endl;
  std::cout << "Target Length:    "  << TgtLength_                << std::endl;
  std::cout << "**********Physics process**********" << std::endl;
  std::cout << "EM process:       "  << (( PhysFlag_ >> 0 )&0x1)  << std::endl;
  std::cout << "Decay process:    "  << (( PhysFlag_ >> 1 )&0x1)  << std::endl;
  std::cout << "Hadronic process: "  << (( PhysFlag_ >> 2 )&0x1)  << std::endl;
  std::cout << "**********Stepping Action**********" << std::endl;
  std::cout << "SKS stop:           "  << (( StepFlag_ >> 0 )&0x1)  << std::endl;
  std::cout << "SKS gamma-ray stop: "  << (( StepFlag_ >> 1 )&0x1)  << std::endl;
  std::cout << "Nuetrino stop:      "  << (( StepFlag_ >> 2 )&0x1)  << std::endl;
  std::cout << "Gamma stop:         "  << (( StepFlag_ >> 3 )&0x1)  << std::endl;
  std::cout << "Electron stop:      "  << (( StepFlag_ >> 4 )&0x1)  << std::endl;
  std::cout << "**********Primary Action**********" << std::endl;
  std::cout << "K18 Mom. [GeV/c]:    "  << K18Momentum_              << std::endl;
  std::cout << "Mom 1.0-1.9 [GeV/c]: "  << (( BeamFlag_ >> 0 )&0x1)  << std::endl;
  std::cout << "Straight beam:       "  << (( BeamFlag_ >> 1 )&0x1)  << std::endl;
  std::cout << "Reaction:            "  << ReactionMode_             << std::endl;
  std::cout << "Primary Energy:      "  << PrimaryE_                 << std::endl;
  std::cout << "Generate ThetaCM:    "  << GenerateThetaCM_          << std::endl;

  std::cout << "-----------------------------------------------" << std::endl;

  InitializeParameterFiles();

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{

  DCGeomManager_ = & DCGeomMan::GetInstance();
  if( DCGeomFileName_!="" )
    DCGeomManager_->Initialize(DCGeomFileName_);
  else
    DCGeomManager_->Initialize();

  return true;
}
