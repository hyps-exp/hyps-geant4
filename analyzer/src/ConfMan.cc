/*
  ConfMan.cc

  2012/1/24
*/

#include "ConfMan.hh"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
//#include <boost/lexical_cast.hpp>
#include "lexical_cast.hh"

#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "K18TransMatrix.hh"

#include "EvDispCFT.hh"

ConfMan * ConfMan::confManager_=0;

ConfMan::ConfMan()
  : m_ConfFileName(),
    DCGeomManager_(0),
    DCTdcCalibManager_(0),
    DCDriftParamManager_(0),
    K18Matrix_(),
    K18Momentum_(1.05),
    SKSFieldNMR_(2.2), SKSFieldCalc_(2.2),
    FlagMHTDC_(false),
    evDispCFT_(0), FlagEvDispCFT_(false), MapScale_(0.7)
{
  std::cout << "MapScale00 : " << MapScale_ << std::endl;
}

ConfMan::ConfMan(const std::string& confFile)
  :m_ConfFileName(confFile)
{
  if (confManager_)
    std::exit(-1);

  confManager_ = this;
}

ConfMan::~ConfMan()
{
}

// ConfMan&
// ConfMan::GetInstance()
// {
//   static ConfMan gConfMan;
//   return gConfMan;
// }

//______________________________________________________________________________
bool ConfMan::ResolveConfDirPath( void )
{
  try {
    std::filesystem::path conf_path(m_ConfFileName);
    std::filesystem::path conf_dir = conf_path.parent_path();

    if (std::filesystem::exists(conf_dir)) {
      m_ConfPath = std::filesystem::canonical(conf_dir).string();
    } else {
      m_ConfPath = std::filesystem::absolute(conf_dir).string();
    }
    return true;

  } catch (const std::filesystem::filesystem_error& e) {
    std::cerr << "#E: [ConfMan::" << __func__ << "] " << e.what() << std::endl;
    return false;
  }
}

// //______________________________________________________________________________
// bool ConfMan::RegisterParamFile( const char* target, std::string& dest )
// {
//   std::string buf_str(target);

//   if (!buf_str.empty()
//       && (buf_str.back() == '\n' || buf_str.back() == '\r'))
//     {
//       buf_str.erase(buf_str.size() - 1);
//     }

//   if (!target || target[0] == '\0') {
//     return false;
//   }

//   try {
//     // std::string tgt_str(target);
//     std::string tgt_str = buf_str;
//     tgt_str.erase(tgt_str.find_last_not_of("\n\r") + 1);

//     std::filesystem::path base(ConfPath_);
//     std::filesystem::path combined = base / tgt_str;

//     dest = std::filesystem::weakly_canonical(combined).string();

//     return true;
//   } catch (const std::exception& e) {
//     std::cerr << "Path error: " << e.what() << std::endl;
//     return false;
//   }
// }

//______________________________________________________________________________
bool ConfMan::RegisterParamFile( const std::string& target, std::string& dest )
{
  if (target.empty()) {
    return false;
  }

  try {
    std::filesystem::path base(m_ConfPath);
    std::filesystem::path combined = base / target;

    dest = std::filesystem::weakly_canonical(combined).string();

    return true;
  } catch (const std::exception& e) {
    std::cerr << "Path error: " << e.what() << std::endl;
    return false;
  }
}

bool
ConfMan::Initialize()
{
  bool status = true;

  double val;
  std::ifstream f(m_ConfFileName.c_str());

  if (f.fail()){
    std::cerr << "#E "
	      << std::endl;
    return false;
  }

  if (!ResolveConfDirPath()) return false;

  MapScale_ = 0.7;
  std::cout << "MapScale0 : " << MapScale_ << std::endl;

  while (f.good()){
    std::string line;
    std::getline(f, line);
    if (line.empty())
      continue;
    std::istringstream input_line(line);
    std::istream_iterator<std::string> line_begin(input_line);
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> param_list(line_begin, line_end);
    if (param_list.size()!=2)
      continue;
    std::cout << param_list[0] << " " << param_list[1] << std::endl;
    //Conf
    //Hodo
    if (param_list[0] == "HDPRM:")
      // HodoParamFileName_ = param_list[1];
      RegisterParamFile(param_list[1], HodoParamFileName_);
    if (param_list[0] == "HDPHC:")
      // HodoPHCFileName_ = param_list[1];
      RegisterParamFile(param_list[1], HodoPHCFileName_);
    //DC
    if (param_list[0] == "DCGEO:")
      // DCGeomFileName_ = param_list[1];
      RegisterParamFile(param_list[1], DCGeomFileName_);
    if (param_list[0] == "DCTDC:")
      // DCTdcCalibFileName_ = param_list[1];
      RegisterParamFile(param_list[1], DCTdcCalibFileName_);
    if (param_list[0] == "DCDRFT:")
      // DCDriftParamFileName_ = param_list[1];
      RegisterParamFile(param_list[1], DCDriftParamFileName_);
    if (param_list[0] == "K18TM:")
      // K18MatrixFileName_ = param_list[1];
      RegisterParamFile(param_list[1], K18MatrixFileName_);
    if (param_list[0] == "PK18:")
      //K18Momentum_ = boost::lexical_cast<double>(param_list[1]);
      K18Momentum_ = hddaq::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDMAP:")
      // FieldMapFileName_ = param_list[1];
      RegisterParamFile(param_list[1], FieldMapFileName_);
    if (param_list[0] == "MAP_SCALE:")
      //MapScale_ = boost::lexical_cast<double>(param_list[1]);
      MapScale_ = hddaq::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDNMR:")
      //SKSFieldNMR_ = boost::lexical_cast<double>(param_list[1]);
      SKSFieldNMR_ = hddaq::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "FLDCALC:")
      //SKSFieldCalc_ = boost::lexical_cast<double>(param_list[1]);
      SKSFieldCalc_ = hddaq::lexical_cast<double>(param_list[1]);
    if (param_list[0] == "BH1FLT:")
      // bh1FilterFileName_ = param_list[1];
      RegisterParamFile(param_list[1], bh1FilterFileName_);
    if (param_list[0] == "BH2FLT:")
      // bh2FilterFileName_ = param_list[1];
      RegisterParamFile(param_list[1], bh2FilterFileName_);

    if (param_list[0] == "MTTRIG:")
      // mtTrigFileName_ = param_list[1];
      RegisterParamFile(param_list[1], mtTrigFileName_);
    if (param_list[0] == "SFTTRIG:")
      // sftTrigFileName_ = param_list[1];
      RegisterParamFile(param_list[1], sftTrigFileName_);
    if (param_list[0] == "MASSTRIG:")
      // massTrigFileName_ = param_list[1];
      RegisterParamFile(param_list[1], massTrigFileName_);
    if (param_list[0] == "CFTEFF:")
      // CftEffFileName_ = param_list[1];
      RegisterParamFile(param_list[1], CftEffFileName_);

    //MHTDC event suppression flag
    if( param_list[0] == "MHTDC:" )
      if(param_list[1]=="1") FlagMHTDC_=true;
      else         FlagMHTDC_=false;

    if( param_list[0] == "EVDISP_CFT:" )
      if(param_list[1]=="1") FlagEvDispCFT_=true;
      else         FlagEvDispCFT_=false;

  }

  std::cout << "DCGEO  file   = " << DCGeomFileName_ << std::endl;
  std::cout << "FLDMAP file   = " << FieldMapFileName_ << std::endl;
  std::cout << "MTTRIG file   = " << mtTrigFileName_ << std::endl;
  std::cout << "SFTTRIG file  = " << sftTrigFileName_ << std::endl;
  std::cout << "MASSTRIG file = " << massTrigFileName_ << std::endl;
  std::cout << "CFTEFF file   = " << CftEffFileName_ << std::endl;

  status = InitializeParameterFiles();

  return status;
}

bool ConfMan::InitializeEvDispCFT( int RunNum )
{
  static const std::string funcname = "[ConfMan::InitializeEvDispCFT]";
  evDispCFT_ = & EvDispCFT::GetInstance();
  evDispCFT_->Initialize(RunNum);

  return true;
}
