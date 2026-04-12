/*
  BH2Filter.cc

  2012/3/12
*/

// -*- C++ -*-

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iterator>

#include "DetectorID.hh"
#include "HodoAnalyzer.hh"
#include "DCHit.hh" 
#include "DCAnalyzer.hh"
#include "BH2Hit.hh"
#include "HodoCluster.hh"
#include "BH2Filter.hh"

//______________________________________________________________________________
// struct BH2Filter::Param
//______________________________________________________________________________
BH2Filter::
Param::Param()
  : m_xmin(NumOfLayersBcOut),
    m_xmax(NumOfLayersBcOut)
{
}

//______________________________________________________________________________
BH2Filter::
Param::~Param()
{
}

//______________________________________________________________________________
void
BH2Filter::
Param::Print(const std::string& arg) const
{

  std::cout << "\n " << arg << " (xmin, xmax) = \n";
  for (int i=0, n=m_xmin.size(); i<n; ++i)
    {
      std::cout << " iplane " << std::setw(3) << i 
		<< " (" << std::setw(6) << m_xmin[i] << " "
		<< std::setw(6) << m_xmax[i] << ")\n";
    }
  std::cout << std::endl;
  return;
}

//______________________________________________________________________________
// class BH2Filter
//______________________________________________________________________________
BH2Filter::BH2Filter()
  : m_param(NumOfSegBH2)
{
}

//______________________________________________________________________________
BH2Filter::~BH2Filter()
{
}

//______________________________________________________________________________
void
BH2Filter::Apply(const HodoAnalyzer& hodo,
		 const DCAnalyzer& dc,
		 std::vector<std::vector<DCHitContainer> >& candidates)
{
  // +++++++++++++++++++++++++++++++++++++++++++
  // candidates [segment id] [plane id] [hit id]
  // +++++++++++++++++++++++++++++++++++++++++++

  m_dc   = &dc;
  m_hodo = &hodo;

  std::set<int> seg;
  for (int i=0, n=hodo.GetNHitsBH2(); i<n; ++i)
    {
      const BH2Hit* const h = hodo.GetHitBH2(i);
      if (!h)
	continue;
      seg.insert(h->SegmentId());
    }

  BuildCandidates(seg, candidates);

  return;
}

//______________________________________________________________________________
void BH2Filter::Apply_withT0Seg(int T0Seg,
		      const DCAnalyzer& dc,
		      std::vector<std::vector<DCHitContainer> >& candidates)
{
  // +++++++++++++++++++++++++++++++++++++++++++
  // candidates [segment id] [plane id] [hit id]
  // +++++++++++++++++++++++++++++++++++++++++++

  m_dc   = &dc;

  std::set<int> seg;
  seg.insert(T0Seg);

  BuildCandidates(seg, candidates);

  return;
}

//______________________________________________________________________________
void BH2Filter::BuildCandidates(std::set<int> & seg,
				std::vector<std::vector<DCHitContainer> >& candidates)
{
  candidates.resize(seg.size());
  std::vector<std::vector<DCHitContainer> >::iterator itCont = candidates.begin();
  for (std::set<int>::const_iterator itSeg = seg.begin(), itSegEnd = seg.end();
       itSeg!=itSegEnd; ++itSeg, ++itCont)
    {
      const int iSeg = *itSeg;
      std::vector<DCHitContainer>& c = *itCont;
      c.resize(NumOfLayersBcOut+1);

//       std::cout << "  BH2 seg = " << iSeg << "\n";
      for (int iplane=0; iplane<NumOfLayersBcOut; ++iplane)
	{
	  int iLayer = iplane + 1;
	  DCHitContainer& after = c[iLayer];
	  const double xmin = m_param[iSeg].m_xmin[iplane];
	  const double xmax = m_param[iSeg].m_xmax[iplane];
	  const DCHitContainer& before = m_dc->GetBcOutHC(iLayer);
	  for (int ih=0, nh=before.size(); ih<nh; ++ih)
	    {
	      const DCHit* const h = before[ih];
	      if (!h)
		continue;
	      const double wpos = h->GetWirePosition();
	      const int layer   = h->GetLayer();

// 	      std::cout << " layer = " << iplane 
// 			<< "(" << layer << ") : " << wpos
// 			<< " (" << xmin << ", " << xmax << ")";
	      if (wpos<xmin || xmax<wpos)
		{
// 		  std::cout << std::endl;
		  continue;
		}
	      // 	      std::cout << " good " << std::endl;
	      after.push_back(const_cast<DCHit*>(h));
	    }
// 	  std::cout << __FILE__ << ":" << __LINE__
// 		    << " " << after.size() << std::endl;
	}
    }
//   std::cout << __FILE__ << ":" << __LINE__ 
// 	    << " " << candidates.size() << std::endl;
  
  return;
}

//______________________________________________________________________________
BH2Filter&
BH2Filter::GetInstance()
{
  static BH2Filter s_instance;
  return s_instance;
}


//______________________________________________________________________________
const std::vector<double>&
BH2Filter::GetXmax(int seg) const
{
  return m_param[seg].m_xmax;
}

//______________________________________________________________________________
const std::vector<double>&
BH2Filter::GetXmin(int seg) const
{
  return m_param[seg].m_xmin;
}

//______________________________________________________________________________
void
BH2Filter::Initialize(const std::string& filename)
{
  static const std::string funcname
    = std::string("[BH2Filter::") + __func__ + "]"; 

  std::ifstream f(filename.c_str());
  if (f.fail())
    {
      std::cerr << "#E " << funcname << " file open fail " 
		<<  filename << std::endl;
      std::exit(-1);
    }
  
  while (f.good())
    {
      std::string l;
      std::getline(f, l);
      if (l.empty() || l.find("#")!=std::string::npos)
	continue;
      std::istringstream iss(l);
      std::istream_iterator<double> issBegin(iss);
      std::istream_iterator<double> issEnd;
      std::vector<double> v(issBegin ,issEnd);
      if (v.size()<kNParam)
	{
// 	  std::cout << "#W " << funcname
// 		    << " number of parameters = " << v.size()
// 		    << ": required = " << kNParam
// 		    << std::endl;
	  continue;
	}
      const int bh2Seg  = static_cast<int>(v[kBH2Segment]);
      const int bcPlane = static_cast<int>(v[kLayerID]);
      const double xmin = v[kXMin];
      const double xmax = v[kXMax];

      int iplane = bcPlane - (PlOffsBc+13);
//       std::cout << " seg = "    << std::setw(2) << bh2Seg
// 		<< ", plane = " << std::setw(3) << bcPlane
// 		<< "(" << iplane << ")"
// 		<< ", xmin = "  << std::setw(5) << xmin
// 		<< ", xmax = "  << std::setw(5) << xmax
// 		<< std::endl;
      m_param[bh2Seg].m_xmin[iplane] = xmin;
      m_param[bh2Seg].m_xmax[iplane] = xmax;

    }
//   Print();
  std::cout << funcname << ": Initialization finished" << std::endl;
  return;
}

//______________________________________________________________________________
void
BH2Filter::Print(const std::string& arg) const
{
  std::cout << "#D BH2Filter::print" << std::endl;
  for (int i=0, n=m_param.size(); i<n; ++i)
    {
      std::stringstream ss;
      ss << "isegment " << i;
      m_param[i].Print(ss.str());
    }
  return;
}
