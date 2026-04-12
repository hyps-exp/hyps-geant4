/*
  RawData.cc

  2012/1/24
*/

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include "DetectorID.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"

#include "ConfMan.hh"
#include "DeleteObject.hh"
#include "TemplateLib.hh"

#ifdef MemoryLeak
debug::Counter RawData::sm_counter("RawData");
#endif

enum EDCDataType
  {
    kDCDataLeading,
    kDCDataTrailing,
    kNDCDataType
  };

RawData::RawData():     
  GCRawHC(0),
  BH1RawHC(0),
  BFTRawHC(),
  SFTRawHC(),
  CFTRawHC(),
  BGORawHC(0),
  T54CounterRawHC(),
  BH2RawHC(0),
  BACRawHC(0),
  
  TOFRawHC(0),
  ACRawHC(0),
  T0RawHC(0),
  AC1RawHC(0),
  LCRawHC(0),

  BcInRawHC(),
  BcOutRawHC(),
  SdcInRawHC(),
  SdcOutRawHC(),

  ScalerRawHC(0),
  MiscRawHC(0),
  MatrixRawHC(0)
{
#ifdef MemoryLeak
  ++sm_counter;
#endif
}

RawData::~RawData()
{
  clearAll();
#ifdef MemoryLeak
  --sm_counter;
#endif
}

bool RawData::AddHodoRawHit(HodoRHitContainer& cont,
			    int DetId,
			    int Plane,
			    int Seg,
			    int AorT,
			    int UorD,
			    int Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp(Data);
      else          p->SetAdcDown(Data);
    }else{
      if( UorD==0 ) p->SetTdcUp(Data);
      else          p->SetTdcDown(Data);
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddCFTHodoRawHit(int DetId,
			    int Plane,
			    int Seg,
			    int AorT,
			    int UorD,
			    double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = CFTRawHC[Plane];

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddBGOHodoRawHit(int DetId,
			       int Plane,
			       int Seg,
			       int AorT,
			       int UorD,
			       double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = BGORawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}


bool RawData::AddTOFHodoRawHit(int DetId,
			       int Plane,
			       int Seg,
			       int AorT,
			       int UorD,
			       double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = TOFRawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddLCHodoRawHit(int DetId,
			      int Plane,
			      int Seg,
			      int AorT,
			      int UorD,
			      double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = LCRawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddCHHodoRawHit(int DetId,
			      int Plane,
			      int Seg,
			      int AorT,
			      int UorD,
			      double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = CHRawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddACHodoRawHit(int DetId,
			      int Plane,
			      int Seg,
			      int AorT,
			      int UorD,
			      double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = ACRawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddT0HodoRawHit(int DetId,
			      int Plane,
			      int Seg,
			      int AorT,
			      int UorD,
			      double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = T0RawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}

bool RawData::AddPiVHodoRawHit(int DetId,
			       int Plane,
			       int Seg,
			       int AorT,
			       int UorD,
			       double Data )
{
  static const std::string funcname = "[RawData::AddHodoRawHit]";

  HodoRHitContainer& cont = PiVRawHC;

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->PlaneId()==Plane &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Plane, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    if( AorT==0 ){
      if( UorD==0 ) p->SetAdcUp((int)(Data*1000));
      else          p->SetAdcDown((int)(Data*1000));
    }else{
      if( UorD==0 ) p->SetTdcUp((int)(Data*1000));
      else          p->SetTdcDown((int)(Data*1000));
    }
    return true;
  }else{
    std::cerr << funcname << ": new fail. DetId="
              << DetId << " PlaneId=" << Plane << " SegId="
              << Seg << " A/T=" << AorT << " U/D=" << UorD << std::endl;
    return false;
  }
}


bool RawData::AddDCRawHit(DCRHitContainer& cont,
			  int Plane,
			  int Wire,
			  int Tdc,
			  int type)
{
  static const std::string funcname = "[RawData::AddDCRawHit]";

  DCRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    DCRawHit *q=cont[i];
    if( q->PlaneId()==Plane &&
	q->WireId()==Wire ){
       p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit( Plane, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
   switch(type)
      {
      case kDCDataLeading:
	p->SetTdc( Tdc );
	break;
      case kDCDataTrailing:
	p->SetTrailing(Tdc);
	break;
      default:
	break;
      }
   return true;
  }else{
    std::cerr << funcname << ": new fail. PlaneId="
              << Plane << " WireId="
              << Wire << std::endl;
    return false;
  }
}

bool RawData::AddSdcInRawHit( int Plane,
			      int Wire,
			      double dTime)
{
  static const std::string funcname = "[RawData::AddSdcInRawHit]";

  DCRHitContainer& cont = SdcInRawHC[Plane];

  DCRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    DCRawHit *q=cont[i];
    if( q->PlaneId()==Plane &&
	q->WireId()==Wire ){
       p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit( Plane, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetTdc( (int)(dTime*1000) );
    return true;
  }else{
    std::cerr << funcname << ": new fail. PlaneId="
              << Plane << " WireId="
              << Wire << std::endl;
    return false;
  }
}


bool RawData::AddSdcOutRawHit( int Plane,
			       int Wire,
			       double dTime)
{
  static const std::string funcname = "[RawData::AddSdcOutRawHit]";
  //int index = Plane - 30;
  int index = Plane - 30;

  DCRHitContainer& cont = SdcOutRawHC[index];

  DCRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    DCRawHit *q=cont[i];
    if( q->PlaneId()==Plane &&
	q->WireId()==Wire ){
       p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit( Plane, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetTdc( (int)(dTime*1000) );
    return true;
  }else{
    std::cerr << funcname << ": new fail. PlaneId="
              << Plane << " WireId="
              << Wire << std::endl;
    return false;
  }
}

bool RawData::AddBcOutRawHit( int Plane,
			       int Wire,
			       double dTime)
{
  static const std::string funcname = "[RawData::AddBcOutRawHit]";
  int index = Plane - (PlOffsBc + PlMinBcOut - 1);

  DCRHitContainer& cont = BcOutRawHC[index];
  Plane -=  PlOffsBc;

  DCRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    DCRawHit *q=cont[i];
    if( q->PlaneId()==Plane &&
	q->WireId()==Wire ){
       p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit( Plane, Wire );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetTdc( (int)(dTime*1000) );
    return true;
  }else{
    std::cerr << funcname << ": new fail. PlaneId="
              << Plane << " WireId="
              << Wire << std::endl;
    return false;
  }
}

void
RawData::clearAll()
{
  std::for_each(GCRawHC.begin(), GCRawHC.end(), DeleteObject());
  GCRawHC.clear();

  std::for_each(BH1RawHC.begin(), BH1RawHC.end(), DeleteObject());
  BH1RawHC.clear();

  for( int p=0; p<NumOfPlaneBFT; ++p){
    for_each( BFTRawHC[p].begin(),  BFTRawHC[p].end(), DeleteObject());
    BFTRawHC[p].clear();
  }

  for( int p=0; p<NumOfPlaneSFT; ++p){
    for_each( SFTRawHC[p].begin(),  SFTRawHC[p].end(), DeleteObject());
    SFTRawHC[p].clear();
  }

  for( int p=0; p<NumOfPlaneCFT; ++p){
    for_each( CFTRawHC[p].begin(),  CFTRawHC[p].end(), DeleteObject());
    CFTRawHC[p].clear();
  }

  std::for_each(BGORawHC.begin(), BGORawHC.end(), DeleteObject());
  BGORawHC.clear();

  for( int p=0; p<NumOfPlaneT54Counter; ++p){
    for_each( T54CounterRawHC[p].begin(),  T54CounterRawHC[p].end(), DeleteObject());
    T54CounterRawHC[p].clear();
  }

  std::for_each(BH2RawHC.begin(), BH2RawHC.end(), DeleteObject());
  BH2RawHC.clear();

  std::for_each(TOFRawHC.begin(), TOFRawHC.end(), DeleteObject());
  TOFRawHC.clear();

  std::for_each(CHRawHC.begin(), CHRawHC.end(), DeleteObject());
  CHRawHC.clear();

  std::for_each(ACRawHC.begin(), ACRawHC.end(), DeleteObject());
  ACRawHC.clear();

  std::for_each(T0RawHC.begin(), T0RawHC.end(), DeleteObject());
  T0RawHC.clear();

  std::for_each(LCRawHC.begin(), LCRawHC.end(), DeleteObject());
  LCRawHC.clear();

  std::for_each(AC1RawHC.begin(), AC1RawHC.end(), DeleteObject());
  AC1RawHC.clear();

  for( int l=0; l<=NumOfLayersBcIn; ++l){
    for_each( BcInRawHC[l].begin(),  BcInRawHC[l].end(),   DeleteObject());
    BcInRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersBcOut; ++l){
    for_each( BcOutRawHC[l].begin(), BcOutRawHC[l].end(),  DeleteObject());
    BcOutRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSsd; ++l){
    for_each( SsdRawHC[l].begin(), SsdRawHC[l].end(),  DeleteObject());
    SsdRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcIn; ++l){
    for_each( SdcInRawHC[l].begin(),  SdcInRawHC[l].end(),   DeleteObject());
    SdcInRawHC[l].clear();
  }

  for( int l=0; l<=NumOfLayersSdcOut; ++l){
    for_each( SdcOutRawHC[l].begin(), SdcOutRawHC[l].end(),  DeleteObject());
    SdcOutRawHC[l].clear();
  }

  // std::for_each(ScalerRawHC.begin(), ScalerRawHC.end(), DeleteObject());
  // ScalerRawHC.clear();

  // std::for_each(MiscRawHC.begin(), MiscRawHC.end(), DeleteObject());
  // MiscRawHC.clear();

  // std::for_each(MatrixRawHC.begin(), MatrixRawHC.end(), DeleteObject());
  // MatrixRawHC.clear();
    
  return;
}

bool
RawData::DecodeHits()
{
  /*
  hddaq::unpacker::UnpackerManager& gUnpacker
    = hddaq::unpacker::GUnpacker::get_instance();
  clearAll();

  //GC
  for( int AorT=0; AorT<2;  ++AorT ){
    int nhit = gUnpacker.get_entries( DetIdGC, 0, 0, 0, AorT );
    if( nhit>0 ){
      int data = gUnpacker.get( DetIdGC, 0, 0, 0, AorT );
      AddHodoRawHit( GCRawHC, DetIdGC, 0, 0, AorT, 0, data );
    }
    else continue;
  }  
  
  //BH1
  for( int seg=0; seg<NumOfSegBH1; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdBH1, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdBH1, 0, seg, UorD, AorT );
	  AddHodoRawHit( BH1RawHC, DetIdBH1, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }

  //BH2 
  for( int seg=0; seg<NumOfSegBH2; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdBH2, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdBH2, 0, seg, UorD, AorT );
	  AddHodoRawHit( BH2RawHC, DetIdBH2, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }


  //TOF
  for( int seg=0; seg<NumOfSegTOF; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdTOF, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdTOF, 0, seg, UorD, AorT );
	  AddHodoRawHit( TOFRawHC, DetIdTOF, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }
  
  //AC
  for( int seg=0; seg<NumOfSegAC; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdAC, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdAC, 0, seg, 0, AorT );
	AddHodoRawHit( ACRawHC, DetIdAC, 0, seg, AorT, 0, data );

      }
      else continue;
    }
  }

  //AC1
  for( int seg=0; seg<NumOfSegAC1; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdAC1, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdAC1, 0, seg, 0, AorT );
	AddHodoRawHit( AC1RawHC, DetIdAC1, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }

  //LC
  for( int seg=0; seg<NumOfSegLC; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      for( int UorD=0; UorD<2; ++UorD ){
	int nhit = gUnpacker.get_entries( DetIdLC, 0, seg, UorD, AorT );
	if( nhit>0 ){
	  int data = gUnpacker.get( DetIdLC, 0, seg, UorD, AorT );
	  AddHodoRawHit( LCRawHC, DetIdLC, 0, seg, AorT, UorD, data );
	}
	else continue;
      }
    }
  }
  
  // BC1&BC2 MWPC
  for(int plane=0; plane<NumOfLayersBcIn; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC1; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC1, plane, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int leading = gUnpacker.get( DetIdBC1, plane, 0, wire, 0, i )  ;

	   
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 leading, kDCDataLeading );
	    int trailing = gUnpacker.get( DetIdBC1, plane, 0, wire, 1, i );
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
	    	 trailing, kDCDataTrailing );

	  }
	}
	else continue; 
      }
    }
    else{
      for(int wire=0; wire<MaxWireBC2; ++wire){
      	int nhit = gUnpacker.get_entries( DetIdBC2, plane-NumOfLayersBc, 0, wire, 0 );
      	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int leading = gUnpacker.get( DetIdBC2, plane-NumOfLayersBc, 0, wire, 0, i );
	    int trailing = gUnpacker.get( DetIdBC2, plane-NumOfLayersBc, 0, wire, 1, i );
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 leading, kDCDataLeading); 
	    AddDCRawHit( BcInRawHC[plane+1], plane+PlMinBcIn, wire+1, 
			 trailing, kDCDataTrailing); 
	  }
	}	
	else continue;
      }
    }
  }
  
  // BC3&BC4 MWDC
  for(int plane=0; plane<NumOfLayersBcOut; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC3; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC3, plane, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdBC3, plane, 0, wire, 0, i);
	    AddDCRawHit( BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data );
	  }
	}
	else continue; 
      }
    }
    else{
      for(int wire=0; wire<MaxWireBC4; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC4, plane-NumOfLayersBc, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data =  gUnpacker.get( DetIdBC4, plane-NumOfLayersBc, 0, wire, 0, i );
	    AddDCRawHit( BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data );
	  }
	}	
	else continue;
      }
    }
  }

  //BFT
  for( int plane=0; plane<NumOfPlaneBFT; ++plane ){
    for(int seg = 0; seg<NumOfSegBFT; ++seg){
      int nhit = gUnpacker.get_entries( DetIdBFT, plane, 0, seg, 0 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int leading  = gUnpacker.get( DetIdBFT, plane, 0, seg, 0, i )  ;
	  int trailing = gUnpacker.get( DetIdBFT, plane, 0, seg, 1, i )  ;
	  AddHodoRawHit( BFTRawHC[plane], DetIdBFT, plane, seg , 1, 0, leading );
	  AddHodoRawHit( BFTRawHC[plane], DetIdBFT, plane, seg , 1, 1, trailing );
	}
      }
      else continue;
    }
  }

  // SFT
  for(int plane=0; plane<NumOfPlaneSFT; ++plane ){
    int NumOfSeg = 0;
    if(plane == SFT_X1 || plane == SFT_X2){
      NumOfSeg = NumOfSegSFT_X;
    }else{
      NumOfSeg = NumOfSegSFT_UV;
    }

    for(int seg = 0; seg<NumOfSeg; ++seg){
      int nhit = gUnpacker.get_entries( DetIdSFT, plane, 0, seg, 0 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int leading  = gUnpacker.get( DetIdSFT, plane, 0, seg, 0, i )  ;
	  int trailing = gUnpacker.get( DetIdSFT, plane, 0, seg, 1, i )  ;
	  AddHodoRawHit( SFTRawHC[plane], DetIdSFT, plane, seg , 1, 0, leading );
	  AddHodoRawHit( SFTRawHC[plane], DetIdSFT, plane, seg , 1, 1, trailing );
	}
      }
      
      else continue;
    }
  }


  // CFT
  for(int plane=0; plane<NumOfPlaneCFT; ++plane ){
    int NumOfSeg = 0;
    if(plane == CFT_PHI1 || plane == CFT_PHI2){
      NumOfSeg = NumOfSegCFT_PHI;
    }else{
      NumOfSeg = NumOfSegCFT_U;
    }

    for(int seg = 0; seg<NumOfSeg; ++seg){
      // TDC
      int nhit = gUnpacker.get_entries( DetIdCFT, plane, 0, seg, 0 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int leading  = gUnpacker.get( DetIdCFT, plane, 0, seg, 0, i )  ;
	  int trailing = gUnpacker.get( DetIdCFT, plane, 0, seg, 1, i )  ;
	  AddHodoRawHit( CFTRawHC[plane], DetIdCFT, plane, seg , 1, 0, leading );
	  AddHodoRawHit( CFTRawHC[plane], DetIdCFT, plane, seg , 1, 1, trailing );
	}
      }

      // ADC High Gain
      // ADC High Gain is stored as AorT=0(ADC) and UorD=0(Up)
      nhit = gUnpacker.get_entries( DetIdCFT, plane, 0, seg, 2 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int adc  = gUnpacker.get( DetIdCFT, plane, 0, seg, 2, i )  ;
	  AddHodoRawHit( CFTRawHC[plane], DetIdCFT, plane, seg , 0, 0, adc );
	}
      }

      // ADC Low Gain
      // ADC Low Gain is stored as AorT=0(ADC) and UorD=1(Down)
      nhit = gUnpacker.get_entries( DetIdCFT, plane, 0, seg, 3 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int adc  = gUnpacker.get( DetIdCFT, plane, 0, seg, 3, i )  ;
	  AddHodoRawHit( CFTRawHC[plane], DetIdCFT, plane, seg , 0, 1, adc );
	}
      }
    }
  }

  // SSD
  for(int plane=0; plane<NumOfLayersSsd;++plane){
    for(int wire=0;wire<768;++wire){
      int nhit = gUnpacker.get_entries(DetIdSsd, plane, 0, wire, 0);
      if(nhit==0) continue;
      for(int i=0;i<nhit;++i){
	int data = gUnpacker.get(DetIdSsd, plane, 0, wire, 0, i);
	AddDCRawHit(SsdRawHC[plane+1], plane+1, wire+1, data);
      }
    }
  }


  // SDC1&SDC2 MWDC
  for(int plane=0; plane<NumOfLayersSdcIn; ++plane ){
    if( plane<NumOfLayersSFT){
 //      for(int wire=0; wire<MaxWireSDC1; ++wire){
// 	int nhit = gUnpacker.get_entries( DetIdSDC1, plane, 0, wire, 2 );
// 	if( nhit>0 ){
// 	  for(int i=0; i<nhit; i++ ){
// 	    int data = (gUnpacker.get( DetIdSDC1, plane, 0, wire, 2, i ) & 0xffff);
// 	    AddDCRawHit( SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data-offset );
// 	  }
// 	}
// 	else continue; 
//       }
    }
    else{
      for(int wire=0; wire<MaxWireSDC2; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC2, plane-NumOfLayersSFT, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC2, plane-NumOfLayersSFT, 0, wire, 0, i );
	    AddDCRawHit( SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data );
	  }
	}	
	else continue;
      }
    }
  }

  // SDC3&SDC4
  for(int plane=0; plane<NumOfLayersSdcOut; ++plane ){
    if( plane<NumOfLayersSdc ){
      int NumOfWireSDC3;
      if( plane==1 || plane==4 ) NumOfWireSDC3 =108;
      else NumOfWireSDC3 =120;
      for(int wire=0; wire<NumOfWireSDC3; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC3, plane, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC3, plane, 0, wire, 0, i );
 	    AddDCRawHit( SdcOutRawHC[plane+1], plane+PlMinSdcOut, wire+1, data );
	  }
	}
	else continue;
      }
    }
    else{
      int NumOfWireSDC4;
      if( plane==(NumOfLayersSdc+1) || plane==(NumOfLayersSdc+4) ) NumOfWireSDC4 =108;
      else NumOfWireSDC4 =120;
      for(int wire=0; wire<NumOfWireSDC4; ++wire){
	int nhit = gUnpacker.get_entries( DetIdSDC4, plane-NumOfLayersSdc, 0, wire, 0 );
	if( nhit>0 ){
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC4, plane-NumOfLayersSdc, 0, wire, 0 ,i );
	    AddDCRawHit( SdcOutRawHC[plane+1],  plane+PlMinSdcOut, NumOfWireSDC4-wire, data );
	    //	    AddDCRawHit( SdcOutRawHC[plane+1],  plane+PlMinSdcOut, wire+1, data );
	  }
	}
	else continue; 
      }
    }
  }

  //Scaler
  for( int seg=0; seg<NumOfScaler; ++seg ){
    
    int nhit = gUnpacker.get_entries( DetIdScaler, 0, 0, seg, 0 );
    //std::cout<<"nhit_raw = "<<nhit<<"seg ="<<seg<<std::endl; 
    if( nhit>0 ){
      int data = gUnpacker.get( DetIdScaler, 0, 0, seg, 0 );
      AddHodoRawHit( ScalerRawHC, DetIdScaler, 0, 0, seg, 0, data );
    }
    
    else continue;
    
  }
  
  
 //Misc
  for( int seg=0; seg<NumOfMisc; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdMisc, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdMisc, 0, seg, 0, AorT );
	AddHodoRawHit( MiscRawHC, DetIdMisc, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }
  */
  /*
  //Matrix
  for( int seg=0; seg<NumOfMatrix; ++seg ){
    for( int AorT=0; AorT<2;  ++AorT ){
      int nhit = gUnpacker.get_entries( DetIdMatrix, 0, seg, 0, AorT );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdMatrix, 0, seg, 0, AorT );
	AddHodoRawHit( MatrixRawHC, DetIdMatrix, 0, seg, AorT, 0, data );
      }
      else continue;
    }
  }
  */

  return true;
}

const HodoRHitContainer& RawData::GetGCRawHC() const
{
  return GCRawHC;
}

const HodoRHitContainer& RawData::GetBH1RawHC() const
{
  return BH1RawHC;
}

const HodoRHitContainer& RawData::GetBH2RawHC() const
{
  return BH2RawHC;
}

const HodoRHitContainer& RawData::GetBACRawHC() const
{
  return BACRawHC;
}

const HodoRHitContainer& RawData::GetCHRawHC() const
{
  return CHRawHC;
}

const HodoRHitContainer& RawData::GetTOFRawHC() const
{
  return TOFRawHC;
}

const HodoRHitContainer& RawData::GetACRawHC() const
{
  return ACRawHC;
}

const HodoRHitContainer& RawData::GetT0RawHC() const
{
  return T0RawHC;
}

const HodoRHitContainer& RawData::GetAC1RawHC() const
{
  return AC1RawHC;
}

const HodoRHitContainer& RawData::GetLCRawHC() const
{
  return LCRawHC;
}

const DCRHitContainer & RawData::GetBcInRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcIn ) layer=0;
  return BcInRawHC[layer];
}

const DCRHitContainer & RawData::GetBcOutRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersBcOut ) layer=0;
  return BcOutRawHC[layer];
}

const DCRHitContainer & RawData::GetSsdRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSsd ) layer=0;
  return SsdRawHC[layer];
}

const DCRHitContainer & RawData::GetSdcInRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcIn ) layer=0;
  return SdcInRawHC[layer];
}

const DCRHitContainer & RawData::GetSdcOutRawHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSdcOut ) layer=0;
  return SdcOutRawHC[layer];
}

const HodoRHitContainer& RawData::GetScalerRawHC() const
{
  return ScalerRawHC;
}

const HodoRHitContainer& RawData::GetMiscRawHC() const
{
  return MiscRawHC;
}

const HodoRHitContainer& RawData::GetMatrixRawHC() const
{
  return MatrixRawHC;
}

const HodoRHitContainer& RawData::GetBFTRawHC( int plane ) const
{
  if( plane<0 || plane>NumOfPlaneBFT ) plane=0;
  return BFTRawHC[plane];
}

const HodoRHitContainer& RawData::GetSFTRawHC( int plane ) const
{
  if( plane<0 || plane>NumOfPlaneSFT ) plane=0;
  return SFTRawHC[plane];
}

const HodoRHitContainer& RawData::GetCFTRawHC( int plane ) const
{
  if( plane<0 || plane>NumOfPlaneCFT ) plane=0;
  return CFTRawHC[plane];
}

const HodoRHitContainer& RawData::GetBGORawHC() const
{
  return BGORawHC;
}

const HodoRHitContainer& RawData::GetPiVRawHC() const
{
  return PiVRawHC;
}

const HodoRHitContainer& RawData::GetT54CounterRawHC( int plane ) const
{
  if( plane<0 || plane>NumOfPlaneT54Counter ) plane=0;
  return T54CounterRawHC[plane];
}
