/*
  EvDisp.hh

  2012/1/24
*/

#ifndef EvDisp_h
#define EvDisp_h 1

#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeometry.h"
#include "TMixture.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTUBS.h"
#include "TTUBS.h"
#include "TRotMatrix.h"
#include "TNode.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView.h"
#include "TPad.h"
#include "TButton.h"
#include "TMarker3DBox.h"
#include "TPave.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TArc.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TBox.h"

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "ThreeVector.hh"
#include "DetectorID.hh"

#define MaxTrack      10
#define MaxBFTHit     1000
#define MaxBFTRateNum 10000
#define MaxComment    10
#define MaxVertex     10

extern  TApplication *theApp;

class DCLocalTrack;
class BFTCluster;

class EvDisp : public TObject {
private:
  TGeometry *gevdisp_;
  TNode     *node_;
  TPad      *tp_;
  TPad      *tp2_;
  TCanvas   *tc_;

  TGeometry *gevdisp_vtx_;
  TNode     *node_vtx_;
  TPad      *tp_vtx_;
  TCanvas   *tc_vtx_;

  TCanvas   *tc_bcIn_;
  TH2F      *hbase_bcIn_;

  TCanvas   *tc_bcOut_;
  TPad      *tp_bcOut1_;
  TPad      *tp_bcOut2_;
  TH2F      *hbase_bcOut_;
  TH2F      *hbase_bcOut_Y_;

  TCanvas   *tc_bft_;
  TH2F      *hbase_bft_;

  TCanvas   *tc_bft_rate_;

  TBRIK     *world_;
  TBRIK     *world_vtx_;

  /*
  TTUBE     *Target_Tube_;
  TNode     *Target_Node_;
  */
  //TNode     *Target_Node_vtx_;

  TBRIK     *Bh2Wall_Brik_;
  TNode     *Bh2Wall_Node_;

  TBRIK     *Bh2Seg_Brik_[NumOfSegBH2];
  TNode     *Bh2Seg_Node_[NumOfSegBH2];

  TBRIK     *Bh2Wall_Brik_vtx_;
  TNode     *Bh2Wall_Node_vtx_;

  TBRIK     *Bh2Seg_Brik_vtx_[NumOfSegBH2];
  TNode     *Bh2Seg_Node_vtx_[NumOfSegBH2];


  TTUBS     *Yoke1_Tubs_;
  TNode     *Yoke1_Node_;

  TTRD1     *Yoke2_Trd_;
  TNode     *Yoke2_Node_;

  TBRIK     *Yoke3_Brik_;
  TNode     *Yoke3_Node_;

  TTRD1     *Yoke4_Trd_;
  TNode     *Yoke4_Node_;

  TTRD1     *Yoke5_Trd_;
  TNode     *Yoke5_Node_;
  /*
  TNode     *Bdc3x_Node_[MaxWireBDC];
  TNode     *Bdc3xp_Node_[MaxWireBDC];
  TNode     *Bdc3u_Node_[MaxWireBDC];
  TNode     *Bdc3up_Node_[MaxWireBDC];
  TNode     *Bdc3v_Node_[MaxWireBDC];
  TNode     *Bdc3vp_Node_[MaxWireBDC];

  TNode     *Bdc3x_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3xp_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3u_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3up_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3v_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3vp_Node_vtx_[MaxWireBDC];

  TNode     *Bdc4x_Node_[MaxWireBDC];
  TNode     *Bdc4xp_Node_[MaxWireBDC];
  TNode     *Bdc4u_Node_[MaxWireBDC];
  TNode     *Bdc4up_Node_[MaxWireBDC];
  TNode     *Bdc4v_Node_[MaxWireBDC];
  TNode     *Bdc4vp_Node_[MaxWireBDC];

  TNode     *Bdc4x_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4xp_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4u_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4up_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4v_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4vp_Node_vtx_[MaxWireBDC];
  */

  TNode     *Sdc1u1_Node_[MaxWireSDC1];
  TNode     *Sdc1u2_Node_[MaxWireSDC1];
  TNode     *Sdc1v1_Node_[MaxWireSDC1];
  TNode     *Sdc1v2_Node_[MaxWireSDC1];

  TNode     *Sdc1u1_Node_vtx_[MaxWireSDC1];
  TNode     *Sdc1u2_Node_vtx_[MaxWireSDC1];
  TNode     *Sdc1v1_Node_vtx_[MaxWireSDC1];
  TNode     *Sdc1v2_Node_vtx_[MaxWireSDC1];

  TNode     *SftX_Node_[NumOfSegSFT_X];
  TNode     *SftV_Node_[NumOfSegSFT_UV];
  TNode     *SftU_Node_[NumOfSegSFT_UV];

  TNode     *SftX_Node_vtx_[NumOfSegSFT_X];
  TNode     *SftV_Node_vtx_[NumOfSegSFT_UV];
  TNode     *SftU_Node_vtx_[NumOfSegSFT_UV];

  TNode     *Sdc2v1_Node_[MaxWireSDC2];
  TNode     *Sdc2v2_Node_[MaxWireSDC2];
  TNode     *Sdc2u1_Node_[MaxWireSDC2];
  TNode     *Sdc2u2_Node_[MaxWireSDC2];
  TNode     *Sdc2x1_Node_[MaxWireSDC2];
  TNode     *Sdc2x2_Node_[MaxWireSDC2];

  TNode     *Sdc2v1_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2v2_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2u1_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2u2_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2x1_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2x2_Node_vtx_[MaxWireSDC2];
  
  TNode     *Sdc3v1_Node_[MaxWireSDC3V];
  TNode     *Sdc3x1_Node_[MaxWireSDC3X];
  TNode     *Sdc3u1_Node_[MaxWireSDC3U];
  TNode     *Sdc3v2_Node_[MaxWireSDC3V];
  TNode     *Sdc3x2_Node_[MaxWireSDC3X];
  TNode     *Sdc3u2_Node_[MaxWireSDC3U];

  TNode     *Sdc4v1_Node_[MaxWireSDC4V];
  TNode     *Sdc4x1_Node_[MaxWireSDC4X];
  TNode     *Sdc4u1_Node_[MaxWireSDC4U];
  TNode     *Sdc4v2_Node_[MaxWireSDC4V];
  TNode     *Sdc4x2_Node_[MaxWireSDC4X];
  TNode     *Sdc4u2_Node_[MaxWireSDC4U];

  TBRIK     *TofWall_Brik_;
  TNode     *TofWall_Node_;

  TBRIK     *TofSeg_Brik_;
  TNode     *TofSeg_Node_[NumOfSegTOF];
  /*
  TBRIK     *Ac1_Brik_;
  TNode     *Ac1_Node_;

  TBRIK     *Ac2_Brik_;
  TNode     *Ac2_Node_;
  */
  TBRIK     *LcWall_Brik_;
  TNode     *LcWall_Node_;

  TBRIK     *LcSeg_Brik_;
  TNode     *LcSeg_Node_[NumOfSegLC];

  mutable TPolyMarker3D   *InitStepMark_;

  //mutable TPolyLine3D     *LocalTrackBdcOut_[MaxTrack];
  mutable TPolyLine3D     *LocalTrackSdcIn_[MaxTrack];
  mutable TPolyLine3D     *LocalTrackSdcOut_[MaxTrack];
  mutable TLine     *LocalTrackSdcIn2_[MaxTrack];
  mutable TLine     *LocalTrackSdcIn2_Y_[MaxTrack];
  mutable TLine     *LocalTrackBcIn2_[MaxTrack];
  mutable TLine     *LocalTrackBcOut2_[MaxTrack];
  mutable TLine     *LocalTrackBcOut2_Y_[MaxTrack];
  mutable TGraph    *BcOut_YPos_gr_[MaxTrack];
  mutable TGraph    *SdcIn_YPos_gr_[MaxTrack];
  //mutable TPolyMarker3D   *VtxPoint_[MaxTrack];

  mutable TPolyMarker3D   *SksStepMark_;

  //mutable int nTrackBdcOut_;
  mutable int nTrackSdcIn_;
  mutable int nTrackSdcOut_;
  mutable int nTrackBcIn_;
  mutable int nTrackBcOut_;
  mutable int nVertex_;

  TArc *BcIn_Arc_[NumOfLayersBcIn][MaxWireBC1];
  TLine *Bh1_Line_[NumOfSegBH1][4];
  TArc *Bft_Arc_[NumOfPlaneBFT][NumOfSegBFT];

  TArc *Bc3_Arc_[6][MaxWireBC3];
  TArc *Bc4_Arc_[6][MaxWireBC4];
  TLine *Bh2_Line_[NumOfSegBH2][4];

  TBox *Target_Box_;
  TBox *Target_Box_Y_;

  mutable TArc *Vertex_Arc_[MaxVertex];
  mutable TArc *Vertex_Arc_Y_[MaxVertex];

  TArc *Sdc1_Arc_[4][MaxWireSDC1];

  TArc *Sft_Arc_[4][NumOfSegSFT_UV];
  TArc *Sdc2_Arc_[6][MaxWireSDC2];

  mutable int nClusterBFT_;
  mutable int nHitBFT_;
  mutable TArc     *Bft_Cl_Arc_[MaxBFTHit];
  mutable TLatex   *Bft_tex_, *tmp_Bft_tex_;
  TLine *TimeGate_Line_[2];


  mutable TGraph *Bft_rate_gr_;
  mutable double Bft_rate_[MaxBFTRateNum], Bft_event_Num_[MaxBFTRateNum];
  mutable int    Bft_rate_Index_;
  mutable TLine     *Bft_rate_line_;
  TH1F   *hBft_rate_;

  int RunNum_;
  mutable TLatex   *EvNum_tex_, *tmp_EvNum_tex_;
  mutable TLatex   *Sks_tex_, *tmp_Sks_tex_[MaxComment];
  mutable TLatex   *vtx_tex_, *tmp_vtx_tex_[MaxComment];
  mutable TLatex   *K18_tex_, *tmp_K18_tex_[MaxComment];

public:
  EvDisp(void);
  ~EvDisp(void);
  static EvDisp & GetInstance( void );
  void Initialize(int RunNum);
  void ConstructTarget(void);
  void ConstructTargetVtx(void);
  void ConstructBH2(void);
  void ConstructBH2Vtx(void);
  void ConstructSKS(void);
  void ConstructBdcOut(void);
  void ConstructBdcOutVtx(void);
  void ConstructSdcIn(void);
  void ConstructSdcInVtx(void);
  void ConstructBcIn(void);
  void ConstructBcOut(void);
  void ConstructSdcOut(void);
  void ConstructTOF(void);
  void ConstructRangeCounter(void);
  void ConstructRangeCounterVtx(void);
  void DrawInitTrack(int nStep, ThreeVector *StepPoint) const;
  void DrawInitTrack(void) const;
  void DrawHitWire(int lnum, int hit_wire, bool range_check=true, bool tdc_check=true) const;
  void DrawSFTWire(int lnum, int hit_wire) const;
  void DrawTrackWire(int lnum, int hit_wire, int it) const;
  void UpdateBcInCanvas() const;
  void UpdateBcOutCanvas() const;
  void UpdateBftCanvas() const;
  void DrawHitWireBcIn(int lnum, int hit_wire) const;
  void DrawHitWireBcOut(int lnum, int hit_wire, int color=0) const;
  void DrawHitWireSdcIn(int lnum, int hit_wire, int color=0) const;
  void DrawHitWireSft(int lnum, int hit_wire, double time, int color=0) const;
  void DrawHitBft(int lnum, int seg, bool timeflag) const;
  void DrawBFTCluster(BFTCluster *cl) const;
  void DrawBFTCluster(double time, double pos) const;
  void DrawHitBH1(int seg, int timeflag) const;
  void DrawHitBH2(int seg, bool timeflag) const;
  void DrawHitBH2(int seg, int Tu, int Td) const;
  void DrawHitHodoscope(int lnum, int seg, int Tu, int Td) const;
  void DrawHitRangeCounter(int lnum, int seg, int Tu, int Td) const;
  void DrawBdcOutLocalTrack(ThreeVector globalPos0, 
			    ThreeVector globalPos1) const;
  void DrawBdcOutLocalTrack(DCLocalTrack *tp) const;
  void DrawSdcInLocalTrack(ThreeVector globalPos0, 
			    ThreeVector globalPos1) const;
  void DrawSdcInLocalTrack(DCLocalTrack *tp) const;
  void DrawBcInLocalTrack(DCLocalTrack *tp) const;
  void DrawBcInLocalTrack(double xbft, double xbh1, bool flagBh1) const;
  void DrawBcOutLocalTrack(DCLocalTrack *tp) const;
  void DrawSdc1LocalTrack(DCLocalTrack *tp) const;
  void DrawSdcOutLocalTrack(ThreeVector globalPos0, 
			    ThreeVector globalPos1) const;
  void DrawSdcOutLocalTrack(DCLocalTrack *tp) const;
  void DrawVertex(ThreeVector vtxPoint, double cost) const;
  void DrawSksTrack(int nStep, ThreeVector *StepPoint) const;
  void WriteEventNumber(int spill, int evNum) const;
  void WriteSksComment(char *buf) const;
  void WriteVertexComment(char *buf) const;
  void WriteK18Comment(char *buf) const;
  void DrawVertex2D(double x, double y, double z, double cdist) const;
  void EndOfEvent(void) const;
  void ResetVisibility() const;
  void calcRotMatrix(double TA, double RA1, double RA2, double *rotMat);
  void get_command(void) const;

private:
  static EvDisp *evDisp_;
  static TApplication *theApp;
};


#endif
