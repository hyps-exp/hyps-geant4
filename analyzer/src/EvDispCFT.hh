/*
  EvDispCFT.hh

  2012/1/24
*/

#ifndef EvDispCFT_h
#define EvDispCFT_h 1

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
#include "TLine.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TBox.h"
#include "TStyle.h"

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "ThreeVector.hh"
#include "DetectorID.hh"

#include <vector>

extern  TApplication *theApp;

typedef std::vector <TLine *> TLineContainer;

class EvDispCFT : public TObject {
private:
  TPad      *tp_[7];
  TCanvas   *tc_;
  TH2F      *hbase_;
  TH2F      *hbaseU_;
  TH1F      *hp_[NumOfLayersCFT];
  TH2F      *hbaseZX_;
  TH2F      *hbaseZY_;
  
  mutable TGraph    *gr_Vertex0XY_;
  mutable TGraph    *gr_Vertex0ZX_;
  mutable TGraph    *gr_Vertex0ZY_;

  mutable TGraph    *gr_VertexXY_;
  mutable TGraph    *gr_VertexZX_;
  mutable TGraph    *gr_VertexZY_;
  
  std::vector <TArc*> CFT_Arc_[NumOfLayersCFT];
  //TLineContainer BGO_Line_[NumOfBGOUnit*(NumOfBGOInOneUnit+NumOfBGOInOneUnit2)];
  //TLineContainer BGO_Line_[NumOfBGOUnit*NumOfBGOInOneUnit];
  TLineContainer BGO_Line_[NumOfSegBGO];

  TLineContainer PiV_Line_[NumOfPiVLayer][NumOfSegPiV];

  //TLine     *U_Line_[NumOfSegCFT_U][2];

  mutable std::vector <TLine*> TrackXYCont;
  mutable std::vector <TLine*> TrackZXCont;
  mutable std::vector <TLine*> TrackZYCont;

  mutable std::vector <TGraph*> CFTVtxContXY;
  mutable std::vector <TGraph*> CFTVtxContZX;
  mutable std::vector <TGraph*> CFTVtxContZY;

public:
  EvDispCFT(void);
  ~EvDispCFT(void);
  static EvDispCFT & GetInstance( void );
  void Initialize(int RunNum);
  void ConstructCFT(void);
  void ConstructBGO(void);
  void ConstructPiV(void);
  void ShowHitFiber(int layer, int segment, double pe) const;
  void ShowHitFiberProton(int layer, int segment) const;
  void ShowHitFiberPi(int layer, int segment) const;
  void ShowHitBGO(int segment, double de) const;
  void ShowHitBGO_Proton(int segment) const;
  void ShowHitBGO_Pi(int segment) const;
  void ShowHitPiV(int layer, int segment, double de) const;
  void ShowHitPiV_Proton(int layer, int segment) const;
  void ShowHitPiV_Pi(int layer, int segment) const;
  void DrawTrackInXYPlane(double x0, double y0, double x1, double y1, bool flagK=false, bool flagP=false) const;
  void DrawTrackInZXPlane(double z0, double x0, double z1, double x1, bool flagK=false, bool flagP=false) const;
  void DrawTrackInZYPlane(double z0, double y0, double z1, double y1, bool flagK=false, bool flagP=false) const;
  void DrawVertex0(double x, double y, double z) const;
  void DrawVertex(double x, double y, double z) const;
  void DrawCFTVertex(double x, double y, double z) const;
  void UpdateCanvas() const;
  void EndOfEvent() const;
  void FiberPos(int layer, int seg, double *x, double *y) const;
  void BGOPos(int seg, double *x, double *y) const;
  //void FiberPosPhi(int layer, int seg, double *x, double *y) const;
  //void FiberPosU(int layer, int seg, double z, double *x, double *y) const;
  void get_command(void) const;

private:
  static EvDispCFT *evDisp_;
  static TApplication *theApp;
};


#endif
