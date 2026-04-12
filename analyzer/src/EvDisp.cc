/*
  EvDisp.cc

  2012/1/24
*/

//Unit is mm.

#include "EvDisp.hh"
#include "DCGeomMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <string>

#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
//#include "BFTCluster.hh"
//#include "BFTHit.hh"

const int MaxChar = 200;

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

EvDisp *EvDisp::evDisp_ = 0;
TApplication *EvDisp::theApp =0;

const int NumColor = 20;
const int Color[NumColor] = {kBlue, kRed, kCyan,  kGreen, kPink, kOrange, kYellow, kMagenta, kViolet, kSpring, kBlue-10, kRed-10, kCyan-10,  kGreen-10, kPink-10, kOrange-10, kYellow-10, kMagenta-10, kViolet-10, kSpring-10};

EvDisp::EvDisp(void)
{
  InitStepMark_ = NULL;
  /*  
  for (int i=0; i<MaxTrack; i++)
    LocalTrackBdcOut_[i] = NULL;
  */
  for (int i=0; i<MaxTrack; i++) {
    LocalTrackSdcIn_[i] = NULL;
    LocalTrackSdcIn2_[i] = NULL;
    LocalTrackSdcIn2_Y_[i] = NULL;
    LocalTrackBcIn2_[i] = NULL;
    LocalTrackBcOut2_[i] = NULL;
    LocalTrackBcOut2_Y_[i] = NULL;
    BcOut_YPos_gr_[i] = NULL;
    SdcIn_YPos_gr_[i] = NULL;
  }

  for (int i=0; i<MaxBFTHit; i++)
    Bft_Cl_Arc_[i] = NULL;

  Bft_tex_ = NULL;
  tmp_Bft_tex_ = NULL;

  Bft_rate_gr_ = NULL;
  Bft_rate_line_ = NULL;

  EvNum_tex_ = NULL;
  tmp_EvNum_tex_ = NULL;

  Sks_tex_ = NULL;
  for (int i=0; i<MaxComment; i++)
    tmp_Sks_tex_[i] = NULL;

  vtx_tex_ = NULL;
  for (int i=0; i<MaxComment; i++)
    tmp_vtx_tex_[i] = NULL;

  K18_tex_ = NULL;
  for (int i=0; i<MaxComment; i++)
    tmp_K18_tex_[i] = NULL;

  for (int i=0; i<MaxTrack; i++)
    LocalTrackSdcOut_[i] = NULL;

  for (int i=0; i<MaxVertex; i++) {
    Vertex_Arc_[i] = NULL;
    Vertex_Arc_Y_[i] = NULL;
  }
  /*
  for (int i=0; i<MaxTrack; i++)
    VtxPoint_[i] = NULL;
  */
  SksStepMark_ = NULL;
}

EvDisp::~EvDisp(void)
{
  if (InitStepMark_)
    delete InitStepMark_;
  /*
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBdcOut_[i])
      delete LocalTrackBdcOut_[i];
  */
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn_[i])
      delete LocalTrackSdcIn_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBcIn2_[i])
      delete LocalTrackBcIn2_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBcOut2_[i])
      delete LocalTrackBcOut2_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBcOut2_Y_[i])
      delete LocalTrackBcOut2_Y_[i];

  for (int i=0; i<MaxTrack; i++)
    if (BcOut_YPos_gr_[i])
      delete BcOut_YPos_gr_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn2_[i])
      delete LocalTrackSdcIn2_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn2_Y_[i])
      delete LocalTrackSdcIn2_Y_[i];

  for (int i=0; i<MaxTrack; i++)
    if (SdcIn_YPos_gr_[i])
      delete SdcIn_YPos_gr_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcOut_[i])
      delete LocalTrackSdcOut_[i];
  /*
  for (int i=0; i<MaxTrack; i++)
    if (VtxPoint_[i])
      delete VtxPoint_[i];
  */
  if (SksStepMark_)
    delete SksStepMark_;

  for (int i=0; i<MaxBFTHit; i++)
    if (Bft_Cl_Arc_[i])
      delete Bft_Cl_Arc_[i];

  if (Bft_tex_)
    delete Bft_tex_;

  if (tmp_Bft_tex_)
    delete tmp_Bft_tex_;

  if (EvNum_tex_)
    delete EvNum_tex_;

  if (tmp_EvNum_tex_)
    delete tmp_EvNum_tex_;

  if (Sks_tex_)
    delete Sks_tex_;

  for (int i=0; i<MaxComment; i++)
    if (tmp_Sks_tex_[i])
      delete tmp_Sks_tex_[i];

  if (vtx_tex_)
    delete vtx_tex_;

  for (int i=0; i<MaxComment; i++)
    if (tmp_vtx_tex_[i])
      delete tmp_vtx_tex_[i];

  if (K18_tex_)
    delete K18_tex_;

  for (int i=0; i<MaxComment; i++)
    if (tmp_K18_tex_[i])
      delete tmp_K18_tex_[i];

  for (int i=0; i<MaxVertex; i++)
    if (Vertex_Arc_[i])
      delete Vertex_Arc_[i];

  for (int i=0; i<MaxVertex; i++)
    if (Vertex_Arc_Y_[i])
      delete Vertex_Arc_Y_[i];


  if (Bft_rate_gr_)
    delete Bft_rate_gr_;

  if (Bft_rate_line_)
    delete Bft_rate_line_;

  if (hBft_rate_)
    delete hBft_rate_;

  delete node_;
  /*
  delete Target_Node_;
  delete Target_Tube_;
  */
  delete Yoke1_Node_;
  delete Yoke1_Tubs_;

  delete Yoke2_Node_;
  delete Yoke2_Trd_;

  delete Yoke3_Node_;
  delete Yoke3_Brik_;

  delete Yoke4_Node_;
  delete Yoke4_Trd_;

  delete Yoke5_Node_;
  delete Yoke5_Trd_;
  /*
  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3x_Node_[i];
    delete Bdc3x_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3xp_Node_[i];
    delete Bdc3xp_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3u_Node_[i];
    delete Bdc3u_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3up_Node_[i];
    delete Bdc3up_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3v_Node_[i];
    delete Bdc3v_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) { 
    delete Bdc3vp_Node_[i];
    delete Bdc3vp_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4x_Node_[i];
    delete Bdc4x_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4xp_Node_[i];
    delete Bdc4xp_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4u_Node_[i];
    delete Bdc4u_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4up_Node_[i];
    delete Bdc4up_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4v_Node_[i];
    delete Bdc4v_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4vp_Node_[i];
    delete Bdc4vp_Node_vtx_[i];
  }
  */

  for (int i=0; i<MaxWireSDC1; i++) {
    delete Sdc1u1_Node_[i];
    delete Sdc1u2_Node_[i];
    delete Sdc1v1_Node_[i];
    delete Sdc1v2_Node_[i];

    delete Sdc1u1_Node_vtx_[i];
    delete Sdc1u2_Node_vtx_[i];
    delete Sdc1v1_Node_vtx_[i];
    delete Sdc1v2_Node_vtx_[i];
  }

  for (int i=0; i<NumOfSegSFT_X; i++) {
    delete SftX_Node_[i];
    delete SftX_Node_vtx_[i];
  }
  for (int i=0; i<NumOfSegSFT_UV; i++) {
    delete SftU_Node_[i];
    delete SftU_Node_vtx_[i];

    delete SftV_Node_[i];
    delete SftV_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireSDC2; i++) {
    delete Sdc2u1_Node_[i];
    delete Sdc2u2_Node_[i];
    delete Sdc2v1_Node_[i];
    delete Sdc2v2_Node_[i];
    delete Sdc2x1_Node_[i];
    delete Sdc2x2_Node_[i];

    delete Sdc2u1_Node_vtx_[i];
    delete Sdc2u2_Node_vtx_[i];
    delete Sdc2v1_Node_vtx_[i];
    delete Sdc2v2_Node_vtx_[i];
    delete Sdc2x1_Node_vtx_[i];
    delete Sdc2x2_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireSDC3V; i++) {
    delete Sdc3v1_Node_[i];
    delete Sdc3v2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC3X; i++) {
    delete Sdc3x1_Node_[i];
    delete Sdc3x2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC3U; i++) {
    delete Sdc3u1_Node_[i];
    delete Sdc3u2_Node_[i];
  }

  for (int i=0; i<MaxWireSDC4V; i++) {
    delete Sdc4v1_Node_[i];
    delete Sdc4v2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC4X; i++) {
    delete Sdc4x1_Node_[i];
    delete Sdc4x2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC4U; i++) {
    delete Sdc4u1_Node_[i];
    delete Sdc4u2_Node_[i];
  }


  for (int i=0; i<NumOfSegBH2; i++) {
    delete Bh2Seg_Node_[i];
    delete Bh2Seg_Brik_[i];
    delete Bh2Seg_Node_vtx_[i];
    delete Bh2Seg_Brik_vtx_[i];
  }
  delete Bh2Wall_Brik_;
  delete Bh2Wall_Node_;
  delete Bh2Wall_Brik_vtx_;
  delete Bh2Wall_Node_vtx_;

  for (int i=0; i<NumOfSegTOF; i++) 
    delete TofSeg_Node_[i];
  delete TofSeg_Brik_;

  delete TofWall_Brik_;
  delete TofWall_Node_;
  /*
  delete Ac1_Brik_;
  delete Ac1_Node_;

  delete Ac2_Brik_;
  delete Ac2_Node_;
  */
  for (int i=0; i<NumOfSegLC; i++) 
    delete LcSeg_Node_[i];
  delete LcSeg_Brik_;

  delete LcWall_Brik_;
  delete LcWall_Node_;


  delete gevdisp_;
  delete tp_;
  delete tp2_;
  delete tc_;


  delete gevdisp_vtx_;
  delete tp_vtx_;
  delete tc_vtx_;

  delete hbase_bcIn_;
  delete tc_bcIn_;

  delete hbase_bcOut_;
  delete hbase_bcOut_Y_;
  delete tp_bcOut1_;
  delete tp_bcOut2_;
  delete tc_bcOut_;

  delete hbase_bft_;
  delete tc_bft_;
  delete tc_bft_rate_;

}

EvDisp & EvDisp::GetInstance( void )
{
  if( !evDisp_ ){
    evDisp_ = new EvDisp();
  }
  if( !theApp ){
    theApp=new TApplication( "App", 0, 0 );
  }
  return *evDisp_;
}

void EvDisp::Initialize(int RunNum)
{
  RunNum_=RunNum;

  //nTrackBdcOut_=0;
  nTrackBcIn_=0;
  nTrackBcOut_=0;
  nTrackSdcIn_=0;
  nTrackSdcOut_=0;
  nVertex_=0;

  gevdisp_ = new TGeometry("evdisp","K1.8 Event Display");

  ThreeVector worldSize(1000.0, 1000.0, 1000.0); /*mm*/
  world_ = new TBRIK("world_","world","void",
		      worldSize.x(), worldSize.y(), worldSize.z());

  node_ = new TNode("node_","node","world_",0.0,0.0,0.0);
  gevdisp_->GetNode("node_")->SetVisibility(0);

  //ConstructBH2();
  std::cout << "Finish Construction of BH2" << std::endl;
  //ConstructTarget();
  //std::cout << "Finish Construction of Target" << std::endl;
  ConstructSKS();
  std::cout << "Finish Construction of SKS" << std::endl;
  /*
  ConstructBdcOut();
  std::cout << "Finish Construction of BdcOut" << std::endl;
  */
  ConstructSdcIn();
  std::cout << "Finish Construction of SdcIn" << std::endl;
  ConstructSdcOut();
  std::cout << "Finish Construction of SdcOut" << std::endl;

  ConstructTOF();
  std::cout << "Finish Construction of Tof" << std::endl;


  tc_ = new TCanvas("canvas","K1.8 Event Display",700,700);
  tp_ = new TPad("pad","K1.8 Event",0.0,0.0,1.0,0.799,10);
  tp_->Draw();
  tp2_ = new TPad("pad2","K1.8 Event",0.0,0.801,1.0,1.0,10);
  tp2_->Draw();

  tp_->cd();
  gevdisp_->Draw();
  //tp_->GetView()->SetParallel();

  tc_->cd();

  tc_->Update();

  gevdisp_vtx_ = new TGeometry("evdisp_vertex","K1.8 Event Display");

  world_vtx_ = new TBRIK("world_vtx_","world","void",
		      worldSize.x(), worldSize.y(), worldSize.z());

  node_vtx_ = new TNode("node_vtx_","node_vtx","world_vtx_",0.0,0.0,0.0);
  gevdisp_vtx_->GetNode("node_vtx_")->SetVisibility(0);


  //ConstructBH2Vtx();
  std::cout << "Finish Construction of BH2Vtx" << std::endl;

  /*
  ConstructTargetVtx();
  std::cout << "Finish Construction of Target at vertex display" << std::endl;
  ConstructBdcOutVtx();
  std::cout << "Finish Construction of BdcOut at vertex display" << std::endl;
  */
  ConstructSdcInVtx();
  std::cout << "Finish Construction of SdcIn at vertex display" << std::endl;

  tc_vtx_ = new TCanvas("canvas_vtx","E559 Event Display",700,700);
  tp_vtx_ = new TPad("pad_vtx","E559 Event",0.0,0.0,1.0,1.0,10);
  tp_vtx_->Draw();

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tp_vtx_->GetView()->ZoomIn();
  //tp_vtx_->GetView()->SetParallel();
  tc_vtx_->cd();

  tc_vtx_->Update();

  tc_bcIn_ = new TCanvas("tc_bdIn_","BcIn tracking",500,500);
  tc_bcIn_->SetFrameFillColor(kGray);
  tc_bcIn_->cd();
  hbase_bcIn_ = new TH2F("hbase_bcIn", "BcIn tracking", 300, -150, 150, 300, -700, 0);
  hbase_bcIn_->Draw();

  ConstructBcIn();

  tc_bcIn_->Update();

  tc_bcOut_ = new TCanvas("tc_bdOut_","BcOut tracking",500,500);
  tp_bcOut1_ = new TPad("tp_bcOut1_","XZ plane",0.0,0.0,0.49,0.99,10);
  tp_bcOut2_ = new TPad("tp_bcOut2_","YZ plane",0.51,0.0,0.99,0.99,10);
  //tc_bcOut_->SetFrameFillColor(kGray);
  tp_bcOut1_->SetFrameFillColor(kGray);
  tp_bcOut2_->SetFrameFillColor(kGray);
  tp_bcOut1_->Draw();
  tp_bcOut2_->Draw();

  //tc_bcOut_->cd();

  tp_bcOut1_->cd();
  hbase_bcOut_ = new TH2F("hbase_bcOut", "BcOut tracking (XZ plane)", 200, -300, 300, 300, -2000, 700);
  hbase_bcOut_->Draw();

  tp_bcOut2_->cd();
  hbase_bcOut_Y_ = new TH2F("hbase_bcOut_Y", "BcOut tracking (YZ plane)", 200, -300, 300, 300, -2000, 700);
  hbase_bcOut_Y_->Draw();

  ConstructBcOut();

  tc_bcOut_->Update();


  nClusterBFT_=0;
  nHitBFT_=0;
  Bft_rate_Index_=0;

  tc_bft_ = new TCanvas("tc_bft_","BFT Hit information",500,500);
  tc_bft_->SetFrameFillColor(kBlack);
  tc_bft_->cd();
  hbase_bft_ = new TH2F("hbase_bft", "BFT Hit information", 200, -100, 100, 1000, -700, 300);
  hbase_bft_->Draw();

  TimeGate_Line_[0] = new TLine(-100, -10, 100, -10);
  TimeGate_Line_[0]->SetLineColor(kWhite);
  TimeGate_Line_[0]->Draw("same");

  TimeGate_Line_[1] = new TLine(-100, 10, 100, 10);
  TimeGate_Line_[1]->SetLineColor(kWhite);
  TimeGate_Line_[1]->Draw("same");

  tc_bft_rate_ = new TCanvas("tc_bft_rate_","BFT Rate Dependence",500,500);
  tc_bft_rate_->Divide(1,2);

  tc_bft_rate_->cd(2);
  hBft_rate_ = new TH1F("hBft_rate_","BFT rate", 500, 0, 500);
  hBft_rate_->GetXaxis()->SetTitle("Rate (MHz)");
  tc_bft_rate_->Draw();

  ResetVisibility();
}
#if 0
void EvDisp::ConstructTarget(void)
{
  static const std::string funcname = "EvDisp::ConstructSKS";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //-----Target
  double TargetRmin = 0.0; 
  double TargetRmax = 33.9; 
  double TargetZ    = 110.0/2.0; 

  Target_Tube_ = new TTUBE("Target_Tube_","Target_Tube_", "void", 
			   TargetRmin, TargetRmax, TargetZ);

  double rotMatTarget[9];
  calcRotMatrix(90.0, 90.0, -50.0, rotMatTarget);
  TRotMatrix *rotTarget = new TRotMatrix("rotTarget","rotTarget",rotMatTarget);

  int lnum = geomMan.GetDetectorId("Target");
  //ThreeVector GlobalPos = geomMan.GetGlobalPosition(lnum);
  ThreeVector GlobalPos = geomMan.Local2GlobalPos(lnum, 
	  ThreeVector(0.0, 0.0, geomMan.GetLocalZ(lnum)+30.0));

  Target_Node_ = new TNode("Target_Node_", "Target_Node", "Target_Tube_",
				   GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
				   "rotTarget", "void");
}

void EvDisp::ConstructTargetVtx(void)
{
  static const std::string funcname = "EvDisp::ConstructSKS";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //-----Target
  double TargetRmin = 0.0; 
  double TargetRmax = 33.9; 
  double TargetZ    = 110.0/2.0; 

  TTUBE *Target_Tube_vtx_ = new TTUBE("Target_Tube_vtx_","Target_Tube_vtx_", "void", 
			   TargetRmin, TargetRmax, TargetZ);

  //-----Target
  double rotMatTarget[9];
  calcRotMatrix(90.0, 90.0, -50.0, rotMatTarget);
  TRotMatrix *rotTarget = new TRotMatrix("rotTarget","rotTarget",rotMatTarget);

  int lnum = geomMan.GetDetectorId("Target");
  //ThreeVector GlobalPos = geomMan.GetGlobalPosition(lnum);
  ThreeVector GlobalPos = geomMan.Local2GlobalPos(lnum, 
	  ThreeVector(0.0, 0.0, geomMan.GetLocalZ(lnum)+30.0));

  Target_Node_vtx_ = new TNode("Target_Node_vtx_", "Target_Node_vtx", "Target_Tube_vtx_",
				   GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
				   "rotTarget", "void");
}
#endif

void EvDisp::ConstructBH2(void)
{
  static const std::string funcname = "EvDisp::ConstructBH2";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BH2
  double rotMatBh2[9];
  double Bh2WallX = 130.0/2.0;
  double Bh2WallY = 14.0/2.0;
  double Bh2WallZ = 60.0/2.0;

  double Bh2SegX[8] = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  double Bh2SegY[8] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};
  double Bh2SegZ[8] = {40./2., 40./2., 40./2., 40./2., 40./2., 40./2., 40./2., 40./2.};

  double localPosX[8] = {-41.5, -19.0, -10.5, -3.5, 3.5, 10.5, 19.0, 41.5};
  double localPosZ[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

  lnum = geomMan.GetDetectorId("BH2(global)");
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBh2);

  TRotMatrix *rotBh2 = new TRotMatrix("rotBh2","rotBh2",rotMatBh2);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector bh2WallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector bh2WallGlobalPos = geomMan.Local2GlobalPos(lnum, bh2WallLocalPos);
  sprintf(object_name, "Bh2Wall_Brik_");
  Bh2Wall_Brik_ = new TBRIK(object_name, object_name, "void", 
			    Bh2WallX, Bh2WallY, Bh2WallZ);    
  sprintf(node_name, "Bh2Wall_Node_");
  Bh2Wall_Node_ = new TNode(node_name, node_name, object_name,
			    bh2WallGlobalPos.x(),
			    bh2WallGlobalPos.y(), 
			    bh2WallGlobalPos.z(),
			    "rotBh2", "void");
  Bh2Wall_Node_->SetVisibility(0);

  Bh2Wall_Node_->cd();

  for (int i=0; i<NumOfSegBH2; i++) {
    sprintf(object_name, "Bh2Seg_Brik_%d", i);
    Bh2Seg_Brik_[i] = new TBRIK(object_name, object_name, "void", 
				Bh2SegX[i], Bh2SegY[i], Bh2SegZ[i]);    

    ThreeVector bh2SegLocalPos = 
      ThreeVector(localPosX[i], 0.0, localPosZ[i]);
    sprintf(node_name, "Bh2Seg_Node_%d", i);
    Bh2Seg_Node_[i] = new TNode(node_name, node_name, object_name,
		bh2SegLocalPos.x(), bh2SegLocalPos.y(), bh2SegLocalPos.z());
  }
  node_->cd();

}

void EvDisp::ConstructBH2Vtx(void)
{
  static const std::string funcname = "EvDisp::ConstructBH2Vtx";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BH2
  double rotMatBh2[9];
  double Bh2WallX = 130.0/2.0;
  double Bh2WallY = 14.0/2.0;
  double Bh2WallZ = 60.0/2.0;

  double Bh2SegX[8] = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  double Bh2SegY[8] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};
  double Bh2SegZ[8] = {40./2., 40./2., 40./2., 40./2., 40./2., 40./2., 40./2., 40./2.};

  double localPosX[8] = {-41.5, -19.0, -10.5, -3.5, 3.5, 10.5, 19.0, 41.5};
  double localPosZ[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

  lnum = geomMan.GetDetectorId("BH2(global)");
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBh2);

  TRotMatrix *rotBh2 = new TRotMatrix("rotBh2","rotBh2",rotMatBh2);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector bh2WallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector bh2WallGlobalPos = geomMan.Local2GlobalPos(lnum, bh2WallLocalPos);
  sprintf(object_name, "Bh2Wall_Brik_");
  Bh2Wall_Brik_vtx_ = new TBRIK(object_name, object_name, "void", 
			    Bh2WallX, Bh2WallY, Bh2WallZ);    
  sprintf(node_name, "Bh2Wall_Node_");
  Bh2Wall_Node_vtx_ = new TNode(node_name, node_name, object_name,
			    bh2WallGlobalPos.x(),
			    bh2WallGlobalPos.y(), 
			    bh2WallGlobalPos.z(),
			    "rotBh2", "void");
  Bh2Wall_Node_vtx_->SetVisibility(0);

  Bh2Wall_Node_vtx_->cd();

  for (int i=0; i<NumOfSegBH2; i++) {
    sprintf(object_name, "Bh2Seg_Brik_vtx_%d", i);
    Bh2Seg_Brik_vtx_[i] = new TBRIK(object_name, object_name, "void", 
				Bh2SegX[i], Bh2SegY[i], Bh2SegZ[i]);    

    ThreeVector bh2SegLocalPos = 
      ThreeVector(localPosX[i], 0.0, localPosZ[i]);
    sprintf(node_name, "Bh2Seg_Node_vtx_%d", i);
    Bh2Seg_Node_vtx_[i] = new TNode(node_name, node_name, object_name,
		bh2SegLocalPos.x(), bh2SegLocalPos.y(), bh2SegLocalPos.z());
  }
  node_vtx_->cd();

}

void EvDisp::ConstructSKS(void)
{

  static const std::string funcname = "EvDisp::ConstructSKS";  

  //-----RotationMatrix
  TRotMatrix *rot = new TRotMatrix("rot","rot",90.0,0.0,0.0,0.0,90.0,-90.0);

  double Yoffset=1024.52;
  //-----Yoke1
  ThreeVector Yoke1Pos(0.0, -1000.0+Yoffset, 0.0); /*mm*/
  double Yoke1Rmin = 2005.0; 
  double Yoke1Rmax = 3405.0; 
  double Yoke1Z    = 500.0/2.0; 
  double Yoke1Phi1 = 44.488;  /*deg*/
  double Yoke1Phi2 = 91.024;  /*deg*/
  
  Yoke1_Tubs_ = new TTUBS("Yoke1_Tubs_", "SKS Yoke1", "void",
			 Yoke1Rmin, Yoke1Rmax, Yoke1Z,
			 Yoke1Phi1, Yoke1Phi1+Yoke1Phi2);

  node_->cd();
  Yoke1_Node_ = new TNode("Yoke1_Node_", "Yoke1_Node", "Yoke1_Tubs_",
			  Yoke1Pos.x(), Yoke1Pos.y(), Yoke1Pos.z());

  //-----Yoke2
  ThreeVector Yoke2Pos(0.0, -1058.0+Yoffset, 0.0); /*mm*/
  double Yoke2dX1 = 1180.0/2.0;
  double Yoke2dX2 = 1655.0/2.0;
  double Yoke2dY  = 500.0/2.0; 
  double Yoke2dZ  = 390.0/2.0;
  
  Yoke2_Trd_ = new TTRD1("Yoke2_Trd_", "SKS Yoke2", "void",
			 Yoke2dX1, Yoke2dX2, Yoke2dY, Yoke2dZ);

  node_->cd();
  Yoke2_Node_ = new TNode("Yoke2_Node_", "Yoke2_Node", "Yoke2_Trd_",
			  Yoke2Pos.x(), Yoke2Pos.y(), Yoke2Pos.z(),
			  "rot", "void");

  //-----Yoke3
  ThreeVector Yoke3Pos(0.0, -1855.5+Yoffset, 0.0); /*mm*/
  double Yoke3X = 1515.0/2.0;
  double Yoke3Y = 1205.0/2.0;
  double Yoke3Z  = 500.0/2.0; 

  Yoke3_Brik_ = new TBRIK("Yoke3_Brik_", "SKS Yoke3", "void",
			 Yoke3X, Yoke3Y, Yoke3Z);

  node_->cd();
  Yoke3_Node_ = new TNode("Yoke3_Node_", "Yoke3_Node", "Yoke3_Brik_",
			  Yoke3Pos.x(), Yoke3Pos.y(), Yoke3Pos.z());

  //-----Yoke4
  ThreeVector Yoke4Pos(0.0, -8.0+Yoffset, -504.0); /*mm*/
  double Yoke4dX1 = 3296.64/2.0;
  double Yoke4dX2 = 1656.0/2.0;
  double Yoke4dY  = 500.0/2.0; 
  double Yoke4dZ  = 2258.0/2.0;
  
  Yoke4_Trd_ = new TTRD1("Yoke4_Trd_", "SKS Yoke4", "void",
			 Yoke4dX1, Yoke4dX2, Yoke4dY, Yoke4dZ);

  node_->cd();
  Yoke4_Node_ = new TNode("Yoke4_Node_", "Yoke4_Node", "Yoke4_Trd_",
			  Yoke4Pos.x(), Yoke4Pos.y(), Yoke4Pos.z(),
			  "rot", "void");

  //-----Yoke5
  ThreeVector Yoke5Pos(0.0, -8.0+Yoffset, 504.0); /*mm*/
  double Yoke5dX1 = 3296.64/2.0;
  double Yoke5dX2 = 1656.0/2.0;
  double Yoke5dY  = 500.0/2.0; 
  double Yoke5dZ  = 2258.0/2.0;
  
  Yoke5_Trd_ = new TTRD1("Yoke5_Trd_", "SKS Yoke5", "void",
			 Yoke5dX1, Yoke5dX2, Yoke5dY, Yoke5dZ);

  node_->cd();
  Yoke5_Node_ = new TNode("Yoke5_Node_", "Yoke5_Node", "Yoke5_Trd_",
			  Yoke5Pos.x(), Yoke5Pos.y(), Yoke5Pos.z(),
			  "rot", "void");

  return;

}
#if 0
void EvDisp::ConstructBdcOut(void)
{
  static const std::string funcname = "EvDisp::ConstructBdcOut";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  const int OffsetGlobal2Local = 30;

  int lnum;
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BDC3X
  double Bdc3xRmin = 0.0; 
  double Bdc3xRmax = 0.01; 
  double Bdc3xZ    = 400.0/2.0; 
  //-- x
  lnum = geomMan.GetDetectorId("BDC3-x-1(global)");
  sprintf(object_name, "Bdc3x_Tube");
  TTUBE *Bdc3x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3x_Node_%d", wire);
    Bdc3x_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- xp
  lnum = geomMan.GetDetectorId("BDC3-x-2(global)");
  sprintf(object_name, "Bdc3xp_Tube");
  TTUBE *Bdc3xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3xp_Node_%d", wire);
    Bdc3xp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-----BDC3V
  double Bdc3vRmin = 0.0; 
  double Bdc3vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  double Bdc3vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3v);
  TRotMatrix *rotBdc3v = new TRotMatrix("rotBdc3v","rotBdc3v",rotMatBdc3v);

  //-- vp
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  sprintf(object_name, "Bdc3vp_Tube");
  TTUBE *Bdc3vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3vp_Node_%d", wire);
    Bdc3vp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3v", "void");
  }

  //-- v
  lnum = geomMan.GetDetectorId("BDC3-v-2(global)");
  sprintf(object_name, "Bdc3v_Tube");
  TTUBE *Bdc3v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3v_Node_%d", wire);
    Bdc3v_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3v", "void");
  }

  //-----BDC3U
  double Bdc3uRmin = 0.0; 
  double Bdc3uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  double Bdc3uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3u);
  TRotMatrix *rotBdc3u = new TRotMatrix("rotBdc3u","rotBdc3u",rotMatBdc3u);

  //-- up
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  sprintf(object_name, "Bdc3up_Tube");
  TTUBE *Bdc3up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3up_Node_%d", wire);
    Bdc3up_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3u", "void");
  }

  //-- u
  lnum = geomMan.GetDetectorId("BDC3-u-2(global)");
  sprintf(object_name, "Bdc3u_Tube");
  TTUBE *Bdc3u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3u_Node_%d", wire);
    Bdc3u_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3u", "void");
  }


  //-----BDC4U
  double Bdc4uRmin = 0.0; 
  double Bdc4uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  double Bdc4uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4u);
  TRotMatrix *rotBdc4u = new TRotMatrix("rotBdc4u","rotBdc4u",rotMatBdc4u);

  //-- u
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  sprintf(object_name, "Bdc4u_Tube");
  TTUBE *Bdc4u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4u_Node_%d", wire);
    Bdc4u_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4u", "void");
  }

  //-- up
  lnum = geomMan.GetDetectorId("BDC4-u-2(global)");
  sprintf(object_name, "Bdc4up_Tube");
  TTUBE *Bdc4up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4up_Node_%d", wire);
    Bdc4up_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4u", "void");
  }

  //-----BDC4V
  double Bdc4vRmin = 0.0; 
  double Bdc4vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  double Bdc4vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4v);
  TRotMatrix *rotBdc4v = new TRotMatrix("rotBdc4v","rotBdc4v",rotMatBdc4v);

  //-- v
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  sprintf(object_name, "Bdc4v_Tube");
  TTUBE *Bdc4v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4v_Node_%d", wire);
    Bdc4v_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4v", "void");
  }

  //-- vp
  lnum = geomMan.GetDetectorId("BDC4-v-2(global)");
  sprintf(object_name, "Bdc4vp_Tube");
  TTUBE *Bdc4vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4vp_Node_%d", wire);
    Bdc4vp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4v", "void");
  }

  //-----BDC4X
  double Bdc4xRmin = 0.0; 
  double Bdc4xRmax = 0.01; 
  double Bdc4xZ    = 400.0/2.0; 

  //-- xp
  lnum = geomMan.GetDetectorId("BDC4-x-1(global)");
  sprintf(object_name, "Bdc4xp_Tube");
  TTUBE *Bdc4xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4xp_Node_%d", wire);
    Bdc4xp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- x
  lnum = geomMan.GetDetectorId("BDC4-x-2(global)");
  sprintf(object_name, "Bdc4x_Tube");
  TTUBE *Bdc4x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4x_Node_%d", wire);
    Bdc4x_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }
}
#endif

#if 0
void EvDisp::ConstructBdcOutVtx(void)
{
  static const std::string funcname = "EvDisp::ConstructBdcOutVtx";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  const int OffsetGlobal2Local = 30;
  
  int lnum;
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BDC3X
  double Bdc3xRmin = 0.0; 
  double Bdc3xRmax = 0.01; 
  double Bdc3xZ    = 400.0/2.0; 
  //-- x
  lnum = geomMan.GetDetectorId("BDC3-x-1(global)");
  sprintf(object_name, "Bdc3x_Tube");
  TTUBE *Bdc3x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3x_Node_vtx_%d", wire);
    Bdc3x_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- xp
  lnum = geomMan.GetDetectorId("BDC3-x-2(global)");
  sprintf(object_name, "Bdc3xp_Tube");
  TTUBE *Bdc3xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3xp_Node_vtx_%d", wire);
    Bdc3xp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-----BDC3V
  double Bdc3vRmin = 0.0; 
  double Bdc3vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  double Bdc3vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3v);
  TRotMatrix *rotBdc3v = new TRotMatrix("rotBdc3v","rotBdc3v",rotMatBdc3v);

  //-- vp
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  sprintf(object_name, "Bdc3vp_Tube");
  TTUBE *Bdc3vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3vp_Node_vtx_%d", wire);
    Bdc3vp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3v", "void");
  }

  //-- v
  lnum = geomMan.GetDetectorId("BDC3-v-2(global)");
  sprintf(object_name, "Bdc3v_Tube");
  TTUBE *Bdc3v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3v_Node_vtx_%d", wire);
    Bdc3v_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3v", "void");
  }

  //-----BDC3U
  double Bdc3uRmin = 0.0; 
  double Bdc3uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  double Bdc3uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3u);
  TRotMatrix *rotBdc3u = new TRotMatrix("rotBdc3u","rotBdc3u",rotMatBdc3u);

  //-- up
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  sprintf(object_name, "Bdc3up_Tube");
  TTUBE *Bdc3up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3up_Node_vtx_%d", wire);
    Bdc3up_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3u", "void");
  }

  //-- u
  lnum = geomMan.GetDetectorId("BDC3-u-2(global)");
  sprintf(object_name, "Bdc3u_Tube");
  TTUBE *Bdc3u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3u_Node_vtx_%d", wire);
    Bdc3u_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3u", "void");
  }


  //-----BDC4U
  double Bdc4uRmin = 0.0; 
  double Bdc4uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  double Bdc4uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4u);
  TRotMatrix *rotBdc4u = new TRotMatrix("rotBdc4u","rotBdc4u",rotMatBdc4u);

  //-- u
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  sprintf(object_name, "Bdc4u_Tube");
  TTUBE *Bdc4u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4u_Node_vtx_%d", wire);
    Bdc4u_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4u", "void");
  }

  //-- up
  lnum = geomMan.GetDetectorId("BDC4-u-2(global)");
  sprintf(object_name, "Bdc4up_Tube");
  TTUBE *Bdc4up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4up_Node_vtx_%d", wire);
    Bdc4up_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4u", "void");
  }

  //-----BDC4V
  double Bdc4vRmin = 0.0; 
  double Bdc4vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  double Bdc4vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4v);
  TRotMatrix *rotBdc4v = new TRotMatrix("rotBdc4v","rotBdc4v",rotMatBdc4v);

  //-- v
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  sprintf(object_name, "Bdc4v_Tube");
  TTUBE *Bdc4v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4v_Node_vtx_%d", wire);
    Bdc4v_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4v", "void");
  }

  //-- vp
  lnum = geomMan.GetDetectorId("BDC4-v-2(global)");
  sprintf(object_name, "Bdc4vp_Tube");
  TTUBE *Bdc4vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4vp_Node_vtx_%d", wire);
    Bdc4vp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4v", "void");
  }

  //-----BDC4X
  double Bdc4xRmin = 0.0; 
  double Bdc4xRmax = 0.01; 
  double Bdc4xZ    = 400.0/2.0; 

  //-- xp
  lnum = geomMan.GetDetectorId("BDC4-x-1(global)");
  sprintf(object_name, "Bdc4xp_Tube");
  TTUBE *Bdc4xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4xp_Node_vtx_%d", wire);
    Bdc4xp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- x
  lnum = geomMan.GetDetectorId("BDC4-x-2(global)");
  sprintf(object_name, "Bdc4x_Tube");
  TTUBE *Bdc4x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4x_Node_vtx_%d", wire);
    Bdc4x_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }
}
#endif

#if 1
void EvDisp::ConstructSdcIn(void)
{
  static const std::string funcname = "EvDisp::ConstructSdcIn";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  
  int lnum; /* layer number */
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

#if 0
  //-----SDC1
  //-- u
  lnum = 1;
  double Sdc1uRmin = 0.0; 
  double Sdc1uRmax = 0.01; 
  double Sdc1uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u1);
  TRotMatrix *rotSdc1u1 = new TRotMatrix("rotSdc1u1","rotSdc1u1",rotMatSdc1u1);
  
  sprintf(object_name, "Sdc1u_Tube");
  TTUBE *Sdc1u_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1uRmin, Sdc1uRmax, Sdc1uZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u1_Node_%d", wire);
    Sdc1u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				     wireGlobalPos.x(), 
				     wireGlobalPos.y(), 
				     wireGlobalPos.z(),
				     "rotSdc1u1", "void");
  }
  //-- u2
  lnum = 2;
  double rotMatSdc1u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u2);
  TRotMatrix *rotSdc1u2 = new TRotMatrix("rotSdc1u2","rotSdc1u2",rotMatSdc1u2);
  
  for (wire=1; wire <= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u2_Node_%d", wire);
    Sdc1u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc1u2", "void");
  }

  //-- v1
  lnum = 3;
  double Sdc1vRmin = 0.0; 
  double Sdc1vRmax = 0.01; 
  double Sdc1vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v1);
  TRotMatrix *rotSdc1v1 = new TRotMatrix("rotSdc1v1","rotSdc1v1",rotMatSdc1v1);
  
  sprintf(object_name, "Sdc1v_Tube");
  TTUBE *Sdc1v_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1vRmin, Sdc1vRmax, Sdc1vZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v1_Node_%d", wire);
    Sdc1v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				     wireGlobalPos.x(), 
				     wireGlobalPos.y(), 
				     wireGlobalPos.z(),
				     "rotSdc1v1", "void");
  }
  //-- v2
  lnum = 4;
  double rotMatSdc1v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v2);
  TRotMatrix *rotSdc1v2 = new TRotMatrix("rotSdc1v2","rotSdc1v2",rotMatSdc1v2);
  
  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v2_Node_%d", wire);
    Sdc1v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc1v2", "void");
  }
#endif

  // SFT
  //-- X
  double SftXRmin = 0.0; 
  double SftXRmax = 0.5; 
  double SftXZ    = 200.0/2.0; 

  sprintf(object_name, "SftX_Tube");
  TTUBE *SftX_Tube = new TTUBE(object_name, object_name, "void", 
			  SftXRmin, SftXRmax, SftXZ);

  lnum = 2;
  double rotMatSftX[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSftX);
  TRotMatrix *rotSftX = new TRotMatrix("rotSftX","rotSftX",rotMatSftX);
  
  for (wire=1; wire <= NumOfSegSFT_X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "SftX_Node_%d", wire);
    SftX_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSftX", "void");
  }

  //-- V
  lnum = 3;
  double SftVRmin = 0.0; 
  double SftVRmax = 0.375; 
  double SftVZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSftV[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSftV);
  TRotMatrix *rotSftV = new TRotMatrix("rotSftV","rotSftV",rotMatSftV);
  
  sprintf(object_name, "SftV_Tube");
  TTUBE *SftV_Tube = new TTUBE(object_name, object_name, "void", 
				 SftVRmin, SftVRmax, SftVZ);

  for (wire=1; wire<= NumOfSegSFT_UV; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "SftV_Node_%d", wire);
    SftV_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				     wireGlobalPos.x(), 
				     wireGlobalPos.y(), 
				     wireGlobalPos.z(),
				     "rotSftV", "void");
  }
  //-- v2
  lnum = 4;
  double rotMatSftU[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSftU);
  TRotMatrix *rotSftU = new TRotMatrix("rotSftU","rotSftU",rotMatSftU);
  
  for (wire=1; wire<= NumOfSegSFT_UV; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "SftU_Node_%d", wire);
    SftU_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSftU", "void");
  }

  //-----SDC2
  //-- v
  lnum = 5;
  double Sdc2vRmin = 0.0; 
  double Sdc2vRmax = 0.01; 
  double Sdc2vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v1);
  TRotMatrix *rotSdc2v1 = new TRotMatrix("rotSdc2v1","rotSdc2v1",rotMatSdc2v1);

  sprintf(object_name, "Sdc2v1_Tube");
  TTUBE *Sdc2v_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2vRmin, Sdc2vRmax, Sdc2vZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v1_Node_%d", wire);
    Sdc2v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2v1", "void");
  }
  // v2
  lnum = 6;
  double rotMatSdc2v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v2);
  TRotMatrix *rotSdc2v2 = new TRotMatrix("rotSdc2v2","rotSdc2v2",rotMatSdc2v2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v2_Node_%d", wire);
    Sdc2v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2v2", "void");
  }

  //-- u1
  lnum = 7;
  double Sdc2uRmin = 0.0; 
  double Sdc2uRmax = 0.01; 
  double Sdc2uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u1);
  TRotMatrix *rotSdc2u1 = new TRotMatrix("rotSdc2u1","rotSdc2u1",rotMatSdc2u1);

  sprintf(object_name, "Sdc2u1_Tube");
  TTUBE *Sdc2u_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2uRmin, Sdc2uRmax, Sdc2uZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u1_Node_%d", wire);
    Sdc2u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2u1", "void");
  }
  // u2
  lnum = 8;
  double rotMatSdc2u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u2);
  TRotMatrix *rotSdc2u2 = new TRotMatrix("rotSdc2u2","rotSdc2u2",rotMatSdc2u2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u2_Node_%d", wire);
    Sdc2u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2u2", "void");
  }

  double Sdc2xRmin = 0.0; 
  double Sdc2xRmax = 0.01; 
  double Sdc2xZ    = 200.0/2.0; 
  //-- x1
  lnum = 9;

  sprintf(object_name, "Sdc2x_Tube");
  TTUBE *Sdc2x_Tube = new TTUBE(object_name, object_name, "void", 
			  Sdc2xRmin, Sdc2xRmax, Sdc2xZ);
  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x1_Node_%d", wire);
    Sdc2x1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z());
  }

  //-- x2
  lnum = 10;

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x2_Node_%d", wire);
    Sdc2x2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

}
#endif 



 void EvDisp::ConstructSdcInVtx(void)
{
  static const std::string funcname = "EvDisp::ConstructSdcInVtx";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  
  int lnum; /* layer number */
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

#if 0
  //-----SDC1
  //-- u
  lnum = 1;
  double Sdc1uRmin = 0.0; 
  double Sdc1uRmax = 0.01; 
  double Sdc1uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u1);
  TRotMatrix *rotSdc1u1 = new TRotMatrix("rotSdc1u1","rotSdc1u1",rotMatSdc1u1);
  
  sprintf(object_name, "Sdc1u_Tube");
  TTUBE *Sdc1u_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1uRmin, Sdc1uRmax, Sdc1uZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u1_Node_vtx_%d", wire);
    Sdc1u1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1u1", "void");
  }
  //-- u2
  lnum = 2;
  double rotMatSdc1u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u2);
  TRotMatrix *rotSdc1u2 = new TRotMatrix("rotSdc1u2","rotSdc1u2",rotMatSdc1u2);
  
  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u2_Node_vtx_%d", wire);
    Sdc1u2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1u2", "void");
  }

  //-- v1
  lnum = 3;
  double Sdc1vRmin = 0.0; 
  double Sdc1vRmax = 0.01; 
  double Sdc1vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v1);
  TRotMatrix *rotSdc1v1 = new TRotMatrix("rotSdc1v1","rotSdc1v1",rotMatSdc1v1);
  
  sprintf(object_name, "Sdc1v_Tube");
  TTUBE *Sdc1v_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1vRmin, Sdc1vRmax, Sdc1vZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v1_Node_vtx_%d", wire);
    Sdc1v1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1v1", "void");
  }
  //-- v2
  lnum = 4;
  double rotMatSdc1v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v2);
  TRotMatrix *rotSdc1v2 = new TRotMatrix("rotSdc1v2","rotSdc1v2",rotMatSdc1v2);
  
  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v2_Node_vtx_%d", wire);
    Sdc1v2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1v2", "void");
  }
#endif

  // SFT
  //-- X
  double SftXRmin = 0.0; 
  double SftXRmax = 0.5; 
  double SftXZ    = 200.0/2.0; 

  sprintf(object_name, "SftX_Tube");
  TTUBE *SftX_Tube = new TTUBE(object_name, object_name, "void", 
			  SftXRmin, SftXRmax, SftXZ);

  lnum = 2;
  double rotMatSftX[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSftX);
  TRotMatrix *rotSftX = new TRotMatrix("rotSftX","rotSftX",rotMatSftX);
  
  for (wire=1; wire <= NumOfSegSFT_X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "SftX_Node_vtx_%d", wire);
    SftX_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSftX", "void");
  }

  //-- V
  lnum = 3;
  double SftVRmin = 0.0; 
  double SftVRmax = 0.375; 
  double SftVZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSftV[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSftV);
  TRotMatrix *rotSftV = new TRotMatrix("rotSftV","rotSftV",rotMatSftV);
  
  sprintf(object_name, "SftV_Tube");
  TTUBE *SftV_Tube = new TTUBE(object_name, object_name, "void", 
				 SftVRmin, SftVRmax, SftVZ);

  for (wire=1; wire<= NumOfSegSFT_UV; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "SftV_Node_vtx_%d", wire);
    SftV_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				     wireGlobalPos.x(), 
				     wireGlobalPos.y(), 
				     wireGlobalPos.z(),
				     "rotSftV", "void");
  }
  //-- U
  lnum = 4;
  double rotMatSftU[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSftU);
  TRotMatrix *rotSftU = new TRotMatrix("rotSftU","rotSftU",rotMatSftU);
  
  for (wire=1; wire<= NumOfSegSFT_UV; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "SftU_Node_vtx_%d", wire);
    SftU_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSftU", "void");

  }

  //-----SDC2
  //-- v
  lnum = 5;
  double Sdc2vRmin = 0.0; 
  double Sdc2vRmax = 0.01; 
  double Sdc2vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v1);
  TRotMatrix *rotSdc2v1 = new TRotMatrix("rotSdc2v1","rotSdc2v1",rotMatSdc2v1);

  sprintf(object_name, "Sdc2v1_Tube");
  TTUBE *Sdc2v_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2vRmin, Sdc2vRmax, Sdc2vZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v1_Node_vtx_%d", wire);
    Sdc2v1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2v1", "void");
  }
  // v2
  lnum = 6;
  double rotMatSdc2v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v2);
  TRotMatrix *rotSdc2v2 = new TRotMatrix("rotSdc2v2","rotSdc2v2",rotMatSdc2v2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v2_Node_vtx_%d", wire);
    Sdc2v2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2v2", "void");
  }

  //-- u1
  lnum = 7;
  double Sdc2uRmin = 0.0; 
  double Sdc2uRmax = 0.01; 
  double Sdc2uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u1);
  TRotMatrix *rotSdc2u1 = new TRotMatrix("rotSdc2u1","rotSdc2u1",rotMatSdc2u1);

  sprintf(object_name, "Sdc2u1_Tube");
  TTUBE *Sdc2u_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2uRmin, Sdc2uRmax, Sdc2uZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u1_Node_vtx_%d", wire);
    Sdc2u1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2u1", "void");
  }
  // u2
  lnum = 8;
  double rotMatSdc2u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u2);
  TRotMatrix *rotSdc2u2 = new TRotMatrix("rotSdc2u2","rotSdc2u2",rotMatSdc2u2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u2_Node_vtx_%d", wire);
    Sdc2u2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2u2", "void");
  }

  double Sdc2xRmin = 0.0; 
  double Sdc2xRmax = 0.01; 
  double Sdc2xZ    = 200.0/2.0; 
  //-- x1
  lnum = 9;

  sprintf(object_name, "Sdc2x_Tube");
  TTUBE *Sdc2x_Tube = new TTUBE(object_name, object_name, "void", 
			  Sdc2xRmin, Sdc2xRmax, Sdc2xZ);
  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x1_Node_vtx_%d", wire);
    Sdc2x1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z());
  }

  //-- x2
  lnum = 10;

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x2_Node_vtx_%d", wire);
    Sdc2x2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z());
  }

}

void EvDisp::ConstructBcIn(void)
{
  static const std::string funcname = "EvDisp::ConstructBcIn";  
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  int wire;
  double localZ, localPos;
  double tiltAngle;

  for (lnum=101; lnum<=112; lnum++) {
    localZ = geomMan.GetLocalZ(lnum);
    tiltAngle = geomMan.GetTiltAngle(lnum);
    for (wire=1; wire<= MaxWireBC1; wire++) {
      localPos = geomMan.calcWirePosition(lnum, wire)/cos(tiltAngle*Deg2Rad);
      
      BcIn_Arc_[lnum-101][wire-1] = new TArc(localPos, localZ, 0.5);
      //BcIn_Arc_[lnum-101][wire-1]->Draw("same");
    }
  }

  double Bh1SegX[NumOfSegBH1] = {30./2., 20./2., 16./2., 12./2., 8./2., 8./2., 8./2., 12./2., 16./2., 20./2., 30./2.};
  double Bh1SegY[NumOfSegBH1] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  //double localPosBh1Z = -515.;
  double localPosBh1Z = DCGeomMan::GetInstance().GetLocalZ( IdBH1 );
  double localPosBh1X_dX = 0.;
  double localPosBh1X[NumOfSegBH1] = {-70. + localPosBh1X_dX,
				      -46. + localPosBh1X_dX,
				      -29. + localPosBh1X_dX,
				      -16. + localPosBh1X_dX,
				      -7. + localPosBh1X_dX,
				      0. + localPosBh1X_dX, 
				      7. + localPosBh1X_dX,
				      16. + localPosBh1X_dX,
				      29. + localPosBh1X_dX,
				      46. + localPosBh1X_dX,
				      70. + localPosBh1X_dX};
  double localPosBh1_dZ[NumOfSegBH1] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};


  for (int i=0; i<NumOfSegBH1; i++) {
    Bh1_Line_[i][0] = new TLine(localPosBh1X[i]-Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i],
				localPosBh1X[i]+Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i]);
    Bh1_Line_[i][1] = new TLine(localPosBh1X[i]+Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i],
				localPosBh1X[i]+Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i]);
    Bh1_Line_[i][2] = new TLine(localPosBh1X[i]+Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i],
				localPosBh1X[i]-Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i]);
    Bh1_Line_[i][3] = new TLine(localPosBh1X[i]+Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]+Bh1SegY[i],
				localPosBh1X[i]+Bh1SegX[i], 
				localPosBh1Z+localPosBh1_dZ[i]-Bh1SegY[i]);
    for (int j=0; j<4; j++)
      Bh1_Line_[i][j]->Draw("same");

  }


  //double localPosBftZ = -160.;  
  double localPosBftZ=DCGeomMan::GetInstance().GetLocalZ( IdBFT );
  double localPos_dX = 0;
  double localPos_dZ = sqrt(3.)/2./2.;
  for (int i=1; i<=NumOfSegBFT; i++) {
    localPos = ((double)i-80.5)*(1.);

    Bft_Arc_[0][i-1] = new TArc(localPos+localPos_dX, localPosBftZ-localPos_dZ, 0.5);
    Bft_Arc_[0][i-1]->Draw("same");

    localPos = ((double)i-80)*(1.);

    Bft_Arc_[1][i-1] = new TArc(localPos+localPos_dX, localPosBftZ+localPos_dZ, 0.5);
    Bft_Arc_[1][i-1]->Draw("same");

    Bft_Arc_[0][i-1]->SetLineColor(kWhite);
    Bft_Arc_[1][i-1]->SetLineColor(kWhite);
  }
  Bft_Arc_[0][80]->SetLineColor(kBlack);
  Bft_Arc_[0][80]->SetFillColor(kRed);

}


void EvDisp::ConstructBcOut(void)
{
  static const std::string funcname = "EvDisp::ConstructBcOut";  
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  int wire;
  double localZ, localPos;
  double tiltAngle;
  double localZ_K18=geomMan.GetLocalZ(130);

  tp_bcOut1_->cd();

  for (lnum=113; lnum<=118; lnum++) {
    localZ = geomMan.GetLocalZ(lnum) - localZ_K18;
    tiltAngle = geomMan.GetTiltAngle(lnum);
    for (wire=1; wire<= MaxWireBC3; wire++) {
      localPos = geomMan.calcWirePosition(lnum, wire)/cos(tiltAngle*Deg2Rad);
      
      Bc3_Arc_[lnum-113][wire-1] = new TArc(localPos, localZ, 1.5);
      Bc3_Arc_[lnum-113][wire-1]->Draw("same");
    }
  }

  for (lnum=119; lnum<=124; lnum++) {
    localZ = geomMan.GetLocalZ(lnum) - localZ_K18;
    tiltAngle = geomMan.GetTiltAngle(lnum);
    for (wire=1; wire<= MaxWireBC4; wire++) {
      localPos = geomMan.calcWirePosition(lnum, wire)/cos(tiltAngle*Deg2Rad);
      
      Bc4_Arc_[lnum-119][wire-1] = new TArc(localPos, localZ, 1.5);
      Bc4_Arc_[lnum-119][wire-1]->Draw("same");
    }
  }

  /*
  for (lnum=1; lnum<=4; lnum++) {
    localZ = geomMan.GetLocalZ(lnum);
    tiltAngle = geomMan.GetTiltAngle(lnum);
    for (wire=1; wire<= MaxWireSDC1; wire++) {
      localPos = geomMan.calcWirePosition(lnum, wire)/cos(tiltAngle*Deg2Rad);
      
      Sdc1_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ, 0.5);
      Sdc1_Arc_[lnum-1][wire-1]->Draw("same");
    }
  }
  */

  for (lnum=1; lnum<=2; lnum++) {
    localZ = geomMan.GetLocalZ(2);
    tiltAngle = geomMan.GetTiltAngle(2);

    for (wire=1; wire<= NumOfSegSFT_X; wire++) {
      if (lnum==1) {
	localPos = geomMan.calcWirePosition(2, wire)/cos(tiltAngle*Deg2Rad);
	//Sft_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ-0.5, 0.5);
	Sft_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ-0.5, 1.5);
	Sft_Arc_[lnum-1][wire-1]->Draw("same");
      } else {
	localPos = geomMan.calcWirePosition(2, wire)/cos(tiltAngle*Deg2Rad)+0.5;
	//Sft_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ+0.5, 0.5);
	Sft_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ+0.5, 1.5);
	Sft_Arc_[lnum-1][wire-1]->Draw("same");
      }
    }
  }

  for (lnum=2; lnum<=4; lnum++) {
    localZ = geomMan.GetLocalZ(lnum);
    tiltAngle = geomMan.GetTiltAngle(lnum);

    for (wire=1; wire<= NumOfSegSFT_UV; wire++) {
      localPos = geomMan.calcWirePosition(lnum, wire)/cos(tiltAngle*Deg2Rad);
      //Sft_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ, 0.375);
      Sft_Arc_[lnum-1][wire-1] = new TArc(localPos, localZ, 1.375);
      Sft_Arc_[lnum-1][wire-1]->Draw("same");
    }
  }

  for (lnum=5; lnum<=10; lnum++) {
    localZ = geomMan.GetLocalZ(lnum);
    tiltAngle = geomMan.GetTiltAngle(lnum);
    for (wire=1; wire<= MaxWireSDC2; wire++) {
      localPos = geomMan.calcWirePosition(lnum, wire)/cos(tiltAngle*Deg2Rad);
      
      Sdc2_Arc_[lnum-5][wire-1] = new TArc(localPos, localZ, 1.5);
      Sdc2_Arc_[lnum-5][wire-1]->Draw("same");
    }
  }

  double Bh2SegX[NumOfSegBH2] = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
  double Bh2SegY[NumOfSegBH2] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};

  double localPosBh2Z = 604.63-1800.;
  double localPosBh2X_dX = 0.;

  double localPosBh2X[NumOfSegBH2] = {-41.5 + localPosBh2X_dX,
				      -19.0 + localPosBh2X_dX,
				      -10.5 + localPosBh2X_dX,
				      -3.5 + localPosBh2X_dX,
				      3.5 + localPosBh2X_dX,
				      10.5 + localPosBh2X_dX,
				      19.0 + localPosBh2X_dX,
                                      41.5 + localPosBh2X_dX};

  double localPosBh2_dZ[NumOfSegBH2] = {0., 0., 0., 0., 0., 0., 0., 0.};


  for (int i=0; i<NumOfSegBH2; i++) {
    Bh2_Line_[i][0] = new TLine(localPosBh2X[i]-Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]-Bh2SegY[i],
				localPosBh2X[i]+Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]-Bh2SegY[i]);
    Bh2_Line_[i][1] = new TLine(localPosBh2X[i]+Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]-Bh2SegY[i],
				localPosBh2X[i]+Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]+Bh2SegY[i]);
    Bh2_Line_[i][2] = new TLine(localPosBh2X[i]+Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]+Bh2SegY[i],
				localPosBh2X[i]-Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]+Bh2SegY[i]);
    Bh2_Line_[i][3] = new TLine(localPosBh2X[i]+Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]+Bh2SegY[i],
				localPosBh2X[i]+Bh2SegX[i], 
				localPosBh2Z+localPosBh2_dZ[i]-Bh2SegY[i]);
    for (int j=0; j<4; j++)
      Bh2_Line_[i][j]->Draw("same");

  }

  double target_width=70.;
  double target_height=40.;
  double target_thickness=77.;
  Target_Box_ = new TBox(-target_width/2., -target_thickness/2,
			 target_width/2., target_thickness/2);
  Target_Box_->SetFillStyle(0);
  Target_Box_->SetLineColor(kBlack);
  Target_Box_->Draw("same");

  tp_bcOut2_->cd();

  Target_Box_Y_ = new TBox(-target_height/2., -target_thickness/2,
			 target_height/2., target_thickness/2);
  Target_Box_Y_->SetFillStyle(0);
  Target_Box_Y_->SetLineColor(kBlack);
  Target_Box_Y_->Draw("same");

}

#if 1
void EvDisp::ConstructSdcOut(void)
{
  static const std::string funcname = "EvDisp::ConstructSdcOut";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  
  int lnum; /* layer number */
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----SDC3
  //-- v1
  double Sdc3vRmin = 0.0; 
  double Sdc3vRmax = 0.01; 
  double Sdc3vZ    = 1240.0/2.0; 
  double rotMatSdc3v1[9];
  lnum = 31;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3v1);
  TRotMatrix *rotSdc3v1 = new TRotMatrix("rotSdc3v1","rotSdc3v1",rotMatSdc3v1);

  sprintf(object_name, "Sdc3v1_Tube");
  TTUBE *Sdc3v1_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc3vRmin, Sdc3vRmax, Sdc3vZ);

  for (wire=1; wire<= MaxWireSDC3V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3v1_Node_%d", wire);
    Sdc3v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3v1", "void");
  }

  //-- x1
  double Sdc3xRmin = 0.0; 
  double Sdc3xRmax = 0.01; 
  double Sdc3xZ    = 1240.0/2.0; 
  lnum = 32;
  sprintf(object_name, "Sdc3x1_Tube");
  TTUBE *Sdc3x1_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc3xRmin, Sdc3xRmax, Sdc3xZ);

  for (wire=1; wire<= MaxWireSDC3X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3x1_Node_%d", wire);
    Sdc3x1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u1
  double Sdc3uRmin = 0.0; 
  double Sdc3uRmax = 0.01; 
  double Sdc3uZ    = 1240.0/2.0; 
  double rotMatSdc3u1[9];
  lnum = 33;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3u1);
  TRotMatrix *rotSdc3u1 = new TRotMatrix("rotSdc3u1","rotSdc3u1",rotMatSdc3u1);

  sprintf(object_name, "Sdc3u1_Tube");
  TTUBE *Sdc3u1_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc3uRmin, Sdc3uRmax, Sdc3uZ);

  for (wire=1; wire<= MaxWireSDC3U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3u1_Node_%d", wire);
    Sdc3u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3u1", "void");
  }

  //-- v2
  double rotMatSdc3v2[9];
  lnum = 34;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3v2);
  TRotMatrix *rotSdc3v2 = new TRotMatrix("rotSdc3v2","rotSdc3v2",rotMatSdc3v2);

  sprintf(object_name, "Sdc3v2_Tube");
  TTUBE *Sdc3v2_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc3vRmin, Sdc3vRmax, Sdc3vZ);

  for (wire=1; wire<= MaxWireSDC3V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3v2_Node_%d", wire);
    Sdc3v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3v2", "void");
  }

  //-- x2
  lnum = 35;
  sprintf(object_name, "Sdc3x2_Tube");
  TTUBE *Sdc3x2_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc3xRmin, Sdc3xRmax, Sdc3xZ);

  for (wire=1; wire<= MaxWireSDC3X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3x2_Node_%d", wire);
    Sdc3x2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u2
  double rotMatSdc3u2[9];
  lnum = 36;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3u2);
  TRotMatrix *rotSdc3u2 = new TRotMatrix("rotSdc3u2","rotSdc3u2",rotMatSdc3u2);
  sprintf(object_name, "Sdc3u2_Tube");
  TTUBE *Sdc3u2_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc3uRmin, Sdc3uRmax, Sdc3uZ);

  for (wire=1; wire<= MaxWireSDC3U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3u2_Node_%d", wire);
    Sdc3u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3u2", "void");
  }

  //-----SDC4
  //-- v1
  double Sdc4vRmin = 0.0; 
  double Sdc4vRmax = 0.01; 
  double Sdc4vZ    = 1240.0/2.0; 
  double rotMatSdc4v1[9];
  lnum = 37;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4v1);
  TRotMatrix *rotSdc4v1 = new TRotMatrix("rotSdc4v1","rotSdc4v1",rotMatSdc4v1);

  sprintf(object_name, "Sdc4v1_Tube");
  TTUBE *Sdc4v1_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc4vRmin, Sdc4vRmax, Sdc4vZ);

  for (wire=1; wire<= MaxWireSDC4V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4v1_Node_%d", wire);
    Sdc4v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4v1", "void");
  }

  //-- x1
  double Sdc4xRmin = 0.0; 
  double Sdc4xRmax = 0.01; 
  double Sdc4xZ    = 1240.0/2.0; 
  lnum = 38;
  sprintf(object_name, "Sdc4x1_Tube");
  TTUBE *Sdc4x1_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc4xRmin, Sdc4xRmax, Sdc4xZ);

  for (wire=1; wire<= MaxWireSDC4X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4x1_Node_%d", wire);
    Sdc4x1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u1
  double Sdc4uRmin = 0.0; 
  double Sdc4uRmax = 0.01; 
  double Sdc4uZ    = 1240.0/2.0; 
  double rotMatSdc4u1[9];
  lnum = 39;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4u1);
  TRotMatrix *rotSdc4u1 = new TRotMatrix("rotSdc4u1","rotSdc4u1",rotMatSdc4u1);

  sprintf(object_name, "Sdc4u1_Tube");
  TTUBE *Sdc4u1_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc4uRmin, Sdc4uRmax, Sdc4uZ);

  for (wire=1; wire<= MaxWireSDC4U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4u1_Node_%d", wire);
    Sdc4u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4u1", "void");
  }

  //-- v2
  double rotMatSdc4v2[9];
  lnum = 40;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4v2);
  TRotMatrix *rotSdc4v2 = new TRotMatrix("rotSdc4v2","rotSdc4v2",rotMatSdc4v2);

  sprintf(object_name, "Sdc4v2_Tube");
  TTUBE *Sdc4v2_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc4vRmin, Sdc4vRmax, Sdc4vZ);

  for (wire=1; wire<= MaxWireSDC4V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4v2_Node_%d", wire);
    Sdc4v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4v2", "void");
  }

  //-- x2
  lnum = 41;
  sprintf(object_name, "Sdc4x2_Tube");
  TTUBE *Sdc4x2_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc4xRmin, Sdc4xRmax, Sdc4xZ);

  for (wire=1; wire<= MaxWireSDC4X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4x2_Node_%d", wire);
    Sdc4x2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u1
  double rotMatSdc4u2[9];
  lnum = 42;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4u2);
  TRotMatrix *rotSdc4u2 = new TRotMatrix("rotSdc4u2","rotSdc4u2",rotMatSdc4u2);
  sprintf(object_name, "Sdc4u2_Tube");
  TTUBE *Sdc4u2_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc4uRmin, Sdc4uRmax, Sdc4uZ);

  for (wire=1; wire<= MaxWireSDC4U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4u2_Node_%d", wire);
    Sdc4u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4u2", "void");
  }

}
#endif

#if 1
void EvDisp::ConstructTOF(void)
{
  static const std::string funcname = "EvDisp::ConstructTOF";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----TOF
  double rotMatTof[9];
  double TofWallX = 70.0*32.0/2.0;
  double TofWallY = 30.0/2.0;
  double TofWallZ = 1000.0/2.0;

  double TofSegX = 70.0/2.0;
  double TofSegY = 30.0/2.0;
  double TofSegZ = 1000.0/2.0;

  lnum = 51;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatTof);
  TRotMatrix *rotTof = new TRotMatrix("rotTof","rotTof",rotMatTof);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector tofWallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector tofWallGlobalPos = geomMan.Local2GlobalPos(lnum, tofWallLocalPos);
  sprintf(object_name, "TofWall_Brik_");
  TofWall_Brik_ = new TBRIK(object_name, object_name, "void", 
			    TofWallX, TofWallY, TofWallZ);    
  sprintf(node_name, "TofWall_Node_");
  TofWall_Node_ = new TNode(node_name, node_name, object_name,
	    tofWallGlobalPos.x(), tofWallGlobalPos.y(), tofWallGlobalPos.z(),
	    "rotTof", "void");

  TofWall_Node_->cd();

  sprintf(object_name, "TofSeg_Brik_");
  TofSeg_Brik_ = new TBRIK(object_name, object_name, "void", 
			   TofSegX, TofSegY, TofSegZ);    
  for (int i=0; i<NumOfSegTOF; i++) {
    ThreeVector tofSegLocalPos = 
      ThreeVector((double)(i-14.5)*TofSegX*2.0, 0.0, 0.0);
    sprintf(node_name, "TofSeg_Node_%d", i);
    TofSeg_Node_[i] = new TNode(node_name, node_name, object_name,
		tofSegLocalPos.x(), tofSegLocalPos.y(), tofSegLocalPos.z());
  }
  node_->cd();
#if 0
  //-----AC1
  double rotMatAc1[9];
  double Ac1X = 1050.0/2.0;
  double Ac1Y = 400.0/2.0;
  double Ac1Z = 1200.0/2.0;

  lnum = 52;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatAc1);
  TRotMatrix *rotAc1 = new TRotMatrix("rotAc1","rotAc1",rotMatAc1);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector Ac1LocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector Ac1GlobalPos = geomMan.Local2GlobalPos(lnum, Ac1LocalPos);

  sprintf(object_name, "Ac1_Brik_");
  Ac1_Brik_ = new TBRIK(object_name, object_name, "void", Ac1X, Ac1Y, Ac1Z);
  sprintf(node_name, "Ac1_Node_");
  Ac1_Node_ = new TNode(node_name, node_name, object_name,
	    Ac1GlobalPos.x(), Ac1GlobalPos.y(), Ac1GlobalPos.z(),
	    "rotAc1", "void");

  //-----AC2
  double rotMatAc2[9];
  double Ac2X = 1400.0/2.0;
  double Ac2Y = 400.0/2.0;
  double Ac2Z = 1400.0/2.0;

  lnum = 53;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatAc2);
  TRotMatrix *rotAc2 = new TRotMatrix("rotAc2","rotAc2",rotMatAc2);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector Ac2LocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector Ac2GlobalPos = geomMan.Local2GlobalPos(lnum, Ac2LocalPos);

  sprintf(object_name, "Ac2_Brik_");
  Ac2_Brik_ = new TBRIK(object_name, object_name, "void", Ac2X, Ac2Y, Ac2Z);
  sprintf(node_name, "Ac2_Node_");
  Ac2_Node_ = new TNode(node_name, node_name, object_name,
	    Ac2GlobalPos.x(), Ac2GlobalPos.y(), Ac2GlobalPos.z(),
	    "rotAc2", "void");

#endif
  //-----LC
  double rotMatLc[9];
  double LcWallX = 100.0*28.0/2.0;
  double LcWallY = 40.0/2.0;
  double LcWallZ = 1400.0/2.0;

  double LcSegX = 100.0/2.0;
  double LcSegY = 40.0/2.0;
  double LcSegZ = 1400.0/2.0;

  lnum = 54;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatLc);
  TRotMatrix *rotLc = new TRotMatrix("rotLc","rotLc",rotMatLc);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector lcWallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector lcWallGlobalPos = geomMan.Local2GlobalPos(lnum, lcWallLocalPos);

  sprintf(object_name, "LcWall_Brik_");
  LcWall_Brik_ = new TBRIK(object_name, object_name, "void", 
			    LcWallX, LcWallY, LcWallZ);    
  sprintf(node_name, "LcWall_Node_");
  LcWall_Node_ = new TNode(node_name, node_name, object_name,
	    lcWallGlobalPos.x(), lcWallGlobalPos.y(), lcWallGlobalPos.z(),
	    "rotLc", "void");

  LcWall_Node_->cd();

  sprintf(object_name, "LcSeg_Brik_");
  LcSeg_Brik_ = new TBRIK(object_name, object_name, "void", 
			  LcSegX, LcSegY, LcSegZ);    
  for (int i=0; i<NumOfSegLC; i++) {
    ThreeVector lcSegLocalPos = 
      ThreeVector((double)(i-13.5)*LcSegX*2.0, 0.0, 0.0);

    sprintf(node_name, "LcSeg_Node_%d", i);
    LcSeg_Node_[i] = new TNode(node_name, node_name, object_name,
		lcSegLocalPos.x(), lcSegLocalPos.y(), lcSegLocalPos.z());
  }
  node_->cd();

}
#endif

#if 0
void EvDisp::DrawInitTrack(int nStep, ThreeVector *StepPoint) const
{
  if ( InitStepMark_ )
    delete InitStepMark_;

  InitStepMark_ = new TPolyMarker3D(nStep);
  for (int i=0; i<nStep; i++) {
    InitStepMark_->SetPoint(i, StepPoint[i].x(), StepPoint[i].y(),  StepPoint[i].z());
  }
  InitStepMark_->SetMarkerSize(1);
  InitStepMark_->SetMarkerColor(kCyan);
  InitStepMark_->SetMarkerStyle(6);

  tp_->cd();
  InitStepMark_->Draw();

  tc_->Update();
}

void EvDisp::DrawInitTrack(void) const
{
  tp_->cd();
  if (InitStepMark_)
    InitStepMark_->Draw();

  tc_->Update();

}
#endif

void EvDisp::DrawHitWire(int lnum, int hit_wire, bool range_check, bool tdc_check) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

    switch (lnum) {
    // SDC1,2 1-10
    case 1:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1u1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1u1_Node_vtx_%d", hit_wire);
      break;
    case 2:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1u2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1u2_Node_vtx_%d", hit_wire);
      break;
    case 3:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1v1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1v1_Node_vtx_%d", hit_wire);
      break;
    case 4:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1v2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1v2_Node_vtx_%d", hit_wire);
      break;
    case 5:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2v1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2v1_Node_vtx_%d", hit_wire);
      break;
    case 6:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2v2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2v2_Node_vtx_%d", hit_wire);
      break;
    case 7:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2u1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2u1_Node_vtx_%d", hit_wire);
      break;
    case 8:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2u2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2u2_Node_vtx_%d", hit_wire);
      break;
    case 9:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2x1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2x1_Node_vtx_%d", hit_wire);
      break;
    case 10:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2x2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2x2_Node_vtx_%d", hit_wire);
      break;
    // SDC3,4 31-42
    case 31:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3V)
	return;
      sprintf(node_name, "Sdc3v1_Node_%d", hit_wire);
      break;
    case 32:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3X)
	return;
      sprintf(node_name, "Sdc3x1_Node_%d", hit_wire);
      break;
    case 33:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3U)
	return;
      sprintf(node_name, "Sdc3u1_Node_%d", hit_wire);
      break;
    case 34:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3V)
	return;
      sprintf(node_name, "Sdc3v2_Node_%d", hit_wire);
      break;
    case 35:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3X)
	return;
      sprintf(node_name, "Sdc3x2_Node_%d", hit_wire);
      break;
    case 36:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3U)
	return;
      sprintf(node_name, "Sdc3u2_Node_%d", hit_wire);
      break;
    case 37:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4V)
	return;
      sprintf(node_name, "Sdc4v1_Node_%d", hit_wire);
      break;
    case 38:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4X)
	return;
      sprintf(node_name, "Sdc4x1_Node_%d", hit_wire);
      break;
    case 39:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4U)
	return;
      sprintf(node_name, "Sdc4u1_Node_%d", hit_wire);
      break;
    case 40:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4V)
	return;
      sprintf(node_name, "Sdc4v2_Node_%d", hit_wire);
      break;
    case 41:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4X)
	return;
      sprintf(node_name, "Sdc4x2_Node_%d", hit_wire);
      break;
    case 42:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4U)
	return;
      sprintf(node_name, "Sdc4u2_Node_%d", hit_wire);
      break;
      /*
    // BDC3,4 113-124
    case 113:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3x_Node_vtx_%d", hit_wire);
      break;
    case 114:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3xp_Node_vtx_%d", hit_wire);
      break;
    case 115:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3vp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3vp_Node_vtx_%d", hit_wire);
      break;
    case 116:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3v_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3v_Node_vtx_%d", hit_wire);
      break;
    case 117:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3up_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3up_Node_vtx_%d", hit_wire);
      break;
    case 118:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3u_Node_vtx_%d", hit_wire);
      break;
    case 119:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4u_Node_vtx_%d", hit_wire);
      break;
    case 120:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4up_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4up_Node_vtx_%d", hit_wire);
      break;
    case 121:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4v_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4v_Node_vtx_%d", hit_wire);
      break;
    case 122:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4vp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4vp_Node_vtx_%d", hit_wire);
      break;
    case 123:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4xp_Node_vtx_%d", hit_wire);
      break;
    case 124:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4x_Node_vtx_%d", hit_wire);
      break;
      */
    default:
      std::cerr << "EvDisp::DrawHitWire No such plane ID " << lnum << std::endl;
      return;
    }

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  if (range_check && tdc_check) 
    gevdisp_->GetNode(node_name)->SetLineColor(kBlue);
  else if (range_check && !tdc_check) 
    gevdisp_->GetNode(node_name)->SetLineColor(28);
  else 
    gevdisp_->GetNode(node_name)->SetLineColor(kBlack);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();


  if ((lnum>=1 && lnum<=10) || (lnum>=113 && lnum<=124 )) {
    gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
    if (range_check && tdc_check) 
      gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlue);
    else if (range_check && !tdc_check) 
      gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(28);
    else 
      gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlack);

    tp_vtx_->cd();
    gevdisp_vtx_->Draw();
    tc_vtx_->Update();
  }

}


void EvDisp::DrawSFTWire(int lnum, int hit_wire) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

  switch (lnum) {
    // SFT
  case 1:
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_X)
      return;
    sprintf(node_name, "SftX_Node_%d", hit_wire);
    sprintf(node_name_vtx, "SftX_Node_vtx_%d", hit_wire);
    break;
  case 2:
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_UV)
      return;
    sprintf(node_name, "SftV_Node_%d", hit_wire);
    sprintf(node_name_vtx, "SftV_Node_vtx_%d", hit_wire);
    break;
  case 3:
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_UV)
      return;
    sprintf(node_name, "SftU_Node_%d", hit_wire);
    sprintf(node_name_vtx, "SftU_Node_vtx_%d", hit_wire);
    break;
  default:
    std::cerr << "EvDisp::DrawSFTWire No such plane ID " << lnum << std::endl;
    return;
  }
  
  gevdisp_->GetNode(node_name)->SetVisibility(1);
  gevdisp_->GetNode(node_name)->SetLineColor(kBlue);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

  gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
  gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlue);

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tc_vtx_->Update();
}

void EvDisp::UpdateBcInCanvas() const
{
  tc_bcIn_->cd();
  tc_bcIn_->Update();

  //tc_bcIn_->Print("tc_bcIn.gif");
}

void EvDisp::UpdateBcOutCanvas() const
{
  tc_bcOut_->cd();
  tc_bcOut_->Update();

  //tc_bcOut_->Print("tc_bcOut.gif");
}

void EvDisp::UpdateBftCanvas() const
{
  static int evNum=0;

  tc_bft_->cd();

  Bft_tex_ = new TLatex();
  Bft_tex_->SetTextSize(0.07);
  Bft_tex_->SetTextColor(kWhite);
  char buf[100];
  sprintf(buf, "%d MHz", nClusterBFT_);
  tmp_Bft_tex_ = Bft_tex_->DrawLatex(40, -600, buf);

  tc_bft_->Update();

  sprintf(buf, "tc_bft%d.gif", evNum);
  tc_bft_->Print(buf);
  evNum++;

  tc_bft_rate_->cd();

  tc_bft_rate_->cd(1);

  if (Bft_rate_Index_ > MaxBFTRateNum)
    Bft_rate_Index_ = 0;

  Bft_event_Num_[Bft_rate_Index_] = Bft_rate_Index_;
  Bft_rate_[Bft_rate_Index_] = nClusterBFT_;
  Bft_rate_Index_++;

  Bft_rate_gr_ = new TGraph(Bft_rate_Index_, Bft_event_Num_, Bft_rate_);
  Bft_rate_gr_->SetName("Bft_rate_gr_");
  Bft_rate_gr_->SetTitle("BFT rate dependence");
  Bft_rate_gr_->SetMarkerColor(kRed);
  Bft_rate_gr_->SetMarkerStyle(20);
  Bft_rate_gr_->SetMarkerSize(0.5);
  Bft_rate_gr_->GetXaxis()->SetTitle("Event #");
  Bft_rate_gr_->GetYaxis()->SetTitle("Rate (MHz)");
  Bft_rate_gr_->Draw("ap");

  double mean_rate=5. * 1.3 ; // Hz, 1.3 is BH1/BH2 ratio
  Bft_rate_line_ = new TLine(0, mean_rate, Bft_rate_Index_, mean_rate);
  Bft_rate_line_->Draw("same");

  tc_bft_rate_->cd(2);
  hBft_rate_->Fill(nClusterBFT_);
  hBft_rate_->Draw();

  tc_bft_rate_->Update();

  //tc_bft_rate_->Print("tc_bft_rate.gif");
}

void EvDisp::DrawHitWireBcIn(int lnum, int hit_wire) const
{
  if (lnum<=0 || lnum>NumOfLayersBcIn) {
    std::cout << "EvDisp::DrawHitWireBcIn invalid layer : " << lnum 
	      << std::endl;
    return;
  }

  if (hit_wire<=0 || hit_wire>MaxWireBC1) {
    std::cout << "EvDisp::DrawHitWireBcIn layer : " << lnum 
	      << " invalid wire : " << hit_wire << std::endl;
    return;
  }

  tc_bcIn_->cd();

  BcIn_Arc_[lnum-1][hit_wire-1]->SetFillColor(kBlue);
  BcIn_Arc_[lnum-1][hit_wire-1]->SetLineColor(kBlue);
  BcIn_Arc_[lnum-1][hit_wire-1]->Draw("same");

  //tc_bcIn_->Update();
}

void EvDisp::DrawHitWireBcOut(int lnum, int hit_wire, int color) const
{
  if (lnum<=0 || lnum>NumOfLayersBcOut) {
    std::cout << "EvDisp::DrawHitWireBcOut invalid layer : " << lnum 
	      << std::endl;
    return;
  }

  if (lnum>=1 && lnum<=6) {
    if (hit_wire<=0 || hit_wire>MaxWireBC3) {
      std::cout << "EvDisp::DrawHitWireBcOut layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } else if (lnum>=7 && lnum<=12) {
    if (hit_wire<=0 || hit_wire>MaxWireBC4) {
      std::cout << "EvDisp::DrawHitWireBcOut layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } 

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  int colorNum;
  if (color < 7)
    colorNum = Color[color];
  else
    colorNum = kBlack;

  if (lnum>=1 && lnum<=6) {
    Bc3_Arc_[lnum-1][hit_wire-1]->SetFillColor(colorNum);
    Bc3_Arc_[lnum-1][hit_wire-1]->SetLineColor(colorNum);
    Bc3_Arc_[lnum-1][hit_wire-1]->Draw("same");
  } else if (lnum>=7 && lnum<=12) {
    Bc4_Arc_[lnum-7][hit_wire-1]->SetFillColor(colorNum);
    Bc4_Arc_[lnum-7][hit_wire-1]->SetLineColor(colorNum);
    Bc4_Arc_[lnum-7][hit_wire-1]->Draw("same");
  }

  //tc_bcIn_->Update();
}

void EvDisp::DrawHitWireSdcIn(int lnum, int hit_wire, int color) const
{
  if (lnum<=0 || lnum>NumOfLayersSdcIn+1) {
    std::cout << "EvDisp::DrawHitWireSdcIn invalid layer : " << lnum 
	      << std::endl;
    return;
  }

  if (lnum>=1 && lnum<=4) {
    if (hit_wire<=0 || hit_wire>MaxWireSDC1) {
      std::cout << "EvDisp::DrawHitWireSdcIn layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } else if (lnum>=5 && lnum<=10) {
    if (hit_wire<=0 || hit_wire>MaxWireSDC2) {
      std::cout << "EvDisp::DrawHitWireSdcIn layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } 

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  int colorNum;
  if (color < 7)
    colorNum = Color[color];
  else
    colorNum = kBlack;

  if (lnum>=1 && lnum<=4) {
    Sdc1_Arc_[lnum-1][hit_wire-1]->SetFillColor(colorNum);
    Sdc1_Arc_[lnum-1][hit_wire-1]->SetLineColor(colorNum);
    Sdc1_Arc_[lnum-1][hit_wire-1]->Draw("same");
  } else if (lnum>=5 && lnum<=10) {
    Sdc2_Arc_[lnum-5][hit_wire-1]->SetFillColor(colorNum);
    Sdc2_Arc_[lnum-5][hit_wire-1]->SetLineColor(colorNum);
    Sdc2_Arc_[lnum-5][hit_wire-1]->Draw("same");
  }

  //tc_bcIn_->Update();
}
#if 0
void EvDisp::DrawHitWireSft(int lnum, int hit_wire, double time) const
{
  if (lnum<=0 || lnum>NumOfLayersSFT) {
    std::cout << "EvDisp::DrawHitWireSft invalid layer : " << lnum 
	      << std::endl;
    return;
  }

  if (lnum==1) {
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_X * 2) {
      std::cout << "EvDisp::DrawHitWireSFT layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } else if (lnum>=2 && lnum<=3) {
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_UV) {
      std::cout << "EvDisp::DrawHitWireSFT layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } 

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  int colorNum = kRed;
  if (time<-5. || time>5.)
    colorNum = kBlack;

  if (lnum==1) {
    std::cout << hit_wire << ", " << (hit_wire+1)/2 << std::endl;
    if ((hit_wire+1)%2 == 0) {
      Sft_Arc_[0][(hit_wire+1)/2]->SetFillColor(colorNum);
      Sft_Arc_[0][(hit_wire+1)/2]->SetLineColor(colorNum);
      Sft_Arc_[0][(hit_wire+1)/2]->Draw("same");
    } else {
      Sft_Arc_[1][(hit_wire+1)/2]->SetFillColor(colorNum);
      Sft_Arc_[1][(hit_wire+1)/2]->SetLineColor(colorNum);
      Sft_Arc_[1][(hit_wire+1)/2]->Draw("same");
    }
  } else if (lnum>=2 && lnum<=3) {
    Sft_Arc_[lnum][hit_wire-1]->SetFillColor(colorNum);
    Sft_Arc_[lnum][hit_wire-1]->SetLineColor(colorNum);
    Sft_Arc_[lnum][hit_wire-1]->Draw("same");
  }

  //tc_bcIn_->Update();
}
#endif

void EvDisp::DrawHitWireSft(int lnum, int hit_wire, double time, int color) const
{
  if (lnum<=0 || lnum>NumOfLayersSFT) {
    std::cout << "EvDisp::DrawHitWireSft invalid layer : " << lnum 
	      << std::endl;
    return;
  }

  if (lnum==1) {
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_X * 2) {
      std::cout << "EvDisp::DrawHitWireSFT layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } else if (lnum>=2 && lnum<=3) {
    if (hit_wire<=0 || hit_wire>NumOfSegSFT_UV) {
      std::cout << "EvDisp::DrawHitWireSFT layer : " << lnum 
		<< " invalid wire : " << hit_wire << std::endl;
      return;
    }
  } 

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  int colorNum;

  if (color < 7)
    colorNum = Color[color];
  else 
    colorNum = kBlack;

  if (time<-5. || time>5.)
    colorNum = kBlack;

  if (lnum==1) {
    Sft_Arc_[0][hit_wire]->SetFillColor(colorNum);
    Sft_Arc_[0][hit_wire]->SetLineColor(colorNum);
    Sft_Arc_[0][hit_wire]->Draw("same");
  } else if (lnum>=2 && lnum<=3) {
    Sft_Arc_[lnum][hit_wire]->SetFillColor(colorNum);
    Sft_Arc_[lnum][hit_wire]->SetLineColor(colorNum);
    Sft_Arc_[lnum][hit_wire]->Draw("same");
  }

  //tc_bcIn_->Update();
}

void EvDisp::DrawHitBft(int lnum, int seg, bool timeflag) const
{
  if (lnum<0 || lnum>=NumOfPlaneBFT) {
    std::cout << "EvDisp::DrawHitBft invalid layer : " << lnum 
	      << std::endl;
    return;
  }

  if (seg<=0 || seg>NumOfSegBFT) {
    std::cout << "EvDisp::DrawHitBFT layer : " << lnum 
	      << " invalid segment : " << seg << std::endl;
    return;
  }

  tc_bcIn_->cd();

  if (timeflag) {
    Bft_Arc_[lnum][seg-1]->SetFillColor(kBlue);
    Bft_Arc_[lnum][seg-1]->SetLineColor(kBlue);
  } else {
    Bft_Arc_[lnum][seg-1]->SetFillColor(kGreen);
    Bft_Arc_[lnum][seg-1]->SetLineColor(kGreen);
  }
  Bft_Arc_[lnum][seg-1]->Draw("same");

  //tc_bcIn_->Update();
}

#if 0
void EvDisp::DrawBFTCluster(BFTCluster *cl) const
{

  if (nHitBFT_ >= MaxBFTHit) {
    std::cerr << "EvDisp::DrawBFTCluster nHitBFT_ is greater than MaxBFTHit"
	      << std::endl;
    return;
  }

  tc_bft_->cd();

  int color_num=nClusterBFT_%NumColor;

  int cs = cl->ClusterSize();
  for (int i=0; i<cs; i++) {
    BFTHit *hit = cl->GetHit(i);
    if (!hit) continue;

    int planeId   = hit->PlaneId();
    int segmentId = hit->SegmentId();
    double tdc       = hit->CTime();


    double localPos_x, localPos_t;
    if (planeId==0) {
      localPos_x = ((double)(segmentId+1)-80.5)*(1.);
      localPos_t = tdc;
    } else if (planeId==1){
      localPos_x = ((double)(segmentId+1)-80)*(1.);
      localPos_t = tdc+0.5;
    }

    Bft_Cl_Arc_[nHitBFT_] = new TArc(localPos_x, localPos_t, 0.5);

    Bft_Cl_Arc_[nHitBFT_]->SetLineColor(Color[color_num]);
    Bft_Cl_Arc_[nHitBFT_]->SetFillColor(Color[color_num]);
    Bft_Cl_Arc_[nHitBFT_]->Draw("same");

    nHitBFT_++;
  }
  nClusterBFT_++;
}
#endif

void EvDisp::DrawBFTCluster(double time, double pos) const
{

  if (nHitBFT_ >= MaxBFTHit) {
    std::cerr << "EvDisp::DrawBFTCluster nHitBFT_ is greater than MaxBFTHit"
	      << std::endl;
    return;
  }

  tc_bft_->cd();

  int color_num=nClusterBFT_%NumColor;

  if (nClusterBFT_ < MaxBFTHit) {
    Bft_Cl_Arc_[nClusterBFT_] = new TArc(pos, time, 0.5);
  
    Bft_Cl_Arc_[nClusterBFT_]->SetLineColor(Color[color_num]);
    Bft_Cl_Arc_[nClusterBFT_]->SetFillColor(Color[color_num]);
    Bft_Cl_Arc_[nClusterBFT_]->Draw("same");
  } else {
    std::cout << "EvDisp::DrawBFTCluster : BFT Cluster Num is over Max : " << nClusterBFT_ << std::endl;
  }

  nClusterBFT_++;
}

#if 0
void EvDisp::DrawTrackWire(int lnum, int hit_wire, int it) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

    switch (lnum) {
    // SDC1,2 1-11
    case 1:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1X)
	return;
      sprintf(node_name, "Sdc1x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1x_Node_vtx_%d", hit_wire);
      break;
    case 2:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1X)
	return;
      sprintf(node_name, "Sdc1xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1xp_Node_vtx_%d", hit_wire);
      break;
    case 3:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1Y)
	return;
      sprintf(node_name, "Sdc1y_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1y_Node_vtx_%d", hit_wire);
      break;
    case 4:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1Y)
	return;
      sprintf(node_name, "Sdc1yp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1yp_Node_vtx_%d", hit_wire);
      break;
    case 5:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1U)
	return;
      sprintf(node_name, "Sdc1u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1u_Node_vtx_%d", hit_wire);
      break;
    case 6:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2vp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2vp_Node_vtx_%d", hit_wire);
      break;
    case 7:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2v_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2v_Node_vtx_%d", hit_wire);
      break;
    case 8:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2up_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2up_Node_vtx_%d", hit_wire);
      break;
    case 9:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2u_Node_vtx_%d", hit_wire);
      break;
    case 10:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2xp_Node_vtx_%d", hit_wire);
      break;
    case 11:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2x_Node_vtx_%d", hit_wire);
      break;
    default:
      std::cerr << "EvDisp::DrawHitWire No such plane ID " << lnum << std::endl;
      return;
    }

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  gevdisp_->GetNode(node_name)->SetLineColor(5+it);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

  gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
  gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(5+it);

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tc_vtx_->Update();

}
#endif

void EvDisp::DrawHitBH1(int seg, int timeflag) const
{
  if (seg<=0 || seg>NumOfSegBH1) {
    std::cout << "EvDisp::DrawHitBH1 invalid BH1 segment : " << seg << std::endl;
    return;
  }

  tc_bcIn_->cd();
  for (int i=0; i<4; i++) {
    if (timeflag==0)
      Bh1_Line_[seg-1][i]->SetLineColor(kRed);
    else if (timeflag==1)
      Bh1_Line_[seg-1][i]->SetLineColor(kGreen);
    else if (timeflag==-1)
      Bh1_Line_[seg-1][i]->SetLineColor(kCyan);
    Bh1_Line_[seg-1][i]->Draw("same");
  }

  //tc_bcIn_->Update();


}

void EvDisp::DrawHitBH2(int seg, bool timeflag) const
{
  if (seg<=0 || seg>NumOfSegBH2) {
    std::cout << "EvDisp::DrawHitBH2 invalid BH2 segment : " << seg << std::endl;
    return;
  }

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  for (int i=0; i<4; i++) {
    if (timeflag)
      Bh2_Line_[seg-1][i]->SetLineColor(kRed);
    else
      Bh2_Line_[seg-1][i]->SetLineColor(kGreen);

    Bh2_Line_[seg-1][i]->Draw("same");
  }

  //tc_bcOut_->Update();


}

void EvDisp::DrawHitBH2(int seg, int Tu, int Td) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

  if (seg<0 || seg>=NumOfSegBH2)
    return;

  sprintf(node_name, "Bh2Seg_Node_%d", seg);
  sprintf(node_name_vtx, "Bh2Seg_Node_vtx_%d", seg);

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  if (Tu>0 && Td>0)
    gevdisp_->GetNode(node_name)->SetLineColor(kBlue);
  else
    gevdisp_->GetNode(node_name)->SetLineColor(kGreen);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

  gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
  if (Tu>0 && Td>0)
    gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlue);
  else
    gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kGreen);

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tc_vtx_->Update();

}

void EvDisp::DrawHitHodoscope(int lnum, int seg, int Tu, int Td) const
{
  char node_name[MaxChar];

  switch (lnum) {
    // Tof 51
  case 51:
    if (seg<0 || seg>=NumOfSegTOF)
      return;
    sprintf(node_name, "TofSeg_Node_%d", seg);
    break;
    // Tof 54
  case 54:
    if (seg<0 || seg>=NumOfSegLC)
      return;
    sprintf(node_name, "LcSeg_Node_%d", seg);
    break;
  default:
    std::cerr << "EvDisp::DrawHitHodoscope No such plane ID " << lnum << std::endl;
    return;
  }

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  if (Tu>0 && Td>0)
    gevdisp_->GetNode(node_name)->SetLineColor(kBlue);
  else
    gevdisp_->GetNode(node_name)->SetLineColor(kGreen);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

}


#if 0
void EvDisp::DrawBdcOutLocalTrack(ThreeVector globalPos0, 
				  ThreeVector globalPos1) const
{
  if (nTrackBdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBdcOutLocalTrack nTrackBdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackBdcOut_[nTrackBdcOut_] = new TPolyLine3D(2);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineColor(kRed);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineWidth(1);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackBdcOut_[nTrackBdcOut_]->Draw();
  tc_->Update();
  nTrackBdcOut_++;

}

void EvDisp::DrawBdcOutLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackBdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBdcOutLocalTrack nTrackBdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  int IdBdc3x = DCGeomMan::GetInstance().GetDetectorId("BDC3-x-1");
  double zBdc3x = DCGeomMan::GetInstance().GetLocalZ(IdBdc3x);
  double x0=tp->GetX(zBdc3x), y0=tp->GetY(zBdc3x);

  int IdK18Target = DCGeomMan::GetInstance().GetDetectorId("K18Target");
  double zK18Target = DCGeomMan::GetInstance().GetLocalZ(IdK18Target);
  double x1=tp->GetX(zK18Target), y1=tp->GetY(zK18Target);

  ThreeVector globalPosTgt = DCGeomMan::GetInstance().GetGlobalPosition(0);
  
  ThreeVector pos0 = ThreeVector(zBdc3x-zK18Target, x0, y0);
  ThreeVector globalPos0 = pos0.rotateZ(130.0*Deg2Rad)+globalPosTgt;
  
  ThreeVector pos1 = ThreeVector(0.0, x1, y1);
  ThreeVector globalPos1 = pos1.rotateZ(130.0*Deg2Rad)+globalPosTgt;
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackBdcOut_[nTrackBdcOut_] = new TPolyLine3D(2);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineColor(kRed);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineWidth(1);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackBdcOut_[nTrackBdcOut_]->Draw();
  tc_->Update();

  tp_vtx_->cd();
  LocalTrackBdcOut_[nTrackBdcOut_]->Draw();
  tc_vtx_->Update();

  nTrackBdcOut_++;

}
#endif


void EvDisp::DrawSdcInLocalTrack(ThreeVector globalPos0, 
				 ThreeVector globalPos1) const
{
  if (nTrackSdcIn_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcIn_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackSdcIn_[nTrackSdcIn_] = new TPolyLine3D(2);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineColor(kRed);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineWidth(1);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcIn_[nTrackSdcIn_]->Draw();
  tc_->Update();
  nTrackSdcIn_++;

}

void EvDisp::DrawBcInLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackBcIn_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBcInLocalTrack nTrackBcIn_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  double x0=tp->GetX0();
  double x1=tp->GetX(-600.);
  double chisqr=tp->GetChiSquare();
  //std::cout << "BcIn : chi2 " << chisqr << std::endl;

  tc_bcIn_->cd();
  LocalTrackBcIn2_[nTrackBcIn_] = new TLine(x0, 0, x1, -600);
  if (chisqr < 10.)
    LocalTrackBcIn2_[nTrackBcIn_]->SetLineColor(kRed);
  else
    LocalTrackBcIn2_[nTrackBcIn_]->SetLineColor(kGreen);

  LocalTrackBcIn2_[nTrackBcIn_]->Draw("same");
  tc_bcIn_->Update();

  nTrackBcIn_++;

}

void EvDisp::DrawBcInLocalTrack(double xbft, double xbh1, bool flagBh1) const
{
  if (nTrackBcIn_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBcInLocalTrack nTrackBcIn_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  double zbft=DCGeomMan::GetInstance().GetLocalZ( IdBFT );
  double zbh1=DCGeomMan::GetInstance().GetLocalZ( IdBH1 );

  tc_bcIn_->cd();
  LocalTrackBcIn2_[nTrackBcIn_] = new TLine(xbft, zbft, xbh1, zbh1);
  if (flagBh1)
    LocalTrackBcIn2_[nTrackBcIn_]->SetLineColor(kRed);
  else
    LocalTrackBcIn2_[nTrackBcIn_]->SetLineColor(kGreen);

  LocalTrackBcIn2_[nTrackBcIn_]->Draw("same");
  tc_bcIn_->Update();

  nTrackBcIn_++;

}

void EvDisp::DrawBcOutLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackBcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBcOutLocalTrack nTrackBcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
  double zK18Tgt = DCGeomMan::GetInstance().GetLocalZ(130);
  double xoffset=3.5;

  double x0=tp->GetX0()+xoffset;
  double x1=tp->GetX(zK18Tgt)+xoffset;

  double yoffset=3.0;
  double y0=tp->GetY0()+yoffset;
  double y1=tp->GetY(zK18Tgt)+yoffset;

  double chisqr=tp->GetChiSquare();
  //std::cout << "BcOut : chi2 " << chisqr << std::endl;

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  int nh=tp->GetNHit();
  int index_gr=0;
  double ypos_gr[12], zpos_gr[12];
  for( int ih=0; ih<nh; ++ih ){
    DCLTrackHit *hit=tp->GetHit(ih);
    int layerId=hit->GetLayer()-112; 
    int wire=(int)hit->GetWire();
    DrawHitWireBcOut(layerId, wire, nTrackBcOut_+1);

    double localPos = hit->GetLocalHitPos();
    double xcal     = hit->GetXcal();
    double tilt     = hit->GetTiltAngle();
    if (fabs(tilt)>1.) {
      double x0 = localPos*cos(tilt*Deg2Rad);
      double y0 = localPos*sin(tilt*Deg2Rad);
      double y  = -cos(tilt*Deg2Rad)/sin(tilt*Deg2Rad)*(xcal-x0)+y0;
      int layer=hit->GetLayer();
      double z = DCGeomMan::GetInstance().GetLocalZ(layer) - zK18Tgt;
      ypos_gr[index_gr] = y;
      zpos_gr[index_gr] = z;
      index_gr++;
    }
  }    

  int colorNum=Color[nTrackBcOut_+1];

  // ZX plane
  LocalTrackBcOut2_[nTrackBcOut_] = new TLine(x0, -zK18Tgt, x1, 0);
  /*
  if (chisqr < 10.)
    LocalTrackBcOut2_[nTrackBcOut_]->SetLineColor(kRed);
  else
    LocalTrackBcOut2_[nTrackBcOut_]->SetLineColor(kGreen);
  */
  LocalTrackBcOut2_[nTrackBcOut_]->SetLineColor(colorNum);
  LocalTrackBcOut2_[nTrackBcOut_]->Draw("same");

  tp_bcOut2_->cd();
  // ZY plane
  LocalTrackBcOut2_Y_[nTrackBcOut_] = new TLine(y0, -zK18Tgt, y1, 0);
  LocalTrackBcOut2_Y_[nTrackBcOut_]->SetLineColor(colorNum);

  LocalTrackBcOut2_Y_[nTrackBcOut_]->Draw("same");

  // Y hit pos
  BcOut_YPos_gr_[nTrackBcOut_] = new TGraph(index_gr, ypos_gr, zpos_gr);
  BcOut_YPos_gr_[nTrackBcOut_]->SetMarkerColor(colorNum);
  BcOut_YPos_gr_[nTrackBcOut_]->SetMarkerStyle(20);
  BcOut_YPos_gr_[nTrackBcOut_]->SetMarkerSize(0.5);
  BcOut_YPos_gr_[nTrackBcOut_]->Draw("p");


  tc_bcOut_->Update();

  nTrackBcOut_++;

}


void EvDisp::DrawSdcInLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackSdcIn_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcIn_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  double x0=tp->GetX0(), y0=tp->GetY0();
  
  int IdSdc2x2 = 10;
  int IdTgt = 0;
  double zSdc2x2 = DCGeomMan::GetInstance().GetLocalZ(IdSdc2x2);
  double zTgt = DCGeomMan::GetInstance().GetLocalZ(IdTgt);
  double x1=tp->GetX(zSdc2x2), y1=tp->GetY(zSdc2x2);
  
  ThreeVector globalPosTgt = DCGeomMan::GetInstance().GetGlobalPosition(0);
  
  ThreeVector pos0 = ThreeVector(0.0, x0, y0);
  ThreeVector globalPos0 = pos0.rotateZ(145.0*Deg2Rad)+globalPosTgt;
  
  ThreeVector pos1 = ThreeVector(zSdc2x2-zTgt, x1, y1);
  ThreeVector globalPos1 = pos1.rotateZ(145.0*Deg2Rad)+globalPosTgt;

  double chisqr=tp->GetChiSquare();
  //std::cout << "SdcIn : chi2 " << chisqr << std::endl;

#if 0  
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackSdcIn_[nTrackSdcIn_] = new TPolyLine3D(2);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineColor(kRed);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineWidth(1);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcIn_[nTrackSdcIn_]->Draw();
  tc_->Update();

  tp_vtx_->cd();
  LocalTrackSdcIn_[nTrackSdcIn_]->Draw();
  tc_vtx_->Update();

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  int nh=tp->GetNHit();

  int index_gr=0;
  double ypos_gr[12], zpos_gr[12];

  for( int ih=0; ih<nh; ++ih ){
    DCLTrackHit *hit=tp->GetHit(ih);
    int layerId=hit->GetLayer(); 
    //int wire=(int)hit->GetWire();
    int wire;

    if (layerId==2) {
      double pos = hit->GetLocalHitPos();
      wire = (int)pos+128;
      if (wire<0)
	wire = 0;
      else if (wire>=256)
	wire = 255;

      DrawHitWireSft(layerId-1, wire, 0., nTrackSdcIn_+1);

    } else if (layerId==3) {
      double pos = hit->GetLocalHitPos();
      wire = 160+(int)(pos/0.75);
      if (wire<0)
	wire = 0;
      else if (wire>320)
	wire = 319;
      DrawHitWireSft(layerId-1, wire, 0., nTrackSdcIn_+1);

    } else if (layerId==4) { 
      double pos = hit->GetLocalHitPos();
      wire = 160-(int)(pos/0.75);
      if (wire<0)
	wire = 0;
      else if (wire>320)
	wire = 319;
      DrawHitWireSft(layerId-1, wire, 0., nTrackSdcIn_+1);
    } else {
      wire=(int)hit->GetWire();
      DrawHitWireSdcIn(layerId, wire, nTrackSdcIn_+1);
    }


    double localPos = hit->GetLocalHitPos();
    double xcal     = hit->GetXcal();
    double tilt     = hit->GetTiltAngle();
    if (fabs(tilt)>1.) {
      double x0 = localPos*cos(tilt*Deg2Rad);
      double y0 = localPos*sin(tilt*Deg2Rad);
      double y  = -cos(tilt*Deg2Rad)/sin(tilt*Deg2Rad)*(xcal-x0)+y0;
      int layer=hit->GetLayer();
      double z = DCGeomMan::GetInstance().GetLocalZ(layer);
      ypos_gr[index_gr] = y;
      zpos_gr[index_gr] = z;
      index_gr++;
    }


  }    

  int colorNum=Color[nTrackSdcIn_+1];

  LocalTrackSdcIn2_[nTrackSdcIn_] = new TLine(x0, 0, x1, zSdc2x2);
  /*
  if (chisqr < 10.)
    LocalTrackSdcIn2_[nTrackSdcIn_]->SetLineColor(kRed);
  else
    LocalTrackSdcIn2_[nTrackSdcIn_]->SetLineColor(kGreen);
  */

  LocalTrackSdcIn2_[nTrackSdcIn_]->SetLineColor(colorNum);
  LocalTrackSdcIn2_[nTrackSdcIn_]->Draw("same");


  tp_bcOut2_->cd();
  LocalTrackSdcIn2_Y_[nTrackSdcIn_] = new TLine(y0, 0, y1, zSdc2x2);
  LocalTrackSdcIn2_Y_[nTrackSdcIn_]->SetLineColor(colorNum);
  LocalTrackSdcIn2_Y_[nTrackSdcIn_]->Draw("same");


  // Y hit pos
  SdcIn_YPos_gr_[nTrackSdcIn_] = new TGraph(index_gr, ypos_gr, zpos_gr);
  SdcIn_YPos_gr_[nTrackSdcIn_]->SetMarkerColor(colorNum);
  SdcIn_YPos_gr_[nTrackSdcIn_]->SetMarkerStyle(20);
  SdcIn_YPos_gr_[nTrackSdcIn_]->SetMarkerSize(0.5);
  SdcIn_YPos_gr_[nTrackSdcIn_]->Draw("p");


  tc_bcOut_->Update();

  nTrackSdcIn_++;

}


#if 0
void EvDisp::DrawSdcOutLocalTrack(ThreeVector globalPos0, 
				 ThreeVector globalPos1) const
{
  if (nTrackSdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
#if 0 
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif

  LocalTrackSdcOut_[nTrackSdcOut_] = new TPolyLine3D(2);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineColor(kRed);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineWidth(1);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcOut_[nTrackSdcOut_]->Draw();
  tc_->Update();
  nTrackSdcOut_++;
}
#endif

void EvDisp::DrawSdcOutLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackSdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
  /*
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 
  double x0=tp->GetX(zTof), y0=tp->GetY(zTof);
  */
  int IdLc = DCGeomMan::GetInstance().GetLcId();
  double zLc = DCGeomMan::GetInstance().GetLocalZ( IdLc ); 
  double x0=tp->GetX(zLc), y0=tp->GetY(zLc);
  
  int IdSdc3x1 = DCGeomMan::GetInstance().GetDetectorId("SDC3-x-1");
  double zSdc3x1 = DCGeomMan::GetInstance().GetLocalZ(IdSdc3x1);
  double x1=tp->GetX(zSdc3x1), y1=tp->GetY(zSdc3x1);


  ThreeVector pos0 = ThreeVector(x0, y0, 0.0);
  //ThreeVector globalPos0 = DCGeomMan::GetInstance().Local2GlobalPos(IdTof, pos0);
  ThreeVector globalPos0 = DCGeomMan::GetInstance().Local2GlobalPos(IdLc, pos0);

  ThreeVector pos1 = ThreeVector(x1, y1, 0.0);
  ThreeVector globalPos1 = DCGeomMan::GetInstance().Local2GlobalPos(IdSdc3x1, pos1);
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif

  LocalTrackSdcOut_[nTrackSdcOut_] = new TPolyLine3D(2);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineColor(kRed);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineWidth(1);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcOut_[nTrackSdcOut_]->Draw();
  tc_->Update();
  nTrackSdcOut_++;
}

void EvDisp::DrawSksTrack(int nStep, ThreeVector *StepPoint) const
{
  if (SksStepMark_)
    SksStepMark_->Delete();

  SksStepMark_ = new TPolyMarker3D(nStep);
  for (int i=0; i<nStep; i++) {
    SksStepMark_->SetPoint(i, StepPoint[i].x(), StepPoint[i].y(),  StepPoint[i].z());
  }
  SksStepMark_->SetMarkerSize(1);
  SksStepMark_->SetMarkerColor(kBlue);
  SksStepMark_->SetMarkerStyle(6);

  tp_->cd();
  SksStepMark_->Draw();

  tc_->Update();
}

#if 0
void EvDisp::DrawVertex(ThreeVector vtxPoint, double cost) const
{
  if (nVertex_>=MaxTrack) {
    std::cerr << "EvDisp::DrawVertex nVertex_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  VtxPoint_[nVertex_] = new TPolyMarker3D(1);
  VtxPoint_[nVertex_]->SetPoint(0, vtxPoint.x(), vtxPoint.y(), vtxPoint.z());
  VtxPoint_[nVertex_]->SetMarkerSize(1);
  if (cost<0.995)
    VtxPoint_[nVertex_]->SetMarkerColor(6); //Pink
  else
    VtxPoint_[nVertex_]->SetMarkerColor(3); //Green
  VtxPoint_[nVertex_]->SetMarkerStyle(29); // Star   

  tp_->cd();
  VtxPoint_[nVertex_]->Draw();
  if (nVertex_==0) {
    tp_->GetView()->ZoomIn();
    tp_->GetView()->ZoomIn();
    tp_->GetView()->ZoomIn();
  }
  tc_->Update();

  tp_vtx_->cd();
  VtxPoint_[nVertex_]->Draw();
  if (nVertex_==0) {
    tp_vtx_->GetView()->ZoomIn();
    tp_vtx_->GetView()->ZoomIn();
    tp_vtx_->GetView()->ZoomIn();
    tp_vtx_->GetView()->ZoomIn();
  }
  tc_vtx_->Update();

  nVertex_++;
}
#endif


void EvDisp::WriteEventNumber(int spill, int evNum) const
{
  tc_->cd();
  tp2_->cd();

  EvNum_tex_ = new TLatex();
  EvNum_tex_->SetTextSize(0.15);
  EvNum_tex_->SetTextColor(kRed);
  char buf[100];
  sprintf(buf, "Run# %d, Spill# %d, Event# %d", RunNum_, spill, evNum);
  tmp_EvNum_tex_ = EvNum_tex_->DrawLatex(0.02, 0.85, buf);

}

void EvDisp::WriteSksComment(char *buf) const
{
  static int num_comment=0;
  tc_->cd();
  tp2_->cd();

  if (! Sks_tex_) {
    Sks_tex_ = new TLatex();
    Sks_tex_->SetTextSize(0.1);
    Sks_tex_->SetTextColor(kBlack);
    num_comment=0;
  }

  if (num_comment<MaxComment)
    tmp_Sks_tex_[num_comment] = Sks_tex_->DrawLatex(0.02, 0.7-0.1*(double)num_comment, buf);
  num_comment++;
}

void EvDisp::WriteVertexComment(char *buf) const
{
  static int num_comment=0;
  tc_->cd();
  tp2_->cd();

  if (! vtx_tex_) {
    vtx_tex_ = new TLatex();
    vtx_tex_->SetTextSize(0.1);
    vtx_tex_->SetTextColor(kBlack);
    num_comment=0;
  }

  if (num_comment<MaxComment)
    tmp_vtx_tex_[num_comment] = vtx_tex_->DrawLatex(0.5, 0.7-0.1*(double)num_comment, buf);
  num_comment++;
}


void EvDisp::WriteK18Comment(char *buf) const
{
  static int num_comment=0;
  tc_bcIn_->cd();


  if (! K18_tex_) {
    K18_tex_ = new TLatex();
    K18_tex_->SetTextSize(0.04);
    K18_tex_->SetTextColor(kBlack);
    num_comment=0;
  }

  if (num_comment<MaxComment)
    tmp_K18_tex_[num_comment] = K18_tex_->DrawLatex(-120, -600-30*(double)num_comment, buf);
  num_comment++;
}


void EvDisp::DrawVertex2D(double x, double y, double z, double cdist) const
{

  tc_bcOut_->cd();
  tp_bcOut1_->cd();

  Vertex_Arc_[nVertex_] = new TArc(x, z, 5.);
  Vertex_Arc_[nVertex_]->SetFillColor(Color[nVertex_+1]);

  if (cdist<10)
    Vertex_Arc_[nVertex_]->SetLineColor(Color[nVertex_+1]);
  else
    Vertex_Arc_[nVertex_]->SetLineColor(kBlack);

  Vertex_Arc_[nVertex_]->Draw("same");

  tp_bcOut2_->cd();

  Vertex_Arc_Y_[nVertex_] = new TArc(y, z, 5.);
  Vertex_Arc_Y_[nVertex_]->SetFillColor(Color[nVertex_+1]);

  if (cdist<10)
    Vertex_Arc_Y_[nVertex_]->SetLineColor(Color[nVertex_+1]);
  else
    Vertex_Arc_Y_[nVertex_]->SetLineColor(kBlack);

  Vertex_Arc_Y_[nVertex_]->Draw("same");


  nVertex_++;

}


void EvDisp::EndOfEvent(void) const
{

  //tc_->Print("tc.gif");

  //nTrackBdcOut_=0;
  nTrackBcIn_=0;
  nTrackBcOut_=0;
  nTrackSdcIn_=0;
  nTrackSdcOut_=0;
  nClusterBFT_=0;
  nHitBFT_=0;
  nVertex_=0;

  if (InitStepMark_) {
    delete InitStepMark_;
    InitStepMark_ = NULL;
  }
  /*
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBdcOut_[i]) {
      delete LocalTrackBdcOut_[i];
      LocalTrackBdcOut_[i] = NULL;
    }
  */
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn_[i]) {
      delete LocalTrackSdcIn_[i];
      LocalTrackSdcIn_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBcIn2_[i]) {
      delete LocalTrackBcIn2_[i];
      LocalTrackBcIn2_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBcOut2_[i]) {
      delete LocalTrackBcOut2_[i];
      LocalTrackBcOut2_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBcOut2_Y_[i]) {
      delete LocalTrackBcOut2_Y_[i];
      LocalTrackBcOut2_Y_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (BcOut_YPos_gr_[i]) {
      delete BcOut_YPos_gr_[i];
      BcOut_YPos_gr_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn2_[i]) {
      delete LocalTrackSdcIn2_[i];
      LocalTrackSdcIn2_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn2_Y_[i]) {
      delete LocalTrackSdcIn2_Y_[i];
      LocalTrackSdcIn2_Y_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (SdcIn_YPos_gr_[i]) {
      delete SdcIn_YPos_gr_[i];
      SdcIn_YPos_gr_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcOut_[i]) {
      delete LocalTrackSdcOut_[i];
      LocalTrackSdcOut_[i] = NULL;
    }
  /*
  for (int i=0; i<MaxTrack; i++)
    if (VtxPoint_[i]) {
      delete VtxPoint_[i];
      VtxPoint_[i] = NULL;
    }
  */
  if (SksStepMark_) {
    delete SksStepMark_;
    SksStepMark_ = NULL;
  }

  for (int i=0; i<MaxBFTHit; i++)
    if (Bft_Cl_Arc_[i]) {
      delete Bft_Cl_Arc_[i];
      Bft_Cl_Arc_[i] = NULL;
    }

  if (Bft_tex_) {
    delete Bft_tex_;
    Bft_tex_ = NULL;
  }

  if (tmp_Bft_tex_) {
    delete tmp_Bft_tex_;
    tmp_Bft_tex_ = NULL;
  }

  if (EvNum_tex_) {
    delete EvNum_tex_;
    EvNum_tex_ = NULL;
  }

  if (tmp_EvNum_tex_) {
    delete tmp_EvNum_tex_;
    tmp_EvNum_tex_ = NULL;
  }

  if (Sks_tex_) {
    delete Sks_tex_;
    Sks_tex_ = NULL;
  }

  for (int i=0; i<MaxComment; i++)
    if (tmp_Sks_tex_[i]) {
      delete tmp_Sks_tex_[i];
      tmp_Sks_tex_[i] = NULL;
    }

  if (vtx_tex_) {
    delete vtx_tex_;
    vtx_tex_ = NULL;
  }

  for (int i=0; i<MaxComment; i++)
    if (tmp_vtx_tex_[i]) {
      delete tmp_vtx_tex_[i];
      tmp_vtx_tex_[i] = NULL;
    }

  if (K18_tex_) {
    delete K18_tex_;
    K18_tex_ = NULL;
  }

  for (int i=0; i<MaxComment; i++)
    if (tmp_K18_tex_[i]) {
      delete tmp_K18_tex_[i];
      tmp_K18_tex_[i] = NULL;
    }

  for (int i=0; i<MaxVertex; i++)
    if (Vertex_Arc_[i]){
      delete Vertex_Arc_[i];
      Vertex_Arc_[i] = NULL;
    }

  for (int i=0; i<MaxVertex; i++)
    if (Vertex_Arc_Y_[i]){
      delete Vertex_Arc_Y_[i];
      Vertex_Arc_Y_[i] = NULL;
    }

  if (Bft_rate_gr_) {
    delete Bft_rate_gr_;
    Bft_rate_gr_ = NULL;
  }

  if (Bft_rate_line_) {
    delete Bft_rate_line_;
    Bft_rate_line_ = NULL;
  }
    
  ResetVisibility();

}

void EvDisp::ResetVisibility(void) const
{
  /*
  for (int i=0; i<MaxWireBDC; i++) {
    Bdc3x_Node_[i]->SetVisibility(0);
    Bdc3xp_Node_[i]->SetVisibility(0);
    Bdc3u_Node_[i]->SetVisibility(0);
    Bdc3up_Node_[i]->SetVisibility(0);
    Bdc3v_Node_[i]->SetVisibility(0);
    Bdc3vp_Node_[i]->SetVisibility(0);
    Bdc4x_Node_[i]->SetVisibility(0);
    Bdc4xp_Node_[i]->SetVisibility(0);
    Bdc4u_Node_[i]->SetVisibility(0);
    Bdc4up_Node_[i]->SetVisibility(0);
    Bdc4v_Node_[i]->SetVisibility(0);
    Bdc4vp_Node_[i]->SetVisibility(0);
  }
  */

  /*
  for (int i=0; i<MaxWireSDC1; i++) {
    Sdc1u1_Node_[i]->SetVisibility(0);
    Sdc1u2_Node_[i]->SetVisibility(0);
    Sdc1v1_Node_[i]->SetVisibility(0);
    Sdc1v2_Node_[i]->SetVisibility(0);

    Sdc1u1_Node_vtx_[i]->SetVisibility(0);
    Sdc1u2_Node_vtx_[i]->SetVisibility(0);
    Sdc1v1_Node_vtx_[i]->SetVisibility(0);
    Sdc1v2_Node_vtx_[i]->SetVisibility(0);
  }
  */

  for (int i=0; i<NumOfSegSFT_X; i++) {
    SftX_Node_[i]->SetVisibility(0);
    SftX_Node_vtx_[i]->SetVisibility(0);
  }
  for (int i=0; i<NumOfSegSFT_UV; i++) {
    SftU_Node_[i]->SetVisibility(0);
    SftU_Node_vtx_[i]->SetVisibility(0);

    SftV_Node_[i]->SetVisibility(0);
    SftV_Node_vtx_[i]->SetVisibility(0);
  }

  for (int i=0; i<MaxWireSDC2; i++) {
    Sdc2u1_Node_[i]->SetVisibility(0);
    Sdc2u2_Node_[i]->SetVisibility(0);
    Sdc2v1_Node_[i]->SetVisibility(0);
    Sdc2v2_Node_[i]->SetVisibility(0);
    Sdc2x1_Node_[i]->SetVisibility(0);
    Sdc2x2_Node_[i]->SetVisibility(0);

    Sdc2u1_Node_vtx_[i]->SetVisibility(0);
    Sdc2u2_Node_vtx_[i]->SetVisibility(0);
    Sdc2v1_Node_vtx_[i]->SetVisibility(0);
    Sdc2v2_Node_vtx_[i]->SetVisibility(0);
    Sdc2x1_Node_vtx_[i]->SetVisibility(0);
    Sdc2x2_Node_vtx_[i]->SetVisibility(0);
  }

  for (int i=0; i<MaxWireSDC3V; i++) {
    Sdc3v1_Node_[i]->SetVisibility(0);
    Sdc3v2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC3X; i++) {
    Sdc3x1_Node_[i]->SetVisibility(0);
    Sdc3x2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC3U; i++) {
    Sdc3u1_Node_[i]->SetVisibility(0);
    Sdc3u2_Node_[i]->SetVisibility(0);
  }

  for (int i=0; i<MaxWireSDC4V; i++) {
    Sdc4v1_Node_[i]->SetVisibility(0);
    Sdc4v2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC4X; i++) {
    Sdc4x1_Node_[i]->SetVisibility(0);
    Sdc4x2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC4U; i++) {
    Sdc4u1_Node_[i]->SetVisibility(0);
    Sdc4u2_Node_[i]->SetVisibility(0);
  }

  /*
  for (int i=0; i<MaxWireBDC; i++) {
    Bdc3x_Node_vtx_[i]->SetVisibility(0);
    Bdc3xp_Node_vtx_[i]->SetVisibility(0);
    Bdc3u_Node_vtx_[i]->SetVisibility(0);
    Bdc3up_Node_vtx_[i]->SetVisibility(0);
    Bdc3v_Node_vtx_[i]->SetVisibility(0);
    Bdc3vp_Node_vtx_[i]->SetVisibility(0);
    Bdc4x_Node_vtx_[i]->SetVisibility(0);
    Bdc4xp_Node_vtx_[i]->SetVisibility(0);
    Bdc4u_Node_vtx_[i]->SetVisibility(0);
    Bdc4up_Node_vtx_[i]->SetVisibility(0);
    Bdc4v_Node_vtx_[i]->SetVisibility(0);
    Bdc4vp_Node_vtx_[i]->SetVisibility(0);
  }
  */
  /*
  for (int i=0; i<NumOfSegBH2; i++) {
    Bh2Seg_Node_[i]->SetLineColor(kBlack);
    Bh2Seg_Node_vtx_[i]->SetLineColor(kBlack);
  }
  */
  for (int i=0; i<NumOfSegTOF; i++) 
    TofSeg_Node_[i]->SetLineColor(kBlack);

  for (int i=0; i<NumOfSegLC; i++) 
    LcSeg_Node_[i]->SetLineColor(kBlack);


  for (int i=0; i<NumOfLayersBcIn; i++)
    for (int j=0; j<MaxWireBC1; j++){
      BcIn_Arc_[i][j]->SetLineColor(kWhite);
      BcIn_Arc_[i][j]->SetFillColor(kWhite);
    }

  for (int i=0; i<NumOfSegBH1; i++)
    for (int j=0; j<4; j++)
      Bh1_Line_[i][j]->SetLineColor(kBlack);

  for (int i=0; i<NumOfPlaneBFT; i++)
    for (int j=0; j<NumOfSegBFT; j++){
      Bft_Arc_[i][j]->SetLineColor(kWhite);
      Bft_Arc_[i][j]->SetFillColor(kWhite);
    }

  for (int i=0; i<6; i++)
    for (int j=0; j<MaxWireBC3; j++){
      Bc3_Arc_[i][j]->SetLineColor(kWhite);
      Bc3_Arc_[i][j]->SetFillColor(kWhite);
    }

  for (int i=0; i<6; i++)
    for (int j=0; j<MaxWireBC4; j++){
      Bc4_Arc_[i][j]->SetLineColor(kWhite);
      Bc4_Arc_[i][j]->SetFillColor(kWhite);
    }

  for (int i=0; i<NumOfSegBH2; i++)
    for (int j=0; j<4; j++)
      Bh2_Line_[i][j]->SetLineColor(kBlack);

  /*
  for (int i=0; i<4; i++)
    for (int j=0; j<MaxWireSDC1; j++){
      Sdc1_Arc_[i][j]->SetLineColor(kWhite);
      Sdc1_Arc_[i][j]->SetFillColor(kWhite);
    }
  */

  for (int i=0; i<2; i++)
    for (int j=0; j<NumOfSegSFT_X; j++){
      Sft_Arc_[i][j]->SetLineColor(kWhite);
      Sft_Arc_[i][j]->SetFillColor(kWhite);
    }

  for (int i=2; i<4; i++)
    for (int j=0; j<NumOfSegSFT_UV; j++){
      Sft_Arc_[i][j]->SetLineColor(kWhite);
      Sft_Arc_[i][j]->SetFillColor(kWhite);
    }

  for (int i=0; i<6; i++)
    for (int j=0; j<MaxWireSDC2; j++){
      Sdc2_Arc_[i][j]->SetLineColor(kWhite);
      Sdc2_Arc_[i][j]->SetFillColor(kWhite);
    }

}


void EvDisp::calcRotMatrix(double TA, double RA1, double RA2, double *rotMat)
{
  double ct1=cos(RA1*Deg2Rad), st1=sin(RA1*Deg2Rad);
  double ct2=cos(RA2*Deg2Rad), st2=sin(RA2*Deg2Rad);
  double ct0=cos(TA*Deg2Rad), st0=sin(TA*Deg2Rad);

  double rotMat1[3][3], rotMat2[3][3];

  /* rotation matrix which is same as DCGeomRecord.cc*/
  rotMat1[0][0] =  ct0*ct2-st0*ct1*st2;
  rotMat1[0][1] = -st0*ct2-ct0*ct1*st2;
  rotMat1[0][2] =  st1*st2;

  rotMat1[1][0] =  ct0*st2+st0*ct1*ct2;
  rotMat1[1][1] = -st0*st2+ct0*ct1*ct2;
  rotMat1[1][2] = -st1*ct2;

  rotMat1[2][0] =  st0*st1;
  rotMat1[2][1] =  ct0*st1;
  rotMat1[2][2] =  ct1;

  /* rotation matrix which rotate -90 deg at x axis*/
  rotMat2[0][0] =  1.0;
  rotMat2[0][1] =  0.0;
  rotMat2[0][2] =  0.0;

  rotMat2[1][0] =  0.0;
  rotMat2[1][1] =  0.0;
  rotMat2[1][2] =  1.0;

  rotMat2[2][0] =  0.0;
  rotMat2[2][1] = -1.0;
  rotMat2[2][2] =  0.0;

  for (int i=0; i<9; i++)
    rotMat[i]=0.0;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	//rotMat[3*i+j] += rotMat1[i][k]*rotMat2[k][j];
	rotMat[i+3*j] += rotMat1[i][k]*rotMat2[k][j];
      }
    }
  }

}


void EvDisp::get_command(void) const
{
  tc_->Update();

  tc_bcOut_->Update();

  char ch;
  char data[100];
  static int stat=0;
  static int Nevent=0;
  static int ev=0;
  
  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    if (nVertex_>0)
      sleep(3);
    else
      sleep(1);

    ev++;
    return;
  } 
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    printf("q|n|p>");

    /* get command */
    scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');

    switch (ch) {
    case 'q': exit(0);
    case 'n':
      stat = 1;
      do {
	printf("event#>");
	scanf("%s",data);
      } while ((Nevent=atoi(data))<=0);
      std::cout << "Continue " << Nevent << "event" << std::endl;
      break;
    case 'p':
      theApp->Run(kTRUE);
      break;
    }
  }
}
