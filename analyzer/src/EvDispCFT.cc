/*
  EvDispCFT.cc

  2012/1/24
*/

//Unit is mm.

#include "EvDispCFT.hh"
#include "DCGeomMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <string>
#include <algorithm>

#include "TemplateLib.hh"

#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
//#include "BFTCluster.hh"
//#include "BFTHit.hh"

const int MaxChar = 200;

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

EvDispCFT *EvDispCFT::evDisp_ = 0;
TApplication *EvDispCFT::theApp =0;

const int NumColor = 10;
const int Color[NumColor] = {kBlue, kRed, kCyan,  kGreen, kPink, kOrange, kYellow, kMagenta, kViolet, kSpring};

EvDispCFT::EvDispCFT(void)
{
  gr_VertexXY_ = 0;
  gr_VertexZX_ = 0;
  gr_VertexZY_ = 0;
}

EvDispCFT::~EvDispCFT(void)
{

}

EvDispCFT & EvDispCFT::GetInstance( void )
{
  if( !evDisp_ ){
    evDisp_ = new EvDispCFT();
  }
  if( !theApp ){
    theApp=new TApplication( "App", 0, 0 );
  }
  return *evDisp_;
}

void EvDispCFT::Initialize(int RunNum)
{
  gStyle->SetPalette(1);

  tc_ = new TCanvas("canvas","CFT Event Display",1000,1000);

  tp_[0] = new TPad("pad0","PHI",0.001,0.001,0.599,0.599,10);
  tp_[0]->SetFrameFillColor(kGray);
  tp_[0]->Draw();
  tp_[1] = new TPad("pad1","",0.601,0.001,0.999,0.199,10);
  tp_[1]->Draw();
  tp_[2] = new TPad("pad2","",0.601,0.201,0.999,0.399,10);
  tp_[2]->Draw();
  tp_[3] = new TPad("pad3","",0.601,0.401,0.999,0.599,10);
  tp_[3]->Draw();
  tp_[4] = new TPad("pad4","",0.601,0.601,0.999,0.999,10);
  tp_[4]->Draw();
  tp_[5] = new TPad("pad5","",0.001,0.601,0.599,0.799,10);
  tp_[5]->Draw();
  tp_[6] = new TPad("pad6","",0.001,0.801,0.599,0.999,10);
  tp_[6]->Draw();

  tp_[0]->cd();

  //hbase_ = new TH2F("hbase","Event Display XY plane", 360, -180,180, 360, -180, 180);
  hbase_ = new TH2F("hbase","Event Display XY plane", 500, -800, 800, 500, -800, 800);
  hbase_->SetMaximum(50);
  hbase_->SetMinimum(-1);
  hbase_->Draw();

  char name[100], title[100];

  for (int i=0; i<NumOfLayersCFT; i++) {
    int segNum = NumOfSegCFT[i];

    sprintf(name, "hp%d", i);
    if (i==CFT_PHI1)
      sprintf(title, "HitPattern PHI1");
    else if (i==CFT_U1)
      sprintf(title, "HitPattern U1");
    else if (i==CFT_PHI2)
      sprintf(title, "HitPattern PHI2");
    else if (i==CFT_V2)
      sprintf(title, "HitPattern V2");
    else if (i==CFT_PHI3)
      sprintf(title, "HitPattern PHI3");
    else if (i==CFT_U3)
      sprintf(title, "HitPattern U3");
    else if (i==CFT_PHI4)
      sprintf(title, "HitPattern PHI4");
    else if (i==CFT_V4)
      sprintf(title, "HitPattern V4");

    hp_[i] = new TH1F(name, title, segNum, 0, segNum);
    //tp_[i+1]->cd();
    //hp_[i]->Draw();
  }


  tp_[4]->cd();
  hbaseU_ = new TH2F("hbaseU","Event Display Z-#phi plane", 420, -10, 410, 380, -10,370);
  hbaseU_->Draw();

  tp_[5]->cd();
  hbaseZX_ = new TH2F("hbaseZX", "Event Display ZX plane", 550, -150, 400, 140, -120, 120);
  hbaseZX_->Draw();

  tp_[6]->cd();
  hbaseZY_ = new TH2F("hbaseZY", "Event Display ZY plane", 550, -150, 400, 140, -120, 120);
  hbaseZY_->Draw();


  ConstructCFT();
  ConstructBGO();
  ConstructPiV();

  tc_->cd();
  tc_->Update();

}


void EvDispCFT::ConstructCFT(void)
{
  static const std::string funcname = "EvDispCFT::ConstructCFT";  

  tp_[0]->cd();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  for (int i=0; i<NumOfLayersCFT; i++) {
    CFT_Arc_[i].reserve(NumOfSegCFT[i]);
    for (int seg=0; seg<NumOfSegCFT[i]; seg++) {
      double x, y;
      FiberPos(i, seg, &x, &y);
      TArc *arc = new TArc(x, y, 0.75/2);	
      arc->SetLineColor(kBlack);
      arc->SetFillStyle(0);
      arc->Draw("same");
      CFT_Arc_[i].push_back(arc);
    }
  }

  /*
  tp_[4]->cd();
  const double zmax=400.;
  //for (int seg=0; seg<NumOfSegCFT_U; seg++) {
  for (int seg=0; seg<426; seg++) {
    if (seg==0) {
      U_Line_[seg][0] = new TLine(0, 0, zmax, 360);
      U_Line_[seg][0]->SetLineColor(kWhite);
      U_Line_[seg][0]->Draw("same");
      U_Line_[seg][1] = new TLine(0, 0, zmax, 360);
      U_Line_[seg][1]->SetLineColor(kWhite);
      U_Line_[seg][1]->Draw("same");
    } else {
      int lnum=65;
      double Z0 = geomMan.calcWirePosition(lnum, seg);
      printf("seg : %d, Z0 : %f\n", seg, Z0);
      double slope = geomMan.GetTiltAngle(lnum);

      double x1 = (zmax-Z0)/slope;

      U_Line_[seg][0] = new TLine(Z0, 0, zmax, x1);
      U_Line_[seg][0]->SetLineColor(kWhite);
      U_Line_[seg][0]->Draw("same");
      U_Line_[seg][1] = new TLine(0, x1, Z0, 360);
      U_Line_[seg][1]->SetLineColor(kWhite);
      U_Line_[seg][1]->Draw("same");

    }
  }  
  */
}

void EvDispCFT::ConstructBGO(void)
{
  static const std::string funcname = "EvDispCFT::ConstructBGO";  

  tp_[0]->cd();

  int unit=0;

  /*
  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5+(double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-1)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*Deg2Rad) - y1*sin(theta*Deg2Rad)),
		       (x1*sin(theta*Deg2Rad) + y1*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos2((x2*cos(theta*Deg2Rad) - y2*sin(theta*Deg2Rad)),
		       (x2*sin(theta*Deg2Rad) + y2*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos3((x3*cos(theta*Deg2Rad) - y3*sin(theta*Deg2Rad)),
		       (x3*sin(theta*Deg2Rad) + y3*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos4((x4*cos(theta*Deg2Rad) - y4*sin(theta*Deg2Rad)),
		       (x4*sin(theta*Deg2Rad) + y4*cos(theta*Deg2Rad)),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      BGO_Line_[unit].push_back(l1);
      BGO_Line_[unit].push_back(l2);
      BGO_Line_[unit].push_back(l3);
      BGO_Line_[unit].push_back(l4);
      unit++;
    }
  }
  */


  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = (double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-0.5)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*Deg2Rad) - y1*sin(theta*Deg2Rad)),
		       (x1*sin(theta*Deg2Rad) + y1*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos2((x2*cos(theta*Deg2Rad) - y2*sin(theta*Deg2Rad)),
		       (x2*sin(theta*Deg2Rad) + y2*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos3((x3*cos(theta*Deg2Rad) - y3*sin(theta*Deg2Rad)),
		       (x3*sin(theta*Deg2Rad) + y3*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos4((x4*cos(theta*Deg2Rad) - y4*sin(theta*Deg2Rad)),
		       (x4*sin(theta*Deg2Rad) + y4*cos(theta*Deg2Rad)),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = j+3*i;
      BGO_Line_[unit].push_back(l1);
      BGO_Line_[unit].push_back(l2);
      BGO_Line_[unit].push_back(l3);
      BGO_Line_[unit].push_back(l4);

      std::cout << unit << std::endl;
    }
  }

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5 + (double)i*45.;

    for (int j=0; j<NumOfBGOInOneUnit2; j++) {
      double x0 = RadiusOfBGOSurface2+BGO_Y/2;
      double y0 = (double)(j)*BGO_X;

      double x1 = x0+BGO_Y/2;
      double y1 = y0+BGO_X/2;

      double x2 = x0-BGO_Y/2;
      double y2 = y0+BGO_X/2;

      double x3 = x0-BGO_Y/2;
      double y3 = y0-BGO_X/2;

      double x4 = x0+BGO_Y/2;
      double y4 = y0-BGO_X/2;

      ThreeVector pos1((x1*cos(theta*Deg2Rad) - y1*sin(theta*Deg2Rad)),
		       (x1*sin(theta*Deg2Rad) + y1*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos2((x2*cos(theta*Deg2Rad) - y2*sin(theta*Deg2Rad)),
		       (x2*sin(theta*Deg2Rad) + y2*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos3((x3*cos(theta*Deg2Rad) - y3*sin(theta*Deg2Rad)),
		       (x3*sin(theta*Deg2Rad) + y3*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos4((x4*cos(theta*Deg2Rad) - y4*sin(theta*Deg2Rad)),
		       (x4*sin(theta*Deg2Rad) + y4*cos(theta*Deg2Rad)),
		       0);

      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");

      unit = j+NumOfBGOInOneUnit+3*i;
      BGO_Line_[unit].push_back(l1);
      BGO_Line_[unit].push_back(l2);
      BGO_Line_[unit].push_back(l3);
      BGO_Line_[unit].push_back(l4);

      std::cout << unit << std::endl;
    }
  }


}


void EvDispCFT::ConstructPiV(void)
{
  static const std::string funcname = "EvDispCFT::ConstructPiV";  

  tp_[0]->cd();


  for (int l=0; l<NumOfPiVLayer; l++) {
    double PiV_X = (RadiusOfPiVSurface + l*PiV_Y)/RadiusOfPiVSurface*PiV_X0;
    int unit=0;
    
    for (int i=0; i<NumOfPiVUnit; i++) {
      //double theta = 22.5+(double)i*45.;
      double theta = (double)i*45.;
   
      for (int j=0; j<NumOfPiVInOneUnit; j++) {
	double x0 = RadiusOfPiVSurface + l*PiV_Y + PiV_Y/2;
	double y0 = (double)(j-1)*PiV_X;

	double x1 = x0+PiV_Y/2;
	double y1 = y0+PiV_X/2;
	
	double x2 = x0-PiV_Y/2;
	double y2 = y0+PiV_X/2;
	
	double x3 = x0-PiV_Y/2;
	double y3 = y0-PiV_X/2;
	
	double x4 = x0+PiV_Y/2;
	double y4 = y0-PiV_X/2;
	
	ThreeVector pos1((x1*cos(theta*Deg2Rad) - y1*sin(theta*Deg2Rad)),
			 (x1*sin(theta*Deg2Rad) + y1*cos(theta*Deg2Rad)),
			 0);
	ThreeVector pos2((x2*cos(theta*Deg2Rad) - y2*sin(theta*Deg2Rad)),
			 (x2*sin(theta*Deg2Rad) + y2*cos(theta*Deg2Rad)),
			 0);
	ThreeVector pos3((x3*cos(theta*Deg2Rad) - y3*sin(theta*Deg2Rad)),
			 (x3*sin(theta*Deg2Rad) + y3*cos(theta*Deg2Rad)),
			 0);
	ThreeVector pos4((x4*cos(theta*Deg2Rad) - y4*sin(theta*Deg2Rad)),
			 (x4*sin(theta*Deg2Rad) + y4*cos(theta*Deg2Rad)),
			 0);
	
	TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
	TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
	TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
	TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");
	
	PiV_Line_[l][unit].push_back(l1);
	PiV_Line_[l][unit].push_back(l2);
	PiV_Line_[l][unit].push_back(l3);
	PiV_Line_[l][unit].push_back(l4);
	unit++;
      }
    }
  }

  for (int l=0; l<NumOfPiVLayer; l++) {
    double PiV2_X = (RadiusOfPiV2Surface + l*PiV2_Y)/RadiusOfPiV2Surface*PiV2_X0;
    int unit = NumOfPiVUnit * NumOfPiVInOneUnit;

    for (int i=0; i<NumOfPiVUnit; i++) {
      double theta = 22.5+(double)i*45.;

      double x0 = RadiusOfPiV2Surface + l*PiV2_Y + PiV2_Y/2;
      double y0 = 0;

      double x1 = x0+PiV2_Y/2;
      double y1 = y0+PiV2_X/2;
      
      double x2 = x0-PiV2_Y/2;
      double y2 = y0+PiV2_X/2;
      
      double x3 = x0-PiV2_Y/2;
      double y3 = y0-PiV2_X/2;
      
      double x4 = x0+PiV2_Y/2;
      double y4 = y0-PiV2_X/2;
      
      ThreeVector pos1((x1*cos(theta*Deg2Rad) - y1*sin(theta*Deg2Rad)),
		       (x1*sin(theta*Deg2Rad) + y1*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos2((x2*cos(theta*Deg2Rad) - y2*sin(theta*Deg2Rad)),
		       (x2*sin(theta*Deg2Rad) + y2*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos3((x3*cos(theta*Deg2Rad) - y3*sin(theta*Deg2Rad)),
		       (x3*sin(theta*Deg2Rad) + y3*cos(theta*Deg2Rad)),
		       0);
      ThreeVector pos4((x4*cos(theta*Deg2Rad) - y4*sin(theta*Deg2Rad)),
		       (x4*sin(theta*Deg2Rad) + y4*cos(theta*Deg2Rad)),
		       0);
      
      TLine *l1 = new TLine(pos1.x(), pos1.y(), pos2.x(), pos2.y());
      TLine *l2 = new TLine(pos2.x(), pos2.y(), pos3.x(), pos3.y());
      TLine *l3 = new TLine(pos3.x(), pos3.y(), pos4.x(), pos4.y());
      TLine *l4 = new TLine(pos4.x(), pos4.y(), pos1.x(), pos1.y());
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");
      
      PiV_Line_[l][unit].push_back(l1);
      PiV_Line_[l][unit].push_back(l2);
      PiV_Line_[l][unit].push_back(l3);
      PiV_Line_[l][unit].push_back(l4);
      unit++;
      
    }
  }
}


void EvDispCFT::ShowHitFiber(int layer, int segment, double pe) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiber";

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);

  if (segment>=0 && segment<NumOfSegCFT[layer]) {
    double x, y;
    FiberPos(layer, segment, &x, &y);
    hbase_->Fill(x, y, pe);
    
    CFT_Arc_[layer][segment]->SetLineColor(kRed);
    
    hp_[CFT_PHI1]->Fill(segment, pe);
  }
}

void EvDispCFT::ShowHitFiberProton(int layer, int segment) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiberProton";

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);

  if (segment>=0 && segment<NumOfSegCFT[layer]) {

    CFT_Arc_[layer][segment]->SetLineColor(kYellow);
  }
}

void EvDispCFT::ShowHitFiberPi(int layer, int segment) const
{
  static const std::string funcname = "EvDispCFT::ShowHitFiberProton";

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //printf("layer %d, seg %d, pe %f\n", layer, segment, pe);

  if (segment>=0 && segment<NumOfSegCFT[layer]) {

    CFT_Arc_[layer][segment]->SetLineColor(kGreen);
  }
}

void EvDispCFT::ShowHitBGO(int segment, double de) const
{
  static const std::string funcname = "EvDispCFT::ShowHitBGO";

  int size = BGO_Line_[segment].size();
  for (int i=0; i<size; i++)
    BGO_Line_[segment][i]->SetLineColor(kRed);

  double x, y;
  BGOPos(segment, &x, &y);
  hbase_->Fill(x, y, de);
}

void EvDispCFT::ShowHitBGO_Proton(int segment) const
{
  static const std::string funcname = "EvDispCFT::ShowHitBGO_Proton";

  int size = BGO_Line_[segment].size();
  for (int i=0; i<size; i++)
    BGO_Line_[segment][i]->SetLineColor(kYellow);
}

void EvDispCFT::ShowHitBGO_Pi(int segment) const
{
  static const std::string funcname = "EvDispCFT::ShowHitBGO_Pi";

  int size = BGO_Line_[segment].size();
  for (int i=0; i<size; i++)
    BGO_Line_[segment][i]->SetLineColor(kGreen);
}


void EvDispCFT::ShowHitPiV(int layer, int segment, double de) const
{
  static const std::string funcname = "EvDispCFT::ShowHitPiV";

  int size = PiV_Line_[layer][segment].size();
  for (int i=0; i<size; i++)
    PiV_Line_[layer][segment][i]->SetLineColor(kRed);
}

void EvDispCFT::ShowHitPiV_Proton(int layer, int segment) const
{
  static const std::string funcname = "EvDispCFT::ShowHitPiV";

  int size = PiV_Line_[layer][segment].size();
  for (int i=0; i<size; i++)
    PiV_Line_[layer][segment][i]->SetLineColor(kYellow);
}

void EvDispCFT::ShowHitPiV_Pi(int layer, int segment) const
{
  static const std::string funcname = "EvDispCFT::ShowHitPiV";

  int size = PiV_Line_[layer][segment].size();
  for (int i=0; i<size; i++)
    PiV_Line_[layer][segment][i]->SetLineColor(kGreen);
}


void EvDispCFT::DrawTrackInXYPlane(double x0, double y0, double x1, double y1, bool flagK, bool flagP) const
{

  TLine *l = new TLine(x0, y0, x1, y1);
  if (flagK)
    l->SetLineColor(kOrange);
  else if (flagP)
    l->SetLineColor(kYellow);
  else
    l->SetLineColor(kGreen);

  TrackXYCont.push_back(l);

}

void EvDispCFT::DrawTrackInZXPlane(double z0, double x0, double z1, double x1, bool flagK, bool flagP) const
{

  TLine *l = new TLine(z0, x0, z1, x1);
  if (flagK)
    l->SetLineColor(kOrange);
  else if (flagP)
    l->SetLineColor(kYellow);
  else
    l->SetLineColor(kGreen);

  TrackZXCont.push_back(l);

}

void EvDispCFT::DrawTrackInZYPlane(double z0, double y0, double z1, double y1, bool flagK, bool flagP) const
{

  TLine *l = new TLine(z0, y0, z1, y1);

  if (flagK)
    l->SetLineColor(kOrange);
  else if (flagP)
    l->SetLineColor(kYellow);
  else
    l->SetLineColor(kGreen);
  TrackZYCont.push_back(l);

}

void EvDispCFT::DrawVertex0(double x, double y, double z) const
{
  gr_Vertex0XY_ = new TGraph(1, &x, &y);
  gr_Vertex0XY_->SetMarkerStyle(21);
  //gr_Vertex0XY_->SetMarkerSize(0.2);
  gr_Vertex0XY_->SetMarkerColor(kRed);

  gr_Vertex0ZX_ = new TGraph(1, &z, &x);
  gr_Vertex0ZX_->SetMarkerStyle(21);
  //gr_Vertex0ZX_->SetMarkerSize(0.2);
  gr_Vertex0ZX_->SetMarkerColor(kRed);

  gr_Vertex0ZY_ = new TGraph(1, &z, &y);
  gr_Vertex0ZY_->SetMarkerStyle(21);
  //gr_Vertex0ZY_->SetMarkerSize(0.2);
  gr_Vertex0ZY_->SetMarkerColor(kRed);

}

void EvDispCFT::DrawVertex(double x, double y, double z) const
{
  gr_VertexXY_ = new TGraph(1, &x, &y);
  gr_VertexXY_->SetMarkerStyle(21);
  //gr_VertexXY_->SetMarkerSize(0.2);
  gr_VertexXY_->SetMarkerColor(kGreen);

  gr_VertexZX_ = new TGraph(1, &z, &x);
  gr_VertexZX_->SetMarkerStyle(21);
  //gr_VertexZX_->SetMarkerSize(0.2);
  gr_VertexZX_->SetMarkerColor(kGreen);

  gr_VertexZY_ = new TGraph(1, &z, &y);
  gr_VertexZY_->SetMarkerStyle(21);
  //gr_VertexZY_->SetMarkerSize(0.2);
  gr_VertexZY_->SetMarkerColor(kGreen);

}

void EvDispCFT::DrawCFTVertex(double x, double y, double z) const
{
  TGraph *gr = new TGraph(1, &x, &y);
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(kYellow);
  CFTVtxContXY.push_back(gr);

  TGraph *gr2 = new TGraph(1, &z, &x);
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerColor(kYellow);
  CFTVtxContZX.push_back(gr2);

  TGraph *gr3 = new TGraph(1, &z, &y);
  gr3->SetMarkerStyle(22);
  gr3->SetMarkerColor(kYellow);
  CFTVtxContZY.push_back(gr3);

}


void EvDispCFT::UpdateCanvas() const
{
  tp_[0]->cd();
  hbase_->Draw("colz");

  for (int layer=0; layer<NumOfLayersCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      CFT_Arc_[layer][seg]->Draw("same");
    }
  }

  //for (int seg=0; seg<NumOfBGOUnit*NumOfBGOInOneUnit; seg++) {
  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = BGO_Line_[seg].size();
    for (int i=0; i<size; i++)
      BGO_Line_[seg][i]->Draw("same");
  }

  for (int l=0; l<NumOfPiVLayer; l++) {
    for (int seg=0; seg<NumOfPiVUnit*NumOfPiVInOneUnit+NumOfPiVUnit; seg++) {
      int size = PiV_Line_[l][seg].size();
      for (int i=0; i<size; i++)
	PiV_Line_[l][seg][i]->Draw("same");
    }
  }

  

  int nt = TrackXYCont.size();
  for (int i=0; i<nt; i++) {
    TLine *l = TrackXYCont[i];
    l->Draw("same");
  }

  if (gr_Vertex0XY_) {
    gr_Vertex0XY_->Draw("p");
  }

  if (gr_VertexXY_) {
    gr_VertexXY_->Draw("p");
  }

  for (int i=0; i<CFTVtxContXY.size(); i++) {
    TGraph *gr = CFTVtxContXY[i];
    gr->Draw("p");
  }


  /*
  tp_[4]->cd();
  //for (int seg=0; seg<NumOfSegCFT_U; seg++) {
  for (int seg=0; seg<426; seg++) {
    U_Line_[seg][0]->Draw("same");
    U_Line_[seg][1]->Draw("same");
  }
  */
  /*
  for (int i=0; i<3; i++) {
    tp_[i+1]->cd();
    hp_[i]->Draw();
  }
  */

  tp_[5]->cd();
  int ntZX = TrackZXCont.size();
  for (int i=0; i<ntZX; i++) {
    TLine *l = TrackZXCont[i];
    l->Draw("same");
  }
  if (gr_Vertex0ZX_)
    gr_Vertex0ZX_->Draw("p");

  if (gr_VertexZX_)
    gr_VertexZX_->Draw("p");

  for (int i=0; i<CFTVtxContZX.size(); i++) {
    TGraph *gr = CFTVtxContZX[i];
    gr->Draw("p");
  }

  tp_[6]->cd();
  int ntZY = TrackZYCont.size();
  for (int i=0; i<ntZY; i++) {
    TLine *l = TrackZYCont[i];
    l->Draw("same");
  }

  if (gr_Vertex0ZY_)
    gr_Vertex0ZY_->Draw("p");

  if (gr_VertexZY_)
    gr_VertexZY_->Draw("p");

  for (int i=0; i<CFTVtxContZY.size(); i++) {
    TGraph *gr = CFTVtxContZY[i];
    gr->Draw("p");
  }

  tc_->cd();
  tc_->Update();
  tc_->Modified();


}

void EvDispCFT::EndOfEvent() const
{
  for (int layer=0; layer<NumOfLayersCFT; layer++) {
    for (int seg=0; seg<NumOfSegCFT[layer]; seg++) {
      CFT_Arc_[layer][seg]->SetLineColor(kBlack);
    }
  }

  //for (int seg=0; seg<NumOfBGOUnit*NumOfBGOInOneUnit; seg++) {
  for (int seg=0; seg<NumOfSegBGO; seg++) {
    int size = BGO_Line_[seg].size();
    for (int i=0; i<size; i++)
      BGO_Line_[seg][i]->SetLineColor(kBlack);
  }

  for (int l=0; l<NumOfPiVLayer; l++) {
    for (int seg=0; seg<NumOfPiVUnit*NumOfPiVInOneUnit+NumOfPiVUnit; seg++) {
      int size = PiV_Line_[l][seg].size();
      for (int i=0; i<size; i++)
	PiV_Line_[l][seg][i]->SetLineColor(kBlack);
    }
  }

  hbase_->Reset("ICES");
  for (int i=0; i<NumOfLayersCFT; i++)
    hp_[i]->Reset("ICES");

  if (gr_Vertex0XY_) {
    delete gr_Vertex0XY_;
  }

  if (gr_Vertex0ZX_) {
    delete gr_Vertex0ZX_;
  }

  if (gr_Vertex0ZY_) {
    delete gr_Vertex0ZY_;
  }

  if (gr_VertexXY_) {
    delete gr_VertexXY_;
    gr_VertexXY_ = 0;
  }

  if (gr_VertexZX_) {
    delete gr_VertexZX_;
    gr_VertexZX_ = 0;
  }

  if (gr_VertexZY_) {
    delete gr_VertexZY_;
    gr_VertexZY_ = 0;
  }

  std::for_each(TrackXYCont.begin(), TrackXYCont.end(), DeleteObject());
  TrackXYCont.clear();

  std::for_each(TrackZXCont.begin(), TrackZXCont.end(), DeleteObject());
  TrackZXCont.clear();

  std::for_each(TrackZYCont.begin(), TrackZYCont.end(), DeleteObject());
  TrackZYCont.clear();

  std::for_each(CFTVtxContXY.begin(), CFTVtxContXY.end(), DeleteObject());
  CFTVtxContXY.clear();

  std::for_each(CFTVtxContZX.begin(), CFTVtxContZX.end(), DeleteObject());
  CFTVtxContZX.clear();

  std::for_each(CFTVtxContZY.begin(), CFTVtxContZY.end(), DeleteObject());
  CFTVtxContZY.clear();

}

void EvDispCFT::FiberPos(int layer, int seg, double *x, double *y) const
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;

  if (layer == CFT_PHI1) {
    lnum = geomMan.GetDetectorId("CFT_PHI1");
  } else if (layer == CFT_PHI2) {
    lnum = geomMan.GetDetectorId("CFT_PHI2");
  } else if (layer == CFT_PHI3) {
    lnum = geomMan.GetDetectorId("CFT_PHI3");
  } else if (layer == CFT_PHI4) {
    lnum = geomMan.GetDetectorId("CFT_PHI4");
  } else if (layer == CFT_U1) {
    lnum = geomMan.GetDetectorId("CFT_U1");
  } else if (layer == CFT_V2) {
    lnum = geomMan.GetDetectorId("CFT_V2");
  } else if (layer == CFT_U3) {
    lnum = geomMan.GetDetectorId("CFT_U3");
  } else if (layer == CFT_V4) {
    lnum = geomMan.GetDetectorId("CFT_V4");
  } else {
    fprintf(stderr, "EvDispCFT::FiberPosPhi : No PHI Layer %d\n", layer);
    return;
  }

  if (layer == CFT_PHI1 || layer == CFT_PHI2 || layer == CFT_PHI3 || layer == CFT_PHI4) {
    double phi = geomMan.calcWirePosition(lnum, seg);
    double r   = geomMan.GetLocalZ(lnum);
    if (seg%2 == 0)
      r -= 0.4;
    else
      r += 0.4;

    *x = r * cos(phi*Deg2Rad);
    *y = r * sin(phi*Deg2Rad);
  } else if (layer == CFT_U1 || layer == CFT_U3) {
    double z = 200.;

    double SegNumU=NumOfSegCFT[layer];
    double phi = -(360./SegNumU)*(double)seg;
    double offset = 360./400.*z;

    phi += offset;

    double r   = geomMan.GetLocalZ(lnum);
    if (seg%2 == 0)
      r -= 0.4755/2;
    else
      r += 0.4755/2;
    
    *x = r * cos(phi*Deg2Rad);
    *y = r * sin(phi*Deg2Rad);
    
  } else if (layer == CFT_V2 || layer == CFT_V4) {
    double z = 200.;

    double SegNumV=NumOfSegCFT[layer];
    double phi = (360./SegNumV)*(double)seg;
    double offset = -360./400.*z;

    phi += offset;

    double r   = geomMan.GetLocalZ(lnum);
    if (seg%2 == 0)
      r -= 0.4755/2;
    else
      r += 0.4755/2;
    
    *x = r * cos(phi*Deg2Rad);
    *y = r * sin(phi*Deg2Rad);
    
  }
}

void EvDispCFT::BGOPos(int seg, double *x, double *y) const
{
  /*
  int UnitNum = seg/NumOfBGOInOneUnit;
  int SegInUnit = seg%NumOfBGOInOneUnit;

  double theta = 22.5+(double)UnitNum*45.;
  double x0 = RadiusOfBGOSurface+BGO_Y/2;
  double y0 = (double)(SegInUnit-1)*BGO_X;

  *x = x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad);
  *y = x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad);
  */

  /* new BGO position 2014/12/03*/
  int UnitNum = seg/(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);
  int SegInUnit = seg%(NumOfBGOInOneUnit+NumOfBGOInOneUnit2);

  double theta;
  double x0;
  double y0;

  if (SegInUnit==0 || SegInUnit==1 ) {
    theta = (double)UnitNum*45.;
    x0 = RadiusOfBGOSurface+BGO_Y/2;
    y0 = (double)((double)SegInUnit-0.5)*BGO_X;
  } else {
    theta = 22.5+(double)UnitNum*45.;
    x0 = RadiusOfBGOSurface2+BGO_Y/2;
    y0 = 0.;

  }


  *x = x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad);
  *y = x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad);

}


#if 0

void EvDispCFT::FiberPosPhi(int layer, int seg, double *x, double *y) const
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;

  if (layer == CFT_PHI1) {
    if (seg>=0 && seg<= 95)
      lnum = 61;
    else if (seg>=96 && seg<= 191)
      lnum = 62;
    else if (seg>=192 && seg<= 287)
      lnum = 63;
    else if (seg>=288 && seg<= 383)
      lnum = 64;

  } else if (layer == CFT_PHI2) {
    if (seg>=0 && seg<= 95)
      lnum = 66;
    else if (seg>=96 && seg<= 191)
      lnum = 67;
    else if (seg>=192 && seg<= 287)
      lnum = 68;
    else if (seg>=288 && seg<= 383)
      lnum = 69;
  } else {
    fprintf(stderr, "EvDispCFT::FiberPosPhi : No PHI Layer %d\n", layer);
    return;
  }

  double phi = geomMan.calcWirePosition(lnum, seg);
  double r   = geomMan.GetLocalZ(lnum);
  if (seg%2 == 0)
    r -= 0.4;
  else
    r += 0.4;

  *x = r * cos(phi*Deg2Rad);
  *y = r * sin(phi*Deg2Rad);

}

void EvDispCFT::FiberPosU(int layer, int seg, double z, double *x, double *y) const
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double phi;

  if (layer == CFT_U) {
    lnum = 65;
    
    double SegNumU=426.;
    phi = -(360./SegNumU)*(double)seg;
    double offset = 360./400.*z;

    phi += offset;
  } else {
    fprintf(stderr, "EvDispCFT::FiberPosU : No U Layer %d\n", layer);
    return;

  }

  double r   = geomMan.GetLocalZ(lnum);
  if (seg%2 == 0)
    r -= 0.4755/2;
  else
    r += 0.4755/2;

  *x = r * cos(phi*Deg2Rad);
  *y = r * sin(phi*Deg2Rad);

}

#endif

void EvDispCFT::get_command(void) const
{
  tc_->Update();

  char ch;
  char data[100];
  static int stat=0;
  static int Nevent=0;
  static int ev=0;
  
  if (stat == 1 && Nevent > 0 && ev<Nevent) {

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
