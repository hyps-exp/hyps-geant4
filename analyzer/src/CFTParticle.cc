#include "CFTParticle.hh"
#include "CFTLocalTrack.hh"
#include "Hodo2Hit.hh"
#include <TRandom.h>

const int NumOfPIDVal = 84;
double Eval[NumOfPIDVal] = {
  1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0,
  21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 37.0, 39.0,
  41.0, 43.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0, 57.0, 59.0,
  61.0, 63.0, 65.0, 67.0, 69.0, 71.0, 73.0, 75.0, 77.0, 79.0,
  81.0, 83.0, 85.0, 87.0, 89.0, 91.0, 93.0, 95.0, 97.0, 99.0,
  101.0, 103.0, 105.0, 107.0, 109.0, 111.0, 113.0, 115.0, 117.0, 119.0,
  121.0, 123.0, 125.0, 127.0, 129.0, 131.0, 133.0, 135.0, 137.0, 139.0,
  141.0, 143.0, 145.0, 147.0, 149.0, 151.0, 153.0, 155.0, 157.0, 160.0,
  164.0, 168.0, 173.0, 179.0 };

double Range1[NumOfPIDVal] = {
  1.84, 2.04, 2.13, 2.02, 1.89, 1.81, 1.73, 1.63, 1.59, 1.52,
  1.47, 1.38, 1.34, 1.25, 1.25, 1.21, 1.14, 1.10, 1.06, 1.03,
  0.98, 0.96, 0.91, 0.90, 0.90, 0.86, 0.83, 0.80, 0.79, 0.76,
  0.77, 0.74, 0.72, 0.74, 0.70, 0.66, 0.65, 0.62, 0.62, 0.61,
  0.61, 0.60, 0.61, 0.55, 0.56, 0.53, 0.50, 0.53, 0.52, 0.52,
  0.54, 0.49, 0.45, 0.51, 0.50, 0.49, 0.47, 0.46, 0.46, 0.49,
  0.41, 0.43, 0.45, 0.43, 0.41, 0.39, 0.39, 0.42, 0.41, 0.39,
  0.39, 0.38, 0.37, 0.37, 0.36, 0.37, 0.41, 0.37, 0.37, 0.37,
  0.37, 0.37, 0.37, 0.37 };

double Range2[NumOfPIDVal] = {
  5.50, 4.58, 4.11, 3.76, 3.55, 3.31, 3.14, 2.97, 2.76, 2.59,
  2.47, 2.39, 2.26, 2.20, 2.04, 1.97, 1.92, 1.86, 1.79, 1.73,
  1.68, 1.63, 1.60, 1.55, 1.48, 1.45, 1.42, 1.39, 1.34, 1.33,
  1.27, 1.26, 1.25, 1.19, 1.19, 1.18, 1.14, 1.14, 1.11, 1.09,
  1.07, 1.05, 1.02, 1.04, 1.01, 1.03, 1.02, 0.95, 0.96, 0.93,
  0.89, 0.92, 0.94, 0.88, 0.86, 0.85, 0.87, 0.87, 0.85, 0.82,
  0.85, 0.83, 0.79, 0.80, 0.81, 0.83, 0.82, 0.76, 0.77, 0.75,
  0.75, 0.75, 0.75, 0.75, 0.73, 0.73, 0.71, 0.73, 0.73, 0.73,
  0.73, 0.73, 0.73, 0.73 };


CFTParticle::CFTParticle(CFTLocalTrack *track)
  : Track_(track), Mass_(-1), BGO_E_(0.), PiV_E_(0.), FiberTotal_E_(0.), 
    FiberMax_E_(0.),
    NormalizedFiberTotal_E_(0.), NormalizedFiberMax_E_(0.), TotalE_(0.),
    //BGO_Ereso(1.0)/*1 MeV for 80 MeV*/,
    BGO_Ereso(0.5)/*1 MeV for 80 MeV*/,
    BGO_Ereso2(4.0)/* constant */,
    Fiber_Ereso(0.4)/*0.4 MeV for 4 MeV Total energy deposit*/,
    CFTVtx_(-999, -999, -999)
{
  
}

CFTParticle::~CFTParticle()
{

}

Hodo2Hit * CFTParticle::GetBGOHit(int i)
{
  if (i>=0 && i<BGOCont_.size())
    return BGOCont_[i];
  else
    return 0;
}

Hodo2Hit * CFTParticle::GetPiVHit(int i)
{
  if (i>=0 && i<PiVCont_.size())
    return PiVCont_[i];
  else
    return 0;
}

bool CFTParticle::Calculate()
{
  FiberTotal_E_ = Track_->TotalDEHiGain();
  FiberMax_E_   = Track_->MaxDEHiGain();
  PathLength_   = Track_->GetTotalPathLength();
  NormalizedFiberTotal_E_ = Track_->NormalizedTotalDEHiGain();
  NormalizedFiberMax_E_   = Track_->NormalizedMaxDEHiGain();


  double Fiber_sigma = Fiber_Ereso * sqrt(FiberTotal_E_ / 4.);
  FiberTotal_E_ += gRandom->Gaus(0., Fiber_sigma);


  int nc = BGOCont_.size();
  for (int i=0; i<nc; i++) {
    Hodo2Hit *hitp = BGOCont_[i];
    BGO_E_ += hitp->DeltaE();
  }


  double BGO_sigma = BGO_Ereso * sqrt(BGO_E_ / 80.) + BGO_Ereso2;
  BGO_E_ += gRandom->Gaus(0., BGO_sigma);
  if (BGO_E_ <0)
    BGO_E_=0.;

  int ncPiV = PiVCont_.size();
  for (int i=0; i<ncPiV; i++) {
    Hodo2Hit *hitp = PiVCont_[i];
    PiV_E_ += hitp->DeltaE();
  }

  TotalE_ = FiberTotal_E_ + BGO_E_;

  if (nc ==1 || nc>= 3) {
    double delta;
    if (checkProton(BGO_E_, &delta))
      Mass_ = 0.9382720;
    else if (checkPi(BGO_E_))
      Mass_ = 0.1395701;
  } else if (nc == 2) {
    Hodo2Hit *hitp1 = BGOCont_[0];
    double BGO_E1 = hitp1->DeltaE();
    Hodo2Hit *hitp2 = BGOCont_[1];
    double BGO_E2 = hitp2->DeltaE();

    bool flag0 = false,  flag1 = false, flag2 = false;
    double delta0, delta1, delta2;
    flag0 = checkProton(BGO_E_, &delta0);
    flag1 = checkProton(BGO_E1, &delta1);
    flag2 = checkProton(BGO_E2, &delta2);

    if (flag0) {
      Mass_ = 0.9382720;      
    } else if (flag1 && flag2) {
      Mass_ = 0.9382720;      
      if (fabs(delta1) < fabs(delta2))
	BGO_E_ = BGO_E1;
      else
	BGO_E_ = BGO_E2;
    } else if (flag1) {
      Mass_ = 0.9382720;      
      BGO_E_ = BGO_E1;
    } else if (flag2) {
      Mass_ = 0.9382720;      
      BGO_E_ = BGO_E2;
    }  else if (checkPi(BGO_E_)){
      Mass_ = 0.1395701;
    }
  } else if (nc == 0) {
    double delta;
    if (checkProton(BGO_E_, &delta))
      Mass_ = 0.9382720;
    else if (checkPi(BGO_E_) && NormalizedFiberTotal_E_<0.5)
      Mass_ = 0.1395701;

  }

  double shiftE = 32.;

  if (PiV_E_>0.2 && Mass_ < 0. ) {
    double delta;
    if (checkProton(BGO_E_ - shiftE, &delta)) {
      Mass_ = 0.9382720;
      if (nc == 1 || nc >= 3) {
	BGO_E_ -= shiftE;
      } else if (nc == 2) {
	Hodo2Hit *hitp1 = BGOCont_[0];
	double BGO_E1 = hitp1->DeltaE();
	Hodo2Hit *hitp2 = BGOCont_[1];
	double BGO_E2 = hitp2->DeltaE();

	bool flag1 = false, flag2 = false;
	double delta1, delta2;
	flag1 = checkProton(BGO_E1, &delta1);
	flag2 = checkProton(BGO_E2, &delta2);

	if (flag1 && flag2) {
	  if (fabs(delta1)<fabs(delta2))
	    BGO_E_ = BGO_E1;
	  else
	    BGO_E_ = BGO_E2;
	} else if (flag1) {
	  BGO_E_ = BGO_E1;
	} else if (flag2) {
	  BGO_E_ = BGO_E2;
	} else {
	  BGO_E_ -= shiftE;
	}
      }
    }
  }

  // assume proton stop before BGO and BGO was hit by pi
  if (Mass_ < 0. ) {
    double delta;
    double E=0.;
    if (checkProton(E, &delta)) {
      Mass_ = 0.9382720;
      BGO_E_ = 0;
    }
  }



  /*
  if (NormalizedFiberTotal_E_>=0.5)
    Mass_ = 0.9382720;
  else if (NormalizedFiberTotal_E_>=0 && NormalizedFiberTotal_E_<0.5)
    Mass_ = 0.1395701;
  */
  return true;

}

bool CFTParticle::checkProton( double BGO_E, double *delta)
{
  double cut1=0., cut2=0.;
  if (BGO_E >= 0 && BGO_E <= Eval[0]) {
    cut1 = Range1[0];
    cut2 = Range2[0];
  } else if (BGO_E >= Eval[NumOfPIDVal-1]) {
    cut1 = Range1[NumOfPIDVal-1];
    cut2 = Range2[NumOfPIDVal-1];
  } else {
    int index=-1;
    for (int i=0; i<NumOfPIDVal-1; i++) {
      if (BGO_E >= Eval[i] && BGO_E <= Eval[i+1]) {
	index = i;
	break;
      }
    }

    if (index != -1) {
      double v1 = BGO_E - Eval[index];
      double v2 = Eval[index+1] - BGO_E;
      cut1 = (v2*Range1[index]+v1*Range1[index+1])/(Eval[index+1]-Eval[index]);
      cut2 = (v2*Range2[index]+v1*Range2[index+1])/(Eval[index+1]-Eval[index]);
    }
  }

  *delta = BGO_E - (cut1+cut2)/2.;

  if (NormalizedFiberTotal_E_ >= cut1 && NormalizedFiberTotal_E_ <= cut2)
    return true;
  else 
    return false;

}

bool CFTParticle::checkPi( double BGO_E)
{
  double cut1=0.;
  if (BGO_E >= 0 && BGO_E <= Eval[0]) {
    cut1 = Range1[0];
  } else if (BGO_E >= Eval[NumOfPIDVal-1]) {
    cut1 = Range1[NumOfPIDVal-1];
  } else {
    int index=-1;
    for (int i=0; i<NumOfPIDVal-1; i++) {
      if (BGO_E >= Eval[i] && BGO_E <= Eval[i+1]) {
	index = i;
	break;
      }
    }

    if (index != -1) {
      double v1 = BGO_E - Eval[index];
      double v2 = Eval[index+1] - BGO_E;
      cut1 = (v2*Range1[index]+v1*Range1[index+1])/(Eval[index+1]-Eval[index]);
    }
  }

  if (NormalizedFiberTotal_E_ < cut1)
    return true;
  else 
    return false;

}
