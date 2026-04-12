#ifndef __KINEMA3RESONANCE_HH__
#define __KINEMA3RESONANCE_HH__

#include "Kinema2Body.hh"
#include "Kinema2BodyVer2.hh"

struct KINEMA_3RESONANCE{
  double E_1_lab;
  double p_1_lab;
  double M_1;
  double P_1_lab[3];

  double E_2_lab;
  double p_2_lab;
  double M_2;
  double P_2_lab[3];

  double E_3_lab;
  double p_3_lab;
  double M_3;
  double P_3_lab[3];
  double theta3, phi3;

  double E_4_lab;
  double p_4_lab;
  double M_4;
  double P_4_lab[3];
  double theta4, phi4;

  double E_5_lab;
  double p_5_lab;
  double M_5;
  double P_5_lab[3];
  double theta5, phi5;

  double E_res_lab;
  double p_res_lab;
  double M_res;
  double P_res_lab[3];
  double theta_res, phi_res;

  double Theta1CM, Theta2CM;
  double Phi1,Phi2;
};

class Kinema3Resonance {
private:

  Kinema2Body kin1;
  Kinema2Body kin2;
  Kinema2BodyVer2 kin1ver2;
  KINEMA_3RESONANCE kin3;

public:
  Kinema3Resonance();
  Kinema3Resonance(double m1, double m2, double m3, double m4, double m5, double m_res, double width, double p1, double p2, int DistFlag=0);
  Kinema3Resonance(double m1, double m2, double m3, double m4, double m5, double m_res, double width, double p1, double *mom2, int *ierr);
  double p2E(double p,double m);
  void CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi);
  double deg2rad(double theta);
  double rag2deg(double rag);
  double RandSin(void);
  double AngularDistFunction(double x);
  double RandAngularDistribution();
  double AngularDistFunctionPik(double x);
  double RandAngularDistributionPik(double p);
  double AngularDistFunctionPik_1282(double x);
  double AngularDistFunctionPik_1328(double x);
  double AngularDistFunctionPik_1377(double x);
  double AngularDistFunctionPik_1419(double x);
  double AngularDistFunctionPik_1490(double x);
  double AngularDistFunctionPik_1518(double x);
  double AngularDistFunctionPik_1582(double x);
  double AngularDistFunctionPik_1614(double x);
  double AngularDistFunctionPik_1687(double x);
  double AngularDistFunctionPik_1712(double x);
  double AngularDistFunctionPiMinuskPlus1_3(double x);
  double RandAngularDistributionPiMinuskPlus1_3();
  double RandAngularDistributionPiMinusP_Elastic_1_32();
  double RandAngularDistributionPiMinusP_ElasticCEX_1_32();
  double RandAngularDistributionPiPlusP_Elastic_1_45();
  double RandAngularDistributionPiMinusP_K0Lambda_1_05();
  double RandAngularDistributionGammaP_KLambda(int index);
  double AngularDistFunctionPd(double x);
  double RandAngularDistributionPd();
  double AngularDistFunctionPn(double x, int flag);
  double RandAngularDistributionPn(int flag);
  void Dump(void);
  double GetEnergy(int i);
  double GetMomentum(int i);
  void GetMomentum(int i, double *mom);
  double GetTheta(int i);
  double GetPhi(int i);
  double GetThetaCM(int i);
  double GetPhiCM(int i);
  void  RotateMom(int i, double deg, double *mom);
  double GetResMass(void);
  void rotVector(int coord,double theta, double *p);
};

#endif
