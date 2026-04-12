#ifndef __KINEMA2BODYVER2_HH__
#define __KINEMA2BODYVER2_HH__

#include "common.h"

class Kinema2BodyVer2 {
private:
  double M1, M2, M3, M4;  /* 粒子の質量 */
  double p1[XYZ], p2[XYZ], p3[XYZ], p4[XYZ];  /* 粒子の実験室系での運動量 */     
  double P1, P2, P3, P4;   /* 運動量の大きさ */
  double E1, E2, E3, E4;  /* 粒子の実験室系でのエネルギー */ 
  double p_para1, p_para2, p_para3, p_para4; 
                     /* 粒子の重心系の運動方向に水平な成分 */     
  double p_parp1, p_parp2, p_parp3, p_parp4; 
                     /* 粒子の重心系の運動方向に垂直な成分 */     

  double theta_fermi;     /* fermi運動している粒子のビームとなす角 */

  double p_cm1, p_cm2, p_cm3[XYZ], p_cm4[XYZ]; /* 粒子の重心系での運動量 */     
  double P_cm3, P_cm4;   /* 運動量の大きさ */
  double E_cm1, E_cm2, E_cm3, E_cm4; /* 粒子の重心系でのエネルギー */ 
  double theta_cm;		     /* 重心系での散乱角 */           
  double p_cm_para1, p_cm_para2, p_cm_para3, p_cm_para4; 
                     /* 粒子の重心系の運動方向に水平な成分 */     
  double p_cm_parp1, p_cm_parp2, p_cm_parp3, p_cm_parp4; 
                     /* 粒子の重心系の運動方向に垂直な成分 */     
  double theta_axis;  /* 重心系での粒子１.２の運動方向と重心の運動とのなす角 */
  double theta_rot;   /* 粒子1,2の運動方向に関しての回転角 */

  double beta_cm, gamma_cm;          /* 重心系の速度 */
  double alpha_cm, tan_alpha;
  double E_cm;
  double theta_target, phi_target;

public:
  Kinema2BodyVer2();
  Kinema2BodyVer2(double m1, double m2, double m3, double m4,
		  double mom1, double *mom2);
  void   GetMomentum(int i,double *mom);
  double GetMomentum(int i);
  double GetEnergy(int i);
  double calcE(double M, double p);
  double calcP(double M, double E);
  double deg2rad(double deg);
  double rad2deg(double rad);
  void rotVector(int coord,double theta, double *p);
  void LorentzBoost(double *beta, double *p);
  double naiseki(double *x, double *y);
  void CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi);
  double RandSin(void);
};
#endif
