#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Kinema3Resonance.hh"
#include "Legendre.hh"
#include "Randomize.hh"

#define PI 3.141592654

using namespace CLHEP;

Kinema3Resonance::Kinema3Resonance()
{
}

Kinema3Resonance::Kinema3Resonance(double m1, double m2, double m3, double m4, double m5, double m_res, double width, double p1, double p2, int DistFlag)
{
  double ECM;
  double vx_res, vy_res, vz_res; /* unit vector */
  double vx3, vy3, vz3; /* unit vector */
  double vx4, vy4, vz4; /* unit vector */
  double vx5, vy5, vz5; /* unit vector */
  double theta1, theta2; /* tmpolary; theta1 represents 1st kinematics of 
			  m1,m2,m_res, and m5.
			 theta2 represents 2nd kinematics of m_res, m3, and m4*/
  double phi3;         /* 2つめの運動学でのphi(CM system)*/ 
  double theta5,phi5;  /* 乱数をふったもの。log note p.43の座標系 */
  double theta_res;  /* 乱数をふったもの。log note p.43の座標系 */
  double Theta3,Phi3;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double Theta4,Phi4;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double Theta5,Phi5;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double Theta_res,Phi_res;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double theta;
  
  kin3.M_1 = m1;
  kin3.M_2 = m2;
  kin3.M_3 = m3;
  kin3.M_4 = m4;
  kin3.M_5 = m5;

  kin3.p_1_lab = p1;
  kin3.p_2_lab = p2;

  kin3.P_1_lab[XCOORD] = p1;
  kin3.P_1_lab[YCOORD] = 0.0;
  kin3.P_1_lab[ZCOORD] = 0.0;

  kin3.P_2_lab[XCOORD] = p2;
  kin3.P_2_lab[YCOORD] = 0.0;
  kin3.P_2_lab[ZCOORD] = 0.0;

  kin3.E_1_lab = p2E(kin3.p_1_lab, kin3.M_1);
  kin3.E_2_lab = p2E(kin3.p_2_lab, kin3.M_2);
  
  ECM = sqrt((kin3.E_1_lab+kin3.E_2_lab)*(kin3.E_1_lab+kin3.E_2_lab)
	     -(kin3.p_1_lab+kin3.p_2_lab)*(kin3.p_1_lab+kin3.p_2_lab));
  do {
    kin3.M_res = RandBreitWigner::shoot(m_res, width);
  } while (kin3.M_3+kin3.M_4 > kin3.M_res || kin3.M_res > ECM-kin3.M_5);
  
  kin1 = Kinema2Body(m1, m2, kin3.M_res, m5);
  kin1.SetMomentum(1, p1);
  kin1.SetMomentum(2, p2);
  /* to generate Theta+ efficiently */
  //do {
  //theta = (double)RandSin();
  //} while (theta<110.0);

  
  if (DistFlag==0) {
    theta = (double)RandSin();
    kin1.SetThetaCM(theta);
  } else if (DistFlag==1) {
    theta = (double)RandAngularDistribution();
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag==2) {
    theta = (double)RandAngularDistributionPik(p1);
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag==3) {
    theta = (double)RandAngularDistributionPiMinuskPlus1_3();
    // theta is Sigma angle
    kin1.SetThetaCM(theta);
  } else if (DistFlag==4) {
    theta = (double)RandAngularDistributionPd();
    kin1.SetThetaCM(theta);
  } else if (DistFlag==5) {
    theta = (double)RandAngularDistributionPiMinusP_Elastic_1_32();
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag==6) {
    theta = (double)RandAngularDistributionPiPlusP_Elastic_1_45();
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag==7) {
    theta = (double)RandAngularDistributionPiMinusP_ElasticCEX_1_32();
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag==8) {
    theta = (double)RandAngularDistributionPiMinusP_K0Lambda_1_05();
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag==9) {
    int index = 0;
    if (p1<= 1.47)
      index = 0;
    else if (p1>= 1.47 && p1>= 1.57)
      index = 1;
    else if (p1>= 1.57 && p1>= 1.67)
      index = 2;
    else if (p1>= 1.67 && p1>= 1.77)
      index = 3;
    else if (p1>= 1.77 && p1>= 1.87)
      index = 4;
    else if (p1>= 1.87 && p1>= 1.97)
      index = 5;
    else if (p1>= 1.97 && p1>= 2.07)
      index = 6;
    else if (p1>= 2.07 && p1>= 2.17)
      index = 7;
    else if (p1>= 2.17 && p1>= 2.27)
      index = 8;
    else if (p1>= 2.27 && p1>= 2.37)
      index = 9;
    else if (p1>= 2.37)
      index = 10;
    
    theta = RandAngularDistributionGammaP_KLambda(index);
    kin1.SetThetaCM(180.-theta);
  } else if (DistFlag>=100 && DistFlag<=129) {
    theta = (double)RandAngularDistributionPn(DistFlag-100);
    kin1.SetThetaCM(theta);
  } else {
    fprintf(stderr, "No such DistFlag : %d\n", DistFlag);
    exit(-1);
  }

  kin1.calc_kinema();
  phi5 = (360.0*(double)RandFlat::shoot());


  /*
  static int kaisuu=0;
  if (kaisuu%2==0)
    phi5 = 180.0;
  else
    phi5 = 0.0;
  kaisuu++;
  */
  //phi5 = 360.0*(double)rand()/(RAND_MAX+1.0);

  kin3.Theta1CM = kin1.GetThetaCM();
  kin3.Phi1     = phi5;

  /* calculate m_res */
  theta_res = kin1.GetThetaLab();
  vx_res = cos(PI*theta_res/180.0);
  vy_res = sin(PI*theta_res/180.0)*cos(PI*phi5/180.0);
  vz_res = -sin(deg2rad(theta_res))*sin(deg2rad(phi5));
  CalcDistoribution(vx_res, vy_res, vz_res, &Theta_res, &Phi_res);

  kin3.E_res_lab = kin1.GetEnergyLab(3);
  kin3.p_res_lab = kin1.GetMomentumLab(3);
  kin3.P_res_lab[0] = kin3.p_res_lab*vx_res;
  kin3.P_res_lab[1] = kin3.p_res_lab*vy_res;
  kin3.P_res_lab[2] = kin3.p_res_lab*vz_res;
  kin3.theta_res = Theta_res;
  kin3.phi_res = Phi_res;

  /* calculate m5 */
  theta5 = -kin1.GetPhiLab();
  vx5 = cos(PI*theta5/180.0);
  vy5 = sin(PI*theta5/180.0)*cos(PI*phi5/180.0);
  vz5 = -sin(deg2rad(theta5))*sin(deg2rad(phi5));
  CalcDistoribution(vx5, vy5, vz5, &Theta5, &Phi5);

  kin3.E_5_lab = kin1.GetEnergyLab(4);
  kin3.p_5_lab = kin1.GetMomentumLab(4);
  kin3.P_5_lab[0] = kin3.p_5_lab*vx5;
  kin3.P_5_lab[1] = kin3.p_5_lab*vy5;
  kin3.P_5_lab[2] = kin3.p_5_lab*vz5;
  kin3.theta5 = Theta5;
  kin3.phi5 = Phi5;

  /* calculate m3,m4 */
  /* もしM_resとm4が等しく、m3が0ならば
     resonanceの値をそのままm4にわたす。*/
  if (kin3.M_res == m4 && m3 == 0.0) {
    kin3.E_3_lab = 0.0;
    kin3.p_3_lab = 0.0;
    kin3.P_3_lab[0] = 0.0;
    kin3.P_3_lab[1] = 0.0;
    kin3.P_3_lab[2] = 0.0;
    kin3.theta3 = 0.0;
    kin3.phi3 = 0.0;

    kin3.E_4_lab = kin3.E_res_lab;
    kin3.p_4_lab =  kin3.p_res_lab;
    kin3.P_4_lab[0] = kin3.P_res_lab[0];
    kin3.P_4_lab[1] = kin3.P_res_lab[1];
    kin3.P_4_lab[2] = kin3.P_res_lab[2];
    kin3.theta4 = Theta_res;
    kin3.phi4 = Phi_res;
  } else {
    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
    kin2.SetMomentum(1, kin1.GetMomentumLab(3));
    kin2.SetMomentum(2, 0.0);
    kin2.SetThetaCM((double)RandSin());
    kin2.calc_kinema();
    
    phi3 = (360.0*(double)RandFlat::shoot());
    //phi3 = 360.0*(double)rand()/(RAND_MAX+1.0);
    
    kin3.Theta2CM = kin2.GetThetaCM();
    kin3.Phi2     = phi3;
    
    /* m3 */
    theta1 = kin1.GetThetaLab();
    theta2 = kin2.GetThetaLab();

    vx3 = cos(deg2rad(theta2))*cos(deg2rad(theta1)) - 
      sin(deg2rad(theta1))*cos(deg2rad(phi3))*sin(deg2rad(theta2));
    
    vy3 = cos(deg2rad(phi5))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) +
      cos(deg2rad(theta1))*cos(deg2rad(phi5))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
      sin(deg2rad(phi5))*sin(deg2rad(phi3))*sin(deg2rad(theta2));
    
    vz3 = -sin(deg2rad(phi5))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) -
      cos(deg2rad(theta1))*sin(deg2rad(phi5))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
      cos(deg2rad(phi5))*sin(deg2rad(phi3))*sin(deg2rad(theta2));
    
    CalcDistoribution(vx3, vy3, vz3, &Theta3, &Phi3);
    
    kin3.E_3_lab = kin2.GetEnergyLab(3);
    kin3.p_3_lab = kin2.GetMomentumLab(3);
    kin3.P_3_lab[0] = kin3.p_3_lab*vx3;
    kin3.P_3_lab[1] = kin3.p_3_lab*vy3;
    kin3.P_3_lab[2] = kin3.p_3_lab*vz3;
    kin3.theta3 = Theta3;
    kin3.phi3 = Phi3;
    
    /* m4 */
    theta1 = kin1.GetThetaLab();
    theta2 = -kin2.GetPhiLab();
    
    vx4 = cos(deg2rad(theta2))*cos(deg2rad(theta1)) - 
      sin(deg2rad(theta1))*cos(deg2rad(phi3))*sin(deg2rad(theta2));
    
    vy4 = cos(deg2rad(phi5))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) +
      cos(deg2rad(theta1))*cos(deg2rad(phi5))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
      sin(deg2rad(phi5))*sin(deg2rad(phi3))*sin(deg2rad(theta2));
    
    vz4 = -sin(deg2rad(phi5))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) -
      cos(deg2rad(theta1))*sin(deg2rad(phi5))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
    cos(deg2rad(phi5))*sin(deg2rad(phi3))*sin(deg2rad(theta2));
    
    CalcDistoribution(vx4, vy4, vz4, &Theta4, &Phi4);
    
    kin3.E_4_lab = kin2.GetEnergyLab(4);
    kin3.p_4_lab = kin2.GetMomentumLab(4);
    kin3.P_4_lab[0] = kin3.p_4_lab*vx4;
    kin3.P_4_lab[1] = kin3.p_4_lab*vy4;
    kin3.P_4_lab[2] = kin3.p_4_lab*vz4;
    kin3.theta4 = Theta4;
    kin3.phi4 = Phi4;
  }
  //Dump();
}

/* for fermi motion */
Kinema3Resonance::Kinema3Resonance(double m1, double m2, double m3, double m4, double m5, double m_res, double width, double p1, double *mom2, int *ierr)
{
  double ECM;
  double vx_res, vy_res, vz_res; /* unit vector */
  double mom[XYZ], Mom; 
  double vx3, vy3, vz3; /* unit vector */
  double vx4, vy4, vz4; /* unit vector */
  double vx5, vy5, vz5; /* unit vector */
  double theta1, theta2; /* tmpolary; theta1 represents 1st kinematics of 
			  m1,m2,m_res, and m5.
			 theta2 represents 2nd kinematics of m_res, m3, and m4*/
  double phi3;         /* 2つめの運動学でのphi(CM system)*/ 
  double theta5,phi5;  /* 乱数をふったもの。log note p.43の座標系 */
  double theta_res;  /* 乱数をふったもの。log note p.43の座標系 */
  double Theta3,Phi3;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double Theta4,Phi4;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double Theta5,Phi5;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double Theta_res,Phi_res;  /* 座標系をかえてTheta,Phiをあらわしたもの */
  double theta;
  int n_try=0;

  kin3.M_1 = m1;
  kin3.M_2 = m2;
  kin3.M_3 = m3;
  kin3.M_4 = m4;
  kin3.M_5 = m5;

  kin3.p_1_lab = p1;
  kin3.p_2_lab = sqrt(mom2[XCOORD]*mom2[XCOORD] + mom2[YCOORD]*mom2[YCOORD]
		      + mom2[ZCOORD]*mom2[ZCOORD]);
  kin3.P_1_lab[XCOORD] = p1;
  kin3.P_1_lab[YCOORD] = 0.0;
  kin3.P_1_lab[ZCOORD] = 0.0;

  kin3.P_2_lab[XCOORD] = mom2[XCOORD];
  kin3.P_2_lab[YCOORD] = mom2[YCOORD];
  kin3.P_2_lab[ZCOORD] = mom2[ZCOORD];

  kin3.E_1_lab = p2E(kin3.p_1_lab, kin3.M_1);
  kin3.E_2_lab = p2E(kin3.p_2_lab, kin3.M_2);

  ECM = sqrt(m1*m1 + m2*m2 + 2.0*(kin3.E_1_lab*kin3.E_2_lab-p1*mom2[XCOORD]));
  //printf("Ecm = %f\n", ECM);
  do {
    kin3.M_res = RandBreitWigner::shoot(m_res, width);
    n_try++;
    if (n_try>100) {
      *ierr = 1;
      return;
    }
    //printf("%d %f  %f  %f\n",n_try,kin3.M_3+kin3.M_4,kin3.M_res,ECM-kin3.M_5);
  } while (kin3.M_3+kin3.M_4 > kin3.M_res || kin3.M_res > ECM-kin3.M_5);

  kin1ver2 = Kinema2BodyVer2(m1, m2, kin3.M_res, m5, p1, mom2);

  /* calculate m_res */
  kin1ver2.GetMomentum(3, mom);
  Mom = kin1ver2.GetMomentum(3);
  CalcDistoribution(mom[XCOORD]/Mom, mom[YCOORD]/Mom, mom[ZCOORD]/Mom, &Theta_res, &Phi_res);

  kin3.E_res_lab = kin1ver2.GetEnergy(3);
  kin3.p_res_lab = kin1ver2.GetMomentum(3);
  kin3.P_res_lab[0] =mom[XCOORD];
  kin3.P_res_lab[1] =mom[YCOORD];
  kin3.P_res_lab[2] =mom[ZCOORD];
  kin3.theta_res = Theta_res;
  kin3.phi_res = Phi_res;

  /* calculate m5 */
  kin1ver2.GetMomentum(4, mom);
  Mom = kin1ver2.GetMomentum(4);
  CalcDistoribution(mom[XCOORD]/Mom, mom[YCOORD]/Mom, mom[ZCOORD]/Mom, &Theta5, &Phi5);

  kin3.E_5_lab = kin1ver2.GetEnergy(4);
  kin3.p_5_lab = kin1ver2.GetMomentum(4);
  kin3.P_5_lab[0] = mom[XCOORD];
  kin3.P_5_lab[1] = mom[YCOORD];
  kin3.P_5_lab[2] = mom[ZCOORD];
  kin3.theta5 = Theta5;
  kin3.phi5 = Phi5;

  /* calculate m3,m4 */
  /* もしM_resとm4が等しく、m3が0ならば
     resonanceの値をそのままm4にわたす。*/
  if (kin3.M_res == m4 && m3 == 0.0) {
    kin3.E_3_lab = 0.0;
    kin3.p_3_lab = 0.0;
    kin3.P_3_lab[0] = 0.0;
    kin3.P_3_lab[1] = 0.0;
    kin3.P_3_lab[2] = 0.0;
    kin3.theta3 = 0.0;
    kin3.phi3 = 0.0;

    kin3.E_4_lab = kin3.E_res_lab;
    kin3.p_4_lab =  kin3.p_res_lab;
    kin3.P_4_lab[0] = kin3.P_res_lab[0];
    kin3.P_4_lab[1] = kin3.P_res_lab[1];
    kin3.P_4_lab[2] = kin3.P_res_lab[2];
    kin3.theta4 = Theta_res;
    kin3.phi4 = Phi_res;
  } else {
    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
    kin2.SetMomentum(1, kin1ver2.GetMomentum(3));
    kin2.SetMomentum(2, 0.0);
    kin2.SetThetaCM((double)RandSin());
    kin2.calc_kinema();
    
    phi3 = (360.0*(double)RandFlat::shoot());
    //phi3 = 360.0*(double)rand()/(RAND_MAX+1.0);
    
    kin3.Theta2CM = kin2.GetThetaCM();
    kin3.Phi2     = phi3;
    
    /* m3 */
    /* まずresonanceの進行方向をx軸に取った系で考える */
    theta2 = kin2.GetThetaLab();
    /* その系で(Mon, 0, 0)をz軸にtheta2度回転、x軸にphi3度回転する */
    Mom = kin2.GetMomentumLab(3);
    mom[XCOORD] = Mom;
    mom[YCOORD] = 0.0;
    mom[ZCOORD] = 0.0;
    
    rotVector(ZCOORD, theta2, mom);
    rotVector(XCOORD, phi3, mom);

    /* 次にresonanceの進行方向が実際の座標になるように回転させる */
    theta1 = Theta_res;
    phi5 = Phi_res;
    /* z軸にtheta1, x軸にphi5回転させれば良い */
    rotVector(ZCOORD, theta1, mom);
    rotVector(XCOORD, phi5, mom);
    
    CalcDistoribution(mom[XCOORD]/Mom, mom[YCOORD]/Mom, mom[ZCOORD]/Mom, &Theta3, &Phi3);
    
    kin3.E_3_lab = kin2.GetEnergyLab(3);
    kin3.p_3_lab = kin2.GetMomentumLab(3);
    kin3.P_3_lab[0] = mom[XCOORD];
    kin3.P_3_lab[1] = mom[YCOORD];
    kin3.P_3_lab[2] = mom[ZCOORD];
    kin3.theta3 = Theta3;
    kin3.phi3 = Phi3;
    
    /* m4 */
    theta1 = Theta_res;
    theta2 = -kin2.GetPhiLab();
    /* まずresonanceの進行方向をx軸に取った系で考える */
    /* その系で(Mon, 0, 0)をz軸にtheta2(-phi)度回転、x軸にphi3度回転する */
    Mom = kin2.GetMomentumLab(4);
    mom[XCOORD] = Mom;
    mom[YCOORD] = 0.0;
    mom[ZCOORD] = 0.0;

    rotVector(ZCOORD, theta2, mom);
    rotVector(XCOORD, phi3, mom);

    /* 次にresonanceの進行方向が実際の座標になるように回転させる */
    theta1 = Theta_res;
    phi5 = Phi_res;
    /* z軸にtheta1, x軸にphi5回転させれば良い */
    rotVector(ZCOORD, theta1, mom);
    rotVector(XCOORD, phi5, mom);

    CalcDistoribution(mom[XCOORD]/Mom, mom[YCOORD]/Mom, mom[ZCOORD]/Mom, &Theta4, &Phi4);    

    
    kin3.E_4_lab = kin2.GetEnergyLab(4);
    kin3.p_4_lab = kin2.GetMomentumLab(4);
    kin3.P_4_lab[0] = mom[XCOORD];
    kin3.P_4_lab[1] = mom[YCOORD];
    kin3.P_4_lab[2] = mom[ZCOORD];
    kin3.theta4 = Theta4;
    kin3.phi4 = Phi4;
  }

  *ierr = 0;
  return;
  //Dump();
  //getchar();
}

double Kinema3Resonance::p2E(double p,double m)
{
  return sqrt(p*p + m*m);
}

void Kinema3Resonance::CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi)
{
  *theta = rag2deg(acos(unitx));

  if (unity>=0.0 && unitz>0.0) 
    *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<0.0 && unitz>=0.0) 
    *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<=0.0 && unitz<0.0) 
    *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity>0.0 && unitz<=0.0) 
    *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
  else {
    fprintf(stderr,
	  "Kinema3Resonance::CalcDistribution No such reagion unity=%f, unitz=%f\n",
	    unity, unitz);
    Dump();
    exit(1);
  }
  return;
}


double Kinema3Resonance::deg2rad(double theta) {
  return 3.141592654*theta/180.0;
}

double Kinema3Resonance::rag2deg(double rag)
{
  return 360.0 * rag/ (2.0 * 3.141592654);
}

    
double Kinema3Resonance::RandSin(void)
{
  int success=0;
  double x,fx;

  do {
    x = 180.0 * (double)RandFlat::shoot();
    //x = 180.0 * (double)rand()/(RAND_MAX+1.0);
    fx = sin(deg2rad(x));
    if (fx >= (double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  return x;
}

double Kinema3Resonance::AngularDistFunction(double x)
{
  const int index=20;
  double cost[index], val[index];

  /* diff_cross_K-p_Sigma+pi-_1.205.txt */
  cost[0]  = -0.950;  val[0]  = 1.31;
  cost[1]  = -0.850;  val[1]  = 0.99;
  cost[2]  = -0.750;  val[2]  = 0.57;
  cost[3]  = -0.650;  val[3]  = 0.42;
  cost[4]  = -0.550;  val[4]  = 0.20;
  cost[5]  = -0.450;  val[5]  = 0.10;
  cost[6]  = -0.350;  val[6]  = 0.02;
  cost[7]  = -0.250;  val[7]  = 0.05;
  cost[8]  = -0.150;  val[8]  = 0.10;
  cost[9]  = -0.050;  val[9]  = 0.22;
  cost[10] =  0.050;  val[10] = 0.20;
  cost[11] =  0.150;  val[11] = 0.63;
  cost[12] =  0.250;  val[12] = 0.41;
  cost[13] =  0.350;  val[13] = 0.44;
  cost[14] =  0.450;  val[14] = 0.53;
  cost[15] =  0.550;  val[15] = 0.71;
  cost[16] =  0.650;  val[16] = 0.74;
  cost[17] =  0.750;  val[17] = 1.55;
  cost[18] =  0.850;  val[18] = 1.98;
  cost[19] =  0.950;  val[19] = 3.26;

  for (int i=0; i<index; i++) {
    if (x>=cost[i]-0.05 && x<=cost[i]+0.05)
      return val[i];
  }
  fprintf(stderr, "Kinema3Resonance::AngularDistFunction\n");
  fprintf(stderr, "Something is wrong, x is not between -1 to 1: %f\n", x);
  exit(-1);
}

double Kinema3Resonance::RandAngularDistribution()
{
  int success=0;
  double x,fx;
  double theta;

  do {
    x = (double)RandFlat::shoot(-1, 1);
    fx = AngularDistFunction(x);
    if (fx >= 4.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

double Kinema3Resonance::AngularDistFunctionPik(double x)
{
  const int index=40;
  double cost[index], val[index];

  /* diff_cross_pi+p_Sigma+K+_1.687.txt */
  cost[0]  = -0.975;  val[0]  = 0.;
  cost[1]  = -0.925;  val[1]  = 0.;
  cost[2]  = -0.875;  val[2]  = 42.6;
  cost[3]  = -0.825;  val[3]  = 28.6;
  cost[4]  = -0.775;  val[4]  = 32.5;
  cost[5]  = -0.725;  val[5]  = 31.1;
  cost[6]  = -0.675;  val[6]  = 33.7;
  cost[7]  = -0.625;  val[7]  = 34.0;
  cost[8]  = -0.575;  val[8]  = 40.9;
  cost[9]  = -0.525;  val[9]  = 39.2;
  cost[10] = -0.475;  val[10] = 41.4;
  cost[11] = -0.425;  val[11] = 45.4;
  cost[12] = -0.375;  val[12] = 39.1;
  cost[13] = -0.325;  val[13] = 36.3;
  cost[14] = -0.275;  val[14] = 34.8;
  cost[15] = -0.225;  val[15] = 25.5;
  cost[16] = -0.175;  val[16] = 29.8;
  cost[17] = -0.125;  val[17] = 28.1;
  cost[18] = -0.075;  val[18] = 28.0;
  cost[19] = -0.025;  val[19] = 27.7;
  cost[20] = 0.025;  val[20] = 26.5;
  cost[21] = 0.075;  val[21] = 32.4;
  cost[22] = 0.125;  val[22] = 29.6;
  cost[23] = 0.175;  val[23] = 33.6;
  cost[24] = 0.225;  val[24] = 28.6;
  cost[25] = 0.275;  val[25] = 24.7;
  cost[26] = 0.325;  val[26] = 28.8;
  cost[27] = 0.375;  val[27] = 31.2;
  cost[28] = 0.425;  val[28] = 24.2;
  cost[29] = 0.475;  val[29] = 31.2;
  cost[30]  = 0.525;  val[30]  = 18.7;
  cost[31]  = 0.575;  val[31]  = 18.6;
  cost[32]  = 0.625;  val[32]  = 21.5;
  cost[33]  = 0.675;  val[33]  = 18.5;
  cost[34]  = 0.725;  val[34]  = 18.6;
  cost[35]  = 0.775;  val[35]  = 28.5;
  cost[36]  = 0.825;  val[36]  = 34.1;
  cost[37]  = 0.875;  val[37]  = 51.2;
  cost[38]  = 0.925;  val[38]  = 87.4;
  cost[39]  = 0.975;  val[39]  = 114.6;

  for (int i=0; i<index; i++) {
    if (x>=cost[i]-0.025 && x<=cost[i]+0.025)
      return val[i];
  }
  fprintf(stderr, "Kinema3Resonance::AngularDistFunctionPik\n");
  fprintf(stderr, "Something is wrong, x is not between -1 to 1: %f\n", x);
  //exit(-1);
  return -1.0;
}

double Kinema3Resonance::AngularDistFunctionPik_1282(double x)
{
  int success=0;
  double fx;
  double theta;
  
  double A[6] = {40.54, -41.90, 26.59, 15.03, 13.06, -14.20};
  double P[6];

  if( !PLGen(5,x,P) ){
    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1282 Error in PLGen\n");
    exit(-1);
  }

  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1328(double x)
{
  int success=0;
  double fx;
  double theta;
  
  //double A[6] = {33.39, -9.08, 3.11, 26.36, 12.64, -0.06};
  double A[6] = {33.04, -6.43, 1.80, 30.07, 11.13, 3.24};
  double P[6];

  if( !PLGen(5,x,P) ){
    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1328 Error in PLGen\n");
    exit(-1);
  }

  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1377(double x)
{
  int success=0;
  double fx;
  double theta;
  
  //double A[6] = {35.04, 0.03, 1.51, 37.99, 8.02, 3.22};
  double A[6] = {36.46, -4.91, 3.72, 39.99, 6.51, 4.57};
  double P[6];

  if( !PLGen(5,x,P) ){
    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1377 Error in PLGen\n");
    exit(-1);
  }

  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1419(double x)
{
  int success=0;
  double fx;
  double theta;

  //double A[6]={30.01, 28.92, -19.64, 51.71, 11.51, 14.01}; //machigai
  double A[6]={40.52, 0.42, 6.20, 49.12, 19.13, 12.98};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1419 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;
}

double Kinema3Resonance::AngularDistFunctionPik_1490(double x)
{
  int success=0;
  double fx;
  double theta;

  double A[6]={36.12, 14.59, -0.70, 42.03, 12.88, 8.23};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1490 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;
}

double Kinema3Resonance::AngularDistFunctionPik_1518(double x)
{
  int success=0;
  double fx;
  double theta;

  double A[6]={11.63, 92.06, -62.57, 82.91, -11.70, 36.57};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1518 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1582(double x)
{
  int success=0;
  double fx;
  double theta;

  double A[6]={22.40, 44.28, -20.82, 47.44, 2.86, 17.63};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1582 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1614(double x)
{
  int success=0;
  double fx;
  double theta;

  double A[6]={23.13, 45.99, -13.73, 40.99, 10.08, 20.63};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1614 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1687(double x)
{
  int success=0;
  double fx;
  double theta;

  double A[6]={6.45, 84.92, -50.06, 50.98, -5.92, 34.65};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1687 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::AngularDistFunctionPik_1712(double x)
{
  int success=0;
  double fx;
  double theta;

  double A[6]={4.53, 89.46, -50.95, 52.14, -4.57, 35.85};
  double P[6];

  if( !PLGen(5,x,P) ) {

    fprintf(stderr, " Kinema3Resonance::AngularDistFunctionPik_1712 Error in PLGen\n");
    exit(-1);
  }
  fx = 0.;
  for (int i=0; i<=5; i++)
    fx += A[i]*P[i];
  
  return fx;

}

double Kinema3Resonance::RandAngularDistributionPik(double p)
{
  int success=0;
  double x,fx;
  double theta;
  
  do {
    //x = (double)RandFlat::shoot(-0.9, 1);
    x = (double)RandFlat::shoot(-1, 1);


    if (p<1.3)
      fx = AngularDistFunctionPik_1282(x);
    else if (p>=1.3&&p<=1.35)
      fx = AngularDistFunctionPik_1328(x);
    else if (p>=1.35&&p<=1.40)
      fx = AngularDistFunctionPik_1377(x);
    else if (p>=1.40&&p<=1.45)
      fx = AngularDistFunctionPik_1419(x);
    else if (p>=1.45&&p<=1.5)
      fx = AngularDistFunctionPik_1490(x);
    else if (p>=1.5&&p<=1.55)
      fx = AngularDistFunctionPik_1518(x);
    else if (p>=1.55&&p<=1.6)
      fx = AngularDistFunctionPik_1582(x);
    else if (p>=1.6&&p<=1.65)
      fx = AngularDistFunctionPik_1614(x);
    else if (p>=1.65&&p<=1.7)
      fx = AngularDistFunctionPik_1687(x);
    else if (p>=1.7)
      fx = AngularDistFunctionPik_1712(x);

    if (fx >= 130.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

double Kinema3Resonance::AngularDistFunctionPiMinuskPlus1_3(double x)
{
  const int index=21;
  double cost[index], val[index];

  cost[0]  = -1.0;  val[0]  = 53.8;
  cost[1]  = -0.9;  val[1]  = 46.34;
  cost[2]  = -0.8;  val[2]  = 39.9;
  cost[3]  = -0.7;  val[3]  = 34.55;
  cost[4]  = -0.6;  val[4]  = 29.89;
  cost[5]  = -0.5;  val[5]  = 25.85;
  cost[6]  = -0.4;  val[6]  = 22.34;
  cost[7]  = -0.3;  val[7]  = 19.27;
  cost[8]  = -0.2;  val[8]  = 16.57;
  cost[9]  = -0.1;  val[9]  = 14.21;
  cost[10] =  0.0;  val[10] = 12.18;
  cost[11] =  0.1;  val[11] = 10.46;
  cost[12] =  0.2;  val[12] =  9.10;
  cost[13] =  0.3;  val[13] =  8.15;
  cost[14] =  0.4;  val[14] =  7.69;
  cost[15] =  0.5;  val[15] =  7.80;
  cost[16] =  0.6;  val[16] =  8.62;
  cost[17] =  0.7;  val[17] = 10.29;
  cost[18] =  0.8;  val[18] = 12.98;
  cost[19] =  0.9;  val[19] = 16.88;
  cost[20] =  1.0;  val[20] = 22.20;

  for (int i=0; i<index; i++) {
    if (x>=cost[i]-0.05 && x<=cost[i]+0.05)
      return val[i];
  }
  fprintf(stderr, "Kinema3Resonance::AngularDistFunctionPiMinuskPlus1_3\n");
  fprintf(stderr, "Something is wrong, x is not between -1 to 1: %f\n", x);
  //exit(-1);
  return -1.0;
}

double Kinema3Resonance::RandAngularDistributionPiMinuskPlus1_3()
{
  int success=0;
  double x,fx;
  double theta;

  double A[5]={19.5, 17.0, 16.3, -1.2, 2.2};
  double P[5];
  do {
    x = (double)RandFlat::shoot(-1, 1);

    if( !PLGen(4,x,P) ){
      fprintf(stderr, " Kinema3Resonance::RandAngularDistributionPiMinuskPlus1_3 Error in PLGen\n");
      exit(-1);
    }
    fx = 0.;
    for (int i=0; i<=4; i++)
      fx += A[i]*P[i];

    if (fx >= 60.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

double Kinema3Resonance::RandAngularDistributionPiMinusP_Elastic_1_32()
{
  int success=0;
  double x,fx;
  double theta;

  double A[7]={1.04, 1.80, 2.52, 2.28, 1.31, 0.61, -0.06};
  double P[7];
  do {
    x = (double)RandFlat::shoot(-1, 1);

    if( !PLGen(6,x,P) ){
      fprintf(stderr, " Kinema3Resonance::RandAngularDistributionPiMinuskPlus1_3 Error in PLGen\n");
      exit(-1);
    }
    fx = 0.;
    for (int i=0; i<=6; i++)
      fx += A[i]*P[i];

    if (fx >= 12.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

double Kinema3Resonance::RandAngularDistributionPiMinusP_ElasticCEX_1_32()
{
  int success=0;
  double x,fx;
  double theta;
    
  double A[8]={0.172, 0.035, 0.171, 0.005, 0.048, 0.026, 0.155, -0.114};
  double P[8];
  do {
    x = (double)RandFlat::shoot(-1, 1);

    if( !PLGen(7,x,P) ){
      fprintf(stderr, " Kinema3Resonance::RandAngularDistributionPiMinuskPlus1_3 Error in PLGen\n");
      exit(-1);
    }
    fx = 0.;
    for (int i=0; i<=7; i++)
      fx += A[i]*P[i];

    if (fx >= 1.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

double Kinema3Resonance::RandAngularDistributionPiPlusP_Elastic_1_45()
{
  int success=0;
  double x,fx;
  double theta;

  double A[8]={0.337, -1.266, 5.272, -1.235, -17.364, 10.548, 20.41, -2.614};

  do {
    x = (double)RandFlat::shoot(-1, 1);

    fx = 0.;
    for (int i=0; i<=7; i++)
      fx += A[i]*pow(x, (double)i);

    if (fx >= 16.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}


double Kinema3Resonance::RandAngularDistributionPiMinusP_K0Lambda_1_05()
{
  int success=0;
  double x,fx;
  double theta;

  double A[3]={1.0, 0.91, 0.25};
  double P[3];
  do {
    x = (double)RandFlat::shoot(-1, 1);

    if( !PLGen(2,x,P) ){
      fprintf(stderr, " Kinema3Resonance::RandAngularDistributionPiMinusP_K0Lambda_1_05 Error in PLGen\n");
      exit(-1);
    }
    fx = 0.;
    for (int i=0; i<=2; i++)
      fx += A[i]*P[i];

    if (fx >= 3.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

double Kinema3Resonance::RandAngularDistributionGammaP_KLambda(int index)
{
  int success=0;

  double A[11][4]={
		 { 1.10814, 0.220499, 0.166374, -0.212446 },  // Eg = 1.421 GeV
		 { 1.07897, 0.304571, 0.320897, -0.210597 },  // Eg = 1.521 GeV
		 { 0.953987, 0.305151, 0.376219, -0.329168 }, // Eg = 1.621 GeV
		 { 0.874421, 0.332706, 0.451402, -0.402187 }, // Eg = 1.723 GeV
		 { 0.81151, 0.400947, 0.460501, -0.365798 },  // Eg = 1.824 GeV
		 { 0.785012, 0.387206, 0.471763, -0.307264 }, // Eg = 1.925 GeV
		 { 0.762287, 0.359112, 0.564651, -0.240777 }, // Eg = 2.026 GeV
		 { 0.675926, 0.469273, 0.497475, -0.0572796 }, // Eg = 2.126 GeV
		 { 0.635498, 0.466133, 0.539799, 0.0307159 }, // Eg = 2.227 GeV
		 { 0.59785, 0.53025, 0.593498, 0.158359 }, // Eg = 2.328 GeV
		 { 0.585356, 0.562789, 0.557159, 0.108846 }// Eg = 2.430 GeV
  };
  
  double x = -999.;
  double P[3];
  do {
    x = (double)RandFlat::shoot(-1, 1);

    double P0 = 1.;
    double P1 = x;
    double P2 = 1./2.*(3*x*x-1);
    double P3 = 1./2.*(5*x*x*x-3*x);

    double fx = (A[index][0]*P0+A[index][1]*P1+A[index][2]*P2+A[index][3]*P3)*(A[index][0]*P0+A[index][1]*P1+A[index][2]*P2+A[index][3]*P3);
    
    if (fx >= 3.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  double theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}


double Kinema3Resonance::AngularDistFunctionPd(double x)
{
  const int index=47;
  double cost[index], val[index];

  /* diff_cross_dp_elastic scattering */
  cost[0]  =  0.985996;  val[0]   =  12.818000;
  cost[1]  =  0.980271;  val[1]   =  13.494000;
  cost[2]  =  0.946085;  val[2]   =  11.950000;
  cost[3]  =  0.935444;  val[3]   =  10.683000;
  cost[4]  =  0.927184;  val[4]   =   9.877000;
  cost[5]  =  0.905569;  val[5]   =   8.187000;
  cost[6]  =  0.881303;  val[6]   =   6.515000;
  cost[7]  =  0.853551;  val[7]   =   5.128000;
  cost[8]  =  0.817145;  val[8]   =   4.140000;
  cost[9]  =  0.783693;  val[9]   =   3.202000;
  cost[10] =  0.747798;  val[10]  =   2.582000;
  cost[11] =  0.716911;  val[11]  =   2.130000;
  cost[12] =  0.676876;  val[12]  =   1.860000;
  cost[13] =  0.633381;  val[13]  =   1.584000;
  cost[14] =  0.567844;  val[14]  =   1.202000;
  cost[15] =  0.538771;  val[15]  =   1.046000;
  cost[16] =  0.486335;  val[16]  =   0.947000;
  cost[17] =  0.466387;  val[17]  =   0.893000;
  cost[18] =  0.228351;  val[18]  =   0.522000;
  cost[19] = -0.026177;  val[19]  =   0.348000;
  cost[20] = -0.297375;  val[20]  =   0.236000;
  cost[21] = -0.487860;  val[21]  =   0.216000;
  cost[22] = -0.569280;  val[22]  =   0.219000;
  cost[23] = -0.599024;  val[23]  =   0.221000;
  cost[24] = -0.627963;  val[24]  =   0.220000;
  cost[25] = -0.661312;  val[25]  =   0.225000;
  cost[26] = -0.687088;  val[26]  =   0.229000;
  cost[27] = -0.724172;  val[27]  =   0.233000;
  cost[28] = -0.771625;  val[28]  =   0.252000;
  cost[29] = -0.794415;  val[29]  =   0.270000;
  cost[30] = -0.825113;  val[30]  =   0.287000;
  cost[31] = -0.845262;  val[31]  =   0.305000;
  cost[32] = -0.864275;  val[32]  =   0.331000;
  cost[33] = -0.877983;  val[33]  =   0.358000;
  cost[34] = -0.895712;  val[34]  =   0.386000;
  cost[35] = -0.911403;  val[35]  =   0.432000;
  cost[36] = -0.928486;  val[36]  =   0.507000;
  cost[37] = -0.941471;  val[37]  =   0.559000;
  cost[38] = -0.953191;  val[38]  =   0.653000;
  cost[39] = -0.963630;  val[39]  =   0.744000;
  cost[40] = -0.972776;  val[40]  =   0.851000;
  cost[41] = -0.980955;  val[41]  =   0.920000;
  cost[42] = -0.987414;  val[42]  =   0.989000;
  cost[43] = -0.992546;  val[43]  =   1.032000;
  cost[44] = -0.995725;  val[44]  =   1.060000;
  cost[45] = -0.998441;  val[45]  =   1.100000;
  cost[46] = -0.999962;  val[46]  =   1.115000;

  if (x>cost[0] || x<cost[46])
    return 0.0;

  for (int i=0; i<index; i++) {
    if (i==0) {
      if (x<=cost[i] && x>(cost[i]+cost[i+1])/2.)
	return val[i];
    } else if (i==46) {
      if (x<=(cost[i-1]+cost[i])/2. && x>=cost[i])
	return val[i];
    } else {
      if (x<=(cost[i-1]+cost[i])/2. && x>(cost[i]+cost[i+1])/2.)
	return val[i];
    }
  }
  fprintf(stderr, "Kinema3Resonance::AngularDistFunctionPd\n");
  fprintf(stderr, "Something is wrong, x is not between -1 to 1: %f\n", x);
  //exit(-1);
  return -1.0;
}


double Kinema3Resonance::RandAngularDistributionPd()
{
  int success=0;
  double x,fx;
  double theta;

  do {
    x = (double)RandFlat::shoot(-1, 1);
    fx = AngularDistFunctionPd(x);
    if (fx >= 120.*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}


double Kinema3Resonance::AngularDistFunctionPn(double x, int flag)
{
  const int index=17;
  double cost[index], val[30][index];

  /* diff_cross_pn_elastic scattering */
  cost[0]  =  -0.984808;
  cost[1]  =  -0.939693;
  cost[2]  =  -0.866025;
  cost[3]  =  -0.766044;
  cost[4]  =  -0.642788;
  cost[5]  =  -0.500000;
  cost[6]  =  -0.342020;
  cost[7]  =  -0.173648;
  cost[8]  =  0.000000 ;
  cost[9]  =  0.173648 ;
  cost[10] =  0.342020 ;
  cost[11] =  0.500000 ;
  cost[12] =  0.642788 ;
  cost[13] =  0.766044 ;
  cost[14] =  0.866025 ;
  cost[15] =  0.939693 ;
  cost[16] =  0.984808 ;

//E(lab)=10.000000
  val[0][0] = 78.814000;
  val[0][1] = 78.584000;
  val[0][2] = 78.238000;
  val[0][3] = 77.816000;
  val[0][4] = 77.355000;
  val[0][5] = 76.882000;
  val[0][6] = 76.416000;
  val[0][7] = 75.969000;
  val[0][8] = 75.549000;
  val[0][9] = 75.160000;
  val[0][10] = 74.805000;
  val[0][11] = 74.486000;
  val[0][12] = 74.207000;
  val[0][13] = 73.969000;
  val[0][14] = 73.778000;
  val[0][15] = 73.637000;
  val[0][16] = 73.550000;
//E(lab)=20.000000
  val[1][0] = 42.862000;
  val[1][1] = 42.344000;
  val[1][2] = 41.641000;
  val[1][3] = 40.891000;
  val[1][4] = 40.184000;
  val[1][5] = 39.556000;
  val[1][6] = 39.017000;
  val[1][7] = 38.571000;
  val[1][8] = 38.220000;
  val[1][9] = 37.966000;
  val[1][10] = 37.807000;
  val[1][11] = 37.735000;
  val[1][12] = 37.735000;
  val[1][13] = 37.790000;
  val[1][14] = 37.879000;
  val[1][15] = 37.979000;
  val[1][16] = 38.058000;
//E(lab)=30.000000
  val[2][0] = 29.570000;
  val[2][1] = 28.743000;
  val[2][2] = 27.710000;
  val[2][3] = 26.716000;
  val[2][4] = 25.861000;
  val[2][5] = 25.152000;
  val[2][6] = 24.578000;
  val[2][7] = 24.139000;
  val[2][8] = 23.844000;
  val[2][9] = 23.704000;
  val[2][10] = 23.722000;
  val[2][11] = 23.885000;
  val[2][12] = 24.165000;
  val[2][13] = 24.529000;
  val[2][14] = 24.934000;
  val[2][15] = 25.322000;
  val[2][16] = 25.614000;
//E(lab)=40.000000
  val[3][0] = 22.981000;
  val[3][1] = 21.865000;
  val[3][2] = 20.573000;
  val[3][3] = 19.426000;
  val[3][4] = 18.492000;
  val[3][5] = 17.735000;
  val[3][6] = 17.126000;
  val[3][7] = 16.672000;
  val[3][8] = 16.397000;
  val[3][9] = 16.329000;
  val[3][10] = 16.473000;
  val[3][11] = 16.817000;
  val[3][12] = 17.325000;
  val[3][13] = 17.946000;
  val[3][14] = 18.623000;
  val[3][15] = 19.277000;
  val[3][16] = 19.778000;
//E(lab)=50.000000
  val[4][0] = 19.197000;
  val[4][1] = 17.828000;
  val[4][2] = 16.349000;
  val[4][3] = 15.123000;
  val[4][4] = 14.155000;
  val[4][5] = 13.364000;
  val[4][6] = 12.718000;
  val[4][7] = 12.237000;
  val[4][8] = 11.961000;
  val[4][9] = 11.932000;
  val[4][10] = 12.162000;
  val[4][11] = 12.638000;
  val[4][12] = 13.317000;
  val[4][13] = 14.137000;
  val[4][14] = 15.029000;
  val[4][15] = 15.905000;
  val[4][16] = 16.596000;
//E(lab)=60.000000
  val[5][0] = 16.811000;
  val[5][1] = 15.228000;
  val[5][2] = 13.624000;
  val[5][3] = 12.371000;
  val[5][4] = 11.393000;
  val[5][5] = 10.577000;
  val[5][6] = 9.898000;
  val[5][7] = 9.386000;
  val[5][8] = 9.097000;
  val[5][9] = 9.086000;
  val[5][10] = 9.372000;
  val[5][11] = 9.938000;
  val[5][12] = 10.739000;
  val[5][13] = 11.705000;
  val[5][14] = 12.758000;
  val[5][15] = 13.810000;
  val[5][16] = 14.664000;
//E(lab)=70.000000
  val[6][0] = 15.206000;
  val[6][1] = 13.442000;
  val[6][2] = 11.762000;
  val[6][3] = 10.517000;
  val[6][4] = 9.542000;
  val[6][5] = 8.706000;
  val[6][6] = 7.998000;
  val[6][7] = 7.457000;
  val[6][8] = 7.151000;
  val[6][9] = 7.146000;
  val[6][10] = 7.463000;
  val[6][11] = 8.086000;
  val[6][12] = 8.970000;
  val[6][13] = 10.037000;
  val[6][14] = 11.204000;
  val[6][15] = 12.389000;
  val[6][16] = 13.378000;
//E(lab)=80.000000
  val[7][0] = 14.071000;
  val[7][1] = 12.156000;
  val[7][2] = 10.438000;
  val[7][3] = 9.220000;
  val[7][4] = 8.254000;
  val[7][5] = 7.401000;
  val[7][6] = 6.671000;
  val[7][7] = 6.106000;
  val[7][8] = 5.782000;
  val[7][9] = 5.773000;
  val[7][10] = 6.105000;
  val[7][11] = 6.759000;
  val[7][12] = 7.693000;
  val[7][13] = 8.825000;
  val[7][14] = 10.066000;
  val[7][15] = 11.347000;
  val[7][16] = 12.446000;
//E(lab)=90.000000
  val[8][0] = 13.239000;
  val[8][1] = 11.196000;
  val[8][2] = 9.466000;
  val[8][3] = 8.286000;
  val[8][4] = 7.329000;
  val[8][5] = 6.463000;
  val[8][6] = 5.717000;
  val[8][7] = 5.134000;
  val[8][8] = 4.793000;
  val[8][9] = 4.777000;
  val[8][10] = 5.110000;
  val[8][11] = 5.776000;
  val[8][12] = 6.736000;
  val[8][13] = 7.903000;
  val[8][14] = 9.187000;
  val[8][15] = 10.532000;
  val[8][16] = 11.718000;
//E(lab)=100.000000
  val[9][0] = 12.608000;
  val[9][1] = 10.458000;
  val[9][2] = 8.734000;
  val[9][3] = 7.598000;
  val[9][4] = 6.646000;
  val[9][5] = 5.769000;
  val[9][6] = 5.015000;
  val[9][7] = 4.421000;
  val[9][8] = 4.065000;
  val[9][9] = 4.039000;
  val[9][10] = 4.367000;
  val[9][11] = 5.030000;
  val[9][12] = 5.995000;
  val[9][13] = 7.176000;
  val[9][14] = 8.477000;
  val[9][15] = 9.860000;
  val[9][16] = 11.113000;
//E(lab)=110.000000
  val[10][0] = 12.119000;
  val[10][1] = 9.877000;
  val[10][2] = 8.172000;
  val[10][3] = 7.078000;
  val[10][4] = 6.128000;
  val[10][5] = 5.243000;
  val[10][6] = 4.486000;
  val[10][7] = 3.887000;
  val[10][8] = 3.519000;
  val[10][9] = 3.484000;
  val[10][10] = 3.801000;
  val[10][11] = 4.451000;
  val[10][12] = 5.407000;
  val[10][13] = 6.586000;
  val[10][14] = 7.886000;
  val[10][15] = 9.285000;
  val[10][16] = 10.587000;
//E(lab)=120.000000
  val[11][0] = 11.730000;
  val[11][1] = 9.410000;
  val[11][2] = 7.732000;
  val[11][3] = 6.676000;
  val[11][4] = 5.725000;
  val[11][5] = 4.833000;
  val[11][6] = 4.080000;
  val[11][7] = 3.480000;
  val[11][8] = 3.106000;
  val[11][9] = 3.062000;
  val[11][10] = 3.366000;
  val[11][11] = 3.996000;
  val[11][12] = 4.933000;
  val[11][13] = 6.097000;
  val[11][14] = 7.381000;
  val[11][15] = 8.780000;
  val[11][16] = 10.117000;
//E(lab)=130.000000
  val[12][0] = 11.416000;
  val[12][1] = 9.028000;
  val[12][2] = 7.380000;
  val[12][3] = 6.359000;
  val[12][4] = 5.403000;
  val[12][5] = 4.505000;
  val[12][6] = 3.761000;
  val[12][7] = 3.166000;
  val[12][8] = 2.788000;
  val[12][9] = 2.737000;
  val[12][10] = 3.029000;
  val[12][11] = 3.634000;
  val[12][12] = 4.545000;
  val[12][13] = 5.685000;
  val[12][14] = 6.944000;
  val[12][15] = 8.330000;
  val[12][16] = 9.687000;
//E(lab)=140.000000
  val[13][0] = 11.157000;
  val[13][1] = 8.709000;
  val[13][2] = 7.095000;
  val[13][3] = 6.102000;
  val[13][4] = 5.137000;
  val[13][5] = 4.236000;
  val[13][6] = 3.504000;
  val[13][7] = 2.919000;
  val[13][8] = 2.542000;
  val[13][9] = 2.487000;
  val[13][10] = 2.765000;
  val[13][11] = 3.345000;
  val[13][12] = 4.225000;
  val[13][13] = 5.336000;
  val[13][14] = 6.561000;
  val[13][15] = 7.923000;
  val[13][16] = 9.292000;
//E(lab)=150.000000
  val[14][0] = 10.939000;
  val[14][1] = 8.439000;
  val[14][2] = 6.859000;
  val[14][3] = 5.890000;
  val[14][4] = 4.911000;
  val[14][5] = 4.008000;
  val[14][6] = 3.293000;
  val[14][7] = 2.722000;
  val[14][8] = 2.349000;
  val[14][9] = 2.292000;
  val[14][10] = 2.559000;
  val[14][11] = 3.112000;
  val[14][12] = 3.959000;
  val[14][13] = 5.037000;
  val[14][14] = 6.225000;
  val[14][15] = 7.556000;
  val[14][16] = 8.927000;
//E(lab)=160.000000
  val[15][0] = 10.754000;
  val[15][1] = 8.207000;
  val[15][2] = 6.661000;
  val[15][3] = 5.709000;
  val[15][4] = 4.716000;
  val[15][5] = 3.812000;
  val[15][6] = 3.116000;
  val[15][7] = 2.562000;
  val[15][8] = 2.196000;
  val[15][9] = 2.140000;
  val[15][10] = 2.397000;
  val[15][11] = 2.925000;
  val[15][12] = 3.737000;
  val[15][13] = 4.781000;
  val[15][14] = 5.928000;
  val[15][15] = 7.222000;
  val[15][16] = 8.588000;
//E(lab)=170.000000
  val[16][0] = 10.593000;
  val[16][1] = 8.006000;
  val[16][2] = 6.492000;
  val[16][3] = 5.552000;
  val[16][4] = 4.541000;
  val[16][5] = 3.639000;
  val[16][6] = 2.965000;
  val[16][7] = 2.430000;
  val[16][8] = 2.074000;
  val[16][9] = 2.022000;
  val[16][10] = 2.271000;
  val[16][11] = 2.774000;
  val[16][12] = 3.552000;
  val[16][13] = 4.561000;
  val[16][14] = 5.667000;
  val[16][15] = 6.921000;
  val[16][16] = 8.275000;
//E(lab)=180.000000
  val[17][0] = 10.452000;
  val[17][1] = 7.828000;
  val[17][2] = 6.345000;
  val[17][3] = 5.413000;
  val[17][4] = 4.383000;
  val[17][5] = 3.482000;
  val[17][6] = 2.832000;
  val[17][7] = 2.320000;
  val[17][8] = 1.975000;
  val[17][9] = 1.929000;
  val[17][10] = 2.173000;
  val[17][11] = 2.654000;
  val[17][12] = 3.399000;
  val[17][13] = 4.373000;
  val[17][14] = 5.438000;
  val[17][15] = 6.648000;
  val[17][16] = 7.986000;
//E(lab)=190.000000
  val[18][0] = 10.326000;
  val[18][1] = 7.669000;
  val[18][2] = 6.216000;
  val[18][3] = 5.287000;
  val[18][4] = 4.237000;
  val[18][5] = 3.339000;
  val[18][6] = 2.714000;
  val[18][7] = 2.225000;
  val[18][8] = 1.894000;
  val[18][9] = 1.856000;
  val[18][10] = 2.097000;
  val[18][11] = 2.558000;
  val[18][12] = 3.271000;
  val[18][13] = 4.212000;
  val[18][14] = 5.237000;
  val[18][15] = 6.402000;
  val[18][16] = 7.720000;
//E(lab)=200.000000
  val[19][0] = 10.211000;
  val[19][1] = 7.525000;
  val[19][2] = 6.100000;
  val[19][3] = 5.171000;
  val[19][4] = 4.099000;
  val[19][5] = 3.207000;
  val[19][6] = 2.607000;
  val[19][7] = 2.143000;
  val[19][8] = 1.828000;
  val[19][9] = 1.799000;
  val[19][10] = 2.039000;
  val[19][11] = 2.482000;
  val[19][12] = 3.165000;
  val[19][13] = 4.075000;
  val[19][14] = 5.061000;
  val[19][15] = 6.181000;
  val[19][16] = 7.475000;
//E(lab)=210.000000
  val[20][0] = 10.105000;
  val[20][1] = 7.394000;
  val[20][2] = 5.996000;
  val[20][3] = 5.062000;
  val[20][4] = 3.970000;
  val[20][5] = 3.082000;
  val[20][6] = 2.510000;
  val[20][7] = 2.070000;
  val[20][8] = 1.772000;
  val[20][9] = 1.754000;
  val[20][10] = 1.996000;
  val[20][11] = 2.423000;
  val[20][12] = 3.079000;
  val[20][13] = 3.959000;
  val[20][14] = 4.909000;
  val[20][15] = 5.983000;
  val[20][16] = 7.252000;
//E(lab)=220.000000
  val[21][0] = 10.006000;
  val[21][1] = 7.274000;
  val[21][2] = 5.901000;
  val[21][3] = 4.960000;
  val[21][4] = 3.846000;
  val[21][5] = 2.965000;
  val[21][6] = 2.419000;
  val[21][7] = 2.005000;
  val[21][8] = 1.724000;
  val[21][9] = 1.718000;
  val[21][10] = 1.964000;
  val[21][11] = 2.379000;
  val[21][12] = 3.009000;
  val[21][13] = 3.862000;
  val[21][14] = 4.777000;
  val[21][15] = 5.806000;
  val[21][16] = 7.048000;
//E(lab)=230.000000
  val[22][0] = 9.912000;
  val[22][1] = 7.161000;
  val[22][2] = 5.812000;
  val[22][3] = 4.862000;
  val[22][4] = 3.727000;
  val[22][5] = 2.854000;
  val[22][6] = 2.335000;
  val[22][7] = 1.947000;
  val[22][8] = 1.683000;
  val[22][9] = 1.690000;
  val[22][10] = 1.942000;
  val[22][11] = 2.346000;
  val[22][12] = 2.953000;
  val[22][13] = 3.782000;
  val[22][14] = 4.665000;
  val[22][15] = 5.650000;
  val[22][16] = 6.863000;
//E(lab)=240.000000
  val[23][0] = 9.823000;
  val[23][1] = 7.056000;
  val[23][2] = 5.729000;
  val[23][3] = 4.769000;
  val[23][4] = 3.613000;
  val[23][5] = 2.748000;
  val[23][6] = 2.256000;
  val[23][7] = 1.892000;
  val[23][8] = 1.647000;
  val[23][9] = 1.668000;
  val[23][10] = 1.926000;
  val[23][11] = 2.322000;
  val[23][12] = 2.909000;
  val[23][13] = 3.716000;
  val[23][14] = 4.570000;
  val[23][15] = 5.512000;
  val[23][16] = 6.695000;
//E(lab)=250.000000
  val[24][0] = 9.736000;
  val[24][1] = 6.956000;
  val[24][2] = 5.651000;
  val[24][3] = 4.678000;
  val[24][4] = 3.503000;
  val[24][5] = 2.646000;
  val[24][6] = 2.181000;
  val[24][7] = 1.842000;
  val[24][8] = 1.614000;
  val[24][9] = 1.650000;
  val[24][10] = 1.917000;
  val[24][11] = 2.306000;
  val[24][12] = 2.875000;
  val[24][13] = 3.663000;
  val[24][14] = 4.490000;
  val[24][15] = 5.391000;
  val[24][16] = 6.544000;
//E(lab)=260.000000
  val[25][0] = 9.652000;
  val[25][1] = 6.861000;
  val[25][2] = 5.577000;
  val[25][3] = 4.590000;
  val[25][4] = 3.396000;
  val[25][5] = 2.549000;
  val[25][6] = 2.110000;
  val[25][7] = 1.795000;
  val[25][8] = 1.585000;
  val[25][9] = 1.635000;
  val[25][10] = 1.911000;
  val[25][11] = 2.297000;
  val[25][12] = 2.851000;
  val[25][13] = 3.621000;
  val[25][14] = 4.425000;
  val[25][15] = 5.286000;
  val[25][16] = 6.408000;
//E(lab)=270.000000
  val[26][0] = 9.569000;
  val[26][1] = 6.770000;
  val[26][2] = 5.506000;
  val[26][3] = 4.505000;
  val[26][4] = 3.293000;
  val[26][5] = 2.455000;
  val[26][6] = 2.042000;
  val[26][7] = 1.750000;
  val[26][8] = 1.557000;
  val[26][9] = 1.623000;
  val[26][10] = 1.909000;
  val[26][11] = 2.293000;
  val[26][12] = 2.833000;
  val[26][13] = 3.590000;
  val[26][14] = 4.372000;
  val[26][15] = 5.195000;
  val[26][16] = 6.286000;
//E(lab)=280.000000
  val[27][0] = 9.487000;
  val[27][1] = 6.682000;
  val[27][2] = 5.437000;
  val[27][3] = 4.422000;
  val[27][4] = 3.192000;
  val[27][5] = 2.366000;
  val[27][6] = 1.977000;
  val[27][7] = 1.708000;
  val[27][8] = 1.532000;
  val[27][9] = 1.612000;
  val[27][10] = 1.910000;
  val[27][11] = 2.294000;
  val[27][12] = 2.823000;
  val[27][13] = 3.567000;
  val[27][14] = 4.330000;
  val[27][15] = 5.118000;
  val[27][16] = 6.177000;
//E(lab)=290.000000
  val[28][0] = 9.406000;
  val[28][1] = 6.598000;
  val[28][2] = 5.371000;
  val[28][3] = 4.340000;
  val[28][4] = 3.095000;
  val[28][5] = 2.280000;
  val[28][6] = 1.915000;
  val[28][7] = 1.667000;
  val[28][8] = 1.507000;
  val[28][9] = 1.602000;
  val[28][10] = 1.912000;
  val[28][11] = 2.298000;
  val[28][12] = 2.817000;
  val[28][13] = 3.551000;
  val[28][14] = 4.299000;
  val[28][15] = 5.052000;
  val[28][16] = 6.080000;
//E(lab)=300.000000
  val[29][0] = 9.325000;
  val[29][1] = 6.515000;
  val[29][2] = 5.307000;
  val[29][3] = 4.261000;
  val[29][4] = 3.001000;
  val[29][5] = 2.197000;
  val[29][6] = 1.855000;
  val[29][7] = 1.628000;
  val[29][8] = 1.483000;
  val[29][9] = 1.593000;
  val[29][10] = 1.915000;
  val[29][11] = 2.304000;
  val[29][12] = 2.817000;
  val[29][13] = 3.543000;
  val[29][14] = 4.277000;
  val[29][15] = 4.997000;
  val[29][16] = 5.995000;
  

  for (int i=0; i<index; i++) {
    if (i==0) {
      if (x < cost[i] )
	return val[flag][i];
      else if ( x >= cost[i] && x < cost[i+1])
	return  (val[flag][i+1]*(x-cost[i]) + val[flag][i]*(cost[i+1]-x))/(cost[i+1]-cost[i]);
    }  else if (i>=1 && i<=index-2) {
      if ( x >= cost[i] && x < cost[i+1])
	return  (val[flag][i+1]*(x-cost[i]) + val[flag][i]*(cost[i+1]-x))/(cost[i+1]-cost[i]);
    } else if (i==index-1) {
      if (x >= cost[i])
	return val[flag][i];
    }
  }
  fprintf(stderr, "Kinema3Resonance::AngularDistFunctionPn\n");
  fprintf(stderr, "Something is wrong, x is not between -1 to 1: %f\n", x);
  //exit(-1);
  return -1.0;
}


double Kinema3Resonance::RandAngularDistributionPn(int flag)
{
  int success=0;
  double x,fx;
  double theta;

  double max_cs[30] = { 80., 43., 30., 23., 20., 17., 16., 15., 14., 13.,
			13., 12., 12., 12., 11., 11., 11., 11., 11., 11., 
			11., 11., 10., 10., 10., 10., 10., 10., 10., 10.};
  do {
    x = (double)RandFlat::shoot(-1, 1);
    fx = AngularDistFunctionPn(x, flag);
    if (fx >= max_cs[flag]*(double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*180./PI;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

void Kinema3Resonance::Dump(void)
{
  printf("======Kinema3Resonance Dump======\n");
  printf("--Particle1--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_1, kin3.p_1_lab, kin3.E_1_lab);
  printf("--Particle2--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_2, kin3.p_2_lab, kin3.E_2_lab);
  printf("--Resonance--\n");
  printf("mass=%f\n",kin3.M_res);
  printf("--Particle3--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_3, kin3.p_3_lab, kin3.E_3_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_3_lab[0], kin3.P_3_lab[1], kin3.P_3_lab[2],
	 kin3.theta3, kin3.phi3);
  printf("--Particle4--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_4, kin3.p_4_lab, kin3.E_4_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_4_lab[0], kin3.P_4_lab[1], kin3.P_4_lab[2],
	 kin3.theta4, kin3.phi4);
  printf("--Particle5--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_5, kin3.p_5_lab, kin3.E_5_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_5_lab[0], kin3.P_5_lab[1], kin3.P_5_lab[2],
	 kin3.theta5, kin3.phi5);

  printf("Energy:E1+E2=%f, E1+E2+E3=%f\n", kin3.E_1_lab+kin3.E_2_lab,
	 kin3.E_3_lab+kin3.E_4_lab+kin3.E_5_lab);
  printf("Momentum: x-> p1+p2=%f, p3+p4+p5=%f\n", 
	 kin3.P_1_lab[0]+kin3.P_2_lab[0],
	 kin3.P_3_lab[0]+kin3.P_4_lab[0]+kin3.P_5_lab[0]);
  printf("          y-> p1+p2=%f p3+p4+p5=%f\n",
	 kin3.P_1_lab[1]+kin3.P_2_lab[1],
	 kin3.P_3_lab[1]+kin3.P_4_lab[1]+kin3.P_5_lab[1]);
  printf("          z-> p1+p2=%f p3+p4+p5=%f\n",
	 kin3.P_1_lab[2]+kin3.P_2_lab[2],
	 kin3.P_3_lab[2]+kin3.P_4_lab[2]+kin3.P_5_lab[2]);

  return;
}

double Kinema3Resonance::GetEnergy(int i)
{
  switch (i) {
  case 1:
    return kin3.E_1_lab;
    break;
  case 2:
    return kin3.E_2_lab;
    break;
  case 3:
    return kin3.E_3_lab;
    break;
  case 4:
    return kin3.E_4_lab;
    break;
  case 5:
    return kin3.E_5_lab;
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetEnergy No such particle %d\n", i);
    exit(1);
  }
}

double Kinema3Resonance::GetMomentum(int i)
{
  switch (i) {
  case 1:
    return kin3.p_1_lab;
    break;
  case 2:
    return kin3.p_2_lab;
    break;
  case 3:
    return kin3.p_3_lab;
    break;
  case 4:
    return kin3.p_4_lab;
    break;
  case 5:
    return kin3.p_5_lab;
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

void Kinema3Resonance::GetMomentum(int i, double *mom)
{
  switch (i) {
  case 1:
    mom[0] = kin3.p_1_lab;
    mom[1] = 0.0;
    mom[2] = 0.0;
    break;
  case 2:
    mom[0] = kin3.p_2_lab;
    mom[1] = 0.0;
    mom[2] = 0.0;
    break;
  case 3:
    mom[0] = kin3.P_3_lab[0];
    mom[1] = kin3.P_3_lab[1];
    mom[2] = kin3.P_3_lab[2];
    break;
  case 4:
    mom[0] = kin3.P_4_lab[0];
    mom[1] = kin3.P_4_lab[1];
    mom[2] = kin3.P_4_lab[2];
    break;
  case 5:
    mom[0] = kin3.P_5_lab[0];
    mom[1] = kin3.P_5_lab[1];
    mom[2] = kin3.P_5_lab[2];
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

double Kinema3Resonance::GetTheta(int i)
{
  switch (i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return kin3.theta3;
    break;
  case 4:
    return kin3.theta4;
    break;
  case 5:
    return kin3.theta5;
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetTheta No such particle %d\n", i);
    exit(1);
  }
}

double Kinema3Resonance::GetPhi(int i)
{
  switch (i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return kin3.phi3;
    break;
  case 4:
    return kin3.phi4;
    break;
  case 5:
    return kin3.phi5;
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetPhi No such particle %d\n", i);
    exit(1);
  }
}

double Kinema3Resonance::GetThetaCM(int i)
{
  switch (i) {
  case 1:
    return kin3.Theta1CM;
    break;
  case 2:
    return kin3.Theta2CM;
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetThetaCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

double Kinema3Resonance::GetPhiCM(int i)
{
  switch (i) {
  case 1:
    return kin3.Phi1;
    break;
  case 2:
    return kin3.Phi2;
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::GetPhiCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}


void Kinema3Resonance::RotateMom(int i, double deg, double *mom)
{
  double Sin,Cos;

  Sin=sin(deg2rad(deg));
  Cos=cos(deg2rad(deg));
  switch (i) {
  case 3:
    mom[0] = Cos*kin3.P_3_lab[0] - Sin*kin3.P_3_lab[1];
    mom[1] = Sin*kin3.P_3_lab[0] + Cos*kin3.P_3_lab[1];
    mom[2] = kin3.P_3_lab[2];
    break;
  case 4:
    mom[0] = Cos*kin3.P_4_lab[0] - Sin*kin3.P_4_lab[1];
    mom[1] = Sin*kin3.P_4_lab[0] + Cos*kin3.P_4_lab[1];
    mom[2] = kin3.P_4_lab[2];
    break;
  case 5:
    mom[0] = Cos*kin3.P_5_lab[0] - Sin*kin3.P_5_lab[1];
    mom[1] = Sin*kin3.P_5_lab[0] + Cos*kin3.P_5_lab[1];
    mom[2] = kin3.P_5_lab[2];
    break;
  default:
    fprintf(stderr, "Kinema3Resonance::RotateMom should be 3,4,5 ->%d\n",i);
    exit(1);
  }

}

double Kinema3Resonance::GetResMass(void)
{
  return kin3.M_res;
}

void Kinema3Resonance::rotVector(int coord,double theta, double *p)
{
  double matrix[3][3];
  double tmp_vec[XYZ];
  int i,j;

  switch (coord) {
  case XCOORD:
    matrix[0][0] = 1.0; matrix[0][1] = 0.0; matrix[0][2] = 0.0;
    matrix[1][0] = 0.0; matrix[1][1] = cos(deg2rad(theta)); matrix[1][2] = -sin(deg2rad(theta));
    matrix[2][0] = 0.0; matrix[2][1] = sin(deg2rad(theta)); matrix[2][2] = cos(deg2rad(theta));
    break;
  case YCOORD:
    matrix[0][0] = cos(deg2rad(theta)); matrix[0][1] = 0.0; matrix[0][2] = sin(deg2rad(theta));    
    matrix[1][0] = 0.0; matrix[1][1] = 1.0; matrix[1][2] = 0.0;
    matrix[2][0] = -sin(deg2rad(theta)); matrix[2][1] = 0.0; matrix[2][2] = cos(deg2rad(theta));    
    break;
  case ZCOORD:
    matrix[0][0] = cos(deg2rad(theta)); matrix[0][1] = -sin(deg2rad(theta)); matrix[0][2] = 0.0;    
    matrix[1][0] = sin(deg2rad(theta)); matrix[1][1] = cos(deg2rad(theta)); matrix[1][2] = 0.0;    
    matrix[2][0] = 0.0; matrix[2][1] = 0.0; matrix[2][2] = 1.0;
    break;
  default:
    fprintf(stderr,"No such coordinate-%d\n", coord);
    break;
  }
  /*
  for (i=0; i<XYZ; i++) {
    for (j=0; j<XYZ; j++) {
      printf("%f  ", matrix[i][j]);
    }
    printf("\n");
  }
  */
  for (i=0; i<XYZ; i++) {
    tmp_vec[i] = 0.0;
    for (j=0; j<XYZ; j++) {
      tmp_vec[i] += matrix[i][j]*p[j];
    }
  }

  for (i=0; i<XYZ; i++) {
    p[i] = tmp_vec[i];
  }

}
