#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Kinema2BodyVer2.hh"
#include "Randomize.hh"

using namespace CLHEP;

Kinema2BodyVer2::Kinema2BodyVer2()
{
}

Kinema2BodyVer2::Kinema2BodyVer2(double m1, double m2, double m3, double m4,
				 double mom1, double *mom2)
{

  double val;
  double beta_vec[XYZ];
  double fourvec[4];

  /* input */
  M1 =   m1;
  M2 =   m2;
  M3 =   m3;
  M4 =   m4;

  P1 = p1[XCOORD] = mom1;
  p1[YCOORD] = p1[ZCOORD] = 0.0;

  p2[XCOORD] = mom2[XCOORD];
  p2[YCOORD] = mom2[YCOORD];
  p2[ZCOORD] = mom2[ZCOORD];
  P2 = sqrt(p2[XCOORD]*p2[XCOORD] + p2[YCOORD]*p2[YCOORD]
	    + p2[ZCOORD]*p2[ZCOORD]);

  if (P2 == 0.0) {
    theta_target = 0.0;
    phi_target = 0.0;
  } else
    CalcDistoribution(p2[XCOORD]/P2, p2[YCOORD]/P2, p2[ZCOORD]/P2, 
		      &theta_target, &phi_target);

  theta_fermi = 180.0-theta_target;

  E1 = calcE(M1, P1);  /* 質量と運動量からエネルギーを求める関数 */
  E2 = calcE(M2, P2);   

  /* まずbeta, gammaを求める */
  /* ターゲットの運動方向を一般的にしたので注意 */
  tan_alpha = (P1+P2)/(P1-P2)*tan(deg2rad(theta_fermi/2.0));
  if (tan_alpha<0)
    alpha_cm = 180.0 + rad2deg(atan(tan_alpha));
  else
    alpha_cm = rad2deg(atan(tan_alpha));

  beta_cm = ((P1+P2)/(E1+E2))*(sin(deg2rad(theta_fermi/2.0))/sin(deg2rad(alpha_cm)));
  gamma_cm = 1.0/sqrt(1.0-beta_cm*beta_cm);

  /* 粒子１、２の重心の運動方向に水平、垂直で分ける */
  p_para1 = P1*cos(deg2rad(alpha_cm-theta_fermi/2.0));
  p_parp1 = P1*sin(deg2rad(alpha_cm-theta_fermi/2.0));
  p_para2 = -P2*cos(deg2rad(alpha_cm+theta_fermi/2.0));
  p_parp2 = P2*sin(deg2rad(alpha_cm+theta_fermi/2.0));

  /* 粒子１、２の重心系での運動量、エネルギーを求める */
  E_cm1 = gamma_cm*E1 - beta_cm*gamma_cm*p_para1;
  p_cm_para1 = -beta_cm*gamma_cm*E1 + gamma_cm*p_para1;
  p_cm_parp1 = p_parp1;
  p_cm1 = sqrt(p_cm_para1*p_cm_para1 + p_cm_parp1*p_cm_parp1);

  E_cm2 = gamma_cm*E2 - beta_cm*gamma_cm*p_para2;
  p_cm_para2 = -beta_cm*gamma_cm*E2 + gamma_cm*p_para2;
  p_cm_parp2 = p_parp2;
  p_cm2 = sqrt(p_cm_para2*p_cm_para2 + p_cm_parp2*p_cm_parp2);

  val = p_cm_parp1/p_cm_para1;
  if (val<0)
    theta_axis = 180.0 + rad2deg(atan(val));
  else
    theta_axis = rad2deg(atan(val));
  
  theta_rot = 360.0*(double)RandFlat::shoot();

  /* 重心系のエネルギーを求める */
  E_cm = E_cm1 + E_cm2;

  /* 散乱粒子の重心系でのエネルギーを求める */
  E_cm3 = (E_cm + (M3*M3 - M4*M4)/E_cm)/2.0;
  E_cm4 = (E_cm - (M3*M3 - M4*M4)/E_cm)/2.0;

  P_cm3 = P_cm4 = calcP(M3, E_cm3); /* 質量とエネルギーから運動量を求める関数 */

  theta_cm = (double)RandSin();
  /* 粒子３、４の重心の運動方向に垂直、水平方向に分ける */
  p_cm_para3 = P_cm3*cos(deg2rad(theta_cm));
  p_cm_parp3 = P_cm3*sin(deg2rad(theta_cm));
  p_cm_para4 = P_cm4*cos(deg2rad(180.0-theta_cm));
  p_cm_parp4 = -P_cm4*sin(deg2rad(180.0-theta_cm));
  
  
  /* 衝突軸に対して回転する */
  p_cm3[XCOORD] = p_cm_para3;
  p_cm3[YCOORD] = p_cm_parp3;
  p_cm3[ZCOORD] = 0.0;
  /* Z軸に対して-theta_axis回転させる */
  rotVector(ZCOORD,-theta_axis, p_cm3);
  /* X軸に対してtheta_rot回転 */
  rotVector(XCOORD,theta_rot, p_cm3);
  /* Z軸に対してtheta_axis回転させる */
  rotVector(ZCOORD,theta_axis, p_cm3);
  
  p_cm4[XCOORD] = p_cm_para4;
  p_cm4[YCOORD] = p_cm_parp4;
  p_cm4[ZCOORD] = 0.0;
  /* Z軸に対して-theta_axis回転させる */
  rotVector(ZCOORD,-theta_axis, p_cm4);
  /* X軸に対してtheta_rot回転 */
  rotVector(XCOORD,theta_rot, p_cm4);
  /* Z軸に対してtheta_axis回転させる */
  rotVector(ZCOORD,theta_axis, p_cm4);

  /* lorentz変換して実験室系に戻す */
  beta_vec[XCOORD] = -beta_cm;
  beta_vec[YCOORD] = 0.0;
  beta_vec[ZCOORD] = 0.0;
  
  fourvec[0] = E_cm3;
  fourvec[1] = p_cm3[XCOORD];
  fourvec[2] = p_cm3[YCOORD];
  fourvec[3] = p_cm3[ZCOORD];
  
  LorentzBoost(beta_vec, fourvec);
  
  E3 = fourvec[0];
  p3[XCOORD] = fourvec[1];
  p3[YCOORD] = fourvec[2];
  p3[ZCOORD] = fourvec[3];
  P3 = sqrt(p3[XCOORD]*p3[XCOORD] + p3[YCOORD]*p3[YCOORD] +
	    p3[ZCOORD]*p3[ZCOORD]);
  /* betaの方向がx軸であった系から本当の系に戻る */
  rotVector(ZCOORD,alpha_cm-theta_fermi/2.0, p3);    
  
  fourvec[0] = E_cm4;
  fourvec[1] = p_cm4[XCOORD];
  fourvec[2] = p_cm4[YCOORD];
  fourvec[3] = p_cm4[ZCOORD];
  
  LorentzBoost(beta_vec, fourvec);
  
  E4 = fourvec[0];
  p4[XCOORD] = fourvec[1];
  p4[YCOORD] = fourvec[2];
  p4[ZCOORD] = fourvec[3];
  P4 = sqrt(p4[XCOORD]*p4[XCOORD] + p4[YCOORD]*p4[YCOORD] +
	    p4[ZCOORD]*p4[ZCOORD]);

  /* betaの方向がx軸であった系から本当の系に戻る */
  rotVector(ZCOORD,alpha_cm-theta_fermi/2.0, p4);    

  /* 最後にターゲットのphiだけ回転する */
  rotVector(XCOORD,phi_target, p3);    
  rotVector(XCOORD,phi_target, p4);    
  /*  
  printf("---LAB---\n");
  printf("p1:%f E1:%f\n", P1, E1);
  printf("p2:%f E2:%f\n", P2, E2);
  printf("p3:%f E3:%f\n", P3, E3);
  printf("p4:%f E4:%f\n", P4, E4);
  printf("p1 (x,y,z)=(%f, 0.0, 0.0)\n", P1);
  printf("p2 (x,y,z)=(%f, %f, %f)\n", p2[XCOORD], p2[YCOORD], p2[ZCOORD]);
  printf("p3 (x,y,z)=(%f, %f, %f)\n", p3[XCOORD], p3[YCOORD], p3[ZCOORD]);
  printf("p4 (x,y,z)=(%f, %f, %f)\n", p4[XCOORD], p4[YCOORD], p4[ZCOORD]);
  printf("---CM---\n");    
  printf("alpha_cm:%f theta_axis:%f\n", alpha_cm, theta_axis);
  printf("p1:%f E1:%f\n", p_cm1, E_cm1);
  printf("p2:%f E2:%f\n", p_cm2, E_cm2);
  printf("p3:%f E3:%f\n", P_cm3, E_cm3);
  printf("p4:%f E4:%f\n", P_cm4, E_cm4);
  printf("p_cm3 (x,y,z)=(%f, %f, %f)\n",p_cm3[XCOORD], p_cm3[YCOORD], p_cm3[ZCOORD]);
  printf("p_cm4 (x,y,z)=(%f, %f, %f)\n",p_cm4[XCOORD], p_cm4[YCOORD], p_cm4[ZCOORD]);
  getchar();
  */
}

void Kinema2BodyVer2::GetMomentum(int i,double *mom)
{
  switch (i) {
  case 1:
    mom[XCOORD] = p1[XCOORD];
    mom[YCOORD] = p1[YCOORD];
    mom[ZCOORD] = p1[ZCOORD];
    return;
    break;
  case 2:
    mom[XCOORD] = p2[XCOORD];
    mom[YCOORD] = p2[YCOORD];
    mom[ZCOORD] = p2[ZCOORD];
    return;
    break;
  case 3:
    mom[XCOORD] = p3[XCOORD];
    mom[YCOORD] = p3[YCOORD];
    mom[ZCOORD] = p3[ZCOORD];
    return;
    break;
  case 4:
    mom[XCOORD] = p4[XCOORD];
    mom[YCOORD] = p4[YCOORD];
    mom[ZCOORD] = p4[ZCOORD];
    return;
    break;
  default:
    fprintf(stderr,"Kinema2BodyVer2::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

double Kinema2BodyVer2::GetMomentum(int i)
{
  switch (i) {
  case 1:
    return P1;
    break;
  case 2:
    return P2;
    break;
  case 3:
    return P3;
    break;
  case 4:
    return P4;
    break;
  default:
    fprintf(stderr,"Kinema2BodyVer2::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

double Kinema2BodyVer2::GetEnergy(int i)
{
  switch (i) {
  case 1:
    return E1;
    break;
  case 2:
    return E2;
    break;
  case 3:
    return E3;
    break;
  case 4:
    return E4;
    break;
  default:
    fprintf(stderr,"Kinema2BodyVer2::GetEnergy No such particle %d\n", i);
    exit(1);
  }
}

#define  PI 3.141592653

double Kinema2BodyVer2::calcE(double M, double p)
{
  return sqrt(M*M + p*p);
}

double Kinema2BodyVer2::calcP(double M, double E)
{
  return sqrt(E*E -M*M);
}

double Kinema2BodyVer2::deg2rad(double deg)
{
  return deg*PI/180.0;
}

double Kinema2BodyVer2::rad2deg(double rad)
{
  return rad*180.0/PI;
}

void Kinema2BodyVer2::rotVector(int coord,double theta, double *p)
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

void Kinema2BodyVer2::LorentzBoost(double *beta, double *p)
{
  int i;
  double Beta, Gamma;
  double tmp_vec[4];

  Beta = sqrt(beta[XCOORD]*beta[XCOORD]+ beta[YCOORD]*beta[YCOORD]+
	      beta[ZCOORD]*beta[ZCOORD]);
  Gamma = 1.0/sqrt(1.0-Beta*Beta);

  tmp_vec[0] = Gamma*(p[0] - naiseki(beta, &p[1]));
  for (i=0; i<3; i++)
    tmp_vec[i+1] = p[i+1] + (Gamma-1.0)/(Beta*Beta)*naiseki(beta, &p[1])*beta[i]
      - Gamma*beta[i]*p[0];

  for (i=0; i<4; i++)
    p[i] = tmp_vec[i];

}

double Kinema2BodyVer2::naiseki(double *x, double *y)
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

void Kinema2BodyVer2::CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi)
{
  *theta = rad2deg(acos(unitx));
  
  if (*theta==0.0) {
    *phi = 0.0;
    return;
  }

  if (unity>=0.0 && unitz>0.0) 
    *phi = rad2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<0.0 && unitz>=0.0) 
    *phi = rad2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<=0.0 && unitz<0.0) 
    *phi = 360.0-rad2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity>0.0 && unitz<=0.0) 
    *phi = 360.0-rad2deg(acos(unity/sin(deg2rad(*theta))));
  else {
    fprintf(stderr,
	  "Kinema2BodyVer2::CalcDistribution No such reagion unity=%f, unitz=%f\n",
	    unity, unitz);
    exit(1);
  }
  
  return;
}

double Kinema2BodyVer2::RandSin(void)
{
  int success=0;
  double x,fx;

  do {
    x = 180.0 * (double)RandFlat::shoot();
    fx = sin(deg2rad(x));
    if (fx >= (double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  return x;
}
