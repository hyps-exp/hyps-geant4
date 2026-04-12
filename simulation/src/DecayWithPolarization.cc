#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "globals.hh"
#include "RadDeg.hh"
#include "Randomize.hh"
#include "DecayWithPolarization.hh"
#include "DetectorID.hh"
#include "G4SystemOfUnits.hh"

using namespace CLHEP;

DecayWithPolarization::DecayWithPolarization(double m1, double m2, double m3, double pol, double asym, G4ThreeVector polVec, G4ThreeVector p1vec, G4ThreeVector posvec, bool phiFlag, double ctau, bool chargeFlag)
  : M1(m1), M2(m2), M3(m3), Pol(pol), Asym(asym), PolVec(polVec), P1Vec(p1vec), PosVec(posvec), CTau(ctau), ChargeFlag(chargeFlag)
{
  double A = M1*M1-(M2+M3)*(M2+M3);
  double B = M1*M1-(M2-M3)*(M2-M3);
  PCM = 1./(2.*M1)*sqrt(A*B);

  ThetaCM = RandAngularDistribution();
  PhiCM   = (double)RandFlat::shoot(0., 360.);

  //std::cout << "ThetaCM = " << ThetaCM << std::endl;
  double cost = cos(ThetaCM*Deg2Rad);
  double sint = sin(ThetaCM*Deg2Rad);
  double cosf = cos(PhiCM*Deg2Rad);
  double sinf = sin(PhiCM*Deg2Rad);

  P2vec_cm = G4ThreeVector(sint*sinf, cost, sint*cosf) * PCM;
  P3vec_cm = -1 * P2vec_cm;

  //std::cout << "P2cec_cm = " << P2vec_cm << std::endl;
  //std::cout << "P3cec_cm = " << P3vec_cm << std::endl;

  Lv2 = G4LorentzVector(P2vec_cm, sqrt(M2*M2+PCM*PCM));
  Lv3 = G4LorentzVector(P3vec_cm, sqrt(M3*M3+PCM*PCM));


  DecayCheck();

  G4ThreeVector boostVec(0, 0, P1Vec.mag()/sqrt(P1Vec.mag2()+M1*M1));
  G4LorentzVector Lv2_lab = Lv2.boost(boostVec);
  G4LorentzVector Lv3_lab = Lv3.boost(boostVec);

  double theta1, phi1;  // theta direction is obtained from z axis

  calcThetaPhi(P1Vec, &theta1, &phi1);
  if (!phiFlag) {
    phi1 -= 180.;
    theta1 *= -1.;
  }
  //std::cout << "theta1 = " << theta1 << ", phi1 = " << phi1 << std::endl;

  P2vec_lab = Lv2_lab.vect();
  P2vec_lab.rotateY(theta1*degree);
  P2vec_lab.rotateZ(phi1*degree);

  P3vec_lab = Lv3_lab.vect();
  P3vec_lab.rotateY(theta1*degree);
  P3vec_lab.rotateZ(phi1*degree);

  // Pol direction adjustment
  G4ThreeVector Xdir(1., 0., 0.);
  G4ThreeVector Ydir(0., 1., 0.);
  G4ThreeVector Zdir(0., 0., 1.);
  Xdir.rotateY(theta1*degree);
  Xdir.rotateZ(phi1*degree);
  Ydir.rotateY(theta1*degree);
  Ydir.rotateZ(phi1*degree);
  Zdir.rotateY(theta1*degree);
  Zdir.rotateZ(phi1*degree);

  double cosPhiY = Ydir*PolVec/(Ydir.mag()*PolVec.mag());
  double PhiY    = acos(cosPhiY)*Rad2Deg;
  if (std::abs(cosPhiY-1.) < 0.000001)
    PhiY    = 0.;

  Xdir.rotate(-PhiY*degree, Zdir);
  Ydir.rotate(-PhiY*degree, Zdir);

  if (std::abs(Ydir*PolVec/(Ydir.mag()*PolVec.mag()) - 1)<0.01) {
    P2vec_lab.rotate(-PhiY*degree, Zdir);    
    P3vec_lab.rotate(-PhiY*degree, Zdir);    
  } else {
    P2vec_lab.rotate(PhiY*degree, Zdir);    
    P3vec_lab.rotate(PhiY*degree, Zdir);    
  }

  double E2 = sqrt(P2vec_lab.mag2()+M2*M2);
  Ekin2 = E2-M2;
  double E3 = sqrt(P3vec_lab.mag2()+M3*M3);
  Ekin3 = E3-M3;

}

double DecayWithPolarization::RandAngularDistribution()
{
  int success=0;
  double x,fx;
  double theta;

  do {
    x = (double)RandFlat::shoot(-1., 1.);
    fx = 0.5*(1+Asym*Pol*x);
    if (fx >= (double)RandFlat::shoot())
      success = 1;
  } while (success==0);

  theta = acos(x)*Rad2Deg;
  //printf("cost = %f, theta = %f\n", x, theta);

  return theta;
}

void DecayWithPolarization::calcThetaPhi(G4ThreeVector vec, double *theta, double *phi)
{
  *theta = acos(vec.z()/vec.mag())*Rad2Deg;

  if (*theta<0.)
    *theta *= -1.;

  if (vec.x()>=0. && vec.y()>=0. ) {
    //std::cout << "1: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = atan(vec.y()/vec.x())*Rad2Deg;
  } else if (vec.x()<0. && vec.y()>=0. ) {
    //std::cout << "2: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = 180.+atan(vec.y()/vec.x())*Rad2Deg;
  } else if (vec.x()<0. && vec.y()<0. ) {
    //std::cout << "3: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = 360.-acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = 180.+atan(vec.y()/vec.x())*Rad2Deg;
  } else if (vec.x()>=0. && vec.y()<0. ) {
    //std::cout << "4: " << atan(vec.y()/vec.x())*Rad2Deg << std::endl;
    *phi = 360.-acos(vec.x()/(vec.mag()*sin(*theta*Deg2Rad)))*Rad2Deg;
    //*phi = 360.+atan(vec.y()/vec.x())*Rad2Deg;
  } else {
    fprintf(stderr, 
	    "PrimaryGeneratorAction::calcThetaPhi strange vector (%f, %f, %f)",
	    vec.x(), vec.y(), vec.z());
    *phi = -1.;
  }
}


void DecayWithPolarization::DecayCheck()
{

  G4double ctau=CTau; /*mm*/
  G4double p = P1Vec.mag();
  G4double m = M1;
  G4double E;

  G4double dx = 0.1*mm; //mm
  G4double totalx=0.0;  //mm

  G4ThreeVector localPos = PosVec;
  G4bool flagDecay=false;

  int nIteration=0;
  while (1) {
    int flagTgtType = -1;

    //check of the position
    flagTgtType = getTargetFlag(localPos);
    totalx += dx;
    localPos += P1Vec*(dx/P1Vec.mag());

    flagDecay=false;

    flagDecay = decayCheck(ctau, p, m, dx/mm );

    //check of the position
    flagTgtType = getTargetFlag(localPos);

    if (ChargeFlag && flagTgtType == 0 ) {
      p = calcEnergyDeposit(p, m, dx/mm, localPos, flagTgtType);
      E = sqrt(p*p + m*m);
    }

    if (fabs(p)<0.000001) {
      //return ; // temporary
      flagDecay = true;
    }

    if (flagDecay )
      break;

    nIteration++;
    if (nIteration>5000)
      return;


  }

  if (ChargeFlag )
    P1Vec *= p/P1Vec.mag();

  DecayPosVec = localPos;

}

G4bool DecayWithPolarization::decayCheck(G4double ctau, G4double momentum, G4double mass, G4double dx)
{
  G4double lambda=ctau*momentum/mass;
  G4double value;

  value = exp(-dx/lambda);

  G4double x;
  x=(G4double)RandFlat::shoot();
  if (x > value)
    return true;
  else 
    return false;  

}


int DecayWithPolarization::getTargetFlag(G4ThreeVector pos)
{
  int flagTgtType = -1;

  double r = sqrt((pos.x()/mm)*(pos.x()/mm)+(pos.y()/mm)*(pos.y()/mm));
  if ( r <= LH2TgtR/2. 
       && fabs(pos.z()/mm) <= Target_Length/2. ) {
    // inside the LH2 target;
    flagTgtType = 0;
    return flagTgtType;
  }

  return flagTgtType;	

}


G4double DecayWithPolarization::calcEnergyDeposit(G4double momentum, G4double mass, G4double dx, G4ThreeVector pos, G4int flagTgtType)
{

  G4double E0 = sqrt(momentum*momentum+mass*mass) * 1000.; //GeV --> MeV
  G4double beta = momentum/sqrt(momentum*momentum+mass*mass);
  
  G4double dEdx = calc_dE_dx(beta, flagTgtType);    // MeV/cm
  if (dEdx < 0.) {
    return -1.;
  }
  G4double dE = dEdx * dx / 10.; // dx (mm), MeV
  G4double E = E0 - dE;          // MeV
  if (E-mass*1000.<0.005) {
    printf("#particle stops in material at (%f, %f, %f)\n",pos.x(), pos.y(), pos.z());
    return 0.0;
  }

  G4double p = sqrt((E/1000.)*(E/1000.) -mass*mass); //GeV/c
  
  return p;
}

G4double DecayWithPolarization::calc_dE_dx(double beta, int flagTgtType)
{
  double value;
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511;
  double logterm;
  double gamma_2;
  double W_max;

  int z = 1;

  double rho;   /*g/cm^3*/
  double I;     /*eV*/
  double Z_A;

  if ( flagTgtType == 0 ) {
    // inside the LH2 target;
    rho = 0.0708; /* g/cm^3 */
    I   = 21.8;   /* eV */
    Z_A = 0.99212;
  } else if (flagTgtType == 1) {
    // inside the LD2 target;
    rho = 0.169; /* g/cm^3 */
    I   = 21.8;  /* eV */
    Z_A = 0.49652;

  } else {
    //fprintf(stderr, "Out of the target : (%f, %f, %f)\n", pos.x(), pos.y(), pos.z());
    return 0.0;
  }

  gamma_2=1/(1-pow(beta,2.0));
  //printf("sqrt(gamma_2):%f\n",sqrt(gamma_2));

  W_max=2.0*m_e*pow(beta,2.0)*gamma_2;

  logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0);

  value=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);

  return value;

}


