const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

void calcBeamDCGeometry()
{
  double Theta = 20.;


  double deltaBDC3[6] = {-22., -18., -2.,  2., 18., 22.};
  double deltaBDC4[6] = {-22., -18., -2.,  2., 18., 22.};

  // distance from KURAMA center
  double targetPosOfs = 100.;

  double distTarget = -574.-420-400+targetPosOfs;
  double distBDC4= 500.;
  double distBDC3= 399+distBDC4;

  double distCFT= 50.;

  char buf[100];

  printf("target   %.1lf    -150.0     0.0   %.1lf\n", distTarget, distTarget-distTarget);

  double gTgt_X = distTarget;
  double gTgt_Y = -150.;
  double gTgt_Z = 0.;

  double L, gx, gy, gz;
  L = distCFT;
  gx = gTgt_X+L*cos(Theta*Deg2Rad);
  gy = gTgt_Y-L*sin(Theta*Deg2Rad);
  gz = gTgt_Z;

  printf("CFT    %.3lf    %.3f    %.3f   \n", gx, gy, gz);

  for (int i=0; i<6; i++) {
    L = -distBDC3+deltaBDC3[i];
    gx = gTgt_X+L*cos(Theta*Deg2Rad);
    gy = gTgt_Y-L*sin(Theta*Deg2Rad);
    gz = gTgt_Z;

    sprintf(buf, "BC3-%d", i);
    printf("%s    %.3lf    %.3f    %.3f   %.1lf\n", buf, gx, gy, gz, L);
  }

  for (int i=0; i<6; i++) {
    L = -distBDC4+deltaBDC4[i];
    gx = gTgt_X+L*cos(Theta*Deg2Rad);
    gy = gTgt_Y-L*sin(Theta*Deg2Rad);
    gz = gTgt_Z;

    sprintf(buf, "BC4-%d", i);
    printf("%s    %.3lf    %.3f    %.3f   %.1lf\n", buf, gx, gy, gz, L);
  }

}
