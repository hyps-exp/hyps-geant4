const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

void calcBeamDCGeometryShift3()
{
  double Theta = 0.;


  double deltaBDC3[6] = {-22., -18., -2.,  2., 18., 22.};
  double deltaBDC4[6] = {-22., -18., -2.,  2., 18., 22.};

  // distance from KURAMA center
  double targetPosOfs = 0.;

  double distTarget = -1387.;
  double distBDC4= 732.;
  double distBDC3= 397.+distBDC4;

  double distCFT= 50.;

  double gTgt_X = distTarget;
  //double gTgt_Y = -320;   // End Gurad Position 
  //double gTgt_Y = -260;   // E07End Gurad Position 
  double gTgt_Y = -240;   // E07End Gurad Position 
  double gTgt_Z = 0.;

  char buf[100];

  printf("target   %.1lf    %.1lf     0.0   %.1lf\n", distTarget, gTgt_Y, distTarget-distTarget);


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
