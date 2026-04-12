void calcDCGeometry()
{
  double deltaSFT_UV[2] = {-16., 16.};
  double deltaSDC2[6] = {-27.2, -22., -2.6, 2.6, 22., 27.2};
  double deltaSDC3[6] = {-30., -18., -6.,  6., 18., 30.};
  double deltaSDC4[6] = {-30., -18., -6.,  6., 18., 30.};

  // distance from KURAMA center
  double targetPosOfs = 0.;

  double distTarget = -574.-420-400 + targetPosOfs;
  double distCFT = -524.-420-400 + targetPosOfs;
  double distSFTX=-956. + targetPosOfs;
  double distSFTUV=-873. + targetPosOfs;
  // KEK version
  //double distSDC2=-545.6;

  // E07 version
  double distSDC2=-630.0;
  double distSDC3= 934.;
  double distSDC4= 1467.5;

  char buf[100];

  printf("target   %.1lf    -150.0     0.0   %.1lf\n", distTarget, distTarget-distTarget);
  printf("CFT      %.1lf    -150.0     0.0   %.1lf\n", distCFT, distCFT-distTarget);
  printf("SFT-X     %.1lf    -150.0     0.0   %.1lf\n", distSFTX, distSFTX-distTarget);
  for (int i=0; i<2; i++) {
    sprintf(buf, "SFT-UV%d", i);
    printf("%s   %.1lf    -150.0     0.0   %.1lf\n", buf, distSFTUV+deltaSFT_UV[i], distSFTUV+deltaSFT_UV[i]-distTarget);
  }

  for (int i=0; i<6; i++) {
    sprintf(buf, "SDC2-%d", i);
    printf("%s    %.1lf    -150.0     0.0   %.1lf\n", buf, distSDC2+deltaSDC2[i], distSDC2+deltaSDC2[i]-distTarget);
  }

  for (int i=0; i<6; i++) {
    sprintf(buf, "SDC3-%d", i);
    printf("%s    %.1lf    60.0     0.0   %.1lf\n", buf, distSDC3+deltaSDC3[i], distSDC3+deltaSDC3[i]-distTarget);
  }

  for (int i=0; i<6; i++) {
    sprintf(buf, "SDC4-%d", i);
    printf("%s    %.1lf    60.0     0.0   %.1lf\n", buf, distSDC4+deltaSDC4[i], distSDC4+deltaSDC4[i]-distTarget);
  }
}
