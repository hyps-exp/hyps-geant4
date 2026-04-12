void calcDCGeometryShift()
{
  const double Deg2Rad = acos(-1)/180.;

  double beam_pos = -320.;

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


  double xposSFT_X = -180.;
  double xposSFT_UV = -180.;
  double xposSDC2 = -120.;
  double xposSDC3 = 60.;
  double xposSDC4 = 60.;
  double tiltSFT_X = 0.;
  double tiltSFT_UV[2] = { 45., -45.};
  double tiltSDC2[6] = {15., 15., -15., -15, 0., 0.};
  double tiltSDC3[6] = {-30., 0., 30., -30, 0., 30.};
  double tiltSDC4[6] = {-30., 0., 30., -30, 0., 30.};

  char buf[100];

  printf("target   %.1lf    %.1lf     0.0   %.1lf\n", distTarget, beam_pos, distTarget-distTarget);
  printf("CFT      %.1lf    %.1lf     0.0   %.1lf\n", distCFT,  beam_pos, distCFT-distTarget);

  
  double offset;
  offset = (xposSFT_X-beam_pos)*cos(tiltSFT_X*Deg2Rad);

  printf("SFT-X     %.1lf    %.1lf     0.0   %.1lf    %.2lf\n", distSFTX,  beam_pos, distSFTX-distTarget, offset);
  for (int i=0; i<2; i++) {
    offset = (xposSFT_UV-beam_pos)*cos(tiltSFT_UV[i]*Deg2Rad);

    sprintf(buf, "SFT-UV%d", i);
    printf("%s   %.1lf    %.1lf     0.0   %.1lf     %.2lf\n", buf, distSFTUV+deltaSFT_UV[i], beam_pos, distSFTUV+deltaSFT_UV[i]-distTarget, offset);
  }

  for (int i=0; i<6; i++) {
    offset = (xposSDC2-beam_pos)*cos(tiltSDC2[i]*Deg2Rad);

    sprintf(buf, "SDC2-%d", i);
    printf("%s    %.1lf    %.1lf     0.0   %.1lf    %.2lf\n", buf, distSDC2+deltaSDC2[i],  beam_pos, distSDC2+deltaSDC2[i]-distTarget, offset);
  }

  for (int i=0; i<6; i++) {
    offset = (xposSDC3-beam_pos)*cos(tiltSDC3[i]*Deg2Rad);

    sprintf(buf, "SDC3-%d", i);
    printf("%s    %.1lf    %.1lf     0.0   %.1lf     %.2lf\n", buf, distSDC3+deltaSDC3[i],  beam_pos, distSDC3+deltaSDC3[i]-distTarget, offset);
  }

  for (int i=0; i<6; i++) {
    offset = (xposSDC4-beam_pos)*cos(tiltSDC4[i]*Deg2Rad);

    sprintf(buf, "SDC4-%d", i);
    printf("%s    %.1lf    %.1f     0.0   %.1lf    %.2lf\n", buf, distSDC4+deltaSDC4[i], beam_pos, distSDC4+deltaSDC4[i]-distTarget, offset);
  }
}
