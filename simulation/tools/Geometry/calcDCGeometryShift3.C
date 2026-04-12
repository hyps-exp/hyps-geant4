void calcDCGeometryShift3()
{
  const double Deg2Rad = acos(-1)/180.;

  double beam_pos = -240.; // E40 beam pos

  double deltaSFT_UV[2] = {-16., 16.};
  double deltaSDC2[6] = {-27.2, -22., -2.6, 2.6, 22., 27.2};
  double deltaSDC3[6] = {-30., -18., -6.,  6., 18., 30.};
  double deltaSDC4[6] = {-30., -18., -6.,  6., 18., 30.};
  double deltaSH1[4] = {-22.5, -7.5, 7.5, 22.5};
  double deltaSH2[4] = {-22.5, -7.5, 7.5, 22.5};

  // distance from KURAMA center

  double distTarget = -1387.;
  double distCFT = distTarget + 50.;
  double distSFTUV = -971;
  double distSFTX  = -888;
  // KEK version
  //double distSDC2=-545.6;
  // E07 version
  double distSDC2=-633.0;


  double distSDC3=  980.;
  double distSDC4= 1230.;

  double distSH1=  distSDC3 - 100.;
  double distSH2=  distSDC4 + 100.;


  /* fixed value from beam position*/
  double xposSFT_UV = 100.;
  double xposSFT_X = 140.;
  double xposSDC2 = 220.;
  /* fixed value from KURAMA center */
  double xposSDC3 = 60.;
  double xposSDC4 = 60.;
  //double xposSH1 = xposSDC3-450-48.;
  double xposSH1 = xposSDC3-450-24.;
  double xposSH2 = xposSDC3-610+24.;


  double tiltSFT_X = 0.;
  double tiltSFT_UV[2] = { 45., -45.};
  double tiltSDC2[6] = {15., 15., -15., -15, 0., 0.};
  double tiltSDC3[6] = {-30., 0., 30., -30, 0., 30.};
  double tiltSDC4[6] = {-30., 0., 30., -30, 0., 30.};
  double tiltSH1[4] = {0., 0., 90., 90.};
  double tiltSH2[4] = {90., 90., 0., 0.};

  char buf[100];

  printf("target   %.1lf    %.1lf     0.0   %.1lf\n", distTarget, beam_pos, distTarget-distTarget);
  printf("CFT      %.1lf    %.1lf     0.0   %.1lf\n", distCFT,  beam_pos, distCFT-distTarget);

  
  double offset;

  for (int i=0; i<2; i++) {
    //offset = (xposSFT_UV-beam_pos)*cos(tiltSFT_UV[i]*Deg2Rad);
    offset = xposSFT_UV*cos(tiltSFT_UV[i]*Deg2Rad);

    sprintf(buf, "SFT-UV%d", i);
    printf("%s   %.1lf    %.1lf     0.0   %.1lf     %.2lf\n", buf, distSFTUV+deltaSFT_UV[i], beam_pos, distSFTUV+deltaSFT_UV[i]-distTarget, offset);
  }

  //offset = (xposSFT_X-beam_pos)*cos(tiltSFT_X*Deg2Rad);
  offset = xposSFT_X*cos(tiltSFT_X*Deg2Rad);

  printf("SFT-X     %.1lf    %.1lf     0.0   %.1lf    %.2lf\n", distSFTX,  beam_pos, distSFTX-distTarget, offset);

  for (int i=0; i<6; i++) {
    //offset = (xposSDC2-beam_pos)*cos(tiltSDC2[i]*Deg2Rad);
    offset = xposSDC2*cos(tiltSDC2[i]*Deg2Rad);

    sprintf(buf, "SDC2-%d", i);
    printf("%s    %.1lf    %.1lf     0.0   %.1lf    %.2lf\n", buf, distSDC2+deltaSDC2[i],  beam_pos, distSDC2+deltaSDC2[i]-distTarget, offset);
  }

  for (int i=0; i<4; i++) {
    offset = (xposSH1-beam_pos)*cos(tiltSH1[i]*Deg2Rad);

    sprintf(buf, "SH1-%d", i);
    printf("%s    %.1lf    %.1lf     0.0   %.1lf     %.2lf\n", buf, distSH1+deltaSH1[i],  beam_pos, distSH1+deltaSH1[i]-distTarget, offset);
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

  for (int i=0; i<4; i++) {
    offset = (xposSH2-beam_pos)*cos(tiltSH2[i]*Deg2Rad);

    sprintf(buf, "SH2-%d", i);
    printf("%s    %.1lf    %.1lf     0.0   %.1lf     %.2lf\n", buf, distSH2+deltaSH2[i],  beam_pos, distSH2+deltaSH2[i]-distTarget, offset);
  }

}
