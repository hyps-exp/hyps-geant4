Double_t fitFunction(Double_t *x, Double_t *par)
{
  double P0 = 1.;
  double P1 = x[0];
  double P2 = 1./2.*(3*x[0]*x[0]-1);
  double P3 = 1./2.*(5*x[0]*x[0]*x[0]-3*x[0]);

  return (par[0]*P0+par[1]*P1+par[2]*P2+par[3]*P3)*(par[0]*P0+par[1]*P1+par[2]*P2+par[3]*P3);
  
}


void fitDiffCrossSection()
{
  gStyle->SetOptStat(0);
  
  char buf[200];

  //FILE *fp = fopen("dcs_1.421.txt","r");
  //FILE *fp = fopen("dcs_1.521.txt","r");
  //FILE *fp = fopen("dcs_1.621.txt","r");
  //FILE *fp = fopen("dcs_1.723.txt","r");
  //FILE *fp = fopen("dcs_1.824.txt","r");
  //FILE *fp = fopen("dcs_1.925.txt","r");
  //FILE *fp = fopen("dcs_2.026.txt","r");
  //FILE *fp = fopen("dcs_2.126.txt","r");
  //FILE *fp = fopen("dcs_2.227.txt","r");
  //FILE *fp = fopen("dcs_2.328.txt","r");
  FILE *fp = fopen("dcs_2.430.txt","r");    
  if (!fp) {
    std::cout << "cannot open file" << std::endl;
    return;
  }

  int MaxIndex = 20;
  double x[MaxIndex], err_x[MaxIndex];
  double y[MaxIndex], err_y[MaxIndex];

  int index=0;
  
  while (1) {
    if (fgets(buf, sizeof(buf), fp) == NULL)
      break;

    if (buf[0] == '\n' || buf[0] == '#')
      continue;

    double cost, val, err;
    sscanf(buf, "%lf %lf %lf", &cost, &val, &err);

    x[index] = cost;
    err_x[index] = 0.05;
    y[index] = val;
    err_y[index] = err;
    index++;
  }
  

  TH2F *hbase = new TH2F("hbase","", 100, -1, 1, 100, 0, 3);
  hbase->Draw();
  
  TGraphErrors *gr = new TGraphErrors(index, x, y, err_x, err_y);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kRed);
  gr->Draw("p");

  TF1 *func = new TF1("func", fitFunction, -1, 1, 4);
  func->SetParameter(0, 1);
  func->SetParameter(1, 0.5);
  func->SetParameter(2, 0.3);
  func->SetParameter(3, 0.);

  gr->Fit("func", "", "", -1, 1);


  double par[4];
  std::cout << " { ";
  for (int i=0; i<4; i++) {
    par[i] = func->GetParameter(i);
    std::cout << par[i];
    if (i!=3)
      std::cout << ", ";
  }
  std::cout << " }" << std::endl;
  
}
