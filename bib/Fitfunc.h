#ifndef Fitfunc
#define Fitfunc

const string lnames[4]={"G+Pol","CB+Pol","G+Exp","CB+Exp"};

TF1* g43(TH1* h, double x1, double x2, double x3, double x4, double x5, double x6, double low, double high){

  TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + [3]*(x-[1])*(x-[1]) - [4]*(x-[1]) + [5]",
  low, high);

  TF1* back = new TF1(((TString)"background_fit")," [1]*(x-[0])*(x-[0]) - [2]*(x-[0]) + [3]",low, high);

  TF1* gauss = new TF1(((TString)"peak"),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
  low, high);


  f->SetParameters(x1, x2, x3, x4,x5,x6);
  f->SetParLimits (0,x1*0.1,x1*15);
  f->SetParLimits (1,3.08,3.12);
  f->SetParLimits (2,0.025,0.15);
  f->SetParLimits (3,0,100000);
  f->SetParLimits (5,0,100000);

  f->SetParNames("X1", "Mean", "sigma", "X4","X5","X6");

  h->Fit(f, "MEQ", "same", low, high);

  double par[6];
  f->GetParameters(par);
  back->SetParameters(par[1], par[3], par[4],par[5]);
  gauss->SetParameters(par[0], par[1], par[2]);


  back->SetLineColor(kBlue);
  gauss->SetLineColor(kGreen);

  back->Draw("same");
  gauss->Draw("same");

  double minB,maxB;
  minB=f->GetParameter(1)-3*f->GetParameter(2);
  maxB=f->GetParameter(1)+3*f->GetParameter(2);

  //number of events under the peak (3sigmas)
  double n3=h->Integral(h->FindFixBin(minB), h->FindFixBin(maxB), "");
  //number of events in the 2.6-minB
  double nbgr=h->Integral(h->FindFixBin(low), h->FindFixBin(minB), "");

  //Nh, Nb, C, Nfb,Cfb as defined in Eq (27)-(31) of "J/Psi analysis steps" 2021 note
  //Nh is the number of events in the 2.6-minB range, we have to take into account bin size
  //Nh should be equal to nbgr
  auto Nh=back->Integral(low,minB)*100;
  //Nb is the number of events under the J/psi peak (3sigmas), we have to take into account bin size.
  auto Nb=back->Integral(minB,maxB)*100;
  //C the ratio of Nb/Nh
  auto C=Nb/Nh;
  //Nfb the total #background in the full spectrum
  auto Nfb=back->Integral(low,high);
  //Scaling factor
  auto Cfb=Nfb/Nh;

  //#of Jpsi events as defined in Eq (32) of note
  auto Njpsi=n3-(C*Nh);
  cout<<"Njpsi is:"<<Njpsi<<" and from fit:"<<f->GetParameter(0)<<endl;
  cout<<endl;


  return f;
}

TF1* g44(TH1* h, double x1, double x2, double x3, double x4, double x5, double low, double high){


  TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + TMath::Exp([3]+[4]*x)",
  low, high);

  TF1* back = new TF1(((TString)"background_fit")," TMath::Exp([0]+[1]*x)",
  low, high);

  TF1* gauss = new TF1(((TString)"peak"),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
  low, high);


  f->SetParameters(x1, x2, x3, x4, x5);

  f->SetParLimits (0,x1*0.1,x1*10);
  f->SetParLimits (1,3.08,3.12);
  f->SetParLimits (2,0.025,0.15);
  f->SetParLimits (3,0.0,1.0E5);
  f->SetParLimits (4,-1.0E5,0.0);

  f->SetParNames("X1", "Mean", "sigma", "X4","X5");

  h->Fit(f, "MEQ", "same", low, high);

  double par[5];
  f->GetParameters(par);
  back->SetParameters( par[3], par[4]);
  gauss->SetParameters(par[0], par[1], par[2]);


  back->SetLineColor(kBlue);
  gauss->SetLineColor(kGreen);

  back->Draw("same");
  gauss->Draw("same");


  return f;
}


double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
    double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
    double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
    double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return N*std::exp(- 0.5 * z * z);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return N*AA * std::pow(arg,n);
  }
}
//-------------------------------------------------------

TF1* CB2(TH1* h, double x1, double x2, double x3, double x4,  double x5, double x6, double x7, double x8, double low, double high){
  TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4]) + [5]*(x-[4])*(x-[4]) - [6]*(x-[4]) + [7]", low, high);

  TF1* CB = new TF1(h->GetName(),"[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4])",low,high);

  TF1* back = new TF1(((TString)"background_fit")," [1]*(x-[0])*(x-[0]) - [2]*(x-[0]) + [3]",
    low, high);



  f->SetParameters(x1, x2, x3, x4, x5,x6,x7,x8);
  f->SetParLimits (0,x1*0.1,x1*10);
  f->SetParLimits (1,0.10,1.5);
  f->SetParLimits (2,1,500);
  f->SetParLimits (3,0.025,0.15);
  f->SetParLimits (4,3.08,3.12);

  f->SetParLimits (6,0.0,1E5);
  //f->SetParLimits (5,0,50);
  f->SetParLimits (7,0.0,1E5);

  f->SetParNames("A","alpha", "n", "sigma", "mean","x6","x7","x8");


  h->Fit(f, "MEQ", "same",low, high);

  double par[8];
  f->GetParameters(par);
  back->SetParameters(par[4], par[5], par[6],par[7]);
  CB->SetParameters(par[0], par[1], par[2],par[3],par[4]);


  back->SetLineColor(kBlue);
  CB->SetLineColor(kGreen);

  back->Draw("same");
  CB->Draw("same");

  return f;
}

TF1* CB3(TH1* h, double x1, double x2, double x3, double x4,  double x5, double x6, double x7, double low, double high){


  TF1* CB = new TF1(h->GetName(),"[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4])",low,high);

  TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4]) + TMath::Exp([5]+[6]*x)", low, high);

  TF1* back = new TF1(((TString)"background_fit")," TMath::Exp([0]+[1]*x)",
    low, high);



  f->SetParameters(x1, x2, x3, x4, x5,x6,x7);
  f->SetParLimits (0,x1*0.1,x1*10);
  f->SetParLimits (1,0.10,1.5);
  f->SetParLimits (2,1,500);
  f->SetParLimits (3,0.025,0.15);
  f->SetParLimits (4,3.08,3.12);

  f->SetParLimits (5,0.0,1.0E5);
  f->SetParLimits (6,-1.0E5,0.0);

  f->SetParNames("A","alpha", "n", "sigma", "mean","x6","x7");


  h->Fit(f, "MEQ", "same",low, high);

  double par[7];
  f->GetParameters(par);
  back->SetParameters(par[5], par[6]);
  CB->SetParameters(par[0], par[1], par[2],par[3],par[4]);


  back->SetLineColor(kBlue);
  CB->SetLineColor(kGreen);

  back->Draw("same");
  CB->Draw("same");

  return f;
}


void plotgraph(double_t* variable, double_t* error, TString name, double_t OValue=0, double_t OValueE=0, int color=600){
  //double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  //double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  double Bins[4] ={8.625, 9.225, 9.73, 10.3};
  double BinsE[4]={0.425, 0.205, 0.27, 0.3};


  TGraphErrors *gr = new TGraphErrors(4,Bins,variable,BinsE,error);
  gr->SetMarkerStyle(8);
  gr->SetTitle(name+"; E_{#gamma} ; # Events");
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  gr->Draw("AEP");

  if(OValue!=0 && OValueE!=0){
    TLine *line = new TLine(8.2,OValue,10.6,OValue);
    TBox *box = new TBox(8.2,OValue-OValueE,10.6,OValue+OValueE);
    line->SetLineColor(color);
    box->SetFillColorAlpha(color, 0.15);
    line->Draw("same");
    box->Draw("same");
  }

}

void plot2graph(double_t variable[4][10], double_t error[4][10], TString name, int g1=0, int g2=1, double_t OValue=0, double_t OValueE=0, double_t OValue2=0, double_t OValueE2=0){
  //double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  //double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  double Bins[4] ={8.625, 9.225, 9.73, 10.3};
  double BinsE[4]={0.425, 0.205, 0.27, 0.3};

  TGraphErrors *gr = new TGraphErrors(4,Bins,variable[g1],BinsE,error[g1]);
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(kBlue);
  gr->SetTitle(name+"; E_{#gamma} ; "+ name);
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  if(name=="Mean")
    gr->GetYaxis()->SetRangeUser(3.06,3.12);
  else if(name=="Sigma")
    gr->GetYaxis()->SetRangeUser(0.0, 0.10);
  
  gr->Draw("AEP");

  TGraphErrors *gr2 = new TGraphErrors(4,Bins,variable[g2],BinsE,error[g2]);
  gr2->SetMarkerStyle(23);
  gr2->SetMarkerColor(kRed);
  gr2->Draw("P");



  auto legend = new TLegend(0.25,0.10);
  legend->AddEntry(gr, lnames[g1].c_str(), "p");
  legend->AddEntry(gr2,lnames[g2].c_str(), "p");
  legend->Draw("same");

  if(OValue!=0 && OValueE!=0){
    TLine *line = new TLine(8.2,OValue,10.6,OValue);
    TBox *box = new TBox(8.2,OValue-OValueE,10.6,OValue+OValueE);
    line->SetLineColor(kBlue);
    box->SetFillColorAlpha(kBlue, 0.15);
    line->Draw("same");
    box->Draw("same");
  }

  if(OValue2!=0 && OValueE2!=0){
    TLine *line2 = new TLine(8.2,OValue2,10.6,OValue2);
    TBox *box2 = new TBox(8.2,OValue2-OValueE2,10.6,OValue2+OValueE2);
    line2->SetLineColor(kRed);
    box2->SetFillColorAlpha(kRed, 0.15);
    line2->Draw("same");
    box2->Draw("same");
  }

}






#endif
