#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "reader.h"
#include "clas12reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "bib/Fitfunc.h"


// QADB header and namespace
#include "QADB.h"
using namespace QA;


struct cartesian {
  double x;
  double y;
  double z;
};


struct response {
  cartesian pos;
  cartesian tpos;
  double time;
  double energy;
  double path;
  int sector;
  int layer;
  int index;
  int component;
  double du;
  double dv;
  double dw;
  double m2u;
  double m2v;
  double m2w;
  double m3u;
  double m3v;
  double m3w;
  double u;
  double v;
  double w;
  double widthu;
  double widthv;
  double widthw;
  double x;
  double y;
  double z;
  double quality;
  int degree;
};


struct particl {
  TLorentzVector lorentz;
  cartesian vertexinfo;
  map<int,response> responses;
  int index = -1;
  double beta;
  double chi2pid;
  int status;
  int pid;
  double vtime;
  double E;
};


/*

TF1* g43(TH1* h, double x1, double x2, double x3, double x4, double x5, double x6, double low, double high){

  TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + [3]*(x-[1])*(x-[1]) - [4]*(x-[1]) + [5]",
  low, high);

  TF1* back = new TF1(((TString)"background_fit")," [1]*(x-[0])*(x-[0]) - [2]*(x-[0]) + [3]",low, high);

  TF1* gauss = new TF1(((TString)"peak"),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
  low, high);


  f->SetParameters(x1, x2, x3, x4,x5,x6);

  f->SetParLimits (0,1,x1);
  f->SetParLimits (1,3.08,3.11);
  f->SetParLimits (2,0.030,0.1);
  //f->SetParLimits (3,0,50);
  f->SetParLimits (5,0.5,x6);

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

  f->SetParLimits (0,1,x1);
  f->SetParLimits (1,3.08,3.11);
  f->SetParLimits (2,0.030,0.09);
  //f->SetParLimits (3,0,50);
  //f->SetParLimits (4,-3,-0.01);

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
  f->SetParLimits (0,1,x1);
  //f->SetParLimits (1,0,1);
  //f->SetParLimits (2,-x3,x3);
  f->SetParLimits (3,0.02,0.09);
  f->SetParLimits (4,3.08,3.11);

  //f->SetParLimits (5,0,50);
  f->SetParLimits (7,0.5,x8);

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
  f->SetParLimits (0,1,x1);
  f->SetParLimits (1,0.50,20);
  //f->SetParLimits (2,-x3,x3);
  f->SetParLimits (3,0.03,0.09);
  f->SetParLimits (4,3.08,3.11);
  //f->SetParLimits (6,-3.025,-0.01);

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
  double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  TGraphErrors *gr = new TGraphErrors(10,Bins,variable,BinsE,error);
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
  double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  TGraphErrors *gr = new TGraphErrors(10,Bins,variable[g1],BinsE,error[g1]);
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(kBlue);
  gr->SetTitle(name+"; E_{#gamma} ; # Events");
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  gr->Draw("AEP");

  TGraphErrors *gr2 = new TGraphErrors(10,Bins,variable[g2],BinsE,error[g2]);
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


void Bin_invariant(TH1F* h[], double_t Epho, double_t invariantMass){
  if(Epho>=8.2 && Epho<8.65)
    h[0]->Fill(invariantMass);
  else if(Epho<8.9)
    h[1]->Fill(invariantMass);
  else if(Epho<9.05)
    h[2]->Fill(invariantMass);
  else if(Epho<9.2)
    h[3]->Fill(invariantMass);
  else if(Epho<9.46)
    h[4]->Fill(invariantMass);
  else if(Epho<9.7)
    h[5]->Fill(invariantMass);
  else if(Epho<10)
    h[6]->Fill(invariantMass);
  else if(Epho<10.2)
    h[7]->Fill(invariantMass);
  else if(Epho<10.4)
    h[8]->Fill(invariantMass);
  else if(Epho<10.6)
    h[9]->Fill(invariantMass);

}*/


int FitResults(string pdfname="Full_Invariant", double EP=-0.030 , double Q2cut=0.5, double MMcut=0.4) {
// Record start time
  auto start = std::chrono::high_resolution_clock::now();

    //RESULTS RESULTS
    TH1F* h_Invariant= new TH1F("h_im",";M(e+e-),GeV ",100,2.5,3.5);

    TH1F* h_invariant_Egamma[10];
    for(int j = 0; j < 10; j++){
      ostringstream name;
      name << "Egamma_bin" << j+1;
      ostringstream sstr;
      sstr<<"Egamma_bin" << j+1<<";M(e+e-),GeV ";
      h_invariant_Egamma[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),100,2.5,3.5);;
      name.str("");
      sstr.str("");
      
    }

  
    //---------------SPRING 2019----------

    TFile *file1 = new TFile("/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/InvariantMassS19.root","READ");
    TTree *results1 = (TTree*)file1->Get("results");
    Double_t invariantMass,t,Q2, EPcut, MM, Egamma;

    results1->SetMakeClass(1);
    results1->SetBranchAddress("invariantMass",&invariantMass);
    results1->SetBranchAddress("t",&t);
    results1->SetBranchAddress("Q2",&Q2);
    results1->SetBranchAddress("EPcut",&EPcut);
    results1->SetBranchAddress("MM",&MM);
    results1->SetBranchAddress("Egamma",&Egamma);

    for(int fc=0;fc<results1->GetEntries();fc++) {//Run 5032 to 5419 // 6616 6783
      results1->GetEntry(fc);
      if(EPcut<EP)
        continue;

      //if(Q2>Q2cut)
        //continue;

      if(MM>MMcut)
        continue;
      h_Invariant->Fill(invariantMass);

      Bin_invariant(h_invariant_Egamma,Egamma,invariantMass);


    }

    //---------------FALL 2018 INBENDING----------
    TFile *file2 = new TFile("/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/InvariantMassF18in.root","READ");
    TTree *results2 = (TTree*)file2->Get("results");

    results2->SetMakeClass(1);
    results2->SetBranchAddress("invariantMass",&invariantMass);
    results2->SetBranchAddress("t",&t);
    results2->SetBranchAddress("Q2",&Q2);
    results2->SetBranchAddress("EPcut",&EPcut);
    results2->SetBranchAddress("MM",&MM);
    results2->SetBranchAddress("Egamma",&Egamma);

    for(int fc=0;fc<results2->GetEntries();fc++) {//Run 5032 to 5419 // 6616 6783
      results2->GetEntry(fc);
      if(EPcut<EP)
        continue;

      //if(Q2>Q2cut)
        //continue;

      if(MM>MMcut)
        continue;
      h_Invariant->Fill(invariantMass);

      Bin_invariant(h_invariant_Egamma,Egamma,invariantMass);
    }

    //---------------FALL 2018  OUTBENDING----------
    TFile *file3 = new TFile("/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/InvariantMassF18out.root","READ");
    TTree *results3 = (TTree*)file3->Get("results");

    results3->SetMakeClass(1);
    results3->SetBranchAddress("invariantMass",&invariantMass);
    results3->SetBranchAddress("t",&t);
    results3->SetBranchAddress("Q2",&Q2);
    results3->SetBranchAddress("EPcut",&EPcut);
    results3->SetBranchAddress("MM",&MM);
    results3->SetBranchAddress("Egamma",&Egamma);

    for(int fc=0;fc<results3->GetEntries();fc++) {//Run 5032 to 5419 // 6616 6783
      results3->GetEntry(fc);
      if(EPcut<EP)
        continue;

      //if(Q2>Q2cut)
        //continue;

      if(MM>MMcut)
        continue;
      h_Invariant->Fill(invariantMass);

      Bin_invariant(h_invariant_Egamma,Egamma,invariantMass);

    }



    TCanvas *can = new TCanvas("can","canvas",200,10,1000,700);
    string pdf_original=pdfname+".pdf";



    gStyle->SetOptFit(1111);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);

    double min,max;
    min=2.6;
    max=3.4;

    double Events[10]={0};
    double Amplitude[4][10]={0};
    double Mean[4][10]={0};
    double Sigma[4][10]={0};

    double EventsE[10]={0};
    double AmplitudeE[4][10]={0};
    double MeanE[4][10]={0};
    double SigmaE[4][10]={0};

    double sum=0;

    //-----------------------------------------
    //-----------------------------------------
    cout<<"Fitting full data sets..."<<endl;

    h_Invariant->SetTitle("Gauss + Pol(2)");
    h_Invariant->Draw();
    TF1 *gp=g43(h_Invariant,2E3,3.09,0.08,80,70,100,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();
    
    h_Invariant->SetTitle("Gauss + Exp");
    h_Invariant->Draw();
    TF1 *ge=g44(h_Invariant,2E3,3.09,0.08,20,-10,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    h_Invariant->SetTitle("CristalBall + Pol(2)");
    h_Invariant->Draw();
    TF1 *cp3=CB2(h_Invariant,2E3,0.75,1E10,0.08,3.1,280,270,100,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    h_Invariant->SetTitle("CristalBall + Exp");
    h_Invariant->Draw();
    TF1 *ce3=CB3(h_Invariant,2E3,20,1,0.08,3.1,15,-5,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    //-----------------------------------------
    //-----------------------------------------


    can->Divide(4,3);
    cout<<"Getting Gauss+ Pol(2) fits..."<<endl;

    for(int i=0;i<10;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();

        //if(i<2)
         // min=3.2;
        
        TF1 *g=g43(h_invariant_Egamma[i],250,3.09,0.08,80,70,100,min,max);

        sum=sum+g->GetParameter(0);
        //Get parameters
        Events[i]=h_invariant_Egamma[i]->GetEntries();
        Amplitude[0][i]=g->GetParameter(0);
        Mean[0][i]=g->GetParameter(1);
        Sigma[0][i]=g->GetParameter(2);
        //Get error
        AmplitudeE[0][i]=g->GetParError(0);
        MeanE[0][i]=g->GetParError(1);
        SigmaE[0][i]=g->GetParError(2);
    }
    can->cd(11);

    TPaveText *pt = new TPaveText(.05,.1,.95,.8);
    pt->AddText("Gauss + Pol(2) ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();


    can->Print( (pdf_original+"(").c_str());
    can->Clear();
    sum=0;
    
    can->Divide(4,3);
    cout<<"Getting Cristal Ball + Pol(2) fits..."<<endl;

    for(int i=0;i<10;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        //if(i<2)
        //  min=3.2;
        TF1 *c=CB2(h_invariant_Egamma[i],250,0.75,1E10,0.08,3.1,280,270,100,min,max);//e-e+ untagged
        //Get parameters
        sum=sum+c->GetParameter(0);
        Amplitude[1][i]=c->GetParameter(0);
        Mean[1][i]=c->GetParameter(4);
        Sigma[1][i]=c->GetParameter(3);
        //Get error
        AmplitudeE[1][i]=c->GetParError(0);
        MeanE[1][i]=c->GetParError(4);
        SigmaE[1][i]=c->GetParError(3);
    }
    can->cd(11);
    pt->Clear();
    pt->AddText("Crystal Ball + Pol(2) ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();


    can->Print( (pdf_original+"(").c_str());
    can->Clear();
    sum=0;

    can->Divide(4,3);
    cout<<"Getting gauss + Exp fits..."<<endl;
    for(int i=0;i<10;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        TF1 *g=g44(h_invariant_Egamma[i],250,3.09,0.08,15,-5,min,max);
        sum=sum+g->GetParameter(0);
        //Get parameters
        Amplitude[2][i]=g->GetParameter(0);
        Mean[2][i]=g->GetParameter(1);
        Sigma[2][i]=g->GetParameter(2);
        //Get error
        AmplitudeE[2][i]=g->GetParError(0);
        MeanE[2][i]=g->GetParError(1);
        SigmaE[2][i]=g->GetParError(2);
    }
    can->cd(11);
    pt->Clear();
    pt->AddText("Gauss + Exp ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();
    can->Print( (pdf_original+"(").c_str());
    can->Clear();
    
    can->Divide(4,3);
    sum=0;
    cout<<"Getting CB + Exp fits..."<<endl;
    for(int i=0;i<10;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        TF1 *c=CB3(h_invariant_Egamma[i],250,0.75,1E10,0.08,3.1,15,-5,min,max);//e-e+ untagged
        sum=sum+c->GetParameter(0);
        //Get parameters
        Amplitude[3][i]=c->GetParameter(0);
        Mean[3][i]=c->GetParameter(4);
        Sigma[3][i]=c->GetParameter(3);
        //Get error
        AmplitudeE[3][i]=c->GetParError(0);
        MeanE[3][i]=c->GetParError(4);
        SigmaE[3][i]=c->GetParError(3);
    }
    can->cd(11);
    pt->Clear();
    pt->AddText("Crystal Ball + Exp ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();

    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    can->Divide(2,2);
    can->cd(1);
    plotgraph(Events,EventsE, "Number Events polinomial");

    can->cd(2);
    plot2graph(Amplitude,AmplitudeE, "N_{J/#psi}",0,1);

    can->cd(3);
    plot2graph(Mean, MeanE, "Mean ",0,1, gp->GetParameter(1),gp->GetParError(1),cp3->GetParameter(4),cp3->GetParError(4));
    

    can->cd(4);
    plot2graph(Sigma, SigmaE, "Sigma ",0,1, gp->GetParameter(2),gp->GetParError(2), cp3->GetParameter(3),cp3->GetParError(3));

    can->Print( (pdf_original+"(").c_str());


    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    plotgraph(Events,EventsE, "Number Events exponential");

    can->cd(2);
    plot2graph(Amplitude,AmplitudeE, "N_{J/#psi}",2,3);

    can->cd(3);
    plot2graph(Mean, MeanE, "Mean ",2,3,ge->GetParameter(1),ge->GetParError(1),ce3->GetParameter(4),ce3->GetParError(4));

    can->cd(4);
    plot2graph(Sigma, SigmaE, "Sigma ",2,3,ge->GetParameter(2),ge->GetParError(2), ce3->GetParameter(3),ce3->GetParError(3));


    can->Print( (pdf_original+")").c_str());




    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;
    

  }//




