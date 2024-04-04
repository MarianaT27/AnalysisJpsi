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

    for(int i=0;i<4;i++){
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

    for(int i=0;i<4;i++){
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
    for(int i=0;i<4;i++){
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
    for(int i=0;i<4;i++){
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




