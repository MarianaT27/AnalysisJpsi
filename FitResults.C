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
#include "bib/func.h"


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

    Get_histos(h_Invariant, h_invariant_Egamma, "InvariantMassS19", EP);
    Get_histos(h_Invariant, h_invariant_Egamma, "InvariantMassF18in", EP);
    Get_histos(h_Invariant, h_invariant_Egamma, "InvariantMassF18out", EP);



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

     TPaveText *pt = new TPaveText(.05,.1,.95,.8);

    double sum=0;
    int init_amp_fit;
    //-----------------------------------------
    //-----------------------------------------
    // GAUSS + POLINOMIAL
    //-----------------------------------------
    //-----------------------------------------


    can->Divide(3,2);
    cout<<"Getting Gauss+ Pol(2) fits..."<<endl;

    for(int i=0;i<4;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        init_amp_fit = (h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) > 0.0) ? h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) : 5;
        TF1 *g=g43(h_invariant_Egamma[i],init_amp_fit,3.096,0.04,7,3,7,min,max);

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
    can->cd(5);
    pt->AddText("Gauss + Pol(2) ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();

    can->cd(6);
    h_Invariant->SetTitle("Full Set");
    h_Invariant->Draw();
    init_amp_fit = (h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) > 0.0) ? h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) : 5;
    TF1 *gp=g43(h_Invariant,init_amp_fit,3.096,0.04,80,70,100,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    //-----------------------------------------
    //-----------------------------------------
    // CRYSTAL BALL + POLINOMIAL
    //-----------------------------------------
    //-----------------------------------------

    sum=0;
    can->Divide(3,2);
    cout<<"Getting Crystal Ball + Pol(2) fits..."<<endl;
    for(int i=0;i<4;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        //if(i<2)
        //  min=3.2;
        init_amp_fit = (h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) > 0.0) ? h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) : 5;
        TF1 *c=CB2(h_invariant_Egamma[i],init_amp_fit,0.75,150,0.04,3.096,7,2,7,min,max);//e-e+ untagged
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
    can->cd(5);
    pt->Clear();
    pt->AddText("Crystal Ball + Pol(2) ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();

    can->cd(6);
    h_Invariant->SetTitle("Full Set");
    h_Invariant->Draw();
    init_amp_fit = (h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) > 0.0) ? h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) : 5;
    TF1 *cp3=CB2(h_Invariant,init_amp_fit,0.75,150,0.04,3.096,70,2,70,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    //-----------------------------------------
    //-----------------------------------------
    // GAUSS + EXPONENTIAL
    //-----------------------------------------
    //-----------------------------------------
    sum=0;
    can->Divide(3,2);
    cout<<"Getting gauss + Exp fits..."<<endl;
    for(int i=0;i<4;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        init_amp_fit = (h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) > 0.0) ? h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) : 5;
        TF1 *g=g44(h_invariant_Egamma[i],init_amp_fit,3.096,0.04,7,-2,min,max);
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
    can->cd(5);
    pt->Clear();
    pt->AddText("Gauss + Exp ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();

    can->cd(6);
    h_Invariant->SetTitle("Full Set");
    h_Invariant->Draw();
    init_amp_fit = (h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) > 0.0) ? h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) : 5;
    TF1 *ge=g44(h_Invariant,init_amp_fit,3.096,0.04,7,-2,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    //-----------------------------------------
    //-----------------------------------------
    // CRYSTAL BALL + EXPONENTIAL
    //-----------------------------------------
    //-----------------------------------------
    can->Divide(3,2);
    sum=0;
    cout<<"Getting CB + Exp fits..."<<endl;
    for(int i=0;i<4;i++){
        can->cd(i+1);
        h_invariant_Egamma[i]->Draw();
        init_amp_fit = (h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) > 0.0) ? h_invariant_Egamma[i]->GetBinContent(h_invariant_Egamma[i]->FindBin(3.096)) : 5;
        TF1 *c=CB3(h_invariant_Egamma[i],init_amp_fit,0.75,150,0.04,3.096,7,-2,min,max);//e-e+ untagged
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
    can->cd(5);
    pt->Clear();
    pt->AddText("Crystal Ball + Exp ");
    pt->AddText(Form("%f",sum));
    pt->SetTextSize(0.1);
    pt->Draw();

    can->cd(6);
    h_Invariant->SetTitle("Full Set");
    h_Invariant->Draw();
    init_amp_fit = (h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) > 0.0) ? h_Invariant->GetBinContent(h_Invariant->FindBin(3.096)) : 5;
    TF1 *ce3=CB3(h_Invariant,init_amp_fit,0.75,150,0.04,3.096,7,-2,min,max);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();


    //-----------------------------------------
    //-----------------------------------------
    // RESULTS
    //-----------------------------------------
    //-----------------------------------------
    can->Divide(2,2);
    can->cd(1);
    plotgraph(Events,EventsE, "Number Events polinomial");

    can->cd(2);
    plot2graph(Amplitude,AmplitudeE, "N_{J/#psi}",0,1);

    can->cd(3);
    plot2graph(Mean, MeanE, "Mean",0,1, gp->GetParameter(1),gp->GetParError(1),cp3->GetParameter(4),cp3->GetParError(4));
    

    can->cd(4);
    plot2graph(Sigma, SigmaE, "Sigma",0,1, gp->GetParameter(2),gp->GetParError(2), cp3->GetParameter(3),cp3->GetParError(3));

    can->Print( (pdf_original+"(").c_str());


    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    plotgraph(Events,EventsE, "Number Events exponential");

    can->cd(2);
    plot2graph(Amplitude,AmplitudeE, "N_{J/#psi}",2,3);

    can->cd(3);
    plot2graph(Mean, MeanE, "Mean",2,3,ge->GetParameter(1),ge->GetParError(1),ce3->GetParameter(4),ce3->GetParError(4));

    can->cd(4);
    plot2graph(Sigma, SigmaE, "Sigma",2,3,ge->GetParameter(2),ge->GetParError(2), ce3->GetParameter(3),ce3->GetParError(3));


    can->Print( (pdf_original+")").c_str());




    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;
    

  }//




