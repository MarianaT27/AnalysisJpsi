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


const string fitnames[4]={"Gauss+Pol","Crystal_Ball+Pol","Gauss+Exp","Crystal_Ball+Exp"};


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
    //gStyle->SetStatW(0.15);
    //gStyle->SetStatH(0.15);

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

    double sum;
    int init_amp_fit;
    //-----------------------------------------
    //-----------------------------------------
    // Fitting
    //-----------------------------------------
    //-----------------------------------------
    TF1* rFit_F[4];
    for(int fit=0;fit<4;fit++){
      can->Clear();
      sum=0;
      can->Divide(3,2);
      Fit_Function Fit_func;
    
      TF1* rFit;

      for(int i=0;i<4;i++){
          can->cd(i+1);
          if((fit==0||fit==1)&&i==0){
            min=2.6;
            max=3.2;
          }
          else if((fit==0||fit==1)&&i==1){
            min=2.6;
            max=3.2;
          }
          else{
            min=2.6;
            max=3.4;
          }
          Fit_func.Set_Limits(min,max);
          h_invariant_Egamma[i]->Draw();
          switch (fit){
            case 0:
              rFit=Fit_func.Fit_Gauss_Pol2(h_invariant_Egamma[i]);
              break;
            case 1:
              rFit=Fit_func.Fit_CrystallBall_Pol2(h_invariant_Egamma[i]);
              break;
            case 2:
              rFit=Fit_func.Fit_Gauss_Exp(h_invariant_Egamma[i]);
              break;
            case 3:
              rFit=Fit_func.Fit_CrystallBall_Exp(h_invariant_Egamma[i]);
              break;
            default:
              break;
          }
          
          sum=sum+rFit->GetParameter(0);
          //Get parameters
          Events[i]=h_invariant_Egamma[i]->GetEntries();
          Amplitude[fit][i]=rFit->GetParameter(0);
          Mean[fit][i]=rFit->GetParameter(1);
          Sigma[fit][i]=rFit->GetParameter(2);
          //Get error
          AmplitudeE[fit][i]=rFit->GetParError(0);
          MeanE[fit][i]=rFit->GetParError(1);
          SigmaE[fit][i]=rFit->GetParError(2);
      }
      can->cd(5);
      pt->AddText((fitnames[fit]).c_str());
      pt->AddText(Form("%f",sum));
      pt->SetTextSize(0.1);
      pt->Draw();

      can->cd(6);
      h_Invariant->SetTitle("Full Set");
      h_Invariant->Draw();

      switch (fit){
        case 0:
          rFit_F[fit]=Fit_func.Fit_Gauss_Pol2(h_Invariant);
          break;
        case 1:
          rFit_F[fit]=Fit_func.Fit_CrystallBall_Pol2(h_Invariant);
          break;
        case 2:
          rFit_F[fit]=Fit_func.Fit_Gauss_Exp(h_Invariant);
          break;
        case 3:
          rFit_F[fit]=Fit_func.Fit_CrystallBall_Exp(h_Invariant);
          break;
        default:
          break;
      }
      can->Print( (pdf_original+"(").c_str());
      can->Clear();
      pt->Clear();

    }



    //-----------------------------------------
    //-----------------------------------------
    // RESULTS
    //-----------------------------------------
    //-----------------------------------------
   can->Divide(2,2);
    can->cd(1);
    plotgraph(Events,EventsE, "Number Events polinomial");

    can->cd(2);
    plot2graph(Amplitude,AmplitudeE, "N_{J/#psi}");

    can->cd(3);
    plot2graph(Mean, MeanE, "Mean",0,1,rFit_F[0]->GetParameter(1),rFit_F[0]->GetParError(1),rFit_F[1]->GetParameter(1),rFit_F[1]->GetParError(1));
    

    can->cd(4);
    plot2graph(Sigma, SigmaE, "Sigma",0,1,rFit_F[0]->GetParameter(2),rFit_F[0]->GetParError(2),rFit_F[1]->GetParameter(2),rFit_F[1]->GetParError(2));

    can->Print( (pdf_original+"(").c_str());


    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    plotgraph(Events,EventsE, "Number Events exponential");

    can->cd(2);
    plot2graph(Amplitude,AmplitudeE, "N_{J/#psi}",2,3);

    can->cd(3);
    plot2graph(Mean, MeanE, "Mean",2,3,rFit_F[2]->GetParameter(1),rFit_F[2]->GetParError(1),rFit_F[3]->GetParameter(1),rFit_F[3]->GetParError(1));

    can->cd(4);
    plot2graph(Sigma, SigmaE, "Sigma",2,3,rFit_F[2]->GetParameter(2),rFit_F[2]->GetParError(2),rFit_F[3]->GetParameter(2),rFit_F[3]->GetParError(2));


    can->Print( (pdf_original+")").c_str());


    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;
    

  }//




