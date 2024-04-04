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


//x1=100---[0]
//x2=3.1---[1]
//x3=0.05--[2]
//x4=15----[3]
//x5=24----[4]
//x6=4.4---[5]


TF1* g43(TH1* h, double x1, double x2, double x3, double x4, double x5, double x6, double low, double high){

  TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + [3]*(x-[1])*(x-[1]) - [4]*(x-[1]) + [5]",
  low, high);

  TF1* back = new TF1(((TString)"background_fit")," [1]*(x-[0])*(x-[0]) - [2]*(x-[0]) + [3]",low, high);

  TF1* gauss = new TF1(((TString)"peak"),"[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
  low, high);


  f->SetParameters(x1, x2, x3, x4,x5,x6);

  f->SetParLimits (0,1,1E3);
  f->SetParLimits (1,3.09,3.11);
  f->SetParLimits (2,-0.5,0.5);
  //f->SetParLimits (3,0,50);
  //f->SetParLimits (4,-5,-2);

  f->SetParNames("X1", "Mean", "sigma", "X4","X5","X6");

  h->Fit(f, "ME", "same", low, high);

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

  f->SetParLimits (0,1,5E2);
  f->SetParLimits (1,3.09,3.11);
  f->SetParLimits (2,-0.5,0.5);
  //f->SetParLimits (3,0,50);
  f->SetParLimits (4,-3,-0.01);

  f->SetParNames("X1", "Mean", "sigma", "X4","X5");

  h->Fit(f, "ME", "same", low, high);

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
  f->SetParLimits (0,1,1E3);
  f->SetParLimits (1,0.5,0.6);
  f->SetParLimits (2,1,1E6);
  f->SetParLimits (3,-0.5,0.5);
  f->SetParLimits (4,3.09,3.11);

  //f->SetParLimits (5,0,50);
  //f->SetParLimits (7,0,700);

  f->SetParNames("A","alpha", "n", "sigma", "mean","x6","x7","x8");


  h->Fit(f, "ME", "same",low, high);

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
  f->SetParLimits (0,1,1E3);
  f->SetParLimits (1,0.5,0.6);
  f->SetParLimits (2,1,1E6);
  f->SetParLimits (3,-0.5,0.5);
  f->SetParLimits (4,3.09,3.11);
  f->SetParLimits (6,-3.025,-0.01);

  f->SetParNames("A","alpha", "n", "sigma", "mean","x6","x7");


  h->Fit(f, "ME", "same",low, high);

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


int untaggedfromroot(string nameFile="PASS1_FALL", int version=-19, double Beam_E=10.2, int PASS=1,string outFile="out", int MC=0, bool QA=true) {
// Record start time
  auto start = std::chrono::high_resolution_clock::now();
  int max=0;


  //********************************
  //DECLARATION OF VARIABLES
  //********************************
  //MISSING MOMENTUM AND INVARIANT MASS
  TLorentzVector miss;
  TLorentzVector invariant_ee;
  TLorentzVector invariant_eeFT;
  TLorentzVector electronFT_vec;  //FT e-
  Double_t Q2, W;

  Double_t electron_E_recon;
  Double_t positron_E_recon;

  Double_t electron_px_recon, electron_py_recon, electron_pz_recon;
  Double_t positron_px_recon, positron_py_recon, positron_pz_recon;

  Double_t px,py,pz,E;
  Double_t p;

  //Testing
  TLorentzVector transverse;
  Double_t Pt;


  //CONFIG VARIABLES
  //Double_t Beam_E;
  //FALL 10.6
  //SPRING 10.2
  TLorentzVector beam(0,0,Beam_E,Beam_E);
  TLorentzVector target(0,0,0,0.938);
  Int_t run_number;
  Int_t event_number;

    //Electron variables
  Double_t electron_p, electron_theta, electron_phi, electron_px, electron_py, electron_pz;
  Double_t electron_vx, electron_vy, electron_vz, electron_vt;
  Double_t electron_pcal_v, electron_pcal_w;
  Double_t electron_pcal_energy, electron_ecin_energy, electron_ecout_energy; //Electron Energies
  Double_t electron_energy;
  Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
  Double_t electron_sfpcal, electron_sfecin, electron_sfecout;


  //Positron variables
  Double_t positron_p, positron_theta, positron_phi, positron_px, positron_py, positron_pz;
  Double_t positron_vx, positron_vy, positron_vz, positron_vt;
  Double_t positron_pcal_v, positron_pcal_w;
    Double_t positron_pcal_energy, positron_ecin_energy, positron_ecout_energy; //Positron Energies
    Double_t positron_energy;
    Double_t positron_chi2pid;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;

    //Proton varianbles
    Double_t proton_p, proton_theta, proton_phi, proton_px, proton_py, proton_pz, proton_status, proton_beta, proton_chi2pid;
    Double_t proton_vx, proton_vy, proton_vz, proton_vt, proton_pid;
    Double_t proton_energy;



    //Photon variables
    Double_t photon_p, photon_theta, photon_phi, photon_px, photon_py, photon_pz;
    Double_t photon_vx, photon_vy, photon_vz, photon_vt,photon_energy;
    int electron_photon;
    int positron_photon;
    Double_t electron_photonE;
    Double_t positron_photonE;
    Double_t E_photon=0;

    //Forward Tagger Electron variables
    Double_t FT_E;
    Double_t FT_Px;
    Double_t FT_Py;
    Double_t FT_Pz;
    Double_t FT_E_new;
    Double_t FT_Px_new;
    Double_t FT_Py_new;
    Double_t FT_Pz_new;
    Double_t FT_vt;
    Double_t FT_x,FT_y,FT_z, R;
    Double_t FT_vtC;
    //Testing
    Double_t FT_vtFcal,FT_xFcal,FT_yFcal,FT_zFcal,RFcal;

    int N_FT;
    int number_of_electrons;
    int number_of_positrons;
    int number_protons;
    int JPSI=0;
    bool electronAccepted, positronAccepted, protonAccepted;

    Double_t score_pos,score_ele;


    //KINEMATICS
    TH1F* h_E_photon = new TH1F("h_E_photon","Photon Energy",120,0,0);
    TH1F* h_Q2 = new TH1F("h_Q2","Q^2;Q^2;Counts",120,0,0);
    TH1F* h_W = new TH1F("h_W","Hadronic Mass;W;Counts",120,4,4.6);
    TH1F* h_MM= new TH1F("h_MM","Missing Mass; MM; Counts ",90,0,2);


    TH2F* h_mm_vs_invariantmass_t= new TH2F("h_mm_vs_invariantmass_t","MM vs M(e^+e^-) for t+-2;M(e^+e^-), GeV;Missing Mass, GeV",150,1,3.5,150,-1,1);
    TH2F* h_qq_vs_invariantmass= new TH2F("h_qq_vs_invariantmass","Q2 vs M(e^+e^-) outside time cut;M(e^+e^-), GeV;Q2",150,2,3.5,150,0,0.5);

    //RESULTS RESULTS
    TH1F* h_Invariant= new TH1F("h_im",";M(e+e-),GeV ",100,2.5,3.5);
    TH1F* h_Invariant_not= new TH1F("h_Invariant_not","Outside time range;M(e+e-),GeV",90,2.6,3.5);


    //MOMENTUM FT ELECTRON
    TH2F* h_sum_cut= new TH2F("h_sum_cut",";Invariant Mass ; #Sum E - #Sum p_{z} -m ",100,2.5,3.5,100,-0.5,0.5);

    
    TString root_file;
    if(MC==1)
      root_file= "/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/"+nameFile+".root";
    else
      root_file = "/w/hallb-scshelf2102/clas12/mtenorio/Analysis/"+nameFile+".root";


    TFile *file = new TFile(root_file,"READ");
    TTree *tree = (TTree*)file->Get("analysis");

    TString out_file="/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/InvariantMass"+nameFile+".root";

    TFile *file2 = new TFile(out_file,"RECREATE");
    TTree *results = new TTree("results",out_file);

    Double_t invariantMass, EPcut,MM, Egamma, t;
    results->Branch("invariantMass",&invariantMass,"invariantMass/d");
    results->Branch("t",&t,"t/d");
    results->Branch("Q2",&Q2,"Q2/d");
    results->Branch("EPcut",&EPcut,"EPcut/d");
    results->Branch("MM",&MM,"MM/d");
    results->Branch("Egamma",&Egamma,"Egamma/d");

    tree->SetMakeClass(1);

    tree->SetBranchAddress("Beam_E",&Beam_E);

    tree->SetBranchAddress("run_number",&run_number);
    tree->SetBranchAddress("event_number",&event_number);

    tree->SetBranchAddress("electron_p",&electron_p);
    tree->SetBranchAddress("electron_theta",&electron_theta);
    tree->SetBranchAddress("electron_phi",&electron_phi);
    tree->SetBranchAddress("electron_px",&electron_px);
    tree->SetBranchAddress("electron_py",&electron_py);
    tree->SetBranchAddress("electron_pz",&electron_pz);
    tree->SetBranchAddress("electron_vx",&electron_vx);
    tree->SetBranchAddress("electron_vy",&electron_vy);
    tree->SetBranchAddress("electron_vz",&electron_vz);
    tree->SetBranchAddress("electron_vt",&electron_vt);
    tree->SetBranchAddress("electron_pcal_v",&electron_pcal_v);
    tree->SetBranchAddress("electron_pcal_w",&electron_pcal_w);
    tree->SetBranchAddress("electron_pcal_energy",&electron_pcal_energy);
    tree->SetBranchAddress("electron_ecin_energy",&electron_ecin_energy);
    tree->SetBranchAddress("electron_ecout_energy",&electron_ecout_energy);
    tree->SetBranchAddress("electron_energy",&electron_energy);
    tree->SetBranchAddress("electron_m2ecin",&electron_m2ecin);
    tree->SetBranchAddress("electron_m2ecout",&electron_m2ecout);
    tree->SetBranchAddress("electron_m2pcal",&electron_m2pcal);
    tree->SetBranchAddress("electron_sfecin",&electron_sfecin);
    tree->SetBranchAddress("electron_sfpcal",&electron_sfpcal);
    tree->SetBranchAddress("electron_sfecout",&electron_sfecout);

    tree->SetBranchAddress("score_ele",&score_ele);
    tree->SetBranchAddress("score_pos",&score_pos);

    
    tree->SetBranchAddress("positron_p",&positron_p);
    tree->SetBranchAddress("positron_theta",&positron_theta);
    tree->SetBranchAddress("positron_phi",&positron_phi);
    tree->SetBranchAddress("positron_px",&positron_px);
    tree->SetBranchAddress("positron_py",&positron_py);
    tree->SetBranchAddress("positron_pz",&positron_pz);
    tree->SetBranchAddress("positron_vx",&positron_vx);
    tree->SetBranchAddress("positron_vy",&positron_vy);
    tree->SetBranchAddress("positron_vz",&positron_vz);
    tree->SetBranchAddress("positron_vt",&positron_vt);
    tree->SetBranchAddress("positron_pcal_v",&positron_pcal_v);
    tree->SetBranchAddress("positron_pcal_w",&positron_pcal_w);
    tree->SetBranchAddress("positron_pcal_energy",&positron_pcal_energy);
    tree->SetBranchAddress("positron_ecin_energy",&positron_ecin_energy);
    tree->SetBranchAddress("positron_ecout_energy",&positron_ecout_energy);
    tree->SetBranchAddress("positron_chi2pid",&positron_chi2pid);
    tree->SetBranchAddress("positron_energy",&positron_energy);
    tree->SetBranchAddress("positron_m2ecin",&positron_m2ecin);
    tree->SetBranchAddress("positron_m2ecout",&positron_m2ecout);
    tree->SetBranchAddress("positron_m2pcal",&positron_m2pcal);
    tree->SetBranchAddress("positron_sfecin",&positron_sfecin);
    tree->SetBranchAddress("positron_sfpcal",&positron_sfpcal);
    tree->SetBranchAddress("positron_sfecout",&positron_sfecout);

    tree->SetBranchAddress("number_protons",&number_protons);
    tree->SetBranchAddress("proton_p",&proton_p);
    tree->SetBranchAddress("proton_theta",&proton_theta);
    tree->SetBranchAddress("proton_phi",&proton_phi);
    tree->SetBranchAddress("proton_px",&proton_px);
    tree->SetBranchAddress("proton_py",&proton_py);
    tree->SetBranchAddress("proton_pz",&proton_pz);
    tree->SetBranchAddress("proton_vx",&proton_vx);
    tree->SetBranchAddress("proton_vy",&proton_vy);
    tree->SetBranchAddress("proton_vz",&proton_vz);
    tree->SetBranchAddress("proton_vt",&proton_vt);
    tree->SetBranchAddress("proton_beta",&proton_beta);
    tree->SetBranchAddress("proton_chi2pid",&proton_chi2pid);
    tree->SetBranchAddress("proton_energy",&proton_energy);


    tree->SetBranchAddress("electron_photon",&electron_photon);
    tree->SetBranchAddress("positron_photon",&positron_photon);
    tree->SetBranchAddress("electron_photonE",&electron_photonE);
    tree->SetBranchAddress("positron_photonE",&positron_photonE);

    tree->SetBranchAddress("N_FT",&N_FT);
    tree->SetBranchAddress("FT_E",&FT_E);
    tree->SetBranchAddress("FT_Px",&FT_Px);
    tree->SetBranchAddress("FT_Py",&FT_Py);
    tree->SetBranchAddress("FT_Pz",&FT_Pz);
    tree->SetBranchAddress("FT_xFcal",&FT_xFcal);
    tree->SetBranchAddress("FT_yFcal",&FT_yFcal);
    tree->SetBranchAddress("FT_zFcal",&FT_zFcal);
    tree->SetBranchAddress("R",&R);
    tree->SetBranchAddress("FT_vtFcal",&FT_vtFcal);
    tree->SetBranchAddress("FT_vt",&FT_vt);

    /////////Instanciate QADB///////////
    QADB *qa = new QADB();


    //ANALYSIS
    //FILE *f_results;
    //FILE *f_events;

    //string nameFile="Pass1_v25";
    //f_results = fopen(("/lustre19/expphy/volatile/clas12/mtenorio/"+nameFile+"_dat_untagged.txt").c_str(), "w");
    //f_events = fopen(("/lustre19/expphy/volatile/clas12/mtenorio/"+nameFile+"_bgk.csv").c_str(), "w");

    Float_t cuts[15]={-0.6,-0.5,-0.4, -0.3,-0.2, -0.15,-0.1,-0.06,-0.04,-0.02,0.0,0.05,0.1, 0.2,0.3};
    TH1F* h_invariant_cuts[15];
    for(int j = 0; j < 15; j++){
      ostringstream name;
      name << "invariant_" << j;
      ostringstream sstr;
      sstr<< cuts[j]<<";M(e+e-),GeV ";
      h_invariant_cuts[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),100,2.0,3.5);;
      name.str("");
      sstr.str("");
      
    }

    int count_runs=0;
    int past_run;
    //Start
    for(int fc=0;fc<tree->GetEntries();fc++) {//Run 5032 to 5419 // 6616 6783
      electronAccepted=false;
      positronAccepted=false;
      protonAccepted=false;
      tree->GetEntry(fc);

      if(fc==0){
        past_run=run_number;
      }


      if(QA){
        bool Keep_event = true;
        if(run_number==5442||run_number==6749)
          continue;

        Keep_event = qa->OkForAsymmetry(run_number, event_number);
        int bad_runs[] = {5610, 5615, 6631, 6757};
        bool Additional_bad_runs = false;
        Additional_bad_runs = (std::find(std::begin(bad_runs), std::end(bad_runs), run_number) != std::end(bad_runs));

        if (!Keep_event || Additional_bad_runs)
          continue;
      }
      


      if(fc==0){
        past_run=run_number;
        count_runs++;
      }
      else{
        if(run_number!=past_run){
          past_run=run_number;
          count_runs++;
        }
      }


      if(number_protons==1){

        //---------CUTS FOR ELECTRON---------
        if(electron_pcal_energy>0.07  && electron_p>1.7 && electron_p<beam.E()&& electron_pcal_v>9 && electron_pcal_w>9){
          if(version==-18||version==-19){
            if(electron_vz>-8 && electron_vz<4){
              if((electron_pcal_energy+electron_ecin_energy+electron_ecout_energy)/electron_p>0.15)
                electronAccepted=true;
            }
          }
          else{
            if((electron_pcal_energy+electron_ecin_energy+electron_ecout_energy)/electron_p>0.15)
              electronAccepted=true;
          }
        }

        //----------CUTS FOR POSITRON---------
        if(positron_pcal_energy>0.07  && positron_p>1.7 && abs(positron_chi2pid)<5 && positron_pcal_v>9&& positron_pcal_w>9){//
          if(((positron_ecin_energy+positron_ecout_energy)/positron_p)>=(0.195-positron_pcal_energy/positron_p)){
            if(version==+18){
              if(positron_vz>-8 && positron_vz<4){
                if((positron_pcal_energy+positron_ecin_energy+positron_ecout_energy)/positron_p>0.15)
                  positronAccepted=true;
              }                  }
              else{
                if((positron_pcal_energy+positron_ecin_energy+positron_ecout_energy)/positron_p>0.15)
                  positronAccepted=true;
              }
            }
          }

        //----------CUTS FOR PROTON---------
          if(proton_beta>0.1&&proton_chi2pid<10&& proton_p>0.4){
            protonAccepted=true;
          }


          if(electronAccepted && positronAccepted && protonAccepted){
            if(abs(electron_vt-positron_vt)<=1){
            //------ENERGY CORRECTION ELECTRON------------
              electron_E_recon=electron_energy+electron_photonE;
              electron_px_recon=electron_px*(electron_E_recon/electron_energy);
              electron_py_recon=electron_py*(electron_E_recon/electron_energy);
              electron_pz_recon=electron_pz*(electron_E_recon/electron_energy);
              TLorentzVector ele(electron_px_recon,electron_py_recon,electron_pz_recon,electron_E_recon);

            //------ENERGY CORRECTION ELECTRON------------
              positron_E_recon=positron_energy+positron_photonE;
              positron_px_recon=positron_px*(positron_E_recon/positron_energy);
              positron_py_recon=positron_py*(positron_E_recon/positron_energy);
              positron_pz_recon=positron_pz*(positron_E_recon/positron_energy);
              TLorentzVector pos(positron_px_recon,positron_py_recon,positron_pz_recon,positron_E_recon);

              TLorentzVector pro(proton_px,proton_py,proton_pz,proton_energy);

              miss=(beam+target)-(pos+ele+pro);
              invariant_ee=(pos+ele);
              Q2=2*Beam_E*(miss.P()-miss.Pz());

              h_mm_vs_invariantmass_t->Fill(invariant_ee.M(),miss.M2());

              Double_t sumE,sumP;
              sumE=pos.E()+ele.E()+pro.E();
              sumP=pos.Pz()+ele.Pz()+pro.Pz();

              h_sum_cut->Fill(invariant_ee.M(),sumE-sumP-0.938);


              invariantMass=invariant_ee.M();
              MM=miss.M2();
              t=2*0.938*(pro.E()-0.938);
              EPcut=sumE-sumP-0.938;
              Egamma=sumE-0.938;

              results->Fill();


              if((sumE-sumP-0.938)<-0.030)//above -0.025
                continue;

              if(abs(miss.M2())>0.4)
                continue;

              h_qq_vs_invariantmass->Fill(invariant_ee.M(),Q2);

              //if(Q2>0.5)
                //continue;

              h_Invariant->Fill(invariant_ee.M());

              for(int j=0;j<15;j++){
                if(score_pos>cuts[j]&&score_ele>cuts[j]){
                  h_invariant_cuts[j]->Fill(invariant_ee.M());
                }
              }
          }//abs(electron_vt-positron_vt)<=1
        }//electronAccepted && positronAccepted && protonAccepted
      }//number_protons==
    }//End for "Runs"

    //fclose(f_results);

    file2->Write();

    printf("\n Total runs: %d \n",count_runs);

    TCanvas *can = new TCanvas("can","canvas",200,10,700,700);
    string pdf_original=nameFile+outFile+".pdf";


    gStyle->SetOptFit(1111);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);


    //---------Using Gauss+polynomial-----------

    
    h_Invariant->Draw();
    TF1 *gp=g43(h_Invariant,5E2,3.09,0.08,80,70,15,2.5,3.5);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    h_Invariant->Draw();
    TF1 *ge=g44(h_Invariant,5E2,3.09,0.08,15,-2,2.5,3.5);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    h_Invariant->Draw();
    //TF1 *c=CB2(h_Invariant,1E2,1.3E2,60,0.01,3.1,25,25,5,2.7,3.3);//e-e+ untagged
    TF1 *cp=CB2(h_Invariant,5E2,0.5,1.5E5,0.08,3.09,80,70,15,2.5,3.5);//e-e+ untagged
    can->Print( (pdf_original+"(").c_str());
    can->Clear();

    h_Invariant->Draw();
    TF1 *ce=CB3(h_Invariant,cp->GetParameter(0),cp->GetParameter(1),cp->GetParameter(2),cp->GetParameter(3),cp->GetParameter(4),15,-2,2.55,3.5);//e-e+ untagged
    can->Print( (pdf_original+")").c_str());
    can->Clear();

   /* FILE *f_readme;

    f_readme = fopen("/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Events.txt", "a");
    fprintf(f_readme,"------------------------- \n");
    fprintf(f_readme,"Name pdf: %s \n",pdf_original.c_str());
    fprintf(f_readme,"QADB: %d  \n",QA);
    fprintf(f_readme,"Configuration : %d \n \n", version);
    fprintf(f_readme,"Gauss+pol(3) : %f \n", gp->GetParameter(0));
    fprintf(f_readme,"Gauss+expo : %f \n", ge->GetParameter(0));
    fprintf(f_readme,"CB+pol(3) : %f \n", cp->GetParameter(0));
    fprintf(f_readme,"CB+expo : %f \n", ce->GetParameter(0));
    fprintf(f_readme,"\n------------------------- \n");*/


    //h_sum_cut->Draw("'colz");
    //can->Print( (pdf_original+")").c_str());

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;

  }




