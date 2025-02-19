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
#include "hipo4/reader.h"
#include <iguana/algorithms/clas12/FTEnergyCorrection/Algorithm.h>
#include "clas12reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

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

string Name_file(string nameFile,int PASS=2, bool QA=true){
    string corrections,QA_s;
    if(PASS==0)
        corrections="_FTcorrON";
    else if(PASS==1)
        corrections="_FTpass1";
    else
        corrections="";


    if(!QA)
        QA_s="_QAOFF";
    else
        QA="";
    
    return nameFile+corrections+QA_s+"_t";
}



int taggedfromroot(string nameFile="S19", int version=-19, double Beam_E=10.2, int PASS=2, bool QA=true) {
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    //********************************
    //DECLARATION OF VARIABLES
    //********************************
    int max;
    //MISSING MOMENTUM AND INVARIANT MASS
    TLorentzVector miss;
    TLorentzVector invariant_ee;
    //FT electron vector
    TLorentzVector electronFT_vec; 
    //Kinematic variables
    Double_t Q2, W;

    //Corrected Electron and Positron E and Px, Py, Pz
    Double_t electron_E_recon;
    Double_t positron_E_recon;
    Double_t electron_px_recon, electron_py_recon, electron_pz_recon;
    Double_t positron_px_recon, positron_py_recon, positron_pz_recon;

    //CONFIG VARIABLES
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

    Float_t score;

    int N_FT;
    int number_of_electrons;
    int number_of_positrons;
    int JPSI=0;
    bool electronAccepted, positronAccepted;



     Double_t score_pos,score_ele;
      Double_t score_pos_6,score_ele_6,score_e_6,score_p_6;



    //EVENT SELECTION
    TH2F* h_electron_ECin_vs_PCAL = new TH2F("h_electron_ECin_vs_PCAL","Electron SF-ECin vs. SF-PCAL; SF_{PCAL};SF_{ECIN}",200,0.0,0.3,100,0.0,0.2);
    TH2F* h_positron_EC_vs_PCAL = new TH2F("h_positron_EC_vs_PCAL","Positrons SF-EC vs. SF-PCALSF_{PCAL};SF_{EC}",200,0.0,0.3,100,0.0,0.2);
    TH1F* h_electron_zvertex = new TH1F("h_electron_zvertex","z-vertex distribution electron;z-vertex,cm",200,-15,35);
    TH1F* h_vertex_timediff = new TH1F("h_vertex_timediff","vertex time difference e- e+ ; vt_e^- - vt_e^+, ns",120,-2,2);

    //Photon energy correction
    TH2D* h_delta_theta_vs_phi_positron;
    TH2D* h_delta_theta_vs_phi_electron;
    TH1I* h_electron_photons = new TH1I("h_electron_photons","Photons radiated from electrons",7,0,7);
    TH1I* h_positron_photons = new TH1I("h_positron_photons","Photons radiated from positrons",7,0,7);


    //TAGGED
    TH1F* h_NFTelectron = new TH1F("h_NFTelectron","# FT electrons; # electrons in FT",10,0,10);
    TH1F* h_vertex_timediff_FT_FD = new TH1F("h_vertex_timediff_FT_FD","Vertex time difference between FT_e^- and FD_e^-; #Delta t_v=FD_e^- - FT_e^-, ns;",350,-25,25);

    //KINEMATICS
    TH1F* h_E_photon = new TH1F("h_E_photon","Photon Energy",120,0,0);
    TH1F* h_Q2 = new TH1F("h_Q2","Q^{2};Q^{2};Counts",120,0,0);
    TH1F* h_Q2_jpsievents = new TH1F("h_Q2_jpsievents","Q^{2};Q^{2};Counts",120,0,0);
    TH1F* h_W = new TH1F("h_W","Hadronic Mass;W;Counts",120,4,4.6);
    TH1F* h_MM= new TH1F("h_MM","Missing Mass; MM; Counts ",90,0,2);

    TH2F* h_mm_vs_invariantmass_t= new TH2F("h_mm_vs_invariantmass_t","MM vs M(e^+e^-) for t+-2;M(e^+e^-), GeV;Missing Mass, GeV",150,0,4,150,0,3);
    TH2F* h_mm_vs_invariantmass_no_t= new TH2F("h_mm_vs_invariantmass_no_t","MM vs M(e^+e^-) outside time cut;M(e^+e^-), GeV;Missing Mass, GeV",150,0,4,150,0,3);

    //RESULTS RESULTS
    TH1F* h_Invariant= new TH1F("h_im",";M(e+e-),GeV ",90,2.0,4);
    TH1F* h_Invariant_not= new TH1F("h_Invariant_not","Outside time range;M(e+e-),GeV",90,2.0,4);


   ////weitghs /lustre19/expphy/volatile/clas12/mtenorio/weights

    TString root_file = "/w/hallb-scshelf2102/clas12/mtenorio/Analysis/"+nameFile+".root";
    //TString root_file = "/volatile/clas12/mtenorio/Root/"+nameFile+".root";
    TFile *file = new TFile(root_file,"READ");
    TTree *tree = (TTree*)file->Get("analysis");

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

    tree->SetBranchAddress("score_ele_6",&score_ele_6);
    tree->SetBranchAddress("score_pos_6",&score_pos_6);

    
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
    tree->SetBranchAddress("score_pos",&score_pos);

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
    int electron_sec,positron_sec,electronFT_sec;
    //tree->SetBranchAddress("electron_sec",&electron_sec);
    //tree->SetBranchAddress("positron_sec",&positron_sec);
    //tree->SetBranchAddress("electronFT_sec",&electronFT_sec);

    h_delta_theta_vs_phi_positron=(TH2D*)file->Get("h_delta_theta_vs_phi_positron");
    h_delta_theta_vs_phi_electron=(TH2D*)file->Get("h_delta_theta_vs_phi_electron");

    //--------SET OUT ROOT file---------
    
    TString out_file_name=Name_file(nameFile,PASS,QA);
    TString out_file="/lustre19/expphy/volatile/clas12/mtenorio/Root/"+out_file_name+".root";

    TFile *file2 = new TFile(out_file,"RECREATE");
    TTree *results = new TTree("results",out_file);

    Double_t invariantMass ,MM, Egamma,score_p,score_e,time_ee;
    results->Branch("invariantMass",&invariantMass,"invariantMass/d");
    results->Branch("Q2",&Q2,"Q2/d");
    results->Branch("MM",&MM,"MM/d");
    results->Branch("time_ee",&time_ee,"time_ee/d");
    results->Branch("Egamma",&Egamma,"Egamma/d");
    results->Branch("W",&W,"W/d");
    results->Branch("score_e",&score_e,"score_e/d");
    results->Branch("score_p",&score_p,"score_p/d");
    results->Branch("score_e_6",&score_e_6,"score_e_6/d");
    results->Branch("score_p_6",&score_p_6,"score_p_6/d");

    gSystem->Load("libIguanaAlgorithms");

    iguana::clas12::FTEnergyCorrection algo_correction;

    algo_correction.Start();


    //ANALYSIS
    //FILE *f_results;

    //f_results = fopen(("/lustre19/expphy/volatile/clas12/mtenorio/"+nameFile+"_dat.txt").c_str(), "w");

     /////////Instanciate QADB///////////
    QADB *qa = new QADB();
    //Variables 
    int count_runs=0;
    int past_run;

    //Start
    for(int fc=0;fc<tree->GetEntries();fc++) {//Run 5032 to 5419 // 6616 6783
        electronAccepted=false;
        positronAccepted=false;
        tree->GetEntry(fc);

        //If QA is active, then apply it
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

           
        //---------CUTS FOR ELECTRON---------
        if(electron_pcal_energy>0.07  && electron_p>1.95 && electron_p<beam.E()&& electron_pcal_v>9 && electron_pcal_w>9){
            if(version==-18||version==-19){
                h_electron_zvertex->Fill(electron_vz);
                if(electron_vz>-8 && electron_vz<4){
                    electronAccepted=true;
                    h_electron_ECin_vs_PCAL->Fill(electron_pcal_energy/electron_p,electron_ecin_energy/electron_p);//------PLOT-------
                }
            }
            else{
                electronAccepted=true;
                h_electron_ECin_vs_PCAL->Fill(electron_pcal_energy/electron_p,electron_ecin_energy/electron_p);//------PLOT-------
            }
        }
        //----------CUTS FOR POSITRON---------
        if(positron_pcal_energy>0.07  && positron_p>1.95 && abs(positron_chi2pid)<5 && positron_pcal_v>9&& positron_pcal_w>9){//
            if(((positron_ecin_energy+positron_ecout_energy)/positron_p)>=(0.195-positron_pcal_energy/positron_p)){
                //h_positronsCUTS->Fill(positron_p);
                if(version==+18){
                    h_electron_zvertex->Fill(positron_vz);
                    if(positron_vz>-8 && positron_vz<4){
                        positronAccepted=true;
                        h_positron_EC_vs_PCAL->Fill(positron_pcal_energy/positron_p,(positron_ecin_energy+positron_ecout_energy)/positron_p);//------PLOT-------
                    }
                }
                else{
                    positronAccepted=true;
                    h_positron_EC_vs_PCAL->Fill(positron_pcal_energy/positron_p,(positron_ecin_energy+positron_ecout_energy)/positron_p);//------PLOT-------
                }
            }
        }

        if(electronAccepted && positronAccepted){
            h_vertex_timediff->Fill(electron_vt-positron_vt);//------PLOT-------
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

                h_electron_photons->Fill(electron_photon);//------PLOT-------
                h_positron_photons->Fill(positron_photon);//------PLOT-------
                h_NFTelectron->Fill(N_FT);//------PLOT-------

                

                TLorentzVector eleFT(FT_Px,FT_Py,FT_Pz,FT_E);

                if(N_FT==1){

                    FT_E_new=FT_E;

                    //NOT USE IN PASS2
                    if(PASS==1){
                        FT_E_new=-0.03689 + (1.1412*FT_E) - (0.04316*pow(FT_E,2))+ (0.007046*pow(FT_E,3))- (0.0004055*pow(FT_E,4));
                        FT_Px_new=FT_Px*(FT_E_new/FT_E);
                        FT_Py_new=FT_Py*(FT_E_new/FT_E);
                        FT_Pz_new=FT_Pz*(FT_E_new/FT_E);
                    }
                    else if (PASS==0){
                        if(version==-19){
                            FT_E_new = FT_E + 0.085643 - 0.0288063*FT_E +0.00894691*pow(FT_E,2) - 0.000725449*pow(FT_E,3);
                            FT_Px_new=FT_Px*(FT_E_new/FT_E);
                            FT_Py_new=FT_Py*(FT_E_new/FT_E);
                            FT_Pz_new=FT_Pz*(FT_E_new/FT_E);

                        }
                        else{
                            auto FT_4vect = algo_correction.Transform(FT_Px,FT_Py,FT_Pz,FT_E);
                            FT_Px_new=std::get<0>(FT_4vect);
                            FT_Py_new=std::get<1>(FT_4vect);
                            FT_Pz_new=std::get<2>(FT_4vect);
                            FT_E_new=std::get<3>(FT_4vect);
                        }
                    }
                    else{
                        FT_E_new=FT_E;
                        FT_Px_new=FT_Px*(FT_E_new/FT_E);
                        FT_Py_new=FT_Py*(FT_E_new/FT_E);
                        FT_Pz_new=FT_Pz*(FT_E_new/FT_E);
                    }

                    

                    //set variables
                    electronFT_vec.SetPxPyPzE(FT_Px_new,FT_Py_new,FT_Pz_new,FT_E_new);


                   



                    Double_t eFT_theta=electronFT_vec.Theta()*57.2958;
                    Double_t eFT_phi=electronFT_vec.Phi()*57.2958;
                    Double_t pos_theta=pos.Theta()*57.2958;
                    Double_t pos_phi=pos.Phi()*57.2958;
                    Double_t ele_theta=ele.Theta()*57.2958;
                    Double_t ele_phi=ele.Phi()*57.2958;

                    //MISSING MASS & INVARIANT MASS
                    miss=(beam+target)-(pos+ele+electronFT_vec);
                    invariant_ee=(pos+ele);

                    Q2=2*Beam_E*FT_E_new*(1-cos(electronFT_vec.Theta()));
                    E_photon=Beam_E-FT_E_new;
                    W=sqrt((0.938*0.938)+(2*0.938*E_photon)-Q2);

                    h_E_photon->Fill(E_photon);//------PLOT-------

                    h_vertex_timediff_FT_FD->Fill(electron_vt-FT_vt);//------PLOT-------

                    invariantMass=invariant_ee.M();
                    MM=miss.M();
                    Egamma=E_photon;
                    score_e=score_ele;
                    score_p=score_pos;
                    score_e_6=score_ele_6;
                    score_p_6=score_pos_6;
                    
                    time_ee=electron_vt-FT_vt;
                    results->Fill();


                   


                    //MM vs IM
                    if((time_ee)<=2 && (time_ee)>=-2){
                        h_mm_vs_invariantmass_t->Fill(invariant_ee.M(),miss.M());
                        h_MM->Fill(miss.M());
                        if(0.7 <=miss.M() &&  miss.M()<= 1.3){//1 sigmas 0.65 <=miss.M() &&  miss.M()<= 1.3
                            h_Invariant->Fill(invariant_ee.M());
                            h_Q2->Fill(Q2);//------PLOT-------
                            if(invariant_ee.M()>=3.0 && invariant_ee.M()<=3.2){
                                h_Q2_jpsievents->Fill(Q2);
                                h_W->Fill(W);
                            }

                        }

                    }
                    else{
                        h_mm_vs_invariantmass_no_t->Fill(invariant_ee.M(),miss.M());
                        h_Invariant_not->Fill(invariant_ee.M());
                    }


                }//N_FT==1
            }
        }
        //printf("run = %d\n",fc);
    }//End for "Runs"

    //fclose(f_results);

    algo_correction.Stop();

    file2->Write();
   
    TCanvas *can = new TCanvas("can","canvas",200,10,700,700);
    TString pdf_original="./R_Tagged/"+out_file_name+".pdf";

    can->Divide(2,2);

    can->cd(1);
    h_electron_ECin_vs_PCAL->Draw("colz");

    can->cd(2);
    h_positron_EC_vs_PCAL->Draw("colz");

    can->cd(3);
    h_electron_zvertex->Draw();

    can->cd(4);
    gPad->SetLogy();
    h_vertex_timediff->Draw();

    can->Print(pdf_original + "(");

    can->Clear();
    can->Divide(2,2);
    gPad->SetLogy(0);
    can->cd(1);
    h_delta_theta_vs_phi_electron->Draw("colz");
    can->cd(2);
    h_delta_theta_vs_phi_positron->Draw("colz");
    can->cd(3);
    gPad->SetLogy();
    h_electron_photons->Draw();
    can->cd(4);
    gPad->SetLogy();
    h_positron_photons->Draw();
    can->Print(pdf_original + "(");


    can->Clear();
    can->Divide(1,2);
    can->cd(1);
    //gPad->SetLogy();
    h_NFTelectron->Draw();

    can->cd(2);
    gPad->SetLogy();
    h_vertex_timediff_FT_FD->Draw();
    can->Print(pdf_original + "(");

    can->Clear();
    gPad->SetLogy(0);

    h_MM->Draw();
    can->Print(pdf_original + "(");


    h_mm_vs_invariantmass_t->Draw("colz");
    can->Print(pdf_original + "(");

    h_mm_vs_invariantmass_no_t->Draw("colz");
    can->Print(pdf_original + "(");


    h_Invariant->Draw();
    can->Print(pdf_original + "(");

    h_W->Draw();
    can->Print( pdf_original + "(");

    h_Q2->Draw();
    can->Print(pdf_original + "(");

    h_Q2_jpsievents->Draw();
    can->Print( pdf_original + ")");


    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;

}




