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
#include "bib/h_analysis.h"

using namespace TMVA;

// QADB header and namespace
#include "QADB.h"
using namespace QA;



TString Name_file(string nameFile, int PASS = 2, int version =-18)
{
    TString ver_s = "";
    TString pass_s = "";
    
    // Set version string based on version number
    if (version == -18)
        ver_s = "F18in";
    else if (version == 18)
        ver_s = "F18out";
    else if (version == -19)
        ver_s = "S19";
    
    // Set pass string if PASS is 0
    if (PASS == 0)
        pass_s = "_mc";
    
    // Combine the components into a single TString
    return TString(nameFile) + "_" + ver_s + pass_s;
}


int exclusivefromrootv2(string nameFile = "S19", TString key="BGK_3", int version = -19, int PASS = 2, bool QA = true)
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    ROOT::EnableImplicitMT(); // Enable multithreading, forcing ROOT to generate the dictionary
    double Beam_E;

    if(abs(version)==18)
        Beam_E=10.6;
    else    
        Beam_E=10.2;

    //********************************
    // DECLARATION OF VARIABLES
    //********************************
    int max;
    // MISSING MOMENTUM AND INVARIANT MASS
    TLorentzVector miss;
    TLorentzVector invariant_ee;
    TLorentzVector miss_ep;
    // FT electron vector
    TLorentzVector eleFT;
    // Kinematic variables
    Double_t Q2, W;

    // Corrected Electron and Positron E and Px, Py, Pz
    Double_t electron_E_recon;
    Double_t positron_E_recon;
    Double_t electron_px_recon, electron_py_recon, electron_pz_recon;
    Double_t positron_px_recon, positron_py_recon, positron_pz_recon;

    // CONFIG VARIABLES
    // FALL 10.6
    // SPRING 10.2
    TLorentzVector beam(0, 0, Beam_E, Beam_E);
    TLorentzVector target(0, 0, 0, 0.938);
    Int_t run_number;
    Int_t event_number;

    // Electron variables
    Double_t electron_p, electron_theta, electron_phi, electron_px, electron_py, electron_pz;
    Double_t electron_vx, electron_vy, electron_vz, electron_vt;
    Double_t electron_pcal_v, electron_pcal_w;
    Double_t electron_pcal_energy, electron_ecin_energy, electron_ecout_energy; // Electron Energies
    Double_t electron_energy;
    Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
    Double_t electron_sfpcal, electron_sfecin, electron_sfecout;

    // Positron variables
    Double_t positron_p, positron_theta, positron_phi, positron_px, positron_py, positron_pz;
    Double_t positron_vx, positron_vy, positron_vz, positron_vt;
    Double_t positron_pcal_v, positron_pcal_w;
    Double_t positron_pcal_energy, positron_ecin_energy, positron_ecout_energy; // Positron Energies
    Double_t positron_energy;
    Double_t positron_chi2pid;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;

    // Proton varianbles
    Double_t proton_p, proton_theta, proton_phi, proton_px, proton_py, proton_pz, proton_status, proton_beta, proton_chi2pid;
    Double_t proton_vx, proton_vy, proton_vz, proton_vt, proton_pid;
    Double_t proton_energy;

    // Photon variables
    Double_t photon_p, photon_theta, photon_phi, photon_px, photon_py, photon_pz;
    Double_t photon_vx, photon_vy, photon_vz, photon_vt, photon_energy;
    int electron_photon;
    int positron_photon;
    Double_t electron_photonE;
    Double_t positron_photonE;
    Double_t E_photon = 0;

    // Forward Tagger Electron variables
    Double_t FT_E;
    Double_t FT_Px;
    Double_t FT_Py;
    Double_t FT_Pz;
    Double_t FT_E_new;
    Double_t FT_Px_new;
    Double_t FT_Py_new;
    Double_t FT_Pz_new;
    Double_t FT_vt;
    Double_t FT_x, FT_y, FT_z, R;
    Double_t FT_vtC;
    // Testing
    Double_t FT_vtFcal, FT_xFcal, FT_yFcal, FT_zFcal, RFcal;

    Float_t score;
    Double_t time_ee;


    int N_FT;
    int number_of_electrons;
    int number_of_positrons, number_protons;
    int JPSI = 0;
    bool electronAccepted, positronAccepted, protonAccepted;
    int number_of_piplus, number_of_piminus, number_of_gamma;
    Double_t piplus_px, piplus_py, piplus_pz;
    Double_t piminus_px, piminus_py, piminus_pz;
    Double_t gamma_px, gamma_py, gamma_pz;

    Double_t score_pos, score_ele;
    // Invariant Mass
    TH2F* h_Mxepex= new TH2F("h_Mxepex",";M(e'p'e^{-}X),GeV ",100,-2,2,100, 2.0, 4);
    // Missing Mass plots
    TH1F *h_Mxep = new TH1F("h_Mxep", "M_{X}(e'p'); M_{x}, GeV; Counts ", 100, 2.0, 4);
  
    // RESULTS RESULTS
    TH1F *h_Invariant = new TH1F("h_im", ";M(e+e-),GeV ", 100, 2.0, 4);

    // Declare histograms for electronFT
    TH1F *h_P_electronFT = new TH1F("h_P_electronFT", "P for electronFT", 100, 0, 10);
    TH1F *h_Theta_electronFT = new TH1F("h_Theta_electronFT", "Theta for electronFT", 100, 0, 180);
    TH1F *h_Phi_electronFT = new TH1F("h_Phi_electronFT", "Phi for electronFT", 100, -180, 180);

    // Declare histograms for electronFD
    TH1F *h_P_electronFD = new TH1F("h_P_electronFD", "P for electronFD", 100, 0, 10);
    TH1F *h_Theta_electronFD = new TH1F("h_Theta_electronFD", "Theta for electronFD", 100, 0, 180);
    TH1F *h_Phi_electronFD = new TH1F("h_Phi_electronFD", "Phi for electronFD", 100, -180, 180);

    // Declare histograms for protonFD
    TH1F *h_P_protonFD = new TH1F("h_P_protonFD", "P for protonFD", 100, 0, 10);
    TH1F *h_Theta_protonFD = new TH1F("h_Theta_protonFD", "Theta for protonFD", 100, 0, 180);
    TH1F *h_Phi_protonFD = new TH1F("h_Phi_protonFD", "Phi for protonFD", 100, -180, 180);

    // Declare histograms for missingFD
    TH1F *h_Mxee = new TH1F("h_Mxee", "M_{X}(e'e^{-})", 100, 0, 4);
    TH1F *h_Mxepe = new TH1F("h_Mxepe", "M_{X}(e'p'e^{-})", 100, -2, 2);
    TH1F *h_Mxpe = new TH1F("h_Mxpe", "M_{X}(p'e^{-})", 100, 0, 4);

    ////weitghs /lustre19/expphy/volatile/clas12/mtenorio/weights

    TString root_file = "/w/hallb-scshelf2102/clas12/mtenorio/Analysis/"+nameFile+".root";
    //TString root_file = "/volatile/clas12/mtenorio/Root/" + nameFile + "_Full_FTJPsi.root";
    TFile *file = new TFile(root_file, "READ");
    TTree *tree = (TTree *)file->Get("analysis");

    tree->SetMakeClass(1);

    tree->SetBranchAddress("Beam_E", &Beam_E);

    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchAddress("event_number", &event_number);

    tree->SetBranchAddress("number_of_positrons", &number_of_positrons);
    tree->SetBranchAddress("number_of_electrons", &number_of_electrons);

    tree->SetBranchAddress("electron_p", &electron_p);
    tree->SetBranchAddress("electron_theta", &electron_theta);
    tree->SetBranchAddress("electron_phi", &electron_phi);
    tree->SetBranchAddress("electron_px", &electron_px);
    tree->SetBranchAddress("electron_py", &electron_py);
    tree->SetBranchAddress("electron_pz", &electron_pz);
    tree->SetBranchAddress("electron_vx", &electron_vx);
    tree->SetBranchAddress("electron_vy", &electron_vy);
    tree->SetBranchAddress("electron_vz", &electron_vz);
    tree->SetBranchAddress("electron_vt", &electron_vt);
    tree->SetBranchAddress("electron_pcal_v", &electron_pcal_v);
    tree->SetBranchAddress("electron_pcal_w", &electron_pcal_w);
    tree->SetBranchAddress("electron_pcal_energy", &electron_pcal_energy);
    tree->SetBranchAddress("electron_ecin_energy", &electron_ecin_energy);
    tree->SetBranchAddress("electron_ecout_energy", &electron_ecout_energy);
    tree->SetBranchAddress("electron_energy", &electron_energy);
    tree->SetBranchAddress("electron_m2ecin", &electron_m2ecin);
    tree->SetBranchAddress("electron_m2ecout", &electron_m2ecout);
    tree->SetBranchAddress("electron_m2pcal", &electron_m2pcal);
    tree->SetBranchAddress("electron_sfecin", &electron_sfecin);
    tree->SetBranchAddress("electron_sfpcal", &electron_sfpcal);
    tree->SetBranchAddress("electron_sfecout", &electron_sfecout);

    tree->SetBranchAddress("positron_p", &positron_p);
    tree->SetBranchAddress("positron_theta", &positron_theta);
    tree->SetBranchAddress("positron_phi", &positron_phi);
    tree->SetBranchAddress("positron_px", &positron_px);
    tree->SetBranchAddress("positron_py", &positron_py);
    tree->SetBranchAddress("positron_pz", &positron_pz);
    tree->SetBranchAddress("positron_vx", &positron_vx);
    tree->SetBranchAddress("positron_vy", &positron_vy);
    tree->SetBranchAddress("positron_vz", &positron_vz);
    tree->SetBranchAddress("positron_vt", &positron_vt);
    tree->SetBranchAddress("positron_pcal_v", &positron_pcal_v);
    tree->SetBranchAddress("positron_pcal_w", &positron_pcal_w);
    tree->SetBranchAddress("positron_pcal_energy", &positron_pcal_energy);
    tree->SetBranchAddress("positron_ecin_energy", &positron_ecin_energy);
    tree->SetBranchAddress("positron_ecout_energy", &positron_ecout_energy);
    tree->SetBranchAddress("positron_chi2pid", &positron_chi2pid);
    tree->SetBranchAddress("positron_energy", &positron_energy);
    tree->SetBranchAddress("positron_m2ecin", &positron_m2ecin);
    tree->SetBranchAddress("positron_m2ecout", &positron_m2ecout);
    tree->SetBranchAddress("positron_m2pcal", &positron_m2pcal);
    tree->SetBranchAddress("positron_sfecin", &positron_sfecin);
    tree->SetBranchAddress("positron_sfpcal", &positron_sfpcal);
    tree->SetBranchAddress("positron_sfecout", &positron_sfecout);

    tree->SetBranchAddress("number_protons", &number_protons);
    tree->SetBranchAddress("proton_p", &proton_p);
    tree->SetBranchAddress("proton_theta", &proton_theta);
    tree->SetBranchAddress("proton_phi", &proton_phi);
    tree->SetBranchAddress("proton_px", &proton_px);
    tree->SetBranchAddress("proton_py", &proton_py);
    tree->SetBranchAddress("proton_pz", &proton_pz);
    tree->SetBranchAddress("proton_vx", &proton_vx);
    tree->SetBranchAddress("proton_vy", &proton_vy);
    tree->SetBranchAddress("proton_vz", &proton_vz);
    tree->SetBranchAddress("proton_vt", &proton_vt);
    tree->SetBranchAddress("proton_beta", &proton_beta);
    tree->SetBranchAddress("proton_chi2pid", &proton_chi2pid);
    tree->SetBranchAddress("proton_energy", &proton_energy);

    tree->SetBranchAddress("electron_photon", &electron_photon);
    tree->SetBranchAddress("positron_photon", &positron_photon);
    tree->SetBranchAddress("electron_photonE", &electron_photonE);
    tree->SetBranchAddress("positron_photonE", &positron_photonE);
    
    tree->SetBranchAddress("number_of_piplus", &number_of_piplus);
    tree->SetBranchAddress("number_of_piminus", &number_of_piminus);
    tree->SetBranchAddress("number_of_gamma", &number_of_gamma);

    tree->SetBranchAddress("piplus_px", &piplus_px);
    tree->SetBranchAddress("piplus_py", &piplus_py);
    tree->SetBranchAddress("piplus_pz", &piplus_pz);

    tree->SetBranchAddress("piminus_px", &piminus_px);
    tree->SetBranchAddress("piminus_py", &piminus_py);
    tree->SetBranchAddress("piminus_pz", &piminus_pz);

    tree->SetBranchAddress("gamma_px", &gamma_px);
    tree->SetBranchAddress("gamma_py", &gamma_py);
    tree->SetBranchAddress("gamma_pz", &gamma_pz);

    tree->SetBranchAddress("FT_E", &FT_E);
    tree->SetBranchAddress("FT_Px", &FT_Px);
    tree->SetBranchAddress("FT_Py", &FT_Py);
    tree->SetBranchAddress("FT_Pz", &FT_Pz);
    tree->SetBranchAddress("FT_xFcal", &FT_xFcal);
    tree->SetBranchAddress("FT_yFcal", &FT_yFcal);
    tree->SetBranchAddress("FT_zFcal", &FT_zFcal);
    tree->SetBranchAddress("R", &R);
    tree->SetBranchAddress("FT_vtFcal", &FT_vtFcal);
    tree->SetBranchAddress("FT_vt", &FT_vt);
    tree->SetBranchAddress("N_FT", &N_FT);
    int electron_sec, positron_sec, electronFT_sec;

    
    TString out_file = "/lustre24/expphy/volatile/clas12/mtenorio/Root/"+nameFile+"_exclusive.root";

    TFile *file2 = new TFile(out_file, "RECREATE");
    TTree *results = new TTree("results", out_file);

    // Declare the 12 variables
    Float_t Mxep, Mxepe, Mxee, Mxpe;
    Float_t electronFT_P, electronFT_Theta;
    Float_t protonFD_P, protonFD_Theta;
    Float_t electronFD_P, electronFD_Theta;
    Float_t D_Phi2,D_Phi1,D_Phi3;
    Float_t electronFT_Phi,electronFD_Phi,protonFD_Phi;


    // Create branches for the variables
    results->Branch("electronFT_P", &electronFT_P, "electronFT_P/F");
    results->Branch("electronFT_Theta", &electronFT_Theta, "electronFT_Theta/F");
    
    results->Branch("protonFD_P", &protonFD_P, "protonFD_P/F");
    results->Branch("protonFD_Theta", &protonFD_Theta, "protonFD_Theta/F");
    
    results->Branch("electronFD_P", &electronFD_P, "electronFD_P/F");
    results->Branch("electronFD_Theta", &electronFD_Theta, "electronFD_Theta/F");

    results->Branch("D_Phi1", &D_Phi1, "D_Phi1/F");
    results->Branch("D_Phi2", &D_Phi2, "D_Phi2/F");
    results->Branch("D_Phi3", &D_Phi3, "D_Phi3/F");

    results->Branch("electronFD_Phi", &electronFD_Phi, "electronFD_Phi/F");
    results->Branch("electronFT_Phi", &electronFT_Phi, "electronFT_Phi/F");
    results->Branch("protonFD_Phi", &protonFD_Phi, "protonFD_Phi/F");
    

    results->Branch("Mxepe", &Mxepe, "Mxepe/F");
    results->Branch("Mxee", &Mxee, "Mxee/F");
    results->Branch("Mxpe", &Mxpe, "Mxpe/F");

    results->Branch("Mxep", &Mxep, "Mxep/F");
    

    // Variables
    int count_runs = 0;
    int past_run;
    int percentageStep = 5;
    int step = tree->GetEntries() * percentageStep / 100;


    int ep_eminus_piplus=0;
    // Start
    for (int fc = 0; fc < tree->GetEntries(); fc++)
    { // Run 5032 to 5419 // 6616 6783
        
        

        if (fc % step == 0)
        {
            double percentage = (fc * 100.0) / tree->GetEntries();
            std::cout << "Progress: " << percentage << "%" << std::endl;
            std::cout << "Events found so far: " << ep_eminus_piplus  << std::endl;
        }
        
        tree->GetEntry(fc);

        if (N_FT !=1 || number_of_electrons!=1||number_protons!=1){
            continue;
             //cout<<"HERE"<<endl;
        }
        if(key=="BGK_1_FTJPsi"&&number_of_piplus==0)
            continue;
        if(key=="BGK_2_FTJPsi"&&number_of_piminus==0)
            continue;
        if(key=="BGK_3_FTJPsi"&&number_of_gamma==0)
            continue;
        if(key=="SIG"&&number_of_positrons!=1)
            continue;

    


        



        
        bool electronAccepted = Accept_electron(electron_pcal_energy, electron_p, electron_pcal_v, electron_pcal_w, electron_vz, version);
        if(proton_beta<=0.1||proton_chi2pid>10 || proton_p<0.4||!electronAccepted){
            continue;
        }
        
        //------ENERGY CORRECTION ELECTRON------------
        electron_E_recon=electron_energy+electron_photonE;
        electron_px_recon=electron_px*(electron_E_recon/electron_energy);
        electron_py_recon=electron_py*(electron_E_recon/electron_energy);
        electron_pz_recon=electron_pz*(electron_E_recon/electron_energy);
        TLorentzVector ele(electron_px_recon,electron_py_recon,electron_pz_recon,electron_E_recon);

        //TLorentzVector ele(electron_px, electron_py, electron_pz, electron_energy);
        //TLorentzVector pos(piplus_px,piplus_py,piplus_pz,sqrt(piplus_px * piplus_px + piplus_py * piplus_py + piplus_pz * piplus_pz + 0.0005 * 0.0005));
        
        //TLorentzVector pos(piminus_px,piminus_py,piminus_pz,sqrt(piminus_px * piminus_px + piminus_py * piminus_py + piminus_pz * piminus_pz + 0.0005 * 0.0005));
        //TLorentzVector pos(positron_px, positron_py, positron_pz, positron_energy);
        

        if(key=="BGK_3"){
            TLorentzVector pos(gamma_px,gamma_py,gamma_pz,sqrt(gamma_px * gamma_px + gamma_py * gamma_py + gamma_pz * gamma_pz + 0.0 * 0.0));
            if(pos.E()<0.7)
                continue;
        }

            FT_E_new = FT_E;

            // NOT USE IN PASS2
            if (PASS == 2)
            {
                if (version == -19)
                {
                    FT_E_new = FT_E + 0.085643 - 0.0288063 * FT_E + 0.00894691 * pow(FT_E, 2) - 0.000725449 * pow(FT_E, 3);
                    FT_Px_new = FT_Px * (FT_E_new / FT_E);
                    FT_Py_new = FT_Py * (FT_E_new / FT_E);
                    FT_Pz_new = FT_Pz * (FT_E_new / FT_E);
                }
                else
                {
                   FT_E_new = FT_E + +  0.0208922 + 0.050158*FT_E- 0.0181107*pow(FT_E,2) + 0.00305671*pow(FT_E,3) - 0.000178235*pow(FT_E,4);
                    FT_Px_new = FT_Px * (FT_E_new / FT_E);
                    FT_Py_new = FT_Py * (FT_E_new / FT_E);
                    FT_Pz_new = FT_Pz * (FT_E_new / FT_E);
                }
            }
            else
            {
                FT_E_new = FT_E;
                FT_Px_new = FT_Px * (FT_E_new / FT_E);
                FT_Py_new = FT_Py * (FT_E_new / FT_E);
                FT_Pz_new = FT_Pz * (FT_E_new / FT_E);
            }

            // set variables
            eleFT.SetPxPyPzE(FT_Px_new, FT_Py_new, FT_Pz_new, FT_E_new);

            TLorentzVector pro(proton_px, proton_py, proton_pz, proton_energy);



            //Kinematic variables
            E_photon=Beam_E-FT_E_new;
            time_ee=electron_vt-FT_vt;
            if (E_photon < 8.1)
                continue;

            if (abs(time_ee) < 2.0)
                continue;

            TLorentzVector BeamTarget = (beam + target);
            // INVARIANT MASS
            //invariant_ee = (pos + ele);
            // MISSING MASS
            miss = BeamTarget - ( ele + eleFT + pro);
            miss_ep = BeamTarget - (eleFT + pro);

            electronFT_P = eleFT.P();
            electronFT_Theta = eleFT.Theta() * 57.2958;
            electronFT_Phi=eleFT.Phi() * 57.2958;
            D_Phi1 = fabs(eleFT.Phi() * 57.2958 -ele.Phi() * 57.2958);
            protonFD_P = pro.P();
            protonFD_Theta = pro.Theta() * 57.2958;
            protonFD_Phi=pro.Phi() * 57.2958;
            D_Phi2 = fabs(eleFT.Phi() * 57.2958 -pro.Phi() * 57.2958);
            electronFD_P = ele.P();
            electronFD_Theta = ele.Theta() * 57.2958;
            electronFD_Phi=ele.Phi() * 57.2958;
            D_Phi3 = fabs(ele.Phi() * 57.2958 -pro.Phi() * 57.2958);

            Mxepe=miss.M2();
            Mxee=(BeamTarget-eleFT-ele).M2();
            Mxpe=(BeamTarget-pro-ele).M2();
            Mxep=(BeamTarget-eleFT-pro).M();

            if (Mxee >= 2.0 && Mxee <= 8.0 &&
                Mxpe >= 0.0 && Mxpe <= 2.0 &&
                electronFT_P <= 2.5 &&
                electronFT_Theta <= 5.5 &&
                electronFD_P <= 7.0 &&
                electronFD_Theta >= 10.0 &&
                protonFD_P <= 4.0 &&
                protonFD_Theta <= 25.0) {
                ep_eminus_piplus++;
                results->Fill();
            }
        
        
        // printf("run = %d\n",fc);
    } // End for "Runs"
    cout<<ep_eminus_piplus<<endl;
    // fclose(f_results);
    file2->Write();




    


    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;
}
