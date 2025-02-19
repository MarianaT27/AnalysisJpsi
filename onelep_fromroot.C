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

int onelep_fromroot(string nameFile = "S19", int version = -19, int top = 3, int PASS = 2,bool Restricted=false,bool AI = false,TString AIModel="12", bool QA = false)
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    // Set Beam energy
    double Beam_E;

    if (abs(version) == 18)
        Beam_E = 10.6;
    else
        Beam_E = 10.2;
    //********************************
    // DECLARATION OF VARIABLES
    //********************************
    int max;
    // MISSING MOMENTUM AND INVARIANT MASS
    TLorentzVector miss_eep;
    TLorentzVector miss_ep;
    // FT electron vector
    TLorentzVector eleFT;
    // Kinematic variables
    Double_t Q2, W;

    // Corrected Electron and Positron E and Px, Py, Pz
    Double_t electron_E_recon;
    Double_t lepton_E_recon;
    Double_t lepton_px_recon, lepton_py_recon, lepton_pz_recon;
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

    int N_FT;
    int number_of_electrons;
    int number_of_positrons, number_protons;
    int JPSI = 0;
    int number_of_piplus=0;
    int number_of_piminus=0;
    int number_of_gamma=0;
    int number_others=0;
    bool electronAccepted, positronAccepted, protonAccepted, leptonAccepted;

    Double_t score_pos, score_ele;

    // 2D M(e+e-) vs MX(e'p')
    // TH2F* h_IMee_vs_Mxep= new TH2F("h_IMee_vs_Mxep","M(e^{+}e^{-}) vs M_{X}(e'p');M(e^{+}e^{-}), GeV;M_{X}(e'p'), GeV",150,0,4,150,0,4);
    // Invariant Mass
    // TH1F* h_IMee= new TH1F("h_IMee",";M(e^{+}e^{-}),GeV ",100,0,4);
    // Missing Mass plots
    TH1F *h_Mxep = new TH1F("h_Mxep", "M_{X}(e'p'); M_{x}, GeV; Counts ", 100, 2.0, 4);

    TH1F *h_MXepe = new TH1F("h_MXepe", "M_{X}(e'p'e); M_{x}, GeV; Counts ", 100, -2, 2);
    TH1F *h_MXee = new TH1F("h_MXee", "M_{X}(e'e); M_{x}, GeV; Counts ", 100, 0, 0);
    TH1F *h_MXpe = new TH1F("h_MXpe", "M_{X}(p'e); M_{x}, GeV; Counts ", 100, 0, 0);

    TH2F *h2_Mxep = new TH2F("h2_Mxep", "M_{X}(e'p'); M_{x}, GeV; Counts ", 100, 2.0, 4, 100, 2.0, 4);
    TH2F *h2_MXepe = new TH2F("h2_MXepe", "M_{X}(e'p'e); M_{x}, GeV;  ", 100, -2, 2, 100, 2.0, 4);
    TH2F *h2_MXee = new TH2F("h2_MXee", "M_{X}(e'e); M_{x}, GeV;  ", 100, 0, 0, 100, 2.0, 4);
    TH2F *h2_MXpe = new TH2F("h2_MXpe", "M_{X}(p'e); M_{x}, GeV;  ", 100, 0, 0, 100, 2.0, 4);

    // KINEMATICS
    TH1F *h_E_photon = new TH1F("h_E_photon", "Photon Energy", 120, 0, 0);
    TH1F *h_Q2 = new TH1F("h_Q2", "Q^{2};Q^{2};Counts", 120, 0, 0);
    TH1F *h_W = new TH1F("h_W", "Hadronic Mass;W;Counts", 120, 4, 4.6);
    TH2F *h_beta_vs_p = new TH2F("h_beta_vs_p", ";P, GeV;#beta", 150, 0, 9, 150, 0, 0);

    TString config;
    if(version==18){
        config="F18out";
    }
    else if(version==-18){
        config="F18in";
    }
    else
        config="S19";

    TString root_file = "/w/hallb-scshelf2102/clas12/mtenorio/Analysis/" + nameFile + ".root";
    // TString root_file = "/volatile/clas12/mtenorio/Root/"+nameFile+".root";
    TFile *file = new TFile(root_file, "READ");
    TTree *tree = (TTree *)file->Get("analysis");

    cout<<"Opening: "<<root_file<<endl;

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

    tree->SetBranchAddress("number_others", &number_others);
    tree->SetBranchAddress("number_of_piplus", &number_of_piplus);
    tree->SetBranchAddress("number_of_piminus", &number_of_piminus);
    tree->SetBranchAddress("number_of_gamma", &number_of_gamma);

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

    tree->SetBranchAddress("N_FT", &N_FT);
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

    //--------SET OUT ROOT file---------

    TString out_file_name;

    out_file_name = Name_file(nameFile, top);

    TString out_file;
    if(AI)
        out_file = "/lustre24/expphy/volatile/clas12/mtenorio/Root/" + out_file_name + "_"+ AIModel+".root";
    else
        out_file = "/lustre24/expphy/volatile/clas12/mtenorio/Root/" + out_file_name + ".root";

    TFile *file2 = new TFile(out_file, "RECREATE");
    TTree *results = new TTree("results", out_file);

    Double_t Egamma, score_p, score_e, time_ee, time_ep, time_efdp, time_posfdp;
    /*Float_t Mxep, Mxepe, Mxee, Mxpe;
    Double_t electronFT_P, electronFT_Theta;
    Double_t protonFD_P, protonFD_Theta;
    Double_t electronFD_P, electronFD_Theta;
    Double_t D_Phi2,D_Phi1,D_Phi3;*/
    Double_t scoreBDT, scoreBDTG;
    Float_t Mxep, Mxepe, Mxee, Mxpe;
    Float_t electronFT_P, electronFT_Theta;
    Float_t protonFD_P, protonFD_Theta;
    Float_t electronFD_P, electronFD_Theta;
    Float_t D_Phi2,D_Phi1,D_Phi3;
    Float_t electronFT_Phi,electronFD_Phi,protonFD_Phi;

    /*
    results->Branch("Q2", &Q2, "Q2/d");
    results->Branch("time_ee", &time_ee, "time_ee/d");
    results->Branch("Egamma", &Egamma, "Egamma/d");
    results->Branch("W", &W, "W/d");
    results->Branch("time_ep", &time_ep, "time_ep/d");

    results->Branch("Mxep", &Mxep, "Mxep/d");
    
    if (top == 3)
    {
        results->Branch("score_p", &score_p, "score_p/d");
        results->Branch("time_posfdp", &time_posfdp, "time_posfdp/d");
    }
    else
    {
        results->Branch("score_e", &score_e, "score_e/d");
        results->Branch("time_efdp", &time_efdp, "time_efdp/d");
    }
    
    */

    

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
    if(AI){
        results->Branch("scoreBDT", &scoreBDT, "scoreBDT/d");
        results->Branch("scoreBDTG", &scoreBDTG, "scoreBDTG/d");
        }

    

    

    // Create a set of variables and declare them to the reader
    
    
    TMVA::Reader *reader = new TMVA::Reader("!Color:Silent");
    
    // Create the Reader object
    if (AI)
    {
        reader->AddVariable("electronFT_P", &electronFT_P);
        reader->AddVariable("electronFT_Theta", &electronFT_Theta);
        if(AIModel=="12")
            reader->AddVariable("D_Phi1", &D_Phi1);

        reader->AddVariable("electronFD_P", &electronFD_P);
        reader->AddVariable("electronFD_Theta", &electronFD_Theta);
        reader->AddVariable("D_Phi3", &D_Phi3);

        reader->AddVariable("protonFD_P", &protonFD_P);
        reader->AddVariable("protonFD_Theta", &protonFD_Theta);
        if(AIModel=="12")
            reader->AddVariable("D_Phi2", &D_Phi2);

        reader->AddVariable("Mxepe", &Mxepe);
        reader->AddVariable("Mxpe", &Mxpe);
        reader->AddVariable("Mxee", &Mxee);
        

        TString dir = "/volatile/clas12/mtenorio/Weights_AI/TMVA_Training_Var"+AIModel+"/weights/";
        TString prefix = "TMVAClassification";

        // Book Methods
        TString weightfile;

        weightfile = dir + prefix + "_BDT" + ".weights.xml";
        reader->BookMVA("BDT method", weightfile);

        weightfile = dir + prefix + "_BDTG" + ".weights.xml";
        reader->BookMVA("BDTG method", weightfile);
    }


    TMVA::Reader *readerTMVA = new TMVA::Reader("!Color:Silent");

    int model = 6;
    // Create a set of variables and declare them to the reader
    Float_t P, Theta, Phi;
    Float_t Nphe;

    Float_t PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;

    readerTMVA->AddVariable("SFPCAL", &PCAL);
    readerTMVA->AddVariable("SFECIN", &ECIN);
    readerTMVA->AddVariable("SFECOUT", &ECOUT);
    readerTMVA->AddVariable("m2PCAL", &m2PCAL);
    readerTMVA->AddVariable("m2ECIN", &m2ECIN);
    readerTMVA->AddVariable("m2ECOUT", &m2ECOUT);

    // Book Methods
    TString weightfile_pos;
    TString weightfile_ele;
    if (version == -19)
    {
        weightfile_ele = "/work/clas12/mtenorio/ML_weights_pass2/S19neg/TMVAClassification_BDT_6.weights.xml";
        weightfile_pos = "/work/clas12/mtenorio/ML_weights_pass2/S19pos/TMVAClassification_BDT_6.weights.xml";
    }
    if (version == -18)
    {
        weightfile_ele = "/work/clas12/mtenorio/ML_weights_pass2/F18inneg/TMVAClassification_BDT_6.weights.xml";
        weightfile_pos = "/work/clas12/mtenorio/ML_weights_pass2/F18inpos/TMVAClassification_BDT_6.weights.xml";
    }
    if (version == +18)
    {
        weightfile_ele = "/work/clas12/mtenorio/ML_weights_pass2/F18outneg/TMVAClassification_BDT_6.weights.xml";
        weightfile_pos = "/work/clas12/mtenorio/ML_weights_pass2/F18outpos/TMVAClassification_BDT_6.weights.xml";
    }
    readerTMVA->BookMVA("BDT pos method", weightfile_pos);
    readerTMVA->BookMVA("BDT ele method", weightfile_ele);

    /////////Instanciate QADB///////////
    //QADB *qa = new QADB();
    // Variables
    int count_runs = 0;
    int past_run;
    int number_misslep = 0;

    int percentageStep = 5;
    int step = tree->GetEntries() * percentageStep / 100;

    // Start
    for (int fc = 0; fc < tree->GetEntries(); fc++)
    { // Run 5032 to 5419 // 6616 6783
        leptonAccepted = false;
        protonAccepted = false;
        tree->GetEntry(fc);

        if (fc % step == 0)
        {
            double percentage = (fc * 100.0) / tree->GetEntries();
            std::cout << "Progress: " << percentage << "%" << std::endl;
            cout<<"Events found: "<<count_runs<<endl;
        }

        // If QA is active, then apply it
        if (QA)
        {
           /* bool Keep_event = true;
            if (run_number == 5442 || run_number == 6749)
                continue;

            Keep_event = qa->OkForAsymmetry(run_number, event_number);
            int bad_runs[] = {5610, 5615, 6631, 6757};
            bool Additional_bad_runs = false;
            Additional_bad_runs = (std::find(std::begin(bad_runs), std::end(bad_runs), run_number) != std::end(bad_runs));

            if (!Keep_event || Additional_bad_runs)
                continue;*/
        }

        // If the missing particle is a electron, then we need 1 e+ in FD
        if (top == 3 && (number_of_positrons != 1 || number_protons != 1 || N_FT != 1))
        {
            continue;
        }

        if (top == 4 && (number_of_electrons != 1 || number_protons != 1 || N_FT != 1))
        {
            continue;
        }

        if (Restricted && (number_others > 0 || number_of_piminus > 0 || number_of_piplus > 0))
            continue;

        if (top == 3)
            leptonAccepted = Accept_positron(positron_ecin_energy, positron_ecout_energy, positron_pcal_energy, positron_p, positron_chi2pid, positron_pcal_v, positron_pcal_w, positron_vz, version);
        else
            leptonAccepted = Accept_electron(electron_pcal_energy, electron_p, electron_pcal_v, electron_pcal_w, electron_vz, version);

        if (proton_beta > 0.1 && proton_chi2pid < 10 && proton_p > 0.4)
            protonAccepted = true;

        if (leptonAccepted && protonAccepted)
        {
            Double_t SF_corr, m2_corr;
            if (top == 3)
            {
                if (version == -19)
                {
                    SF_corr = 0.01;
                    m2_corr = 0.8;
                }
                else if (version == 18)
                {
                    SF_corr = 0.03;
                    m2_corr = 1.0;
                }
                else
                {
                    SF_corr = 0.01;
                    m2_corr = 0.8;
                }
                PCAL = positron_sfpcal;
                ECIN = positron_sfecin + SF_corr;
                ;
                ECOUT = positron_sfecout;
                m2PCAL = positron_m2pcal;
                m2ECIN = positron_m2ecin * m2_corr;
                m2ECOUT = positron_m2ecout;
                score_pos = readerTMVA->EvaluateMVA("BDT pos method");

                lepton_E_recon = positron_energy + positron_photonE;
                lepton_px_recon = positron_px * (lepton_E_recon / positron_energy);
                lepton_py_recon = positron_py * (lepton_E_recon / positron_energy);
                lepton_pz_recon = positron_pz * (lepton_E_recon / positron_energy);
            }
            else if (top == 4)
            {
                if (version == -19)
                {
                    SF_corr = 0.03;
                    m2_corr = 0.8;
                }
                else if (version == 18)
                {
                    SF_corr = 0.05;
                    m2_corr = 1.1;
                }
                else
                {
                    SF_corr = 0.02;
                    m2_corr = 0.8;
                }
                PCAL = electron_sfpcal;
                ECIN = electron_sfecin + SF_corr;
                ECOUT = electron_sfecout;
                m2PCAL = electron_m2pcal;
                m2ECIN = electron_m2ecin * m2_corr;
                m2ECOUT = electron_m2ecout;
                score_ele = readerTMVA->EvaluateMVA("BDT ele method");

                lepton_E_recon = electron_energy + electron_photonE;
                lepton_px_recon = electron_px * (lepton_E_recon / electron_energy);
                lepton_py_recon = electron_py * (lepton_E_recon / electron_energy);
                lepton_pz_recon = electron_pz * (lepton_E_recon / electron_energy);
            }
            TLorentzVector lepton(lepton_px_recon, lepton_py_recon, lepton_pz_recon, lepton_E_recon);
            if (PASS == 1)
            {
                FT_E_new = -0.03689 + (1.1412 * FT_E) - (0.04316 * pow(FT_E, 2)) + (0.007046 * pow(FT_E, 3)) - (0.0004055 * pow(FT_E, 4));
                FT_Px_new = FT_Px * (FT_E_new / FT_E);
                FT_Py_new = FT_Py * (FT_E_new / FT_E);
                FT_Pz_new = FT_Pz * (FT_E_new / FT_E);
            }
            else if (PASS == 2)
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
            double Ppro_ec = proton_pcorr(proton_theta, proton_p);
            double scale_factor = Ppro_ec / proton_p;

            double Ppx_new = proton_px * scale_factor;
            double Ppy_new = proton_py * scale_factor;
            double Ppz_new = proton_pz * scale_factor;

            TLorentzVector pro(Ppx_new, Ppy_new, Ppz_new, proton_energy);

            miss_eep = (beam + target) - (lepton + eleFT + pro); // Other lepton
            miss_ep = (beam + target) - (eleFT + pro);

            TLorentzVector miss_leptonFDp = (beam + target) - (lepton + pro);

            Mxepe = miss_eep.M2();
            Mxep = miss_ep.M();
            Mxee = ((beam + target) - (eleFT + lepton)).M2();
            Mxpe = miss_leptonFDp.M2();

            if (top == 3)
            {
                time_ee = positron_vt - FT_vt;
                time_posfdp = positron_vt - proton_vt;
            }
            else
            {
                time_ee = electron_vt - FT_vt;
                time_efdp = electron_vt - proton_vt;
            }

            time_ep = proton_vt - FT_vt;

            

            

            Q2 = 2 * Beam_E * FT_E_new * (1 - cos(eleFT.Theta()));
            E_photon = Beam_E - FT_E_new;
            W = sqrt((0.938 * 0.938) + (2 * 0.938 * E_photon) - Q2);
            

            if (E_photon < 8.1)
                continue;


            if (abs(time_ee) > 2.5)
                continue;
            


            

            count_runs++;

            electronFT_P = eleFT.P();
            electronFT_Theta = eleFT.Theta() * 57.2958;
            electronFT_Phi=eleFT.Phi() * 57.2958;
            D_Phi1 = fabs(eleFT.Phi() * 57.2958 -lepton.Phi() * 57.2958);
            protonFD_P = pro.P();
            protonFD_Theta = pro.Theta() * 57.2958;
            protonFD_Phi=pro.Phi() * 57.2958;
            D_Phi2 = fabs(eleFT.Phi() * 57.2958 -pro.Phi() * 57.2958);
            electronFD_P = lepton.P();
            electronFD_Theta = lepton.Theta() * 57.2958;
            electronFD_Phi=lepton.Phi() * 57.2958;
            D_Phi3 = fabs(lepton.Phi() * 57.2958 -pro.Phi() * 57.2958);


            scoreBDT = 99;
            scoreBDTG = 99;
            //cout<<"HERE"<<endl;

            //scoreBDT = reader->EvaluateMVA("BDT method");
            //scoreBDTG = reader->EvaluateMVA("BDTG method");
            //cout<<"HERE1"<<endl;
            results->Fill();
            //cout<<"HERE2"<<endl;

        }
    } // End for "Runs"

    cout<<count_runs<<endl;


    file2->Write();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;


}
