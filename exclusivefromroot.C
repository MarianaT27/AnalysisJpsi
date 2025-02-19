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

int exclusivefromroot(string nameFile = "S19", int version = -19, int PASS = 2,bool Restricted=false,bool AI = false,TString AIModel="", bool QA = true)
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

    int N_FT;
    int number_of_electrons;
    int number_of_positrons, number_protons;
    int JPSI = 0;
    bool electronAccepted, positronAccepted, protonAccepted;
    int number_of_piplus=0;
    int number_of_piminus=0;
    int number_of_gamma=0;
    int number_others=0;

    Double_t score_pos, score_ele;

    // EVENT SELECTION
    TH2F *h_electron_ECin_vs_PCAL = new TH2F("h_electron_ECin_vs_PCAL", "Electron SF-ECin vs. SF-PCAL; SF_{PCAL};SF_{ECIN}", 200, 0.0, 0.3, 100, 0.0, 0.2);
    TH2F *h_positron_EC_vs_PCAL = new TH2F("h_positron_EC_vs_PCAL", "Positrons SF-EC vs. SF-PCALSF_{PCAL};SF_{EC}", 200, 0.0, 0.3, 100, 0.0, 0.2);
    TH1F *h_electron_zvertex = new TH1F("h_electron_zvertex", "z-vertex distribution electron;z-vertex,cm", 200, -15, 35);
    TH1F *h_vertex_timediff = new TH1F("h_vertex_timediff", "vertex time difference e- e+ ; vt_e^- - vt_e^+, ns", 120, -2, 2);
    TH2F *h_beta_vs_p = new TH2F("h_beta_vs_p", ";P, GeV;#beta", 150, 0, 9, 150, 0, 0);

    // Photon energy correction
    TH2D *h_delta_theta_vs_phi_positron;
    TH2D *h_delta_theta_vs_phi_electron;
    TH1I *h_electron_photons = new TH1I("h_electron_photons", "Photons radiated from electrons", 7, 0, 7);
    TH1I *h_positron_photons = new TH1I("h_positron_photons", "Photons radiated from positrons", 7, 0, 7);

    // TAGGED
    TH1F *h_p = new TH1F("h_p", "; P", 100, 0, 10);
    TH1F *h_p_corr = new TH1F("h_p_corr", "; P", 100, 0, 10);

    TH1F *h_NFTelectron = new TH1F("h_NFTelectron", "# FT electrons; # electrons in FT", 10, 0, 10);
    TH1F *h_vertex_timediff_FT_FD = new TH1F("h_vertex_timediff_FT_FD", "Vertex time difference between FT_e^- and FD_e^-; #Delta t_v=FD_e^- - FT_e^-, ns;", 350, -25, 25);
    TH1F *h_vertex_timediff_FT_eFD = new TH1F("h_vertex_timediff_FT_eFD", "Vertex time difference between FT_e^- and FD_e^+; #Delta t_v=FD_e^+ - FT_e^-, ns;", 350, -25, 25);
    TH1F *h_vertex_timediff_FT_pFD = new TH1F("h_vertex_timediff_FT_pFD", "Vertex time difference between FT_e^- and FD_p; #Delta t_v=FD_p - FT_e^-, ns;", 350, -25, 25);

    // KINEMATICS

    // 2D M(e+e-) vs MX(e'p')
    TH2F *h_IMee_vs_Mxep = new TH2F("h_IMee_vs_Mxep", "M(e^{+}e^{-}) vs M_{X}(e'p');M(e^{+}e^{-}), GeV;M_{X}(e'p'), GeV", 150, 0, 4, 150, 0, 4);
    // Invariant Mass
    // TH1F* h_IMee= new TH1F("h_IMee",";M(e^{+}e^{-}),GeV ",100,0,4);
    // Missing Mass plots
    TH1F *h_Mxep = new TH1F("h_Mxep", "M_{X}(e'p'); M_{x}, GeV; Counts ", 100, 2.0, 4);
    TH1F *h_MXepee = new TH1F("h_MXepee", "M_{X}(e'p'e^{+}e^{-}); M_{x}, GeV; Counts ", 100, -1, 1);
    TH1F *h_MXeee = new TH1F("h_MXeee", "M_{X}(e'e^{+}e^{-}); M_{x}, GeV; Counts ", 100, 0, 2);
    TH1F *h_MXepe1 = new TH1F("h_MXepe1", "M_{X}(e'p'e^{+}); M_{x}, GeV; Counts ", 100, -1, 1);
    TH1F *h_MXepe2 = new TH1F("h_MXepe2", "M_{X}(e'p'e^{-}); M_{x}, GeV; Counts ", 100, -1, 1);
    TH1F *h_MXee1 = new TH1F("h_MXee1", "M_{X}(e'e^{+}); M_{x}, GeV; Counts ", 100, 0, 4);
    TH1F *h_MXee2 = new TH1F("h_MXee2", "M_{X}(e'e^{-}); M_{x}, GeV; Counts ", 100, 0, 4);
    TH1F *h_MXpe1 = new TH1F("h_MXpe1", "M_{X}(p'e^{+}); M_{x}, GeV; Counts ", 100, -1, 3);
    TH1F *h_MXpe2 = new TH1F("h_MXpe2", "M_{X}(p'e^{-}); M_{x}, GeV; Counts ", 100, -1, 3);

    TH2F *h2_Invariant = new TH2F("h2_Invariant", ";M(e+e-),GeV ", 100, 2.0, 4, 100, 2.0, 4);
    TH2F *h2_MXep = new TH2F("h2_MXep", "M_{X}(e'p'); M_{x}, GeV;  ", 100, 2.0, 4, 100, 2.0, 4);
    TH2F *h2_MXeee = new TH2F("h2_MXeee", "M_{X}(e'e^{+}e^{-}); M_{x}, GeV;  ", 100, 0, 2, 100, 2.0, 4);
    TH2F *h2_MXepe1 = new TH2F("h2_MXepe1", "M_{X}(e'p'e^{+}); M_{x}, GeV;  ", 100, -1, 1, 100, 2.0, 4);
    TH2F *h2_MXepe2 = new TH2F("h2_MXepe2", "M_{X}(e'p'e^{-}); M_{x}, GeV;  ", 100, -1, 1, 100, 2.0, 4);
    TH2F *h2_MXee1 = new TH2F("h2_MXee1", "M_{X}(e'e^{+}); M_{x}, GeV;  ", 100, 0, 4, 100, 2.0, 4);
    TH2F *h2_MXee2 = new TH2F("h2_MXee2", "M_{X}(e'e^{-}); M_{x}, GeV;  ", 100, 0, 4, 100, 2.0, 4);
    TH2F *h2_MXpe1 = new TH2F("h2_MXpe1", "M_{X}(p'e^{+}); M_{x}, GeV;  ", 100, -1, 3, 100, 2.0, 4);
    TH2F *h2_MXpe2 = new TH2F("h2_MXpe2", "M_{X}(p'e^{-}); M_{x}, GeV;  ", 100, -1, 3, 100, 2.0, 4);

    // KINEMATICS
    TH1F *h_E_photon = new TH1F("h_E_photon", "Photon Energy", 120, 0, 0);
    TH1F *h_Q2 = new TH1F("h_Q2", "Q^{2};Q^{2};Counts", 120, 0, 0);
    TH1F *h_Q2_jpsievents = new TH1F("h_Q2_jpsievents", "Q^{2};Q^{2};Counts", 120, 0, 0);
    TH1F *h_W = new TH1F("h_W", "Hadronic Mass;W;Counts", 120, 4, 4.6);

    // RESULTS RESULTS
    TH1F *h_Invariant = new TH1F("h_im", ";M(e+e-),GeV ", 100, 2.0, 4);

    ////weitghs /lustre19/expphy/volatile/clas12/mtenorio/weights

    TString config;
    if(version==18){
        config="F18out";
    }
    else if(version==-18){
        config="F18in";
    }
    else
        config="S19";

    TString root_file = "/w/hallb-scshelf2102/clas12/mtenorio/Analysis/" + config + "_Full.root";
    // TString root_file = "/volatile/clas12/mtenorio/Root/"+nameFile+".root";
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

    tree->SetBranchAddress("number_of_piplus", &number_of_piplus);
    tree->SetBranchAddress("number_of_piminus", &number_of_piminus);
    tree->SetBranchAddress("number_of_gamma", &number_of_gamma);
    tree->SetBranchAddress("number_others", &number_others);

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
    int electron_sec, positron_sec, electronFT_sec;

    h_delta_theta_vs_phi_positron = (TH2D *)file->Get("h_delta_theta_vs_phi_positron");
    h_delta_theta_vs_phi_electron = (TH2D *)file->Get("h_delta_theta_vs_phi_electron");

    //--------SET OUT ROOT file---------

    TString out_file_name = Name_file(nameFile, 2);
    TString out_file = "/lustre24/expphy/volatile/clas12/mtenorio/Root/AI_Dec24/" + out_file_name + ".root";

    TFile *file2 = new TFile(out_file, "RECREATE");
    TTree *results = new TTree("results", out_file);

    Double_t invariantMass, Egamma, score_p, score_e, time_ee, time_epos, time_ep, time_efdp, time_posfdp;
    Double_t Mxep, MXepee, MXeee, MXepe1, MXepe2, MXee1, MXee2, MXpe1, MXpe2;
    Double_t scoreBDT, scoreBDTG;
    results->Branch("time_ee", &time_ee, "time_ee/d");
    results->Branch("time_ep", &time_ep, "time_ep/d");
    results->Branch("time_efdp", &time_efdp, "time_efdp/d");
    results->Branch("time_posfdp", &time_posfdp, "time_posfdp/d");
    results->Branch("time_epos", &time_epos, "time_epos/d");
    results->Branch("Egamma", &Egamma, "Egamma/d");
    results->Branch("W", &W, "W/d");
    results->Branch("Q2", &Q2, "Q2/d");
    results->Branch("invariantMass", &invariantMass, "invariantMass/d");

    results->Branch("Mxep", &Mxep, "Mxep/d");
    results->Branch("MXepee", &MXepee, "MXepee/d");
    results->Branch("MXeee", &MXeee, "MXeee/d");
    results->Branch("MXepe1", &MXepe1, "MXepe1/d");
    results->Branch("MXepe2", &MXepe2, "MXepe2/d");
    results->Branch("MXee1", &MXee1, "MXee1/d");
    results->Branch("MXee2", &MXee2, "MXee2/d");
    results->Branch("MXpe1", &MXpe1, "MXpe1/d");
    results->Branch("MXpe2", &MXpe2, "MXpe2/d");

    results->Branch("scoreBDT", &scoreBDT, "scoreBDT/d");
    results->Branch("scoreBDTG", &scoreBDTG, "scoreBDTG/d");


    results->Branch("score_e", &score_e, "score_e/d");
    results->Branch("score_p", &score_p, "score_p/d");

    gSystem->Load("libIguanaAlgorithms");

    iguana::clas12::FTEnergyCorrection algo_correction;

    algo_correction.Start();
        // Create a set of variables and declare them to the reader
    Float_t Mxepe, Mxee, Mxpe;
    Float_t electronFT_P, electronFT_Theta, electronFT_Phi;
    Float_t protonFD_P, protonFD_Theta, protonFD_Phi;
    Float_t electronFD_P, electronFD_Theta, electronFD_Phi;
    TMVA::Reader *reader = new TMVA::Reader("!Color:Silent");
    // Create the Reader object
    if (AI)
    {
        reader->AddVariable("electronFT_P", &electronFT_P);
        reader->AddVariable("electronFT_Theta", &electronFT_Theta);
        reader->AddVariable("electronFT_Phi", &electronFT_Phi);

        reader->AddVariable("protonFD_P", &protonFD_P);
        reader->AddVariable("protonFD_Theta", &protonFD_Theta);
        reader->AddVariable("protonFD_Phi", &protonFD_Phi);

        reader->AddVariable("electronFD_P", &electronFD_P);
        reader->AddVariable("electronFD_Theta", &electronFD_Theta);
        reader->AddVariable("electronFD_Phi", &electronFD_Phi);

        reader->AddVariable("Mxepe", &Mxepe);
        
        if(AIModel=="12"||AIModel=="12_MIX"||AIModel=="12_SIDIS"){
            reader->AddVariable("Mxpe", &Mxpe);
            reader->AddVariable("Mxee", &Mxee);
        }

        TString dir = "/volatile/clas12/mtenorio/weights_TaggedAI/"+AIModel+"/";
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
    Float_t P, Theta, Phi, PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;
    Float_t Nphe;

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
    QADB *qa = new QADB();
    // Variables
    int count_runs = 0;
    int past_run;
    int percentageStep = 5;
    int step = tree->GetEntries() * percentageStep / 100;

    // Start
    for (int fc = 0; fc < tree->GetEntries(); fc++)
    { // Run 5032 to 5419 // 6616 6783
        electronAccepted = false;
        positronAccepted = false;
        protonAccepted = false;

        

        tree->GetEntry(fc);

        if (fc % step == 0)
        {
            double percentage = (fc * 100.0) / tree->GetEntries();
            std::cout << "Progress: " << percentage << "%" << std::endl;
        }

        // If QA is active, then apply it
        if (QA)
        {
            bool Keep_event = true;
            if (run_number == 5442 || run_number == 6749)
                continue;

            Keep_event = qa->OkForAsymmetry(run_number, event_number);
            int bad_runs[] = {5610, 5615, 6631, 6757};
            bool Additional_bad_runs = false;
            Additional_bad_runs = (std::find(std::begin(bad_runs), std::end(bad_runs), run_number) != std::end(bad_runs));

            if (!Keep_event || Additional_bad_runs)
                continue;
        }

        if (number_of_electrons != 1||number_of_positrons != 1)
            continue;

        electronAccepted = Accept_electron(electron_pcal_energy, electron_p, electron_pcal_v, electron_pcal_w, electron_vz, version);

        positronAccepted = Accept_positron(positron_ecin_energy, positron_ecout_energy, positron_pcal_energy, positron_p, positron_chi2pid, positron_pcal_v, positron_pcal_w, positron_vz, version);

        if (number_protons != 1||N_FT != 1)
            continue;

        if (Restricted && (number_others > 0 || number_of_piminus > 0 || number_of_piplus > 0))
            continue;

        if (proton_beta > 0.1 && proton_chi2pid < 10 && proton_p > 0.4)
        {
            protonAccepted = true;
        }

        if (electronAccepted && positronAccepted && protonAccepted)
        {
            h_vertex_timediff->Fill(electron_vt - positron_vt); //------PLOT-------
            if (abs(electron_vt - positron_vt) > 1)
                continue;
            //------ENERGY CORRECTION ELECTRON------------
            electron_E_recon = electron_energy + electron_photonE;
            electron_px_recon = electron_px * (electron_E_recon / electron_energy);
            electron_py_recon = electron_py * (electron_E_recon / electron_energy);
            electron_pz_recon = electron_pz * (electron_E_recon / electron_energy);
            TLorentzVector ele(electron_px_recon, electron_py_recon, electron_pz_recon, electron_E_recon);

            //------ENERGY CORRECTION ELECTRON------------
            positron_E_recon = positron_energy + positron_photonE;
            positron_px_recon = positron_px * (positron_E_recon / positron_energy);
            positron_py_recon = positron_py * (positron_E_recon / positron_energy);
            positron_pz_recon = positron_pz * (positron_E_recon / positron_energy);
            TLorentzVector pos(positron_px_recon, positron_py_recon, positron_pz_recon, positron_E_recon);

            h_electron_photons->Fill(electron_photon); //------PLOT-------
            h_positron_photons->Fill(positron_photon); //------PLOT-------
            h_NFTelectron->Fill(N_FT);                 //------PLOT-------

            FT_E_new = FT_E;

            // NOT USE IN PASS2
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
                    auto FT_4vect = algo_correction.Transform(FT_Px, FT_Py, FT_Pz, FT_E);
                    FT_Px_new = std::get<0>(FT_4vect);
                    FT_Py_new = std::get<1>(FT_4vect);
                    FT_Pz_new = std::get<2>(FT_4vect);
                    FT_E_new = std::get<3>(FT_4vect);
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
            h_p->Fill(proton_p);
            double Ppro_ec = proton_pcorr(proton_theta, proton_p);
            h_p_corr->Fill(Ppro_ec);
            double scale_factor = Ppro_ec / proton_p;

            double Ppx_new = proton_px * scale_factor;
            double Ppy_new = proton_py * scale_factor;
            double Ppz_new = proton_pz * scale_factor;

            TLorentzVector pro(Ppx_new, Ppy_new, Ppz_new, proton_energy);

            Double_t SF_corr, m2_corr;

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

            Egamma = E_photon;
            PCAL = positron_sfpcal;
            ECIN = positron_sfecin + SF_corr;
            ECOUT = positron_sfecout;
            m2PCAL = positron_m2pcal;
            m2ECIN = positron_m2ecin * m2_corr;
            m2ECOUT = positron_m2ecout;
            score_pos = readerTMVA->EvaluateMVA("BDT pos method");

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

            score_e = score_ele;
            score_p = score_pos;

            TLorentzVector BeamTarget = (beam + target);
            // INVARIANT MASS
            invariant_ee = (pos + ele);
            // MISSING MASS
            miss = BeamTarget - (pos + ele + eleFT + pro);
            miss_ep = BeamTarget - (eleFT + pro);

            Q2 = 2 * Beam_E * FT_E_new * (1 - cos(eleFT.Theta()));
            E_photon = Beam_E - FT_E_new;
            W = sqrt((0.938 * 0.938) + (2 * 0.938 * E_photon) - Q2);

            h_E_photon->Fill(E_photon); //------PLOT-------

            h_vertex_timediff_FT_FD->Fill(electron_vt - FT_vt);  //------PLOT-------
            h_vertex_timediff_FT_eFD->Fill(positron_vt - FT_vt); //------PLOT-------
            h_vertex_timediff_FT_pFD->Fill(proton_vt - FT_vt);   //------PLOT-------

            Mxep = miss_ep.M();
            MXepee = miss.M2();
            MXeee = (BeamTarget - (eleFT + ele + pos)).M();
            MXepe1 = (BeamTarget - (eleFT + pro + pos)).M2();
            MXepe2 = (BeamTarget - (eleFT + pro + ele)).M2();
            MXee1 = (BeamTarget - (eleFT + pos)).M();
            MXee2 = (BeamTarget - (eleFT + ele)).M();
            MXpe1 = (BeamTarget - (pro + pos)).M();
            MXpe2 = (BeamTarget - (pro + ele)).M();

            if (AI)
            {
                Mxepe = static_cast<Float_t>(MXepe2);
                Mxee = static_cast<Float_t>(MXee2);
                Mxpe = static_cast<Float_t>(MXpe2);
                electronFT_P = eleFT.P();
                electronFT_Theta = eleFT.Theta() * 57.2958;
                electronFT_Phi = eleFT.Phi() * 57.2958;
                protonFD_P = pro.P();
                protonFD_Theta = pro.Theta() * 57.2958;
                protonFD_Phi = pro.Phi() * 57.2958;
                electronFD_P = ele.P();
                electronFD_Theta = ele.Theta() * 57.2958;
                electronFD_Phi = ele.Phi() * 57.2958;

                scoreBDT = reader->EvaluateMVA("BDT method");
                scoreBDTG = reader->EvaluateMVA("BDTG method");
            }

            invariantMass = invariant_ee.M();

            time_ee = electron_vt - FT_vt;
            time_epos = positron_vt - FT_vt;
            time_ep = proton_vt - FT_vt;
            time_efdp = electron_vt - proton_vt;
            time_posfdp = positron_vt - proton_vt;
            results->Fill();

            h_beta_vs_p->Fill(proton_p, proton_beta);

            if (score_ele < 0.05 || score_pos < 0.05)
                continue;

            if (E_photon < 8.1)
                continue;

            h_MXepee->Fill(MXepee);
            if (abs(time_ee) <= 2.0 && abs(MXepee) < 0.5)
            {
                h_Q2->Fill(Q2);
                h_IMee_vs_Mxep->Fill(invariantMass, Mxep);

                h_Invariant->Fill(invariantMass);
                h_Mxep->Fill(Mxep);
                h_MXeee->Fill(MXeee);
                h_MXepe1->Fill(MXepe1);
                h_MXepe2->Fill(MXepe2);
                h_MXee1->Fill(MXee1);
                h_MXee2->Fill(MXee2);
                h_MXpe1->Fill(MXpe1);
                h_MXpe2->Fill(MXpe2);

                h2_Invariant->Fill(invariantMass, invariantMass);
                h2_MXep->Fill(Mxep, invariantMass);
                h2_MXeee->Fill(MXeee, invariantMass);
                h2_MXepe1->Fill(MXepe1, invariantMass);
                h2_MXepe2->Fill(MXepe2, invariantMass);
                h2_MXee1->Fill(MXee1, invariantMass);
                h2_MXee2->Fill(MXee2, invariantMass);
                h2_MXpe1->Fill(MXpe1, invariantMass);
                h2_MXpe2->Fill(MXpe2, invariantMass);

                if (invariantMass >= 3.0 && invariantMass <= 3.2)
                {
                    h_Q2_jpsievents->Fill(Q2);
                    h_W->Fill(W);
                }
            }
        }
        // printf("run = %d\n",fc);
    } // End for "Runs"

    // fclose(f_results);

    algo_correction.Stop();

    file2->Write();

    TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
    TString pdf_original = "./R_others/December2024/" + out_file_name + ".pdf";
    // Fidutial Cuts and vertex time Plots
    // Page 1
    can->Divide(1, 2);
    can->cd(1);
    h_p->Draw();
    can->cd(2);
    h_p_corr->Draw();
    can->Print(pdf_original + "(");
    can->Clear();

    // Page 1
    can->Divide(2, 2);
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
    // Page 2
    can->Divide(2, 2);
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
    // Page 3
    can->Divide(1, 4);
    can->cd(1);
    h_NFTelectron->Draw();
    can->cd(2);
    gPad->SetLogy();
    h_vertex_timediff_FT_FD->Draw();
    can->cd(3);
    gPad->SetLogy();
    h_vertex_timediff_FT_eFD->Draw();
    can->cd(4);
    gPad->SetLogy();
    h_vertex_timediff_FT_pFD->Draw();
    can->Print(pdf_original + "(");
    can->Clear();
    gPad->SetLogy(0);
    // Variables
    // Page 4
    can->Divide(2, 2);
    can->cd(1);
    h_E_photon->Draw();
    can->cd(2);
    h_W->Draw();
    can->cd(3);
    h_Q2->Draw();
    can->cd(4);
    h_Q2_jpsievents->Draw();
    can->Print(pdf_original + "(");
    can->Clear();
    // Missing mass cut, Invariant mass
    // Page 5
    h_MXepee->Draw();
    can->Print(pdf_original + "(");
    h_IMee_vs_Mxep->Draw("colz");
    can->Print(pdf_original + "(");
    can->Clear();
    // Missing mass
    // Page 6
    can->Divide(3, 3);
    can->cd(1);
    h_Invariant->Draw();
    can->cd(2);
    h_Mxep->Draw();
    can->cd(3);
    h_MXeee->Draw();
    can->cd(4);
    h_MXepe1->Draw();
    can->cd(5);
    h_MXepe2->Draw();
    can->cd(6);
    h_MXee1->Draw();
    can->cd(7);
    h_MXee2->Draw();
    can->cd(8);
    h_MXpe1->Draw();
    can->cd(9);
    h_MXpe2->Draw();
    can->Print(pdf_original + "(");
    can->Clear();
    // Page 7
    can->Divide(3, 3);
    can->cd(1);
    h2_Invariant->Draw("colz");
    can->cd(2);
    h2_MXep->Draw("colz");
    can->cd(3);
    h2_MXeee->Draw("colz");
    can->cd(4);
    h2_MXepe1->Draw("colz");
    can->cd(5);
    h2_MXepe2->Draw("colz");
    can->cd(6);
    h2_MXee1->Draw("colz");
    can->cd(7);
    h2_MXee2->Draw("colz");
    can->cd(8);
    h2_MXpe1->Draw("colz");
    can->cd(9);
    h2_MXpe2->Draw("colz");
    can->Print(pdf_original + ")");

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;
}
