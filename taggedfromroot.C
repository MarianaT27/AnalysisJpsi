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

int taggedfromroot(string nameFile = "S19", int version = -19, int PASS = 2, bool Restricted=false, bool QA = false)
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    // Set Beam energy
    double Beam_E;

    if (abs(version) == 18)
        Beam_E = 10.6;
    else if(abs(version)==19)
        Beam_E = 10.2;
    else if(abs(version)==10)
        Beam_E = 10.6;
    else   
        Beam_E = 6.4; 

        cout<<"Beam Energy is :"<<Beam_E<<endl;

    
    //********************************
    // DECLARATION OF VARIABLES
    //********************************
    int max;
    // MISSING MOMENTUM AND INVARIANT MASS
    TLorentzVector miss;
    TLorentzVector invariant_ee;
    // FT electron vector
    TLorentzVector electronFT_vec;
    TLorentzVector electronFD_vec;
    TLorentzVector positronFD_vec;
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
    Double_t FT_vtFcal, FT_xFcal, FT_yFcal, FT_zFcal, RFcal,beamCharge_mC;

    Float_t score;

    int N_FT;
    int number_of_electrons;
    int number_of_positrons;
    int JPSI = 0;
    bool electronAccepted, positronAccepted;

    int number_of_piplus = 0;
    int number_of_piminus = 0;
    int number_of_gamma = 0;
    int number_others = 0;

    Double_t score_pos, score_ele;
    Double_t score_pos_6, score_ele_6, score_e_6, score_p_6;

    // EVENT SELECTION
    TH2F *h_electron_ECin_vs_PCAL = new TH2F("h_electron_ECin_vs_PCAL", "Electron SF-ECin vs. SF-PCAL; SF_{PCAL};SF_{ECIN}", 200, 0.0, 0.3, 100, 0.0, 0.2);
    TH2F *h_positron_EC_vs_PCAL = new TH2F("h_positron_EC_vs_PCAL", "Positrons SF-EC vs. SF-PCALSF_{PCAL};SF_{EC}", 200, 0.0, 0.3, 100, 0.0, 0.2);
    TH1F *h_electron_zvertex = new TH1F("h_electron_zvertex", "z-vertex distribution lepton;z-vertex,cm", 200, -15, 35);
    TH1F *h_vertex_timediff = new TH1F("h_vertex_timediff", "vertex time difference e- e+ ; vt_e^- - vt_e^+, ns", 120, -2, 2);

    // Photon energy correction
    TH1I *h_electron_photons = new TH1I("h_electron_photons", "Photons radiated from electrons", 7, 0, 7);
    TH1I *h_positron_photons = new TH1I("h_positron_photons", "Photons radiated from positrons", 7, 0, 7);

    // TAGGED
    TH1F *h_NFTelectron = new TH1F("h_NFTelectron", "# FT electrons; # electrons in FT", 10, 0, 10);
    TH1F *h_vertex_timediff_FT_FD = new TH1F("h_vertex_timediff_FT_FD", "Vertex time difference between FT_e^- and FD_e^-; #Delta t_v=FD_e^- - FT_e^-, ns;", 350, -25, 25);

    // KINEMATICS
    TH1F *h_E_photon = new TH1F("h_E_photon", "Photon Energy", 120, 0, 0);
    TH1F *h_Q2 = new TH1F("h_Q2", "Q^{2};Q^{2};Counts", 120, 0, 0);
    TH1F *h_Q2_jpsievents = new TH1F("h_Q2_jpsievents", "Q^{2};Q^{2};Counts", 120, 0, 0);
    TH1F *h_W = new TH1F("h_W", "Hadronic Mass;W;Counts", 120, 4, 4.6);
    TH1F *h_MM = new TH1F("h_MM", "Missing Mass; MM; Counts ", 90, 0, 2);

    TH2F *h_mm_vs_invariantmass_t = new TH2F("h_mm_vs_invariantmass_t", "MM vs M(e^+e^-) for t+-2;M(e^+e^-), GeV;Missing Mass, GeV", 150, 0, 4, 150, 0, 3);
    TH2F *h_mm_vs_invariantmass_no_t = new TH2F("h_mm_vs_invariantmass_no_t", "MM vs M(e^+e^-) outside time cut;M(e^+e^-), GeV;Missing Mass, GeV", 150, 0, 4, 150, 0, 3);

    // RESULTS RESULTS
    TH1F *h_Invariant = new TH1F("h_im", ";M(e+e-),GeV ", 90, 2.0, 4);
    TH1F *h_Invariant_not = new TH1F("h_Invariant_not", "Outside time range;M(e+e-),GeV", 90, 2.0, 4);


    ////weitghs /lustre19/expphy/volatile/clas12/mtenorio/weights

    //TString root_file = "/w/hallb-scshelf2102/clas12/mtenorio/Analysis/" + nameFile + ".root";
     TString root_file = "/volatile/clas12/mtenorio/Root/S18/"+nameFile+".root";
     TChain *tree = new TChain("analysis"); // "analysis" is the tree name in the files
    
     if(nameFile=="S18_10GeV_In")
        root_file = "/volatile/clas12/mtenorio/Root/S18/"+nameFile+"1.root";

    tree->Add(root_file);
    if(nameFile=="S18_10GeV_In"){
        root_file = "/volatile/clas12/mtenorio/Root/S18/"+nameFile+"2.root";
        tree->Add(root_file);
    }

    tree->SetMakeClass(1);

    //tree->SetBranchAddress("Beam_E", &Beam_E);

    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchAddress("event_number", &event_number);
    tree->SetBranchAddress("beamCharge_mC", &beamCharge_mC);


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

    tree->SetBranchAddress("electron_photon", &electron_photon);
    tree->SetBranchAddress("positron_photon", &positron_photon);
    tree->SetBranchAddress("electron_photonE", &electron_photonE);
    tree->SetBranchAddress("positron_photonE", &positron_photonE);

    tree->SetBranchAddress("number_of_piplus", &number_of_piplus);
    tree->SetBranchAddress("number_of_piminus", &number_of_piminus);
    tree->SetBranchAddress("number_of_gamma", &number_of_gamma);
    tree->SetBranchAddress("number_others", &number_others);

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
    // tree->SetBranchAddress("electron_sec",&electron_sec);
    // tree->SetBranchAddress("positron_sec",&positron_sec);
    // tree->SetBranchAddress("electronFT_sec",&electronFT_sec);
    TString root_file1, root_file2;
    bool useTwoFiles = false;
    
    if (nameFile == "S18_10GeV_In") {
        root_file1 = "/volatile/clas12/mtenorio/Root/S18/S18_10GeV_In1.root";
        root_file2 = "/volatile/clas12/mtenorio/Root/S18/S18_10GeV_In2.root";
        cout<<"Using 2 files"<<endl;
        useTwoFiles = true;
    } else {
        root_file1 = root_file;
    }
    
    // Open first file and retrieve histograms
    TFile *file1 = new TFile(root_file1, "READ");
    if (!file1 || file1->IsZombie()) {
        cout << "Error: Could not open file " << root_file1 << std::endl;
       return 0;
    }


    TH2F *h1_positron = (TH2F *)file1->Get("h_delta_theta_vs_phi_positron");
    TH2F *h1_electron = (TH2F *)file1->Get("h_delta_theta_vs_phi_electron");
    
    if (!h1_positron || !h1_electron) {
        cout << "Error: One or more histograms missing in " << root_file1 << std::endl;
        file1->Close();
       return 0;
    }
    if (!h1_electron) {
        std::cerr << "Error: h1_electron is NULL! Histogram missing in " << root_file1 << std::endl;
        file1->ls();  // List objects in the ROOT file
        file1->Close();
        return 0;
    } else {
        std::cout << "Successfully retrieved h1_electron from " << root_file1 << std::endl;
    }
    gROOT->cd(); 
    // Clone the first histogram to create a combined histogram
    TH2F *h_delta_theta_vs_phi_positron = (TH2F *)h1_positron->Clone("h_delta_theta_vs_phi_positron_combined");
    TH2F *h_delta_theta_vs_phi_electron = (TH2F *)h1_electron->Clone("h_delta_theta_vs_phi_electron_combined");
    if (!h_delta_theta_vs_phi_positron || !h_delta_theta_vs_phi_electron) {
        cout << "Error: One or more histograms missing in " << root_file1 << std::endl;
        file1->Close();
       return 0;
    }


    if (useTwoFiles) {
        // Open second file
        TFile *file2 = new TFile(root_file2, "READ");
        if (!file2 || file2->IsZombie()) {
            cout << "Error: Could not open file " << root_file2 << std::endl;
            file1->Close();
           return 0;
        }
    
        TH2F *h2_positron = (TH2F *)file2->Get("h_delta_theta_vs_phi_positron");
        TH2F *h2_electron = (TH2F *)file2->Get("h_delta_theta_vs_phi_electron");

        if (!h2_positron || !h2_electron) {
            cout << "Error: One or more histograms missing in " << root_file2 << std::endl;
        } else {
            // Add histograms from the second file
            cout<<"adding histograms from "<<root_file2<<endl;
            h_delta_theta_vs_phi_positron->Add(h2_positron);
            h_delta_theta_vs_phi_electron->Add(h2_electron);
        }
    
        file2->Close();
    }
    
    // Close first file
    file1->Close();

    if (!h_delta_theta_vs_phi_electron) {
        std::cerr << "Error: h_delta_theta_vs_phi_electron is NULL!" << std::endl;
    } else {
        std::cout << "Histogram successfully created." << std::endl;
    }


    //--------SET OUT ROOT file---------

    TString out_file_name = Name_file(nameFile, 1);
    TString out_file = "/lustre24/expphy/volatile/clas12/mtenorio/Root/S18/" + out_file_name + ".root";

    TFile *file2 = new TFile(out_file, "RECREATE");
    TTree *results = new TTree("results", out_file);

    Double_t invariantMass, MM, Egamma, score_p, score_e, time_ee;
    results->Branch("beamCharge_mC", &beamCharge_mC, "beamCharge_mC/d");
    results->Branch("electron_p", &electron_p, "electron_p/d");
    results->Branch("positron_p", &positron_p, "positron_p/d");
    results->Branch("invariantMass", &invariantMass, "invariantMass/d");
    results->Branch("Q2", &Q2, "Q2/d");
    results->Branch("MM", &MM, "MM/d");
    results->Branch("time_ee", &time_ee, "time_ee/d");
    results->Branch("Egamma", &Egamma, "Egamma/d");
    results->Branch("W", &W, "W/d");

  
    //TMVA::Reader *readerTMVA = new TMVA::Reader("!Color:Silent");

    int model = 6;
    // Create a set of variables and declare them to the reader
    Float_t P, Theta, Phi, PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;
    Float_t Nphe;
    /*
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
    */

    // ANALYSIS
    /////////Instanciate QADB///////////
    //QADB *qa = new QADB();
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

        tree->GetEntry(fc);

        if (fc % step == 0)
        {
            double percentage = (fc * 100.0) / tree->GetEntries();
            std::cout << "Progress: " << percentage << "%" << std::endl;
        }


        if (number_of_electrons != 1 || number_of_positrons != 1)
            continue;
       
       
        if(version<0)
            h_electron_zvertex->Fill(electron_vz);
        else
            h_electron_zvertex->Fill(positron_vz);
    

        electronAccepted = Accept_electron(electron_pcal_energy, electron_p, electron_pcal_v, electron_pcal_w, electron_vz, version);

        positronAccepted = Accept_positron(positron_ecin_energy, positron_ecout_energy, positron_pcal_energy, positron_p, positron_chi2pid, positron_pcal_v, positron_pcal_w, positron_vz, version);

        
        if (Restricted && (number_others > 0 || number_of_piminus > 0 || number_of_piplus > 0))
            continue;

        if (electronAccepted && positronAccepted)
        {
            h_electron_ECin_vs_PCAL->Fill(electron_sfpcal,electron_sfecin);
            h_positron_EC_vs_PCAL->Fill(positron_sfpcal,(positron_ecin_energy+positron_ecout_energy)/positron_p);
            h_vertex_timediff->Fill(electron_vt - positron_vt); //------PLOT-------
            if (abs(electron_vt - positron_vt) <= 1)
            {
                //------ENERGY CORRECTION ELECTRON------------
                electron_E_recon = electron_energy + electron_photonE;
                electron_px_recon = electron_px * (electron_E_recon / electron_energy);
                electron_py_recon = electron_py * (electron_E_recon / electron_energy);
                electron_pz_recon = electron_pz * (electron_E_recon / electron_energy);
                electronFD_vec.SetPxPyPzE(electron_px_recon, electron_py_recon, electron_pz_recon, electron_E_recon);

                //------ENERGY CORRECTION ELECTRON------------
                positron_E_recon = positron_energy + positron_photonE;
                positron_px_recon = positron_px * (positron_E_recon / positron_energy);
                positron_py_recon = positron_py * (positron_E_recon / positron_energy);
                positron_pz_recon = positron_pz * (positron_E_recon / positron_energy);
                positronFD_vec.SetPxPyPzE(positron_px_recon, positron_py_recon, positron_pz_recon, positron_E_recon);

                h_electron_photons->Fill(electron_photon); //------PLOT-------
                h_positron_photons->Fill(positron_photon); //------PLOT-------
                h_NFTelectron->Fill(N_FT);                 //------PLOT-------

                TLorentzVector eleFT(FT_Px, FT_Py, FT_Pz, FT_E);

                if (N_FT == 1)
                {

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
                    electronFT_vec.SetPxPyPzE(FT_Px_new, FT_Py_new, FT_Pz_new, FT_E_new);

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
                    //score_pos = readerTMVA->EvaluateMVA("BDT pos method");

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
                    //score_ele = readerTMVA->EvaluateMVA("BDT ele method");

                    //score_e = score_ele;
                    //score_p = score_pos;

                    // MISSING MASS & INVARIANT MASS
                    miss = (beam + target) - (positronFD_vec + electronFD_vec + electronFT_vec);
                    invariant_ee = (positronFD_vec + electronFD_vec);

                    Q2 = 2 * Beam_E * FT_E_new * (1 - cos(electronFT_vec.Theta()));
                    E_photon = Beam_E - FT_E_new;
                    W = sqrt((0.938 * 0.938) + (2 * 0.938 * E_photon) - Q2);
                    invariantMass = invariant_ee.M();
                    MM = miss.M();
                    Egamma = E_photon;

                    time_ee = electron_vt - FT_vt;
                    results->Fill();

                    h_E_photon->Fill(E_photon); //------PLOT-------
                    h_vertex_timediff_FT_FD->Fill(time_ee); //------PLOT-------


                    //if (score_ele < 0.05 || score_pos < 0.05)
                      //  continue;
                    if (E_photon < 8.1)
                        continue;

                    // MM vs IM
                    if (abs(time_ee) <= 2.5)
                    {
                        h_mm_vs_invariantmass_t->Fill(invariant_ee.M(), miss.M());
                        h_MM->Fill(miss.M());
                        if (0.7 <= miss.M() && miss.M() <= 1.3)
                        { // 1 sigmas 0.65 <=miss.M() &&  miss.M()<= 1.3
                            h_Invariant->Fill(invariant_ee.M());
                            h_Q2->Fill(Q2); //------PLOT-------
                            if (invariant_ee.M() >= 3.0 && invariant_ee.M() <= 3.2)
                            {
                                h_Q2_jpsievents->Fill(Q2);
                                h_W->Fill(W);
                            }
                        }
                    }
                    else
                    {
                        h_mm_vs_invariantmass_no_t->Fill(invariant_ee.M(), miss.M());
                        h_Invariant_not->Fill(invariant_ee.M());
                    }

                } // N_FT==1
            }
        }
        // printf("run = %d\n",fc);
    } // End for "Runs"

    // fclose(f_results);

    file2->Write();

    TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
    TString pdf_original = "./R_Tagged/" + out_file_name + ".pdf";

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
    can->Divide(1, 2);
    can->cd(1);
    // gPad->SetLogy();
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
    can->Print(pdf_original + "(");

    h_Q2->Draw();
    can->Print(pdf_original + "(");

    h_Q2_jpsievents->Draw();
    can->Print(pdf_original + ")");

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;
}
