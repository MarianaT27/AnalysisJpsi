

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <TLorentzVector.h>
#include <cmath>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
namespace fs = std::filesystem;

int evaluate_TMVA(TString outname = "Var9_v1", TString model = "Var9_v1", TString inputFile="F18in_Full_epe-_12") {
    TMVA::Tools::Instance();

    std::cout << std::endl;
    std::cout << "==> Start TMVA_Tagged_Test" << std::endl;

    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

    // Variables declaration
    // Variables to read
    Float_t Mxep, Mxepe, Mxee, Mxpe;
    Float_t electronFT_P, electronFT_Theta;
    Float_t protonFD_P, protonFD_Theta;
    Float_t electronFD_P, electronFD_Theta;
    Float_t D_Phi2,D_Phi1,D_Phi3;
    Float_t electronFT_Phi,electronFD_Phi,protonFD_Phi;
    Double_t scoreBDT,scoreBDTG,scoreMLP;

    reader->AddVariable("electronFT_P", &electronFT_P);
    reader->AddVariable("electronFT_Theta", &electronFT_Theta);
    reader->AddVariable("electronFD_P", &electronFD_P);
    reader->AddVariable("electronFD_Theta", &electronFD_Theta);
    reader->AddVariable("D_Phi3", &D_Phi3);
    reader->AddVariable("protonFD_P", &protonFD_P);
    reader->AddVariable("protonFD_Theta", &protonFD_Theta);
    if(model=="Var10")
        reader->AddVariable("Mxepe", &Mxepe);
    reader->AddVariable("Mxpe", &Mxpe);
    reader->AddVariable("Mxee", &Mxee);
    


    
    TString dir = "/volatile/clas12/mtenorio/Weights_AI/TMVA_Training_" + model + "/weights/";
    TString prefix = "TMVAClassification";

    // Book MVA methods
    reader->BookMVA("BDT method", dir + prefix + "_BDT.weights.xml");
    reader->BookMVA("BDTG method", dir + prefix + "_BDTG.weights.xml");
    
    bool MLP=false;
    std::string filename = (dir + prefix + "_MLP.weights.xml").Data();;
    if (fs::exists(filename)) {
        MLP=true;
    } 

    if(MLP){
        reader->BookMVA("MLP method", dir + prefix + "_MLP.weights.xml");
    }

    // Input ROOT file
    TString file1 = "/volatile/clas12/mtenorio/Root/"+inputFile+".root";
    TFile *input1 = new TFile(file1);
    std::cout << "--- TMVA evaluation : Using input file: " << input1->GetName() << std::endl;

    TTree *tree = (TTree *)input1->Get("results");

    // Declare variables for existing tree branches
    
    tree->SetBranchAddress("Mxepe", &Mxepe);
    tree->SetBranchAddress("Mxee", &Mxee);
    tree->SetBranchAddress("Mxpe", &Mxpe);
    tree->SetBranchAddress("Mxep", &Mxep);
    tree->SetBranchAddress("electronFT_P", &electronFT_P);
    tree->SetBranchAddress("electronFT_Theta", &electronFT_Theta);
    tree->SetBranchAddress("electronFT_Phi", &electronFT_Phi);
    tree->SetBranchAddress("protonFD_P", &protonFD_P);
    tree->SetBranchAddress("protonFD_Theta", &protonFD_Theta);
    tree->SetBranchAddress("protonFD_Phi", &protonFD_Phi);
    tree->SetBranchAddress("electronFD_P", &electronFD_P);
    tree->SetBranchAddress("electronFD_Theta", &electronFD_Theta);
    tree->SetBranchAddress("electronFD_Phi", &electronFD_Phi);
    tree->SetBranchAddress("D_Phi1", &D_Phi1);
    tree->SetBranchAddress("D_Phi2", &D_Phi2);
    tree->SetBranchAddress("D_Phi3", &D_Phi3);

    // Create a new file to save the tree with additional branches
    TFile *outputFile = new TFile("/volatile/clas12/mtenorio/Root/TMVA_test/"+outname+".root", "RECREATE");
    TTree *newTree = new TTree("results", "/volatile/clas12/mtenorio/Root/TMVA_test/"+outname+".root");

    // Add new branches for scores
    newTree->Branch("electronFT_P", &electronFT_P, "electronFT_P/F");
    newTree->Branch("electronFT_Theta", &electronFT_Theta, "electronFT_Theta/F");
    
    newTree->Branch("protonFD_P", &protonFD_P, "protonFD_P/F");
    newTree->Branch("protonFD_Theta", &protonFD_Theta, "protonFD_Theta/F");
    
    newTree->Branch("electronFD_P", &electronFD_P, "electronFD_P/F");
    newTree->Branch("electronFD_Theta", &electronFD_Theta, "electronFD_Theta/F");

    newTree->Branch("D_Phi1", &D_Phi1, "D_Phi1/F");
    newTree->Branch("D_Phi2", &D_Phi2, "D_Phi2/F");
    newTree->Branch("D_Phi3", &D_Phi3, "D_Phi3/F");

    newTree->Branch("electronFD_Phi", &electronFD_Phi, "electronFD_Phi/F");
    newTree->Branch("electronFT_Phi", &electronFT_Phi, "electronFT_Phi/F");
    newTree->Branch("protonFD_Phi", &protonFD_Phi, "protonFD_Phi/F");
    

    newTree->Branch("Mxepe", &Mxepe, "Mxepe/F");
    newTree->Branch("Mxee", &Mxee, "Mxee/F");
    newTree->Branch("Mxpe", &Mxpe, "Mxpe/F");

    newTree->Branch("Mxep", &Mxep, "Mxep/F");
    newTree->Branch("scoreBDT", &scoreBDT, "scoreBDT/D");
    newTree->Branch("scoreBDTG", &scoreBDTG, "scoreBDTG/D");
    newTree->Branch("scoreMLP", &scoreMLP, "scoreMLP/D");

    std::cout << "--- Processing: " << tree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();
    int percentageStep = 5;
    int step = tree->GetEntries() * percentageStep / 100;
    scoreBDT=-99;
    scoreBDTG=-99;
    scoreMLP=-99;

    for (Long64_t ievt = 0; ievt < tree->GetEntries(); ievt++) {
        if (ievt % step == 0){
            double percentage = (ievt * 100.0) / tree->GetEntries();
            std::cout << "Progress: " << percentage << std::endl;
        }
        tree->GetEntry(ievt);
        // Compute MVA scores
        scoreBDT = reader->EvaluateMVA("BDT method");
        scoreBDTG = reader->EvaluateMVA("BDTG method");
        if(MLP){
            scoreMLP= reader->EvaluateMVA("MLP method");
        }

        // Fill the new tree
        newTree->Fill();
    }

    // Save the new tree and close files
    newTree->Write();
    outputFile->Close();
    input1->Close();

    sw.Stop();
    std::cout << "--- End of event loop: ";
    sw.Print();

    std::cout << "--- Created root file: "<<outname<<".root containing the MVA output scores" << std::endl;

    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
    return 0;
}


