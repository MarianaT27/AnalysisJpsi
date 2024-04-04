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

void upd() { 

	//Electron variables
    Double_t electron_p, electron_theta, electron_phi;
    Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
    Double_t electron_sfpcal, electron_sfecin, electron_sfecout;

    //Positron variables
    Double_t positron_p, positron_theta, positron_phi;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;

    cout<<"Updating tree"<<endl;


	TFile *f = new TFile("/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/S19.root","update"); 
	TTree *tree = (TTree*)f->Get("analysis");
	TMVA::Reader *readerTMVA = new TMVA::Reader( "!Color:!Silent" );
	int model=9;
    // Create a set of variables and declare them to the reader
	Float_t P, Theta, Phi, PCAL,ECIN,ECOUT;

	Double_t score_pos,score_ele;
	Float_t m2PCAL=-1;
	Float_t m2ECIN=-1;
	Float_t m2ECOUT=-1;
	Float_t Nphe;


	readerTMVA->AddVariable( "P",&P );
	readerTMVA->AddVariable( "Theta",&Theta);
	readerTMVA->AddVariable( "Phi",&Phi);
     //readerTMVA->AddVariable( "Nphe",&Nphe);
	readerTMVA->AddVariable( "SFPCAL",&PCAL);
	readerTMVA->AddVariable( "SFECIN",&ECIN);
	readerTMVA->AddVariable( "SFECOUT",&ECOUT );
	readerTMVA->AddVariable( "m2PCAL",&m2PCAL);
	readerTMVA->AddVariable( "m2ECIN",&m2ECIN);
	readerTMVA->AddVariable( "m2ECOUT",&m2ECOUT);

    //Book Methods
	TString weightfile_pos; 
	TString weightfile_ele;

	weightfile_ele= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19neg/TMVAClassification_BDT.weights.xml";
	weightfile_pos= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19pos/TMVAClassification_BDT.weights.xml";

	readerTMVA->BookMVA( "BDT pos method", weightfile_pos );
	readerTMVA->BookMVA( "BDT ele method", weightfile_ele );

	TBranch *score_p = tree->Branch("score_pos",&score_pos,"score_pos/d");
	TBranch *score_e = tree->Branch("score_ele",&score_ele,"score_ele/d");

	cout<<"Adding new branches"<<endl;

	tree->SetBranchAddress("electron_p",&electron_p);
	tree->SetBranchAddress("electron_theta",&electron_theta);
	tree->SetBranchAddress("electron_phi",&electron_phi);
	tree->SetBranchAddress("electron_m2ecin",&electron_m2ecin);
	tree->SetBranchAddress("electron_m2ecout",&electron_m2ecout);
	tree->SetBranchAddress("electron_m2pcal",&electron_m2pcal);
	tree->SetBranchAddress("electron_sfecin",&electron_sfecin);
	tree->SetBranchAddress("electron_sfpcal",&electron_sfpcal);
	tree->SetBranchAddress("electron_sfecout",&electron_sfecout);

	tree->SetBranchAddress("positron_p",&positron_p);
	tree->SetBranchAddress("positron_theta",&positron_theta);
	tree->SetBranchAddress("positron_phi",&positron_phi);
	tree->SetBranchAddress("positron_m2ecin",&positron_m2ecin);
	tree->SetBranchAddress("positron_m2ecout",&positron_m2ecout);
	tree->SetBranchAddress("positron_m2pcal",&positron_m2pcal);
	tree->SetBranchAddress("positron_sfecin",&positron_sfecin);
	tree->SetBranchAddress("positron_sfpcal",&positron_sfpcal);
	tree->SetBranchAddress("positron_sfecout",&positron_sfecout);


	Long64_t nentries = tree->GetEntries(); 
	cout<<"Starting loop"<<endl;
	for (Long64_t i=0;i<nentries;i++) {
		tree->GetEntry(i); 
		P=positron_p;
		Theta=positron_theta/57.2958;
		Phi=positron_phi/57.2958;
		PCAL=positron_sfpcal;
		ECIN=positron_sfecin;
		ECOUT=positron_sfecout;
		m2PCAL=positron_m2pcal;
		m2ECIN=positron_m2ecin;
		m2ECOUT=positron_m2ecout;
        score_pos=readerTMVA->EvaluateMVA("BDT pos method");

		P=electron_p;
		Theta=electron_theta/57.2958;
		Phi=electron_phi/57.2958;
		PCAL=electron_sfpcal;
		ECIN=electron_sfecin;
		ECOUT=electron_sfecout;
		m2PCAL=electron_m2pcal;
		m2ECIN=electron_m2ecin;
		m2ECOUT=electron_m2ecout;
        score_ele=readerTMVA->EvaluateMVA("BDT ele method");

		score_p->Fill(); 
		score_e->Fill();
	} 
	tree->Print(); 
	tree->Write(); 
	delete f; 
}
