#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "THStack.h"
#include <iostream>
#include <fstream>
using namespace std;

int Variable_Plots_Tagged(TString name_pdf="10_vs_12_distributions",TString name_file1="F18in_Full_epe-_10",Double_t score_cut=0.0,TString name_file2="Var9_accidentals", TString BGK_name="BGK_1")
{
	gROOT->SetBatch(kTRUE);

	gStyle->SetPaintTextFormat("4.1f");
	gStyle->SetPalette(kLightTemperature);
	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(.045, "xyz");
	gStyle->SetTitleSize(.050, "xyz");
	gStyle->SetTitleSize(.07, "t");
	gStyle->SetFrameLineWidth(1);
	gStyle->SetLineWidth(1);
	gStyle->SetHistLineWidth(1);
	gStyle->SetMarkerStyle(13);
	gStyle->SetTitleW(0.8); // per cent of the pad width
	gStyle->SetTitleH(0.1); // per cent of the pad height

	//TString name = "Tagged_12";

	TFile *Data1_file = new TFile("/lustre24/expphy/volatile/clas12/mtenorio/Root/TMVA_test/"+name_file1+".root");
	TTree *Data1_tree = (TTree *)Data1_file->Get("results");
	//int entries1_sig = Data1_tree->GetEntries();
	TString label1_sig = name_file1;
	TCut cut_0 = "";
	

	TFile *Data2_file = new TFile("/lustre24/expphy/volatile/clas12/mtenorio/Root/TMVA_test/Var9_0-100.root");
	//TFile *Data2_file = new TFile("/lustre24/expphy/volatile/clas12/mtenorio/Root/SIG.root");
	TTree *Data2_tree = (TTree *)Data2_file->Get("results");
	//int entries2_sig = Data2_tree->GetEntries();
	TString label2_sig = "epe-";
	TCut cut_1 = "";

	TFile *Data3_file = new TFile("/lustre24/expphy/volatile/clas12/mtenorio/Root/SIG_Processed.root");
	TTree *Data3_tree = (TTree *)Data3_file->Get("results");
	//int entries2_sig = Data2_tree->GetEntries();
	TString label3_sig = "MC Signal";
	TCut cut_2 = "";

	TFile *Data4_file = new TFile("/lustre24/expphy/volatile/clas12/mtenorio/Root/"+BGK_name+".root");
	TTree *Data4_tree = (TTree *)Data4_file->Get("results");
	//int entries2_sig = Data2_tree->GetEntries();
	TString label4_sig = "BGK";
	

	


	std::vector<std::vector<TString>> labels1D{

		// ElectronFT
		{"electronFT_P", "P for electronFT", "custom_range", "0.", "2.5",
		 "custom_binning", "100", "", "legend", "P_electronFT"},

		{"electronFT_Theta", "#theta for electronFT", "custom_range", "2.5", "5.5",
		 "custom_binning", "100", "", "legend", "Theta_electronFT"},

		{"D_Phi1", "#Delta #phi = #phi eFT-#phi eFD", "custom_range", "0.", "360.",
		 "custom_binning", "100", "", "legend", "D_Phi1"},

		// ElectronFD Variables
		{"electronFD_P", "P for electronFD", "custom_range", "2.", "7.",
		 "custom_binning", "100", "", "legend", "P_electronFD"},

		{"electronFD_Theta", "#theta for electronFD", "custom_range", "10.", "37.",
		 "custom_binning", "100", "", "legend", "Theta_electronFD"},

		{"D_Phi2", "#Delta #phi = #phi eFT-#phi pFD", "custom_range", "0.", "360.",
		 "custom_binning", "100", "", "legend", "D_Phi2"},

		// ProtonFD Variables
		{"protonFD_P", "P for protonFD", "custom_range", "0.5", "4.",
		 "custom_binning", "100", "", "legend", "P_protonFD"},

		{"protonFD_Theta", "#theta for protonFD", "custom_range", "5.", "25.",
		 "custom_binning", "100", "", "legend", "Theta_protonFD"},

		{"D_Phi3", "#Delta #phi = #phi eFD-#phi pFD", "custom_range", "0.", "360.",
		 "custom_binning", "100", "", "legend", "D_Phi3"},

		// MissingFD
		{"Mxepe", "M_{X}(e'p'e^{-})", "custom_range", "-1.", "1.",
		 "custom_binning", "100", "", "legend", "Mxepe"},

		{"Mxee", "M_{X}(e'e^{-})", "custom_range", "2.", "8.",
		 "custom_binning", "100", "", "legend", "Mxee"},

		{"Mxpe", "M_{X}(p'e^{-})", "custom_range", "0.", "2.",
		 "custom_binning", "100", "", "legend", "Mxpe"}};


	std::vector<std::vector<TString>> Invariant{
		{"Mxep", "M_{X}(e'p')", "custom_range", "2.0", "3.5",
		 "custom_binning", "100", "", "legend", "Mxep"}};

	
	/*
	TCanvas *cancG0 = new TCanvas("", "can0", 1200, 1000);
	cancG0->Divide(3,4,0.01,0.01,0);
	for (int i = 0; i < labels1D.size(); i++)
	{
		cancG0->cd(i+1);
		//////////////////////////////////////////////////
		// Set options for each label
		TString label1 = labels1D[i][0];
		TString xAxis_label = labels1D[i][1];
		TString range_x_option = labels1D[i][2];
		TString min_x_option = labels1D[i][3];
		TString max_x_option = labels1D[i][4];
		TString binning_x_option = labels1D[i][5];
		TString nb_x_bins = labels1D[i][6];
		TString string_cut = labels1D[i][7];
		TString legend_option = labels1D[i][8];
		TString output_string = labels1D[i][9];
		//////////////////////////////////////////////////

		cout << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << "Doing 1D " << label1 << " plot" << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << endl;

		TCut cut = string_cut.Data();

		TString cut_string = cut.GetTitle();

		float min_X_histo_ini;
		float max_X_histo_ini;

		int nBins_X = 20;

		if (range_x_option == "custom_range")
		{
			//cout << "here " << endl;
			min_X_histo_ini = stof((string)min_x_option.Data());
			max_X_histo_ini = stof((string)max_x_option.Data());
		}

		if (binning_x_option == "custom_binning")
		{
			//cout << "here " << nb_x_bins.Data() << endl;
			nBins_X = stoi((string)nb_x_bins.Data());
		}

		TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		Data1_tree->Draw(label1 + ">>Data_hist",cut_0 * cut);
		// Data1_tree->Draw("Q2>>Data_hist");
		// Styling for Signal (Blue Solid Fill)
		Data_hist->SetFillColorAlpha( kCyan-10  ,0.50);
		Data_hist->SetFillStyle(1001); // Solid fill
		Data_hist->SetLineColor(kCyan-8);
		Data_hist->SetLineWidth(2);
		Data_hist->SetTitle(";" + xAxis_label+";");
		
		
		TH1F *BGK_hist = new TH1F("BGK_hist", "BGK_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		//TCut weight = Form("%i/%i", entries_sig, entries_bgk);

		Data2_tree->Draw(label1 + ">>BGK_hist",cut_0 * cut);


		

		// Styling for Background (Red Hatched Fill)
		BGK_hist->SetFillColor(kRed-7);
		BGK_hist->SetFillStyle(3004); // Hatched pattern
		BGK_hist->SetLineColor(kRed-7);
		BGK_hist->SetLineWidth(2);
        BGK_hist->SetTitle(";" + xAxis_label+";");



		TH1F *SIG_hist = new TH1F("SIG_hist", "SIG_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		//TCut weight = Form("%i/%i", entries_sig, entries_bgk);

		Data3_tree->Draw(label1 + ">>SIG_hist",cut_0 * cut);
		SIG_hist->SetFillColor(kGreen+3);
		SIG_hist->SetFillStyle(3003); // DOt pattern
		SIG_hist->SetLineColor(kGreen+3);
		SIG_hist->SetLineWidth(2);
        SIG_hist->SetTitle(";" + xAxis_label+";");

		
		BGK_hist->Scale(1./BGK_hist->Integral(),"nosw2");
		Data_hist->Scale(1./Data_hist->Integral(),"nosw2");
		SIG_hist->Scale(1./SIG_hist->Integral(),"nosw2");
		//gStyle->SetLegendTextSize(0.06);

		Data_hist->SetMaximum(std::max({Data_hist->GetMaximum(), SIG_hist->GetMaximum(),BGK_hist->GetMaximum()}) * 1.2);

		int hist = Data_hist->GetEntries();
		int hist_1 = BGK_hist->GetEntries();
		int hist_2 = SIG_hist->GetEntries();

		
		

		Data_hist->Draw("hist");
		BGK_hist->SetMarkerStyle(20);
		BGK_hist->Draw("hist same");
		SIG_hist->Draw("hist same");
		 if(i==0){
			auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
			legend->AddEntry(Data_hist, Form("%s (%d)", label1_sig.Data(),hist), "f");
			legend->AddEntry(BGK_hist, Form("%s (%d)", label2_sig.Data(),hist_1), "f");
			legend->AddEntry(SIG_hist, Form("%s (%d)", label3_sig.Data(),hist_2), "f");

			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");
			// Labels
			double x_top_label = 0.90;

			TPaveText *CLAS12_Internal = new TPaveText(0.2, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
			CLAS12_Internal->SetFillStyle(4050);
			CLAS12_Internal->SetLineColor(0);
			CLAS12_Internal->SetTextFont(42);
			CLAS12_Internal->SetTextSize(0.0599401);
			CLAS12_Internal->SetBorderSize(0);
			CLAS12_Internal->AddText(Form("Tagged AI Variables"));
			CLAS12_Internal->Draw();
		 }

		

	}
	

	cancG0->SaveAs(name_pdf+".pdf(");

	cancG0->Clear();
	cancG0->Divide(3,4,0.01,0.01,0);
	for (int i = 0; i < labels1D.size(); i++)
	{
		cancG0->cd(i+1);
		//////////////////////////////////////////////////
		// Set options for each label
		TString label1 = labels1D[i][0];
		TString xAxis_label = labels1D[i][1];
		TString range_x_option = labels1D[i][2];
		TString min_x_option = labels1D[i][3];
		TString max_x_option = labels1D[i][4];
		TString binning_x_option = labels1D[i][5];
		TString nb_x_bins = labels1D[i][6];
		TString string_cut = labels1D[i][7];
		TString legend_option = labels1D[i][8];
		TString output_string = labels1D[i][9];
		//////////////////////////////////////////////////

		cout << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << "Doing 1D " << label1 << " plot" << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << endl;

		TCut cut = string_cut.Data();

		TString cut_string = cut.GetTitle();

		float min_X_histo_ini;
		float max_X_histo_ini;

		int nBins_X = 20;

		if (range_x_option == "custom_range")
		{
			//cout << "here " << endl;
			min_X_histo_ini = stof((string)min_x_option.Data());
			max_X_histo_ini = stof((string)max_x_option.Data());
		}

		if (binning_x_option == "custom_binning")
		{
			//cout << "here " << nb_x_bins.Data() << endl;
			nBins_X = stoi((string)nb_x_bins.Data());
		}

		TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		Data1_tree->Draw(label1 + ">>Data_hist",cut_0 * cut);
		// Data1_tree->Draw("Q2>>Data_hist");
		// Styling for Signal (Blue Solid Fill)
		Data_hist->SetFillColorAlpha( kCyan-10  ,0.50);
		Data_hist->SetFillStyle(1001); // Solid fill
		Data_hist->SetLineColor(kCyan-8);
		Data_hist->SetLineWidth(2);
		Data_hist->SetTitle(";" + xAxis_label+";");
		
		
		TH1F *BGK_hist = new TH1F("BGK_hist", "BGK_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		//TCut weight = Form("%i/%i", entries_sig, entries_bgk);

		Data4_tree->Draw(label1 + ">>BGK_hist",cut_0 * cut);


		

		// Styling for Background (Red Hatched Fill)
		BGK_hist->SetFillColor(kRed-7);
		BGK_hist->SetFillStyle(3004); // Hatched pattern
		BGK_hist->SetLineColor(kRed-7);
		BGK_hist->SetLineWidth(2);
        BGK_hist->SetTitle(";" + xAxis_label+";");



		TH1F *SIG_hist = new TH1F("SIG_hist", "SIG_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		//TCut weight = Form("%i/%i", entries_sig, entries_bgk);

		Data3_tree->Draw(label1 + ">>SIG_hist",cut_0 * cut);
		SIG_hist->SetFillColor(kGreen+3);
		SIG_hist->SetFillStyle(3003); // DOt pattern
		SIG_hist->SetLineColor(kGreen+3);
		SIG_hist->SetLineWidth(2);
        SIG_hist->SetTitle(";" + xAxis_label+";");

		
		BGK_hist->Scale(1./BGK_hist->Integral(),"nosw2");
		Data_hist->Scale(1./Data_hist->Integral(),"nosw2");
		SIG_hist->Scale(1./SIG_hist->Integral(),"nosw2");
		//gStyle->SetLegendTextSize(0.06);

		Data_hist->SetMaximum(std::max({Data_hist->GetMaximum(), SIG_hist->GetMaximum(),BGK_hist->GetMaximum()}) * 1.2);

		int hist = Data_hist->GetEntries();
		int hist_1 = BGK_hist->GetEntries();
		int hist_2 = SIG_hist->GetEntries();

		
		

		Data_hist->Draw("hist");
		BGK_hist->SetMarkerStyle(20);
		BGK_hist->Draw("hist same");
		SIG_hist->Draw("hist same");
		 if(i==0){
			auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
			legend->AddEntry(Data_hist, Form("%s (%d)", label1_sig.Data(),hist), "f");
			legend->AddEntry(BGK_hist, Form("%s (%d)", label4_sig.Data(),hist_1), "f");
			legend->AddEntry(SIG_hist, Form("%s (%d)", label3_sig.Data(),hist_2), "f");

			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");
			// Labels
			double x_top_label = 0.90;

			TPaveText *CLAS12_Internal = new TPaveText(0.2, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
			CLAS12_Internal->SetFillStyle(4050);
			CLAS12_Internal->SetLineColor(0);
			CLAS12_Internal->SetTextFont(42);
			CLAS12_Internal->SetTextSize(0.0599401);
			CLAS12_Internal->SetBorderSize(0);
			CLAS12_Internal->AddText(Form("Tagged AI Variables"));
			CLAS12_Internal->Draw();
		 }

		

	}
	cancG0->SaveAs(name_pdf+".pdf");
	

	*/

	TCanvas *cancG1 = new TCanvas("", "can1", 1200, 1000);
	cancG1->Divide(3,4,0.01,0.01,0);
	for (int i = 0; i < labels1D.size(); i++)
	{
		cancG1->cd(i+1);
		//////////////////////////////////////////////////
		// Set options for each label
		TString label1 = labels1D[i][0];
		TString xAxis_label = labels1D[i][1];
		TString range_x_option = labels1D[i][2];
		TString min_x_option = labels1D[i][3];
		TString max_x_option = labels1D[i][4];
		TString binning_x_option = labels1D[i][5];
		TString nb_x_bins = labels1D[i][6];
		TString string_cut = labels1D[i][7];
		TString legend_option = labels1D[i][8];
		TString output_string = labels1D[i][9];
		//////////////////////////////////////////////////

		cout << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << "Doing 1D " << label1 << " plot" << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << endl;

		TCut cut = string_cut.Data();

		TString cut_string = cut.GetTitle();

		float min_X_histo_ini;
		float max_X_histo_ini;

		int nBins_X = 20;

		if (range_x_option == "custom_range")
		{
			//cout << "here " << endl;
			min_X_histo_ini = stof((string)min_x_option.Data());
			max_X_histo_ini = stof((string)max_x_option.Data());
		}

		if (binning_x_option == "custom_binning")
		{
			//cout << "here " << nb_x_bins.Data() << endl;
			nBins_X = stoi((string)nb_x_bins.Data());
		}

		TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		Data1_tree->Draw(label1 + ">>Data_hist",cut_0 * cut);
		// Data1_tree->Draw("Q2>>Data_hist");


		// Styling for Signal (Blue Solid Fill)
		Data_hist->SetFillColorAlpha( kCyan+2  ,0.50);
		Data_hist->SetFillStyle(1001); // Solid fill
		Data_hist->SetLineColor(kCyan-3);
		Data_hist->SetLineWidth(2);
		 Data_hist->SetTitle(";" + xAxis_label+";");

		
		
		TH1F *Data_aftercut_hist = new TH1F("Data_aftercut_hist", "Data_aftercut_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		//TCut weight = Form("%i/%i", entries_sig, entries_bgk);

		Data1_tree->Draw(label1 + ">>Data_aftercut_hist", Form(cut_0 * cut&&"scoreBDT>=%g",score_cut));

		// Styling for Background (Red Hatched Fill)
		Data_aftercut_hist->SetFillColor(kGreen+3);
		Data_aftercut_hist->SetFillStyle(3004); // Hatched pattern
		Data_aftercut_hist->SetLineColor(kGreen+3);
		Data_aftercut_hist->SetLineWidth(2);
        Data_aftercut_hist->SetTitle(";" + xAxis_label+";");
		

		
		Data_aftercut_hist->Scale(1./Data_aftercut_hist->Integral(),"nosw2");
		Data_hist->Scale(1./Data_hist->Integral(),"nosw2");
		//gStyle->SetLegendTextSize(0.06);

		Data_hist->SetMaximum(std::max(Data_hist->GetMaximum(), Data_aftercut_hist->GetMaximum()) * 1.2);

		int hist = Data_hist->GetEntries();
		int hist_1 = Data_aftercut_hist->GetEntries();

		
		

		Data_hist->Draw("hist");
		Data_aftercut_hist->SetMarkerStyle(20);
		Data_aftercut_hist->Draw("hist same");
		 if(i==0){
			auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
			legend->AddEntry(Data_hist, Form("%s (%d)", label1_sig.Data(),hist), "f");
			legend->AddEntry(Data_aftercut_hist, Form("%s After %g Cut(%d)", label1_sig.Data(),score_cut,hist_1), "f");

			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");
			// Labels
			double x_top_label = 0.90;

			TPaveText *CLAS12_Internal = new TPaveText(0.2, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
			CLAS12_Internal->SetFillStyle(4050);
			CLAS12_Internal->SetLineColor(0);
			CLAS12_Internal->SetTextFont(42);
			CLAS12_Internal->SetTextSize(0.0599401);
			CLAS12_Internal->SetBorderSize(0);
			CLAS12_Internal->AddText(Form("Tagged AI Variables"));
			CLAS12_Internal->Draw();
		 }

		

	}
	cancG1->SaveAs(name_pdf + ".pdf(");
	cancG1->Clear();
	cancG1->Divide(3,4,0.01,0.01,0);
	for (int i = 0; i < labels1D.size(); i++)
	{
		cancG1->cd(i+1);
		//////////////////////////////////////////////////
		// Set options for each label
		TString label1 = labels1D[i][0];
		TString xAxis_label = labels1D[i][1];
		TString range_x_option = labels1D[i][2];
		TString min_x_option = labels1D[i][3];
		TString max_x_option = labels1D[i][4];
		TString binning_x_option = labels1D[i][5];
		TString nb_x_bins = labels1D[i][6];
		TString string_cut = labels1D[i][7];
		TString legend_option = labels1D[i][8];
		TString output_string = labels1D[i][9];
		//////////////////////////////////////////////////

		cout << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << "Doing 1D " << label1 << " plot" << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << endl;

		TCut cut = string_cut.Data();

		TString cut_string = cut.GetTitle();

		float min_X_histo_ini;
		float max_X_histo_ini;

		int nBins_X = 20;

		if (range_x_option == "custom_range")
		{
			//cout << "here " << endl;
			min_X_histo_ini = stof((string)min_x_option.Data());
			max_X_histo_ini = stof((string)max_x_option.Data());
		}

		if (binning_x_option == "custom_binning")
		{
			//cout << "here " << nb_x_bins.Data() << endl;
			nBins_X = stoi((string)nb_x_bins.Data());
		}

		TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		Data2_tree->Draw(label1 + ">>Data_hist",cut_0 * cut);
		// Data1_tree->Draw("Q2>>Data_hist");


		// Styling for Signal (Blue Solid Fill)
		Data_hist->SetFillColorAlpha( kCyan+2  ,0.50);
		Data_hist->SetFillStyle(1001); // Solid fill
		Data_hist->SetLineColor(kCyan-3);
		Data_hist->SetLineWidth(2);
		 Data_hist->SetTitle(";" + xAxis_label+";");

		
		
		TH1F *Data_aftercut_hist = new TH1F("Data_aftercut_hist", "Data_aftercut_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
		//TCut weight = Form("%i/%i", entries_sig, entries_bgk);

		Data2_tree->Draw(label1 + ">>Data_aftercut_hist", Form(cut_0 * cut&&"scoreBDT>=%g",score_cut));

		// Styling for Background (Red Hatched Fill)
		Data_aftercut_hist->SetFillColor(kGreen+3);
		Data_aftercut_hist->SetFillStyle(3004); // Hatched pattern
		Data_aftercut_hist->SetLineColor(kGreen+3);
		Data_aftercut_hist->SetLineWidth(2);
        Data_aftercut_hist->SetTitle(";" + xAxis_label+";");
		

		
		Data_aftercut_hist->Scale(1./Data_aftercut_hist->Integral(),"nosw2");
		Data_hist->Scale(1./Data_hist->Integral(),"nosw2");
		//gStyle->SetLegendTextSize(0.06);

		Data_hist->SetMaximum(std::max(Data_hist->GetMaximum(), Data_aftercut_hist->GetMaximum()) * 1.2);

		int hist = Data_hist->GetEntries();
		int hist_1 = Data_aftercut_hist->GetEntries();

		
		

		Data_hist->Draw("hist");
		Data_aftercut_hist->SetMarkerStyle(20);
		Data_aftercut_hist->Draw("hist same");
		 if(i==0){
			auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
			legend->AddEntry(Data_hist, Form("%s (%d)", label2_sig.Data(),hist), "f");
			legend->AddEntry(Data_aftercut_hist, Form("%s After %g Cut(%d)", label2_sig.Data(),score_cut,hist_1), "f");
			legend->SetTextSize(0.02);
			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");
			// Labels
			double x_top_label = 0.90;

			TPaveText *CLAS12_Internal = new TPaveText(0.2, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
			CLAS12_Internal->SetFillStyle(4050);
			CLAS12_Internal->SetLineColor(0);
			CLAS12_Internal->SetTextFont(42);
			CLAS12_Internal->SetTextSize(0.0599401);
			CLAS12_Internal->SetBorderSize(0);
			CLAS12_Internal->AddText(Form("Tagged AI Variables"));
			CLAS12_Internal->Draw();
		 }

		

	}
	cancG1->SaveAs(name_pdf + ".pdf");


	TCanvas *cancG2 = new TCanvas("", "can2", 1200, 1000);
	//////////////////////////////////////////////////
	// Set options for each label
	TString label1 = Invariant[0][0];
	TString xAxis_label = Invariant[0][1];
	TString range_x_option = Invariant[0][2];
	TString min_x_option = Invariant[0][3];
	TString max_x_option = Invariant[0][4];
	TString binning_x_option = Invariant[0][5];
	TString nb_x_bins = Invariant[0][6];
	TString string_cut = Invariant[0][7];
	TString legend_option = Invariant[0][8];
	TString output_string = Invariant[0][9];
	//////////////////////////////////////////////////

	cout << endl;
	cout << "//////////////////////////////////////////////////" << endl;
	cout << "Doing 1D " << label1 << " plot" << endl;
	cout << "//////////////////////////////////////////////////" << endl;
	cout << endl;

	TCut cut = string_cut.Data();

	TString cut_string = cut.GetTitle();

	float min_X_histo_ini;
	float max_X_histo_ini;

	int nBins_X = 20;

	if (range_x_option == "custom_range"){
		min_X_histo_ini = stof((string)min_x_option.Data());
		max_X_histo_ini = stof((string)max_x_option.Data());
	}

	if (binning_x_option == "custom_binning"){
		nBins_X = stoi((string)nb_x_bins.Data());
	}

	TH1F *Data_aftercut_hist = new TH1F("Data_aftercut_hist", "Data_aftercut_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
	Data1_tree->Draw(label1 + ">>Data_aftercut_hist",  Form(cut_0 * cut&&" scoreBDT >= %g && abs(Mxepe) < %g", score_cut, 0.1));

	// Styling for Background (Red Hatched Fill)
	Data_aftercut_hist->SetFillColor(kGreen+3);
	Data_aftercut_hist->SetFillStyle(3004); // Hatched pattern
	Data_aftercut_hist->SetLineColor(kGreen+3);
	Data_aftercut_hist->SetLineWidth(2);
    Data_aftercut_hist->SetTitle(";" + xAxis_label+";");
		

	//int hist = Data_hist->GetEntries();
	int hist_1 = Data_aftercut_hist->GetEntries();


	//Data_hist->Draw("hist");
	//Data_aftercut_hist->SetMaximum(10);
	Data_aftercut_hist->Draw("hist");
	
	auto legend = new TLegend();
	//legend->AddEntry(Data_hist, Form("%s (%d)", label1_sig.Data(),hist), "f");
	legend->AddEntry(Data_aftercut_hist, Form("%s After %g Cut(%d)", label1_sig.Data(),score_cut,hist_1), "l");

	legend->SetFillStyle(0);
	legend->SetLineWidth(0);
	legend->Draw("same ");

	cancG2->SaveAs(name_pdf + ".pdf(");

	TCanvas *cancG3 = new TCanvas("", "can3", 1200, 1000);
	

	TH1F *Data2_aftercut_hist = new TH1F("Data2_aftercut_hist", "Data2_aftercut_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
	Data2_tree->Draw(label1 + ">>Data2_aftercut_hist",  Form(cut_0 * cut&&"scoreBDT >= %g && abs(Mxepe) < %g", score_cut, 0.1));

	// Styling for Background (Red Hatched Fill)
	Data2_aftercut_hist->SetFillColor(kGreen+3);
	Data2_aftercut_hist->SetFillStyle(3004); // Hatched pattern
	Data2_aftercut_hist->SetLineColor(kGreen+3);
	Data2_aftercut_hist->SetLineWidth(2);
    Data2_aftercut_hist->SetTitle(";" + xAxis_label+";");
		

	 //hist = Data2_hist->GetEntries();
	 hist_1 = Data2_aftercut_hist->GetEntries();


	//Data2_hist->Draw("hist");
	//Data2_aftercut_hist->SetMarkerStyle(20);
	Data2_aftercut_hist->Draw("hist");
	
	legend->Clear();
	//legend->AddEntry(Data2_hist, Form("%s (%d)", label2_sig.Data(),hist), "f");
	legend->AddEntry(Data2_aftercut_hist, Form("%s After %g Cut(%d)", label2_sig.Data(),score_cut,hist_1), "l");

	legend->SetFillStyle(0);
	legend->SetLineWidth(0);
	legend->Draw("same ");

	


	cancG3->SaveAs(name_pdf + ".pdf)");


	

	gApplication->Terminate();

	return 0;
}
