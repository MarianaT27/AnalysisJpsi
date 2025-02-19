#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

void compareRootFiles(TString file1, TString file2,TString pdfname) {
    gStyle->SetOptStat(0);
    // Open the ROOT files
    TFile *f1 = TFile::Open("/lustre24/expphy/volatile/clas12/mtenorio/Root/S18/"+file1+".root");
    TFile *f2 = TFile::Open("/lustre24/expphy/volatile/clas12/mtenorio/Root/S18/"+file2+".root");

    if (!f1 || !f2 || f1->IsZombie() || f2->IsZombie()) {
        std::cerr << "Error: One or both ROOT files could not be opened.\n";
        return;
    }

    // Get the TTree named "results" from both files
    TTree *tree1 = (TTree*)f1->Get("results");
    TTree *tree2 = (TTree*)f2->Get("results");

    if (!tree1 || !tree2) {
        std::cerr << "Error: TTree 'results' not found in one or both files.\n";
        return;
    }

    // Variables to read from trees
    double beamCharge_mC, electron_p, positron_p, invariantMass, Q2, MM, time_ee, Egamma, W;

    // Set branch addresses
    for (TTree* tree : {tree1, tree2}) {
        tree->SetBranchAddress("beamCharge_mC", &beamCharge_mC);
        tree->SetBranchAddress("electron_p", &electron_p);
        tree->SetBranchAddress("positron_p", &positron_p);
        tree->SetBranchAddress("invariantMass", &invariantMass);
        tree->SetBranchAddress("Q2", &Q2);
        tree->SetBranchAddress("MM", &MM);
        tree->SetBranchAddress("time_ee", &time_ee);
        tree->SetBranchAddress("Egamma", &Egamma);
        tree->SetBranchAddress("W", &W);
    }

    
     // Define binning for each variable: {nBins, xMin, xMax}
     std::map<std::string, std::tuple<int, double, double>> binning = {
        {"electron_p", {100, 0, 10}},
        {"positron_p", {100, 0, 10}},
        {"time_ee", {100, -5, 5}},
        {"Egamma", {100, 0, 0}},
        {"MM", {100, 0, 3.5}},
        {"invariantMass", {80, 0.5, 3.5}},
        {"Q2", {100, 0, 0}},
        {"W", {120, 4, 4.6}}
    };
    // Define histograms
    std::map<std::string, TH1D*> h1, h2;

    for (const auto& [var, binParams] : binning) {
        int nBins;
        double xMin, xMax;
        std::tie(nBins, xMin, xMax) = binParams;

        h1[var] = new TH1D(("h1_" + var).c_str(), (var + " Comparison; " + var + "; Events").c_str(), nBins, xMin, xMax);
        h2[var] = new TH1D(("h2_" + var).c_str(), (var + " Comparison; " + var + "; Events").c_str(), nBins, xMin, xMax);
    }

    double totalBeamCharge1 = 0, totalBeamCharge2 = 0;

    // Fill histograms for tree1
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; i++) {
        tree1->GetEntry(i);
        if (beamCharge_mC > 0 ) {  // Avoid division by zero
            double weight = 1.0;
            h1["electron_p"]->Fill(electron_p, weight);
            h1["positron_p"]->Fill(positron_p, weight);
            h1["time_ee"]->Fill(time_ee, weight);
            h1["Egamma"]->Fill(Egamma, weight);
            if(fabs(time_ee)<2.5&& Egamma>8.1){
                h1["Q2"]->Fill(Q2, weight);
                h1["MM"]->Fill(MM, weight);
                h1["W"]->Fill(W, weight);
                //if(MM>0.7&&MM<1.3)
                    h1["invariantMass"]->Fill(invariantMass, weight);
            }
            totalBeamCharge1 += beamCharge_mC;
        }
    }

    // Fill histograms for tree2
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; i++) {
        tree2->GetEntry(i);
        if (beamCharge_mC > 0) {  // Avoid division by zero
            double weight = 1.0;
            h2["electron_p"]->Fill(electron_p, weight);
            h2["positron_p"]->Fill(positron_p, weight);
            h2["time_ee"]->Fill(time_ee, weight);
            h2["Egamma"]->Fill(Egamma, weight);

            if (fabs(time_ee) < 2.5 && Egamma > 8.1) {
                h2["Q2"]->Fill(Q2, weight);
                h2["MM"]->Fill(MM, weight);
                h2["W"]->Fill(W, weight);

                //if (MM > 0.7 && MM < 1.3)
                    h2["invariantMass"]->Fill(invariantMass, weight);
            }
            totalBeamCharge2 += beamCharge_mC;
        }
    }

    
    // Normalize histograms by total beam charge
    for (const auto& [var, _] : binning) {
        if (totalBeamCharge1 > 0) h1[var]->Scale(1.0 / totalBeamCharge1);
        if (totalBeamCharge2 > 0) h2[var]->Scale(1.0 / totalBeamCharge2);

        h1[var]->SetMaximum(std::max(h1[var]->GetMaximum(), h2[var]->GetMaximum()) * 1.2);

		
    }

    

    // Create a canvas to plot the histograms
    TCanvas *c = new TCanvas("c", "Comparison Plots", 800, 1200);
    c->Divide(2, 4); // Divide canvas into 3x3 subplots

    // Define the desired order explicitly
    std::vector<std::string> plot_order = {
        "electron_p", "positron_p", "time_ee", "Egamma",
        "MM", "invariantMass", "Q2", "W"
    };

    int padNumber = 1;
    for (const auto& var : plot_order) {
        c->cd(padNumber);
        cout<<"Plotting "<<var<<endl;
        int hist_1 = h1[var]->GetEntries();
		int hist_2= h2[var]->GetEntries();

        // Style for h1 (red)
        h1[var]->SetLineColor(kRed + 1);
        h1[var]->SetLineWidth(1);
        h1[var]->SetFillColorAlpha(kRed, 0.3);  // Light red transparency
        h1[var]->Draw("HIST");

        // Style for h2 (blue)
        h2[var]->SetLineColor(kBlue + 1);
        h2[var]->SetLineWidth(1);
        //h2[var]->SetFillColorAlpha(kBlue, 0.3);  // Light blue transparency
        h2[var]->Draw("HIST SAME");

        TLegend *legend = new TLegend();
        legend->SetTextSize(0.03);
        legend->SetFillStyle(0);
		legend->SetLineWidth(0);
        legend->AddEntry(h1[var],Form("%s (%d)", file1.Data(),hist_1), "l");
        legend->AddEntry(h2[var], Form("%s (%d)", file2.Data(),hist_2), "l");
        legend->Draw();

        padNumber++;
    }
    c->SaveAs(pdfname+".pdf");

    std::cout << "Comparison plots saved to '"<<pdfname<<".pdf'\n";

}