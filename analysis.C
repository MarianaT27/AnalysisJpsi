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
#include "clas12reader.h"
#include "bib/Utils.h"
#include "bib/Functions.h"
#include "bib/h_Fit.h"
#include "bib/Plot.h"
#include "bib/Basic.h"
#include "bib/Update.h"

       

void analysis() { 


    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    cout<<"Running analysis..."<<endl;
    //Read 2 different root files 
    TH1F *h1_Mee1_case1_in;
    TH1F *h1_Mee2_case1_in;
    TH1F *h1_Mee1_case2_in;
    TH1F *h1_Mee2_case2_in;
    TH1F *h1_Mee1_case3_in;
    TH1F *h1_Mee2_case3_in;

    TH1F *h1_Mee1_case1_out;
    TH1F *h1_Mee2_case1_out;
    TH1F *h1_Mee1_case2_out;
    TH1F *h1_Mee2_case2_out;
    TH1F *h1_Mee1_case3_out;
    TH1F *h1_Mee2_case3_out;

    TFile *file_in = new TFile("/lustre19/expphy/volatile/clas12/mtenorio/Root/F18in_Bg_est.root", "READ");
    h1_Mee1_case1_in = (TH1F*)file_in->Get("h1_Mee1_case1");
    h1_Mee2_case1_in = (TH1F*)file_in->Get("h1_Mee2_case1");
    h1_Mee1_case2_in = (TH1F*)file_in->Get("h1_Mee1_case2");
    h1_Mee2_case2_in = (TH1F*)file_in->Get("h1_Mee2_case2");
    h1_Mee1_case3_in = (TH1F*)file_in->Get("h1_Mee1_case3");
    h1_Mee2_case3_in = (TH1F*)file_in->Get("h1_Mee2_case3");

    TFile *file_out = new TFile("/lustre19/expphy/volatile/clas12/mtenorio/Root/F18out_Bg_est.root", "READ");
    h1_Mee1_case1_out = (TH1F*)file_out->Get("h1_Mee1_case1");
    h1_Mee2_case1_out = (TH1F*)file_out->Get("h1_Mee2_case1");
    h1_Mee1_case2_out = (TH1F*)file_out->Get("h1_Mee1_case2");
    h1_Mee2_case2_out = (TH1F*)file_out->Get("h1_Mee2_case2");
    h1_Mee1_case3_out = (TH1F*)file_out->Get("h1_Mee1_case3");
    h1_Mee2_case3_out = (TH1F*)file_out->Get("h1_Mee2_case3");

    TCanvas *can = new TCanvas("can", "canvas", 1500, 1000);
    TString outname = "./R_bg_est/Test_analysis.pdf";
    gStyle->SetOptStat(0);
    can->Divide(3,2);
    can->cd(1);
    gPad->SetLogy();
    h1_Mee2_case1_in->SetLineColor(kRed);
    h1_Mee1_case1_in->SetMaximum(std::max(h1_Mee1_case1_in->GetMaximum(),h1_Mee2_case1_in->GetMaximum()) * 1.2);
    h1_Mee1_case1_in->Draw();
    h1_Mee2_case1_in->Draw("same");
    TLegend *legend_gr = new TLegend();
    legend_gr->AddEntry(h1_Mee1_case1_in, "(e^{-}e^{+}p)", "l");
    legend_gr->AddEntry(h1_Mee2_case1_in, "(e^{-}e^{-}p)", "l");
    legend_gr->SetLineWidth(0);
    legend_gr->Draw();
    can->cd(2);
    gPad->SetLogy();
    h1_Mee2_case2_in->SetLineColor(kRed);
    h1_Mee1_case2_in->SetMaximum(std::max(h1_Mee1_case2_in->GetMaximum(),h1_Mee2_case2_in->GetMaximum()) * 1.2);
    h1_Mee1_case2_in->Draw();
    h1_Mee2_case2_in->Draw("same");
    can->cd(3);
    gPad->SetLogy();
    h1_Mee2_case3_in->SetLineColor(kRed);
    h1_Mee1_case3_in->SetMaximum(std::max(h1_Mee1_case3_in->GetMaximum(),h1_Mee2_case3_in->GetMaximum()) * 1.2);
    h1_Mee1_case3_in->Draw();
    h1_Mee2_case3_in->Draw("same");

    can->cd(4);
    gPad->SetLogy();
    h1_Mee2_case1_out->SetLineColor(kRed);
    h1_Mee1_case1_out->SetMaximum(std::max(h1_Mee1_case1_out->GetMaximum(),h1_Mee2_case1_out->GetMaximum()) * 1.2);
    h1_Mee1_case1_out->Draw();
    h1_Mee2_case1_out->Draw("same");
    can->cd(5);
    gPad->SetLogy();
    h1_Mee2_case2_out->SetLineColor(kRed);
    h1_Mee1_case2_out->SetMaximum(std::max(h1_Mee1_case2_out->GetMaximum(),h1_Mee2_case2_out->GetMaximum()) * 1.2);
    h1_Mee1_case2_out->Draw();
    h1_Mee2_case2_out->Draw("same");
    can->cd(6);
    gPad->SetLogy();
    h1_Mee2_case3_out->SetLineColor(kRed);
    h1_Mee1_case3_out->SetMaximum(std::max(h1_Mee1_case3_out->GetMaximum(),h1_Mee2_case3_out->GetMaximum()) * 1.2);
    h1_Mee1_case3_out->Draw();
    h1_Mee2_case3_out->Draw("same");

    can->Print(outname+"(");


    TH1F *h1_ratio_case1_in = new TH1F(*h1_Mee2_case1_in);
    TH1F *h1_ratio_case2_in = new TH1F(*h1_Mee2_case2_in);
    TH1F *h1_ratio_case3_in = new TH1F(*h1_Mee2_case3_in);
    TH1F *h1_ratio_case1_out = new TH1F(*h1_Mee2_case1_out);
    TH1F *h1_ratio_case2_out = new TH1F(*h1_Mee2_case2_out);
    TH1F *h1_ratio_case3_out = new TH1F(*h1_Mee2_case3_out);

    h1_ratio_case1_in->Divide(h1_Mee1_case1_in);
    h1_ratio_case2_in->Divide(h1_Mee1_case2_in);
    h1_ratio_case3_in->Divide(h1_Mee1_case3_in);
    h1_ratio_case1_out->Divide(h1_Mee1_case1_out);
    h1_ratio_case2_out->Divide(h1_Mee1_case2_out);
    h1_ratio_case3_out->Divide(h1_Mee1_case3_out);

    TH1F *h1_R_case1_in = new TH1F(*h1_ratio_case1_in);
    TH1F *h1_R_case2_in = new TH1F(*h1_ratio_case2_in);
    TH1F *h1_R_case3_in = new TH1F(*h1_ratio_case3_in);

    h1_R_case1_in->Divide(h1_ratio_case1_out);
    h1_R_case2_in->Divide(h1_ratio_case2_out);
    h1_R_case3_in->Divide(h1_ratio_case3_out);


    can->Clear();
    TCanvas *can2 = new TCanvas("can", "canvas", 1500, 500);
    can2->Divide(3,1);
    can2->cd(1);
    gPad->SetLogy();
    h1_ratio_case1_in->SetTitle("Inbending Ratio=(e^{-}e^{-}p)/(e^{-}e^{+}p)");
    h1_ratio_case1_in->SetLineColor(kRed);
    h1_ratio_case2_in->SetLineColor(kBlue);
    h1_ratio_case3_in->SetLineColor(kBlack);
    h1_ratio_case1_in->Draw();
    h1_ratio_case2_in->Draw("same");
    h1_ratio_case3_in->Draw("same");
    legend_gr->Clear();
    legend_gr->AddEntry(h1_ratio_case1_in, "pt/p<0.05", "l");
    legend_gr->AddEntry(h1_ratio_case2_in, "pt/p<0.1", "l");
    legend_gr->AddEntry(h1_ratio_case3_in, "Q2<0.2", "l");
    legend_gr->SetLineWidth(0);
    legend_gr->Draw();
    can2->cd(2);
    gPad->SetLogy();
    h1_ratio_case1_out->SetTitle("Outbending Ratio=(e^{-}e^{-}p)/(e^{-}e^{+}p)");
    h1_ratio_case1_out->SetLineColor(kRed);
    h1_ratio_case2_out->SetLineColor(kBlue);
    h1_ratio_case3_out->SetLineColor(kBlack);
    h1_ratio_case1_out->Draw();
    h1_ratio_case2_out->Draw("same");
    h1_ratio_case3_out->Draw("same");

    can2->cd(3);
    gPad->SetLogy();
    h1_R_case1_in->SetTitle("R_{In}/R_{Out}");
    h1_R_case1_in->SetLineColor(kRed);
    h1_R_case2_in->SetLineColor(kBlue);
    h1_R_case3_in->SetLineColor(kBlack);
    h1_R_case1_in->Draw();
    h1_R_case2_in->Draw("same");
    h1_R_case3_in->Draw("same");
    

    can2->Print(outname+")");


    //Can I limit the histograms to the region of 1.2 to 3.2?
    //Calculate ratios R_in/out Mee1/Mee2 for case 1, 2 and 3 plot in same canvas (2 plots from this with 3 histograms each)
    //Calculate ratio r2=R_in/R_out for each case (1 plot with 3 histograms)
    //Calculate r as r(Mee)=sqrt(r2(Mee)) for each mass point
    //Calculate Errors as error_r=error_r2/sqrt(2)



    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";
}
