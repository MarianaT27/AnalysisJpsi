#ifndef Plot
#define Plot

void plot_ratio(TString name, double *Var)
{
  double Bins[10] = {8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  double BinsE[10] = {0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  double error[10] = {0};

  TGraphErrors *gr = new TGraphErrors(10, Bins, Var, BinsE, error);
  gr->SetMarkerStyle(8);
  gr->SetTitle(name + "; E_{#gamma} Bin; Ratio");
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  gr->GetYaxis()->SetRangeUser(0, 1);
  gr->Draw("AEP");
}

void plot_ratio_same(double *Var1, double *Var2, double *Var3, double *Var4, TString name_Y = "Ratio NJpsi/Fit")
{
  double Bins[10] = {8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  double BinsE[10] = {0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};
  // double BinsE[10]={0};

  double error[10] = {0};
  TString name[4] = {"Gauss", "CrystallBall", "GaussExp", "dNdMee"};

  TGraphErrors *gr = new TGraphErrors(10, Bins, Var1, BinsE, error);
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  gr->GetYaxis()->SetRangeUser(0.0, 2.0);
  gr->SetMarkerStyle(8);
  gr->SetMarkerColor(kRed);
  gr->SetTitle("; E_{#gamma} Bin; " + name_Y);
  gr->Draw("AEP");

  TGraphErrors *gr2 = new TGraphErrors(10, Bins, Var2, BinsE, error);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerColor(kAzure - 3);
  gr2->Draw("P");

  TGraphErrors *gr3 = new TGraphErrors(10, Bins, Var3, BinsE, error);
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerColor(kMagenta + 1);
  gr3->Draw("P");

  TGraphErrors *gr4 = new TGraphErrors(10, Bins, Var4, BinsE, error);
  gr4->SetMarkerStyle(8);
  gr4->SetMarkerColor(kOrange - 2);
  gr4->Draw("P");

  auto legend = new TLegend(0.54, 0.87, 0.9, 0.75); // 0.54, 0.87, 0.90, 0.60
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(gr, name[0], "p");
  legend->AddEntry(gr2, name[1], "p");
  legend->AddEntry(gr3, name[2], "p");
  legend->AddEntry(gr4, name[3], "p");
  legend->Draw("same");
}

void plot_function(int n, double *X, double *Y, int color)
{
  TGraph *gr = new TGraph(n, X, Y);
  // gr->GetXaxis()->SetRangeUser(2.5,3.5);
  // gr->GetYaxis()->SetRangeUser(0.5,2.5);
  gr->SetMarkerStyle(8);
  gr->SetLineColor(color);
  gr->SetLineWidth(3);
  // gr->SetTitle("; E_{#gamma} Bin; Ratio NJpsi/Fit");
  gr->Draw("AL");
}

void plot_function_same(int n, double *X, double *Y, int color)
{
  TGraph *gr = new TGraph(n, X, Y);
  // gr->GetXaxis()->SetRangeUser(2.5,3.5);
  // gr->GetYaxis()->SetRangeUser(0.5,2.5);
  gr->SetMarkerStyle(8);
  gr->SetLineColor(color);
  gr->SetLineWidth(3);
  // gr->SetTitle("; E_{#gamma} Bin; Ratio NJpsi/Fit");
  gr->Draw("L");
}

void plot_results(string pdfname)
{
  TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
  string pdf_original = pdfname + ".pdf";
  gStyle->SetOptFit(1011);

  TPaveText *pt = new TPaveText();
  auto legend = new TLegend(); // 0.54, 0.87, 0.90, 0.60

  // Histograms to print
  Fit_Function Fit_func;
  Fit_func.Set_Limits(2.6, 3.3);
  double nb_JPsi;

  // h_Invariant->GetYaxis()->SetRangeUser(0,270);
  h_Invariant->Draw();
  Fit_func.Fit_dNdMee(h_Invariant);
  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);

  pt->AddText("Fit_dNdMee");
  pt->SetTextSize(0.03);
  pt->Draw();
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.1f) ", nb_JPsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  h_mm_vs_invariantmass->Draw("colz");
  can->Print((pdf_original + "(").c_str());

  h_qq_vs_invariantmass->Draw("colz");
  can->Print((pdf_original + "(").c_str());

  h_MM->Draw();
  can->Print((pdf_original + "(").c_str());
  // Fit_func.Set_Limits(-0.04, 0.04);
  // Fit_func.Fit_Gauss(h_MM);

  h_Q2->Draw();
  can->Print((pdf_original + ")").c_str());
}

void plot_results_t(string pdfname)
{

  TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
  string pdf_original = pdfname + ".pdf";
  gStyle->SetOptFit(1011);
  gStyle->SetStatFontSize(0.03);

  TPaveText *pt = new TPaveText();
  auto legend = new TLegend(); // 0.54, 0.87, 0.90, 0.60

  // Histograms to print
  Fit_Function Fit_func;
  Fit_func.Set_Limits(2.0, 3.4);
  double nb_JPsi;

  // h_Invariant->GetYaxis()->SetRangeUser(0,270);
  h_Invariant->Draw();
  Fit_func.Fit_Gauss_Pol2(h_Invariant);
  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);

  pt->AddText("Fit_Gauss_Pol2");
  pt->SetTextSize(0.03);
  pt->Draw();
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.1f) ", nb_JPsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  h_Invariant->Draw();
  Fit_func.Fit_Gauss_Exp(h_Invariant);
  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);

  pt->AddText("Fit_Gauss_Exp");
  pt->SetTextSize(0.03);
  pt->Draw();
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.1f) ", nb_JPsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  h_Invariant->Draw();
  Fit_func.Fit_GaussExp(h_Invariant);
  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);

  pt->AddText("Fit_GaussExp");
  pt->SetTextSize(0.03);
  pt->Draw();
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.1f) ", nb_JPsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  h_Invariant->Draw();
  Fit_func.Fit_dNdMee(h_Invariant);
  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);

  pt->AddText("Fit_dNdMee");
  pt->SetTextSize(0.03);
  pt->Draw();
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.1f) ", nb_JPsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  h_mm_vs_invariantmass_t->Draw("colz");
  can->Print((pdf_original + "(").c_str());

  h_mm_vs_invariantmass_no_t->Draw("colz");
  can->Print((pdf_original + "(").c_str());

  h_MX->Draw();
  can->Print((pdf_original + "(").c_str());

  h_W->Draw();
  can->Print((pdf_original + "(").c_str());

  h_Q2->Draw();
  can->Print((pdf_original + ")").c_str());
}

void plotgraph(double_t *variable, double_t *error, TString name, double_t OValue = 0, double_t OValueE = 0, int color = 600)
{
  // double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  // double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  double Bins[4] = {8.625, 9.225, 9.73, 10.3};
  double BinsE[4] = {0.425, 0.205, 0.27, 0.3};

  TGraphErrors *gr = new TGraphErrors(4, Bins, variable, BinsE, error);
  gr->SetMarkerStyle(8);
  gr->SetTitle(name + "; E_{#gamma} ; # Events");
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  gr->Draw("AEP");

  if (OValue != 0 && OValueE != 0)
  {
    TLine *line = new TLine(8.2, OValue, 10.6, OValue);
    TBox *box = new TBox(8.2, OValue - OValueE, 10.6, OValue + OValueE);
    line->SetLineColor(color);
    box->SetFillColorAlpha(color, 0.15);
    line->Draw("same");
    box->Draw("same");
  }
}

void plot(double_t *variable, TString name = "", int key = 0)
{
  // double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  // double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  // double Bins[4] ={0.375, 1.125, 2.25};
  // double BinsE[4]={0.375, 0.375, 0.75};

  double Bins[4] = {0.5, 1.5, 2.75};
  double BinsE[4] = {0.5, 0.5, 0.75};

  double error[4] = {0};

  double Var_3[3] = {0};
  TString out;

  if (key == 0)
  {
    Var_3[0] = variable[0];
    Var_3[1] = variable[2];
    Var_3[2] = variable[4];
    out = name + " E_{#gamma}=[8.2,9.8]";
  }
  else
  {
    Var_3[0] = variable[1];
    Var_3[1] = variable[3];
    Var_3[2] = variable[5];
    out = name + " E_{#gamma}=[9.8,10.6]";
  }

  TGraphErrors *gr = new TGraphErrors(3, Bins, Var_3, BinsE, error);
  gr->SetMarkerStyle(8);
  gr->SetTitle(out + "; -t ; BSA");
  gr->GetXaxis()->SetRangeUser(0, 3.5);
  gr->GetYaxis()->SetRangeUser(-1, 1);
  gr->Draw("AEP");
}

void plot2graph(double_t variable[4][10], double_t error[4][10], TString name, int g1 = 0, int g2 = 1, double_t OValue = 0, double_t OValueE = 0, double_t OValue2 = 0, double_t OValueE2 = 0)
{
  // double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  // double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  double Bins[4] = {8.625, 9.225, 9.73, 10.3};
  double BinsE[4] = {0.425, 0.205, 0.27, 0.3};

  TGraphErrors *gr = new TGraphErrors(4, Bins, variable[g1], BinsE, error[g1]);
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(kBlue);
  gr->SetTitle(name + "; E_{#gamma} ; " + name);
  gr->GetXaxis()->SetRangeUser(8.2, 10.6);
  if (name == "Mean")
  {
    gr->GetYaxis()->SetRangeUser(3.01, 3.12);
  }
  else if (name == "Sigma")
  {
    gr->GetYaxis()->SetRangeUser(0.0, 0.10);
  }

  gr->Draw("AEP");

  TGraphErrors *gr2 = new TGraphErrors(4, Bins, variable[g2], BinsE, error[g2]);
  gr2->SetMarkerStyle(23);
  gr2->SetMarkerColor(kRed);
  gr2->Draw("P");

  auto legend = new TLegend(0.25, 0.10);
  legend->AddEntry(gr, lnames[g1].c_str(), "p");
  legend->AddEntry(gr2, lnames[g2].c_str(), "p");
  legend->Draw("same");

  if (OValue != 0 && OValueE != 0)
  {
    TLine *line = new TLine(8.2, OValue, 10.6, OValue);
    TBox *box = new TBox(8.2, OValue - OValueE, 10.6, OValue + OValueE);
    line->SetLineColor(kBlue);
    box->SetFillColorAlpha(kBlue, 0.15);
    line->Draw("same");
    box->Draw("same");
  }

  if (OValue2 != 0 && OValueE2 != 0)
  {
    TLine *line2 = new TLine(8.2, OValue2, 10.6, OValue2);
    TBox *box2 = new TBox(8.2, OValue2 - OValueE2, 10.6, OValue2 + OValueE2);
    line2->SetLineColor(kRed);
    box2->SetFillColorAlpha(kRed, 0.15);
    line2->Draw("same");
    box2->Draw("same");
  }
}

#endif