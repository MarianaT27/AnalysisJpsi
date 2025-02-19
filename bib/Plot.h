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
  gStyle->SetStatFontSize(0.020);
  double_t par_Jpsi;
  double_t err_Jpsi;

  TPaveText *pt = new TPaveText();
  auto legend = new TLegend(); // 0.54, 0.87, 0.90, 0.60

  // Histograms to print
  Fit_Function Fit_func;
  Fit_func.Set_Limits(2.0, 3.5);
  double nb_JPsi;

  h_Invariant->Draw();
  Fit_func.Fit_dNdMee(h_Invariant);

  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);
  par_Jpsi = Fit_func.f->GetParameter(0);
  err_Jpsi = Fit_func.f->GetParError(0);
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.2f +/- %3.2f) ", par_Jpsi, err_Jpsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  gStyle->SetStatFontSize(0.03);

  h_mm_vs_invariantmass_t->Draw("colz");
  can->Print((pdf_original + "(").c_str());

  h_mm_vs_invariantmass_no_t->Draw("colz");
  can->Print((pdf_original + "(").c_str());

  h_W->Draw();
  can->Print((pdf_original + "(").c_str());

  h_Q2->Draw();
  can->Print((pdf_original + "(").c_str());

  h_MX->Draw();
  // Fit_Function Fit_func_Mx;
  // Fit_func_Mx.Set_Limits(0.5,1.5);
  // Fit_func_Mx.Fit_GaussMx_Pol2(h_MX);
  can->Print((pdf_original + ")").c_str());

  /*TCanvas *can2 = new TCanvas("can2", "canvas",2100,1400);

  can2->Divide(3,2);
  for(int i=0;i<6;i++){
    can2->cd(i+1);
    h_Data_Var[i]->Draw();
    Fit_func_Mx.Fit_GaussMx_Pol2(h_Data_Var[i]);
    if(pdfname.find("FTcorrON")!= std::string::npos){
      mean_FT[i]=Fit_func_Mx.f->GetParameter(1);
      mean_FT_Error[i]=Fit_func_Mx.f->GetParError(1);
      sigma_FT[i]=Fit_func_Mx.f->GetParameter(2);
      sigma_FT_Error[i]=Fit_func_Mx.f->GetParError(2);
     // cout<<"HEre"<<endl;
    }
    else{
      mean_NoFT[i]=Fit_func_Mx.f->GetParameter(1);
      mean_NoFT_Error[i]=Fit_func_Mx.f->GetParError(1);
      sigma_NoFT[i]=Fit_func_Mx.f->GetParameter(2);
      sigma_NoFT_Error[i]=Fit_func_Mx.f->GetParError(2);
      //cout<<"HEre2"<<endl;
    }

  }*/
  // can2->Print((pdf_original + ")").c_str());
}

void plot_results_exc(string pdfname)
{

  TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
  string pdf_original = pdfname + ".pdf";
  gStyle->SetOptFit(1011);
  gStyle->SetStatFontSize(0.020); // 0.015
  double_t par_Jpsi;
  double_t err_Jpsi;

  TPaveText *pt = new TPaveText();
  auto legend = new TLegend(); // 0.54, 0.87, 0.90, 0.60

  // Histograms to print
  Fit_Function Fit_func;
  Fit_func.Set_Limits(2.0, 3.5);

  double nb_JPsi;

  h_Invariant->Draw();
  h_Invariant->SetMinimum(0); // Set the lower limit
  h_Invariant->SetMaximum(21); // Set the upper limit
  Fit_func.Fit_dNdMee(h_Invariant);

  nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);
  par_Jpsi = Fit_func.f->GetParameter(0);
  err_Jpsi = Fit_func.f->GetParError(0);
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.2f +/- %3.2f) ", par_Jpsi, err_Jpsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  h_Mxep->Draw();
  h_Mxep->SetMinimum(0); // Set the lower limit
  h_Mxep->SetMaximum(18); // Set the upper limit
  Fit_func.Set_Limits(2.5, 3.4);
  Fit_func.Fit_dNdMee(h_Mxep);
  nb_JPsi = Fit_func.Get_Integral_Signal(h_Mxep);
  par_Jpsi = Fit_func.f->GetParameter(0);
  err_Jpsi = Fit_func.f->GetParError(0);
  // pt->AddText("Fit_dNdMee");
  // pt->SetTextSize(0.03);
  // pt->Draw();
  legend->SetFillStyle(0);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.2f +/- %3.2f) ", par_Jpsi, err_Jpsi), "l");
  legend->Draw("same");
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  gStyle->SetStatFontSize(0.03);

  can->Divide(1, 5);
  can->cd(1);
  gPad->SetLogy();
  h_vertex_timediff_FT_FD->Draw();
  can->cd(2);
  gPad->SetLogy();
  h_vertex_timediff_FT_eFD->Draw();
  can->cd(3);
  gPad->SetLogy();
  h_vertex_timediff_FT_pFD->Draw();
  can->cd(4);
  gPad->SetLogy();
  h_vertex_timediff_eFD_pFD->Draw();
  can->cd(5);
  gPad->SetLogy();
  h_vertex_timediff_FD_pFD->Draw();

  can->Print((pdf_original + "(").c_str());
  can->Clear();
  gPad->SetLogy(0);

  can->Divide(2, 2);
  can->cd(1);
  h_W->Draw();
  can->cd(2);
  h_W_ep->Draw();
  can->cd(3);
  h_Q2->Draw();
  can->cd(4);
  h_MXepee->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  can->Divide(2, 2);
  can->cd(1);
  h_IMee_vs_Mxep->Draw("colz");
  can->cd(2);
  h2_Invariant->Draw("colz");
  can->cd(3);
  h_MXeee->Draw();
  can->cd(4);
  h2_MXeee->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  // Missing mass
  // Page 6

  can->Divide(2, 2);
  can->cd(1);
  h_Mxep->Draw();
  can->cd(2);
  h_MXepe1->Draw();
  can->cd(3);
  h_MXee1->Draw();
  can->cd(4);
  h_MXpe1->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  can->Divide(2, 2);
  can->cd(1);
  h2_MXep->Draw("colz");
  can->cd(2);
  h2_MXepe1->Draw("colz");
  can->cd(3);
  h2_MXee1->Draw("colz");
  can->cd(4);
  h2_MXpe1->Draw("colz");
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  can->Divide(2, 2);
  can->cd(1);
  h_Mxep->Draw();
  can->cd(2);
  h_MXepe2->Draw();
  can->cd(3);
  h_MXee2->Draw();
  can->cd(4);
  h_MXpe2->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  can->Divide(2, 2);
  can->cd(1);
  h2_MXep->Draw("colz");
  can->cd(2);
  h2_MXepe2->Draw("colz");
  can->cd(3);
  h2_MXee2->Draw("colz");
  can->cd(4);
  h2_MXpe2->Draw("colz");
  can->Print((pdf_original + ")").c_str());
  can->Clear();
}

void plot_results_onelp(string pdfname, int top = 3)
{

  TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
  string pdf_original = pdfname + ".pdf";
  gStyle->SetOptFit(1011);
  gStyle->SetStatFontSize(0.020);
  double_t par_Jpsi;
  double_t err_Jpsi;

  TPaveText *pt = new TPaveText();
  auto legend = new TLegend(); // 0.54, 0.87, 0.90, 0.60

  // Histograms to print
  Fit_Function Fit_func;
  Fit_func.Set_Limits(2.0, 3.5);
  double nb_JPsi;

  cout << "Fitting Fit_dNdMee " << endl;
  h_Mxep->Draw();
  if (top == 3)
  {
    h_Mxep->SetMinimum(0); // Set the lower limit
    h_Mxep->SetMaximum(48.3); // Set the upper limit
    Fit_func.Fit_dNdMee(h_Mxep);
    cout << "after fit" << endl;

    nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);
    par_Jpsi = Fit_func.f->GetParameter(0);
    err_Jpsi = Fit_func.f->GetParError(0);
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.2f +/- %3.2f) ", par_Jpsi, err_Jpsi), "l");
    legend->Draw("same");
  }
  else{
    Fit_func.Fit_dNdMee_pol(h_Mxep);
    cout << "after fit" << endl;

    nb_JPsi = Fit_func.Get_Integral_Signal(h_Invariant);
    par_Jpsi = Fit_func.f->GetParameter(0);
    err_Jpsi = Fit_func.f->GetParError(0);
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(Fit_func.f_signal, Form("J#psi fit (%3.2f +/- %3.2f) ", par_Jpsi, err_Jpsi), "l");
    legend->Draw("same");
  }

  can->Print((pdf_original + "(").c_str());
  can->Clear();
  pt->Clear();
  legend->Clear();

  can->Divide(1, 2);
  can->cd(1);
  gPad->SetLogy();
  h_vertex_timediff_FT_FD->Draw();
  can->cd(2);
  gPad->SetLogy();
  h_vertex_timediff_FT_pFD->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  gPad->SetLogy(0);

  h_W->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();
  h_Q2->Draw();
  can->Print((pdf_original + "(").c_str());
  can->Clear();

  if (top == 3)
  {

    can->Divide(2, 2);
    can->cd(1);
    h_Mxep->Draw();
    can->cd(2);
    h_MXepe1->Draw();
    can->cd(3);
    h_MXee1->Draw();
    can->cd(4);
    h_MXpe1->Draw();
    can->Print((pdf_original + "(").c_str());
    can->Clear();

    can->Divide(2, 2);
    can->cd(1);
    h2_MXep->Draw("colz");
    can->cd(2);
    h2_MXepe1->Draw("colz");
    can->cd(3);
    h2_MXee1->Draw("colz");
    can->cd(4);
    h2_MXpe1->Draw("colz");
    can->Print((pdf_original + ")").c_str());
  }
  else
  {
    can->Divide(2, 2);
    can->cd(1);
    h_Mxep->Draw();
    can->cd(2);
    h_MXepe2->Draw();
    can->cd(3);
    h_MXee2->Draw();
    can->cd(4);
    h_MXpe2->Draw();
    can->Print((pdf_original + "(").c_str());
    can->Clear();

    can->Divide(2, 2);
    can->cd(1);
    h2_MXep->Draw("colz");
    can->cd(2);
    h2_MXepe2->Draw("colz");
    can->cd(3);
    h2_MXee2->Draw("colz");
    can->cd(4);
    h2_MXpe2->Draw("colz");
    can->Print((pdf_original + ")").c_str());
  }
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

void plot2graph(double_t *variable1, double_t *error1, double_t *variable2, double_t *error2, TString name, TString config)
{
  // double Bins[10] ={8.425, 8.775, 8.975, 9.125, 9.33, 9.58, 9.85, 10.1, 10.3, 10.5};
  // double BinsE[10]={0.225, 0.125, 0.075, 0.075, 0.13, 0.12, 0.15, 0.10, 0.10, 0.10};

  // double Bins[4] = {8.625, 9.225, 9.73, 10.3};
  // double BinsE[4] = {0.425, 0.205, 0.27, 0.3};

  cout << "plot2graph " << endl;

  double Bins[6] = {0.5, 1.1, 1.3, 1.5, 1.7, 2.65};
  // double BinsE[6] = {0.5,0.1,0.1,0.1,0.1,0.85};
  double BinsE[6] = {0};

  TGraphErrors *gr = new TGraphErrors(6, Bins, variable1, BinsE, error1);
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(kBlue);
  gr->SetTitle(name + " " + config + "; M(e^{-}e^{+}) ; M_{X} " + name);
  gr->GetXaxis()->SetRangeUser(0.0, 3.5);
  if (name == "Mean")
  {
    gr->GetYaxis()->SetRangeUser(0.8, 1.2);
  }
  else if (name == "Sigma")
  {
    gr->GetYaxis()->SetRangeUser(0.0, 0.2);
  }

  gr->Draw("AEP");

  TGraphErrors *gr2 = new TGraphErrors(6, Bins, variable2, BinsE, error2);
  gr2->SetMarkerStyle(23);
  gr2->SetMarkerColor(kRed);
  gr2->Draw("P");

  auto legend = new TLegend(0.25, 0.10);
  legend->AddEntry(gr, "With FT correction", "p");
  legend->AddEntry(gr2, "No correction", "p");
  legend->Draw("same");

  if (name == "Mean")
  {
    TLine *line = new TLine(0.3, 0.938, 2.9, 0.938);
    line->SetLineColor(kViolet);
    line->Draw("same");
  }
}

#endif