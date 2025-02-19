#ifndef Basic
#define Basic

bool MM_bool = true;
Double_t MMcut = 0.1;

bool time_bool = true;
Double_t timecut = 2.0;

bool Egamma_bool = true;
Double_t Egammacut = 8.1;

void Get_histos(string name, Double_t score_e_cut = -0.06, Double_t score_p_cut = -0.06)
{
  TFile *file = new TFile(("/lustre24/expphy/volatile/clas12/mtenorio/Root/Results/" + name + ".root").c_str(), "READ");
  TTree *results = (TTree *)file->Get("results");
  Double_t invariantMass, t, Q2, EP, MM, Egamma, score_e, score_p, score_e_6, score_p_6;
  int helicity;

  results->SetMakeClass(1);
  results->SetBranchAddress("invariantMass", &invariantMass);
  results->SetBranchAddress("t", &t);
  results->SetBranchAddress("Q2", &Q2);
  results->SetBranchAddress("EPcut", &EP);
  results->SetBranchAddress("MM", &MM);
  results->SetBranchAddress("Egamma", &Egamma);
  results->SetBranchAddress("score_e", &score_e);
  results->SetBranchAddress("score_p", &score_p);
  results->SetBranchAddress("score_e_6", &score_e_6);
  results->SetBranchAddress("score_p_6", &score_p_6);

  bool Q2_bool = true;
  Double_t Q2cut = 0.5;

  bool EP_bool = false;
  Double_t EPcut = -0.3;

  int count = 0;

  for (int fc = 0; fc < results->GetEntries(); fc++)
  { // Run 5032 to 5419 // 6616 6783
    results->GetEntry(fc);

    if (score_e_6 < score_e_cut)
      continue;
    if (score_p_6 < score_p_cut)
      continue;

    h_Q2->Fill(Q2);

    h_mm_vs_invariantmass->Fill(invariantMass, MM);
    h_qq_vs_invariantmass->Fill(invariantMass, Q2);

    if (Q2_bool && Q2 > Q2cut)
      continue;

    h_MM->Fill(MM);

    if (MM_bool && MM > MMcut)
      continue;
    if (EP_bool && EP < EPcut)
      continue;

    h_Invariant->Fill(invariantMass);
    h_E_photon->Fill(Egamma);

    // Bin_Egamma(Egamma,invariantMass);
    // Bin_Egamma_t(Egamma,t,helicity,invariantMass);
  }
}

void Reset_histos()
{
  h_Invariant->Reset("ICES");
  h_Invariant->Reset("ICES");
  h_IMee_vs_Mxep->Reset("ICES");
  h_Invariant->Reset("ICES");
  h_Mxep->Reset("ICES");
  h_MXeee->Reset("ICES");
  h_MXepe1->Reset("ICES");
  h_MXepe2->Reset("ICES");
  h_MXee1->Reset("ICES");
  h_MXee2->Reset("ICES");
  h_MXpe1->Reset("ICES");
  h_MXpe2->Reset("ICES");

  h2_Invariant->Reset("ICES");
  h2_MXep->Reset("ICES");
  h2_MXeee->Reset("ICES");
  h2_MXepe1->Reset("ICES");
  h2_MXepe2->Reset("ICES");
  h2_MXee1->Reset("ICES");
  h2_MXee2->Reset("ICES");
  h2_MXpe1->Reset("ICES");
  h2_MXpe2->Reset("ICES");

  h_In_M_ep->Reset("ICES");
  h_In_M_posFDp->Reset("ICES");
  h_In_M_eleFDp->Reset("ICES");
  h_mm_vs_invariantmass->Reset("ICES");
  h_qq_vs_invariantmass->Reset("ICES");

  h_Mee->Reset("ICES");
  h_MM->Reset("ICES");
  h_Q2->Reset("ICES");
  h_E_photon->Reset("ICES");

  h_im_onelep->Reset("ICES");

  h_Invariant_noCuts->Reset("ICES");
  h_Invariant_tCut->Reset("ICES");
  h_MX->Reset("ICES");
  h_W->Reset("ICES");
  h_W_ep->Reset("ICES");
  h_mm_vs_invariantmass_t->Reset("ICES");
  h_mm_vs_invariantmass_no_t->Reset("ICES");
  h_dt_vs_invariantmass->Reset("ICES");
  h_MM_vs_dt->Reset("ICES");

  h_vertex_timediff_FT_FD->Reset("ICES");
  h_vertex_timediff_FT_eFD->Reset("ICES");
  h_vertex_timediff_FT_pFD->Reset("ICES");
  h_vertex_timediff_eFD_pFD->Reset("ICES");
  h_vertex_timediff_FD_pFD->Reset("ICES");

  for (int i = 0; i < 6; i++)
    h_Data_Var[i]->Reset("ICES");
}

void Get_histos_t(string name, Double_t score_e_cut = -0.6, Double_t score_p_cut = -0.6)
{
  TFile *file = new TFile(("/lustre24/expphy/volatile/clas12/mtenorio/Root/Results/" + name + ".root").c_str(), "READ");
  TTree *results = (TTree *)file->Get("results");
  Double_t invariantMass, Q2, W, MM, Egamma, score_e, score_p, time_ee, score_e_6, score_p_6;

  results->SetMakeClass(1);
  results->SetBranchAddress("invariantMass", &invariantMass);
  results->SetBranchAddress("Q2", &Q2);
  results->SetBranchAddress("MM", &MM);
  results->SetBranchAddress("time_ee", &time_ee);
  results->SetBranchAddress("Egamma", &Egamma);
  results->SetBranchAddress("W", &W);
  results->SetBranchAddress("score_e", &score_e);
  results->SetBranchAddress("score_p", &score_p);
  results->SetBranchAddress("score_e_6", &score_e_6);
  results->SetBranchAddress("score_p_6", &score_p_6);

  double_t Mx_min, Mx_max;
  if (name.find("FTcorrON") != std::string::npos)
  {
    // Mx_min=0.9593-(0.0942*3);
    // Mx_max=0.9593+(0.0942*3);
    Mx_min = 0.963 - (0.08452 * 3);
    Mx_max = 0.963 + (0.08452 * 3);
  }
  else
  {
    // Mx_min=1.015-(0.09222*3);
    // Mx_max=1.015+(0.09222*3);
    Mx_min = 0.999 - (0.0756 * 3);
    Mx_max = 0.999 + (0.0756 * 3);
  }

  Mx_min = 0.963 - (0.08452 * 3);
  Mx_max = 0.963 + (0.08452 * 3);

  int count = 0;

  for (int fc = 0; fc < results->GetEntries(); fc++)
  {
    results->GetEntry(fc);

    if (score_e < score_e_cut)
      continue;
    if (score_p < score_p_cut)
      continue;

    if (Egamma_bool && Egamma < Egammacut)
    {
      continue;
    }

    h_Q2->Fill(Q2);

    if (time_bool && abs(time_ee) > timecut)
    {
      h_mm_vs_invariantmass_no_t->Fill(invariantMass, MM);
      h_Invariant_tCut->Fill(invariantMass);
      if (MM < Mx_min || MM > Mx_max)
        h_Invariant_noCuts->Fill(invariantMass);
      continue;
    }

    h_mm_vs_invariantmass_t->Fill(invariantMass, MM);

    h_MX->Fill(MM);
    Bin_oneVar(invariantMass, MM);

    if (MM_bool && (MM < Mx_min || MM > Mx_max))
    {
      continue;
    }

    h_Invariant->Fill(invariantMass);

    if (invariantMass >= 3.0 && invariantMass <= 3.2)
    {
      h_W->Fill(W);
    }
  }
}

void Get_histos_exclusive(string name, Double_t score_e_cut = -0.6, Double_t score_p_cut = -0.6)
{
  TFile *file = new TFile(("/lustre24/expphy/volatile/clas12/mtenorio/Root/AI_Dec24/" + name + ".root").c_str(), "READ");
  TTree *results = (TTree *)file->Get("results");
  Double_t invariantMass, Q2, W, Egamma, score_e, score_p, time_ee, time_epos, time_ep;
  Double_t time_efdp, time_posfdp;
  Double_t Mxep, MXepee, MXeee, MXepe1, MXepe2, MXee1, MXee2, MXpe1, MXpe2;
  Double_t scoreBDT, scoreMLP, scoreBDTG; 
  results->SetMakeClass(1);
  results->SetBranchAddress("invariantMass", &invariantMass);
  results->SetBranchAddress("Q2", &Q2);
  results->SetBranchAddress("time_ee", &time_ee);
  results->SetBranchAddress("time_epos", &time_epos);
  results->SetBranchAddress("time_ep", &time_ep);
  results->SetBranchAddress("time_efdp", &time_efdp);
  results->SetBranchAddress("time_posfdp", &time_posfdp);
  results->SetBranchAddress("Egamma", &Egamma);
  results->SetBranchAddress("W", &W);
  results->SetBranchAddress("score_e", &score_e);
  results->SetBranchAddress("score_p", &score_p);
  results->SetBranchAddress("scoreBDT", &scoreBDT);
  results->SetBranchAddress("scoreBDTG", &scoreBDTG);

  results->SetBranchAddress("Mxep", &Mxep);
  results->SetBranchAddress("MXepee", &MXepee);
  results->SetBranchAddress("MXeee", &MXeee);
  results->SetBranchAddress("MXepe1", &MXepe1);
  results->SetBranchAddress("MXepe2", &MXepe2);
  results->SetBranchAddress("MXee1", &MXee1);
  results->SetBranchAddress("MXee2", &MXee2);
  results->SetBranchAddress("MXpe1", &MXpe1);
  results->SetBranchAddress("MXpe2", &MXpe2);

  int count = 0;

  for (int fc = 0; fc < results->GetEntries(); fc++)
  {
    results->GetEntry(fc);

    if (score_e < score_e_cut|| scoreBDT<-0.2)
    {
      continue;
    }
    if (score_p < score_p_cut)
    {
      continue;
    }

    h_vertex_timediff_FT_FD->Fill(time_ee);
    h_vertex_timediff_FT_eFD->Fill(time_epos);
    h_vertex_timediff_FT_pFD->Fill(time_ep);
    h_vertex_timediff_eFD_pFD->Fill(time_posfdp);
    h_vertex_timediff_FD_pFD->Fill(time_efdp);

    if (Egamma_bool && Egamma < Egammacut)
    {
      continue;
    }

    h_Q2->Fill(Q2);

    if (time_bool && (abs(time_ee) > timecut && abs(time_epos) > timecut && abs(time_ep) > timecut))
    {
      continue;
    }

    h_MXepee->Fill(MXepee);
    h_IMee_vs_Mxep->Fill(invariantMass, Mxep);
    if (MM_bool && abs(MXepee) > 0.1)
    {
      continue;
    }

    
    double_t Mx_min, Mx_max;
    Mx_min = 0.963 - (0.08452 * 3);
     Mx_max = 0.963 + (0.08452 * 3);

     if (MM_bool && (MXeee < Mx_min || MXeee > Mx_max))
    {
      //continue;
    }
    

    h_Invariant->Fill(invariantMass);
    h_Invariant->Fill(invariantMass);
    h_Mxep->Fill(Mxep);
    h_MXeee->Fill(MXeee);
    h_MXepe1->Fill(MXepe1);
    h_MXepe2->Fill(MXepe2);
    h_MXee1->Fill(MXee1);
    h_MXee2->Fill(MXee2);
    h_MXpe1->Fill(MXpe1);
    h_MXpe2->Fill(MXpe2);

    h2_Invariant->Fill(invariantMass, Mxep);
    h2_MXep->Fill(Mxep, Mxep);
    h2_MXeee->Fill(MXeee, Mxep);
    h2_MXepe1->Fill(MXepe1, Mxep);
    h2_MXepe2->Fill(MXepe2, Mxep);
    h2_MXee1->Fill(MXee1, Mxep);
    h2_MXee2->Fill(MXee2, Mxep);
    h2_MXpe1->Fill(MXpe1, Mxep);
    h2_MXpe2->Fill(MXpe2, Mxep);

    if (invariantMass >= 3.0 && invariantMass <= 3.2)
    {
      h_W->Fill(W);
    }
    if (Mxep >= 3.0 && Mxep <= 3.2)
    {
      h_W_ep->Fill(W);
    }
  }
}

void Get_histos_onelp(string name, int top=3, Double_t score_e_cut = -0.6, Double_t score_p_cut = -0.6, Double_t AIcut=0.04)
{
  TFile *file = new TFile(("/lustre24/expphy/volatile/clas12/mtenorio/Root/AI_Dec24/" + name + ".root").c_str(), "READ");
  cout << "Reading file " << name << ".root" << endl;
  TTree *results = (TTree *)file->Get("results");
  bool isEle, isPos;
  isEle = false;
  isPos = false;

  Double_t Q2, W, Egamma, score_p, score_e, time_ee, time_ep;
  Double_t Mxep, MXepe, MXee, MXpe;
  Double_t scoreBDT, scoreMLP, scoreBDTG;

  results->SetMakeClass(1);
  results->SetBranchAddress("Q2", &Q2);
  results->SetBranchAddress("time_ee", &time_ee);
  results->SetBranchAddress("time_ep", &time_ep);
  results->SetBranchAddress("Egamma", &Egamma);
  results->SetBranchAddress("W", &W);
  results->SetBranchAddress("Mxep", &Mxep);
  results->SetBranchAddress("MXepe", &MXepe);
  results->SetBranchAddress("MXee", &MXee);
  results->SetBranchAddress("MXpe", &MXpe);
  results->SetBranchAddress("scoreBDT", &scoreBDT);
  results->SetBranchAddress("scoreBDTG", &scoreBDTG);

  if (top==3){
    cout<<"Reading e'p'e+"<<endl;
    isPos = true;
    results->SetBranchAddress("score_p", &score_p);
  }
  else{
    isEle = true;
    results->SetBranchAddress("score_e", &score_e);
  }

  int count = 0;

  for (int fc = 0; fc < results->GetEntries(); fc++)
  {
    results->GetEntry(fc);

    if (top==3 && (score_p < score_p_cut||scoreBDT<-0.2 ))
    {
      continue;
    }

    if (top==4  && (score_e < score_e_cut||scoreBDT<AIcut))//0.0481 
    {
      continue;
    }

    h_vertex_timediff_FT_FD->Fill(time_ee);
    h_vertex_timediff_FT_pFD->Fill(time_ep);

    if (time_bool && (abs(time_ee) > timecut  && abs(time_ep) > timecut))
    {
      continue;
    }

    if (Egamma_bool && Egamma < Egammacut)
    {
      continue;
    }

    h_Q2->Fill(Q2);

    

    if (top==3){
      h_MXepe1->Fill(MXepe);
    }
    else{
      h_MXepe2->Fill(MXepe);
    }

    if (MM_bool && abs(MXepe) > MMcut)
    {
      continue;
    }

/*
    if(MXee<3.45|| MXee>6.47)
      continue;

    if(MXpe<0.18 || MXpe>1.44)
      continue;*/
      

      

    h_Mxep->Fill(Mxep);
    h2_MXep->Fill(Mxep, Mxep);

    if (top==3)
    {
      h2_MXepe1->Fill(MXepe, Mxep);
      h_MXee1->Fill(MXee);
      h_MXpe1->Fill(MXpe);
      h2_MXee1->Fill(MXee, Mxep);
      h2_MXpe1->Fill(MXpe, Mxep);
    }
    else
    {
      h2_MXepe2->Fill(MXepe, Mxep);
      h_MXee2->Fill(MXee);
      h_MXpe2->Fill(MXpe);
      h2_MXee2->Fill(MXee, Mxep);
      h2_MXpe2->Fill(MXpe, Mxep);
    }

    if (Mxep >= 3.0 && Mxep <= 3.2)
    {
      h_W->Fill(W);
    }
  }
}

/*
void Fit_Egamma_all(string pdfname)
{
  TCanvas *can = new TCanvas("can", "canvas", 200, 10, 1000, 700);
  string pdf_original = pdfname + ".pdf";

  gStyle->SetOptFit(1111);

  double min, max;
  min = 2.6;
  max = 3.4;

  double Events[10] = {0};
  double Amplitude[5][10] = {0};
  double Mean[5][10] = {0};
  double Sigma[5][10] = {0};

  double EventsE[10] = {0};
  double AmplitudeE[5][10] = {0};
  double MeanE[5][10] = {0};
  double SigmaE[5][10] = {0};

  TPaveText *pt = new TPaveText(.05, .1, .95, .8);

  double sum;
  int init_amp_fit;
  //-----------------------------------------
  //-----------------------------------------
  // Fitting
  //-----------------------------------------
  //-----------------------------------------
  TF1 *rFit_F[5];
  for (int fit = 0; fit < 5; fit++)
  {
    can->Clear();
    sum = 0;
    can->Divide(3, 2);
    Fit_Function Fit_func;

    TF1 *rFit;
    cout << fitnames[fit] << endl;

    for (int i = 0; i < 4; i++)
    {
      can->cd(i + 1);
      if ((fit == 0 || fit == 1) && i == 0)
      {
        min = 2.6;
        max = 3.2;
      }
      else if ((fit == 0 || fit == 1) && i == 1)
      {
        min = 2.6;
        max = 3.2;
      }
      else
      {
        min = 2.6;
        max = 3.4;
      }
      Fit_func.Set_Limits(min, max);
      h_invariant_Egamma[i]->Draw();
      switch (fit)
      {
      case 0:
        rFit = Fit_func.Fit_Gauss_Pol2(h_invariant_Egamma[i]);
        break;
      case 1:
        rFit = Fit_func.Fit_CrystallBall_Pol2(h_invariant_Egamma[i]);
        break;
      case 2:
        rFit = Fit_func.Fit_Gauss_Exp(h_invariant_Egamma[i]);
        break;
      case 3:
        rFit = Fit_func.Fit_CrystallBall_Exp(h_invariant_Egamma[i]);
        break;
      case 4:
        rFit = Fit_func.Fit_GaussExp(h_invariant_Egamma[i]);
      default:
        break;
      }

      sum = sum + Fit_func.f_signal->Integral(0., 10.) / (h_invariant_Egamma[i]->GetXaxis()->GetBinWidth(2));
      // cout<<rFit->GetParameter(0)<<"and"<<Fit_func.f_signal->Integral(0.,10.)/(h_invariant_Egamma[i]->GetXaxis()->GetBinWidth(2))<<endl;
      // Get parameters
      Events[i] = h_invariant_Egamma[i]->GetEntries();
      Amplitude[fit][i] = rFit->GetParameter(0);
      Mean[fit][i] = rFit->GetParameter(1);
      Sigma[fit][i] = rFit->GetParameter(2);
      // Get error
      AmplitudeE[fit][i] = rFit->GetParError(0);
      MeanE[fit][i] = rFit->GetParError(1);
      SigmaE[fit][i] = rFit->GetParError(2);
      cout << "Bin " << i << endl;
      Fit_func.GetNJpsi(h_invariant_Egamma[i]);
    }
    can->cd(5);
    pt->AddText((fitnames[fit]).c_str());
    // if(fit!=4)
    pt->AddText(Form("%f", sum));

    pt->SetTextSize(0.1);
    pt->Draw();

    can->cd(6);
    h_Invariant->SetTitle("Full Set");
    h_Invariant->Draw();

    switch (fit)
    {
    case 0:
      rFit_F[fit] = Fit_func.Fit_Gauss_Pol2(h_Invariant);
      break;
    case 1:
      rFit_F[fit] = Fit_func.Fit_CrystallBall_Pol2(h_Invariant);
      break;
    case 2:
      rFit_F[fit] = Fit_func.Fit_Gauss_Exp(h_Invariant);
      break;
    case 3:
      rFit_F[fit] = Fit_func.Fit_CrystallBall_Exp(h_Invariant);
      break;
    case 4:
      rFit_F[fit] = Fit_func.Fit_GaussExp(h_Invariant);
    default:
      break;
    }
    cout << "Total " << endl;
    Fit_func.GetNJpsi(h_Invariant);

    can->Print((pdf_original + "(").c_str());
    can->Clear();
    pt->Clear();
  }

  //-----------------------------------------
  //-----------------------------------------
  // RESULTS
  //-----------------------------------------
  //-----------------------------------------
  can->Divide(2, 2);
  can->cd(1);
  plotgraph(Events, EventsE, "Number Events polinomial");

  can->cd(2);
  //plot2graph(Amplitude, AmplitudeE, "N_{J/#psi}");

  can->cd(3);
  plot2graph(Mean, MeanE, "Mean", 0, 1, rFit_F[0]->GetParameter(1), rFit_F[0]->GetParError(1), rFit_F[1]->GetParameter(1), rFit_F[1]->GetParError(1));

  can->cd(4);
  plot2graph(Sigma, SigmaE, "Sigma", 0, 1, rFit_F[0]->GetParameter(2), rFit_F[0]->GetParError(2), rFit_F[1]->GetParameter(2), rFit_F[1]->GetParError(2));

  can->Print((pdf_original + "(").c_str());

  can->Clear();
  can->Divide(2, 2);
  can->cd(1);
  plotgraph(Events, EventsE, "Number Events exponential");

  can->cd(2);
  //plot2graph(Amplitude, AmplitudeE, "N_{J/#psi}", 2, 3);

  can->cd(3);
  plot2graph(Mean, MeanE, "Mean", 2, 3, rFit_F[2]->GetParameter(1), rFit_F[2]->GetParError(1), rFit_F[3]->GetParameter(1), rFit_F[3]->GetParError(1));

  can->cd(4);
  plot2graph(Sigma, SigmaE, "Sigma", 2, 3, rFit_F[2]->GetParameter(2), rFit_F[2]->GetParError(2), rFit_F[3]->GetParameter(2), rFit_F[3]->GetParError(2));

  can->Print((pdf_original + ")").c_str());
}

void Fit_Egamma_t_all(string pdfname)
{
  TCanvas *can = new TCanvas("can", "canvas", 700, 700);
  string pdf_original = pdfname + ".pdf";

  gStyle->SetOptFit(1111);

  double min, max;
  min = 2.6;
  max = 3.4;

  double Events[10] = {0};
  double Amplitude[5][10] = {0};
  double Mean[5][10] = {0};
  double Sigma[5][10] = {0};

  double EventsE[10] = {0};
  double AmplitudeE[5][10] = {0};
  double MeanE[5][10] = {0};
  double SigmaE[5][10] = {0};

  TPaveText *pt = new TPaveText(.05, .1, .95, .8);

  double sum;
  int init_amp_fit;
  //-----------------------------------------
  //-----------------------------------------
  // Fitting
  //-----------------------------------------
  //-----------------------------------------
  TF1 *rFit_F[5];
  for (int fit = 0; fit < 5; fit++)
  {
    can->Clear();
    sum = 0;
    can->Divide(2, 2);
    Fit_Function Fit_func;

    TF1 *rFit;
    cout << fitnames[fit] << endl;

    for (int i = 0; i < 4; i++)
    {
      can->cd(i + 1);
      if ((fit == 0 || fit == 1) && i == 0)
      {
        min = 2.6;
        max = 3.2;
      }
      else if ((fit == 0 || fit == 1) && i == 1)
      {
        min = 2.6;
        max = 3.2;
      }
      else
      {
        min = 2.6;
        max = 3.4;
      }
      Fit_func.Set_Limits(min, max);
      h_invariant_Egamma_t[i]->Draw();
      switch (fit)
      {
      case 0:
        rFit = Fit_func.Fit_Gauss_Pol2(h_invariant_Egamma_t[i]);
        break;
      case 1:
        rFit = Fit_func.Fit_CrystallBall_Pol2(h_invariant_Egamma_t[i]);
        break;
      case 2:
        rFit = Fit_func.Fit_Gauss_Exp(h_invariant_Egamma_t[i]);
        break;
      case 3:
        rFit = Fit_func.Fit_CrystallBall_Exp(h_invariant_Egamma_t[i]);
        break;
      case 4:
        rFit = Fit_func.Fit_GaussExp(h_invariant_Egamma_t[i]);
      default:
        break;
      }

      sum = sum + Fit_func.f_signal->Integral(0., 10.) / (h_invariant_Egamma_t[i]->GetXaxis()->GetBinWidth(2));
      // cout<<rFit->GetParameter(0)<<"and"<<Fit_func.f_signal->Integral(0.,10.)/(h_invariant_Egamma[i]->GetXaxis()->GetBinWidth(2))<<endl;
      // Get parameters
      Events[i] = h_invariant_Egamma_t[i]->GetEntries();
      Amplitude[fit][i] = rFit->GetParameter(0);
      Mean[fit][i] = rFit->GetParameter(1);
      Sigma[fit][i] = rFit->GetParameter(2);
      // Get error
      AmplitudeE[fit][i] = rFit->GetParError(0);
      MeanE[fit][i] = rFit->GetParError(1);
      SigmaE[fit][i] = rFit->GetParError(2);
      cout << "Bin " << i << endl;
      Fit_func.GetNJpsi(h_invariant_Egamma_t[i]);
    }
    if (fit != 4)
      can->Print((pdf_original + "(").c_str());
    else
      can->Print((pdf_original + ")").c_str());
    can->Clear();
  }
}

void Fit_Egamma_t_hel(string pdfname, int fit)
{
  TCanvas *canplus = new TCanvas("canplus", "canvas", 1000, 700);
  TCanvas *canminus = new TCanvas("canminus", "canvas", 1000, 700);
  string pdf_original = pdfname + ".pdf";
  gStyle->SetOptFit(1111);

  double min, max;
  min = 2.6;
  max = 3.4;

  TPaveText *pt = new TPaveText(.05, .1, .95, .8);

  //-----------------------------------------
  // Fitting
  //-----------------------------------------

  canplus->Divide(3, 4);
  // canminus->Divide(3,1);
  Fit_Function Fit_func_p;
  Fit_Function Fit_func_m;

  double BSA_NJ[6];
  double BSA_Fit[6];

  TF1 *rFit;
  cout << fitnames[fit] << endl;

  for (int i = 0; i < 6; i++)
  {
    cout << "Bin " << i << endl;
    canplus->cd(i + 1);
    if ((fit == 0 || fit == 1) && i == 0)
      max = 3.2;
    else if ((fit == 0 || fit == 1) && i == 1)
      max = 3.2;
    else
      max = 3.4;

    Fit_func_p.Set_Limits(min, max);
    Fit_func_m.Set_Limits(min, max);
    h_invariant_Egamma_t_plus[i]->Draw();
    switch (fit)
    {
    case 0:
      rFit = Fit_func_p.Fit_Gauss_Pol2(h_invariant_Egamma_t_plus[i]);
      break;
    case 1:
      rFit = Fit_func_p.Fit_CrystallBall_Pol2(h_invariant_Egamma_t_plus[i]);
      break;
    case 2:
      rFit = Fit_func_p.Fit_Gauss_Exp(h_invariant_Egamma_t_plus[i]);
      break;
    case 3:
      rFit = Fit_func_p.Fit_CrystallBall_Exp(h_invariant_Egamma_t_plus[i]);
      break;
    case 4:
      rFit = Fit_func_p.Fit_GaussExp(h_invariant_Egamma_t_plus[i]);
      break;
    default:
    {
      cout << "Fit function not found/invalid! Cancel Fitting." << endl;
      break;
    }
    }
    // cout<<"Helicity +1: ";
    auto NJ_p = Fit_func_p.GetNJpsi(h_invariant_Egamma_t_plus[i]);

    canplus->cd(i + 7);
    h_invariant_Egamma_t_minus[i]->Draw();
    switch (fit)
    {
    case 0:
      rFit = Fit_func_m.Fit_Gauss_Pol2(h_invariant_Egamma_t_minus[i]);
      break;
    case 1:
      rFit = Fit_func_m.Fit_CrystallBall_Pol2(h_invariant_Egamma_t_minus[i]);
      break;
    case 2:
      rFit = Fit_func_m.Fit_Gauss_Exp(h_invariant_Egamma_t_minus[i]);
      break;
    case 3:
      rFit = Fit_func_m.Fit_CrystallBall_Exp(h_invariant_Egamma_t_minus[i]);
      break;
    case 4:
      rFit = Fit_func_m.Fit_GaussExp(h_invariant_Egamma_t_minus[i]);
      break;
    default:
    {
      cout << "Fit function not found/invalid! Cancel Fitting." << endl;
      break;
    }
    }
    // cout<<"Helicity -1: ";
    auto NJ_m = Fit_func_m.GetNJpsi(h_invariant_Egamma_t_minus[i]);

    auto NJ_fit_p = Fit_func_p.f_signal->Integral(0., 10.) / (h_invariant_Egamma_t_plus[i]->GetXaxis()->GetBinWidth(2));
    auto NJ_fit_m = Fit_func_m.f_signal->Integral(0., 10.) / (h_invariant_Egamma_t_minus[i]->GetXaxis()->GetBinWidth(2));

    BSA_NJ[i] = GetBSA(NJ_p, NJ_m);
    BSA_Fit[i] = GetBSA(NJ_fit_p, NJ_fit_m);
    cout << "BSA from GetNJpsi is:" << BSA_NJ[i] << endl;
    cout << "BSA from Fit is:" << BSA_Fit[i] << endl;
    cout << "------" << endl;
  }
  canplus->Print((pdf_original + "(").c_str());

  canplus->Clear();
  canplus->Divide(2, 2);
  canplus->cd(1);
  plot(BSA_NJ, "NJpsi ", 0);
  canplus->cd(2);
  plot(BSA_NJ, "From NJpsi", 1);

  canplus->cd(3);
  plot(BSA_Fit, "From Fit ", 0);
  canplus->cd(4);
  plot(BSA_Fit, "From Fit ", 1);
  canplus->Print((pdf_original + ")").c_str());
}
*/

#endif
