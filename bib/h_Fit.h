#ifndef h_Fit
#define h_Fit


class Fit_Function
{
public:
  TF1 *f;
  TF1 *f_signal;
  TF1 *f_BG;
  float min_fit;
  float max_fit;
  int bin = 100;
  float range;

  Fit_Function() {}


  void Set_Limits(float in_min_fit, float in_max_fit)
  {
    min_fit = in_min_fit;
    max_fit = in_max_fit;
  }

  void Set_Bin_Range(int in_bin, float in_min, float in_max)
  {
    bin = in_bin;
    range = in_max - in_min;
  }

  double GetNJpsi(TH1 *h)
  {
    double minB, maxB;
    minB = f->GetParameter(1) - 3 * f->GetParameter(2);
    maxB = f->GetParameter(1) + 3 * f->GetParameter(2);
    // cout<<"Min-max:"<<minB<<"-"<<maxB;

    // number of events under the peak (3sigmas)
    double n3 = h->Integral(h->FindFixBin(minB), h->FindFixBin(maxB), "");
    // number of events in the 2.6-minB
    double nbgr = h->Integral(h->FindFixBin(min_fit), h->FindFixBin(minB), "");

    // Nh, Nb, C, Nfb,Cfb as defined in Eq (27)-(31) of "J/Psi analysis steps" 2021 note
    // Nh is the number of events in the 2.6-minB range, we have to take into account bin size
    // Nh should be equal to nbgr
    auto Nh = f_BG->Integral(min_fit, minB) * bin;
    // Nb is the number of events under the J/psi peak (3sigmas), we have to take into account bin size.
    auto Nb = f_BG->Integral(minB, maxB) * bin;
    // C the ratio of Nb/Nh
    auto C = Nb / Nh;
    // Nfb the total #background in the full spectrum
    auto Nfb = f_BG->Integral(min_fit, max_fit) * bin;
    // Scaling factor
    auto Cfb = Nfb / Nh;

    // #of Jpsi events as defined in Eq (32) of note
    auto Njpsi = n3 - (C * Nh);
    // cout<<"Njpsi is:"<<Njpsi<<" and from fit:"<<f_signal->Integral(0.,10.)/(h->GetXaxis()->GetBinWidth(2))<<endl;

    return Njpsi;
  }

  float Get_Integral_Signal(TH1 *h)
  {
    double integral = 0.0;
    integral = (f_signal->Integral(0., 10.)) / (h->GetXaxis()->GetBinWidth(2));
    return integral;
  }

  TF1 *Fit_Gauss(TH1 *h)
  {
    f = new TF1(h->GetName(), "[0]*0.398942*(1.0/90)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
                 min_fit, max_fit);

    int init_amp_fit = (h->GetBinContent(h->FindBin(0.0)) > 0.0) ? h->GetBinContent(h->FindBin(0.0)) : 5;

    f->SetParameters(init_amp_fit, 0.0, 0.01);
    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 10000);
    // f->SetParLimits (1,-0.05,0.05);
    f->SetParLimits(2, 0.00, 0.1);

    f->SetParNames("A", "Mean", "Sigma");

    h->Fit(f, "ME", "same", min_fit, max_fit);

    return f;
  }

  TF1 *Fit_Gauss_IM(TH1 *h)
  {
    f = new TF1(h->GetName(), "[0]*0.398942*(1.0/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
                min_fit, max_fit);

    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 5;

    f->SetParameters(init_amp_fit, 3.096, 0.01);
    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 10000);
    // f->SetParLimits (1,-0.05,0.05);
    // f->SetParLimits (2,0.00,0.1);

    f->SetParNames("A", "Mean", "Sigma");

    h->Fit(f, "ME", "same", min_fit, max_fit);

    return f;
  }

  TF1 *Fit_Gauss_Pol2(TH1 *h)
  {
    f = new TF1(((TString) "background_fit") + h->GetName(), "[0]*0.398942*(1.5/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + [3]*(x-[1])*(x-[1]) - [4]*(x-[1]) + [5]",
                min_fit, max_fit);

    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 1;

    f->SetParameters(init_amp_fit, 3.096, 0.04, 7, 3, 7);
    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100);
    f->SetParLimits(1, 3.01, 3.12);
    // f->SetParLimits (2,0.025,0.15);
    f->SetParLimits(3, 0, 100000);
    f->SetParLimits(5, 1, 100000);

    f->SetParNames("A", "Mean", "Sigma", "a", "b", "c");

    h->Fit(f, "MEQ", "same", min_fit, max_fit);

    double par[6];
    f->GetParameters(par);

    f_BG = new TF1(((TString) "background_fit"), " [1]*(x-[0])*(x-[0]) - [2]*(x-[0]) + [3]", min_fit, max_fit);
    f_signal = new TF1(((TString) "peak"), "[0]*0.398942*(1.5/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
                       min_fit, max_fit);

    f_BG->SetParameters(par[1], par[3], par[4], par[5]);
    f_signal->SetParameters(par[0], par[1], par[2]);

    f_BG->SetLineColor(kBlue);
    f_signal->SetLineColor(kGreen);

    f_BG->Draw("same");
    f_signal->Draw("same");

    return f;
  }

  TF1 *Fit_Gauss_Exp(TH1 *h)
  {

    f = new TF1(((TString) "background_fit") + h->GetName(), "[0]*0.398942*(1.5/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + TMath::Exp([3]+[4]*x)",
                min_fit, max_fit);

    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 5;

    f->SetParameters(init_amp_fit, 3.096, 0.04, 7, -2);

    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 1000);
    f->SetParLimits(1, 3.01, 3.12);
    f->SetParLimits(2, 0.025, 0.15);
    f->SetParLimits(3, 0.0, 1.0E5);
    f->SetParLimits(4, -1.0E5, 0.0);

    f->SetParNames("A", "Mean", "sigma", "X4", "X5");

    h->Fit(f, "MEQ", "same", min_fit, max_fit);

    double par[5];

    f_BG = new TF1(((TString) "background_fit"), " TMath::Exp([0]+[1]*x)", min_fit, max_fit);

    f_signal = new TF1(((TString) "peak"), "[0]*0.398942*(1.5/100)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
                       min_fit, max_fit);

    f->GetParameters(par);
    f_BG->SetParameters(par[3], par[4]);
    f_signal->SetParameters(par[0], par[1], par[2]);

    f_BG->SetLineColor(kBlue);
    f_signal->SetLineColor(kGreen);

    f_BG->Draw("same");
    f_signal->Draw("same");

    return f;
  }

  TF1 *Fit_CrystallBall_Pol2(TH1 *h)
  {
    f = new TF1(((TString) "background_fit") + h->GetName(), "[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4]) + [5]*(x-[1])*(x-[1]) - [6]*(x-[1]) + [7]", min_fit, max_fit);

    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 5;

    // f->SetParameters(init_amp_fit,3.096,0.04,0.85,1.1,7,3,7);
    f->SetParameters(init_amp_fit, 3.096, 0.04, 1.05, 1.1, 7, 3, 7);

    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100);
    f->SetParLimits(1, 3.01, 3.12);
    f->SetParLimits(2, 0.025, 0.15);
    f->SetParLimits(3, 1, 1.25);

    f->SetParLimits(4, 1, 1.2);

    f->SetParLimits(5, 0, 100000);
    // f->SetParLimits (6,0,250);
    f->SetParLimits(7, 1, 100000);

    f->SetParNames("A", "mean", "sigma", "alpha", "n", "a", "b", "c");

    h->Fit(f, "MEQ", "same", min_fit, max_fit);

    double par[8];
    f_signal = new TF1(h->GetName(), "[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4])", min_fit, max_fit);
    f_BG = new TF1(((TString) "background_fit"), " [1]*(x-[0])*(x-[0]) - [2]*(x-[0]) + [3]", min_fit, max_fit);

    f->GetParameters(par);
    f_signal->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    f_BG->SetParameters(par[1], par[5], par[6], par[7]);

    f_BG->SetLineColor(kBlue);
    f_signal->SetLineColor(kGreen);

    f_BG->Draw("same");
    f_signal->Draw("same");
    gPad->Modified();
    gPad->Update();

    TPaveStats *st = (TPaveStats *)h->FindObject("stats");

    st->SetName("mystats");
    h->SetStats(0);
    TList *list = st->GetListOfLines();
    TText *tconst = st->GetLineWith("A");
    list->Remove(tconst);

    return f;
  }

  TF1 *Fit_CrystallBall_Exp(TH1 *h)
  {
    f = new TF1(((TString) "background_fit") + h->GetName(), "[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4]) + TMath::Exp([5]+[6]*x)", min_fit, max_fit);
    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 5;
    f->SetParameters(init_amp_fit, 3.096, 0.04, 1.05, 1.1, 7, -2);
    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100);
    f->SetParLimits(1, 3.01, 3.12);
    f->SetParLimits(2, 0.025, 0.15);
    //
    f->SetParLimits(3, 1, 1.25);
    f->SetParLimits(4, 1, 1.2);

    f->SetParLimits(5, 0.0, 1.0E5);
    f->SetParLimits(6, -1.0E5, 0.0);

    f->SetParNames("A", "mean", "sigma", "alpha", "n", "x6", "x7");

    h->Fit(f, "MEQ", "same", min_fit, max_fit);

    double par[7];
    f_BG = new TF1(((TString) "background_fit"), " TMath::Exp([0]+[1]*x)", min_fit, max_fit);

    f_signal = new TF1(h->GetName(), "[0]*(1.0/100)*crystalball_function(x, [1], [2], [3], [4])", min_fit, max_fit);

    f->GetParameters(par);
    f_BG->SetParameters(par[5], par[6]);
    f_signal->SetParameters(par[0], par[1], par[2], par[3], par[4]);

    f_BG->SetLineColor(kBlue);
    f_signal->SetLineColor(kGreen);

    f_BG->Draw("same");
    f_signal->Draw("same");

    return f;
  }

  TF1 *Fit_GaussExp(TH1 *h)
  {
    f = new TF1(((TString) "background_fit") + h->GetName(), "[0]*gaussexp_function(x, [1], [2], [3]) + TMath::Exp([4]+[5]*x)", min_fit, max_fit);
    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 5;
    f->SetParameters(init_amp_fit, 3.09, 0.04, 0.5, 7, -2);

    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100);
    f->SetParLimits(1, 3.01, 3.10);
    f->SetParLimits(2, 0.025, 0.15);

    f->SetParLimits(3, 0.30, 1.5);

    f->SetParLimits(4, 0.0, 1.0E5);
    f->SetParLimits(5, -1.0E5, 0.0);

    f->SetParNames("A", "mean", "sigma", "k", "x6", "x7");

    h->Fit(f, "MEQ", "same", min_fit, max_fit);

    double par[7];
    f_BG = new TF1(((TString) "background_fit"), " TMath::Exp([0]+[1]*x)", min_fit, max_fit);

    f_signal = new TF1(h->GetName(), "[0]*gaussexp_function(x, [1], [2], [3])", min_fit, max_fit);

    f->GetParameters(par);
    f_BG->SetParameters(par[4], par[5]);
    f_signal->SetParameters(par[0], par[1], par[2], par[3]);

    f_BG->SetLineColor(kBlue);
    f_signal->SetLineColor(kGreen);

    f_BG->Draw("same");
    f_signal->Draw("same");
    /*
    gPad->Modified(); gPad->Update();

    TPaveStats *st = (TPaveStats*)h->FindObject("stats");


    st->SetName("mystats");
    h->SetStats(0);
    TList *list = st->GetListOfLines();
    TText *tconst = st->GetLineWith("A");
    list->Remove(tconst);
    TLatex *myt = new TLatex(0, 0, Form("J/psi = %f",(f_signal->Integral(0.,10.))/(h->GetXaxis()->GetBinWidth(2))));
    //myt ->SetTextFont(10);
    myt ->SetTextSize(0.025);
    //myt ->SetTextColor(kRed);
    list->Add(myt);

    gPad->Modified();
    gPad->Update(); */

    return f;
  }

  TF1 *Fit_dNdMee(TH1 *h)
  {
    f = new TF1(((TString) "background_fit") + h->GetName(), "[0]*dNdMee(x, [1], [2], [3]) + TMath::Exp([4]+[5]*x)", min_fit, max_fit);
    int init_amp_fit = (h->GetBinContent(h->FindBin(3.096)) > 0.0) ? h->GetBinContent(h->FindBin(3.096)) : 5;
    f->SetParameters(init_amp_fit, 3.09, 0.04, 0.3, 7, -2);
    f->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100);
    f->SetParLimits(1, 3.01, 3.10);
    f->SetParLimits(2, 0.025, 0.15);
    f->SetParLimits(3, 0.1, 0.5);

    f->SetParLimits(4, 0.0, 1.0E5);
    f->SetParLimits(5, -1.0E5, 0.0);

    f->SetParNames("A", "mean", "sigma", "br", "x6", "x7");

    h->Fit(f, "MEQ", "same", min_fit, max_fit);

    double par[7];
    f_BG = new TF1(((TString) "background_fit"), " TMath::Exp([0]+[1]*x)", min_fit, max_fit);

    f_signal = new TF1(h->GetName(), "[0]*dNdMee(x, [1], [2], [3])", min_fit, max_fit);

    f->GetParameters(par);
    f_BG->SetParameters(par[4], par[5]);
    f_signal->SetParameters(par[0], par[1], par[2], par[3]);

    f_BG->SetLineColor(kBlue);
    f_signal->SetLineColor(kGreen);

    f_BG->Draw("same");
    f_signal->Draw("same");

    /*gPad->Modified(); gPad->Update();

    TPaveStats *st = (TPaveStats*)h->FindObject("stats");


    st->SetName("mystats");
    h->SetStats(0);
    TList *list = st->GetListOfLines();
    TText *tconst = st->GetLineWith("A");
    list->Remove(tconst);
    TLatex *myt = new TLatex(0, 0, Form("J/psi = %f",(f_signal->Integral(0.,10.))/(h->GetXaxis()->GetBinWidth(2))));
    //myt ->SetTextFont(10);
    myt ->SetTextSize(0.025);
    //myt ->SetTextColor(kRed);
    list->Add(myt);

    gPad->Modified();
    gPad->Update(); */

    return f;
  }
};


#endif
