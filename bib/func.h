#ifndef func
#define func

void Bin_invariant(TH1F* h[], double_t Epho, double_t invariantMass){

  if(Epho>=8.2 && Epho<9.05)
    h[0]->Fill(invariantMass);
  else if(Epho<9.46)
    h[1]->Fill(invariantMass);
  else if(Epho<10)
    h[2]->Fill(invariantMass);
  else if(Epho<10.6)
    h[3]->Fill(invariantMass);

}

void Get_histos(TH1F* h_In, TH1F* h_In_bin[],string name, Double_t EPcut=-0.030){
    TFile *file = new TFile(("/lustre19/expphy/volatile/clas12/mtenorio/Root/"+name+".root").c_str(),"READ");
    TTree *results = (TTree*)file->Get("results");
    Double_t invariantMass,t,Q2, EP, MM, Egamma;

    results->SetMakeClass(1);
    results->SetBranchAddress("invariantMass",&invariantMass);
    results->SetBranchAddress("t",&t);
    results->SetBranchAddress("Q2",&Q2);
    results->SetBranchAddress("EPcut",&EP);
    results->SetBranchAddress("MM",&MM);
    results->SetBranchAddress("Egamma",&Egamma);

    bool MM_bool=true;
    Double_t MMcut=0.4;

    bool Q2_bool=true;
    Double_t Q2cut=0.5;

    bool EP_bool=false;
    

    for(int fc=0;fc<results->GetEntries();fc++) {//Run 5032 to 5419 // 6616 6783
      results->GetEntry(fc);
      if(MM_bool&& MM>MMcut)
        continue;
      if(EP_bool && EP<EPcut)
        continue;
      if(Q2_bool && Q2>Q2cut)
        continue;
      
      h_In->Fill(invariantMass);

      Bin_invariant(h_In_bin,Egamma,invariantMass);
    }
}



#endif
