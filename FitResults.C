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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "bib/Utils.h"
#include "bib/Functions.h"
#include "bib/h_Fit.h"
#include "bib/Plot.h"
#include "bib/Basic.h"
#include "bib/Update.h"


// QADB header and namespace
#include "QADB.h"
using namespace QA;

int FitResults(string pdfname="Full_Invariant", double EP=-0.030 , double Q2cut=0.5, double MMcut=0.4) {
// Record start time
  auto start = std::chrono::high_resolution_clock::now();
  //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);


    //RESULTS RESULTS
    for(int j = 0; j < 10; j++){
      ostringstream name;
      name << "Egamma_bin" << j+1;
      ostringstream sstr;
      sstr<<"Egamma_bin" << j+1<<";M(e+e-),GeV ";
      h_invariant_Egamma[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),100,2.5,3.5);;
      name.str("");
      sstr.str("");    
    }

    for(int j = 0; j < 4; j++){
      ostringstream name;
      name << "Egamma_-t_bin" << j+1;
      ostringstream sstr;
      sstr<<"Egamma_-t_bin" << j+1<<";M(e+e-),GeV ";
      h_invariant_Egamma_t[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),100,2.5,3.5);;
      name.str("");
      sstr.str("");    
    }

    for(int j = 0; j < 6; j++){
      ostringstream name;
      name << "Egamma_-t_bin_+1_" << j+1;
      ostringstream sstr;
      sstr<<"Egamma_-t_bin_+1_" << j+1<<";M(e+e-),GeV ";
      h_invariant_Egamma_t_plus[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),100,2.5,3.5);;
      name.str("");
      sstr.str("");    
    }

    for(int j = 0; j < 6; j++){
      ostringstream name;
      name << "Egamma_-t_bin_-1_" << j+1;
      ostringstream sstr;
      sstr<<"Egamma_-t_bin_-1_" << j+1<<";M(e+e-),GeV ";
      h_invariant_Egamma_t_minus[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),100,2.5,3.5);;
      name.str("");
      sstr.str("");    
    }


    

    Get_histos("InvariantMassS19", EP);
    Get_histos("InvariantMassF18in", EP);
    Get_histos("InvariantMassF18out", EP);


    Fit_Egamma_t_hel("Gauss+Pol",0);
    Fit_Egamma_t_hel("CB+Pol",1);
    Fit_Egamma_t_hel("Gauss+Exp",2);
    Fit_Egamma_t_hel("CB+Exp",3);
    Fit_Egamma_t_hel("GaussExp+Exp",4);




    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;
    

  }//




