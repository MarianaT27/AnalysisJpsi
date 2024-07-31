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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "bib/Utils.h"
#include "bib/Functions.h"
#include "bib/h_Fit.h"
#include "bib/Plot.h"
#include "bib/Basic.h"
#include "bib/Update.h"
#include "bib/Analysis.h"


// QADB header and namespace
#include "QADB.h"
using namespace QA;

int Results(string filename = "", string outfilename = "")
{
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  // ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);
  cout<<"----------------------"<<endl;
  cout<<"Tagged Analysis"<<endl;
  cout<<"----------------------"<<endl;
  Analysis_tagged("All","Tagged_Analysis_LeptonID");
  Reset_histos();
  Analysis_tagged();
  /*cout<<"----------------------"<<endl;
  cout<<"Untagged Analysis"<<endl;
  cout<<"----------------------"<<endl;
  Analysis_untagged();*/

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;

} //
