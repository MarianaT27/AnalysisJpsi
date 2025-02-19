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
#include <filesystem>
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

int Results(int function=0,string filename = "_All_", string outfilename = "",Double_t AIcut=0.04)
{
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
   ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);

  for(int j = 0; j < 6; j++){
      ostringstream name;
      name << "Mee_bin" << j+1;
      ostringstream sstr;
      sstr<<"Mee_bin" << j+1<<";M_{X},GeV ";
      h_Data_Var[j] = new TH1F(name.str().c_str(), sstr.str().c_str(),90,0.5,1.5);
      name.str("");
      sstr.str("");
      
    }

    //LocateFiles(-18,"FTJPsi");
    //LocateFiles(18,"InclusivePositron");
    //LocateFiles(-19,"FTJPsi");

  switch (function)
  {
  case 1:
    Analysis_tagged(filename+"ee+e-","Tagged_Analysis"+outfilename,true,true);
    break;
  case 2:
    Analysis_exclusive(filename+"ee+e-p","All_Analysis"+outfilename,true,true);
    break;
  case 3:
     Analysis_onelp(filename+"epe+_12",3,"epe+_Analysis"+outfilename,AIcut,true,true);
    break;
  case 4:
    Analysis_onelp(filename+"epe-_12",4,"epe-_Analysis"+outfilename,AIcut,true,true);
    break;

  case 5:
    Analysis_exclusive(filename,outfilename,true,false);
    break;

  case 6:
    Analysis_onelp(filename,4,outfilename,true,false);
    break;

  default:
    cout<<"Non-valid option"<<endl;
    break;
  }


  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";

  return 0;

} //
