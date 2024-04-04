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

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))


struct cartesian {

  double x;
  double y;
  double z;

};


struct response {

  cartesian pos;
  cartesian tpos;
  double time;
  double energy;
  double path;
  int sector;
  int layer;
  int index;
  int component;
  double du;
  double dv;
  double dw;
  double m2u;
  double m2v;
  double m2w;
  double m3u;
  double m3v;
  double m3w;
  double u;
  double v;
  double w;
  double widthu;
  double widthv;
  double widthw;
  double x;
  double y;
  double z;
  double quality;
  int degree;


};


struct particl {

  TLorentzVector lorentz;
  cartesian vertexinfo;
  map<int,response> responses;
  int index = -1;
  double beta;
  double chi2pid;
  int status;
  int pid;
  double vtime;
  double E;

};
//job,name,job,fc


int toroot_MC(string nameFile="S19_noB", int job=7326, int version=-19, double Beam_E=10.2) {
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    //********************************
    //FALL 10.6
    //SPRING 10.2
    TLorentzVector beam(0,0,Beam_E,Beam_E);
    TLorentzVector target(0,0,0,0.938);
    Int_t run_number;
    Int_t event_number;

    //Electron variables
    Double_t electron_p, electron_theta, electron_phi, electron_px, electron_py, electron_pz;
    Double_t electron_vx, electron_vy, electron_vz, electron_vt;
    Double_t electron_pcal_v, electron_pcal_w;
    Double_t electron_pcal_energy, electron_ecin_energy, electron_ecout_energy; //Electron Energies
    Double_t electron_energy;
    Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
    Double_t electron_sfpcal, electron_sfecin, electron_sfecout;

    //Positron variables
    Double_t positron_p, positron_theta, positron_phi, positron_px, positron_py, positron_pz;
    Double_t positron_vx, positron_vy, positron_vz, positron_vt;
    Double_t positron_pcal_v, positron_pcal_w;
    Double_t positron_pcal_energy, positron_ecin_energy, positron_ecout_energy; //Positron Energies
    Double_t positron_energy;
    Double_t positron_chi2pid;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;

    //Proton varianbles
    Double_t proton_p, proton_theta, proton_phi, proton_px, proton_py, proton_pz, proton_status, proton_beta, proton_chi2pid;
    Double_t proton_vx, proton_vy, proton_vz, proton_vt, proton_pid;
    Double_t proton_energy;

    //Photon variables
    Double_t photon_p, photon_theta, photon_phi, photon_px, photon_py, photon_pz;
    Double_t photon_vx, photon_vy, photon_vz, photon_vt,photon_energy;
    int electron_photon;
    int positron_photon;
    Double_t electron_photonE;
    Double_t positron_photonE;
    TH2F* h_delta_theta_vs_phi_positron = new TH2F("h_delta_theta_vs_phi_positron","delta_theta vs delta_phi positron",200,-10,10,150,-30,30);
    TH2F* h_delta_theta_vs_phi_electron = new TH2F("h_delta_theta_vs_phi_electron","delta_theta vs delta_phi electron",200,-10,10,150,-30,30);


    //Forward Tagger Electron variables
    Double_t FT_E;
    Double_t FT_p;
    Double_t FT_Px;
    Double_t FT_Py;
    Double_t FT_Pz;
    Double_t FT_vt;
    Double_t FT_x,FT_y,FT_z, R;
    Double_t FT_vtC;
    Double_t FT_vtFcal,FT_xFcal,FT_yFcal,FT_zFcal;

    TLorentzVector ele, pos, ele_FT;

     Double_t score_pos,score_ele;

    

    int N_FT;
    int number_of_electrons;
    int number_of_positrons;
    int number_protons;
    int JPSI=0;


    TString root_file = "/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/"+nameFile+"_MC.root";
    TFile *file = new TFile(root_file,"RECREATE");
    TTree *analysis = new TTree("analysis",root_file);

    analysis->Branch("Beam_E",&Beam_E,"Beam_E/d");
    analysis->Branch("run_number",&run_number,"run_number/i");
    analysis->Branch("event_number",&event_number,"event_number/i");

    analysis->Branch("number_of_electrons",&number_of_electrons,"number_of_electrons/i");
    analysis->Branch("electron_p",&electron_p,"electron_p/d");
    analysis->Branch("electron_theta",&electron_theta,"electron_theta/d");
    analysis->Branch("electron_phi",&electron_phi,"electron_phi/d");
    analysis->Branch("electron_px",&electron_px,"electron_px/d");
    analysis->Branch("electron_py",&electron_py,"electron_py/d");
    analysis->Branch("electron_pz",&electron_pz,"electron_pz/d");
    analysis->Branch("electron_vx",&electron_vx,"electron_vx/d");
    analysis->Branch("electron_vy",&electron_vy,"electron_vy/d");
    analysis->Branch("electron_vz",&electron_vz,"electron_vz/d");
    analysis->Branch("electron_vt",&electron_vt,"electron_vt/d");
    analysis->Branch("electron_pcal_v",&electron_pcal_v,"electron_pcal_v/d");
    analysis->Branch("electron_pcal_w",&electron_pcal_w,"electron_pcal_w/d");
    analysis->Branch("electron_pcal_energy",&electron_pcal_energy,"electron_pcal_energy/d");
    analysis->Branch("electron_ecin_energy",&electron_ecin_energy,"electron_ecin_energy/d");
    analysis->Branch("electron_ecout_energy",&electron_ecout_energy,"electron_ecout_energy/d");
    analysis->Branch("electron_energy",&electron_energy,"electron_energy/d");
    analysis->Branch("electron_m2ecin",&electron_m2ecin,"electron_m2ecin/d");
    analysis->Branch("electron_m2ecout",&electron_m2ecout,"electron_m2ecout/d");
    analysis->Branch("electron_m2pcal",&electron_m2pcal,"electron_m2pcal/d");
    analysis->Branch("electron_sfecin",&electron_sfecin,"electron_sfecin/d");
    analysis->Branch("electron_sfpcal",&electron_sfpcal,"electron_sfpcal/d");
    analysis->Branch("electron_sfecout",&electron_sfecout,"electron_sfecout/d");


    analysis->Branch("number_of_positrons",&number_of_positrons,"number_of_positrons/i");    
    analysis->Branch("positron_p",&positron_p,"positron_p/d");
    analysis->Branch("positron_theta",&positron_theta,"positron_theta/d");
    analysis->Branch("positron_phi",&positron_phi,"positron_phi/d");
    analysis->Branch("positron_px",&positron_px,"positron_px/d");
    analysis->Branch("positron_py",&positron_py,"positron_py/d");
    analysis->Branch("positron_pz",&positron_pz,"positron_pz/d");
    analysis->Branch("positron_vx",&positron_vx,"positron_vx/d");
    analysis->Branch("positron_vy",&positron_vy,"positron_vy/d");
    analysis->Branch("positron_vz",&positron_vz,"positron_vz/d");
    analysis->Branch("positron_vt",&positron_vt,"positron_vt/d");
    analysis->Branch("positron_pcal_v",&positron_pcal_v,"positron_pcal_v/d");
    analysis->Branch("positron_pcal_w",&positron_pcal_w,"positron_pcal_w/d");
    analysis->Branch("positron_pcal_energy",&positron_pcal_energy,"positron_pcal_energy/d");
    analysis->Branch("positron_ecin_energy",&positron_ecin_energy,"positron_ecin_energy/d");
    analysis->Branch("positron_ecout_energy",&positron_ecout_energy,"positron_ecout_energy/d");
    analysis->Branch("positron_chi2pid",&positron_chi2pid,"positron_chi2pid/d");
    analysis->Branch("positron_energy",&positron_energy,"positron_energy/d");
    analysis->Branch("positron_m2ecin",&positron_m2ecin,"positron_m2ecin/d");
    analysis->Branch("positron_m2ecout",&positron_m2ecout,"positron_m2ecout/d");
    analysis->Branch("positron_m2pcal",&positron_m2pcal,"positron_m2pcal/d");
    analysis->Branch("positron_sfecin",&positron_sfecin,"positron_sfecin/d");
    analysis->Branch("positron_sfpcal",&positron_sfpcal,"positron_sfpcal/d");
    analysis->Branch("positron_sfecout",&positron_sfecout,"positron_sfecout/d");

    analysis->Branch("number_protons",&number_protons,"number_protons/i");
    analysis->Branch("proton_p",&proton_p,"proton_p/d");
    analysis->Branch("proton_theta",&proton_theta,"proton_theta/d");
    analysis->Branch("proton_phi",&proton_phi,"proton_phi/d");
    analysis->Branch("proton_px",&proton_px,"proton_px/d");
    analysis->Branch("proton_py",&proton_py,"proton_py/d");
    analysis->Branch("proton_pz",&proton_pz,"proton_pz/d");
    analysis->Branch("proton_vx",&proton_vx,"proton_vx/d");
    analysis->Branch("proton_vy",&proton_vy,"proton_vy/d");
    analysis->Branch("proton_vz",&proton_vz,"proton_vz/d");
    analysis->Branch("proton_vt",&proton_vt,"proton_vt/d");
    analysis->Branch("proton_beta",&proton_beta,"proton_beta/d");
    analysis->Branch("proton_chi2pid",&proton_chi2pid,"proton_chi2pid/d");
    analysis->Branch("proton_energy",&proton_energy,"proton_energy/d");

    analysis->Branch("electron_photon",&electron_photon,"electron_photon/i");
    analysis->Branch("positron_photon",&positron_photon,"positron_photon/i");
    analysis->Branch("electron_photonE",&electron_photonE,"electron_photonE/d");
    analysis->Branch("positron_photonE",&positron_photonE,"positron_photonE/d");

    analysis->Branch("N_FT",&N_FT,"N_FT/i");
    analysis->Branch("FT_p",&FT_p,"FT_p/d");
    analysis->Branch("FT_E",&FT_E,"FT_E/d");
    analysis->Branch("FT_Px",&FT_Px,"FT_Px/d");
    analysis->Branch("FT_Py",&FT_Py,"FT_Py/d");
    analysis->Branch("FT_Pz",&FT_Pz,"FT_Pz/d");
    analysis->Branch("FT_xFcal",&FT_xFcal,"FT_xFcal/d");
    analysis->Branch("FT_yFcal",&FT_yFcal,"FT_yFcal/d");
    analysis->Branch("FT_zFcal",&FT_zFcal,"FT_zFcal/d");
    analysis->Branch("R",&R,"R/d");
    analysis->Branch("FT_vtFcal",&FT_vtFcal,"FT_vtFcal/d");
    analysis->Branch("FT_vt",&FT_vt,"FT_vt/d");

    analysis->Branch("score_pos",&score_pos,"score_pos/d");
    analysis->Branch("score_ele",&score_ele,"score_ele/d");
                                                    
                
      //********************************************
    //:::::::::::::::::::TMVA::::::::::::::::::::::
    //********************************************
     TMVA::Reader *readerTMVA = new TMVA::Reader( "!Color:!Silent" );
     int model=9;
    // Create a set of variables and declare them to the reader
     Float_t P, Theta, Phi, PCAL,ECIN,ECOUT;
     Float_t m2PCAL=-1;
     Float_t m2ECIN=-1;
     Float_t m2ECOUT=-1;
     Float_t Nphe;

     
    readerTMVA->AddVariable( "P",&P );
    readerTMVA->AddVariable( "Theta",&Theta);
    readerTMVA->AddVariable( "Phi",&Phi);
     //readerTMVA->AddVariable( "Nphe",&Nphe);
     readerTMVA->AddVariable( "SFPCAL",&PCAL);
     readerTMVA->AddVariable( "SFECIN",&ECIN);
    readerTMVA->AddVariable( "SFECOUT",&ECOUT );
    readerTMVA->AddVariable( "m2PCAL",&m2PCAL);
    readerTMVA->AddVariable( "m2ECIN",&m2ECIN);
    readerTMVA->AddVariable( "m2ECOUT",&m2ECOUT);

    //Book Methods
    TString weightfile_pos; 
    TString weightfile_ele;
    if(version==-19){
        weightfile_ele= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19neg/TMVAClassification_BDT.weights.xml";
        weightfile_pos= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19pos/TMVAClassification_BDT.weights.xml";
    
    }
    if(version==-18){
        weightfile_ele= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18inneg/TMVAClassification_BDT.weights.xml";
        weightfile_pos= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18inpos/TMVAClassification_BDT.weights.xml";
    }
    if(version==+18){
        weightfile_ele= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18outneg/TMVAClassification_BDT.weights.xml";
        weightfile_pos= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18outpos/TMVAClassification_BDT.weights.xml";

    }
    readerTMVA->BookMVA( "BDT pos method", weightfile_pos );
    readerTMVA->BookMVA( "BDT ele method", weightfile_ele );



    //Start
    for(int fc =0; fc<=100; fc++) {//Run 5032 to 5419 // 6616 6783

        char filename1[500];

        sprintf(filename1,"/volatile/clas12/osg/marianat/job_%d/output/%s-%d-%d.hipo",job,nameFile.c_str(),job,fc);


        hipo::reader  reader;
        reader.open(filename1);

        hipo::dictionary  factory;

        reader.readDictionary(factory);

        factory.show();
        hipo::structure  particles;
        hipo::structure  detectors;
        hipo::event      event;

        hipo::bank  dataPART;
        hipo::bank  dataFTPART;
        hipo::bank  dataCALO;
        hipo::bank  dataHEADER;
        hipo::bank  dataFT;
        hipo::bank  dataEVENT;



        hipo::bank EVENT(factory.getSchema("REC::Event"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank FTPART(factory.getSchema("RECFT::Particle"));
        hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
        hipo::bank FT(factory.getSchema("REC::ForwardTagger"));
        hipo::bank HEADER(factory.getSchema("RUN::config"));

        int counter = 0;
        while(reader.next()==true ){//Loops all events

            counter++;

            particl electron;
            particl positron;
            particl electronFT;
            particl photon;
            particl proton;

            reader.read(event);

            event.getStructure(EVENT);
            event.getStructure(PART);
            event.getStructure(FTPART);
            event.getStructure(CALO);
            event.getStructure(FT);
            event.getStructure(HEADER);




            double partNumber = 100;

            int rn = 0;
            int en = 0;


            if(PART.getSize()<1) {
              partNumber = 200;
              continue;
            }



            if(HEADER.getRows()==1) {
                for(int i = 0; i < HEADER.getRows(); i++) {
                    rn = HEADER.getInt("run",i);
                    en = HEADER.getInt("event",i);
                }
            }



            int event_start_time;

            for(int i = 0; i < EVENT.getRows(); i++) {
                event_start_time = EVENT.getFloat("startTime",i);

            }

            number_of_electrons = 0;
            number_of_positrons = 0;
            number_protons = 0;
            N_FT= 0;



            int nrows = PART.getRows();



            bool electronAccepted;
            bool positronAccepted;


            for(int i = 0; i < nrows; i++){
                int   pid = PART.getInt("pid",i);
                int charge = PART.getByte("charge",i);
                float  px = PART.getFloat("px",i);
                float  py = PART.getFloat("py",i);
                float  pz = PART.getFloat("pz",i);
                float  vx = PART.getFloat("vx",i);
                float  vy = PART.getFloat("vy",i);
                float  vz = PART.getFloat("vz",i);
                float  vt = PART.getFloat("vt",i);
                float beta = PART.getFloat("beta",i);
                float chi2pid = PART.getFloat("chi2pid",i);
                int    status = PART.getInt("status",i);
                int    index = i;
                TLorentzVector particle_vector;
                particle_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0005*0.0005));

                //Status Between 2000 and 4000 Means Forward Detector
                if(abs(status)>=2000 && abs(status)<4000) {

                    //PID 11 Means Electron
                    if(pid==11) {
                        number_of_electrons = number_of_electrons + 1;
                        TLorentzVector electron_vector;
                        electron_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0005*0.0005));
                        electron.lorentz = electron_vector;
                        electron.vertexinfo.x = vx;
                        electron.vertexinfo.y = vy;
                        electron.vertexinfo.z = vz;
                        electron.vtime = vt;
                        electron.status = status;
                        electron.index = index;
                        electron.chi2pid = chi2pid;
                    }
                    //PID -11 Means Positron
                    else if(pid==-11) {
                        number_of_positrons = number_of_positrons + 1;
                        TLorentzVector positron_vector;
                        positron_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0005*0.0005));
                        positron.lorentz = positron_vector;
                        positron.vertexinfo.x = vx;
                        positron.vertexinfo.y = vy;
                        positron.vertexinfo.z = vz;
                        positron.vtime = vt;
                        positron.status = status;
                        positron.index = index;
                        positron.chi2pid = chi2pid;

                    }
                    else if(pid==2212) {
                        number_protons = number_protons + 1;
                        TLorentzVector proton_vector;
                        proton_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.938*0.938));
                        proton.lorentz = proton_vector;
                        proton.vertexinfo.x = vx;
                        proton.vertexinfo.y = vy;
                        proton.vertexinfo.z = vz;
                        proton.vtime = vt;
                        proton.status = status;
                        proton.index = index;
                        proton.pid = pid;
                        proton.chi2pid = chi2pid;
                        proton.beta = beta;
                    }
                }


                if(pid==11 && abs(status)>=1000 && abs(status)<2000) {//PID 11 Means Electron, Status Between 1000 and 2000 Means Forward Tagger
                    N_FT = N_FT + 1;
                    TLorentzVector electronFT_vector;
                    electronFT_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0005*0.0005));
                    electronFT.lorentz = electronFT_vector;
                    electronFT.vertexinfo.x = FT.getFloat("x",i);
                    electronFT.vertexinfo.y = FT.getFloat("y",i);
                    electronFT.vertexinfo.z = FT.getFloat("z",i);
                    electronFT.status = status;
                    electronFT.index = index;
                    electronFT.chi2pid = chi2pid;
                }

            }


            //FT
            for(int i1 = 0; i1 < FT.getRows(); i1++){
                response ftresponse;
                int detectorID = FT.getByte("detector",i1);
                int pindex = FT.getShort("pindex",i1);
                int layerID = FT.getByte("layer",i1);

                ftresponse.time = FT.getFloat("time",i1);
                ftresponse.path = FT.getFloat("path",i1);

                if(pindex==electronFT.index && detectorID==10 && layerID==1 ) { //layer 1 works layer 2 is zero THIS WORKS
                    FT_xFcal=FT.getFloat("x",i1);
                    FT_yFcal=FT.getFloat("y",i1);
                    FT_zFcal=FT.getFloat("z",i1);
                    FT_vtFcal=FT.getFloat("time",i1);
                }
            }



            //REC::Calorimeter
            for (int i1 = 0; i1 < CALO.getRows(); i1++) {
                int pindex = CALO.getShort("pindex",i1);
                response ecalresponse;
                ecalresponse.energy = CALO.getFloat("energy",i1);
                ecalresponse.time = CALO.getFloat("time",i1);
                ecalresponse.index = CALO.getFloat("index",i1);
                ecalresponse.sector = CALO.getByte("sector",i1);
                ecalresponse.layer = CALO.getByte("layer",i1);
                ecalresponse.u = CALO.getFloat("lu",i1);
                ecalresponse.v = CALO.getFloat("lv",i1);
                ecalresponse.w = CALO.getFloat("lw",i1);
                ecalresponse.du = CALO.getFloat("du",i1);
                ecalresponse.dv = CALO.getFloat("dv",i1);
                ecalresponse.dw = CALO.getFloat("dw",i1);
                ecalresponse.m2u = CALO.getFloat("m2u",i1);
                ecalresponse.m2v = CALO.getFloat("m2v",i1);
                ecalresponse.m2w = CALO.getFloat("m2w",i1);

                ecalresponse.path = CALO.getFloat("path",i1);


                if(pindex==positron.index && ecalresponse.layer==1) {
                  positron.responses.insert(pair<int,response>(109,ecalresponse));
                }

                if(pindex==positron.index && ecalresponse.layer==4) {
                  positron.responses.insert(pair<int, response>(110,ecalresponse));
                }
                if(pindex==positron.index && ecalresponse.layer==7) {
                  positron.responses.insert(pair<int, response>(111,ecalresponse));
                }

                if(pindex==electron.index && ecalresponse.layer==1) {//PCAL
                  electron.responses.insert(pair<int,response>(109,ecalresponse));
                }
                if(pindex==electron.index && ecalresponse.layer==4) {//ECIN
                  electron.responses.insert(pair<int, response>(110,ecalresponse));
                }
                if(pindex==electron.index && ecalresponse.layer==7) {//ECOUT
                  electron.responses.insert(pair<int, response>(111,ecalresponse));
                }


            }


            if(partNumber==100) {

               //tagged e+e-
                if(number_of_electrons==1 &&number_of_positrons==1) {

                    event_number = en;
                    run_number = rn;

                    electron_p = electron.lorentz.P();
                    electron_theta = electron.lorentz.Theta()*57.2958;
                    electron_phi = electron.lorentz.Phi()*57.2958;
                    electron_px= electron.lorentz.Px();
                    electron_py= electron.lorentz.Py();
                    electron_pz= electron.lorentz.Pz();
                    electron_vx = electron.vertexinfo.x;
                    electron_vy = electron.vertexinfo.y;
                    electron_vz = electron.vertexinfo.z;
                    electron_vt = electron.vtime;
                    electron_pcal_v=electron.responses[109].v;
                    electron_pcal_w=electron.responses[109].w;
                    electron_pcal_energy = electron.responses[109].energy;
                    electron_ecin_energy = electron.responses[110].energy;
                    electron_ecout_energy = electron.responses[111].energy;
                    electron_energy = electron.lorentz.E();
                    electron_m2pcal=(electron.responses[109].m2u+electron.responses[109].m2v+electron.responses[109].m2w)/3;
                    electron_m2ecin=(electron.responses[110].m2u+electron.responses[110].m2v+electron.responses[110].m2w)/3;
                    electron_m2ecout=(electron.responses[111].m2u+electron.responses[111].m2v+electron.responses[111].m2w)/3;
                    electron_sfpcal=electron_pcal_energy/electron_p;
                    electron_sfecin=electron_ecin_energy/electron_p;
                    electron_sfecout=electron_ecout_energy/electron_p;

                    positron_pcal_v=positron.responses[109].v;
                    positron_pcal_w=positron.responses[109].w;
                    positron_pcal_energy = positron.responses[109].energy;
                    positron_ecin_energy = positron.responses[110].energy;
                    positron_ecout_energy = positron.responses[111].energy;
                    positron_vx = positron.vertexinfo.x;
                    positron_vy = positron.vertexinfo.y;
                    positron_vz = positron.vertexinfo.z;
                    positron_vt=positron.vtime;
                    positron_p = positron.lorentz.P();
                    positron_px =positron.lorentz.Px();
                    positron_py =positron.lorentz.Py();
                    positron_pz =positron.lorentz.Pz();
                    positron_theta = positron.lorentz.Theta()*57.2958;
                    positron_phi = positron.lorentz.Phi()*57.2958;
                    positron_energy = positron.lorentz.E();
                    positron_chi2pid= positron.chi2pid;
                    positron_m2pcal=(positron.responses[109].m2u+positron.responses[109].m2v+positron.responses[109].m2w)/3;
                    positron_m2ecin=(positron.responses[110].m2u+positron.responses[110].m2v+positron.responses[110].m2w)/3;
                    positron_m2ecout=(positron.responses[111].m2u+positron.responses[111].m2v+positron.responses[111].m2w)/3;
                    positron_sfpcal=positron_pcal_energy/positron_p;
                    positron_sfecin=positron_ecin_energy/positron_p;
                    positron_sfecout=positron_ecout_energy/positron_p;

                    P=positron_p;
                    Theta=positron_theta/57.2958;
                    Phi=positron_phi/57.2958;
                    PCAL=positron_sfpcal;
                    ECIN=positron_sfecin;
                    ECOUT=positron_sfecout;
                    m2PCAL=positron_m2pcal;
                    m2ECIN=positron_m2ecin;
                    m2ECOUT=positron_m2ecout;
                    score_pos=readerTMVA->EvaluateMVA("BDT pos method");

                    P=electron_p;
                    Theta=electron_theta/57.2958;
                    Phi=electron_phi/57.2958;
                    PCAL=electron_sfpcal;
                    ECIN=electron_sfecin;
                    ECOUT=electron_sfecout;
                    m2PCAL=electron_m2pcal;
                    m2ECIN=electron_m2ecin;
                    m2ECOUT=electron_m2ecout;
                    score_ele=readerTMVA->EvaluateMVA("BDT ele method");


                    proton_p=proton.lorentz.P();
                    proton_theta=proton.lorentz.Theta()*57.2958;
                    proton_phi=proton.lorentz.Phi()*57.2958;
                    proton_px= proton.lorentz.Px();
                    proton_py= proton.lorentz.Py();
                    proton_pz= proton.lorentz.Pz();
                    proton_vx=proton.vertexinfo.x;
                    proton_vy=proton.vertexinfo.y;
                    proton_vz=proton.vertexinfo.z;
                    proton_vt=proton.vtime;
                    proton_beta=proton.beta;
                    proton_chi2pid=proton.chi2pid;
                    proton_energy = proton.lorentz.E();

                    electron_photon=0;
                    positron_photon=0;
                    electron_photonE=0;
                    positron_photonE=0;

                    for(int i = 0; i < PART.getRows(); i++){
                        float  px = PART.getFloat("px",i);
                        float  py = PART.getFloat("py",i);
                        float  pz = PART.getFloat("pz",i);
                        float beta = PART.getFloat("beta",i);
                        if(PART.getInt("pid",i)==22  && beta>=0.91 && beta<=1.09) {    //){// 0.94 to 1.06
                            TLorentzVector photon_vector;
                            photon_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0*0.0));
                            h_delta_theta_vs_phi_electron->Fill(photon_vector.Theta()*57.2958-electron_theta,photon_vector.Phi()*57.2958-electron_phi);
                            h_delta_theta_vs_phi_positron->Fill(photon_vector.Theta()*57.2958-positron_theta,photon_vector.Phi()*57.2958-positron_phi);
                            if((abs(photon_vector.Theta()*57.2958-electron_theta)<0.7)){//Was 0.7
                                electron_photonE=electron_photonE+photon_vector.E();
                                electron_photon++;
                            }
                            else if((abs(photon_vector.Theta()*57.2958-positron_theta)<0.7)){
                                positron_photonE=positron_photonE+photon_vector.E();
                                positron_photon++;
                            }
                        }
                    }
                    
                    FT_E=electronFT.lorentz.E();
                    FT_p=electronFT.lorentz.P();
                    FT_Px=electronFT.lorentz.Px();
                    FT_Py=electronFT.lorentz.Py();
                    FT_Pz=electronFT.lorentz.Pz();
                    R=sqrt((FT_xFcal*FT_xFcal)+(FT_yFcal*FT_yFcal)+(FT_zFcal*FT_zFcal));
                    FT_vt=FT_vtFcal-(R/29.9792458);
                    analysis->Fill();

                }//SELECT THE ELECTRON-POSITRON PAIR
                
            }//If particle
        }//end while
        printf("run = %d\n",fc);
    }//End for "Runs"

  h_delta_theta_vs_phi_electron->Write();
   h_delta_theta_vs_phi_positron->Write();
   file->Write();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;

}




