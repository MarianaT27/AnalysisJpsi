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
#include <iguana/algorithms/clas12/FTEnergyCorrection/Algorithm.h>
#include "clas12reader.h"

#include "QADB.h"
using namespace QA;

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))

struct cartesian
{

    double x;
    double y;
    double z;
};

struct response
{

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

struct particl
{

    TLorentzVector lorentz;
    cartesian vertexinfo;
    map<int, response> responses;
    int index = -1;
    double beta;
    double chi2pid;
    int status;
    int pid;
    double vtime;
    double E;
};

bool Accept_positron(Double_t positron_ecin_energy, Double_t positron_ecout_energy, Double_t positron_pcal_energy, Double_t positron_p, Double_t positron_chi2pid, Double_t positron_pcal_v, Double_t positron_pcal_w, Double_t positron_vz, int version)
{
    if (positron_pcal_energy > 0.07 && positron_p > 1.95 && abs(positron_chi2pid) < 5 && positron_pcal_v > 9 && positron_pcal_w > 9)
    { //
        if (((positron_ecin_energy + positron_ecout_energy) / positron_p) >= (0.195 - positron_pcal_energy / positron_p))
        {
            if (version == +18)
            {
                if (positron_vz > -8 && positron_vz < 4)
                    return true;
            }
            else
                return true;
        }
    }
    return false;
}

bool FC_bool(Double_t pcal_v, Double_t pcal_w, Double_t pcal_u)
{
    if (pcal_v > 9 && pcal_w > 9 && pcal_u>9)
    { //
        return true;
    }
    return false;
}

int Hipo_Analysis()
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    Int_t argc = gApplication->Argc();
    char **argv = gApplication->Argv();

    TString output_file = TString(argv[argc - 1]);
    cout << "Outputting on " << output_file << "\n";

    Double_t Beam_E = 10.6;

    //********************************
    // FALL 10.6
    // SPRING 10.2
    TLorentzVector beam(0, 0, Beam_E, Beam_E);
    TLorentzVector target(0, 0, 0, 0.938);
    Int_t run_number;
    Int_t event_number;
    Int_t helicity;

    // Electron variables
    Double_t electron_p, electron_theta, electron_phi, electron_px, electron_py, electron_pz;
    Double_t electron_vx, electron_vy, electron_vz, electron_vt;
    Double_t electron_pcal_v, electron_pcal_w,electron_pcal_u;
    Double_t electron_pcal_energy, electron_ecin_energy, electron_ecout_energy; // Electron Energies
    Double_t electron_energy;
    Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
    Double_t electron_sfpcal, electron_sfecin, electron_sfecout;

    // Positron variables
    Double_t positron_p, positron_theta, positron_phi, positron_px, positron_py, positron_pz;
    Double_t positron_vx, positron_vy, positron_vz, positron_vt;
    Double_t positron_pcal_v, positron_pcal_w, positron_pcal_u;
    Double_t positron_pcal_energy, positron_ecin_energy, positron_ecout_energy; // Positron Energies
    Double_t positron_energy;
    Double_t positron_chi2pid;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;

    // Proton varianbles
    Double_t proton_p, proton_theta, proton_phi, proton_px, proton_py, proton_pz, proton_status, proton_beta, proton_chi2pid;
    Double_t proton_vx, proton_vy, proton_vz, proton_vt, proton_pid;
    Double_t proton_energy;

    // Photon variables
    Double_t photon_p, photon_theta, photon_phi, photon_px, photon_py, photon_pz;
    Double_t photon_vx, photon_vy, photon_vz, photon_vt, photon_energy;
    int electron_photon;
    int positron_photon;
    Double_t electron_photonE;
    Double_t positron_photonE;
    
    TH2F* h_Mxep_vs_Mxepe= new TH2F("h_Mxep_vs_Mxepe","M_{X}(e'p') vs M_{X}(e'p'e^{+});M_{X}(e'p'), GeV;M_{X}(e'p'e^{+}), GeV",80,2.5,4,90,-1,8);
    TH1F* h_Mxep= new TH1F("h_Mxep","M_{X}(e'p');M_{X}(e'p'),GeV ",80,0,4);
    TH1F* h_Mxepe= new TH1F("h_Mxepe","M_{X}(e'p'e^{+}); M_{x}; Counts ",90,-1,8);

    TH1F *h_Mxep_cuts= new TH1F("h_MXep_cuts","M_{X}(e'p');M_{X}(e'p'); Counts ",80,2.5,4);
    TH1F *h_Mxepe_cuts = new TH1F("h_Mxepe_cuts", "M_{X}(e'p'e^{+}); M_{X}(e'p'e^{+}); Counts ", 50, -0.5, 0.5);
    TH1I *h_electrons = new TH1I("h_electrons", "# of electrons; ; Counts ", 5, 0,5);

    TH1F *h_e_P = new TH1F("h_e_P", "P", 100, 0, 14);
    TH1F *h_e_theta = new TH1F("h_e_theta", "#Theta", 100, 0, 1);
    TH1F *h_e_phi = new TH1F("h_e_phi", "#Phi", 100, -3.2, 3.2);
    TH1F *h_e_SFPCAL = new TH1F("h_e_SFPCAL", "SFPCAL", 100, 0, 0.3);
    TH1F *h_e_SFECIN = new TH1F("h_e_SFECIN", "SFECIN", 100, 0, 0.25);
    TH1F *h_e_SFECOUT = new TH1F("h_e_SFECOUT", "SFECOUT", 100, 0, 0.1);
    TH1F *h_e_m2PCAL = new TH1F("h_e_m2PCAL", "m2PCAL", 100, 0, 80);
    TH1F *h_e_m2ECIN = new TH1F("h_e_m2ECIN", "m2ECIN", 100, 0, 350);
    TH1F *h_e_m2ECOUT = new TH1F("h_e_m2ECOUT", "m2ECOUT", 100, 0, 200);

    TH1F *h_e_P_cut = new TH1F("h_e_P_cut", "P", 100, 0, 14);
    TH1F *h_e_theta_cut = new TH1F("h_e_theta_cut", "#Theta", 100, 0, 1);
    TH1F *h_e_phi_cut = new TH1F("h_e_phi_cut", "#Phi", 100, -3.2, 3.2);
    TH1F *h_e_SFPCAL_cut = new TH1F("h_e_SFPCAL_cut", "SFPCAL", 100, 0, 0.3);
    TH1F *h_e_SFECIN_cut = new TH1F("h_e_SFECIN_cut", "SFECIN", 100, 0, 0.25);
    TH1F *h_e_SFECOUT_cut = new TH1F("h_e_SFECOUT_cut", "SFECOUT", 100, 0, 0.1);
    TH1F *h_e_m2PCAL_cut = new TH1F("h_e_m2PCAL_cut", "m2PCAL", 100, 0, 80);
    TH1F *h_e_m2ECIN_cut = new TH1F("h_e_m2ECIN_cut", "m2ECIN", 100, 0, 350);
    TH1F *h_e_m2ECOUT_cut = new TH1F("h_e_m2ECOUT_cut", "m2ECOUT", 100, 0, 200);


    // Forward Tagger Electron variables
    Double_t FT_E;
    Double_t FT_p;
    Double_t FT_Px;
    Double_t FT_Py;
    Double_t FT_Pz;
    Double_t FT_vt;
    Double_t FT_x, FT_y, FT_z, R;
    Double_t FT_vtC;
    Double_t FT_vtFcal, FT_xFcal, FT_yFcal, FT_zFcal;

    TLorentzVector ele, pos, ele_FT;

    Double_t score_pos, score_ele;

    int N_FT;
    int number_of_electrons;
    int number_of_positrons;
    int number_protons;
    int JPSI = 0;


   gSystem->Load("libIguanaAlgorithms");

    iguana::clas12::FTEnergyCorrection algo_correction;

    algo_correction.Start();
    //********************************************
    //:::::::::::::::::::TMVA::::::::::::::::::::::
    //********************************************
    TMVA::Reader *readerTMVA = new TMVA::Reader("!Color:Silent");
    int model = 9;
    // Create a set of variables and declare them to the reader
    Float_t P, Theta, Phi, PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;
    Float_t Nphe;

    //readerTMVA->AddVariable("P", &P);
    //readerTMVA->AddVariable("Theta", &Theta);
    //readerTMVA->AddVariable("Phi", &Phi);
    // readerTMVA->AddVariable( "Nphe",&Nphe);
    readerTMVA->AddVariable("SFPCAL", &PCAL);
    readerTMVA->AddVariable("SFECIN", &ECIN);
    readerTMVA->AddVariable("SFECOUT", &ECOUT);
    readerTMVA->AddVariable("m2PCAL", &m2PCAL);
    readerTMVA->AddVariable("m2ECIN", &m2ECIN);
    readerTMVA->AddVariable("m2ECOUT", &m2ECOUT);

    // Book Methods
    TString weightfile_pos;
    TString weightfile_ele;
    weightfile_ele= "/work/clas12/mtenorio/ML_weights_pass2/S19neg/TMVAClassification_BDT_6.weights.xml";
    weightfile_pos= "/work/clas12/mtenorio/ML_weights_pass2/S19pos/TMVAClassification_BDT_6.weights.xml";
    /* if(version==-19){
         weightfile_ele= "/work/clas12/mtenorio/ML_weights_pass2/S19neg/TMVAClassification_BDT.weights.xml";
         weightfile_pos= "/work/clas12/mtenorio/ML_weights_pass2/S19pos/TMVAClassification_BDT.weights.xml";

     }*/
    // if(version==-18){
    //weightfile_ele = "/work/clas12/mtenorio/ML_weights_pass2/F18inneg/TMVAClassification_BDT_6.weights.xml";
    //weightfile_pos = "/work/clas12/mtenorio/ML_weights_pass2/F18inpos/TMVAClassification_BDT_6.weights.xml";
    //}
    // if(version==+18){
    // weightfile_ele= "/work/clas12/mtenorio/ML_weights_pass2/F18outneg/TMVAClassification_BDT.weights.xml";
    // weightfile_pos= "/work/clas12/mtenorio/ML_weights_pass2/F18outpos/TMVAClassification_BDT.weights.xml";

    //}
    readerTMVA->BookMVA("BDT pos method", weightfile_pos);
    readerTMVA->BookMVA("BDT ele method", weightfile_ele);
    int num_files = 0;

    int extra_electrons=0;
    QADB *qa = new QADB();
    // Start
    for (int fc = 5; fc < (argc - 1); fc++)
    { // Run 5032 to 5419 // 6616 6783
        // if(num_files>10)
        //   break;

        TString filename1 = TString(argv[fc]);
        cout << "Analysis running on " << filename1 << endl;
        num_files++;

        hipo::reader reader;
        reader.open(filename1);

        hipo::dictionary factory;

        reader.readDictionary(factory);

        // factory.show();
        // hipo::structure  particles;
        // hipo::structure  detectors;
        hipo::event event;

        hipo::bank EVENT(factory.getSchema("REC::Event"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank FTPART(factory.getSchema("RECFT::Particle"));
        hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
        hipo::bank FT(factory.getSchema("REC::ForwardTagger"));
        hipo::bank HEADER(factory.getSchema("RUN::config"));

        int counter = 0;
        while (reader.next() == true)
        { // Loops all events

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

            if (PART.getSize() < 1)
            {
                partNumber = 200;
                continue;
            }

            if (HEADER.getRows() == 1)
            {
                for (int i = 0; i < HEADER.getRows(); i++)
                {
                    rn = HEADER.getInt("run", i);
                    en = HEADER.getInt("event", i);
                }
            }

            bool Keep_event = true;
            if (rn == 5442 || rn == 6749)
                continue;

            Keep_event = qa->OkForAsymmetry(rn, en);
            int bad_runs[] = {5610, 5615, 6631, 6757};
            bool Additional_bad_runs = false;
            Additional_bad_runs = (std::find(std::begin(bad_runs), std::end(bad_runs), run_number) != std::end(bad_runs));

            if (!Keep_event || Additional_bad_runs)
                continue;

            int event_start_time;

            int hel = 0;

            for (int i = 0; i < EVENT.getRows(); i++)
            {
                event_start_time = EVENT.getFloat("startTime", i);
                hel = EVENT.getByte("helicity", i);
            }
            // if(hel==0)
            //   continue;

            number_of_electrons = 0;
            number_of_positrons = 0;
            number_protons = 0;
            N_FT = 0;

            int nrows = PART.getRows();

            bool electronAccepted;
            bool positronAccepted;

            for (int i = 0; i < nrows; i++)
            {
                int pid = PART.getInt("pid", i);
                int charge = PART.getByte("charge", i);
                float px = PART.getFloat("px", i);
                float py = PART.getFloat("py", i);
                float pz = PART.getFloat("pz", i);
                float vx = PART.getFloat("vx", i);
                float vy = PART.getFloat("vy", i);
                float vz = PART.getFloat("vz", i);
                float vt = PART.getFloat("vt", i);
                float beta = PART.getFloat("beta", i);
                float chi2pid = PART.getFloat("chi2pid", i);
                int status = PART.getInt("status", i);
                int index = i;
                TLorentzVector particle_vector;
                particle_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));

                // Status Between 2000 and 4000 Means Forward Detector
                if (abs(status) >= 2000 && abs(status) < 4000)
                {

                    // PID 11 Means Electron
                    if (pid == 11)
                    {
                        number_of_electrons = number_of_electrons + 1;
                        TLorentzVector electron_vector;
                        electron_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                        electron.lorentz = electron_vector;
                        electron.vertexinfo.x = vx;
                        electron.vertexinfo.y = vy;
                        electron.vertexinfo.z = vz;
                        electron.vtime = vt;
                        electron.status = status;
                        electron.index = index;
                        electron.chi2pid = chi2pid;
                    }
                    // PID -11 Means Positron
                    else if (pid == -11)
                    {
                        number_of_positrons = number_of_positrons + 1;
                        TLorentzVector positron_vector;
                        positron_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                        positron.lorentz = positron_vector;
                        positron.vertexinfo.x = vx;
                        positron.vertexinfo.y = vy;
                        positron.vertexinfo.z = vz;
                        positron.vtime = vt;
                        positron.status = status;
                        positron.index = index;
                        positron.chi2pid = chi2pid;
                    }
                    else if (pid == 2212)
                    {
                        number_protons = number_protons + 1;
                        TLorentzVector proton_vector;
                        proton_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.938 * 0.938));
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

                if (pid == 11 && abs(status) >= 1000 && abs(status) < 2000)
                { // PID 11 Means Electron, Status Between 1000 and 2000 Means Forward Tagger
                    N_FT = N_FT + 1;
                    TLorentzVector electronFT_vector;
                    electronFT_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                    electronFT.lorentz = electronFT_vector;
                    electronFT.vertexinfo.x = FT.getFloat("x", i);
                    electronFT.vertexinfo.y = FT.getFloat("y", i);
                    electronFT.vertexinfo.z = FT.getFloat("z", i);
                    electronFT.status = status;
                    electronFT.index = index;
                    electronFT.chi2pid = chi2pid;
                }
            }

            // FT
            for (int i1 = 0; i1 < FT.getRows(); i1++)
            {
                response ftresponse;
                int detectorID = FT.getByte("detector", i1);
                int pindex = FT.getShort("pindex", i1);
                int layerID = FT.getByte("layer", i1);

                ftresponse.time = FT.getFloat("time", i1);
                ftresponse.path = FT.getFloat("path", i1);

                if (pindex == electronFT.index && detectorID == 10 && layerID == 1)
                { // layer 1 works layer 2 is zero THIS WORKS
                    FT_xFcal = FT.getFloat("x", i1);
                    FT_yFcal = FT.getFloat("y", i1);
                    FT_zFcal = FT.getFloat("z", i1);
                    FT_vtFcal = FT.getFloat("time", i1);
                }
            }

            // REC::Calorimeter
            for (int i1 = 0; i1 < CALO.getRows(); i1++)
            {
                int pindex = CALO.getShort("pindex", i1);
                response ecalresponse;
                ecalresponse.energy = CALO.getFloat("energy", i1);
                ecalresponse.time = CALO.getFloat("time", i1);
                ecalresponse.index = CALO.getFloat("index", i1);
                ecalresponse.sector = CALO.getByte("sector", i1);
                ecalresponse.layer = CALO.getByte("layer", i1);
                ecalresponse.u = CALO.getFloat("lu", i1);
                ecalresponse.v = CALO.getFloat("lv", i1);
                ecalresponse.w = CALO.getFloat("lw", i1);
                ecalresponse.du = CALO.getFloat("du", i1);
                ecalresponse.dv = CALO.getFloat("dv", i1);
                ecalresponse.dw = CALO.getFloat("dw", i1);
                ecalresponse.m2u = CALO.getFloat("m2u", i1);
                ecalresponse.m2v = CALO.getFloat("m2v", i1);
                ecalresponse.m2w = CALO.getFloat("m2w", i1);

                ecalresponse.path = CALO.getFloat("path", i1);

                if (pindex == positron.index && ecalresponse.layer == 1)
                {
                    positron.responses.insert(pair<int, response>(109, ecalresponse));
                }

                if (pindex == positron.index && ecalresponse.layer == 4)
                {
                    positron.responses.insert(pair<int, response>(110, ecalresponse));
                }
                if (pindex == positron.index && ecalresponse.layer == 7)
                {
                    positron.responses.insert(pair<int, response>(111, ecalresponse));
                }

                if (pindex == electron.index && ecalresponse.layer == 1)
                { // PCAL
                    electron.responses.insert(pair<int, response>(109, ecalresponse));
                }
                if (pindex == electron.index && ecalresponse.layer == 4)
                { // ECIN
                    electron.responses.insert(pair<int, response>(110, ecalresponse));
                }
                if (pindex == electron.index && ecalresponse.layer == 7)
                { // ECOUT
                    electron.responses.insert(pair<int, response>(111, ecalresponse));
                }
            }

            // tagged e+e-
            /*
            if (number_of_positrons == 1 && number_protons == 1 && N_FT == 1)
            {

                event_number = en;
                run_number = rn;

                helicity = hel;

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

                positron_pcal_v = positron.responses[109].v;
                positron_pcal_w = positron.responses[109].w;
                positron_pcal_energy = positron.responses[109].energy;
                positron_ecin_energy = positron.responses[110].energy;
                positron_ecout_energy = positron.responses[111].energy;
                positron_vx = positron.vertexinfo.x;
                positron_vy = positron.vertexinfo.y;
                positron_vz = positron.vertexinfo.z;
                positron_vt = positron.vtime;
                positron_p = positron.lorentz.P();
                positron_px = positron.lorentz.Px();
                positron_py = positron.lorentz.Py();
                positron_pz = positron.lorentz.Pz();
                positron_theta = positron.lorentz.Theta() * 57.2958;
                positron_phi = positron.lorentz.Phi() * 57.2958;
                positron_energy = positron.lorentz.E();
                positron_chi2pid = positron.chi2pid;
                positron_m2pcal = (positron.responses[109].m2u + positron.responses[109].m2v + positron.responses[109].m2w) / 3;
                positron_m2ecin = (positron.responses[110].m2u + positron.responses[110].m2v + positron.responses[110].m2w) / 3;
                positron_m2ecout = (positron.responses[111].m2u + positron.responses[111].m2v + positron.responses[111].m2w) / 3;
                positron_sfpcal = positron_pcal_energy / positron_p;
                positron_sfecin = positron_ecin_energy / positron_p;
                positron_sfecout = positron_ecout_energy / positron_p;


                proton_p = proton.lorentz.P();
                proton_theta = proton.lorentz.Theta() * 57.2958;
                proton_phi = proton.lorentz.Phi() * 57.2958;
                proton_px = proton.lorentz.Px();
                proton_py = proton.lorentz.Py();
                proton_pz = proton.lorentz.Pz();
                proton_vx = proton.vertexinfo.x;
                proton_vy = proton.vertexinfo.y;
                proton_vz = proton.vertexinfo.z;
                proton_vt = proton.vtime;
                proton_beta = proton.beta;
                proton_chi2pid = proton.chi2pid;
                proton_energy = proton.lorentz.E();

                electron_photon = 0;
                positron_photon = 0;
                electron_photonE = 0;
                positron_photonE = 0;

                for (int i = 0; i < PART.getRows(); i++)
                {
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    float beta = PART.getFloat("beta", i);
                    if (PART.getInt("pid", i) == 22 && beta >= 0.91 && beta <= 1.09)
                    { //){// 0.94 to 1.06
                        TLorentzVector photon_vector;
                        photon_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0 * 0.0));
                        if ((abs(photon_vector.Theta() * 57.2958 - positron_theta) < 0.7))
                        {
                            positron_photonE = positron_photonE + photon_vector.E();
                            positron_photon++;
                        }
                    }
                }

                FT_E = electronFT.lorentz.E();
                FT_p = electronFT.lorentz.P();
                FT_Px = electronFT.lorentz.Px();
                FT_Py = electronFT.lorentz.Py();
                FT_Pz = electronFT.lorentz.Pz();
                R = sqrt((FT_xFcal * FT_xFcal) + (FT_yFcal * FT_yFcal) + (FT_zFcal * FT_zFcal));
                FT_vt = FT_vtFcal - (R / 29.9792458);

                bool leptonAccepted = Accept_positron(positron_ecin_energy, positron_ecout_energy, positron_pcal_energy, positron_p, positron_chi2pid, positron_pcal_v, positron_pcal_w, positron_vz, -18);
                bool protonAccepted = false;
                if (proton_beta > 0.1 && proton_chi2pid < 10 && proton_p > 0.4)
                    protonAccepted = true;
                if (leptonAccepted && protonAccepted)
                {
                    PCAL = positron_sfpcal;
                    ECIN = positron_sfecin;
                    ECOUT = positron_sfecout;
                    m2PCAL = positron_m2pcal;
                    m2ECIN = positron_m2ecin;
                    m2ECOUT = positron_m2ecout;
                    score_pos = readerTMVA->EvaluateMVA("BDT pos method");

                    Double_t lepton_E_recon = positron_energy + positron_photonE;
                    Double_t lepton_px_recon = positron_px * (lepton_E_recon / positron_energy);
                    Double_t lepton_py_recon = positron_py * (lepton_E_recon / positron_energy);
                    Double_t lepton_pz_recon = positron_pz * (lepton_E_recon / positron_energy);

                    TLorentzVector lepton(lepton_px_recon, lepton_py_recon, lepton_pz_recon, lepton_E_recon);
                    auto FT_4vect = algo_correction.Transform(FT_Px, FT_Py, FT_Pz, FT_E);

                    // set variables
                    TLorentzVector electronFT_vec(std::get<0>(FT_4vect), std::get<1>(FT_4vect), std::get<2>(FT_4vect), std::get<3>(FT_4vect));
                    TLorentzVector pro(proton_px, proton_py, proton_pz, proton_energy);
                    TLorentzVector miss_eep = (beam + target) - (lepton + electronFT_vec + pro); // Other lepton
                    TLorentzVector miss_ep = (beam + target) - (electronFT_vec + pro);

                    Double_t Mx_epe = miss_eep.M2();
                    Double_t Mx_ep = miss_ep.M();

                    
                    Double_t time_ee = positron_vt - FT_vt;
                    Double_t time_ep = proton_vt - FT_vt;

                    Double_t Q2 = 2 * Beam_E * std::get<3>(FT_4vect) * (1 - cos(electronFT_vec.Theta()));
                    Double_t E_photon = Beam_E - std::get<3>(FT_4vect);
                    Double_t W = sqrt((0.938 * 0.938) + (2 * 0.938 * E_photon) - Q2);

                    
                    extra_electrons=extra_electrons+number_of_electrons;
                    if(number_of_electrons>0)
                        h_electrons->Fill(number_of_electrons);
                    

                    if (score_pos < 0.00)
                        continue;


                    //h_E_photon->Fill(E_photon);
                    if(abs(Mx_ep)<2.5)
                        continue;


                    h_Mxep_vs_Mxepe->Fill(Mx_ep, Mx_epe);
                    h_Mxep->Fill(Mx_ep);
                    h_Mxepe->Fill(Mx_epe);

                    if (abs(time_ee) > 2.5)
                        continue;

                    h_Mxepe_cuts->Fill(Mx_epe);

                    if(abs(Mx_epe)>0.1)
                        continue;

                    h_Mxep_cuts->Fill(Mx_ep);
                    cout<<"Run: "<<run_number<<" Event: "<<event_number<<endl;
                }

                // analysis->Fill();

            } // SELECT THE ELECTRON-POSITRON PAIR


            */

            
            if (number_of_positrons == 1)
            {
                extra_electrons++;

                event_number = en;
                run_number = rn;

                helicity = hel;

                /*electron_p = electron.lorentz.P();
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
                electron_pcal_u=electron.responses[109].u;
                electron_pcal_energy = electron.responses[109].energy;
                electron_ecin_energy = electron.responses[110].energy;
                electron_ecout_energy = electron.responses[111].energy;
                electron_energy = electron.lorentz.E();
                electron_m2pcal=(electron.responses[109].m2u+electron.responses[109].m2v+electron.responses[109].m2w)/3;
                electron_m2ecin=(electron.responses[110].m2u+electron.responses[110].m2v+electron.responses[110].m2w)/3;
                electron_m2ecout=(electron.responses[111].m2u+electron.responses[111].m2v+electron.responses[111].m2w)/3;
                electron_sfpcal=electron_pcal_energy/electron_p;
                electron_sfecin=electron_ecin_energy/electron_p;
                electron_sfecout=electron_ecout_energy/electron_p;*/

                positron_pcal_v = positron.responses[109].v;
                positron_pcal_w = positron.responses[109].w;
                positron_pcal_u = positron.responses[109].u;
                positron_pcal_energy = positron.responses[109].energy;
                positron_ecin_energy = positron.responses[110].energy;
                positron_ecout_energy = positron.responses[111].energy;
                positron_vx = positron.vertexinfo.x;
                positron_vy = positron.vertexinfo.y;
                positron_vz = positron.vertexinfo.z;
                positron_vt = positron.vtime;
                positron_p = positron.lorentz.P();
                positron_px = positron.lorentz.Px();
                positron_py = positron.lorentz.Py();
                positron_pz = positron.lorentz.Pz();
                positron_theta = positron.lorentz.Theta();// * 57.2958;
                positron_phi = positron.lorentz.Phi();// * 57.2958;
                positron_energy = positron.lorentz.E();
                positron_chi2pid = positron.chi2pid;
                positron_m2pcal = (positron.responses[109].m2u + positron.responses[109].m2v + positron.responses[109].m2w) / 3;
                positron_m2ecin = (positron.responses[110].m2u + positron.responses[110].m2v + positron.responses[110].m2w) / 3;
                positron_m2ecout = (positron.responses[111].m2u + positron.responses[111].m2v + positron.responses[111].m2w) / 3;
                positron_sfpcal = positron_pcal_energy / positron_p;
                positron_sfecin = positron_ecin_energy / positron_p;
                positron_sfecout = positron_ecout_energy / positron_p;


                bool FC_p=FC_bool(positron_pcal_v, positron_pcal_w,positron_pcal_u);
                bool FC_e=FC_bool(electron_pcal_v, electron_pcal_w,electron_pcal_u);

                PCAL = positron_sfpcal;
                ECIN = positron_sfecin;
                ECOUT = positron_sfecout;
                m2PCAL = positron_m2pcal;
                m2ECIN = positron_m2ecin;
                m2ECOUT = positron_m2ecout;
                score_pos = readerTMVA->EvaluateMVA("BDT pos method");

                h_e_P->Fill(positron_p);
                h_e_theta->Fill(positron_theta);
                h_e_phi->Fill(positron_phi);
                h_e_SFPCAL->Fill(PCAL);
                h_e_SFECIN->Fill(ECIN);
                h_e_SFECOUT->Fill(ECOUT);
                h_e_m2PCAL->Fill(m2PCAL);
                h_e_m2ECIN->Fill(m2ECIN);
                h_e_m2ECOUT->Fill(m2ECOUT);

                if(score_pos>0.0){
                    h_e_P_cut->Fill(positron_p);
                    h_e_theta_cut->Fill(positron_theta);
                    h_e_phi_cut->Fill(positron_phi);
                    h_e_SFPCAL_cut->Fill(PCAL);
                    h_e_SFECIN_cut->Fill(ECIN);
                    h_e_SFECOUT_cut->Fill(ECOUT);
                    h_e_m2PCAL_cut->Fill(m2PCAL);
                    h_e_m2ECIN_cut->Fill(m2ECIN);
                    h_e_m2ECOUT_cut->Fill(m2ECOUT);
                }
            } // SELECT THE ELECTRON-POSITRON PAIR


            


        } // end while
    } // End for "Runs"
    algo_correction.Stop();

    cout<<"extra_electrons: "<<extra_electrons<<endl;

    // gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
    /*can->Divide(2, 3);
    can->cd(1);
    h_Mxep_vs_Mxepe->Draw("colz");
    can->cd(2);
    gPad->SetLogy();
    h_electrons->Draw();
    can->cd(3);
    gPad->SetLogy(0);
    h_Mxepe->Draw();
    can->cd(4);
    h_Mxep->Draw();
    can->cd(5);
    h_Mxepe_cuts->Draw();
    can->cd(6);
    h_Mxep_cuts->Draw();*/
    can->Clear();
    can->Divide(3, 3);
    can->cd(1);

    h_e_P->SetLineColor(kRed);
    h_e_P->SetTitle("P");
    h_e_P->Draw();
    h_e_P_cut->Draw("same");

    can->cd(2);
    h_e_theta->SetLineColor(kRed);
    h_e_theta->SetTitle("Theta");
    h_e_theta->Draw();
    h_e_theta_cut->Draw("same");

    can->cd(3);
    h_e_phi->SetLineColor(kRed);
    h_e_phi->SetTitle("Phi");
    h_e_phi->Draw();
    h_e_phi_cut->Draw("same");

    can->cd(4);
    h_e_SFPCAL->SetLineColor(kRed);
    h_e_SFPCAL->SetTitle("SFPCAL");
    h_e_SFPCAL->Draw();
    h_e_SFPCAL_cut->Draw("same");

    can->cd(5);
    h_e_SFECIN->SetLineColor(kRed);
    h_e_SFECIN->SetTitle("SFECIN");
    h_e_SFECIN->Draw();
    h_e_SFECIN_cut->Draw("same");

    can->cd(6);
    h_e_SFECOUT->SetLineColor(kRed);
    h_e_SFECOUT->SetTitle("SFECOUT");
    h_e_SFECOUT->Draw();
    h_e_SFECOUT_cut->Draw("same");

    can->cd(7);
    h_e_m2PCAL->SetLineColor(kRed);
    h_e_m2PCAL->SetTitle("m2PCAL");
    h_e_m2PCAL->Draw();
    h_e_m2PCAL_cut->Draw("same");

    can->cd(8);
    h_e_m2ECIN->SetLineColor(kRed);
    h_e_m2ECIN->SetTitle("m2ECIN");
    h_e_m2ECIN->Draw();
    h_e_m2ECIN_cut->Draw("same");

    can->cd(9);
    h_e_m2ECOUT->SetLineColor(kRed);
    h_e_m2ECOUT->SetTitle("m2ECOUT");
    h_e_m2ECOUT->Draw();
    h_e_m2ECOUT_cut->Draw("same");

    //can->Print((pdf_original + ")").c_str());
    //can->Clear();

    TString outname = output_file + ".png";
    can->Print(outname);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    gApplication->Terminate();

    return 0;
}
