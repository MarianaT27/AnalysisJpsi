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

namespace fs = std::filesystem;

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
    int sector;
};

int toroot_v3(string nameFile = "PASS1_FALL", int version = -19, double Beam_E = 10.2, int PASS = 1, int max = 0)
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    //********************************
    // FALL 10.6
    // SPRING 10.2
    TLorentzVector beam(0, 0, Beam_E, Beam_E);
    TLorentzVector target(0, 0, 0, 0.938);

    particl electronFD;
    particl positronFD;
    particl electronFT;
    particl photon;
    particl proton;

    Int_t run_number;
    Int_t event_number;
    Int_t helicity;

    // Radiated Photon Energy correction variables
    int electron_photon;
    int positron_photon;
    Double_t electron_photonE;
    Double_t positron_photonE;

    TH1I *h_PID_others = new TH1I("h_PID_others", "", 50, -25, 25);
    TH2F *h_delta_theta_vs_phi_positron = new TH2F("h_delta_theta_vs_phi_positron", "#Delta#phi vs #Delta#theta e^{+}", 200, -10, 10, 150, -30, 30);
    TH2F *h_delta_theta_vs_phi_electron = new TH2F("h_delta_theta_vs_phi_electron", "#Delta#phi v #Delta#theta e^{-}", 200, -10, 10, 150, -30, 30);

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
    int number_others;

    TString root_file = "/lustre24/expphy/volatile/clas12/mtenorio/Root/" + nameFile + ".root";
    TFile *file = new TFile(root_file, "RECREATE");
    TTree *analysis = new TTree("analysis", root_file);

    analysis->Branch("Beam_E", &Beam_E, "Beam_E/d");
    analysis->Branch("run_number", &run_number, "run_number/i");
    analysis->Branch("event_number", &event_number, "event_number/i");
    analysis->Branch("helicity", &helicity, "helicity/I");

    analysis->Branch("number_of_electrons", &number_of_electrons, "number_of_electrons/i");
    analysis->Branch("electronFD", &electronFD);

    analysis->Branch("number_of_positrons", &number_of_positrons, "number_of_positrons/i");
    analysis->Branch("positronFD", &positronFD);

    analysis->Branch("number_protons", &number_protons, "number_protons/i");
    analysis->Branch("proton", &proton);

    analysis->Branch("electron_photon", &electron_photon, "electron_photon/i");
    analysis->Branch("positron_photon", &positron_photon, "positron_photon/i");
    analysis->Branch("electron_photonE", &electron_photonE, "electron_photonE/d");
    analysis->Branch("positron_photonE", &positron_photonE, "positron_photonE/d");

    analysis->Branch("N_FT", &N_FT, "N_FT/i");
    analysis->Branch("electronFT", &electronFT);
    analysis->Branch("R", &R, "R/d");

    // analysis->Branch("score_pos",&score_pos,"score_pos/d");
    // analysis->Branch("score_ele",&score_ele,"score_ele/d");

    //********************************************
    //:::::::::::::::::::TMVA::::::::::::::::::::::
    //********************************************
    // TMVA::Reader *readerTMVA = new TMVA::Reader( "!Color:!Silent" );
    int model = 9;
    // Create a set of variables and declare them to the reader
    Float_t P, Theta, Phi, PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;
    Float_t Nphe;

    if (max == 0)
    {
        if (version == -18)
            max = 170;
        if (version == +18)
            max = 180;
        if (version == -19)
            max = 120;
    }

    // Start
    std::string folder_path;

    if (version == -18)
        folder_path = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/FTJPsi/";
    else if (version == +18)
        folder_path = "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/FTJPsi/";
    else if (version == -19)
        folder_path = "/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/FTJPsi/";
    if (PASS == 0)
        folder_path = "/volatile/clas12/osg/marianat/job_8326/output/";

    int f_count = 0;

    for (const auto &entry : fs::directory_iterator(folder_path))
    {

        char filename1[500];
        if (f_count >= max)
        {
            break;
        }
        f_count++;

        if (entry.path().extension() == ".hipo")
        {
            std::string filename = entry.path().filename().string();
            if (filename == "recon.hipo" || filename == "gemc_denoised.hipo")
            {
                continue; // Skip this iteration
            }
            sprintf(filename1, "%s%s", folder_path.c_str(), filename.c_str());
            // Print or use the full path
            std::cout << filename1 << std::endl;
        }

        cout << "Analysis running on " << filename1 << endl;

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
        hipo::bank SCI(factory.getSchema("REC::Scintillator"));
        hipo::bank TRK(factory.getSchema("REC::Track"));
        hipo::bank FT(factory.getSchema("REC::ForwardTagger"));
        hipo::bank HEADER(factory.getSchema("RUN::config"));

        int counter = 0;
        while (reader.next() == true)
        { // Loops all events

            counter++;

            particl electron;
            particl positron;
            particl electron_FT;
            particl proton;

            reader.read(event);

            event.getStructure(EVENT);
            event.getStructure(PART);
            event.getStructure(FTPART);
            event.getStructure(CALO);
            event.getStructure(SCI);
            event.getStructure(TRK);
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
            number_others = 0;

            int nrows = PART.getRows();

            bool electronAccepted;
            bool positronAccepted;

            // REC::Particle
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
                int sector;
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
                    else
                    {
                        number_others = number_others + 1;
                        //h_PID_others->Fill(pid);
                    }
                }

                if (pid == 11 && abs(status) >= 1000 && abs(status) < 2000)
                { // PID 11 Means Electron, Status Between 1000 and 2000 Means Forward Tagger
                    N_FT = N_FT + 1;
                    TLorentzVector electron_FT_vector;
                    electron_FT_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                    electron_FT.lorentz = electron_FT_vector;
                    electron_FT.vertexinfo.x = FT.getFloat("x", i);
                    electron_FT.vertexinfo.y = FT.getFloat("y", i);
                    electron_FT.vertexinfo.z = FT.getFloat("z", i);
                    electron_FT.status = status;
                    electron_FT.index = index;
                    electron_FT.chi2pid = chi2pid;
                }
            }

            // Check if at least two variables are greater than zero
            int count_non_zero = 0;

            if (number_of_electrons > 0) count_non_zero++;
            if (number_of_positrons > 0) count_non_zero++;
            if (number_protons > 0) count_non_zero++;
            if (N_FT > 0) count_non_zero++;

            // Skip the event if fewer than 2 variables are greater than zero
            if (count_non_zero < 2) continue;

            // FT
            for (int i1 = 0; i1 < FT.getRows(); i1++)
            {
                response ftresponse;
                int detectorID = FT.getByte("detector", i1);
                int pindex = FT.getShort("pindex", i1);
                int layerID = FT.getByte("layer", i1);

                ftresponse.time = FT.getFloat("time", i1);
                ftresponse.path = FT.getFloat("path", i1);

                if (pindex == electron_FT.index && detectorID == 10 && layerID == 1)
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
                ecalresponse.du = CALO.getFloat("du", i1); // add this too
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

            //Store Variables
        


            event_number = en;
            run_number = rn;
            helicity = hel;
            R = sqrt((FT_xFcal * FT_xFcal) + (FT_yFcal * FT_yFcal) + (FT_zFcal * FT_zFcal));
            electron_FT.vtime = FT_vtFcal - (R / 29.9792458);

            electronFD=electron;
            positronFD=positron;
            electronFT=electron_FT;

            electron_photon = 0;
            positron_photon = 0;
            electron_photonE = 0;
            positron_photonE = 0;

            if (number_of_electrons > 0 || number_of_positrons > 0)
            {
                for (int i = 0; i < PART.getRows(); i++)
                {
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    float beta = PART.getFloat("beta", i);
                    if (PART.getInt("pid", i) == 22 && beta >= 0.91 && beta <= 1.09)
                    {
                        TLorentzVector photon_vector;
                        photon_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0 * 0.0));
                        if(number_of_electrons>0)
                            h_delta_theta_vs_phi_electron->Fill(photon_vector.Theta() * 57.2958 - electronFD.lorentz.Theta()* 57.2958, photon_vector.Phi() * 57.2958 - electronFD.lorentz.Phi()* 57.2958);
                        if(number_of_positrons>0)
                            h_delta_theta_vs_phi_positron->Fill(photon_vector.Theta() * 57.2958 - positronFD.lorentz.Theta()* 57.2958, photon_vector.Phi() * 57.2958 - positronFD.lorentz.Phi()* 57.2958);
                        if ((abs(photon_vector.Theta() * 57.2958 - electronFD.lorentz.Theta()* 57.2958) < 0.7)&&number_of_electrons>0)
                        {
                            electron_photonE = electron_photonE + photon_vector.E();
                            electron_photon++;
                        }
                        if ((abs(photon_vector.Theta() * 57.2958 - positronFD.lorentz.Theta()* 57.2958) < 0.7)&&number_of_positrons>0)
                        {
                            positron_photonE = positron_photonE + photon_vector.E();
                            positron_photon++;
                        }
                    }
                }
            }

            
            analysis->Fill();

            //}//SELECT THE ELECTRON-POSITRON PAIR

            //}//If particle
        } // end while

    } // End for "Runs"

    h_delta_theta_vs_phi_electron->Write();
    h_delta_theta_vs_phi_positron->Write();
    file->Write();

    // gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("can", "canvas", 1000, 1500);
    h_PID_others->Draw();

    TString outname = nameFile + "PID.png";
    can->Print(outname);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;
}
