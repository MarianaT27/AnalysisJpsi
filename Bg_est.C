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

int Bg_est()
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
    Double_t lepton1_p, lepton1_theta, lepton1_phi, lepton1_px, lepton1_py, lepton1_pz;
    Double_t lepton1_pcal_energy, lepton1_ecin_energy, lepton1_ecout_energy; // Electron Energies
    Double_t lepton1_energy;
    Double_t lepton1_sfpcal, lepton1_sfecin, lepton1_sfecout;

    // Positron variables
    Double_t lepton2_p, lepton2_theta, lepton2_phi, lepton2_px, lepton2_py, lepton2_pz;
    Double_t lepton2_pcal_energy, lepton2_ecin_energy, lepton2_ecout_energy; // Positron Energies
    Double_t lepton2_energy;
    Double_t lepton2_sfpcal, lepton2_sfecin, lepton2_sfecout;

    // Proton varianbles
    Double_t proton_p, proton_theta, proton_phi, proton_px, proton_py, proton_pz, proton_status, proton_beta, proton_chi2pid;
    Double_t proton_vx, proton_vy, proton_vz, proton_vt, proton_pid;
    Double_t proton_energy;

    // Photon variables
    Double_t photon_p, photon_theta, photon_phi, photon_px, photon_py, photon_pz;
    Double_t photon_vx, photon_vy, photon_vz, photon_vt, photon_energy;
    int lepton1_photon;
    int lepton2_photon;
    Double_t lepton1_photonE;
    Double_t lepton2_photonE;

    TH1F *h1_Mee1_before = new TH1F("h1_Mee1_before", "eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee2_before = new TH1F("h1_Mee2_before", "eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee1_case1 = new TH1F("h1_Mee1_case1", "pt/p<0.05 eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee2_case1 = new TH1F("h1_Mee2_case1", "pt/p<0.05 eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee1_case2 = new TH1F("h1_Mee1_case2", "pt/p<0.1 eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee2_case2 = new TH1F("h1_Mee2_case2", "pt/p<0.1 eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee1_case3 = new TH1F("h1_Mee1_case3", "Q2<0.2 eep; M(ee),GeV;", 100, 0, 4);
    TH1F *h1_Mee2_case3 = new TH1F("h1_Mee2_case3", "Q2<0.2eep; M(ee),GeV;", 100, 0, 4);
    
    TH2F *h2_Ptp_Mx2 = new TH2F("h2_Ptp_Mx2", "Pt_{X}/P_{X} vs MM^2; ;", 100, -3, 3, 100, 0, 0.5);

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
    int number_protonsCD;
    int JPSI = 0;
    int num_files = 0;

    // Start
    for (int fc = 5; fc < (argc - 1); fc++)
    {

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

            particl electron[3];
            particl positron;
            particl electronFT;
            particl photon;
            particl proton;
            particl protonCD;

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
            number_protonsCD = 0;
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
                        if (number_of_electrons < 3)
                        {
                            TLorentzVector electron_vector;
                            electron_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                            electron[number_of_electrons].lorentz = electron_vector;
                            electron[number_of_electrons].vertexinfo.x = vx;
                            electron[number_of_electrons].vertexinfo.y = vy;
                            electron[number_of_electrons].vertexinfo.z = vz;
                            electron[number_of_electrons].vtime = vt;
                            electron[number_of_electrons].status = status;
                            electron[number_of_electrons].index = index;
                            electron[number_of_electrons].chi2pid = chi2pid;
                        }
                        number_of_electrons = number_of_electrons + 1;
                    }
                    // PID -11 Means Positron
                    else if (pid == -11)
                    {
                        number_of_positrons = number_of_positrons + 1;
                        TLorentzVector lepton2_vector;
                        lepton2_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                        positron.lorentz = lepton2_vector;
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

                if (pid == 2212 && abs(status) >= 4000 && abs(status) < 8000)
                { // PID 2212 Means Proton, Status Between 4000 and 8000 Means Central Detector
                    number_protonsCD = number_protonsCD + 1;
                    TLorentzVector protonCD_vector;
                    protonCD_vector.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.938 * 0.938));
                    protonCD.lorentz = protonCD_vector;
                    protonCD.vertexinfo.x = vx;
                    protonCD.vertexinfo.y = vy;
                    protonCD.vertexinfo.z = vz;
                    protonCD.vtime = vt;
                    protonCD.status = status;
                    protonCD.index = index;
                    protonCD.pid = pid;
                    protonCD.chi2pid = chi2pid;
                    protonCD.beta = beta;
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

                for (int el = 0; el < 3; el++)
                {
                    if (pindex == electron[el].index && ecalresponse.layer == 1)
                    { // PCAL
                        electron[el].responses.insert(pair<int, response>(109, ecalresponse));
                    }
                    if (pindex == electron[el].index && ecalresponse.layer == 4)
                    { // ECIN
                        electron[el].responses.insert(pair<int, response>(110, ecalresponse));
                    }
                    if (pindex == electron[el].index && ecalresponse.layer == 7)
                    { // ECOUT
                        electron[el].responses.insert(pair<int, response>(111, ecalresponse));
                    }
                }
            }

            TLorentzVector pro;
            if (number_protons == 1)
            {
                pro=proton.lorentz;
                proton_p = proton.lorentz.P();
                proton_chi2pid = proton.chi2pid;
                proton_energy = proton.lorentz.E();
                if (abs(proton_chi2pid) > 5 || proton_p < 0.4)
                    continue;
            }
            else if (number_protons == 0 && number_protonsCD == 1)
            {
                pro=protonCD.lorentz;
                proton_p = protonCD.lorentz.P();
                proton_chi2pid = protonCD.chi2pid;
                proton_energy = protonCD.lorentz.E();
                if (abs(proton_chi2pid) > 4)
                    continue;
            }
            else
            {
                continue;
            }

            if (number_of_electrons > 3)
                cout << number_of_electrons << " e- at run: " << rn << endl;

            // e+e-
            if (number_of_electrons == 1 && number_of_positrons == 1)
            {
                lepton1_sfpcal = electron[0].responses[109].energy / electron[0].lorentz.P();
                lepton1_sfecin = electron[0].responses[110].energy / electron[0].lorentz.P();

                if (electron[0].lorentz.P() < 1.9 || electron[0].lorentz.P() > 10.6 || positron.lorentz.P() < 1.9)
                    continue;
                if (electron[0].chi2pid < -3 || positron.chi2pid < -3)
                    continue;
                if (electron[0].responses[109].energy < 0.07 || positron.responses[109].energy < 0.07)
                    continue;
                if (electron[0].responses[109].v < 9 || electron[0].responses[109].w < 9)
                    continue;
                if (positron.responses[109].v < 9 || positron.responses[109].w < 9)
                    continue;
                if ((0.9 * lepton1_sfpcal + lepton1_sfecin) < 0.17 || (electron[0].responses[110].energy == 0 && electron[0].responses[109].energy < 0.15))
                    continue;
                if ((0.9 * positron.responses[109].energy + positron.responses[110].energy) < 0.17 || (positron.responses[110].energy == 0 && positron.responses[109].energy < 0.15))
                    continue;

                lepton1_p = electron[0].lorentz.P();
                lepton1_px = electron[0].lorentz.Px();
                lepton1_py = electron[0].lorentz.Py();
                lepton1_pz = electron[0].lorentz.Pz();
                lepton1_theta = electron[0].lorentz.Theta() * 57.2958;
                lepton1_phi = electron[0].lorentz.Phi() * 57.2958;
                lepton1_energy = electron[0].lorentz.E();

                lepton2_p = positron.lorentz.P();
                lepton2_px = positron.lorentz.Px();
                lepton2_py = positron.lorentz.Py();
                lepton2_pz = positron.lorentz.Pz();
                lepton2_theta = positron.lorentz.Theta() * 57.2958;
                lepton2_phi = positron.lorentz.Phi() * 57.2958;
                lepton2_energy = positron.lorentz.E();
            }
            // e-e-
            else if (number_of_electrons == 2)
            {
                lepton1_sfpcal = electron[0].responses[109].energy / electron[0].lorentz.P();
                lepton1_sfecin = electron[0].responses[110].energy / electron[0].lorentz.P();

                lepton2_sfpcal = electron[1].responses[109].energy / electron[1].lorentz.P();
                lepton2_sfecin = electron[1].responses[110].energy / electron[1].lorentz.P();

                if (electron[0].lorentz.P() < 1.9 || electron[0].lorentz.P() > 10.6 ||electron[1].lorentz.P() < 1.9 || electron[1].lorentz.P() > 10.6 )
                    continue;
                if (electron[0].chi2pid < -3 || electron[1].chi2pid < -3 )
                    continue;
                if (electron[0].responses[109].energy < 0.07 || electron[1].responses[109].energy < 0.07)
                    continue;
                if (electron[0].responses[109].v < 9 || electron[0].responses[109].w < 9)
                    continue;
                if (electron[1].responses[109].v < 9 || electron[1].responses[109].w < 9)
                    continue;
                if ((0.9 * lepton1_sfpcal + lepton1_sfecin) < 0.17 || (electron[0].responses[110].energy == 0 && electron[0].responses[109].energy < 0.15))
                    continue;
                if ((0.9 * lepton2_sfpcal + lepton2_sfecin) < 0.17 || (electron[1].responses[110].energy == 0 && electron[1].responses[109].energy < 0.15))
                    continue;

                lepton1_p = electron[0].lorentz.P();
                lepton1_px = electron[0].lorentz.Px();
                lepton1_py = electron[0].lorentz.Py();
                lepton1_pz = electron[0].lorentz.Pz();
                lepton1_theta = electron[0].lorentz.Theta() * 57.2958;
                lepton1_phi = electron[0].lorentz.Phi() * 57.2958;
                lepton1_energy = electron[0].lorentz.E();

                lepton2_p = electron[1].lorentz.P();
                lepton2_px = electron[1].lorentz.Px();
                lepton2_py = electron[1].lorentz.Py();
                lepton2_pz = electron[1].lorentz.Pz();
                lepton2_theta = electron[1].lorentz.Theta() * 57.2958;
                lepton2_phi = electron[1].lorentz.Phi() * 57.2958;
                lepton2_energy = electron[1].lorentz.E();
            }
            else
                continue;

            lepton1_photon = 0;
            lepton2_photon = 0;
            lepton1_photonE = 0;
            lepton2_photonE = 0;

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

                    if ((abs(photon_vector.Theta() * 57.2958 - lepton1_theta) < 0.7))
                    {
                        lepton1_photonE = lepton1_photonE + photon_vector.E();
                        lepton1_photon++;
                    }
                    else if ((abs(photon_vector.Theta() * 57.2958 - lepton2_theta) < 0.7))
                    {
                        lepton2_photonE = lepton2_photonE + photon_vector.E();
                        lepton2_photon++;
                    }
                }
            }
            //------ENERGY CORRECTION ELECTRON------------
            auto lepton1_E_recon = lepton1_energy + lepton1_photonE;
            auto lepton1_px_recon = lepton1_px * (lepton1_E_recon / lepton1_energy);
            auto lepton1_py_recon = lepton1_py * (lepton1_E_recon / lepton1_energy);
            auto lepton1_pz_recon = lepton1_pz * (lepton1_E_recon / lepton1_energy);
            TLorentzVector e1(lepton1_px_recon, lepton1_py_recon, lepton1_pz_recon, lepton1_E_recon);

            //------ENERGY CORRECTION POSITRON------------
            auto lepton2_E_recon = lepton2_energy + lepton2_photonE;
            auto lepton2_px_recon = lepton2_px * (lepton2_E_recon / lepton2_energy);
            auto lepton2_py_recon = lepton2_py * (lepton2_E_recon / lepton2_energy);
            auto lepton2_pz_recon = lepton2_pz * (lepton2_E_recon / lepton2_energy);
            TLorentzVector e2(lepton2_px_recon, lepton2_py_recon, lepton2_pz_recon, lepton2_E_recon);

            auto Mee = e1 + e2;

            if(number_of_electrons == 1 && number_of_positrons == 1)
                h1_Mee1_before->Fill(Mee.M());
            if(number_of_electrons == 2 )
                h1_Mee2_before->Fill(Mee.M());

            auto Mx=(beam+target)-(e1+e2+pro);
            if(abs(Mx.M2())>0.4)
                continue;
            auto pt= Mx.Perp();
            auto Q2= 2*Beam_E*(Mx.P()-Mx.Pz());
            if(pt/Mx.P()<0.05){
                if(number_of_electrons == 1 && number_of_positrons == 1)
                    h1_Mee1_case1->Fill(Mee.M());
                if(number_of_electrons == 2 )
                    h1_Mee2_case1->Fill(Mee.M());    
            }

            if(pt/Mx.P()<0.1){
                if(number_of_electrons == 1 && number_of_positrons == 1)
                    h1_Mee1_case2->Fill(Mee.M());
                if(number_of_electrons == 2 )
                    h1_Mee2_case2->Fill(Mee.M());
            }

            if(Q2<0.2){
                if(number_of_electrons == 1 && number_of_positrons == 1)
                    h1_Mee1_case3->Fill(Mee.M());
                if(number_of_electrons == 2 )
                    h1_Mee2_case3->Fill(Mee.M());
            }
            
            // SELECT THE ELECTRON-POSITRON PAIR
        } // end while
        printf("num_files = %d\n", num_files);
    } // End for "Runs"

    // gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("can", "canvas", 1000, 1000);
    TString outname = output_file + ".pdf";
    can->Divide(2,2);
    can->cd(1);
    gPad->SetLogy();
    h1_Mee2_before->SetLineColor(kRed);
    h1_Mee1_before->SetMaximum(std::max(h1_Mee1_before->GetMaximum(),h1_Mee2_before->GetMaximum()) * 1.2);
    h1_Mee1_before->Draw();
    h1_Mee2_before->Draw("same");
    can->cd(2);
    gPad->SetLogy();
    h1_Mee2_case1->SetLineColor(kRed);
    h1_Mee1_case1->SetMaximum(std::max(h1_Mee1_case1->GetMaximum(),h1_Mee2_case1->GetMaximum()) * 1.2);
    h1_Mee1_case1->Draw();
    h1_Mee2_case1->Draw("same");
    can->cd(3);
    gPad->SetLogy();
    h1_Mee2_case2->SetLineColor(kRed);
    h1_Mee1_case2->SetMaximum(std::max(h1_Mee1_case2->GetMaximum(),h1_Mee2_case2->GetMaximum()) * 1.2);
    h1_Mee1_case2->Draw();
    h1_Mee2_case2->Draw("same");
    can->cd(4);
    gPad->SetLogy();
    h1_Mee2_case3->SetLineColor(kRed);
    h1_Mee1_case3->SetMaximum(std::max(h1_Mee1_case3->GetMaximum(),h1_Mee2_case3->GetMaximum()) * 1.2);
    h1_Mee1_case3->Draw();
    h1_Mee2_case3->Draw("same");

    
    can->Print(outname);


    TFile *outputfile = new TFile("./R_bg_est/"output_file+"_Bg_est.root","RECREATE");
    outputfile->cd();
    h1_Mee1_before->Write();
    h1_Mee2_before->Write();
    h1_Mee1_case1->Write();
    h1_Mee2_case1->Write();
    h1_Mee1_case2->Write();
    h1_Mee2_case2->Write();
    h1_Mee1_case3->Write();
    h1_Mee2_case3->Write();
    
    outputfile->Close();


    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    gApplication->Terminate();

    return 0;
}
