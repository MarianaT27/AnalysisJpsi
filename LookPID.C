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

void LookPID() {
    std::string folder_path;
    folder_path = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/";
    int f_count = 0;
    int n_jpsi = 0;

    std::cout << "Folder path: " << folder_path << std::endl;

    for (const auto& entry : fs::directory_iterator(folder_path)) {
        char filename1[500];
        f_count++;

        if (entry.path().extension() == ".hipo") {
            std::string filename = entry.path().filename().string();
            if (filename == "recon.hipo" || filename == "gemc_denoised.hipo" || filename == "gemc.hipo" || filename == "00083.hipo" || filename == "00048.hipo") {
                continue; // Skip this iteration
            }
            sprintf(filename1, "%s%s", folder_path.c_str(), filename.c_str());
            std::cout << filename1 << std::endl;
        }

        std::cout << "Analysis running on " << filename1 << std::endl;

        hipo::reader reader;
        reader.open(filename1);

        hipo::dictionary factory;
        reader.readDictionary(factory);

        hipo::event event;

        hipo::bank EVENT(factory.getSchema("REC::Event"));
        hipo::bank MCPART(factory.getSchema("MC::Particle"));

        while (reader.next() == true) { // Loop through all events
            reader.read(event);

            event.getStructure(EVENT);
            event.getStructure(MCPART);

            int nParticles = MCPART.getRows();
            for (int i = 0; i < nParticles; i++) {
                // Extract particle data from the MC::Particle bank
                int pid = MCPART.getInt("pid",i);
 
                // Check if the particle is a J/ψ meson
                if (pid == 443) {
                    n_jpsi++;
                    std::cout << "Particle " << i << " is a J/ψ meson!" << std::endl;
                }
            }
        } // End while
    } // End for loop over runs

    std::cout << "Total number of J/ψ particles found: " << n_jpsi << std::endl;
}