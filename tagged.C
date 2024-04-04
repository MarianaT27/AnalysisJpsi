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



int tagged(string nameFile="PASS1_FALL", int version=-19, double Beam_E=10.2, int PASS=1, int max=0, int train=0) {
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();


    //********************************
    //DECLARATION OF VARIABLES
    //********************************

    //MISSING MOMENTUM AND INVARIANT MASS
    //Acting
    TLorentzVector miss;
    TLorentzVector invariant_ee;
    TLorentzVector invariant_eeFT;
    TLorentzVector electronFT_vec;  //FT e-
    Double_t Q2, W;

    Double_t electron_E_recon;
    Double_t positron_E_recon;

    Double_t electron_px_recon, electron_py_recon, electron_pz_recon;
    Double_t positron_px_recon, positron_py_recon, positron_pz_recon;

    Double_t px,py,pz,E;
    Double_t p;

    //Testing
    TLorentzVector transverse;
    Double_t Pt;


    //CONFIG VARIABLES
    //Double_t Beam_E=E;//
    //FALL 10.6
    //SPRING 10.2
    TLorentzVector beam(0,0,Beam_E,Beam_E);
    TLorentzVector target(0,0,0,0.938);
    Int_t run_number;
    Int_t event_number;

    //Electron variables
    Double_t electron_p, electron_theta, electron_phi, electron_px, electron_py, electron_pz;
    Double_t electron_vx, electron_vy, electron_vz, electron_vt;
    Double_t electron_sf;
    Double_t electron_pcal_v, electron_pcal_w;
    Double_t electron_pcal_energy, electron_ecin_energy, electron_ecout_energy; //Electron Energies
    Double_t electron_SCI_time,electron_SCI_path;
    Double_t electron_energy;

    //Positron variables
    Double_t positron_p, positron_theta, positron_phi, positron_px, positron_py, positron_pz;
    Double_t positron_vx, positron_vy, positron_vz, positron_vt;
    Double_t positron_sf;
    Double_t positron_pcal_v, positron_pcal_w;
    Double_t positron_pcal_energy, positron_ecin_energy, positron_ecout_energy; //Positron Energies
    Double_t positron_energy;
    Double_t positron_SCI_time,positron_SCI_path, positron_chi2pid;

    //Photon variables
    Double_t photon_p, photon_theta, photon_phi, photon_px, photon_py, photon_pz;
    Double_t photon_vx, photon_vy, photon_vz, photon_vt,photon_energy;

    //Forward Tagger Electron variables
    Double_t FT_E_old,FT_E_new;
    Double_t FT_Px_old,FT_Px_new;
    Double_t FT_Py_old,FT_Py_new;
    Double_t FT_Pz_old,FT_Pz_new;
    Double_t FT_vt;
    Double_t FT_x,FT_y,FT_z, R;
    Double_t FT_vtC;
    //Testing
    Double_t FT_vtFcal,FT_xFcal,FT_yFcal,FT_zFcal,RFcal;

    int N_FT;
    int number_of_electrons;
    int number_of_positrons;
    int number_protons;
    int JPSI=0;
                    //  0    1    2   3     4    5    6    7    8    9
    int RunList[108]={5036,5038,5043,5046,5051,5053,5117,5124,5125,5126,//0
                      5127,5128,5130,5139,5153,5159,5160,5163,5164,5165,//10
                      5166,5167,5168,5169,5180,5181,5182,5183,5189,5190,//20
                      5191,5193,5196,5197,5198,5199,5200,5202,5203,5204,//30
                      5205,5206,5208,5211,5212,5215,5216,5219,5220,5221,//40
                      5222,5229,5230,5231,5232,5233,5234,5235,5238,5239,//50
                      5247,5248,5249,5252,5253,5257,5258,5259,5261,5262,//60
                      5303,5304,5305,5306,5317,5318,5319,5320,5323,5324,//70
                      5333,5339,5340,5341,5344,5345,5346,5347,5349,5354,//80
                      5355,5356,5357,5358,5359,5361,5362,5366,5367,5368,//90
                      5369,5372,5373,5375,5376,5379,5380,5381};          //100

                         //  1    2     3     4     5    6      7     8    9     10    11    12    13    14   15
    int Inbending18[146] ={5036, 5038, 5039, 5040, 5041, 5043, 5045, 5046, 5047, 5051, 5052, 5053, 5116, 5117,//3
                          5119, 5120, 5124, 5125, 5126, 5127, 5128, 5129, 5130, 5137, 5139, 5153, 5158, 5159, 5162,//4
                          5163, 5165, 5166, 5167, 5168, 5169, 5180, 5181, 5182, 5183, 5190, 5191, 5193, 5194, 5195,//5
                          5196, 5197, 5198, 5199, 5200, 5201, 5202, 5204, 5205, 5206, 5208, 5211, 5212, 5215, 5216,//6
                          5219, 5220, 5221, 5222, 5223, 5230, 5231, 5232, 5233, 5234, 5235, 5237, 5238, 5247, 5248,//7
                          5249, 5250, 5252, 5253, 5257, 5258, 5259, 5261, 5262, 5303, 5304, 5305, 5306, 5307, 5310,//8
                          5311, 5315, 5317, 5318, 5319, 5320, 5324, 5325, 5333, 5334, 5339, 5341, 5342, 5343, 5344,//9
                          5345, 5346, 5347, 5349, 5351, 5354, 5355, 5356, 5357, 5358, 5359, 5360, 5361, 5362, 5366,//10
                          5367, 5368, 5369, 5372, 5373, 5374, 5375, 5376, 5377, 5378, 5379, 5380, 5381, 5383, 5386,//11
                          5390, 5391, 5392, 5393, 5394, 5398, 5400, 5401, 5403, 5404, 5406, 5407};

                          int Outbending18[157]={5666, 5665, 5664, 5663, 5662, 5656, 5655, 5654, 5652, 5651, 5650, 5649, 5648, 5647, 5646,
                          5645, 5644, 5643, 5641, 5639, 5638, 5637, 5635, 5633, 5632, 5631, 5630, 5629, 5628, 5626,
                          5625, 5624, 5623, 5621, 5619, 5618, 5616, 5615, 5614, 5613, 5612, 5611, 5607, 5606, 5603,
                          5602, 5601, 5598, 5597, 5594, 5592, 5591, 5578, 5577, 5574, 5573, 5572, 5571, 5570, 5569,
                          5567, 5562, 5561, 5559, 5558, 5557, 5556, 5555, 5554, 5552, 5551, 5550, 5549, 5548, 5547,
                          5546, 5545, 5544, 5543, 5541, 5540, 5538, 5537, 5536, 5535, 5534, 5533, 5532, 5530, 5528,
                          5527, 5526, 5525, 5524, 5523, 5522, 5521, 5520, 5519, 5518, 5517, 5516, 5507, 5505, 5500,
                          5499, 5498, 5497, 5487, 5486, 5485, 5483, 5482, 5481, 5480, 5479, 5478, 5476, 5475, 5474,
                          5473, 5472, 5471, 5470, 5469, 5468, 5467, 5466, 5465, 5464, 5460, 5456, 5455, 5454, 5453,
                          5452, 5451, 5450, 5449, 5448, 5447, 5445, 5441, 5440, 5438, 5437, 5436, 5435, 5434, 5432,
                          5430, 5429, 5426, 5425, 5424, 5423};

                          int Spring19[114] = {6619, 6620, 6631, 6632, 6636, 6637, 6638, 6639, 6640, 6642, 6645, 6647, 6648, 6650, 6651,
                          6652, 6654, 6655, 6656, 6657, 6658, 6660, 6661, 6662, 6663, 6664, 6665, 6666, 6667, 6668,
                          6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678, 6680, 6682, 6683, 6684, 6685, 6687, 6688,
                          6689, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6704, 6705, 6706, 6707, 6708,
                          6709, 6710, 6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6728, 6729, 6730, 6731,
                          6732, 6733, 6734, 6736, 6737, 6738, 6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747, 6748,
                          6749, 6750, 6753, 6754, 6755, 6756, 6757, 6759, 6760, 6762, 6763, 6764, 6765, 6767, 6768,
                          6769, 6775, 6776, 6777, 6778, 6779, 6780, 6781, 6783};

    int S19_P2[120]  =   {6616, 6618, 6619, 6620, 6631, 6632, 6636, 6637, 6638, 6639, 6640, 6642, 6645, 6647, 6648, //15
                          6650, 6651, 6652, 6654, 6655, 6656, 6657, 6658, 6660, 6661, 6662, 6663, 6664, 6665, 6666, //30
                          6667, 6668, 6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678, 6680, 6682, 6683, 6684, 6685, //45
                          6687, 6688, 6689, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6704, 6705, 6706, //60
                          6707, 6708, 6709, 6710, 6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6722, 6723, //75
                          6724, 6725, 6728, 6729, 6730, 6731, 6732, 6733, 6734, 6736, 6737, 6738, 6739, 6740, 6741, //90
                          6742, 6743, 6744, 6746, 6747, 6748, 6749, 6750,       6753, 6754, 6755, 6756, 6757, 6759, //104
                          6760, 6762, 6763, 6764, 6765, 6767, 6768, 6769, 6775, 6776, 6777, 6778, 6779, 6780, 6781, //119
                          6783};                                                                                    //120

                          //6748 6751 6657 missing
                          


                          //GOLD: 6645,6648,6667,6669,6673,6683,6718,6719,6731,6762,6764
                          //SILVER: 6650,6658,6682

                       /* 6616, 6618, 6619, 6620, 6631, 6632, 6636, 6637, 6638, 6639, 6640, 6642, 6645, 6647, 6648, 
                          6650, 6651, 6652, 6654, 6655, 6656, 6657, 6658, 6660, 6661, 6662, 6663, 6664, 6665, 6666, 
                          6667, 6668, 6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678, 6680, 6682, 6683, 6684, 6685, 
                          6687, 6688, 6689, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6704, 6705, 6706, 
                          6707, 6708, 6709, 6710, 6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6722, 6723, 
                          6724, 6725, 6728, 6729, 6730, 6731, 6732, 6733, 6734, 6736, 6737, 6738, 6739, 6740, 6741, 
                          6742, 6743, 6744, 6746, 6747, 6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756, 6757, 6759, 
                          6760, 6762, 6763, 6764, 6765, 6767, 6768, 6769, 6775, 6776, 6777, 6778, 6779, 6780, 6781, 
                          6783*/


                          int Gold[11]  =   {6645,6648,6667,6669,6673,6683,6718,6719,6731,6762,6764};
                          int Silver[3]  = {6650,6658,6682};
   // int S19_P2[11]  =   {6712, 6716, 6769};



    //EVENT SELECTION
                          TH2F* h_electron_ECin_vs_PCAL = new TH2F("h_electron_ECin_vs_PCAL","Electron SF-ECin vs. SF-PCAL",200,0.0,0.3,100,0.0,0.2);
                          TH2F* h_positron_EC_vs_PCAL = new TH2F("h_positron_EC_vs_PCAL","Positrons SF-EC vs. SF-PCAL",200,0.0,0.3,100,0.0,0.2);
                          TH1F* h_electron_zvertex = new TH1F("h_electron_zvertex","z-vertex distribution electron",200,-15,35);
                          TH1F* h_vertex_timediff = new TH1F("h_vertex_timediff","vertex time difference e- e+ ",120,-2,2);


                          TH1F* h_p_electron_before = new TH1F("h_p_electron_before","M(e+e-) Before Radiative Loss Correction; [GeV]; ",120,0,4);
                          TH1F* h_p_electron_after = new TH1F("h_p_electron_after","M(e+e-) After Radiative Loss Correction; [GeV]; ",120,2,4);


    //ENERGY CORRECTION
                          TH2F* h_delta_theta_vs_phi_positron = new TH2F("h_delta_theta_vs_phi_positron","delta_theta vs delta_phi positron",200,-10,10,150,-30,30);
                          TH2F* h_beta_momentum = new TH2F("h_beta_momentum","beta vs momentum",27,6616,6755,7,0,7);
                          TH2F* h_delta_theta_vs_phi_electron = new TH2F("h_delta_theta_vs_phi_electron","delta_theta vs delta_phi electron",200,-10,10,150,-30,30);

                          TH1I* h_electron_photons = new TH1I("h_electron_photons","Photons radiated from electrons",7,0,7);
                          TH1I* h_positron_photons = new TH1I("h_positron_photons","Photons radiated from positrons",7,0,7);


    //TAGGED
                          TH1F* h_NFTelectron = new TH1F("h_NFTelectron","# FT electrons",10,0,10);
                          TH1F* h_vertex_timediff_FT_FD = new TH1F("h_vertex_timediff_FT_FD","Vertex time difference between FT_e^- and FT_e^-; #Deltat_v=FD_e^- - FT_e^-, ns;",350,-25,25);



    //KINEMATICS
                          TH1F* h_E_photon = new TH1F("h_E_photon","Photon Energy",120,0,0);
                          TH1F* h_Q2 = new TH1F("h_Q2","Q^2;Q^2;Counts",120,0,0);
                          TH1F* h_W = new TH1F("h_W","Hadronic Mass;W;Counts",120,4,4.6);
                          TH1F* h_MM= new TH1F("h_MM","Missing Mass; MM; Counts ",90,0,2);
                          TH2F* h_MM_vs_Deltat= new TH2F("h_MM_vs_Deltat","MM vs #Delta t; #Delta t_{e^{-}-e'^{-}};Missing Mass",150,-2,2,150,0,2);
                          TH2F* h_IM_vs_Deltat= new TH2F("h_IM_vs_Deltat","M(e+e-) vs #Delta t; #Delta t_{e^{-}-e'^{-}};Invariant Mass",150,-2,2,150,2,4);


                          TH2F* h_mm_vs_invariantmass_t= new TH2F("h_mm_vs_invariantmass_t","MM vs M(e^+e^-) for t+-2;M(e^+e^-), GeV;Missing Mass, GeV",150,0,4,150,0,3);
                          TH2F* h_mm_vs_invariantmass_no_t= new TH2F("h_mm_vs_invariantmass_no_t","MM vs M(e^+e^-)",150,0,4,150,0,3);

    //RESULTS RESULTS
                          TH1F* h_Invariant= new TH1F("h_im","M(e+e-) ",90,2.0,4);
                          TH1F* h_Invariant_not= new TH1F("h_Invariant_not","M(e+e-) outside time range",90,2.0,4);
                          TH1F* h_Invariant_noMcut= new TH1F("h_Invariant_noPcut","M(e+e-) no MM cut",90,2.0,4);


    //ANALYSIS
                          FILE *f_results;
                          FILE *f_events;

    //string nameFile="Pass1_v25";
                          f_results = fopen(("/lustre19/expphy/volatile/clas12/mtenorio/"+nameFile+"_dat.txt").c_str(), "w");
    //f_events = fopen(("/lustre19/expphy/volatile/clas12/mtenorio/"+nameFile+"_bgk.csv").c_str(), "w");

                          if(max==0){
                            if(version==-18)
                                max=146;
                            if(version==+18)
                                max=157;
                            if(version==-19)
                                max=114;
                        }

    //Start
    for(int fc =1; fc<2; fc++) {//Run 5032 to 5419 // 6616 6783

        char filename1[500];

        if(PASS==1){
            if(version==-18)
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",Inbending18[fc]);  //Inbending FALL
                //sprintf(filename1,"/lustre19/expphy/volatile/clas12/mtenorio/MC/epee/out.hipo");  //Inbending FALL
            else if(version==+18)
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",Outbending18[fc]);  //Outbending FALL
            else if(version==-19)
                //sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",fc);  //Inbending SPRING
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",Silver[fc]);  //Inbending SPRING

        }
        else if(PASS==2){
            //if(train==1)//train=1 is old v30
                //sprintf(filename1,"/volatile/clas12/rg-a/production/pass0/Spring19/sp19_dst_v1_32/dst/train/jpsitcs/jpsitcs_00%d.hipo",S19_P2[fc]);
            //if(train==2)
            //sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/jpsitcs/jpsitcs_00%d.hipo",Silver[fc]);
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",Silver[fc]);
            //if(train==0)
                //sprintf(filename1,"/volatile/clas12/mtenorio/DC_smearing/MCPhysics/JPsi/recon/ep_epjp_00%d.hipo ",fc);
            
        }


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
        hipo::bank  dataCHE;
        hipo::bank  dataFT;
        hipo::bank  dataSCI;
        hipo::bank  dataEVENT;



        hipo::bank EVENT(factory.getSchema("REC::Event"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank FTPART(factory.getSchema("RECFT::Particle"));
        hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
        hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
        hipo::bank FT(factory.getSchema("REC::ForwardTagger"));
        hipo::bank HEADER(factory.getSchema("RUN::config"));
        hipo::bank SCI(factory.getSchema("REC::Scintillator"));




        int counter = 0;


        while(reader.next()==true ){//Loops all events




            counter++;

            particl electron;
            particl positron;
            particl electronFT;
            particl photon;
            particl muon;
            particl muonp;
            particl proton;

            particl electron_photon1;
            particl electron_photon2;
            particl positron_photon1;
            particl positron_photon2;

            reader.read(event);

            event.getStructure(EVENT);
            event.getStructure(PART);
            event.getStructure(FTPART);
            event.getStructure(CALO);
            event.getStructure(CHE);
            event.getStructure(FT);
            event.getStructure(HEADER);
            event.getStructure(SCI);





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

                //PID 11 Means Electron, Status Between 2000 and 4000 Means Forward Detector
            if(pid==11 && abs(status)>=2000 && abs(status)<4000) {
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
            if(pid==-11 && abs(status)>=2000 && abs(status)<4000 ) {
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
            if(pid==2212 && abs(status)>2000 && abs(status)<4000) {
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

            //REC::Scintillator
            for (int i1 = 0; i1 < SCI.getRows(); i1++) {
                int pindex = SCI.getShort("pindex",i1);
                response ftofresponse;
                ftofresponse.energy = SCI.getFloat("energy",i1);
                ftofresponse.time = SCI.getFloat("time",i1);
                ftofresponse.index = SCI.getShort("index",i1);
                ftofresponse.layer = SCI.getByte("layer",i1);
                ftofresponse.sector = SCI.getByte("sector",i1);
                ftofresponse.path = SCI.getFloat("path",i1);


                if(pindex==electron.index && ftofresponse.layer==1) {//Layer 1 is FTOF1A
                  electron.responses.insert(pair<int, response>(10,ftofresponse));
              }
                if(pindex==electron.index && ftofresponse.layer==2) {//Layer 2 is FTOF1B
                  electron.responses.insert(pair<int, response>(11,ftofresponse));
              }
              if(pindex==positron.index && ftofresponse.layer==1) {
                  positron.responses.insert(pair<int, response>(10,ftofresponse));
              }
              if(pindex==positron.index && ftofresponse.layer==2) {
                  positron.responses.insert(pair<int, response>(11,ftofresponse));
              }

          }

          for(int i = 0; i < CHE.getRows(); i++){
            response htccresponse;
            int pindex = CHE.getShort("pindex",i);
            double nphe = CHE.getFloat("nphe",i);
            double time = CHE.getFloat("time",i);
            int detectorID = CHE.getByte("detector",i);
            htccresponse.energy = nphe;
            htccresponse.time = time;

            if(pindex==electron.index && detectorID==15) {
              electron.responses.insert(pair<int, response>(6,htccresponse));
          }
          if(pindex==positron.index && detectorID==15) {
            positron.responses.insert(pair<int, response>(6,htccresponse));
        }
        if(pindex==photon.index && detectorID==15) {
            photon.responses.insert(pair<int, response>(6,htccresponse));
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

            //if(event_number==108736645 && partNumber==100){
          if(partNumber==100) {

            //if(number_of_electrons==1 &&number_of_positrons==1) {
            if(number_of_electrons==1 && number_of_positrons==1) {



                electronAccepted=false;
                positronAccepted=false;


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
                electron_sf = (electron.responses[109].energy + electron.responses[110].energy + electron.responses[111].energy)/electron.lorentz.P();
                electron_pcal_v=electron.responses[109].v;
                electron_pcal_w=electron.responses[109].w;
                electron_pcal_energy = electron.responses[109].energy;
                electron_ecin_energy = electron.responses[110].energy;
                electron_ecout_energy = electron.responses[111].energy;
                electron_SCI_time=electron.responses[10].time;
                electron_SCI_path=electron.responses[10].path;
                electron_energy = electron.lorentz.E();

                positron_pcal_v=positron.responses[109].v;
                positron_pcal_w=positron.responses[109].w;


                positron_pcal_energy = positron.responses[109].energy;
                positron_ecin_energy = positron.responses[110].energy;
                positron_ecout_energy = positron.responses[111].energy;
                positron_vz = positron.vertexinfo.z;
                positron_vt=positron.vtime;
                positron_p = positron.lorentz.P();
                positron_px =positron.lorentz.Px();
                positron_py =positron.lorentz.Py();
                positron_pz =positron.lorentz.Pz();
                positron_theta = positron.lorentz.Theta()*57.2958;
                positron_phi = positron.lorentz.Phi()*57.2958;
                positron_sf = (positron.responses[109].energy + positron.responses[110].energy + positron.responses[111].energy)/positron.lorentz.P();
                positron_SCI_time=positron.responses[11].time;//Try 11 11
                positron_SCI_path=positron.responses[11].path;//Try 11
                positron_energy = positron.lorentz.E();
                positron_chi2pid= positron.chi2pid;



                    

                    if(electron_pcal_energy>0.07  && electron_p>1.7 && electron_p<beam.E()&& electron_pcal_v>9 && electron_pcal_w>9){//
                        if(version==-18||version==-19){
                            h_electron_zvertex->Fill(electron_vz);
                            if(electron_vz>-8 && electron_vz<4){
                                electronAccepted=true;
                            }
                        }
                        else
                            electronAccepted=true;
                    }
                    //CUTS FOR POSITRON
                    if(positron_pcal_energy>0.07  && positron_p>1.7 && abs(positron_chi2pid)<5 && positron_pcal_v>9&& positron_pcal_w>9){//
                        if(((positron_ecin_energy+positron_ecout_energy)/positron_p)>=(0.195-positron_pcal_energy/positron_p)){
                            if(version==+18){
                                h_electron_zvertex->Fill(positron_vz);
                                if(positron_vz>-8 && positron_vz<4){
                                    positronAccepted=true;
                                }
                            }
                            else
                                positronAccepted=true;
                        }
                    }

                    
                    int electron_photon=0;
                    int positron_photon=0;


                    //CUTS FOR TIME
                    Double_t T_vertex=abs(positron_SCI_time)-positron_SCI_path/(29.9792458);


                    if(electronAccepted){
                        h_electron_ECin_vs_PCAL->Fill(electron_pcal_energy/electron_p,electron_ecin_energy/electron_p);

                    }

                    if(positronAccepted){
                        h_positron_EC_vs_PCAL->Fill(positron_pcal_energy/positron_p,(positron_ecin_energy+positron_ecout_energy)/positron_p);
                    }



                    if(electronAccepted && positronAccepted){
                        h_vertex_timediff->Fill(electron_vt-positron_vt);

                        //CUT TO TIME AND VERTEX VERTEX
                        if(abs(electron_vt-positron_vt)<=1){


                            Double_t NocorrectedIM=0;

                            Double_t electron_photonE=0;
                            Double_t positron_photonE=0;
                            Double_t E_photon=0;

                            NocorrectedIM=(electron.lorentz+positron.lorentz).M();




                            for(int i = 0; i < PART.getRows(); i++){
                                float  px = PART.getFloat("px",i);
                                float  py = PART.getFloat("py",i);
                                float  pz = PART.getFloat("pz",i);
                                float beta = PART.getFloat("beta",i);
                                if(PART.getInt("pid",i)==22  && beta>=0.91 && beta<=1.09) {    //){// 0.94 to 1.06
                                    TLorentzVector photon_vector;
                                    photon_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0*0.0));
                                    //-----PLOT ENERGY CORRECTION -----
                                    h_delta_theta_vs_phi_electron->Fill(photon_vector.Theta()*57.2958-electron_theta,photon_vector.Phi()*57.2958-electron_phi);
                                    h_delta_theta_vs_phi_positron->Fill(photon_vector.Theta()*57.2958-positron_theta,photon_vector.Phi()*57.2958-positron_phi);
                                        //h_p_electron_before->Fill(photon_vector.E());
                                    //-----------


                                    if(abs(photon_vector.Theta()*57.2958-electron_theta)<0.7){//&&electron_photon<=2){//Was 0.7
                                        h_p_electron_before->Fill(photon_vector.E());
                                        electron_photonE=electron_photonE+photon_vector.E();
                                        electron_photon++;
                                    }
                                    else if((abs(photon_vector.Theta()*57.2958-positron_theta)<0.7)){//&&positron_photon<=2){
                                        h_p_electron_before->Fill(photon_vector.E());
                                        positron_photonE=positron_photonE+photon_vector.E();
                                        positron_photon++;
                                    }
                                }
                            }
                            //-------PLOT NUMBER PHOTONS AND # FT ELECTRONS
                            h_electron_photons->Fill(electron_photon);
                            h_positron_photons->Fill(positron_photon);
                            h_beta_momentum->Fill(rn,electron_photon);
                            //------------


                            //ENERGY CORRECTION
                            electron_E_recon=electron_energy+electron_photonE;
                            positron_E_recon=positron_energy+positron_photonE;

                            electron_px_recon=electron_px*(electron_E_recon/electron_energy);
                            electron_py_recon=electron_py*(electron_E_recon/electron_energy);
                            electron_pz_recon=electron_pz*(electron_E_recon/electron_energy);


                            positron_px_recon=positron_px*(positron_E_recon/positron_energy);
                            positron_py_recon=positron_py*(positron_E_recon/positron_energy);
                            positron_pz_recon=positron_pz*(positron_E_recon/positron_energy);




                            //MISSING MASS
                            TLorentzVector ele(electron_px_recon,electron_py_recon,electron_pz_recon,electron_E_recon);
                            TLorentzVector pos(positron_px_recon,positron_py_recon,positron_pz_recon,positron_E_recon);


                            h_NFTelectron->Fill(N_FT);

                            if(N_FT==1){

                                R=sqrt((FT_xFcal*FT_xFcal)+(FT_yFcal*FT_yFcal)+(FT_zFcal*FT_zFcal));
                                FT_vt=FT_vtFcal-(R/29.9792458);


                                //Do correction of energy
                                FT_E_old=electronFT.lorentz.E();


                                //NOT USE IN PASS2
                                if(PASS==1){
                                    FT_E_new=-0.03689 + (1.1412*FT_E_old) - (0.04316*pow(FT_E_old,2))+ (0.007046*pow(FT_E_old,3))- (0.0004055*pow(FT_E_old,4));
                                }
                                else{
                                    FT_E_new=electronFT.lorentz.E();
                                }

                                FT_Px_new=electronFT.lorentz.Px()*(FT_E_new/FT_E_old);
                                FT_Py_new=electronFT.lorentz.Py()*(FT_E_new/FT_E_old);
                                FT_Pz_new=electronFT.lorentz.Pz()*(FT_E_new/FT_E_old);

                                //set variables
                                electronFT_vec.SetPxPyPzE(FT_Px_new,FT_Py_new,FT_Pz_new,FT_E_new);
                                Double_t eFT_theta=electronFT_vec.Theta()*57.2958;
                                Double_t eFT_phi=electronFT_vec.Phi()*57.2958;

                                Double_t pos_theta=pos.Theta()*57.2958;
                                Double_t pos_phi=pos.Phi()*57.2958;
                                Double_t ele_theta=ele.Theta()*57.2958;
                                Double_t ele_phi=ele.Phi()*57.2958;

                                //MISSING MASS & INVARIANT MASS
                                miss=(beam+target)-(pos+ele+electronFT_vec);
                                invariant_ee=(pos+ele);

                                //MISSING MASS & INVARIANT MASS With NO RADIATIVE CORRECTION
                                //miss=(beam+target)-(positron.lorentz+electron.lorentz+electronFT_vec);
                                //invariant_ee=(positron.lorentz+electron.lorentz);

                                Q2=2*Beam_E*FT_E_new*(1-cos(electronFT_vec.Theta()));
                                E_photon=Beam_E-FT_E_new;
                                W=sqrt((0.938*0.938)+(2*0.938*E_photon)-Q2);

                                h_E_photon->Fill(E_photon);
                                h_Q2->Fill(Q2);

                                h_vertex_timediff_FT_FD->Fill(electron_vt-FT_vt);

                                h_MM_vs_Deltat->Fill(electron_vt-FT_vt,miss.M());
                                h_IM_vs_Deltat->Fill(electron_vt-FT_vt,invariant_ee.M());


                                fprintf(f_results,"%d %d %7.6f %7.6f %7.6f \n",rn,en,electron.vtime-FT_vt,  miss.M(),invariant_ee.M());
                                

                                //MM vs IM
                                if((electron.vtime-FT_vt)<=1.5 && (electron.vtime-FT_vt)>=-2.5){
                                    h_mm_vs_invariantmass_t->Fill(invariant_ee.M(),miss.M());
                                    h_MM->Fill(miss.M());


                                    if(0.65 <=miss.M() &&  miss.M()<= 1.3){//1 sigmas 0.65 <=miss.M() &&  miss.M()<= 1.3
                                        h_Invariant->Fill(invariant_ee.M());
                                        //if(invariant_ee.M()>=3.0 && invariant_ee.M()<=3.2){
                                            //fprintf(f_results,"%d, %d, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f,%7.6f, %7.6f, %7.6f, 1 \n",rn,en, invariant_ee.M(),
                                             //   ele.P(), ele_theta, ele_phi, pos.P(),pos_theta, pos_phi,electronFT_vec.P(),eFT_theta,eFT_phi );
                                            //fprintf(f_results,"%d %d %7.6f \n",rn,en, invariant_ee.M(), miss.M());
                                        if(invariant_ee.M()>=3.0 && invariant_ee.M()<=3.2){
                                            h_W->Fill(W);
                                        }
                                        
                                    }
                                    else{
                                        h_Invariant_noMcut->Fill(invariant_ee.M());
                                    }

                                }
                                else{
                                    h_mm_vs_invariantmass_no_t->Fill(invariant_ee.M(),miss.M());
                                    h_Invariant_not->Fill(invariant_ee.M());
                                    //fprintf(f_events,"%d, %d, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f,%7.6f, %7.6f, %7.6f, 0 \n",rn,en, invariant_ee.M(),
                                                //ele.P(), ele_theta, ele_phi, pos.P(),pos_theta, pos_phi,electronFT_vec.P(),eFT_theta,eFT_phi );
                                }


                            }//N_FT==1
                            


                        }//CUT TO TIME AND VERTEX

                    }//iF A ELECTRON-POSITRON PAIR IS ACCEPTED

                    
                }//SELECT THE ELECTRON-POSITRON PAIR


            }//If particle


            counter++;
        }//end while
        printf("processed events = %d\n",counter);
        printf("run = %d\n",fc);
    }//End for "Runs"

    fclose(f_results);


    TCanvas *can = new TCanvas("can","canvas",200,10,700,700);
    string pdf_original=nameFile+".pdf";

    

    can->Divide(2,2);

    can->cd(1);
    h_electron_ECin_vs_PCAL->Draw("colz");

    can->cd(2);
    h_positron_EC_vs_PCAL->Draw("colz");

    can->cd(3);
    h_electron_zvertex->Draw();

    can->cd(4);
    gPad->SetLogy();
    h_vertex_timediff->Draw();

    can->Print( (pdf_original + "(").c_str());


    can->cd(1);
    gPad->SetLogy(0);
    h_delta_theta_vs_phi_electron->Draw("colz");

    can->cd(2);
   h_delta_theta_vs_phi_positron->Draw("colz");

    can->cd(3);
    gPad->SetLogy();//log
    h_electron_photons->Draw();

    can->cd(4);
    h_positron_photons->Draw();
   can->Print( (pdf_original + "(").c_str());

    can->Clear();
    can->Divide(1,2);
    can->cd(1);
    gPad->SetLogy();
    h_NFTelectron->Draw();

    can->cd(2);
    h_vertex_timediff_FT_FD->Draw();
    can->Print( (pdf_original + "(").c_str());

    can->Clear();
    gPad->SetLogy(0);

    h_MM->Draw();
    can->Print( (pdf_original + "(").c_str());


    h_mm_vs_invariantmass_t->Draw("colz");
    can->Print( (pdf_original + "(").c_str());

    h_mm_vs_invariantmass_no_t->Draw("colz");
    can->Print( (pdf_original + "(").c_str());


    h_Invariant->Draw();
    can->Print( (pdf_original + "(").c_str());

    h_W->Draw();
    can->Print( (pdf_original + "(").c_str());

    h_E_photon->Draw();
    can->Print( (pdf_original + "(").c_str());

    h_Invariant_not->Draw();
    can->Print( (pdf_original + ")").c_str());




    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;

}




