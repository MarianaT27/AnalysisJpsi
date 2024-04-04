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



int toroot(string nameFile="PASS1_FALL", int version=-19, double Beam_E=10.2, int PASS=1, int max=0) {
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


    TString root_file = "/lustre19/expphy/volatile/clas12/mtenorio/AnalysisJpsi/Root/"+nameFile+".root";
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

    //----------------PASS 1 run list-----------------------------------------

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

    //----------------PASS 2 run list-----------------------------------------

    int F18in_P2[170] =  {                        5032, 5036, 5038, 5039, 5040, 5041, 5043, 5045, 5046, 5047, 5051, //11
                          5052, 5053, 5116, 5117, 5119, 5120, 5124, 5125, 5126, 5127, 5128, 5129, 5130, 5137, 5138, //26
                          5139, 5153, 5158, 5159, 5160, 5162, 5163, 5164, 5165, 5166, 5167, 5168, 5169, 5180, 5181, //41
                          5182, 5183,       5190, 5191, 5193, 5194, 5195, 5196, 5197, 5198, 5199, 5200, 5201, 5202, //55
                          5203, 5204, 5205, 5206, 5208, 5211, 5212, 5215, 5216, 5219, 5220, 5221, 5222, 5223, 5225, //70
                          5229, 5230, 5231, 5232, 5233, 5234, 5235,       5237, 5238, 5239, 5247, 5248, 5249, 5250, //84
                          5252, 5253, 5257, 5258, 5259, 5261, 5262, 5300, 5301, 5302, 5303, 5304, 5305, 5306, 5307, //99
                          5310, 5311, 5315, 5316, 5317, 5318, 5319, 5320, 5323, 5324, 5325, 5333, 5334, 5335, 5336, //115
                          5339, 5340, 5341, 5342, 5343, 5344, 5345, 5346, 5347, 5349, 5351, 5354, 5355, 5356, 5357, //129
                          5358, 5359, 5360, 5361, 5362, 5366, 5367, 5368, 5369, 5370, 5371, 5372, 5373, 5374, 5375, //144
                          5376, 5377, 5378,       5380, 5381, 5382, 5383, 5386, 5390, 5391, 5392, 5393, 5398,       //157
                          5400, 5401, 5402, 5403, 5404, 5406, 5407, 5414, 5415, 5416, 5417, 5418, 5419};            //170
                          //MISSING: 5026, 5027, 5030, 5031, 5189, 5236, 5379, 5399


    int F18out_P2[181]  = {      5423, 5424, 5425, 5426, 5428, 5429, 5430, 5431, 5432, 5434, 5435, 5436, 5437, 5438, //14
                          5439, 5440, 5441, 5442, 5443, 5444, 5445, 5447, 5448, 5449, 5450, 5451, 5452, 5453, 5454, //29
                          5455, 5456, 5457, 5460, 5462, 5464, 5465, 5466, 5467, 5468, 5469, 5470, 5471, 5472, 5473, //44
                          5474, 5475, 5476, 5478, 5479, 5481, 5482, 5483, 5485, 5486, 5487, 5495, 5496, 5497, 5498, //59
                          5499, 5500, 5504, 5505, 5507, 5516, 5517, 5518, 5519, 5520, 5521, 5522, 5523, 5524, 5525, //74
                          5526, 5527, 5528, 5530, 5532, 5533, 5534, 5535, 5536, 5537, 5538, 5540, 5541, 5542, 5543, //89
                          5544, 5545, 5546, 5547, 5548, 5549, 5550, 5551, 5552, 5554, 5555, 5556, 5557, 5558, 5559, //104
                          5561, 5562, 5564, 5565, 5566, 5567, 5569, 5570, 5571, 5572, 5573, 5574, 5577, 5578,       //118
                                5586, 5589, 5590, 5591, 5592, 5594, 5595, 5597, 5598, 5600, 5601, 5602, 5603, 5604, //132
                          5606, 5607,       5610, 5611, 5612, 5613, 5614, 5615, 5616, 5617, 5618, 5619, 5620, 5621, //146
                          5623, 5624, 5625, 5626, 5627, 5628, 5629, 5630, 5631, 5632, 5633, 5634, 5635, 5637, 5638, //161
                          5639, 5641, 5643, 5644, 5645, 5646, 5647, 5648, 5649, 5650, 5651, 5652, 5654, 5655, 5656, //176
                          5662, 5663, 5664, 5665, 5666};                                                            //181
                          //MISSING: 5422, 5581, 5584, 5609

     

    int S19_P2[120]  =   {6616, 6618, 6619, 6620, 6631, 6632, 6636, 6637, 6638, 6639, 6640, 6642, 6645, 6647, 6648, //15
                          6650, 6651, 6652, 6654, 6655, 6656, 6657, 6658, 6660, 6661, 6662, 6663, 6664, 6665, 6666, //30
                          6667, 6668, 6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678, 6680, 6682, 6683, 6684, 6685, //45
                          6687, 6688, 6689, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6704, 6705, 6706, //60
                          6707, 6708, 6709, 6710, 6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6722, 6723, //75
                          6724, 6725, 6728, 6729, 6730, 6731, 6732, 6733, 6734, 6736, 6737, 6738, 6739, 6740, 6741, //90
                          6742, 6743, 6744, 6746, 6747, 6748, 6749, 6750,       6753, 6754, 6755, 6756, 6757, 6759, //104
                          6760, 6762, 6763, 6764, 6765, 6767, 6768, 6769, 6775, 6776, 6777, 6778, 6779, 6780, 6781, //119
                          6783};                                                                                    //120
                          //MISSING: 6751                                                          
                
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


    if(max==0){
        if(version==-18)
            max=170;
        if(version==+18)
            max=181;
        if(version==-19)
            max=120;
    }

    //Start
    for(int fc =0; fc<=max; fc++) {//Run 5032 to 5419 // 6616 6783

        char filename1[500];

        if(PASS==1){
            if(version==-18)
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",F18in_P2[fc]);  //Inbending FALL
            else if(version==+18)
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",F18out_P2[fc]);  //Outbending FALL
            else if(version==-19)
                sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/jpsitcs/jpsitcs_00%d.hipo",S19_P2[fc]);  //Inbending SPRING

        }
        else if(PASS==2){
          if(version==-18)
              sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/jpsitcs/jpsitcs_00%d.hipo",F18in_P2[fc]);
          else if(version==+18)
              sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/jpsitcs/jpsitcs_00%d.hipo",F18out_P2[fc]); 
          else if(version==-19)
            sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/jpsitcs/jpsitcs_00%d.hipo",S19_P2[fc]);
        }
        else{
            //sprintf(filename1,"/volatile/clas12/mtenorio/DC_smearing/MCPhysics/Jpsi/recon/ep_epjp_00%d.hipo",fc);

            sprintf(filename1,"/volatile/clas12/osg/marianat/job_7203/output/Jpsi-ep_epjp_00%d-7203-%d.hipo",fc+1,fc);
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




