#ifndef Utils
#define Utils

int F18in_P2[170] = {5032, 5036, 5038, 5039, 5040, 5041, 5043, 5045, 5046, 5047, 5051,                         // 11
                     5052, 5053, 5116, 5117, 5119, 5120, 5124, 5125, 5126, 5127, 5128, 5129, 5130, 5137, 5138, // 26
                     5139, 5153, 5158, 5159, 5160, 5162, 5163, 5164, 5165, 5166, 5167, 5168, 5169, 5180, 5181, // 41
                     5182, 5183, 5190, 5191, 5193, 5194, 5195, 5196, 5197, 5198, 5199, 5200, 5201, 5202,       // 55
                     5203, 5204, 5205, 5206, 5208, 5211, 5212, 5215, 5216, 5219, 5220, 5221, 5222, 5223, 5225, // 70
                     5229, 5230, 5231, 5232, 5233, 5234, 5235, 5237, 5238, 5239, 5247, 5248, 5249, 5250,       // 84
                     5252, 5253, 5257, 5258, 5259, 5261, 5262, 5300, 5301, 5302, 5303, 5304, 5305, 5306, 5307, // 99
                     5310, 5311, 5315, 5316, 5317, 5318, 5319, 5320, 5323, 5324, 5325, 5333, 5334, 5335, 5336, // 115
                     5339, 5340, 5341, 5342, 5343, 5344, 5345, 5346, 5347, 5349, 5351, 5354, 5355, 5356, 5357, // 129
                     5358, 5359, 5360, 5361, 5362, 5366, 5367, 5368, 5369, 5370, 5371, 5372, 5373, 5374, 5375, // 144
                     5376, 5377, 5378, 5380, 5381, 5382, 5383, 5386, 5390, 5391, 5392, 5393, 5398,             // 157
                     5400, 5401, 5402, 5403, 5404, 5406, 5407, 5414, 5415, 5416, 5417, 5418, 5419};            // 170
// MISSING: 5026, 5027, 5030, 5031, 5189, 5236, 5379, 5399

int F18out_P2[180] = {5423, 5424, 5425, 5426, 5428, 5429, 5430, 5431, 5432, 5434, 5435, 5436, 5437, 5438,       // 14
                      5439, 5440, 5441, 5442, 5443, 5444, 5445, 5447, 5448, 5449, 5450, 5451, 5452, 5453, 5454, // 29
                      5455, 5456, 5457, 5460, 5462, 5464, 5465, 5466, 5467, 5468, 5469, 5470, 5471, 5472, 5473, // 44
                      5474, 5475, 5476, 5478, 5479, 5481, 5482, 5483, 5485, 5486, 5487, 5495, 5496, 5497, 5498, // 59
                      5499, 5500, 5504, 5505, 5507, 5516, 5517, 5518, 5519, 5520, 5521, 5522, 5523, 5524, 5525, // 74
                      5526, 5527, 5528, 5530, 5532, 5533, 5534, 5535, 5536, 5537, 5538, 5540, 5541, 5542, 5543, // 89
                      5544, 5545, 5546, 5547, 5548, 5549, 5550, 5551, 5552, 5554, 5555, 5556, 5557, 5558, 5559, // 104
                      5561, 5562, 5564, 5565, 5566, 5567, 5569, 5570, 5571, 5572, 5573, 5574, 5577, 5578,       // 118
                      5586, 5589, 5590, 5591, 5592, 5594, 5595, 5597, 5598, 5600, 5601, 5602, 5603, 5604,       // 132
                      5606, 5607, 5610, 5611, 5612, 5613, 5614, 5615, 5616, 5617, 5618, 5619, 5620, 5621,       // 146
                      5623, 5624, 5625, 5626, 5627, 5628, 5629, 5630, 5631, 5632, 5633, 5635, 5637, 5638,       // 161
                      5639, 5641, 5643, 5644, 5645, 5646, 5647, 5648, 5649, 5650, 5651, 5652, 5654, 5655, 5656, // 176
                      5662, 5663, 5664, 5665, 5666};                                                            // 181
// MISSING: 5422, 5581, 5584, 5609, ADD: 5480

int S19_P2[120] = {6616, 6618, 6619, 6620, 6631, 6632, 6636, 6637, 6638, 6639, 6640, 6642, 6645, 6647, 6648, // 15
                   6650, 6651, 6652, 6654, 6655, 6656, 6657, 6658, 6660, 6661, 6662, 6663, 6664, 6665, 6666, // 30
                   6667, 6668, 6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678, 6680, 6682, 6683, 6684, 6685, // 45
                   6687, 6688, 6689, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6704, 6705, 6706, // 60
                   6707, 6708, 6709, 6710, 6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6722, 6723, // 75
                   6724, 6725, 6728, 6729, 6730, 6731, 6732, 6733, 6734, 6736, 6737, 6738, 6739, 6740, 6741, // 90
                   6742, 6743, 6744, 6746, 6747, 6748, 6749, 6750, 6753, 6754, 6755, 6756, 6757, 6759,       // 104
                   6760, 6762, 6763, 6764, 6765, 6767, 6768, 6769, 6775, 6776, 6777, 6778, 6779, 6780, 6781, // 119
                   6783};

const string fitnames[5] = {"Gauss+Pol", "Crystal_Ball+Pol", "Gauss+Exp", "Crystal_Ball+Exp", "GaussExp+Exp"};
const string lnames[5] = {"G+Pol", "CB+Pol", "G+Exp", "CB+Exp", "GaussExp+Exp"};

TH1F *h_Invariant = new TH1F("h_im", ";M(e+e-),GeV ", 100, 2.0, 3.5);
TH2F *h_mm_vs_invariantmass = new TH2F("h_mm_vs_invariantmass", "MM vs M(e^+e^-);M(e^+e^-), GeV;Missing Mass, GeV", 150, 1, 3.5, 150, -1, 1);
TH2F *h_qq_vs_invariantmass = new TH2F("h_qq_vs_invariantmass", "Q^{2} vs M(e^+e^-);M(e^+e^-), GeV; Q^{2}", 150, 1, 3.5, 150, 0, 0.5);

TH1F *h_Mee = new TH1F("h_Mee", ";M(e+e-),GeV ", 90, 2.0, 4);
TH1F *h_MM = new TH1F("h_MM", ";M_{X},GeV ", 90, -0.5, 0.5);
TH1F *h_Q2 = new TH1F("h_Q2", ";Q^{2} ", 150, 0, 1);
TH1F *h_E_photon = new TH1F("h_E_photon", "Photon Energy", 120, 0, 0);

TH1F *h_Invariant_noCuts = new TH1F("h_im_noCuts", ";M(e+e-),GeV ", 100, 2.5, 3.5);
TH1F *h_Invariant_tCut = new TH1F("h_im_tCut", ";M(e+e-),GeV ", 100, 2.5, 3.5);
TH1F *h_MX = new TH1F("h_MX", ";M_{X},GeV ", 90, 0.5, 1.5);
TH1F *h_W = new TH1F("h_W", ";W,GeV ", 90, 4, 4.6);
TH2F* h_mm_vs_invariantmass_t= new TH2F("h_mm_vs_invariantmass_t","M_{X} vs M(e^+e^-) for t<=2.5ns;M(e^+e^-), GeV;Missing Mass, GeV",150,0,4,150,0,3);
TH2F* h_mm_vs_invariantmass_no_t= new TH2F("h_mm_vs_invariantmass_no_t","M_{X} vs M(e^+e^-) for t>2.5ns;M(e^+e^-), GeV;Missing Mass, GeV",150,0,4,150,0,3);


TH1F *h_invariant_Egamma[4];
TH1F *h_invariant_Egamma_t[4];

TH1F *h_invariant_Egamma_t_plus[6];
TH1F *h_invariant_Egamma_t_minus[6];

double GetBSA(double Njpsi_p, double Njpsi_n)
{
  return (Njpsi_p - Njpsi_n) / (Njpsi_p + Njpsi_n);
}

void Bin_Egamma(double_t Epho, double_t invariantMass)
{

  if (Epho >= 8.2 && Epho < 9.05)
    h_invariant_Egamma[0]->Fill(invariantMass);
  else if (Epho < 9.46)
    h_invariant_Egamma[1]->Fill(invariantMass);
  else if (Epho < 10)
    h_invariant_Egamma[2]->Fill(invariantMass);
  else if (Epho < 10.6)
    h_invariant_Egamma[3]->Fill(invariantMass);
}

void Bin_oneVar(double_t Var, double_t Data)
{

  if (Var >= 8.2 && Var < 9.05)
    h_invariant_Egamma[0]->Fill(Data);
  else if (Var < 9.46)
    h_invariant_Egamma[1]->Fill(Data);
  else if (Var < 10)
    h_invariant_Egamma[2]->Fill(Data);
  else if (Var < 10.6)
    h_invariant_Egamma[3]->Fill(Data);
}

void Bin_Egamma_t(double_t Epho, double_t t, int hel, double_t invariantMass)
{

  if (t < 1)
  {
    if (Epho >= 8.2 && Epho < 9.46)
      h_invariant_Egamma_t[0]->Fill(invariantMass);
    else if (Epho < 10.6)
      h_invariant_Egamma_t[1]->Fill(invariantMass);
  }
  else if (t < 2.5)
  {
    if (Epho >= 8.2 && Epho < 9.46)
      h_invariant_Egamma_t[2]->Fill(invariantMass);
    else if (Epho < 10.6)
      h_invariant_Egamma_t[3]->Fill(invariantMass);
  }

  if (hel == +1)
  {
    if (t < 1)
    {
      if (Epho >= 8.2 && Epho < 9.8)
        h_invariant_Egamma_t_plus[0]->Fill(invariantMass);
      else if (Epho < 10.6)
        h_invariant_Egamma_t_plus[1]->Fill(invariantMass);
    }
    else if (t < 2)
    {
      if (Epho >= 8.2 && Epho < 9.8)
        h_invariant_Egamma_t_plus[2]->Fill(invariantMass);
      else if (Epho < 10.6)
        h_invariant_Egamma_t_plus[3]->Fill(invariantMass);
    }
    else if (t < 3.5)
    {
      if (Epho >= 8.2 && Epho < 9.8)
        h_invariant_Egamma_t_plus[4]->Fill(invariantMass);
      else if (Epho < 10.6)
        h_invariant_Egamma_t_plus[5]->Fill(invariantMass);
    }
  }
  if (hel == -1)
  {
    if (t < 1)
    {
      if (Epho >= 8.2 && Epho < 9.8)
        h_invariant_Egamma_t_minus[0]->Fill(invariantMass);
      else if (Epho < 10.6)
        h_invariant_Egamma_t_minus[1]->Fill(invariantMass);
    }
    else if (t < 2)
    {
      if (Epho >= 8.2 && Epho < 9.8)
        h_invariant_Egamma_t_minus[2]->Fill(invariantMass);
      else if (Epho < 10.6)
        h_invariant_Egamma_t_minus[3]->Fill(invariantMass);
    }
    else if (t < 3.5)
    {
      if (Epho >= 8.2 && Epho < 9.8)
        h_invariant_Egamma_t_minus[4]->Fill(invariantMass);
      else if (Epho < 10.6)
        h_invariant_Egamma_t_minus[5]->Fill(invariantMass);
    }
  }
}

#endif
