#ifndef h_analysis
#define h_analysis
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


string Name_file(string nameFile, int topology)
{
    string tops, miss_s;
    switch (topology) {
        case 1:
            tops="_ee+e-";
            break;
        case 2:
            tops="_ee+e-p";
            break;
        case 3:
            tops="_epe+";
            break;
        case 4:
            tops="_epe-";
            break;
        default:
            tops="_TEST";
            break;
    }

    return nameFile + tops;
}

bool Accept_electron(Double_t electron_pcal_energy, Double_t electron_p, Double_t electron_pcal_v, Double_t electron_pcal_w, Double_t electron_vz, int version)
{
    Double_t Beam_E;
    if (abs(version) == 18)
        Beam_E = 10.6;
    else if(abs(version)==19)
        Beam_E = 10.2;
    else if(abs(version)==10)
        Beam_E = 10.0;
    else   
        Beam_E = 6.0; 

    if (electron_pcal_energy > 0.07 && electron_p > 1.95 && electron_p < Beam_E && electron_pcal_v > 9 && electron_pcal_w > 9)
    {
        if (version<0)
        {
            if (electron_vz > -8 && electron_vz < 4)
                return true;
        }
        else
            return true;
    }
    return false;
}

bool Accept_positron(Double_t positron_ecin_energy, Double_t positron_ecout_energy, Double_t positron_pcal_energy, Double_t positron_p, Double_t positron_chi2pid, Double_t positron_pcal_v, Double_t positron_pcal_w, Double_t positron_vz, int version)
{
    if (positron_pcal_energy > 0.07 && positron_p > 1.95 && abs(positron_chi2pid) < 5 && positron_pcal_v > 9 && positron_pcal_w > 9)
    { //
        if (((positron_ecin_energy + positron_ecout_energy) / positron_p) >= (0.195 - positron_pcal_energy / positron_p))
        {
            if (version > 0)
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


double proton_pcorr(double theta, double e) {
    double ag[5][15] = {
        {5.00, 7.00, 9.00, 11.00, 13.00, 15.00, 17.00, 19.00, 21.00, 23.00, 25.00, 27.00, 29.00, 31.00, 33.00},
        {0.1039E-02, -0.6952E-02, -0.9509E-02, -0.9879E-02, -0.1279E-01, -0.1157E-01, -0.1018E-01, -0.9222E-02, -0.1355E-01, -0.1207E-01, -0.9474E-02, -0.2216E-01, -0.2105E-01, -0.2118E-01, -0.2360E-01},
        {-0.6922E-03, 0.9763E-03, 0.1482E-02, 0.1530E-02, 0.2187E-02, 0.1953E-02, 0.1688E-02, 0.1668E-02, 0.2849E-02, 0.2495E-02, 0.1508E-02, 0.4215E-02, 0.3911E-02, 0.3948E-02, 0.4634E-02},
        {0.9806E-03, 0.1157E-01, 0.1485E-01, 0.1588E-01, 0.1945E-01, 0.1736E-01, 0.1551E-01, 0.1383E-01, 0.1926E-01, 0.1720E-01, 0.1464E-01, 0.3250E-01, 0.3231E-01, 0.3296E-01, 0.3608E-01},
        {-0.8024E-02, -0.1035E-01, -0.1240E-01, -0.1361E-01, -0.1518E-01, -0.1432E-01, -0.1341E-01, -0.1255E-01, -0.1462E-01, -0.1388E-01, -0.1574E-01, -0.2646E-01, -0.2820E-01, -0.3000E-01, -0.3259E-01}
    };
    double ec;

    if (theta >= ag[0][14]) {
        ec = e - e * (ag[1][14] + ag[2][14] * e + ag[3][14] / e + ag[4][14] / (e * e));
    } else if (theta < ag[0][0]) {
        ec = e - e * (ag[1][0] + ag[2][0] * e + ag[3][0] / e + ag[4][0] / (e * e));
    } else {
        for (int i = 0; i < 15; i++) {
            if (theta >= ag[0][i] - 1 && theta < ag[0][i] + 1) {
                ec = e - e * (ag[1][i] + ag[2][i] * e + ag[3][i] / e + ag[4][i] / (e * e));
                break;
            }
        }
    }
    return ec;
}



#endif