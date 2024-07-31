#ifndef Analysis
#define Analysis
void Analysis_tagged(string filename="", string outfilename="Tagged_Analysis", bool All=true, bool LeptonID=true)
{
    if (All)
    {
        if(LeptonID){
            Get_histos_t("S19_t");
            Get_histos_t("F18in_t");
            Get_histos_t("F18out_t");
        }
        else{
            Get_histos_t("S19_t", 0.0, 0.05);
            Get_histos_t("F18in_t", 0.0, 0.05);
            Get_histos_t("F18out_t", 0.0, 0.05);
        }
    }
    else
        Get_histos_t(filename);

    plot_results_t(outfilename);
}

void Analysis_untagged(string filename="All", string outfilename="Untagged_Analysis")
{
    if (filename == "All")
    {
        Get_histos("S19_ut", 0.0, 0.05);
        Get_histos("F18in_ut", 0.0, 0.05);
        Get_histos("F18out_ut", 0.0, 0.05);
    }
    else
        Get_histos(filename);

    plot_results(outfilename);
}

#endif