#ifndef Analysis
#define Analysis
void Analysis_tagged(string filename = "", string outname = "Tagged_Analysis", bool LeptonID = true, bool All = true)
{
    string outfilename;
    if (All)
    {
        if (LeptonID)
        {
            //lep_s = "_LeptonID";
            Get_histos_t("S19" + filename, 0.0, 0.05);
            Get_histos_t("F18in" + filename, 0.0, 0.05);
            Get_histos_t("F18out" + filename, 0.0, 0.05);
        }
        else
        {
            Get_histos_t("S19" + filename);
            Get_histos_t("F18in" + filename);
            Get_histos_t("F18out" + filename);
        }
    }
    else
    {
        if (LeptonID){
            //lep_s = "_LeptonID";
            Get_histos_t(filename, 0.0, 0.05);
        }
        else
            Get_histos_t(filename);
    }

    outfilename=outname;

    plot_results_t("./R_ALL/" + outfilename);
}

void Analysis_exclusive(string filename = "", string outname = "Exclusive_Analysis", bool LeptonID = true, bool All = true)
{
    string outfilename;
    string lep_s = "";
    string folder="";
    if (All)
    {
        if (LeptonID)
        {
            lep_s = "_LeptonID";
            Get_histos_exclusive("S19" + filename + "", 0.0, 0.05);
            Get_histos_exclusive("F18in" + filename + "", 0.0, 0.05);
            Get_histos_exclusive("F18out" + filename + "", 0.0, 0.05);
        }
        else
        {
            Get_histos_exclusive("S19" + filename + "");
            Get_histos_exclusive("F18in" + filename + "");
            Get_histos_exclusive("F18out" + filename + "");
        }
    }
    else
    {
        if (LeptonID){
            lep_s = "_LeptonID";
            Get_histos_exclusive(folder+filename, 0.0, 0.05);
        }
        else
            Get_histos_exclusive(folder+filename);
    }
    outfilename=outname;

    plot_results_exc("./R_others/Jan2025/" + outfilename);
}

void Analysis_onelp(string filename = "", int top=3, string outname = "epe+_Analysis",Double_t AIcut=0.04, bool LeptonID = true, bool All = true)
{
    string outfilename;
    string folder="";
    if (All)
    {
        if (LeptonID)
        {
            Get_histos_onelp("S19" + filename ,top, 0.0,0.05,AIcut);
            Get_histos_onelp("F18in" + filename ,top, 0.0,0.05,AIcut);
            Get_histos_onelp("F18out" + filename,top, 0.0,0.05,AIcut);
        }
        else
        {
            Get_histos_onelp("S19" + filename,top );
            Get_histos_onelp("F18in" + filename,top);
            Get_histos_onelp("F18out" + filename,top);
        }
    }
    else
    {
        if (LeptonID){
            Get_histos_onelp(folder+filename,top, 0.0, 0.05);
        }
        else
            Get_histos_onelp(folder+filename,top);
    }
    outfilename=outname;

    plot_results_onelp("./R_others/Jan2025/" + outfilename,top);
}



void Analysis_untagged(string filename = "All", string outfilename = "Untagged_Analysis")
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