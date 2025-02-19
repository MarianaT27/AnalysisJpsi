#!/bin/bash

#---------------FOR EFFICIENCY STUDIES USING MC SIMULATIONS------------------
#clas12root4 'taggedMC.C("MC_Rec_5830","Eff_5830", 5830,199, 4.0,4.5)' -q -b
#clas12root4 'lund.C("MC_Gen_5830","Eff_5830", 5830,199, 4.0,4.5)' -q -b

#---------------FOR CLASSIC TAGGED AND TAGGED MU ANALYSIS----------
#-------------tagged.C("NameFile",version, beam energy, pass, #max)
#clas12root 'tagged.C("DCSmearing_Data_t",-18,10.6,2, 10)' -q -b 
#clas12root4 'taggedmu.C("S19_p2_muons_G",-19,10.2,2,10)' -q -b 
#clas12root4 'untaggedmu.C("S19_p2_muons",-19,10.2,2,119)' -q -b 


#---------------TO CREATE A ROOT FILE -------------------------

#clas12root -l -b -q Hipo_Analysis.C /volatile/clas12/osg/marianat/job_7898/output/* F18in_JPSIGEN_OSG_v2
#clas12root -l -b -q Bg_est.C /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/jpsitcs/jpsitcs_005* F18in
#clas12root 'toroot_v2.C("F18in_Full_FTJPsi",18,10.6,2)' -q -b
#clas12root 'toroot_v2.C("S19_Full",-19,10.2,2)' -q -b


#---------------TO READ A ROOT FILE -------------------------
#clas12root4 'untaggedfromroot.C("DCSmearing_new",-19,10.2,2,"untagged")' -q -b
#clas12root 'untagged.C("S19",-19,10.2,"_ut")' -q -b
#root 'compareRootFiles.C("S18_10GeV_In_ee+e-","F18in_ee+e-","10GeV_Inbending")' -q
#clas12root 'toroot_v2.C("S18_10GeV_In1",-10,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_inb_10GeV_physics_1/dst/train/jpsitcs/")' -q -b
#clas12root 'toroot_v2.C("S18_10GeV_In2",-10,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_inb_10GeV_physics_2/dst/train/jpsitcs/")' -q -b
#clas12root 'toroot_v2.C("S18_10GeV_Out",10,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_outb_10GeV_physics/dst/train/jpsitcs/")' -q -b
#clas12root 'toroot_v2.C("F18in",10.6,"",2,-18,7)' -q -b
#clas12root 'toroot_v2.C("F18out",10.6,"",2,18,7)' -q -b
#clas12root 'toroot_v2.C("S18_6GeV_In",-6,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_inb_6GeV_physics/dst/train/jpsitcs/")' -q -b
#clas12root 'toroot_v2.C("S18_6GeV_Out",6,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_outb_6GeV_physics/dst/train/jpsitcs/")' -q -b


#clas12root 'taggedfromroot.C("S18_10GeV_In",-10,0)' -q -b
#clas12root 'taggedfromroot.C("F18in",-18,0)' -q -b
root 'compareRootFiles.C("S18_10GeV_In_ee+e-","F18in_ee+e-","10GeV_Inbending")' -q
#clas12root 'taggedfromroot.C("F18out",+18,0)' -q -b
#clas12root 'taggedfromroot.C("S18_10GeV_Out",10,0)' -q -b
root 'compareRootFiles.C("S18_10GeV_Out_ee+e-","F18out_ee+e-","10GeV_Outbending")' -q
#clas12root 'taggedfromroot.C("S18_6GeV_In",-6,0)' -q -b
#clas12root 'taggedfromroot.C("S18_6GeV_Out",6,0)' -q -b

#clas12root 'exclusivefromroot.C("F18in_Full",-18)' -q -b

#string nameFile = "S19", int version = -19, int top = 3, int PASS = 2,bool Restricted=false,bool AI = true,TString AIModel="12", bool QA = false
#top=3 is e+e'p' top=4 is e-e'p'
#Restricted means selecting only 1 lepton, 1 e' and 1 p'
#QA must be false because is not working 
#clas12root 'onelep_fromroot.C("S19_Full",-19,4,2,false,true,"12")' -q -b
#clas12root 'onelep_fromroot.C("S19_Full",-19,4,2,false,true,"10")' -q -b
#root 'Variable_Plots_Tagged.C("F18in_10_vs_12","F18in_Full_epe-_","")' -q
#root 'Variable_Plots_Tagged.C("F18out_10_vs_12","F18out_Full_epe-_","")' -q
#root 'Variable_Plots_Tagged.C("S19_10_vs_12","S19_Full_epe-_","")' -q
#clas12root 'onelep_fromroot.C("F18in_Full",-18,4,2,false,true,"12_MIX")' -q -b


#root 'Variable_Plots_Tagged.C("Var10","Var10")' -q
#root 'Variable_Plots_Tagged.C("Var9_v2","Var9_v2")' -q
#root 'Variable_Plots_Tagged.C("Var9_v4","Var9_v4")' -q
#clas12root 'onelep_fromroot.C("MC_F18",-18,4,2)' -q -b
#clas12root 'onelep_fromroot.C("MC_F18",-18,3,2)' -q -b
#root 'evaluate_TMVA.C("Var9_0-100_MC_epe+","Var9_0-100","MC_F18_epe+")' -q
#root 'evaluate_TMVA.C("Var9_0-100_MC_epe-","Var9_0-100","MC_F18_epe-")' -q
#root 'evaluate_TMVA.C("Var9_0-100_MC_exc","Var9_0-100","MC_F18_exclusive")' -q

#root 'Variable_Plots_Tagged.C("Var9_0-100_MC_epe+","Var9_0-100_MC_epe+",-0.2)' -q
#root 'Variable_Plots_Tagged.C("Var9_0-100_MC_exc","Var9_0-100_MC_exc",-0.2)' -q
#root 'evaluate_TMVA.C("Var9_0-100_MC_epe+","Var9_0-100","MC_F18_epe+")' -q
#root 'evaluate_TMVA.C("Var9_50-50_allSIG_accidentals","Var9_50-50_allSIG","F18in_Full_epe-_accidentals")' -q

#root 'evaluate_TMVA.C("1")'



#clas12root 'Results.C(2,"_")' -q -b -l
#clas12root 'Results.C(2,"_","AI")' -q -b -l
#clas12root 'Results.C(3,"_","AI")' -q -b -l


#clas12root 'Results.C(4,"_","_c_cut18",0.20)' -q -b -l
#clas12root 'Results.C(4,"_Restricted_","_R_c_cut14",0.0)' -q -b -l
#clas12root 'Results.C(4,"_Restricted_","_R_c_cut17",0.15)' -q -b -l




#clas12root 'Results.C(1,"_Restricted_","_R")' -q -b -l
#clas12root 'Results.C(2,"_Restricted_","_R")' -q -b -l
#clas12root 'Results.C(3,"_Restricted_","_R")' -q -b -l


#----------------------- TO RUN DEC/03/2024 ------------------------------
#clas12root 'Results.C(6,"F18in_All_epe-_","epe-")' -q -b -l
#clas12root 'Results.C(6,"F18in_Full_epe-_","epe-_restricted")' -q -b -l3

#clas12root 'Results.C(6,"F18in_All_epe-_10","epe-_10")' -q -b -l
#clas12root 'Results.C(6,"F18in_Full_epe-_10","epe-_restricted_10")' -q -b -l


#clas12root 'Results.C(6,"F18in_All_epe-_10_MIX","epe-_10_MIX")' -q -b -l
#clas12root 'Results.C(6,"F18in_Full_epe-_10_MIX","epe-_restricted_10_MIX")' -q -b -l

#clas12root 'Results.C(6,"F18in_All_epe-_10_SIDIS","epe-_10_SIDIS")' -q -b -l
#clas12root 'Results.C(6,"F18in_Full_epe-_10_SIDIS","epe-_restricted_10_SIDIS")' -q -b -l






#clas12root FitResults.C -q -b



#clas12root 'Results.C("F18in_temp")' -q -b
#clas12root 'Results.C(0)' -q -b -l

#clas12root 'Results.C(1)' -q -b -l
#clas12root 'Results.C(2)' -q -b -l
#clas12root 'Results.C(3)' -q -b -l
#clas12root 'Results.C(4)' -q -b -l




#clas12root 'FitResults.C("Cut2",-0.030)' -q -b
#clas12root 'FitResults.C("Cut3",-0.025)' -q -b
#clas12root 'FitResults.C("Cut4",-0.020)' -q -b
#clas12root 'FitResults.C("Cut5",-0.015)' -q -b
#evince Cut*.pdf


