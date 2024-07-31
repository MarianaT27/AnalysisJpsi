#!/bin/bash

#---------------FOR EFFICIENCY STUDIES USING MC SIMULATIONS------------------
#clas12root4 'taggedMC.C("MC_Rec_5830","Eff_5830", 5830,199, 4.0,4.5)' -q -b
#clas12root4 'lund.C("MC_Gen_5830","Eff_5830", 5830,199, 4.0,4.5)' -q -b


#---------------FOR CLASSIC TAGGED AND TAGGED MU ANALYSIS----------
#-------------tagged.C("NameFile",version, beam energy, pass, #max)
#clas12root 'tagged.C("DCSmearing_Data_t",-18,10.6,2, 10)' -q -b 
#clas12root 'tagged.C("DCSmearing_v5.4_t",-18,10.6,2, 10,0)' -q -b 
#clas12root4 'taggedmu.C("S19_p2_muons_G",-19,10.2,2,10)' -q -b 
#clas12root4 'untaggedmu.C("S19_p2_muons",-19,10.2,2,119)' -q -b 


 
#---------------To create a root file -------------------------
#clas12root4 'toroot_MC.C("F18in_noBg",7324,-18,10.6)' -q -b
#clas12root4 'toroot_MC.C("F18out_noBg",7325,+18,10.6)' -q -b
#clas12root4 'toroot_MC.C("S19_noB",7326,-19,10.2)' -q -b
#clas12root 'toroot.C("DCSmearing_v5.9",-18,10.6,0,10,"_v5.9")' -q -b
#clas12root 'toroot.C("F18in_7772",-18,10.6,0,200,"",7772)' -q -b
#clas12root 'toroot.C("F18in_7773",-18,10.6,0,61,"test_nobg-",7773)' -q -b
#clas12root 'toroot.C("F18in_7878-80",-18,10.6,4,100,"F18",7878)' -q -b
#clas12root 'toroot.C("F18in_SIDIS300",-18,10.6,3,50,"SIDIS",7900)' -q -b

#clas12root -l -b -q Hipo_Analysis.C /volatile/clas12/osg/marianat/job_7898/output/* F18in_JPSIGEN_OSG_v2
#clas12root -l -b -q Hipo_Analysis.C /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/jpsitcs/jpsitcs_0050* F18in_JPSIGEN_Data
#clas12root -l -b -q Hipo_Analysis.C /volatile/clas12/osg/marianat/job_7900/output/* F18in_SIDIS_OSG
#clas12root -l -b -q Hipo_Analysis.C /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_0050* F18in_SIDIS_Data_v2

#clas12root -l -b -q Bg_est.C /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/jpsitcs/jpsitcs_005* F18in
#clas12root -l -b -q Bg_est.C /cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/jpsitcs/jpsitcs_00* F18out
#clas12root -l -b -q Bg_est.C /volatile/clas12/osg/marianat/job_7900/output/* F18in_SIDIS_OSG


#clas12root 'toroot.C("F18in_temp",-18,10.6,2,15)' -q -b


#clas12root 'toroot.C("DCSmearing_Data",-18,10.6,2,10)' -q -b

#clas12root 'untagged.C("F18in_osg",-18,10.6,"_ut",1,0)' -q -b
#clas12root 'untagged.C("DCSmearing_v5.4",-18,10.6,"_ut",1,0)' -q -b
#clas12root 'untagged.C("DCSmearing_Data",-18,10.6,"_ut",1,0)' -q -b

#clas12root 'Results.C("F18in_osg","OSG")' -q -b
#clas12root 'Results.C("DCSmearing_v5.9_ut","New")' -q -b
#clas12root 'Results.C("DCSmearing_v5.4_ut","Old")' -q -b

#clas12root 'toroot.C("F18out",-18,10.6,2)' -q -b
#--------------------------------------------------------------
#	(S)-silver runs, max=2	; (G)-golden runs, max=10  ;   //Spring19
#			Full-all funs, max=119		       //Spring 19


#----------------FOR DOING ANALYSIS FROM ROOT FILE--------------
#	(S)-silver runs (G)-golden runs 
#	PASS1
#	Spring 2019 MC ..........................S19_MC.root 
#	Spring 2019 tagged J/psi->e-e+...........S19_p1_Full.root
#
#	PASS2    
#	Spring 2019 tagged J/psi->e-e+...........S19_p2_Full.root  
#	Spring 2019 tagged J/psi->mu-mu+.........S19_p2_mu.root
#	Spring 2019(S) tagged J/psi->mu-mu+......S19_p2_mu_S.root
#
#

#clas12root4 'taggedfromroot.C("osg_new",-19,10.2,2)' -q -b

#PASS=0 is to test corrections
clas12root 'taggedfromroot.C("F18in",-18,10.6)' -q -b
clas12root 'taggedfromroot.C("F18out",+18,10.6)' -q -b
clas12root 'taggedfromroot.C("S19",-19,10.2)' -q -b
clas12root 'taggedfromroot.C("F18in",-18,10.6,0)' -q -b
clas12root 'taggedfromroot.C("F18out",+18,10.6,0)' -q -b
clas12root 'taggedfromroot.C("S19",-19,10.2,0)' -q -b
#clas12root4 'untaggedfromroot.C("DCSmearing_new",-19,10.2,2,"untagged")' -q -b

#clas12root 'untagged.C("S19",-19,10.2,"_ut")' -q -b
#clas12root 'untagged.C("F18in",-18,10.6,"_TEST",1,0)' -q -b
#clas12root 'untagged.C("F18out",+18,10.6,"_ut")' -q -b




#clas12root FitResults.C -q -b
#clas12root 'Results.C("UT_BDT6_zero_pos",6,0.0)' -q -b
#clas12root 'Results.C("UT_BDT9_zero_pos",9,0.0)' -q -b
#clas12root 'Results.C("UT_BDT6_OutEvents_pos",-6,0.0)' -q -b
#clas12root 'Results.C("UT_BDT9_OutEvents_pos",-9,0.0)' -q -b
#clas12root 'Results.C("UT_noML")' -q -b
#clas12root 'Results.C("noFC_NOScore",0)' -q -b


#clas12root 'Results.C("F18in_temp")' -q -b
#clas12root 'Results.C()' -q -b

#clas12root 'FitResults.C("Cut2",-0.030)' -q -b
#clas12root 'FitResults.C("Cut3",-0.025)' -q -b
#clas12root 'FitResults.C("Cut4",-0.020)' -q -b
#clas12root 'FitResults.C("Cut5",-0.015)' -q -b
#evince Cut*.pdf


