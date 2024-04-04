#!/bin/bash

#---------------FOR EFFICIENCY STUDIES USING MC SIMULATIONS------------------
#clas12root4 'taggedMC.C("MC_Rec_5830","Eff_5830", 5830,199, 4.0,4.5)' -q -b
#clas12root4 'lund.C("MC_Gen_5830","Eff_5830", 5830,199, 4.0,4.5)' -q -b


#---------------FOR CLASSIC TAGGED AND TAGGED MU ANALYSIS----------
#-------------tagged.C("NameFile",version, beam energy, pass, #max)
#clas12root4 'tagged.C("DCSmearing_test_tagged",-19,10.2,2, 10)' -q -b 
#clas12root4 'taggedmu.C("S19_p2_muons_G",-19,10.2,2,10)' -q -b 
#clas12root4 'untaggedmu.C("S19_p2_muons",-19,10.2,2,119)' -q -b 


 
#---------------To create a root file -------------------------
#clas12root4 'toroot_MC.C("F18in_noBg",7324,-18,10.6)' -q -b
#clas12root4 'toroot_MC.C("F18out_noBg",7325,+18,10.6)' -q -b
#clas12root4 'toroot_MC.C("S19_noB",7326,-19,10.2)' -q -b
#clas12root4 'toroot.C("F18in_pass2",-18,10.6,2)' -q -b
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
#clas12root4 'untaggedfromroot.C("DCSmearing_new",-19,10.2,2,"untagged")' -q -b
#clas12root 'untaggedfromroot.C("S19",-19,10.2,2,"_noQ2_QADB")' -q -b


#clas12root FitResults.C -q -b
clas12root 'FitResults.C("Cut1",-0.035)' -q -b
clas12root 'FitResults.C("Cut2",-0.030)' -q -b
clas12root 'FitResults.C("Cut3",-0.025)' -q -b
clas12root 'FitResults.C("Cut4",-0.020)' -q -b
clas12root 'FitResults.C("Cut5",-0.015)' -q -b
evince Cut*.pdf


