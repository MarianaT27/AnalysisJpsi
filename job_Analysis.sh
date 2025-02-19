#!/bin/bash
#SBATCH --job-name=exclusivefromrootv2                   # Job name
#SBATCH --partition=ifarm
#SBATCH --output=/farm_out/marianat/exclusivefromrootv2.out     # Output file (no extra space before path)
#SBATCH --error=/farm_out/marianat/exclusivefromrootv2.err      # Error file (no extra space before path)
#SBATCH --time=4:00:00                              # Max run time (hh:mm:ss)
#SBATCH --mem=4G                                     # Memory request
#SBATCH --ntasks=1                                # Number of tasks (usually 1 for ROOT scripts)
#SBATCH --cpus-per-task=1                            # Number of CPUs per task
#SBATCH --mail-user=marianat@jlab.org
#SBATCH --mail-type=END

# Load necessary modules (if applicable)
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles   
module load clas12 

clas12root 'exclusivefromrootv2.C("F18in","SIG",-18,2,false)' -q -b
#root 'evaluate_TMVA.C("Var9_0-100_exc","Var9_0-100","F18in_exclusive")' -q

clas12root 'exclusivefromrootv2.C("F18out","SIG",18,2,false)' -q -b
#clas12root 'onelep_fromroot.C("MC_F18",-18,4,2)' -q -b
#clas12root 'onelep_fromroot.C("MC_F18",-18,3,2)' -q -b
#root 'evaluate_TMVA.C("Var9_0-100_exc","Var9_0-100","F18in_exclusive")' -q


clas12root 'exclusivefromrootv2.C("S19","SIG",-19,2,false)' -q -b