#!/bin/bash
#SBATCH --job-name=toroot                   # Job name
#SBATCH --partition=ifarm
#SBATCH --output=/farm_out/marianat/toroot.out     # Output file (no extra space before path)
#SBATCH --error=/farm_out/marianat/toroot.err      # Error file (no extra space before path)
#SBATCH --time=4:00:00                              # Max run time (hh:mm:ss)
#SBATCH --mem=4G                                     # Memory request
#SBATCH --ntasks=1                                # Number of tasks (usually 1 for ROOT scripts)
#SBATCH --cpus-per-task=1                            # Number of CPUs per task
#SBATCH --mail-user=marianat@jlab.org
#SBATCH --mail-type=END

# Load necessary modules (if applicable)
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles   
module load clas12 

clas12root 'toroot_v2.C("S18_10GeV_In1",-10,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_inb_10GeV_physics_1/dst/train/jpsitcs/")' -q -b
clas12root 'toroot_v2.C("S18_10GeV_In2",-10,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_inb_10GeV_physics_2/dst/train/jpsitcs/")' -q -b
clas12root 'toroot_v2.C("S18_10GeV_Out",10,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_outb_10GeV_physics/dst/train/jpsitcs/")' -q -b
clas12root 'toroot_v2.C("S18_6GeV_In",-6,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_inb_6GeV_physics/dst/train/jpsitcs/")' -q -b
clas12root 'toroot_v2.C("S18_6GeV_Out",6,"/volatile/clas12/rg-a/production/pass0/Spring18/dst/rga_sp18_outb_6GeV_physics/dst/train/jpsitcs/")' -q -b

clas12root 'taggedfromroot.C("S18_10GeV_In",-10,0)' -q -b
clas12root 'taggedfromroot.C("S18_10GeV_Out",10,0)' -q -b
clas12root 'taggedfromroot.C("S18_6GeV_In",-6,0)' -q -b
clas12root 'taggedfromroot.C("S18_6GeV_Out",6,0)' -q -b

