#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:05
#SBATCH -p priority
#SBATCH --mem-per-cpu=1M
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

cd /n/groups/cbdm-db/dam41/FC_07318

sbatch run_cuttag_alignment.sh LIB055320_CHS00231818_S9 MEC_Hnf4g_fix_R1 mm10 /n/groups/cbdm-db/dam41/FC_07318
sbatch run_cuttag_alignment.sh LIB055320_CHS00231819_S10 MEC_Hnf4g_fix_R2 mm10 /n/groups/cbdm-db/dam41/FC_07318
sbatch run_cuttag_alignment.sh LIB055320_CHS00231820_S11 MEC_Hnf4g_nofix_R3 mm10 /n/groups/cbdm-db/dam41/FC_07318
sbatch run_cuttag_alignment.sh LIB055320_CHS00231821_S12 MEC_Hnf4g_nofix_R4 mm10 /n/groups/cbdm-db/dam41/FC_07318

cd /n/groups/cbdm-db/dam41/FC_07782

sbatch run_cuttag_alignment.sh LIB057739_CHS00253003_S3 IEC_Hnf4g_R1 mm10 /n/groups/cbdm-db/dam41/FC_07782
sbatch run_cuttag_alignment.sh LIB057739_CHS00253004_S4 IEC_Hnf4g_R2 mm10 /n/groups/cbdm-db/dam41/FC_07782
