#$ -cwd
#$ -l h_rt=001:00:00
#                       $ -l h_vmem=16G
#$ -S /bin/sh
#$ -j y
#$ -N CCMpred
#$ -pe thread 8 
 
echo "Running job $JOB_NAME, $JOB_ID on $HOSTNAME"

APP=/projects/mikem/CCMpred/bin/ccmpred
INP=/projects01/arifakhanlabs/CCMPred/ccmpred_aln/COG1925_COG1105_byspecies_paired.afa.ccmp
OUT=/projects01/arifakhanlabs/CCMPred/output/COG1925_COG1105_byspecies_paired.mat


time $APP -t $NSLOTS $INP $OUT 

