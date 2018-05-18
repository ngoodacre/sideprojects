#$ -cwd
#$ -S /bin/sh
#                                   $ -l h_vmem=10G
#$ -l h_rt=40:00:00
#$ -j y
#$ -N Syn3.0_BLASTP_vs_UNIPROTSPROT
#$ -l bigbox

#$ -pe thread 40

NT=8

export PATH=/home/mikem/python/bin:$PATH
export PATH=/home/mikem/blast+2.3.0_fda/bin:$PATH
export PATH=/home/mikem/muscle:$PATH  


echo "Running task $SGE_TASK_ID of job $JOB_NAME (ID of $JOB_ID) on $HOSTNAME"

DIR=/projects01/arifakhanlabs/BLAST_runs

time blastp -query $DIR/Syn3.0.fasta -db outdir/uniprot_trembl -evalue 0.00001 -dbseqnum 63039659 -out $DIR/Syn3.0_BLASTP_vs_UNIPROTTREMBL_output.txt 


# time tblastx -query $QF/VDBv5.1.insertionseqs.fasta -db nr -out HMnrblastxSF9LRassembly-1.fa -evalue 0.001 -outfmt 0 -num_threads $NT

