#! /bin/bash -e
#BSUB -G cellgeni
#BSUB -J selma[2]
#BSUB -o %J.%I.selma.log
#BSUB -e %J.%I.selma.err
#BSUB -q long
#BSUB -n 2
#BSUB -M250000
#BSUB -R "span[hosts=1] select[mem>250000] rusage[mem=250000]"

cd /lustre/scratch117/cellgen/cellgeni/pasham/tmp/atac.bkg.quant/selma

TAGs=(E7.5_rep1 E7.5_rep2 E7.75_rep1 E8.0_rep1 E8.0_rep2 E8.5_CRISPR_T_KO E8.5_CRISPR_T_WT E8.5_rep1 E8.5_rep2 E8.75_rep1 E8.75_rep2 combined)
TAG=${TAGs[$LSB_JOBINDEX-1]}

source activate selma
export PATH=/nfs/cellgeni/pasham/bin/selma.m01/bin:/software/R-4.1.0/bin/:$PATH   
export PYTHONPATH=/nfs/cellgeni/pasham/bin/selma.m01/lib/python3.6/site-packages:$PYTHONPATH


SELMA -m sc \
 --scATAC10x -i  \
 ../${TAG}/atac_fragments.tsv.gz \
 -g mm10 \
 -f PE \
 -o ${TAG} \
 -t ATAC \
 -s mm10.2bit \
 --SCcorrection \
 -p /nfs/team205/vk7/sanger_projects/large_data/gastrulation_multiome_anndata/latest/data/processed/Annotations/all_peaks.bed
