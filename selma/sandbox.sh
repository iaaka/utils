conda activate selma
export PATH=/nfs/cellgeni/pasham/bin/selma.m01/bin:/software/R-4.1.0/bin/:$PATH   
export PYTHONPATH=/nfs/cellgeni/pasham/bin/selma.m01/lib/python3.6/site-packages:$PYTHONPATH

# test on one example
SELMA -m sc --scATAC10x -i  ../E7.5_rep1/atac_fragments.tsv.gz -g mm10 -f PE -o E7.5_rep1_my01 -t ATAC -s mm10.2bit --SCcorrection -p /nfs/team205/vk7/sanger_projects/large_data/gastrulation_multiome_anndata/latest/data/processed/Annotations/all_peaks.bed

# combine fragments
cd ~/nfs/lustre/tmp/atac.bkg.quant
zcat */atac_fragments.tsv.gz | grep -v '#' | gzip > combined_atac_fragments.bed.gz
# bedtools sort -i combined_atac_fragments.bed.gz doesnt work even with 140G RAM
gunzip -c combined_atac_fragments.bed.gz | sort -t $'\t' -k 1,1 -k 2,2n | gzip > combined_atac_fragments.sorted.bed.gz
mkdir combined
mv combined_atac_fragments.sorted.bed.gz combined/atac_fragments.tsv.gz
