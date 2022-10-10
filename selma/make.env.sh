conda create -y -n selma python=3.6
conda activate selma
conda config --env --add channels conda-forge
conda install numpy

# original version
#  unfortuntelly it doesn't work
# git clone https://github.com/Tarela/SELMA.git
# cd SELMA
# python setup.py install --prefix /nfs/cellgeni/pasham/bin/selma
# 
# export PATH=/nfs/cellgeni/pasham/bin/selma/bin:/software/R-4.1.0/bin/:$PATH   
# export PYTHONPATH=/nfs/cellgeni/pasham/bin/selma/lib/python3.6/site-packages:$PYTHONPATH
# SELMA -m sc --scATAC10x -i  ../E7.5_rep1/atac_fragments.bed.gz -g mm10 -f PE -o E7.5_rep1      -t ATAC -s mm10.2bit --SCcorrection -p /nfs/team205/vk7/sanger_projects/large_data/gastrulation_multiome_anndata/latest/data/processed/Annotations/all_peaks.bed



# use my fork
git clone https://github.com/iaaka/SELMA.git
cd SELMA
python setup.py install --prefix /nfs/cellgeni/pasham/bin/selma.m01