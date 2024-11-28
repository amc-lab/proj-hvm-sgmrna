#!/bin/bash

# Process HVM Amplicon Data
# A. M. Chakrabarti
# 24th May 2024

#SBATCH --job-name="hvm_amplicon"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --output=hvm_amplicon-%A.out

ml purge
ml Nextflow/23.10.0
ml Singularity/3.6.4
ml Graphviz/2.47.2-GCCcore-10.3.0

export NXF_SINGULARITY_CACHEDIR=/camp/lab/bauerd/home/shared/singularity
export NXF_HOME=/nemo/lab/ulej/home/users/luscomben/users/chakraa2/.nextflow

cd /camp/lab/bauerd/home/users/chakraa2/projects/harriet/202405

nextflow run main.nf \
-resume \
--input samplesheet_full.csv \
--fasta ref/amplicon.fa \
--outdir results_full