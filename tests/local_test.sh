#!/bin/bash
## activate nf-core conda environment
env_nf=$HOME/miniforge3/envs/m1/env_nf
source $HOME/miniforge3/bin/activate ${env_nf}
## specify params
outdir=$HOME/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/SarcAtlas/APS033_rnaseq_qc/nanoqc_test/filtered_results
pipelinedir=$HOME/VSCodeProjects/shahcompbio-nanoqc
samplesheet=${pipelinedir}/assets/local_samplesheet.csv
refgenome=$HOME/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/code/ref_genomes/hg38p14/GRCh38.primary_assembly.genome.fa
gtf=$HOME/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/code/localpgdata/resources/gencode.v31.primary_assembly.annotation.gtf
mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile arm,docker \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --input ${samplesheet} \
    --fasta ${refgenome} \
    --gtf ${gtf} \
    --filter_reads \
    --protocol 'RNA'