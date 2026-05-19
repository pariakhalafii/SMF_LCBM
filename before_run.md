# to get started:

#install miniforge that gives me conda and mamba 

brew install miniforge
   conda init zsh    # or `bash`, depending on your shell


# close and open terminal again 


# create an env with everything that this pipeline needs

conda create -n smf -c bioconda -c conda-forge \
       flash fastqc trim-galore cutadapt bismark bowtie2 samtools \
       python=3.11 numpy pandas matplotlib pysam pyfaidx pytest pyyaml







# activate it whenever I need to run the pipeline 
conda activate smf


# verify it works

flash --version
fastqc --version
trim_galore --version
bismark --version
samtools --version