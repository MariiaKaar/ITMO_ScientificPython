FROM anaconda/miniconda3


RUN conda install -y -c bioconda -c conda-forge \
    python=3.11 \

    bwa \

    samtools \

    bcftools \

    fastqc=0.12.1 \

    trimmomatic \

    matplotlib \

    && conda clean -a -y
