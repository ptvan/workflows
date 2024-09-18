FROM condaforge/miniforge3:latest
ENV PATH=/opt/conda/bin:$PATH
RUN mamba update -n base -c defaults mamba
RUN mamba install -c bioconda deeptools && \
    mamba install -c bioconda macs3

RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker
RUN apt-get update && apt-get install -y bedtools && \
    apt-get install -y bowtie2 && \
    apt-get install -y samtools && \
    apt-get install -y python3


