FROM ubuntu:20.04
WORKDIR /pipeline/packages/

SHELL ["bash", "-l" ,"-c"]
RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker
RUN apt update && apt install -y wget bzip2 openjdk-21-jre-headless bowtie2 samtools bedtools && \
    wget --no-check-certificate  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh 
RUN wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar -O /pipeline/packages/picard.jar
RUN bash miniforge.sh -b -p /root/miniforge
ENV PATH="$PATH:/root/miniforge/bin"
RUN mamba update -n base -c defaults mamba
RUN mamba install -c bioconda deeptools && \
    mamba install -c bioconda macs3
