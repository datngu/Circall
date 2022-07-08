FROM ubuntu:22.04

MAINTAINER Dat T Nguyen "ndat<at>utexas.edu"
LABEL authors="Dat T Nguyen" \
      description="Docker image containing all requirements for running Circall" 


ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update -y
RUN apt install -yq r-base

RUN R -e 'install.packages("data.table")'
RUN R -e 'install.packages("foreach")'
RUN R -e 'install.packages("doParallel")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'BiocManager::install("Biostrings")'

ADD Circall_v0.2.0 /Circall
ENV LD_LIBRARY_PATH /Circall/linux/lib:$LD_LIBRARY_PATH
ENV PATH /Circall/linux/bin:$PATH

WORKDIR /Circall
RUN bash ./config.sh
WORKDIR /

##conda - not use
#ADD environment.yml /
#RUN conda install mamba -n base -c conda-forge -y
#RUN mamba env create -f /environment.yml
#ENV PATH /opt/conda/envs/circall/bin:$PATH




