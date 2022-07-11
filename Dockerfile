FROM nfcore/base:2.1

MAINTAINER Dat T Nguyen "ndat<at>utexas.edu"
LABEL authors="Dat T Nguyen" \
      description="Docker image containing all requirements for running Circall" 


ADD environment.yml /
ADD Circall_v1.0.1 /Circall


RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/circall/bin:$PATH
# circall
ENV LD_LIBRARY_PATH /Circall/linux/lib:$LD_LIBRARY_PATH
ENV PATH /Circall/linux/bin:$PATH

WORKDIR /Circall
RUN bash ./config.sh
WORKDIR /
