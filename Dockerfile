FROM condaforge/mambaforge:4.9.2-5 AS conda

# Ignore interactive prompts
ARG DEBIAN_FRONTEND=noninteractive

# update packages and install make
RUN apt-get update && apt-get install make

# install open-ssh
RUN apt-get install -y openssh-client

# install cmake and build-essentials
RUN apt-get -y install cmake build-essential

# Add environment lock file
ADD conda-lock.yml /tmp/conda-lock.yml

# install conda-lock
RUN mamba install conda-lock=2.5.7

# update mamba
RUN mamba update -n base mamba

# create a conda env
ENV CONDA_ENV=./conda/bin
RUN conda-lock install --name myenv /tmp/conda-lock.yml

RUN echo "source activate myenv" > ~/.bashrc
ENV PATH=/opt/conda/envs/myenv/bin:$PATH

ADD . /
RUN pip install --upgrade pip
RUN pip install .
