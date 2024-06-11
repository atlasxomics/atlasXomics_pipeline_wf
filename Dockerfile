FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

# install background packages
RUN apt-get update -y && \
    apt-get install -y \
        aptitude \
        curl \
        default-jdk \
        gdebi-core \
        git \
        libjpeg-dev \
        libz-dev \
        unzip \
        tabix \
        vim \
        wget

# install BBMAP
RUN curl -L https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download -o \
    BBMap_39.01.tar.gz && \
    tar -xvzf BBMap_39.01.tar.gz && \
    rm BBMap_39.01.tar.gz

# install Chromap
RUN git clone -b li_dev4 https://github.com/haowenz/chromap.git &&\
    cd /root/chromap && \
    make && \
    cd /root 

# Install python packages
RUN pip install --upgrade pip
COPY requirements.txt /root/requirements.txt
RUN python3 -m pip install -r requirements.txt

# Install modified pycisTopic from commit e9b0e1a
COPY pycisTopic.tar.gz /root/pycisTopic.tar.gz
RUN tar -xzvf pycisTopic.tar.gz && \
    cd /root/pycisTopic && \
    pip install . && \
    cd /root

# Copy support dirs 
COPY barcodes /root/barcodes
COPY blacklist /root/blacklist
COPY chrom_sizes /root/chrom_sizes
COPY scripts /root/scripts

COPY version /root/version

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN pip install marshmallow-enum
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
