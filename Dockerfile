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

# Install pycisTopic, before it required python3.11; pin scipy install
RUN git clone https://github.com/aertslab/pycisTopic.git && \
    cd /root/pycisTopic && \
    git checkout e9b0e1a && \
    sed -i 's/scipy/scipy==1.12/' requirements.txt && \
    pip install . && \
    cd /root

# Install R
RUN apt-get update -y && \
    apt-get install -y \
        r-base \
        r-base-dev \
        apt-transport-https \
        build-essential \
        gfortran \
        libhdf5-dev \
        libatlas-base-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libgit2-dev \
        libglpk-dev \
        libgsl-dev \
        libicu-dev \
        liblzma-dev \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libpcre3-dev \
        libssl-dev \
        libtcl8.6 \
        libtiff5 \
        libtk8.6 \
        libxml2-dev \
        libxt-dev \
        libx11-dev \
        libtiff-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        locales \
        make \
        pandoc \
        tzdata \
        vim \
        wget \
        zlib1g-dev \
        r-cran-rjava \
        libmagick++-dev

# Fix systemd conflict with timedatectl
RUN echo "TZ=$( cat /etc/timezone )" >> /etc/R/Renviron.site

# Install devtools (https://stackoverflow.com/questions/20923209), also cairo
RUN apt-get install -y r-cran-devtools libcairo2-dev

# Upgrade R to version 4.3.0
RUN wget https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
RUN tar zxvf R-4.3.0.tar.gz
RUN cd R-4.3.0 && ./configure --enable-R-shlib
RUN cd R-4.3.0 && make && make install
RUN R CMD javareconf

# Installation of R packages with renv
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/renv_1.0.5.tar.gz', repos = NULL, type = 'source')"
COPY renv.lock /root/renv.lock
COPY .Rprofile /root/.Rprofile
RUN mkdir /root/renv
COPY renv/activate.R /root/renv/activate.R
COPY renv/settings.json /root/renv/settings.json
RUN R -e "renv::restore()"

# Copy files to root dir
COPY bc50.txt /root/bc50.txt
COPY bc50_old.txt /root/bc50_old.txt
COPY bc96.txt /root/bc96.txt
COPY bc96_fg.txt /root/bc96_fg.txt

COPY blacklist /root/blacklist
COPY hg38_chrom_sizes.txt /root/hg38_chrom_sizes.txt
COPY mm10_chrom_sizes.txt /root/mm10_chrom_sizes.txt

COPY rn6_chrom_sizes.txt /root/rn6_chrom_sizes.txt
COPY singlecellsummary.py /root/singlecellsummary.py
COPY pycis.py /root/pycis.py
COPY peak_files.R /root/peak_files.R
COPY version /root/version
COPY bc_process_newbulk.py /root/bc_process_newbulk.py


# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN pip install marshmallow-enum
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
