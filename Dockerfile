FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main

RUN apt-get update -y && \
    apt-get install -y curl unzip

RUN apt-get install -y default-jdk
RUN apt-get install -y vim wget git libz-dev
RUN echo "alias ll='ls -l --color=auto'" >> .bashrc
RUN apt-get install -y default-jdk
RUN curl -L https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download -o \
    BBMap_39.01.tar.gz && \
    tar -xvzf BBMap_39.01.tar.gz && \
    rm BBMap_39.01.tar.gz
RUN python3 -m pip install biopython slims-python-api


# cellranger download link, needs to be updated periodically


# Copy barcode files to root dir
COPY bc50.txt /root/bc50.txt
COPY bc50_old.txt /root/bc50_old.txt
COPY bc96.txt /root/bc96.txt

RUN python3 -m pip install matplotlib
RUN python3 -m pip install pandas



#install chromap

RUN wget https://github.com/haowenz/chromap/archive/refs/heads/li_dev4.zip && \
    unzip li_dev4.zip && \
    mv /root/chromap-li_dev4 /root/chromap && \
    cd /root/chromap && \
    make && \
    cd /root 

RUN apt-get update -y
RUN apt-get install -y gdebi-core
RUN apt install -y aptitude
RUN aptitude install -y libjpeg-dev
RUN apt-get update -y

RUN echo "alias ll='ls -l --color=auto'" >> .bashrc
RUN apt-get install tabix

RUN pip install pyranges
RUN pip install pyrle

RUN git clone https://github.com/aertslab/pycisTopic.git && \
    cd /root/pycisTopic && \
    pip install . && \
    cd /root


COPY blacklist /root/blacklist
COPY hg38_chrom_sizes.txt /root/hg38_chrom_sizes.txt
COPY mm10_chrom_sizes.txt /root/mm10_chrom_sizes.txt

#references

#RUN wget http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz && \
#    tar -xvzf Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz && \
#    rm  Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz

#RUN wget http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz && \
#    tar -xvzf Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz && \
#    rm Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz 



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

RUN echo "alias ll='ls -l --color=auto'" >> .bashrc
# Fix systemd conflict with timedatectl
RUN echo "TZ=$( cat /etc/timezone )" >> /etc/R/Renviron.site
# Have to install devtools like this; see https://stackoverflow.com/questions/20923209, also cairo
RUN apt-get install -y r-cran-devtools libcairo2-dev

# Install packages
RUN R -e "install.packages(c('Cairo', 'BiocManager', 'Matrix', 'Seurat','shiny', 'shinyhelper', 'data.table', 'Matrix', 'DT', 'magrittr','ggplot2','ggrepel','hdf5r','ggdendro','gridExtra', 'ggseqlogo', 'circlize','tidyverse','qdap','topGO'))"
RUN R -e "devtools::install_github('immunogenomics/harmony')"
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
RUN R -e "devtools::install_github('GreenleafLab/chromVARmotifs')"
RUN R -e "library('ArchR'); ArchR::installExtraPackages()"

RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "devtools::install_github('SGDDNB/ShinyCell')"
RUN R -e "BiocManager::install('ComplexHeatmap')"


# Upgrade R to version 4.3.0

RUN wget https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
RUN tar zxvf R-4.3.0.tar.gz
RUN cd R-4.3.0 && ./configure --enable-R-shlib
RUN cd R-4.3.0 && make && make install
RUN apt install -y default-jdk
RUN R CMD javareconf

RUN R -e "install.packages(c('pkgconfig','munsell','zip','zoo','xtable','listenv','lazyeval','bit64','rJava','labeling','magick'),repos = 'http://cran.us.r-project.org')"

RUN R -e "ArchR::installExtraPackages()"
RUN R -e "BiocManager::install(version = '3.17', ask = FALSE)"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "install.packages(c('kableExtra','viridis'),repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages(c('rmdformats'),repos = 'http://cran.us.r-project.org')"

RUN R -e "install.packages(c('Signac','Seurat'),repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('rhdf5')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('BSgenome')"
RUN R -e "BiocManager::install('EnsDb.Mmusculus.v79')"
RUN R -e "BiocManager::install('EnsDb.Hsapiens.v86')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('BSgenome.Rnorvegicus.UCSC.rn6')"
RUN R -e "install.packages(c('Seurat'), dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

COPY rn6_chrom_sizes.txt /root/rn6_chrom_sizes.txt
COPY summary.py /root/summary.py
COPY singlecell.py /root/singlecell.py
COPY peak_files.R /root/peak_files.R
COPY version /root/version

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN echo 'hello'
RUN pip install marshmallow-enum
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root


