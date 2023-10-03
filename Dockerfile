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

# Copy helper script, rat reference genome
#COPY bc_process.py /root/bc_process.py
#COPY extractBarcode.py /root/extractBarcode.py

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


COPY version /root/version
COPY blacklist /root/blacklist
COPY hg38_chrom_sizes.txt /root/hg38_chrom_sizes.txt
COPY mm10_chrom_sizes.txt /root/mm10_chrom_sizes.txt
COPY qc.py /root/qc.py
#references
#RUN wget http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz && \
#    tar -xvzf Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz && \
#    rm  Refdata_scATAC_MAESTRO_GRCh38_1.1.0.tar.gz

#RUN wget http://cistrome.org/~galib/MAESTRO/references/scATAC/Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz && \
#    tar -xvzf Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz && \
#    rm Refdata_scATAC_MAESTRO_GRCm38_1.1.0.tar.gz 


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


