FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_total_wf:16.0.7-549da6


# Copy barcode files to root dir
COPY bc50.txt /root/bc50.txt
COPY bc50_old.txt /root/bc50_old.txt
COPY bc96.txt /root/bc96.txt
COPY bc96_fg.txt /root/bc96_fg.txt
COPY bcFG210v4.txt /root/bcFG210v4.txt

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
RUN echo 'hello'
RUN pip install marshmallow-enum
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root


