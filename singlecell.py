import os
import glob
import gzip
import pickle
import logging
import argparse
import subprocess
import numpy as np
import requests
import pandas as pd
import pyranges as pr
from ast import literal_eval
from pathlib import Path
from typing import List, Tuple, Union




ap = argparse.ArgumentParser()
ap.add_argument('-r2', required=True)
ap.add_argument('-bc1', required=True)
ap.add_argument('-bc2', required=True)
ap.add_argument('-bed', required=True)
ap.add_argument('-w', required=True)
ap.add_argument('-tmp1', required=True)
ap.add_argument('-tmp2', required=True)
ap.add_argument('-fbif', required=True)
ap.add_argument('-cbif', required=True)
ap.add_argument('-sc', required=True)
ap.add_argument('-cis', required=True)

args = vars(ap.parse_args())

r2 = args['r2']
bc1 = args['bc1']
bc2 = args['bc2']
aln = args['bed']
whitelist = args['w']
tmp1 = args['tmp1']
tmp2 = args['tmp2']
fastq_bc_inlst_freq = args['fbif']
chromap_bc_inlst_freq = args['cbif']
singlecell = args['sc']
cistopic = args['cis']




command1 = ['zcat', r2]
command2 = ['awk "NR%4==2"']
command3 = ['cut', '-c', '61-68']
command4 = ['cut', '-c', '23-30']

# Open the output file for writing
with open(str(bc1), 'w') as a, open(str(bc2), 'w') as b:

        process1 = subprocess.Popen(command1, stdout=subprocess.PIPE)

        process2 = subprocess.Popen(
            command2, stdin=process1.stdout, stdout=subprocess.PIPE, shell=True
        )
        process1.stdout.close()

        process3 = subprocess.Popen(command3, stdin=process2.stdout, stdout=a)
        process2.stdout.close()

        process1.wait()
        process2.wait()
        process3.wait()
        process1 = subprocess.Popen(command1, stdout=subprocess.PIPE)

        process2 = subprocess.Popen(
            command2, stdin=process1.stdout, stdout=subprocess.PIPE, shell=True
        )
        process1.stdout.close()

        process4 = subprocess.Popen(command4, stdin=process2.stdout, stdout=b)
        process2.stdout.close()

        process1.wait()
        process2.wait()
        process4.wait()

# Extract genome name from local genome ref dir
#genome_id = (glob.glob(f"{species.local_path}/*.fa")[0].split("/")[-1].split("_")[0])

# Assign genome metadata
#genome_dict = {
#            "GRCh38": ["hs", "hg38_chrom_sizes.txt", "blacklist/hg38-blacklist.v2.bed"],
#            "GRCm38": ["mm", "mm10_chrom_sizes.txt", "blacklist/mm10-blacklist.v2.bed"],
#            "Rnor6": ["rnor6", "rn6_chrom_sizes.txt", None]}


# make data frame from bc files produced of fastq2 

with open(tmp1, 'w') as fw:

    subprocess.run(["paste" ,bc2, bc1], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(["awk", '{print $1, $2, $1$2}',tmp1], stdout=fw
    )

# make sure it is tab delimited

with open(tmp1, 'w') as fw:

    subprocess.run(["sed", 's/ /\t/g',tmp2], stdout=fw
    )


# subset tab delimited file so that only those in white-list remain
awk_command1 = [
    "awk",
    'NR == FNR { criteria[$1] = 1; next } $3 in criteria',
    whitelist,
    tmp1
]

try:
    with open(tmp2, "w") as b:
        subprocess.run(awk_command1, stdout=b, check=True)
    print(f"Subset rows written to {tmp2}")
except subprocess.CalledProcessError:
    print("Error occurred while executing the AWK command1")

# count frequency of each barcode


with open(tmp1, 'w') as fw:

    subprocess.run(["cut", '-f3',tmp2], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(['sort',tmp1], stdout=fw
    )

with open(tmp1, 'w') as fw:

    subprocess.run(["uniq", "-c",tmp2], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(["awk", '{print $2, $1}',tmp1], stdout=fw
    )

with open(fastq_bc_inlst_freq, 'w') as fw:

    subprocess.run(["sed", 's/ /\t/g',tmp2], stdout=fw
    )


# read chromap aln.bed file and count bc frequency

with open(tmp1, 'w') as fw:

    subprocess.run(["cut", '-f4',aln], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(['sort',tmp1], stdout=fw
    )

with open(tmp1, 'w') as fw:

    subprocess.run(["uniq", "-c",tmp2], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(["awk", '{print $2, $1}',tmp1], stdout=fw
    )

with open(chromap_bc_inlst_freq, 'w') as fw:

    subprocess.run(["sed", 's/ /\t/g',tmp2], stdout=fw
    )


subset_df = pd.read_csv(fastq_bc_inlst_freq, sep="\t", header=None)
df_bed = pd.read_csv(chromap_bc_inlst_freq, sep="\t", header=None)
cistopic_obj = pd.read_csv(cistopic)
cistopic_obj.index = cistopic_obj['barcode']
singlecell_df = pd.concat([subset_df[1], df_bed[1]], axis=1).dropna()
singlecell_df.index = subset_df[0]
singlecell_df.columns=['total','passed_filters']
singlecell_df.index= singlecell_df.index + "-1"
singlecell_df= singlecell_df[singlecell_df.index.isin(cistopic_obj.index)].reindex(cistopic_obj.index)
singlecell_df= pd.concat([singlecell_df,cistopic_obj],axis=1)
singlecell_df= singlecell_df.drop(['barcode','Unnamed: 0'], axis=1)

sum_row = singlecell_df.sum().to_frame().T
sum_row.index = ["NO_BARCODE"]
singlecell_df = pd.concat([sum_row, singlecell_df])
singlecell_df.reset_index(drop=False, inplace=True)
singlecell_df.rename(columns={'index': 'barcode'}, inplace=True)
singlecell_df['sample_id'][0] = singlecell_df['sample_id'][1]
singlecell_df.columns= singlecell_df.columns.str.replace('cisTopic_', '')
singlecell_df.iloc[0, 3:] = '-'
singlecell_df.to_csv(singlecell, index= False, header=True)


print(r"""

    _     _    _             __  __                   _            
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___ 
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/
                                                                   
                                                                                          
 
 """)
