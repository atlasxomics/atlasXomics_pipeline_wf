import argparse
import subprocess
import pandas as pd
import numpy as np
import pyranges as pr
from ast import literal_eval

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
ap.add_argument('-i', required=True)
ap.add_argument('-g', required=True)
ap.add_argument('-l', required=True)
ap.add_argument('-v', required=True)



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
run_id = args['i']
genome = args['g']
logfile = args['l']
version = args['v']





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

# make data frame from bc files produced of fastq2
with open(tmp1, 'w') as fw:

    subprocess.run(
        ["paste", bc2, bc1], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(
        ["awk", '{print $1, $2, $1$2}', tmp1], stdout=fw
    )

# make sure it is tab delimited

with open(tmp1, 'w') as fw:

    subprocess.run(
        ["sed", 's/ /\t/g', tmp2], stdout=fw
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

    subprocess.run(
        ["cut", '-f3', tmp2], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(
        ['sort', tmp1], stdout=fw
    )

with open(tmp1, 'w') as fw:

    subprocess.run(
        ["uniq", "-c", tmp2], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(
        ["awk", '{print $2, $1}', tmp1], stdout=fw
    )

with open(fastq_bc_inlst_freq, 'w') as fw:

    subprocess.run(
        ["sed", 's/ /\t/g', tmp2], stdout=fw
    )


# read chromap aln.bed file and count bc frequency

with open(tmp1, 'w') as fw:

    subprocess.run(
        ["cut", '-f4', aln], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(
        ['sort', tmp1], stdout=fw
    )

with open(tmp1, 'w') as fw:

    subprocess.run(
        ["uniq", "-c", tmp2], stdout=fw
    )

with open(tmp2, 'w') as fw:

    subprocess.run(
        ["awk", '{print $2, $1}', tmp1], stdout=fw
    )

with open(chromap_bc_inlst_freq, 'w') as fw:

    subprocess.run(
        ["sed", 's/ /\t/g', tmp2], stdout=fw
    )


subset_df = pd.read_csv(fastq_bc_inlst_freq, sep="\t", header=None)
df_bed = pd.read_csv(chromap_bc_inlst_freq, sep="\t", header=None)
missed_barcodes = list(set(subset_df[0].tolist()) - set(df_bed[0].tolist()))
df_bed = df_bed.append(missed_barcodes, ignore_index=True)
df_bed.replace(np.nan, 0, inplace=True)
df_bed = df_bed.sort_values([0], ascending=[True])
df_bed.index=subset_df.index
cistopic_obj = pd.read_csv("./Statistics/cistopic_cell_data.csv")
cistopic_obj.index = cistopic_obj['barcode']
singlecell_df = pd.concat([subset_df[1], df_bed[1]], axis=1).dropna()
singlecell_df.index = subset_df[0]
singlecell_df.columns = ['total', 'passed_filters']
singlecell_df.index = singlecell_df.index + "-1"
singlecell_df = singlecell_df[
    singlecell_df.index.isin(cistopic_obj.index)
    ].reindex(cistopic_obj.index)
singlecell_df = pd.concat([singlecell_df, cistopic_obj], axis=1)
singlecell_df = singlecell_df.drop(['barcode', 'Unnamed: 0'], axis=1)

sum_row = singlecell_df.sum().to_frame().T
sum_row.index = ["NO_BARCODE"]
singlecell_df = pd.concat([sum_row, singlecell_df])
singlecell_df.reset_index(drop=False, inplace=True)
singlecell_df.rename(columns={'index': 'barcode'}, inplace=True)
singlecell_df['sample_id'][0] = singlecell_df['sample_id'][1]
singlecell_df.columns = singlecell_df.columns.str.replace('cisTopic_', '')
singlecell_df.iloc[0, 3:] = '-'
singlecell_df.to_csv("/root/Statistics/singlecell.csv", header=True, index= False)
print("Output written to singlecell.csv")





def chromap_log_stats(log_file, req_string):
    with open(log_file, 'r') as file:
        for line in file:
            if req_string in line:
                val = literal_eval(line.strip().split(": ")[1])
                return val
                break



summary_df = pd.DataFrame(columns=['Sample ID'])
summary_df.at[0, 'Sample ID'] = f"{run_id}"
summary_df.at[0, 'Genome'] = genome
summary_df.at[0, 'Pipeline version'] = "AtlasXomics-" + version
summary_df.at[0, 'Fraction uniq-aligned reads'] = chromap_log_stats(logfile,"Number of uni-mappings") / chromap_log_stats(logfile,"Number of mappings")
summary_df.at[0, 'Chromap input read pairs'] = singlecell_df['total'][0]
summary_df.at[0, 'Fraction unaligned reads'] = 1-(chromap_log_stats(logfile,"Number of mapped reads") / chromap_log_stats(logfile,"Number of reads"))
summary_df.at[0, 'Fraction reads with valid barcode'] = 1-(chromap_log_stats(logfile,"Number of corrected barcodes") / chromap_log_stats(logfile,"Number of barcodes in whitelist"))

peak_file = pd.read_csv(f"./Statistics/{run_id}_peaks.bed", sep='\t', header=None)
summary_df.at[0, 'Number of peaks'] = len(peak_file.index)
summary_df.at[0, 'TSS_enrichment'] = max([cistopic_obj['TSS_enrichment'].mean(),cistopic_obj['TSS_enrichment'].median()])
summary_df.at[0, 'FRIP'] = max([cistopic_obj['FRIP'].mean(),cistopic_obj['FRIP'].median()])
summary_df.at[0, 'Fraction duplicate reads'] = max([cistopic_obj['Dupl_rate'].mean(),cistopic_obj['Dupl_rate'].median()])
summary_df.to_csv("/root/Statistics/summary.csv", header=True, index= False)
print("Output written to summary.csv")



print(r"""

    _     _    _             __  __                   _            
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___ 
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/



 """)
