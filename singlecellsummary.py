import argparse
import numpy as np
import pandas as pd
import subprocess

from ast import literal_eval


def chromap_log_stats(log_file, req_string):
    with open(log_file, "r") as file:
        for line in file:
            if req_string in line:
                val = literal_eval(line.strip().split(": ")[1])
                return val


ap = argparse.ArgumentParser()
ap.add_argument("-r2", required=True)
ap.add_argument("-bc1", required=True)
ap.add_argument("-bc2", required=True)
ap.add_argument("-bed", required=True)
ap.add_argument("-w", required=True)
ap.add_argument("-tmp1", required=True)
ap.add_argument("-tmp2", required=True)
ap.add_argument("-fbif", required=True)
ap.add_argument("-cbif", required=True)
ap.add_argument("-i", required=True)
ap.add_argument("-g", required=True)
ap.add_argument("-l", required=True)
ap.add_argument("-v", required=True)

args = vars(ap.parse_args())

r2 = args["r2"]
bc1 = args["bc1"]
bc2 = args["bc2"]
aln = args["bed"]
whitelist = args["w"]
tmp1 = args["tmp1"]
tmp2 = args["tmp2"]
fastq_bc_inlst_freq = args["fbif"]
chromap_bc_inlst_freq = args["cbif"]
run_id = args["i"]
genome = args["g"]
logfile = args["l"]
version = args["v"]

# Create singlecell.csv
command1 = ["zcat", r2]
command2 = ["awk 'NR%4==2'"]
command3 = ["cut", "-c", "61-68"]
command4 = ["cut", "-c", "23-30"]

# Parse the barcodes from read2 into new files (bc1.txt, bc2.txt)
with open(str(bc1), "w") as a, open(str(bc2), "w") as b:

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

# Write bc1, bc2 to tmp1.txt
with open(tmp1, "w") as fw:
    subprocess.run(["paste", bc2, bc1], stdout=fw)

# Write bc1, bc2, bc1bc2 to tmp2.txt
with open(tmp2, "w") as fw:
    subprocess.run(["awk", "{print $1, $2, $1$2}", tmp1], stdout=fw)

# Ensure tmp1.txt is tab delimited
with open(tmp1, "w") as fw:
    subprocess.run(["sed", "s/ /\t/g", tmp2], stdout=fw)

# Subset tab delimited file so that only those in whitelist remain
awk_command1 = [
    "awk",
    "NR == FNR { criteria[$1] = 1; next } $3 in criteria",
    whitelist,
    tmp1
]

try:
    with open(tmp2, "w") as b:
        subprocess.run(awk_command1, stdout=b, check=True)
    print(f"Barcodes that match whitelist written to {tmp2}")
except subprocess.CalledProcessError:
    print("Error occurred with 'awk' while filtering barcodes")

# Count the frequency of each barcode
# Extract filtered, concatenated barcodes from tmp2.txt, write to tmp1.txt
with open(tmp1, "w") as fw:
    subprocess.run(["cut", "-f3", tmp2], stdout=fw)

# Sort barcodes from tmp1.txt, write to tmp2.txt
with open(tmp2, "w") as fw:
    subprocess.run(["sort", tmp1], stdout=fw)

# Get barcode counts from tmp2.txt, write to tmp1.txt
with open(tmp1, "w") as fw:
    subprocess.run(["uniq", "-c", tmp2], stdout=fw)

# Write counts from tmp1.txt to tmp2.txt
with open(tmp2, "w") as fw:
    subprocess.run(["awk", "{print $2, $1}", tmp1], stdout=fw)

# Write tab-delimited bc counts to fastq_bc_inlst_freq.txt
with open(fastq_bc_inlst_freq, "w") as fw:
    subprocess.run(["sed", "s/ /\t/g", tmp2], stdout=fw)

# Count barcode frequency from aln.bed, write to chromap_bc_inlst_freq.txt
with open(tmp1, "w") as fw:
    subprocess.run(["cut", "-f4", aln], stdout=fw)

with open(tmp2, "w") as fw:
    subprocess.run(['sort', tmp1], stdout=fw)

with open(tmp1, "w") as fw:
    subprocess.run(["uniq", "-c", tmp2], stdout=fw)

with open(tmp2, 'w') as fw:
    subprocess.run(["awk", "{print $2, $1}", tmp1], stdout=fw)

with open(chromap_bc_inlst_freq, "w") as fw:
    subprocess.run(["sed", "s/ /\t/g", tmp2], stdout=fw)

# Open barcode counts from r2 fastq, aln.bed as DataFrames
fastq_df = pd.read_csv(fastq_bc_inlst_freq, sep="\t", header=None)
chromap_df = pd.read_csv(chromap_bc_inlst_freq, sep="\t", header=None)

# Check if fastq has barcodes that chromap didn't find, add to chromap_df
missed_barcodes = list(set(fastq_df[0].tolist()) - set(chromap_df[0].tolist()))
chromap_df = chromap_df.append(missed_barcodes, ignore_index=True)

# cleanup chromap barcodes
chromap_df.replace(np.nan, 0, inplace=True)
chromap_df = chromap_df.sort_values([0], ascending=[True])

# Check if chromap has barcodes that fastq didn't find, add to fastq_df
missed_barcodes = list(set(chromap_df[0].tolist()) - set(fastq_df[0].tolist()))
fastq_df = fastq_df.append(missed_barcodes, ignore_index=True)

# cleanup fastq barcodes
fastq_df.replace(np.nan, 0, inplace=True)
fastq_df = fastq_df.sort_values([0], ascending=[True])

# Concat fastq (total), chromap (passed filters) barcodes, clean, save
chromap_df.index = fastq_df.index
singlecell_df = pd.concat([fastq_df[1], chromap_df[1]], axis=1).dropna()
singlecell_df.index = fastq_df[0]
singlecell_df.columns = ["total", "passed_filters"]
singlecell_df.index = singlecell_df.index + "-1"
singlecell_df.to_csv(
    "/root/Statistics/singlecell.csv", header=True, index=False
)
print("Output written to singlecell.csv")

# Create summary.csv, add metadata
summary_df = pd.DataFrame(columns=["Sample ID"])
summary_df.at[0, "Sample ID"] = run_id
summary_df.at[0, "Genome"] = genome
summary_df.at[0, "Pipeline version"] = f"AtlasXomics-{version}"

# Extract summary stats from chromap log file
summary_df.at[0, "Fraction uniq-aligned reads"] = (
    chromap_log_stats(logfile, "Number of uni-mappings") /
    chromap_log_stats(logfile, "Number of mappings")
)
summary_df.at[0, "Chromap input read pairs"] = singlecell_df["total"].sum()
summary_df.at[0, "Fraction unaligned reads"] = 1 - (
    chromap_log_stats(logfile, "Number of mapped reads") /
    chromap_log_stats(logfile, "Number of reads")
)
summary_df.at[0, "Fraction reads with valid barcode"] = 1 - (
    chromap_log_stats(logfile, "Number of corrected barcodes") /
    chromap_log_stats(logfile, "Number of barcodes in whitelist")
)

# Extract number for peaks from .bed file
peak_file = pd.read_csv(
    f"./Statistics/{run_id}_peaks.bed", sep="\t", header=None
)
summary_df.at[0, "Number of peaks"] = len(peak_file.index)

# Open cistopic results csv; calculate TSS, FRIP, pct_duplicates
cistopic_obj = pd.read_csv("./Statistics/cistopic_cell_data.csv")
summary_df.at[0, "TSS_enrichment"] = max(
    [
        cistopic_obj["TSS_enrichment"].mean(),
        cistopic_obj["TSS_enrichment"].median()
    ]
)
summary_df.at[0, "FRIP"] = max(
    [
        cistopic_obj["FRIP"].mean(),
        cistopic_obj["FRIP"].median()
    ]
)
summary_df.at[0, "Fraction duplicate reads"] = max(
    [
        cistopic_obj["Dupl_rate"].mean(),
        cistopic_obj["Dupl_rate"].median()
    ]
)
summary_df.to_csv("/root/Statistics/summary.csv", header=True, index=False)
print("Output written to summary.csv")

print(r"""

    _     _    _             __  __                   _
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/


 """)
