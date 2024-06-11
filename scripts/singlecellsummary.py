import argparse
import gzip
import logging
import numpy as np
import pandas as pd

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from pathlib import Path

from ast import literal_eval


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def chromap_log_stats(log_file, req_string):
    with open(log_file, "r") as file:
        for line in file:
            if req_string in line:
                val = literal_eval(line.strip().split(": ")[1])
                return val


ap = argparse.ArgumentParser()
ap.add_argument("-r2", required=True)
ap.add_argument("-bed", required=True)
ap.add_argument("-w", required=True)
ap.add_argument("-i", required=True)
ap.add_argument("-g", required=True)
ap.add_argument("-l", required=True)
ap.add_argument("-v", required=True)

args = vars(ap.parse_args())

r2 = args["r2"]
aln = args["bed"]
whitelist = args["w"]
run_id = args["i"]
genome = args["g"]
logfile = args["l"]
version = args["v"]

tmp = Path("Statistics/tmp1.txt").resolve()

# Extract barcodes from r2 fastq, save to tmp
bc1_s = 60
bc1_e = 67

bc2_s = 22
bc2_e = 29

r2 = gzip.open(r2, "rt")

with r2, open(tmp, "w") as tmp_out:
    for _, seq, _ in FastqGeneralIterator(r2):
        tmp_out.write(f"{seq[bc2_s : bc2_e + 1]}{seq[bc1_s : bc1_e + 1]}\n")

# Read in fastq barcodes as pd.DataFrame
fastq_df = pd.read_csv(tmp, names=["barcodes"], header=None)
logging.info(f"Length fastq barcodes, pre-filter: {len(fastq_df)}")

# Filter raw fastq bcs based on whitelist
whitelist = pd.read_csv(whitelist, header=None)
fastq_df = fastq_df[fastq_df["barcodes"].isin(whitelist[0])]
logging.info(f"Length fastq barcodes, post-filter: {len(fastq_df)}")

fastq_df = fastq_df.sort_values(by=["barcodes"])
fastq_df = fastq_df.value_counts(sort=False).reset_index()
fastq_df.columns = ["barcodes", "count"]


# Count barcode frequency from aln.bed, write to chromap_bc_inlst_freq.txt
chromap_df = pd.read_csv(
    aln, header=None, names=["barcodes"], sep="\t", usecols=[3]
)
chromap_df = chromap_df.value_counts(sort=False).reset_index()
chromap_df = chromap_df.sort_values(by=["barcodes"])

# Check if fastq has barcodes that chromap didn't find, add to chromap_df
missed_barcodes = pd.Series(
    list(
        set(fastq_df["barcodes"].tolist()) -
        set(chromap_df["barcodes"].tolist())
    )
)
missed_barcodes = pd.DataFrame({"barcodes": missed_barcodes})
chromap_df = pd.concat([chromap_df, missed_barcodes], ignore_index=True)

# cleanup chromap barcodes
chromap_df.replace(np.nan, 0, inplace=True)
chromap_df = chromap_df.sort_values(["barcodes"], ascending=[True])
chromap_df.columns = ["barcodes", "count"]

# Check if chromap has barcodes that fastq didn't find, add to fastq_df
missed_barcodes = pd.Series(
    list(
        set(chromap_df["barcodes"].tolist()) -
        set(fastq_df["barcodes"].tolist())
    )
)
missed_barcodes = pd.DataFrame({"barcodes": missed_barcodes})
fastq_df = pd.concat([fastq_df, missed_barcodes], ignore_index=True)

# cleanup fastq barcodes
fastq_df.replace(np.nan, 0, inplace=True)
fastq_df = fastq_df.sort_values(["barcodes"], ascending=[True])
fastq_df.columns = ["barcodes", "count"]

# Concat fastq (total), chromap (passed filters) barcodes, clean, save
chromap_df.index = fastq_df.index
singlecell_df = pd.concat(
    [fastq_df["count"], chromap_df["count"]], axis=1
).dropna()
singlecell_df.index = fastq_df["barcodes"]
singlecell_df.columns = ["total", "passed_filters"]
singlecell_df.index = singlecell_df.index + "-1"
singlecell_df.reset_index(drop=False, inplace=True)

singlecell_df.to_csv(
    "/root/Statistics/singlecell.csv", header=True, index=False
)
logging.info("Output written to singlecell.csv")

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

# Open cistopic results csv; calculate TSS, FRIP, pct_duplicates
cistopic_obj = pd.read_csv("./Statistics/cistopic_cell_data.csv")
summary_df.at[0, "TSS_enrichment"] = max(
    [
        cistopic_obj["TSS_enrichment"].mean(),
        cistopic_obj["TSS_enrichment"].median()
    ]
)
summary_df.at[0, "Fraction duplicate reads"] = max(
    [
        cistopic_obj["Dupl_rate"].mean(),
        cistopic_obj["Dupl_rate"].median()
    ]
)
summary_df.to_csv("/root/Statistics/summary.csv", header=True, index=False)
logging.info("Output written to summary.csv")

print(r"""

    _     _    _             __  __                   _
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/


 """)
