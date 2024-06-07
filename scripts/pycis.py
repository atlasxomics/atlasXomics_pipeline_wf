import argparse
import gzip
import pandas as pd
import pyranges as pr
import pybiomart as pbm

from pathlib import Path
from pycisTopic.qc import compute_qc_stats, plot_sample_metrics

# globals
ap = argparse.ArgumentParser()
ap.add_argument("-f", required=True)
ap.add_argument("-i", required=True)
ap.add_argument("-g", required=True)
ap.add_argument("-c", required=True)
ap.add_argument("-k", required=True)
ap.add_argument("-w", required=True)
ap.add_argument("-p", required=True)

args = vars(ap.parse_args())

fragment = args["f"]
run_id = args["i"]
genome = args["g"]
chrsize = args["c"]
black_lst = args["k"] if genome != "rnor6" else None
whitelist = args["w"]
positions_file = args["p"]

work_dir = Path("Statistics").resolve()
tmp_dir = Path("tmp_dir").resolve()

# Filter fragment file, save path in dict for pyscistopic
print("preparing fragment file starting")

with gzip.open(fragment, "rb") as file:
    df = pd.read_table(file, header=None)

targeted_chromosomes = {
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
}

subset_df = df[df[0].isin(targeted_chromosomes)]
subset_df = subset_df[subset_df[1] < subset_df[2]]
fragment_edited = f"{work_dir}/fragments_edited.tsv.gz"
subset_df.to_csv(
   fragment_edited, index=False, compression="gzip", sep="\t", header=None
)
fragments_dict = {run_id: fragment_edited}

print("Preparing for pycisTopic...")
chromsizes = pd.read_csv(chrsize, sep="\t", header=None)
chromsizes.columns = ["Chromosome", "End"]
chromsizes["Start"] = [0] * chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ["Chromosome", "Start", "End"]]

# make format agree with CellRangerARC annotations
chromsizes["Chromosome"] = chromsizes["Chromosome"].str.replace("v", ".")
chromsizes["Chromosome"] = chromsizes["Chromosome"].apply(
    lambda x: x.split("_")[1] if len(x.split("_")) > 1 else x
)
chromsizes = pr.PyRanges(chromsizes)

# Quality Control
print("Quality control started ...")

# Create annot for tss calculation
if genome == "mm":
    dataset = pbm.Dataset(
       name="mmusculus_gene_ensembl",
       host="http://nov2020.archive.ensembl.org/"
    )
elif genome == "hs":
    dataset = pbm.Dataset(
        name="hsapiens_gene_ensembl", host="http://www.ensembl.org"
    )
elif genome == "rnor6":
    dataset = pbm.Dataset(
        name="rnorvegicus_gene_ensembl", host="http://www.ensembl.org"
    )
else:
    raise Exception("Incorrect genome id.")

annot = dataset.query(
   attributes=[
      "chromosome_name",
      "transcription_start_site",
      "strand",
      "external_gene_name",
      "transcript_biotype"
    ]
)
annot["Chromosome/scaffold name"] = annot[
   "Chromosome/scaffold name"
].to_numpy(dtype=str)

filter = annot["Chromosome/scaffold name"].str.contains("CHR|GL|JH|MT")
annot = annot[~filter]

annot["Chromosome/scaffold name"] = annot[
   "Chromosome/scaffold name"
].str.replace(r"(\b\S)", r"chr\1")

annot.columns = ["Chromosome", "Start", "Strand", "Gene", "Transcript_type"]
annot = annot[annot.Transcript_type == "protein_coding"]

# Compute qc stats, save as .pkl
path_to_regions = {
   run_id: f"{work_dir}/scATAC/consensus_peak_calling/consensus_regions.bed"
}

metadata_bc, profile_data_dict = compute_qc_stats(
    fragments_dict=fragments_dict,
    tss_annotation=annot,
    stats=[
       "barcode_rank_plot",
       "duplicate_rate",
       "insert_size_distribution",
       "profile_tss",
    ],
    label_list=None,
    n_cpu=77,
    valid_bc=None,
    n_frag=100,
    n_bc=None,
    tss_flank_window=1000,
    tss_minimum_signal_window=100,
    tss_rolling_window=10,
    remove_duplicates=True,
    _temp_dir=f"{tmp_dir}/ray_spill",
    use_polars=False  # True gives TypeError
)

# Create standard QC plot
plot_sample_metrics(
    profile_data_dict,
    insert_size_distribution_xlim=[0, 600],
    ncol=5,
    plot=True,
    save=f"{work_dir}/qc_plot.pdf"
)

metadata_bc[run_id].to_csv(
    f"{work_dir}/cistopic_cell_data.csv", sep=",", header=True
)

print(
    r"""

    _     _    _             __  __                   _
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/


"""
)
