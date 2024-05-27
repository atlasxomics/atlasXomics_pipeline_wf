import argparse
import gzip
import os
import pandas as pd
import pathlib
import pickle
import pyranges as pr
import pybiomart as pbm
import sys

from pycisTopic.qc import *
from pycisTopic.cistopic_class import *
from pycisTopic.pseudobulk_peak_calling import *
from pycisTopic.iterative_peak_calling import *


ap = argparse.ArgumentParser()
ap.add_argument("-f", required=True)
ap.add_argument("-i", required=True)
ap.add_argument("-g", required=True)
ap.add_argument("-c", required=True)
ap.add_argument("-k", required=True)
ap.add_argument("-w", required=True)
ap.add_argument("-d", required=True)
ap.add_argument("-t", required=True)
ap.add_argument("-p", required=True)

args = vars(ap.parse_args())

fragment = args["f"]
run_id = args["i"]
genome = args["g"]
chrsize = args["c"]
black_lst = args["k"] if genome != "rnor6" else None
whitelist = args["w"]
work_dir = args["d"]
tmp_dir = args["t"]
positions_file = args["p"]

#########################################
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

#################################pycisTopic###############################
print("pycisTopic Starting ...")

fragments_dict = {run_id: fragment_edited}

cell_data = pd.read_table(whitelist, header=None, index_col=None)
cell_data.rename(columns={0: "barcode"}, inplace=True)
cell_data.index = (
    cell_data["barcode"].astype(str) + "-1" + "___" + list(fragments_dict)[0]
)
cell_data["Sample"] = list(fragments_dict)[0]
cell_data["barcode"] = cell_data["barcode"] + "-1"
cell_data.index.name = None
cell_data["variable"] = cell_data["Sample"]

chromsizes = pd.read_csv(chrsize, sep="\t", header=None)
chromsizes.columns = ["Chromosome", "End"]
chromsizes["Start"] = [0] * chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ["Chromosome", "Start", "End"]]

# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes["Chromosome"] = [
    chromsizes["Chromosome"][x].replace("v", ".")
    for x in range(len(chromsizes["Chromosome"]))
]
chromsizes["Chromosome"] = [
    chromsizes["Chromosome"][x].split("_")[1]
    if len(chromsizes["Chromosome"][x].split("_")) > 1
    else chromsizes["Chromosome"][x]
    for x in range(len(chromsizes["Chromosome"]))
]
chromsizes = pr.PyRanges(chromsizes)

bw_paths, bed_paths = export_pseudobulk(
    input_data=cell_data,
    variable="variable",
    sample_id_col="Sample",
    chromsizes=chromsizes,
    bed_path=work_dir + "/scATAC/consensus_peak_calling/pseudobulk_bed_files/",
    bigwig_path=work_dir + "/scATAC/consensus_peak_calling/pseudobulk_bw_files/",
    path_to_fragments=fragments_dict,
    n_cpu=77,  # sometimes increasing n_cpu > 1 results in  ray memory error
    normalize_bigwig=True,
    remove_duplicates=True,
    _temp_dir=tmp_dir + "ray_spill",  # tmp dir must not be too long
    split_pattern="___",
    use_polars=False
)

with open(
   f"{work_dir}/scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl",
   "wb"
) as f:
    pickle.dump(bed_paths, f)

with open(
   f"{work_dir}/scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl",
   "wb"
) as f:
    pickle.dump(bw_paths, f)

infile = open(
   f"{work_dir}/scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl",
   "rb"
)
bed_paths = pickle.load(infile)
infile.close()

macs_path = "/usr/local/bin/macs2"
outdir = work_dir + "/scATAC/consensus_peak_calling/MACS/"

genome_size = genome if genome != "rnor6" else "2.9e9"

# Run peak calling
narrow_peaks_dict = peak_calling(
    macs_path,
    bed_paths,
    outdir,
    genome_size=genome_size,
    n_cpu=77,
    input_format="BEDPE",
    shift=73,
    ext_size=146,
    keep_dup="all",
    q_value=0.05,
    _temp_dir=tmp_dir + "ray_spill"
)

peak_half_width = 250
consensus_peaks = get_consensus_peaks( # iterative_peak_calling
   narrow_peaks_dict,
   peak_half_width,
   chromsizes=chromsizes,
   path_to_blacklist=black_lst
)

# Write to bed
consensus_peaks.to_bed(
   path=work_dir + "/scATAC/consensus_peak_calling/consensus_regions.bed",
   keep=True,
   compression="infer",
   chain=False
)

print("Quality control started ...")
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
       "frip"
    ],
    label_list=None,
    path_to_regions=path_to_regions,
    n_cpu=77,
    valid_bc=None,
    n_frag=100,
    n_bc=None,
    tss_flank_window=1000,
    tss_window=50,
    tss_minimum_signal_window=100,
    tss_rolling_window=10,
    remove_duplicates=True,
    _temp_dir=os.path.join(tmp_dir + 'ray_spill'),
    use_polars=False  # True gives TypeError: __init__() got an unexpected
)                     # keyword argument 'encoding'.

if not os.path.exists(f"{work_dir}/scATAC/quality_control"):
    os.makedirs(f"{work_dir}/scATAC/quality_control")

pickle.dump(
    metadata_bc,
    open(f"{work_dir}/scATAC/quality_control/metadata_bc.pkl", "wb")
)

pickle.dump(
    profile_data_dict,
    open(f"{work_dir}/scATAC/quality_control/profile_data_dict.pkl", "wb")
)

with open(f"{work_dir}/scATAC/quality_control/metadata_bc.pkl", "wb") as f:
    pickle.dump(metadata_bc, f)

with open(
    f"{work_dir}/scATAC/quality_control/profile_data_dict.pkl", "wb"
) as f:
    pickle.dump(profile_data_dict, f)

plot_sample_metrics(
    profile_data_dict,
    insert_size_distribution_xlim=[0, 600],
    ncol=5,
    plot=True,
    save=f"{work_dir}/scATAC/quality_control/sample_metrics.pdf")

qc_filters = {
    run_id: {
        "Log_unique_nr_frag": [0, None],
        "FRIP": [0.0, None],
        "TSS_enrichment": [0, None],
        "Dupl_rate": [None, None]
    }
}

# Return figure to plot together with other metrics, and cells passing filters.
# Figure will be saved as pdf.
FRIP_NR_FRAG_filterDict = {}
TSS_NR_FRAG_filterDict = {}
FRIP_NR_FRAG_figDict = {}
TSS_NR_FRAG_figDict = {}
DR_NR_FRAG_figDict = {}

for runID in metadata_bc:
    FRIP_NR_FRAG_fig, FRIP_NR_FRAG_filter = plot_barcode_metrics(
        metadata_bc[runID],
        var_x="Log_unique_nr_frag",
        var_y="FRIP",
        min_x=qc_filters[runID]["Log_unique_nr_frag"][0],
        max_x=qc_filters[runID]["Log_unique_nr_frag"][1],
        min_y=qc_filters[runID]["FRIP"][0],
        max_y=qc_filters[runID]["FRIP"][1],
        return_cells=True,
        return_fig=True,
        plot=False
    )

    # Return figure plot together with other metrics, cells passing filters.
    TSS_NR_FRAG_fig, TSS_NR_FRAG_filter = plot_barcode_metrics(
        metadata_bc[runID],
        var_x='Log_unique_nr_frag',
        var_y='TSS_enrichment',
        min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
        max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
        min_y=qc_filters[runID]['TSS_enrichment'][0],
        max_y=qc_filters[runID]['TSS_enrichment'][1],
        return_cells=True,
        return_fig=True,
        plot=False
    )

    # Return figure to plot together with other metrics, but not returning
    # cells (no filter applied for the duplication rate  per barcode)
    DR_NR_FRAG_fig = plot_barcode_metrics(
        metadata_bc[runID],
        var_x="Log_unique_nr_frag",
        var_y="Dupl_rate",
        min_x=qc_filters[runID]["Log_unique_nr_frag"][0],
        max_x=qc_filters[runID]["Log_unique_nr_frag"][1],
        min_y=qc_filters[runID]["Dupl_rate"][0],
        max_y=qc_filters[runID]["Dupl_rate"][1],
        return_cells=False,
        return_fig=True,
        plot=False,
        plot_as_hexbin=True
    )

    # Barcodes passing filters
    FRIP_NR_FRAG_filterDict[runID] = FRIP_NR_FRAG_filter
    TSS_NR_FRAG_filterDict[runID] = TSS_NR_FRAG_filter

    # Figs
    FRIP_NR_FRAG_figDict[runID] = FRIP_NR_FRAG_fig
    TSS_NR_FRAG_figDict[runID] = TSS_NR_FRAG_fig
    DR_NR_FRAG_figDict[runID] = DR_NR_FRAG_fig

bc_passing_filters = dict()
for runID in metadata_bc:
    bc_passing_filters[runID] = list(
        (set(FRIP_NR_FRAG_filterDict[runID]) &
         set(TSS_NR_FRAG_filterDict[runID]))
    )
    bc_len = len(bc_passing_filters[runID])
    print(
        f"{bc_len} barcodes passed filters for sample {runID}"
    )

with open(
    f"{work_dir}/scATAC/quality_control/bc_passing_filters.pkl", "wb"
) as f:
    pickle.dump(bc_passing_filters, f)

metadata_bc = pickle.load(
    open(
        os.path.join(f"{work_dir}/scATAC/quality_control/metadata_bc.pkl"),
        "rb"
    )
)
# n.b. we use the same regions twice
path_to_regions = {
    run_id: os.path.join(
        f"{work_dir}/scATAC/consensus_peak_calling/consensus_regions.bed"
    )
}
metadata_bc = pickle.load(
    open(
        os.path.join(f"{work_dir}/scATAC/quality_control/metadata_bc.pkl"),
        "rb"
    )
)
bc_passing_filters = pickle.load(
    open(
        os.path.join(
            work_dir + "/scATAC/quality_control/bc_passing_filters.pkl"
        ),
        "rb"
    )
)
cistopic_obj_list = [
    create_cistopic_object_from_fragments( # cistopic_class
        path_to_fragments=fragments_dict[key],
        path_to_regions=path_to_regions[key],
        path_to_blacklist=black_lst,
        metrics=metadata_bc[key],
        valid_bc=bc_passing_filters[key],
        n_cpu=77,
        use_polars=False,  # True gives TypeError
        project=key
    )
    for key in fragments_dict.keys()
]

cistopic_obj = merge(cistopic_obj_list)

cistopic_obj.cell_data.to_csv(
    f"{work_dir}/cistopic_cell_data.csv", sep=",", header=True
)

peak_file = pd.read_table(
    str(
        list(
            pathlib.Path(
                os.path.join(
                    work_dir + "/scATAC/consensus_peak_calling/MACS/"
                )
            ).glob("*.narrowPeak")
        )[0]
    ),
    sep="\t",
    header=None
)

peak_file.columns = [
    "chr",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "fold_change",
    "neg_log10pvalue_summit",
    "neg_log10qvalue_summit",
    "relative_summit_position"
]

peak_file = peak_file.drop_duplicates(subset="start", keep="last")

peak_file.to_csv(
    f"{work_dir}/{run_id}_peaks.bed", sep="\t", header=False, index=False
)

print("peaks written to _peaks.bed")

print(
    r"""

    _     _    _             __  __                   _
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/


"""
)
