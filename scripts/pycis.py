import argparse
import gzip
import logging
import pandas as pd
import pyranges as pr
import pybiomart as pbm

from pathlib import Path

from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.qc import compute_qc_stats, plot_sample_metrics


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


# globals
ap = argparse.ArgumentParser()
ap.add_argument("-f", required=True)
ap.add_argument("-i", required=True)
ap.add_argument("-g", required=True)
ap.add_argument("-c", required=True)
ap.add_argument("-k", required=True)
ap.add_argument("-w", required=True)
ap.add_argument("-p", required=True)
ap.add_argument("-b", required=True)
ap.add_argument("-nl", required=True)
ap.add_argument("-cp", required=True)

args = vars(ap.parse_args())

fragment = args["f"]
run_id = args["i"]
genome = args["g"]
chrsize = args["c"]
black_lst = args["k"] if genome != "rnor6" else None
whitelist = args["w"]
positions_file = args["p"]
bulk = eval(args["b"])
noLigation_bulk = eval(args["nl"])
call_peaks = eval(args["cp"])

work_dir = Path("Statistics").resolve()
tmp_dir = Path("tmp_dir").resolve()

# Filter fragment file, save path in dict for pyscistopic
logging.info("Preparing fragments file...")

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

logging.info("Preparing inputs for pycistopic...")
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

if any([bulk, noLigation_bulk, call_peaks]):

    logging.info("Performing peak calling...")

    # Prepare cell_data, chomsizes for pycistopic
    # cell_data: {index=barcode} | barcode | variable | Sample
    cell_data = pd.read_table(whitelist, header=None, index_col=None)
    cell_data.rename(columns={0: "barcode"}, inplace=True)
    cell_data.index = (
        cell_data["barcode"].astype(str) +
        "-1" +
        "___" +
        list(fragments_dict)[0]
    )
    cell_data["Sample"] = list(fragments_dict)[0]
    cell_data["barcode"] = cell_data["barcode"] + "-1"
    cell_data.index.name = None
    cell_data["variable"] = cell_data["Sample"]

    # Run pseudobulk
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data,
        variable="variable",
        sample_id_col="Sample",
        chromsizes=chromsizes,
        bed_path=f"{work_dir}/scATAC/consensus_peak_calling/pseudobulk_bed_files/",
        bigwig_path=f"{work_dir}/scATAC/consensus_peak_calling/pseudobulk_bw_files/",
        path_to_fragments=fragments_dict,
        n_cpu=77,  # sometimes increasing n_cpu > 1 results in ray memory error
        _temp_dir=f"{tmp_dir}/ray_spill",
        use_polars=False
    )

    # Run peak calling, derive consensus peaks, write to .bed
    narrow_peaks_dict = peak_calling(
        "/usr/local/bin/macs2",
        bed_paths,
        f"{work_dir}/scATAC/consensus_peak_calling/MACS/",
        genome_size=genome if genome != "rnor6" else "2.9e9",
        n_cpu=77,
        _temp_dir=f"{tmp_dir}/ray_spill"
    )
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict,
        peak_half_width=250,
        chromsizes=chromsizes,
        path_to_blacklist=black_lst
    )
    consensus_peaks.to_bed(
        path=f"{work_dir}/scATAC/consensus_peak_calling/consensus_regions.bed",
        keep=True,
        compression="infer",
        chain=False
    )
    peak_file = pd.read_table(
        f"{work_dir}/scATAC/consensus_peak_calling/MACS/{run_id}_peaks.narrowPeak",
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
    logging.info("Peaks written to Statistscs/[run_id]_peaks.bed.")

# Quality Control
logging.info("Quality control started ...")

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

stats = [
    "barcode_rank_plot",
    "duplicate_rate",
    "insert_size_distribution",
    "profile_tss"
]
if any([bulk, noLigation_bulk, call_peaks]):
    stats = stats + ["frip"]

path_to_regions = {
       run_id: f"{work_dir}/scATAC/consensus_peak_calling/consensus_regions.bed"
    } if bulk or call_peaks else None

metadata_bc, profile_data_dict = compute_qc_stats(
    fragments_dict=fragments_dict,
    tss_annotation=annot,
    stats=stats,
    label_list=None,
    path_to_regions=path_to_regions,
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
    profile_list=stats,
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
