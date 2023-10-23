import os
import gzip
import pickle
import logging
import argparse
import matplotlib.pyplot as plt
import numpy as np
import requests
import pandas as pd
import pyranges as pr
import pybiomart as pbm
from pycisTopic.qc import *
from ast import literal_eval
from pycisTopic.cistopic_class import *
from pycisTopic.pseudobulk_peak_calling import *
from pycisTopic.pseudobulk_peak_calling import *
from pycisTopic.iterative_peak_calling import *
from pathlib import Path
from typing import List, Tuple, Union


ap = argparse.ArgumentParser()
ap.add_argument('-bc1', required=True)
ap.add_argument('-bc2', required=True)
ap.add_argument('-bcff', required=True)
ap.add_argument('-bed', required=True)
ap.add_argument('-f', required=True)
ap.add_argument('-bedbf', required=True)
ap.add_argument('-i', required=True)
ap.add_argument('-g', required=True)
ap.add_argument('-c', required=True)
ap.add_argument('-k', required=True)
ap.add_argument('-w', required=True)
ap.add_argument('-l', required=True)
ap.add_argument('-d', required=True)
ap.add_argument('-t', required=True)
#ap.add_argument('-otd', required=True)
ap.add_argument('-p', required=True)
ap.add_argument('-v', required=True)
 
args = vars(ap.parse_args())

bc1 = args['bc1']
bc2 = args['bc2']
bc_freq_filtered = args['bcff']
bed = args['bed']
fragment = args['f']
bed_bc_freq = args['bedbf']
run_id = args['i']
genome = args['g']
chrsize = args['c']
black_lst = args['k']
whitelist = args['w']
logfile = args['l']
work_dir = args['d']
tmp_dir = args['t']
#out_dir = args['otd']
positions_file = args['p']
version = args['v']



#########################################
print('preparing fragment file starting')

with gzip.open(fragment, 'rb') as file:
    df = pd.read_table(file, header=None)

targeted_chromosomes = {"chr1", "chr2","chr3", "chr4","chr5", "chr6", "chr7","chr8", "chr9","chr10",
                       "chr11", "chr12","chr13", "chr14","chr15", "chr16", "chr17","chr18", "chr19","chr20","chr21","chrX", "chrY"}

subset_df = df[df[0].isin(targeted_chromosomes)]


subset_df = subset_df[subset_df[1] < subset_df[2]]


fragment_edited = work_dir + "/fragments_edited.tsv.gz"

subset_df.to_csv(fragment_edited, index=False, compression='gzip', sep = "\t", header=None)

print(f"Output written to '{work_dir} + singlecell.csv'")
###########################################################pycisTopic#############################################
print('pycisTopic Starting ...')

fragments_dict = {run_id:fragment_edited}


cell_data = pd.read_table(whitelist, header=None, index_col=None)
cell_data.rename(columns={0: 'barcode'}, inplace=True)
cell_data.index=cell_data['barcode'].astype(str) + '-1' + '___' + list(fragments_dict)[0]
cell_data['Sample'] = list(fragments_dict)[0]
cell_data['barcode']=cell_data['barcode'] + '-1'
cell_data.index.name = None
cell_data['variable'] = cell_data['Sample']

chromsizes=pd.read_csv(chrsize, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)


bw_paths, bed_paths = export_pseudobulk(
    input_data=cell_data,
    variable='variable',
    sample_id_col='Sample',
    chromsizes=chromsizes,
    bed_path=work_dir + '/scATAC/consensus_peak_calling/pseudobulk_bed_files/',
    bigwig_path=work_dir + '/scATAC/consensus_peak_calling/pseudobulk_bw_files/',
    path_to_fragments=fragments_dict,
    n_cpu=77,  # sometimes increasing n_cpu > 1 doesn't produce bed files. ray memory Error happens
    normalize_bigwig=True,
    remove_duplicates=True,
    _temp_dir=tmp_dir + 'ray_spill',  # tmp dir must not be too long
    split_pattern='___',
    use_polars=False
) 

with open(work_dir + '/scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl', 'wb') as f:
  pickle.dump(bed_paths, f)

with open(work_dir + '/scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl', 'wb') as f:
  pickle.dump(bw_paths, f)



infile = open(work_dir + '/scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl', 'rb')
bed_paths = pickle.load(infile)
infile.close()


macs_path = '/usr/local/bin/macs2'

outdir = work_dir + '/scATAC/consensus_peak_calling/MACS/'
# Run peak calling
narrow_peaks_dict = peak_calling(
    macs_path,
    bed_paths,
    outdir,
    genome_size=genome,
    n_cpu=77,
    input_format='BEDPE',
    shift=73,
    ext_size=146,
    keep_dup='all',
    q_value=0.05,
    _temp_dir=tmp_dir + 'ray_spill'
)

peak_half_width = 250
consensus_peaks = get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes, path_to_blacklist=black_lst)


# Write to bed
consensus_peaks.to_bed(path= work_dir + '/scATAC/consensus_peak_calling/consensus_regions.bed'
                       , keep=True
                       , compression='infer'
                       , chain=False)

print('Quality control started ...')


if (genome == 'mm'):
   dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host='http://nov2020.archive.ensembl.org/')
else:
   dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')

annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
annot = annot[~filter]
annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']

path_to_regions = {run_id: work_dir + '/scATAC/consensus_peak_calling/consensus_regions.bed'}

sys.stderr = open(os.devnull, "w")  # silence stderr
metadata_bc, profile_data_dict = compute_qc_stats(
    fragments_dict=fragments_dict,
    tss_annotation=annot,
    stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution',
           'profile_tss', 'frip'],
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
    use_polars=False  # True gives TypeError: __init__() got an unexpected keyword argument 'encoding'
)


sys.stderr = sys.__stderr__  # unsilence stderr

if not os.path.exists(work_dir + '/scATAC/quality_control'):
    os.makedirs(work_dir + '/scATAC/quality_control')

pickle.dump(metadata_bc,
            open(work_dir + '/scATAC/quality_control/metadata_bc.pkl', 'wb'))

pickle.dump(profile_data_dict,
            open(work_dir + '/scATAC/quality_control/profile_data_dict.pkl', 'wb'))



with open(work_dir + '/scATAC/quality_control/metadata_bc.pkl', 'wb') as f:
  pickle.dump(metadata_bc, f)


with open(work_dir + '/scATAC/quality_control/profile_data_dict.pkl', 'wb') as f:
  pickle.dump(profile_data_dict, f)


plot_sample_metrics(profile_data_dict,
           insert_size_distribution_xlim=[0,600],
           ncol=5,
           plot=True,
           save= work_dir + '/scATAC/quality_control/sample_metrics.pdf')

qc_filters = {
    run_id: {
        'Log_unique_nr_frag':   [0, None],
        'FRIP':                 [0.0, None],
        'TSS_enrichment':       [0, None],
        'Dupl_rate':            [None, None]
    }
}


# Return figure to plot together with other metrics, and cells passing filters. Figure will be saved as pdf.
FRIP_NR_FRAG_filterDict = {}
TSS_NR_FRAG_filterDict = {}
FRIP_NR_FRAG_figDict = {}
TSS_NR_FRAG_figDict = {}
DR_NR_FRAG_figDict={}

from pycisTopic.qc import *
for runID in metadata_bc:
    FRIP_NR_FRAG_fig, FRIP_NR_FRAG_filter=plot_barcode_metrics(metadata_bc[runID],
                                           var_x='Log_unique_nr_frag',
                                           var_y='FRIP',
                                           min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
                                           max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
                                           min_y=qc_filters[runID]['FRIP'][0],
                                           max_y=qc_filters[runID]['FRIP'][1],
                                           return_cells=True,
                                           return_fig=True,
                                           plot=False)
    # Return figure to plot together with other metrics, and cells passing filters
    TSS_NR_FRAG_fig, TSS_NR_FRAG_filter=plot_barcode_metrics(metadata_bc[runID],
                                          var_x='Log_unique_nr_frag',
                                          var_y='TSS_enrichment',
                                          min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
                                          max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
                                          min_y=qc_filters[runID]['TSS_enrichment'][0],
                                          max_y=qc_filters[runID]['TSS_enrichment'][1],
                                          return_cells=True,
                                          return_fig=True,
                                          plot=False)
    # Return figure to plot together with other metrics, but not returning cells (no filter applied for the duplication rate  per barcode)
    DR_NR_FRAG_fig=plot_barcode_metrics(metadata_bc[runID],
                                          var_x='Log_unique_nr_frag',
                                          var_y='Dupl_rate',
                                          min_x=qc_filters[runID]['Log_unique_nr_frag'][0],
                                          max_x=qc_filters[runID]['Log_unique_nr_frag'][1],
                                          min_y=qc_filters[runID]['Dupl_rate'][0],
                                          max_y=qc_filters[runID]['Dupl_rate'][1],
                                          return_cells=False,
                                          return_fig=True,
                                          plot=False,
                                          plot_as_hexbin = True)

    # Barcodes passing filters
    FRIP_NR_FRAG_filterDict[runID] = FRIP_NR_FRAG_filter
    TSS_NR_FRAG_filterDict[runID] = TSS_NR_FRAG_filter
    # Figs
    FRIP_NR_FRAG_figDict[runID] = FRIP_NR_FRAG_fig
    TSS_NR_FRAG_figDict[runID] = TSS_NR_FRAG_fig
    DR_NR_FRAG_figDict[runID]=DR_NR_FRAG_fig


bc_passing_filters = dict()
for runID in metadata_bc:
    bc_passing_filters[runID] = list((set(FRIP_NR_FRAG_filterDict[runID]) & set(TSS_NR_FRAG_filterDict[runID])))
    print(f"{len(bc_passing_filters[runID])} barcodes passed filters for sample {runID}")


with open(work_dir +'/scATAC/quality_control/bc_passing_filters.pkl', 'wb') as f:
  pickle.dump(bc_passing_filters, f)


metadata_bc = pickle.load(open(os.path.join(work_dir + '/scATAC/quality_control/metadata_bc.pkl'), 'rb'))


# note we use twice the same regions!
path_to_regions = {
    run_id: os.path.join(work_dir + '/scATAC/consensus_peak_calling/consensus_regions.bed')
}

path_to_blacklist = black_lst

metadata_bc = pickle.load(open(os.path.join(work_dir + '/scATAC/quality_control/metadata_bc.pkl'), 'rb'))
bc_passing_filters = pickle.load(open(os.path.join(work_dir + '/scATAC/quality_control/bc_passing_filters.pkl'), 'rb'))

cistopic_obj_list = [
    create_cistopic_object_from_fragments(
        path_to_fragments=fragments_dict[key],
        path_to_regions=path_to_regions[key],
        path_to_blacklist=black_lst,
        metrics=metadata_bc[key],
        valid_bc=bc_passing_filters[key],
        n_cpu=77,
        use_polars=False,  # True gives TypeError: __init__() got an unexpected keyword argument 'encoding'
        project=key)
    for key in fragments_dict.keys()
]

cistopic_obj = merge(cistopic_obj_list)

cistopic_obj.cell_data.to_csv(work_dir + '/cistopic_cell_data.csv', sep=",", header=True)


df = pd.read_csv(bc2, names=['bc2'])
df['bc1'] = pd.read_csv(bc1, names=['bc1'])
df['bc'] = df['bc2'] + df['bc1']
with open(whitelist, 'r') as file:
     subset_values = [line.strip() for line in file]
subset_df = df[df['bc'].isin(subset_values)]
subset_df['bc'].value_counts().to_csv(bc_freq_filtered, sep='\t', header=False)
df_bed = pd.read_table(bed, sep="\t", header=None)
df_bed.columns = ['chr', 'start', 'end', 'bc','freq']
df_bed['bc'].value_counts().to_csv(bed_bc_freq, sep='\t', header=False)

singlecell_df = pd.concat([subset_df['bc'].value_counts(), df_bed['bc'].value_counts()], axis=1).dropna()
singlecell_df.columns=['total','passed_filters']
singlecell_df.index= singlecell_df.index + "-1"
cistopic_obj.cell_data.index = cistopic_obj.cell_data['barcode']
singlecell_df= singlecell_df[singlecell_df.index.isin(cistopic_obj.cell_data.index)].reindex(cistopic_obj.cell_data.index)
singlecell_df= pd.concat([singlecell_df,cistopic_obj.cell_data],axis=1)
singlecell_df= singlecell_df.drop(['barcode'], axis=1)


sum_row = singlecell_df.sum().to_frame().T
sum_row.index = ["NO_BARCODE"]
singlecell_df = pd.concat([sum_row, singlecell_df])
singlecell_df.reset_index(drop=False, inplace=True)
singlecell_df.rename(columns={'index': 'barcode'}, inplace=True)
singlecell_df['sample_id'][0] = singlecell_df['sample_id'][1]
singlecell_df.columns= singlecell_df.columns.str.replace('cisTopic_', '')
singlecell_df.iloc[0, 3:] = '-'
singlecell_df.to_csv("/root/Statistics/singlecell.csv", index= False, header=True)
print("Output written to singlecell.csv")

peak_file = pd.read_table(open(os.path.join(work_dir + f'/scATAC/consensus_peak_calling/MACS/{run_id}_peaks.narrowPeak'), 'rb'), sep='\t', header=None)
peak_file.columns=["chr","start","end","name","score","strand"
                                 ,"fold_change","neg_log10pvalue_summit","neg_log10qvalue_summit","relative_summit_position"]

peak_file = peak_file.drop_duplicates(subset='start', keep="last")

peak_file.to_csv(f"{work_dir}/{run_id}_peaks.bed", sep='\t', header=False, index=False)

print("peaks written to _peaks.bed")

def chromap_log_stats(file_path, target_string):
    # Open the file for reading
    with open(file_path, 'r') as file:
        for line in file:
            # Check if the target string is present in the line
            if target_string in line:
                # If the string is found, print the line and break out of the loop
                val = literal_eval(line.strip().split(": ")[1])
                return val
                break



summary_df = pd.DataFrame(columns=['Sample ID'])
summary_df.at[0, 'Sample ID'] = f"{run_id}"
summary_df.at[0, 'Genome'] = genome
summary_df.at[0, 'Pipeline version'] = "AtlasXomics-" + version
summary_df.at[0, 'Fraction aligned reads'] = chromap_log_stats(logfile,"Number of uni-mappings") / chromap_log_stats(logfile,"Number of mappings")
summary_df.at[0, 'Chromap input read pairs'] = singlecell_df['total'][0]
summary_df.at[0, 'Fraction unaligned reads'] = 1-(chromap_log_stats(logfile,"Number of mapped reads") / chromap_log_stats(logfile,"Number of reads"))
summary_df.at[0, 'Fraction reads with valid barcode'] = 1-(chromap_log_stats(logfile,"Number of corrected barcodes") / chromap_log_stats(logfile,"Number of barcodes in whitelist"))
summary_df.at[0, 'Number of peaks'] = len(peak_file.index)
summary_df.at[0, 'TSS_enrichment'] = max([cistopic_obj.cell_data['TSS_enrichment'].mean(),cistopic_obj.cell_data['TSS_enrichment'].median()])
summary_df.at[0, 'FRIP'] = max([cistopic_obj.cell_data['FRIP'].mean(),cistopic_obj.cell_data['FRIP'].median()])
summary_df.at[0, 'Fraction duplicate reads'] = max([cistopic_obj.cell_data['Dupl_rate'].mean(),cistopic_obj.cell_data['Dupl_rate'].median()])
summary_df.to_csv("/root/Statistics/summary.csv", header=True, index= False)
print(f"Output written to summary.csv")

def get_axis_avgs(
    singlecell_path: Path,
    positions_path: Path
) -> pd.DataFrame:
    """
    Given singlecell.csv and barcode file, return pd.DataFrame with median
    row and column counts for each lane.
    """
    
    sc = pd.read_csv(singlecell_path)
    sc = sc[sc["barcode"] != "NO_BARCODE"]
    sc["barcode"] = sc["barcode"].apply(lambda x: x.strip("-1"))

    positions = pd.read_csv(positions_path, header=None)
    positions.columns = ["barcode", "on_tissue", "row", "column"]

    merged = pd.merge(sc, positions)

    avgs = [merged.groupby([axis]).median(numeric_only=True)["passed_filters"]
            for axis in ["row", "column"]]

    df = pd.merge(avgs[0], avgs[1], right_index=True, left_index=True)
    df.columns = ["row_avg", "col_avg"]
    df.index.name = None
    
    return df

def get_upper_bounds(
    avgs: Union[pd.Series, List[float]],
    sigma: int=1
) -> Tuple[float, float]:    
    """
    Given pd.DataFrame with median row and column counts for each lane,
    return values above with a lane average is considered an outlier for
    """
        
    mean = np.mean(avgs)
    std = np.std(avgs)

    row, col = mean + sigma * std
    
    return row, col

def get_outliers(
    avgs: Union[pd.Series, List[float]],
    row_bound: float,
    col_bound: float
) -> pd.DataFrame:
    """Add boolean column identifying outliers to lane averages;
    merger with average to return a dataframe."""
    
    row_outliers = avgs["row_avg"] > row_bound
    col_outliers = avgs["col_avg"] > col_bound
    
    avgs = avgs.merge(row_outliers, left_index=True, right_index=True)
    avgs = avgs.merge(col_outliers, left_index=True, right_index=True)
    
    avgs.columns = ["row_avg", "col_avg", "row_outlier", "col_outlier"]
    
    return avgs
def plotting_task(
    singlecell_path: Path,
    positions_path: Path
):
    """Save a barplot.pdf containing row/column medians with outlier
    highlighted at sigma 1 and 2"""

    dfs = []
    for sigma in [1,2]:
        avgs = get_axis_avgs(singlecell_path, positions_path)
        row_bound, col_bound = get_upper_bounds(avgs, sigma=sigma)
        df = get_outliers(avgs, row_bound, col_bound)
        dfs.append((sigma, df))

    plt.rc("figure", figsize=(15,10))

    fig, ax = plt.subplots(4, 1)
    plt.subplots_adjust(wspace=0, hspace=0.7)

    i = 0
    for sigma, df in dfs:

        row_x = df.index
        row_y = df["row_avg"]
        row_colors = df["row_outlier"]

        ax[i].bar(    
            [x + 1 for x in row_x[~row_colors]],
            row_y[~row_colors],
            color="blue",
            edgecolor=(0,0,0),
        )

        ax[i].bar(    
            [x + 1 for x in row_x[row_colors]],
            row_y[row_colors],
            color="red",
            edgecolor=(0,0,0),
            label = "outliers"
        )

        ax[i].set_title(f"row medians (sigma={sigma})")
        ax[i].set_ylabel("median frag counts")
        ax[i].legend()
        
        i += 1

        col_x = df.index
        col_y = df["col_avg"]
        col_colors = df["col_outlier"]

        ax[i].bar(    
            [x + 1 for x in col_x[~col_colors]],
            col_y[~col_colors],
            color="blue",
            edgecolor=(0,0,0),
        )

        ax[i].bar(    
            [x + 1 for x in col_x[col_colors]],
            col_y[col_colors],
            color="red",
            edgecolor=(0,0,0),
            label = "outlier"
        )

        ax[i].set_title(f"column medians (sigma={sigma})")
        ax[i].set_ylabel("median frag counts")
        ax[i].legend()

        xticks = [x + 1 for x in range(len(df.index))]
        font_size = 5 if len(df.index) == 96 else 9
        ax[i].set_xticks(xticks)
        ax[i].set_xticklabels(xticks, fontsize=font_size)

        i += 1
        
    plt.savefig("/root/Statistics/lane_qc.pdf")

plotting_task("/root/Statistics/singlecell.csv", positions_file)

print("guess what! statistcs calculations successfully finished!!")


print(r"""

    _     _    _             __  __                   _            
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___ 
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/
                                                                   
                                                                                          
 
 """)
