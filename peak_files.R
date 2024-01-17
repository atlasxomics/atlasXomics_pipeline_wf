library(Signac)
library(Seurat)
library(rhdf5)
library("ArchR")
library("GenomicRanges")
library('BSgenome')
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Rnorvegicus.UCSC.rn6)

args <- commandArgs(trailingOnly=TRUE)

fragpath <- args[1]
#narrowPeak <- args[2]
narrowPeak <- list.files('/root/Statistics/scATAC/consensus_peak_calling/MACS/'
                         , pattern = '*.narrowPeak'
                         ,full.names = TRUE
                         , include.dirs = TRUE)

species <- args[3]
run_id <- args[4]


counts = CountFragments(fragments = fragpath, verbose = FALSE)
df <- read.table(narrowPeak,sep = "\t")
names(df) <- c("chr","start","end","name","score","strand"
               ,"fold_change","neg_log10pvalue_summit","neg_log10qvalue_summit","relative_summit_position"
              )

peaks <- makeGRangesFromDataFrame(df)


rawdata <- fragpath
dataPath <- "/root/Statistics/"
file.copy(rawdata, dataPath, overwrite = TRUE)

system('cd "/root/Statistics/"; tabix -p bed fragments.tsv.gz -f')


rawdata <- "/root/Statistics/scATAC/quality_control/sample_metrics.pdf"
dataPath <- "/root/Statistics/"
file.copy(rawdata, dataPath, overwrite = TRUE)

file.rename("/root/Statistics/sample_metrics.pdf","/root/Statistics/qc_plot.pdf")

file.remove('/root/Statistics/fragments_edited.tsv.gz')


fragpath = "/root/Statistics/fragments.tsv.gz"

frags = CreateFragmentObject(path = fragpath, cells = counts$CB)

if (species=='mm'){
   genome1=seqlengths(BSgenome.Mmusculus.UCSC.mm10)
}else if (species=='hs') {
   genome1=seqlengths(BSgenome.Hsapiens.UCSC.hg38)
} else{
  genome1=seqlengths(BSgenome.Rnorvegicus.UCSC.rn6)
}


bin_matrix = GenomeBinMatrix(
  fragments = frags,
  genome = genome1, 
  binsize = 5000, 
  verbose = TRUE
)


if (species=='mm'){
   genome2="mm10"
}else if (species=='hs') {
   genome2="hg38"
} else{
  genome2="rn6"
}


chrom = CreateChromatinAssay(
  counts = bin_matrix,
  genome = genome2,
  fragments = fragpath
)


seurat_obj = CreateSeuratObject(
  counts = chrom, 
  assay = 'ATAC'
)
    

# convert peaks to count matrix

# quantify counts in each peak

macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat_obj),
  features = peaks,
    sep = c(":", "-"),
  cells = colnames(seurat_obj))


#' Write a dgCMatrix to an h5 file similar to cellRanger format
#'
#' @param mat a dgCMatrix to write.
#' @param cols_are Whether columns are "gene_names" or "sample_names". If "gene_names", mat will be transposed to match 10X conventions.
#' @param h5_target The target .h5 file for output.
#' @param ref_name Reference name for storing the data
#' @param gene_ids If available, ENSEMBL IDs for each gene
#'
#' @return an .h5 file with these mappings from dgCMatrix -> .h5:
#' \itemize{
#'   \item colnames(mat) -> /ref_name/barcodes (after transposition if cols_are == "gene_names")
#'   \item rownames(mat) -> /ref_name/gene_names (after transposition if cols_are == "gene_names")
#'   \item mat@x -> /ref_name/data
#'   \item mat@i -> /ref_name/indices
#'   \item mat@p -> /ref_name/indptr
#'   \item gene_ids -> /ref_name/gene
#' }
#'
write_dgCMatrix_h5 <- function(mat,
                               cols_are = "gene_names",
                               h5_target,
                               ref_name = "mm10-1.2.0_premrna",
                               gene_ids = NULL) {

  #library(Matrix)

  if(grepl("gene",cols_are)) {
    mat <- Matrix::t(mat)
  }

  # Create target file
  rhdf5::h5createFile(h5_target)
  # Create data group
  rhdf5::h5createGroup(h5_target,
                       ref_name)

  # Store sample ids (barcodes) and gene names
  rhdf5::h5write(colnames(mat),
                 h5_target,
                 paste0("/",ref_name,"/barcodes"))
  rhdf5::h5write(rownames(mat),
                 h5_target,
                 paste0("/",ref_name,"/gene_names"))

  if(is.null(gene_ids)) {
    gene_ids <- rownames(mat)
  }

  rhdf5::h5write(gene_ids,
                 h5_target,
                 paste0("/",ref_name,"/gene"))

  # Store dimensions as shape
  rhdf5::h5write(dim(mat),
                 h5_target,
                 paste0("/",ref_name,"/shape"))

  # Store values from mat@x as data
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/data"),
                         dims = length(mat@x),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@x,
                 h5_target,
                 paste0("/",ref_name,"/data"))

  # Store row indices from mat@i as indices
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/indices"),
                         dims = length(mat@i),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@i,
                 h5_target,
                 paste0("/",ref_name,"/indices"))

  # Store column pointers from mat@p as indptr
  rhdf5::h5write(mat@p,
                 h5_target,
                 paste0("/",ref_name,"/indptr"))

}


h5_target = paste0('/root/Statistics/', run_id,'_raw_peak_bc_matrix.h5')

write_dgCMatrix_h5(macs2_counts, cols_are = "sample_names", h5_target,ref_name = genome2, gene_ids = NULL)

file.remove('/root/Statistics/fragments.tsv.gz')
file.remove('/root/Statistics/fragments.tsv.gz.tbi')
file.remove('/root/Statistics/tmp1.txt')
file.remove('/root/Statistics/tmp2.txt')
file.remove('/root/Statistics/bc1.txt')
file.remove('/root/Statistics/bc2.txt')
file.remove('/root/Statistics/fastq_bc_inlst_freq.txt')
file.remove('/root/Statistics/chromap_bc_inlst_freq.txt')


print(r"""

    _     _    _             __  __                   _            
   / \   | |_ | |  __ _  ___ \ \/ /  ___   _ __ ___  (_)  ___  ___ 
  / _ \  | __|| | / _` |/ __| \  /  / _ \ | '_ ` _ \ | | / __|/ __|
 / ___ \ | |_ | || (_| |\__ \ /  \ | (_) || | | | | || || (__ \__ \
/_/   \_\ \__||_| \__,_||___//_/\_\ \___/ |_| |_| |_||_| \___||___/
                                                                   
                                                                                          
 
 """)
