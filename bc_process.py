########################################################
# Process R2 for cellranger-atac pipeline
# revised Dec 2022 for spatial ATAC v2
########################################################

import argparse
import random
from gzip import open as gzopen

from Bio.SeqIO.QualityIO import FastqGeneralIterator

print("initialzing")
seq_start = 117
bc2_start, bc2_end = 22, 30
bc1_start, bc1_end = 60, 68

ap = argparse.ArgumentParser()
ap.add_argument('-i', required=True)
ap.add_argument('-o2', required=True)
ap.add_argument('-o3', required=True)

args = vars(ap.parse_args())
print(args)

input_file_R2 = args['i']
output_file_R3 = args['o3']
output_file_R2 = args['o2']
with gzopen(input_file_R2, 'rt') as in_handle_R2, \
	   open(output_file_R3, 'w') as out_handle_R3, \
	   open(output_file_R2, 'w') as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle_R2):
        new_seq_R3 = seq[seq_start:]
        new_qual_R3 = qual[seq_start:]
        out_handle_R3.write(f'@{title}\n{new_seq_R3}\n+\n{new_qual_R3}\n'
                    )
        barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end]   # BC2 + BC1 
        new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end] 
        out_handle_R2.write(f'@{title}\n{barcode}\n+\n{new_qual_R2}\n')

print("complete")
