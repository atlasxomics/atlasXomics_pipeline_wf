########################################################
# Process R2 for cellranger-atac or chromap pipeline
# revised Dec 2023 
########################################################

import argparse
import random
import logging
from gzip import open as gzopen

from Bio.SeqIO.QualityIO import FastqGeneralIterator

## define some functions and setup ##
logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)

#print("initialzing")
#seq_start = 117
#bc2_start, bc2_end = 22, 30
#bc1_start, bc1_end = 60, 68


def makeLeader(fullBC, merLen=8, fullLeader=True):
    fullL = f"CAAGCGTTGGCTTCTCGCATCT{fullBC[:(merLen)]}ATCCACGTGCTTGAGAGGCCAGAGCATTCG{fullBC[(merLen):]}GTGGCCGATGTTTCGCATCGGCGTACGACTAGATGTGTATAAGAGACAG"
    noLigL = f"CAAGCGTTGGCTTCTCGCATCT{fullBC[:(merLen)]}ATCCACGTGCTTGAGAGGCCAGAGCATTCG{fullBC[(merLen):]}GTGGCCGATGTTTCG"
    return fullL if fullLeader else noLigL

### end functions ###

logging.info("Initializing Barcode Read Processing")

ap = argparse.ArgumentParser()
ap.add_argument('-i', required=True)
ap.add_argument('-o2', required=True)
ap.add_argument('-o3', required=True)
ap.add_argument('-b', required=False, action='store_true')
ap.add_argument('-nl', required=False, action='store_true', help="Optional flag to process no-ligation bulks")
ap.add_argument('-bcf', required=False, default='/root/bc50.txt.gz', help='Optional flag to specify barcode list')
ap.add_argument('-bcs', required=False, default='22,30,60,68,117', help='Optional specify bc read structure as list [bc2_start, bc2_end, bc1_start, b1_end, genomicSeq_start]')
ap.add_argument('-cm', required=False, action='store_true', help='Optional flag to specify chromap-style bulk barcode reads.')

args = vars(ap.parse_args())
logging.info(f"args: {args}")

bc2_start, bc2_end, bc1_start, bc1_end, seq_start = list(map(int, args['bcs'].split(',')))
logging.info(f"barcode coordinates: {(bc2_start, bc2_end, bc1_start, bc1_end, seq_start)}")

input_file_R2 = args['i']
output_file_R3 = args['o3']
output_file_R2 = args['o2']


if args['b'] or args['nl'] or args['cm']:
    bc_list = gzopen(args['bcf'], 'rt').readlines()

with gzopen(input_file_R2, 'rt') as in_handle_R2, \
    open(output_file_R3, 'w') as out_handle_R3, \
        open(output_file_R2, 'w') as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle_R2):
        new_seq_R3 = seq[seq_start:]
        new_qual_R3 = qual[seq_start:]
        out_handle_R3.write(
            f'@{title}\n{new_seq_R3}\n+\n{new_qual_R3}\n'
        )
        if args['b'] or args['cm'] or args['nl']: # if we're doing bulk or we need a chromap-style output, get a random barcode from list
            barcode = random.choice(bc_list).rstrip()
        else:
            barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end]  #if we're not doing bulk, get the barcode from the read sequence
        if not args['cm']: #if we do not want a chromap style read, do what we used to do and make a read from the bc    
            new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]
            out_handle_R2.write(f'@{title}\n{barcode}\n+\n{new_qual_R2}\n')
        elif args['nl']:  #but if we do want a chromap style read, and it's a noligation one...
            leaderRead2 = makeLeader(barcode, merLen=8, fullLeader=False)  # we make a leader using the chosen random barcode
            outRead2 = leaderRead2+seq  # we pre-pend the read seq w/ the leader
            outRead2Qual = qual[:len(leaderRead2)]+qual  #we do a sloppy hack of copying initial read qualities from read to leader
            out_handle_R2.write(f'@{title}\n{outRead2}\n+\n{outRead2Qual}\n')  # we output the entry
        else: #we want chromap style and reads already have a full leader (ie. nl ==False), we remove leader and replace with synthetic leader
            leaderRead2 = makeLeader(barcode, merLen=8, fullLeader=True)  # we make a leader using the chosen random barcode
            outRead2 = leaderRead2+seq[len(leaderRead2):]  # we substitute real read bases with synth leader read bases
            outRead2Qual = qual  #we we keep the intact quality scores for all the bases
            out_handle_R2.write(f'@{title}\n{outRead2}\n+\n{outRead2Qual}\n')  # we output the entry
            
print("complete")
