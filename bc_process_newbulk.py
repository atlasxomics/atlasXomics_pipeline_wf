"""Process R2 for cellranger-atac or chromap pipeline
revised Dec 2023"""

import argparse
import gzip
import logging
import random

from Bio.SeqIO.QualityIO import FastqGeneralIterator


# define some functions and setup
logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def makeLeader(fullBC, merLen=8, fullLeader=True):
    fullL = f"CAAGCGTTGGCTTCTCGCATCT{fullBC[:(merLen)]}ATCCACGTGCTTGAGAGGCCAGAGCATTCG{fullBC[(merLen):]}GTGGCCGATGTTTCGCATCGGCGTACGACTAGATGTGTATAAGAGACAG"
    noLigL = f"CAAGCGTTGGCTTCTCGCATCT{fullBC[:(merLen)]}ATCCACGTGCTTGAGAGGCCAGAGCATTCG{fullBC[(merLen):]}GTGGCCGATGTTTCG"

    return fullL if fullLeader else noLigL


logging.info("Initializing Barcode Read Processing")

ap = argparse.ArgumentParser()
ap.add_argument("-i", required=True)
ap.add_argument(
    "-o2",
    required=True,
    help="Required path to output: Gzipped FastQ for read2"
)
ap.add_argument("-o3", required=True)
ap.add_argument("-b", required=False, action="store_true")
ap.add_argument(
    "-nl",
    required=False,
    action="store_true",
    help="Optional flag to process no-ligation bulks"
)
ap.add_argument(
    "-bcf",
    required=False,
    default="/root/bc50.txt.gz",
    help="Optional flag to specify barcode list"
)
ap.add_argument(
    "-bcs",
    required=False,
    default="22,30,60,68,117",
    help="Optional specify bc read structure as list [bc2_start, bc2_end, \
            bc1_start, b1_end, genomicSeq_start]"
)
ap.add_argument(
    "-cm",
    required=False,
    action="store_true",
    help="Optional flag to specify chromap-style bulk barcode reads."
)

args = vars(ap.parse_args())
logging.info(f"args: {args}")

bc2_start, bc2_end, bc1_start, bc1_end, seq_start = list(
    map(int, args["bcs"].split(","))
)
logging.info(
    f"barcode coordinates: \
        {(bc2_start, bc2_end, bc1_start, bc1_end, seq_start)}"
)

input_file_R2 = args["i"]
output_file_R3 = args["o3"]
output_file_R2 = args["o2"]

if args["b"] or args["nl"]:
    try:
        bc_list = gzip.open(args["bcf"], "rt").readlines()
    except gzip.BadGzipFile:
        logging.info(
            f"BC file {args['bcf']} is not gzipped. Opening as plain text."
        )
        bc_list = open(args["bcf"], "rt").readlines()

try:
    logging.info(f"Reading input Fastq file {input_file_R2} as gzip.")
    testRead = gzip.open(input_file_R2, "rt")
    pk = testRead.read(1)
    in_handle_R2 = gzip.open(input_file_R2, "rt")
except gzip.BadGzipFile:
    logging.info(
        f"Input Fastq file {input_file_R2} is not gzipped. Opening as plain \
            text."
    )
    in_handle_R2 = open(input_file_R2, "rt")

with (
    in_handle_R2,
    gzip.open(output_file_R3, "wt") as out_handle_R3,
    gzip.open(output_file_R2, "wt") as out_handle_R2
):

    readCount = 0
    for title, seq, qual in FastqGeneralIterator(in_handle_R2):
        readCount += 1

        if not (readCount % 500000):
            logging.info(f"Processed {readCount} reads")
        new_seq_R3 = seq[seq_start:]
        new_qual_R3 = qual[seq_start:]
        out_handle_R3.write(f"@{title}\n{new_seq_R3}\n+\n{new_qual_R3}\n")

        # if we're doing bulk or we need a chromap-style output, get a random
        # barcode from list
        if args["b"] or args["nl"]:
            barcode = random.choice(bc_list).rstrip()
        # if we're not doing bulk, get the barcode from the read sequence
        else:
            barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end]

        # if we do not want a chromap style read, do what we used to do and
        # make a read from the bc
        if not args["cm"]:
            new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]
            out_handle_R2.write(f'@{title}\n{barcode}\n+\n{new_qual_R2}\n')

        # but if we do want a chromap style read, and it's a noligation one...
        elif args["nl"]:
            # we make a leader using the chosen random barcode
            leaderRead2 = makeLeader(barcode, merLen=8, fullLeader=False)
            # we pre-pend the read seq w/ the leader
            outRead2 = leaderRead2+seq
            # we do a sloppy hack of copying initial read qualities from read
            # to leader
            outRead2Qual = qual[:len(leaderRead2)] + qual
            # we output the entry
            out_handle_R2.write(f"@{title}\n{outRead2}\n+\n{outRead2Qual}\n")

        # we want chromap style and reads already have a full leader
        # (ie. nl ==False), we remove leader and replace with synthetic leader
        else:
            # we make a leader using the chosen random barcode
            leaderRead2 = makeLeader(barcode, merLen=8, fullLeader=True)
            # we substitute real read bases with synth leader read bases
            outRead2 = leaderRead2+seq[len(leaderRead2):]
            # we we keep the intact quality scores for all the bases
            outRead2Qual = qual
            # we output the entry
            out_handle_R2.write(f"@{title}\n{outRead2}\n+\n{outRead2Qual}\n")

logging.info("complete")
