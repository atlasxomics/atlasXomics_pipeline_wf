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
ap.add_argument("-i", required=True, help="Path to input")
ap.add_argument(
    "-o2", required=True, help="Path to output: Gzipped FastQ for read2")
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
    default="/root/barcodes/bc50.txt.gz",
    help="Optional flag to specify barcode list"
)

args = vars(ap.parse_args())
logging.info(f"args: {args}")

logging.info("barcode coordinates: 22, 30, 60, 68, 117")

input_file_R2 = args["i"]
output_file_R2 = args["o2"]

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

with in_handle_R2, gzip.open(output_file_R2, "wt") as out_handle_R2:

    readCount = 0

    for title, seq, qual in FastqGeneralIterator(in_handle_R2):

        readCount += 1

        if not (readCount % 500000):
            logging.info(f"Processed {readCount} reads")

        read_len = len(seq)
        if read_len != 151:
            logging.warn(
                f"Read {readCount} is ({read_len}bp, not 151bp. Skipping \
                    bulkification for this read."
            )
            out_handle_R2.write(f"@{title}\n{seq}\n+\n{qual}\n")
            continue

        # get a random barcode from list
        barcode = random.choice(bc_list).rstrip()

        # if we do want a chromap style read, and it's a noligation one...
        if args["nl"]:
            # we make a leader using the chosen random barcode
            leaderRead2 = makeLeader(barcode, merLen=8, fullLeader=False)
            # we pre-pend the read seq w/ the leader
            outRead2 = leaderRead2 + seq
            # we do a sloppy hack of copying initial read qualities from read
            # to leader
            outRead2Qual = qual[:len(leaderRead2)] + qual
            # we output the entry
            out_handle_R2.write(f"@{title}\n{outRead2}\n+\n{outRead2Qual}\n")

        # if nl == False, we remove leader and replace with synthetic leader
        else:
            # we make a leader using the chosen random barcode
            leaderRead2 = makeLeader(barcode, merLen=8, fullLeader=True)
            # we substitute real read bases with synth leader read bases
            outRead2 = leaderRead2 + seq[len(leaderRead2):]
            # we we keep the intact quality scores for all the bases
            outRead2Qual = qual
            # we output the entry
            out_handle_R2.write(f"@{title}\n{outRead2}\n+\n{outRead2Qual}\n")

logging.info("complete")
