import logging
import pandas as pd


chromsize_paths = {
    "GRCh38":  "chrom_sizes/hg38_chrom_sizes.txt",
    "GRCm38": "chrom_sizes/mm10_chrom_sizes.txt",
    "Rnor6": "chrom_sizes/rn6_chrom_sizes.txt"
}


def call_oob(row: pd.Series, chromsizes: dict) -> bool:

    chrom, start, end = row["chrom"], row["start"], row["end"]
    if chrom in chromsizes:
        chrom_length = chromsizes[chrom]
        if start > chrom_length or end > chrom_length:
            logging.warning(
                f"Found one beyond {chrom=}, {chrom_length=}: \
                {(chrom, start, end)}"
            )
            return False
        else:
            return True
    else:
        logging.warning(f"Wrong chromosome found: {(chrom, start, end)}")
        return False


def filter_oob(input_file: str, chromsizes: dict) -> pd.DataFrame:
    """Read in the input file as a DataFrame and apply call_oob() to remove any
    fragments with start/stop greater than its corresponding chromosome sizes.
    """

    columns = ["chrom", "start", "end", "cellBarcode", "duplicates"]
    df = pd.read_csv(input_file, sep=" ", header=None, names=columns)

    condition = df.apply(lambda row: call_oob(row, chromsizes), axis=1)

    return df[condition]


def load_chromsizes(path: str) -> dict:

    sizes = {}
    with open(path, 'r') as f:
        for line in f:
            chrom, length = line.strip().split()
            sizes[chrom] = int(length)

    return sizes
