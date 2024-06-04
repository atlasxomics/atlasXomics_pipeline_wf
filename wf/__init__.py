"""AtlasXomics/Latch workflow for DBiT-seq epigenomic preprocessing
"""

import glob
import logging
import os
import subprocess

from enum import Enum
from pathlib import Path
from typing import List, Optional, Union, Tuple

from latch import custom_task, large_task, medium_task, small_task, workflow
from latch.account import Account
from latch.functions.messages import message
from latch.resources.launch_plan import LaunchPlan
from latch.registry.table import Table
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule
)

import wf.lims as lims

from wf.outliers import plotting_task


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


class BarcodeFile(Enum):
    x50 = "bc50.txt"
    x50_old = "bc50_old.txt"
    x96 = "bc96.txt"
    x96_fg = "bc96_fg.txt"
    x210 = "bcFG210v4.txt"
    x220 = "bc220-20-MAY.txt"


@large_task(retries=0)
def filtering(
    r1: LatchFile,
    r2: LatchFile,
    run_id: str,
    skip1: bool,
    skip2: bool
) -> Tuple[LatchFile, LatchFile, LatchFile, LatchFile]:

    filtered_r1_l1 = Path(f"{run_id}_linker1_R1.fastq.gz").resolve()
    filtered_r2_l1 = Path(f"{run_id}_linker1_R2.fastq.gz").resolve()
    l1_stats = Path(f"{run_id}_l1_stats.txt").resolve()

    _bbduk1_cmd = [
        "bbmap/bbduk.sh",
        f"in1={r1.local_path}",
        f"in2={r2.local_path}",
        f"outm1={str(filtered_r1_l1)}",
        f"outm2={str(filtered_r2_l1)}",
        "skipr1=t",
        "threads=96",
        "-Xmx196g",
        "k=30",
        "mm=f",
        "rcomp=f",
        "restrictleft=103",
        "hdist=3",
        f"stats={l1_stats}",
        "literal=GTGGCCGATGTTTCGCATCGGCGTACGACT"
    ]

    if skip1:
        _bbduk1_cmd = ["echo", str("the first step skipped!!")]
        subprocess.run(["mv", r1.local_path, str(filtered_r1_l1)])
        subprocess.run(["mv", r2.local_path, str(filtered_r2_l1)])
        subprocess.run(["touch", str(l1_stats)])

    subprocess.run(_bbduk1_cmd)

    filtered_r1_l2 = Path(f"{run_id}_linker2_R1.fastq.gz").resolve()
    filtered_r2_l2 = Path(f"{run_id}_linker2_R2.fastq.gz").resolve()
    l2_stats = Path(f"{run_id}_l2_stats.txt").resolve()

    _bbduk2_cmd = [
        "bbmap/bbduk.sh",
        f"in1={str(filtered_r1_l1)}",
        f"in2={str(filtered_r2_l1)}",
        f"outm1={str(filtered_r1_l2)}",
        f"outm2={str(filtered_r2_l2)}",
        "skipr1=t",
        "threads=96",
        "-Xmx196g",
        "k=30",
        "mm=f",
        "rcomp=f",
        "restrictleft=65",
        "hdist=3",
        f"stats={l2_stats}",
        "literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
    ]

    if skip2:
        _bbduk2_cmd = ["echo", str("the second step skipped!!")]
        subprocess.run(["mv", str(filtered_r1_l1), str(filtered_r1_l2)])
        subprocess.run(["mv", str(filtered_r2_l1), str(filtered_r2_l2)])
        subprocess.run(["touch", str(l2_stats)])

    subprocess.run(_bbduk2_cmd)

    return (
        LatchFile(
            str(filtered_r1_l2),
            f"latch:///chromap_outputs/{run_id}/preprocessing/{run_id}_S1_L001_R1_001.fastq.gz"
        ),
        LatchFile(
            str(filtered_r2_l2),
            f"latch:///chromap_outputs/{run_id}/preprocessing/{run_id}_linker2_R2.fastq.gz"
        ),
        LatchFile(
            str(l1_stats),
            f"latch:///chromap_outputs/{run_id}/preprocessing/{l1_stats.name}"
        ),
        LatchFile(
            str(l2_stats),
            f"latch:///chromap_outputs/{run_id}/preprocessing/{l2_stats.name}"
        )
    )


@medium_task(retries=0)
def process_bc_task(
    r2: LatchFile,
    run_id: str,
    bulk: bool,
    noLigation_bulk: bool,
    barcode_file: BarcodeFile
) -> LatchFile:
    """ Process read2: save genomic portion as read3, extract 16 bp
    barcode seqs and save as read3
    """
    if (not (bulk or noLigation_bulk)):
        return LatchFile(r2.local_path, r2.remote_path)

    outdir = Path("chromap_inputs/").resolve()
    os.mkdir(outdir)

    new_r2 = Path(f"{outdir}/{run_id}_S1_L001_R2_001.fastq").resolve()
    r3 = Path(f"{outdir}/{run_id}_S1_L001_R3_001.fastq").resolve()

    _bc_cmd = [
        "python",
        "scripts/bc_process_newbulk.py",
        "-i",
        r2.local_path,
        "-o2",
        f"{str(new_r2)}",
        "-o3",
        f"{str(r3)}",
        "-bcf",
        f"barcodes/{barcode_file.value}"
    ]

    if bulk:
        _bc_cmd.append("-b")
        _bc_cmd.append("-cm")
    if noLigation_bulk:
        _bc_cmd.append("-nl")
        _bc_cmd.append("-cm")

    subprocess.run(_bc_cmd)

    return LatchFile(
            str(new_r2),
            f"latch:///chromap_outputs/{run_id}/preprocessing/{run_id}_bulkd_R2.fastq.gz"
        )


@large_task(retries=0)
def alignment(
    r1: LatchFile,
    r2: LatchFile,
    species: LatchDir,
    run_id: str,
    barcode_file: BarcodeFile
) -> Tuple[LatchFile, LatchFile, LatchFile, LatchFile]:
    """Run chromap alignment
    """

    logfile = Path("chromap_log.txt").resolve()
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=str(logfile)
    )

    fragment = Path("aln.bed").resolve()
    reference = glob.glob(f"{species.local_path}/*.fa")[0]
    index = glob.glob(f"{species.local_path}/*.index")[0]

    _chromap_command = [
        "/root/chromap/chromap",
        "-t",
        "96",
        "--preset",
        "atac",
        "-x",
        str(index),
        "-r",
        str(reference),
        "-1",
        r1.local_path,
        "-2",
        r2.local_path,
        "-o",
        str(fragment),
        "-b",
        r2.local_path,
        "--barcode-whitelist",
        f"barcodes/{barcode_file.value}",
        "--read-format",
        "bc:22:29,bc:60:67,r1:0:-1,r2:117:-1"
    ]

    try:
        # Open the log file for writing
        with open(str(logfile), "w") as log_file:

            # Start the subprocess with stdout redirected to a pipe
            process = subprocess.Popen(
                _chromap_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                shell=False
            )

            # Continuously read and log the output as it becomes available
            while True:
                output_line = process.stdout.readline()
                if output_line == "":
                    break  # End of output
                log_file.write(output_line)
                log_file.flush()
                logging.info(output_line.strip())

            # Wait for the subprocess to finish
            process.wait()

    except Exception as e:
        logging.error(f"An error occurred: {e}")

    outdir = Path("chromap_output/").resolve()
    os.mkdir(outdir)

    temp_file = Path(f"{outdir}/temp.bed").resolve()
    unzip_file = Path(f"{outdir}/fragments.tsv").resolve()
    output_file = Path(f"{outdir}/fragments.tsv.gz").resolve()
    output_file_index = Path(f"{outdir}/fragments.tsv.gz.tbi").resolve()

    subprocess.run(["echo", str("add -1 to barcodes and zip the file!!")])
    with open(str(temp_file), "w") as fw:
        subprocess.run(  # awk expression needs to be ""; don't change...
            ["awk", 'BEGIN{FS=OFS=" "}{$4=$4"-1"}4', str(fragment)],
            stdout=fw
        )

    subprocess.run(
        ["echo", str("make sure the out put file is tab delimited!")]
    )

    with open(str(unzip_file), "w") as fw:
        subprocess.run(["sed", "s/ /\t/g", str(temp_file)], stdout=fw)

    subprocess.run(
        ["echo", "use tabix bgzip to convert bed file into gz file"]
    )

    with open(str(output_file), "w") as fw:
        subprocess.run(["bgzip", "-c", str(unzip_file)], stdout=fw)

    with open(str(output_file_index), "w") as fw:
        subprocess.run(
            ["tabix", "-p", "bed", str(output_file), "-f"],  stdout=fw
        )

    return (
        LatchFile(
            str(fragment),
            f"latch:///chromap_outputs/{run_id}/chromap_output/aln.bed"
        ),
        LatchFile(
            str(logfile),
            f"latch:///chromap_outputs/{run_id}/chromap_output/chromap_log.txt"
        ),
        LatchFile(
            str(output_file),
            f"latch:///chromap_outputs/{run_id}/chromap_output/fragments.tsv.gz"
        ),
        LatchFile(
            str(output_file_index),
            f"latch:///chromap_outputs/{run_id}/chromap_output/fragments.tsv.gz.tbi"
        )
    )


def allocate_mem(
    r2: LatchFile,
    frag: LatchFile,
    bed: LatchFile,
    logfile: LatchFile,
    species: LatchDir,
    run_id: str,
    barcode_file: BarcodeFile
) -> int:
    bigs = ["bcFG210v4.txt", "bc220-20-MAY.txt"]
    return 750 if barcode_file.value in bigs else 192


@custom_task(cpu=30, memory=allocate_mem, storage_gib=500)
def statistics(
    r2: LatchFile,
    frag: LatchFile,
    bed: LatchFile,
    logfile: LatchFile,
    species: LatchDir,
    run_id: str,
    barcode_file: BarcodeFile
) -> LatchDir:

    work_dir = Path("Statistics/").resolve()
    os.mkdir(work_dir)

    tmp_dir = Path("tmp_dir/").resolve()
    os.mkdir(tmp_dir)

    whitelist = Path(f"barcodes/{barcode_file.value}").resolve()
    singlecell = Path(f"{work_dir}/singlecell.csv").resolve()

    positions_paths = {
        "x50": "s3://latch-public/test-data/13502/x50_all_tissue_positions_list.csv",
        "x50_old": "s3://latch-public/test-data/13502/x50-old_tissue_positions_list.csv",
        "x96": "s3://latch-public/test-data/13502/x96_all_tissue_positions_list.csv",
        "x96_fg": "s3://latch-public/test-data/13502/xfg96_11DEC_all_tissue_positions_list.csv",
        "x210": "s3://latch-public/test-data/13502/xFG210v4_all_tissue_positions_list.csv",
        "x220": "s3://latch-public/test-data/13502/xbc220-20-MAY_alltissue_positions_list.csv"
    }
    positions_path = LatchFile(positions_paths[barcode_file.name])
    positions_file = Path(positions_path.local_path).resolve()

    # Extract genome name from local genome ref dir
    genome_id = (
        glob.glob(f"{species.local_path}/*.fa")[0].split("/")[-1].split("_")[0]
    )

    # Assign genome metadata
    genome_dict = {
        "GRCh38": [
            "hs",
            "chrom_sizes/hg38_chrom_sizes.txt",
            "blacklist/hg38-blacklist.v2.bed"
        ],
        "GRCm38": [
            "mm",
            "chrom_sizes/mm10_chrom_sizes.txt",
            "blacklist/mm10-blacklist.v2.bed"
        ],
        "Rnor6": [
            "rnor6",
            "chrom_sizes/rn6_chrom_sizes.txt",
            "na"
        ]
    }

    _pyct_cmd = [
        "python",
        "scripts/pycis.py",
        "-f",
        frag.local_path,
        "-i",
        run_id,
        "-g",
        genome_dict[genome_id][0],
        "-c",
        Path(genome_dict[genome_id][1]).resolve(),
        "-k",
        Path(genome_dict[genome_id][2]).resolve(),
        "-w",
        whitelist,
        "-d",
        work_dir,
        "-t",
        tmp_dir,
        "-p",
        positions_file
    ]

    subprocess.run(_pyct_cmd)

    _sc_cmd = [
        "python",
        "scripts/singlecellsummary.py",
        "-r2",
        r2.local_path,
        "-bed",
        bed.local_path,
        "-w",
        whitelist,
        "-i",
        run_id,
        "-g",
        genome_dict[genome_id][0],
        "-l",
        logfile,
        "-v",
        open(Path("version").resolve(), "r").read()
    ]

    subprocess.run(_sc_cmd)

    logging.info("Creating lane qc plot.")
    plotting_task(singlecell, positions_file)

    logging.info("Cleaning up /Statistics/.")
    _cp_cmd = [
        "cp",
        "/root/Statistics/scATAC/quality_control/sample_metrics.pdf",
        "/root/Statistics/qc_plot.pdf"
    ]
    subprocess.run(_cp_cmd)

    try:  # Put in try/except so not to break if missing
        os.remove("/root/Statistics/fragments_edited.tsv.gz")
        os.remove("/root/Statistics/tmp1.txt")
    except FileNotFoundError:
        pass

    return LatchDir(
        str(work_dir),
        f"latch:///chromap_outputs/{run_id}/Statistics"
    )


@small_task(retries=0)
def lims_task(
    results_dir: LatchDir,
    run_id: str,
    upload: bool,
    ng_id: Optional[str]
) -> LatchDir:

    if upload:

        data = Path(results_dir.local_path + '/summary.csv').resolve()

        slims = lims.slims_init()
        results = lims.csv_to_dict(data)

        payload = {lims.mapping[key]: value for (key, value) in results.items()
                   if key in lims.mapping.keys() and value not in lims.NA}

        if ng_id:
            pk = lims.get_pk(ng_id, slims)
        else:
            try:
                pk = lims.get_pk(run_id.split("_")[-1], slims)
            except IndexError:
                logging.warning("Invalid SLIMS ng_id.")
                message(
                    typ="warning",
                    data={
                        "title": "SLIMS fail", "body": "Invalid SLIMS ng_id."
                    }
                )
                return results_dir

        payload["rslt_fk_content"] = pk
        payload["rslt_fk_test"] = 39
        payload["rslt_value"] = "upload"

        logging.info(f"upload succeeded: {lims.push_result(payload, slims)}")

        return results_dir

    return results_dir


@small_task(retries=0)
def upload_latch_registry(
    run_id: str,
    r1: LatchFile,
    r2: LatchFile,
    chromap_frag: LatchFile,
    results_dir: LatchDir,
    table_id: str = "761"
):

    acc = Account.current()

    if acc.id == "13502":

        table = Table(table_id)

        prefix = f"{results_dir.remote_path}/"

        peaks_bed = f"{prefix}{run_id}_peaks.bed"
        single_cell_file = f"{prefix}/singlecell.csv"

        fragments_file_tbi = f"{chromap_frag.remote_path}.tbi"

        with table.update() as updater:
            try:
                updater.upsert_record(
                    run_id,
                    fastq_read_1=r1,
                    fastq_read_2=r2,
                    fragments_file=chromap_frag,
                    peaks_bed=LatchFile(peaks_bed),
                    fragment_file_tbi=LatchFile(fragments_file_tbi),
                    single_cell_file=LatchFile(single_cell_file)
                )
                return
            except TypeError:
                logging.warning(f"No table with id {table_id} found.")
                return
            finally:
                return
    else:
        return


metadata = LatchMetadata(
    display_name="ATX epigenomic preprocessing",
    author=LatchAuthor(
        name="Noori Sotudeh",
        email="noorisotudeh@gmail.com",
        github="https://github.com/atlasxomics/atlasXomics_pipeline_wf",
    ),
    repository="https://github.com/atlasxomics/atlasXomics_pipeline_wf",
    parameters={
        "species": LatchParameter(
            display_name="Chromap genome directory",
            description="Select reference genome for chromap alignment. Make \
                    sure to use right directory, eg. \
                    '/Chromap_references/Refdata_scATAC_MAESTRO_GRCm38_1.1.0'",
            batch_table_column=True,
        ),
        "r1": LatchParameter(
            display_name="read 1",
            description="Read 1 must contain genomic sequence.",
            batch_table_column=True,
        ),
        "r2": LatchParameter(
            display_name="read 2",
            description="Read 2 must contain barcode sequences \
                        and end with >35bp of genomic sequence.",
            batch_table_column=True,
        ),
        "skip1": LatchParameter(
            display_name="skip linker1 filtering",
            description="If True, will skip first step linker1 filtering!",
            batch_table_column=True,
        ),
        "skip2": LatchParameter(
            display_name="skip linker2 filtering",
            description="If True, will skip first step linker2 filtering!",
            batch_table_column=True,
        ),
        "run_id": LatchParameter(
            display_name="run id",
            description="ATX Run ID with optional prefix, default to \
                        Dxxxxx_NGxxxxx format.",
            batch_table_column=True,
            placeholder="Dxxxxx_NGxxxxx",
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="run id cannot start with a '/'"
                ),
                LatchRule(
                    regex="_NG[0-9]{5}$",
                    message="Provide ng_id in ng_id field if upload to \
                    SLIMS desired."
                )
            ]
        ),
        "barcode_file": LatchParameter(
            display_name="barcode file",
            description="Expected sequences of barcodes used is assay; \
                        bc50.txt for SOP 50x50, bc96.txt for 96x96, \
                        bc50_old.txt for previous version of 50x50.",
            batch_table_column=True,
        ),
        "bulk": LatchParameter(
            display_name="bulk",
            description=" If true, barcode sequences in reads will be \
                        replaced with random barcodes (current: bulk).",
            batch_table_column=True,
        ),
        "noLigation_bulk": LatchParameter(
            display_name="no-ligation primer bulk",
            description="If true, reads will have linker sequences added and \
                        random barcodes assigned (current: no-ligation \
                        primer bulk).",
            batch_table_column=True,
        ),
        "upload_to_slims": LatchParameter(
            display_name="upload to slims",
            description="Select for run metrics to be upload to SLIMS; if \
                        selected provide ng_id.",
            batch_table_column=True,
        ),
        "ng_id": LatchParameter(
            display_name="ng_id",
            description="Provide SLIMS ng_id (ie. NG00001) if pushing to \
                        SLIMS and run_id does not end in '_NG00001'.",
            placeholder="NGxxxxx",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^NG[0-9]{5}$",
                    message="ng_id must match NGxxxxx format."
                ),
                LatchRule(
                    regex="^\S+$",
                    message="run id cannot contain whitespace"
                )
            ]
        ),
        "table_id": LatchParameter(
            display_name="Registry Table ID",
            description="Provide the ID of the Registry table. Files that \
                        will be populated in the table are: singlecell.csv, \
                        fragments.tsv.gz, and summary.csv"
        )
    },
)


@workflow(metadata)
def total_wf(
    r1: LatchFile,
    r2: LatchFile,
    run_id: str,
    skip1: bool,
    skip2: bool,
    species: LatchDir,
    ng_id: Optional[str],
    barcode_file: BarcodeFile = BarcodeFile.x50,
    noLigation_bulk: bool = False,
    bulk: bool = False,
    upload_to_slims: bool = False,
    table_id: str = "761"
) -> List[Union[LatchDir, LatchFile]]:
    """Workflow for processing epigenomic data generated via DBiT-seq

    ATX epigenomic preprocessing
    ----
    The ATX epigenomic preprocessing workflow requires the barcoding schema
    described in [Zhang et al. 2023](https://www.nature.com/articles/s41586-023-05795-1#MOESM1)
    for Illumina short-read sequencing:
    - read1: genomic sequence
    - read2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence

    The workflow is comprised of the following steps:

    1. Filter reads on read2 ligation linker sequences with
    bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/).
        1. reads with >3 mismatches in the ligation linker1 are removed
        from analysis
        2. reads with >3 mismatches in the ligation linker2 are removed
        from analysis

    2. Perform sequence alignment with [Chromap](https://github.com/haowenz/chromap).
        - Chromap will process barcodes from fastq_R2.gz and align reads of both
        fastq_R1.gz and fastq_R2.gz files. The result is a fragment file which
        can be used for downstream analysis by [ArchR](https://www.archrproject.com/)
        or [Signac](https://stuartlab.org/signac/).

    3. Generate quality control metrics (ie. FRIP, TSS score) and peak
    .bed/.h5 files using the [pycisTopic](https://github.com/aertslab/pycisTopic)package.

    Questions? Comments?  Contact support@atlasxomics.com or post in the
    AtlasXomics [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).
    """

    filtered_r1, filtered_r2, _, _ = filtering(
        r1=r1,
        r2=r2,
        run_id=run_id,
        skip1=skip1,
        skip2=skip2
    )

    bulkd_r2 = process_bc_task(
        r2=filtered_r2,
        run_id=run_id,
        bulk=bulk,
        noLigation_bulk=noLigation_bulk,
        barcode_file=barcode_file
    )

    chromap_bed, chromap_log, chromap_frag, chromap_index = alignment(
        r1=filtered_r1,
        r2=bulkd_r2,
        run_id=run_id,
        species=species,
        barcode_file=barcode_file
    )

    reports = statistics(
        r2=bulkd_r2,
        bed=chromap_bed,
        frag=chromap_frag,
        logfile=chromap_log,
        species=species,
        run_id=run_id,
        barcode_file=barcode_file
    )

    lims_task(
        results_dir=reports,
        run_id=run_id,
        upload=upload_to_slims,
        ng_id=ng_id
    )

    upload_latch_registry(
        run_id=run_id,
        r1=r2,
        r2=r1,
        chromap_frag=chromap_frag,
        results_dir=reports,
        table_id=table_id
    )

    return [chromap_bed, chromap_frag, chromap_log, chromap_index, reports]


LaunchPlan(
    total_wf,
    "demo",
    {
        "r1": LatchFile(
            "s3://latch-public/test-data/13502/atx_demo_R1_001.fastq.gz"
        ),
        "r2": LatchFile(
            "s3://latch-public/test-data/13502/atx_demo_R2_001.fastq.gz"
        ),
        "run_id": "demo",
        "skip1": False,
        "skip2": False,
        "species": LatchDir("latch:///Chromap_references/Human")
    }
)

if __name__ == "__main__":

    r2 = LatchFile("latch://13502.account/downsampled/D01033_NG01681/ds_D01033_NG01681_S3_L001_R2_001.fastq.gz")
    species = LatchDir("latch://13502.account/Chromap_references/Human")
    bed = LatchFile("latch://13502.account/chromap_outputs/demo/chromap_output/aln.bed")
    frag = LatchFile("latch://13502.account/chromap_outputs/demo/chromap_output/fragments.tsv.gz")
    logfile = LatchFile("latch://13502.account/chromap_outputs/demo/chromap_output/chromap_log.txt")

    statistics(
        r2=r2,
        species=species,
        run_id="sc_debug",
        barcode_file=BarcodeFile.x50,
        bed=bed,
        frag=frag,
        logfile=logfile
    )
