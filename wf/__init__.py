import glob, os
import subprocess
from enum import Enum
from pathlib import Path
from typing import List, Optional, Union, Tuple
import glob
from latch import large_task, large_task, small_task, workflow
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

class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

class BarcodeFile(Enum):
    x50 = "bc50.txt"
    x50_old = "bc50_old.txt"
    x96 = "bc96.txt"


class ExtractBarcode(Enum):
    extB = "bc_process.py"
#    extB_old = "extractBarcode_old.py"



@large_task(retries=0)
def filter_task(
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
        "k=30",
        "mm=f",
        "rcomp=f",
        "restrictleft=103",
        "hdist=3",
        f"stats={l1_stats}",
        "literal=GTGGCCGATGTTTCGCATCGGCGTACGACT"
        ]
    if skip1:
            _bbduk1_cmd = [
            "echo",
            str("the first step skiped!!")]
            subprocess.run(["mv",
            r1.local_path,
            str(filtered_r1_l1)])
            subprocess.run(["mv",
            r2.local_path,
            str(filtered_r2_l1)])
            subprocess.run(["touch",
            str(l1_stats)]) 
             
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
        "k=30",
        "mm=f",
        "rcomp=f",
        "restrictleft=65",
        "hdist=3",
        f"stats={l2_stats}",
        "literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
    ]

    if skip2:
            _bbduk2_cmd = [
            "echo",
            str("the second step skiped!!")]
            subprocess.run(["mv",
            str(filtered_r1_l1),
            str(filtered_r1_l2)])
            subprocess.run(["mv",
            str(filtered_r2_l1),
            str(filtered_r2_l2)])
            subprocess.run(["touch",
            str(l2_stats)])



    subprocess.run(_bbduk2_cmd)

    return (
            LatchFile(
                str(filtered_r1_l2),
                f"latch:///chromap_outputs/{run_id}/chromap_inputs/{run_id}_S1_L001_R1_001.fastq.gz"
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


@large_task(retries=0)
def process_bc_task(
    r2: LatchFile,
    run_id: str,
    extract_bc_file: ExtractBarcode
) -> Tuple[LatchFile, LatchFile]:
    """ Process original read2: 
        save genomic portion as read3, 
        extract 16 bp barcode seqs and save as new read2
    """

    outdir = Path("chromap_inputs/").resolve()
    os.mkdir(outdir)

    new_r2 = Path(f"{outdir}/{run_id}_S1_L001_R2_001.fastq").resolve()
    r3 = Path(f"{outdir}/{run_id}_S1_L001_R3_001.fastq").resolve()
    extract_bc_path = Path(f"{extract_bc_file.value}").resolve()

    _bc_cmd = [
        "python",
        f"{extract_bc_file.value}",
        "-i",
        r2.local_path,
        "-o2",
        f"{str(new_r2)}",
        "-o3",
        f"{str(r3)}"


    ]


    try:
        subprocess.run(_bc_cmd)
    except subprocess.CalledProcessError as e:
        message(
            "error",
            {
                "title": f"Barcode Processing error",
                "body": f"Error surfaced:\n{str(e)}",
            },
        )
        raise RuntimeError(f"NanoPlot error: {e}")

    return (
        LatchFile(
        str(new_r2),
        f"latch:///chromap_outputs/{run_id}/chromap_inputs/{run_id}_S1_L001_R2_001.fastq"
         ),
        LatchFile(
        str(r3),
        f"latch:///chromap_outputs/{run_id}/chromap_inputs/{run_id}_S1_L001_R3_001.fastq"
         )
          )

        

@large_task(retries=0)
def chromap_alignment(
    r1: LatchFile,
    r2: LatchFile,
    r3: LatchFile,
    species: LatchDir,
    run_id: str,
    barcode_file: BarcodeFile

) -> LatchFile:
    """Run chromap alignment
    """

    barcode_path = Path(f"{barcode_file.value}").resolve()

    fragment = Path("aln.bed").resolve()
    reference = glob.glob(f"{species.local_path}/*.fa")[0]
    index = glob.glob(f"{species.local_path}/*.index")[0]
    
    _chromap_command = [

    "/root/chromap/chromap",
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
    r3.local_path,
    "--barcode-whitelist",
    f"{barcode_file.value}"


    ]

    subprocess.run(_chromap_command)
    

    # Ensure chromap ran correctly

    try:
        subprocess.run(_chromap_command)
    except subprocess.CalledProcessError as e:
        message(
            "error",
            {
                "title": f"chromap running error",
                "body": f"Error surfaced:\n{str(e)}",
            },
        )
        raise RuntimeError(f"NanoPlot error: {e}")


    return LatchFile(str(fragment), f"latch:///chromap_outputs/{run_id}/chromap_output/aln.bed")

@small_task
def bed2fragment(bed: LatchFile, run_id: str) -> LatchFile:

    outdir = Path("chromap_output/").resolve()
    os.mkdir(outdir)
    temp_file = Path(f"{outdir}/temp.bed").resolve()
    unzip_file = Path(f"{outdir}/fragments.tsv").resolve()
    output_file = Path(f"{outdir}/fragments.tsv.gz").resolve()

    subprocess.run([
    "echo",
    str("add -1 to barcodes and zip the file!!")])

    with open(str(temp_file), 'w') as fw:
          subprocess.run(['awk', 'BEGIN{FS=OFS=" "}{$4=$4"-1"}4', bed.local_path], stdout=fw)

    subprocess.run([
    "echo",
    str("make sure the out put file is tab delimited!")])

    with open(str(unzip_file), 'w') as fw:
          subprocess.run(['sed', 's/ /\t/g', str(temp_file)], stdout=fw)
    subprocess.run([
    "echo",
    str("use tabix bgzip to convert bed file into gz file")])

    with open(str(output_file), 'w') as fw:
          subprocess.run(['bgzip', '-c', str(unzip_file)], stdout=fw)
          

    return LatchFile(str(output_file), f"latch:///chromap_outputs/{run_id}/chromap_output/fragments.tsv.gz")




metadata = LatchMetadata(
    display_name="atlasXomics_pipeline_wf",
    author=LatchAuthor(
        name="Noori Sotoudeh",
        email="noorisotudeh@gmail.com",
        github="https://github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics",
    parameters={
        "species": LatchParameter(
            display_name="genome chromap index and fasta files Directory",
            description="Select reference genome for chromap alignment. make sure to use right directory, eg: /Chromap_refernces/Refdata_scATAC_MAESTRO_GRCm38_1.1.0",
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
        "extract_bc_file": LatchParameter(
            display_name="extract barcode file",
            description="this is the code will be use in step2 to extract barcodes from read2; \
                        extractBarcode_coprofile.py for current platform, \
                        extractBarcode_old.py for previous platform.",
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
        )

    },
)

@workflow(metadata)
def chromap_alignment_wf(
    r1: LatchFile,
    r2: LatchFile,
    run_id: str,
    skip1: bool,
    skip2: bool,
    species: LatchDir,
    barcode_file: BarcodeFile = BarcodeFile.x50,
    extract_bc_file: ExtractBarcode = ExtractBarcode.extB
) -> List[Union[LatchDir, LatchFile]]:
    """Pipeline for processing Spatial Whole Transcriptome data generated via DBiT-seq.

    Whole Transcriptome
    ----

    Process data from DBiT-seq experiments for spatially-resolved epigenomics:

    > See Deng, Y. et al 2022.

    # Steps

    * filter read2 on linker 1 identify via bbduk
    * filter read2 on linker 2 identify via bbduk
    * extract barcode (new read2) from read2
    * run chromap alignment and get fragment file
    """
    filtered_r1, filtered_r2, _, _ = filter_task(
        r1=r1,
        r2=r2,
        run_id=run_id,
        skip1=skip1,
        skip2=skip2
    )
    bc_fq, seq_fq = process_bc_task(
        r2=filtered_r2,
        run_id=run_id,
        extract_bc_file=extract_bc_file
    )

    chromap_outs = chromap_alignment(
        r1=filtered_r1,
        r2=seq_fq,
        r3=bc_fq,
        run_id=run_id,
        species=species,
        barcode_file=barcode_file
    )
    reports = bed2fragment(
        bed=chromap_outs,
        run_id=run_id
    )


    return [chromap_outs, reports]


LaunchPlan(
    chromap_alignment_wf,
    "default",
    {
        "r1" : LatchFile("latch:///noori/coProf/cp1.100cov.fq.gz"),
        "r2" : LatchFile("latch:///noori/coProf/cp2.100cov.fq.gz"),
        "run_id" : "ds_D01033_NG01681",
        "skip1" : False,
        "skip2" : False,
        "extract_bc_file" : ExtractBarcode.extB,
        "species" : LatchDir("latch:///Chromap_refernces/Refdata_scATAC_MAESTRO_GRCm38_1.1.0"),
        "barcode_file" : BarcodeFile.x50


    },
)
