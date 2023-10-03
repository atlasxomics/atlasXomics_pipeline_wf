import glob, os
import subprocess
from enum import Enum
from pathlib import Path
from typing import List, Optional, Union, Tuple
import glob
from latch import large_task, medium_task, small_task, workflow
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
import pandas as pd
from ast import literal_eval
import logging



class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

class BarcodeFile(Enum):
    x50 = "bc50.txt"
    x50_old = "bc50_old.txt"
    x96 = "bc96.txt"


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
def chromap_alignment(
    r1: LatchFile,
    r2: LatchFile,
    species: LatchDir,
    run_id: str,
    barcode_file: BarcodeFile

) -> Tuple[LatchFile,LatchFile,LatchFile]:
    """Run chromap alignment
    """
    barcode_path = Path(f"{barcode_file.value}").resolve()

    fragment = Path("aln.bed").resolve()
    logfile = Path("chromap_log.txt").resolve()
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
    f"{barcode_file.value}",
    "--read-format",
    "bc:22:29,bc:60:67,r1:0:-1,r2:117:-1"

    ]
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=str(logfile)
        )
    try:
        # Open the log file for writing
        with open(str(logfile), 'w') as log_file:

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
                if output_line == '':
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

    subprocess.run([
    "echo",
    str("add -1 to barcodes and zip the file!!")])

    with open(str(temp_file), 'w') as fw:
          subprocess.run(['awk', 'BEGIN{FS=OFS=" "}{$4=$4"-1"}4', str(fragment)], stdout=fw)

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
          

    return (
           LatchFile(str(fragment), f"latch:///chromap_outputs/{run_id}/chromap_output/aln.bed"),
           LatchFile(str(logfile), f"latch:///chromap_outputs/{run_id}/chromap_output/chromap_log.txt"),
           LatchFile(str(output_file), f"latch:///chromap_outputs/{run_id}/chromap_output/fragments.tsv.gz")
           )



@large_task
def statistics(r2: LatchFile,frag: LatchFile, bed: LatchFile,logfile: LatchFile, species: LatchDir, run_id: str, barcode_file: BarcodeFile) -> LatchDir:
    work_dir = Path("statistics/").resolve()
    os.mkdir(work_dir)
    tmp_dir = Path("tmp_dir/").resolve()
    os.mkdir(tmp_dir)

    bc1 = Path(f"{work_dir}/bc1.txt").resolve()
    bc2 = Path(f"{work_dir}/bc2.txt").resolve()
    bc_freq_filtered = Path(f"{work_dir}/bc_freq_filtered.txt").resolve()
    bed_bc_freq = Path(f"{work_dir}/bed_bc_freq.txt").resolve()
    stats = Path(f"{work_dir}/stats.csv").resolve()
    whitelist =  Path(f"{barcode_file.value}").resolve()
#    fragment_edited = Path(f"{work_dir}/fragments_edited.tsv.gz").resolve()

    subprocess.run([
    "echo",
    str("find number of read per barcode from read_R2.fq.gz file!!")])
    

    command1 = ['zcat', r2.local_path]
    command2 = ['awk "NR%4==2"']
    command3 = ['cut', '-c', '61-68']
    command4 = ['cut', '-c', '23-30']

    # Open the output file for writing
    
    with open(str(bc1), 'w') as a, open(str(bc2), 'w') as b:
        process1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
        
        process2 = subprocess.Popen(command2, stdin=process1.stdout, stdout=subprocess.PIPE, shell=True)
        process1.stdout.close()
        
        process3 = subprocess.Popen(command3, stdin=process2.stdout, stdout=a)
        process2.stdout.close()

        process1.wait()
        process2.wait()
        process3.wait()

        process1 = subprocess.Popen(command1, stdout=subprocess.PIPE)

        process2 = subprocess.Popen(command2, stdin=process1.stdout, stdout=subprocess.PIPE, shell=True)
        process1.stdout.close()

        process4 = subprocess.Popen(command4, stdin=process2.stdout, stdout=b)
        process2.stdout.close()

        process1.wait()
        process2.wait()
        process4.wait()

        
        if (glob.glob(f"{species.local_path}/*.index")[0].split("/")[-1].split("_")[0]=='GRCm38'):
           genome='mm'
        else:
           genome='hs'

        if (genome=='mm'):
           chrsize=Path("mm10_chrom_sizes.txt").resolve()
        else:
           chrsize=Path("hg38_chrom_sizes.txt").resolve()

        if (genome=='mm'):
           black_lst=Path("blacklist/mm10-blacklist.v2.bed").resolve()
        else:
           black_lst=Path("blacklist/hg38-blacklist.v2.bed").resolve()


        _bc_cmd = [
                  "python",
                  "qc.py",
                  "-bc1",
                  str(bc1),
                  "-bc2",
                  str(bc2),
                  "-bcff",
                  str(bc_freq_filtered),
                  "-bed",
                  bed.local_path,
                  "-f",
                  frag.local_path,
                  "-bedbf",
                  str(bed_bc_freq),
                  "-i",
                  run_id,
                  "-g",
                  genome,
                  "-c",
                  chrsize,
                  "-k",
                  black_lst,
                  "-w",
                  barcode_file.value,
                  "-l",
                  logfile,
                  "-d",
                  work_dir,
                  "-t",
                  tmp_dir,
                  "-v",
                  open(Path("version").resolve(), 'r').read(),
                  "-o",
                  str(stats)
        ]

        subprocess.run(_bc_cmd)


        return LatchDir(str(work_dir), f"latch:///chromap_outputs/{run_id}/statistics")






metadata = LatchMetadata(
    display_name="atlasXomics_pipeline_wf",
    author=LatchAuthor(
        name="Noori Sotoudeh",
        email="noorisotudeh@gmail.com",
        github="https://github.com/atlasxomics/atlasXomics_pipeline_wf",
    ),
    repository="https://github.com/atlasxomics/atlasXomics_pipeline_wf",
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

    chromap_bed, chromap_log, chromap_frag = chromap_alignment(
        r1=filtered_r1,
        r2=filtered_r2,
        run_id=run_id,
        species=species,
        barcode_file=barcode_file
    )

    reports =  statistics(
        r2=r2,
        bed=chromap_bed,
        frag=chromap_frag,
        logfile=chromap_log,
        species=species,
        run_id=run_id,
        barcode_file=barcode_file
    )


    return [chromap_bed, chromap_frag, chromap_log, reports]


LaunchPlan(
    chromap_alignment_wf,
    "default",
    {
        "r1" : LatchFile("latch:///noori/coProf/cp1.100cov.fq.gz"),
        "r2" : LatchFile("latch:///noori/coProf/cp2.100cov.fq.gz"),
        "run_id" : "ds_D01033_NG01681",
        "skip1" : False,
        "skip2" : False,
        "species" : LatchDir("latch:///Chromap_refernces/Refdata_scATAC_MAESTRO_GRCm38_1.1.0")


    },
)

if __name__ == '__main__':
    statistics(
            r2="./chromap_R2.fq.gz",
            bed="aln.bed",
            logfile="chromap_log.txt",
            species="Refdata_scATAC_MAESTRO_GRCm38_1.1.0",
            run_id="test",
            barcode_file=BarcodeFile.x50
              )
