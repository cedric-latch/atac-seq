"""latch/ATACseq"""

import _io
import functools
import os
import shutil
import subprocess
import types
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union

import lgenome
from dataclasses_json import dataclass_json
from flytekit import task
from flytekitplugins.pod import Pod
from kubernetes.client.models import (V1Container, V1PodSpec,
                                      V1ResourceRequirements, V1Toleration)
from latch import (map_task, message, workflow)
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile, file_glob

print = functools.partial(print, flush=True)

# Available latch genomes
class LatchGenome(Enum):
    RefSeq_hg38_p14 = "Homo sapiens (RefSeq hg38.p14)"
    RefSeq_T2T_CHM13v2_0 = "Homo sapiens (RefSeq T2T-CHM13v2.0)"
    RefSeq_R64 = "Saccharomyces cerevisiae (RefSeq R64)"
    RefSeq_GRCm39 = "Mus musculus (RefSeq GRCm39)"

# Handling of single vs paired-end reads
@dataclass_json
@dataclass
class SingleEndReads:
    r1: LatchFile

@dataclass_json
@dataclass
class PairedEndReads:
    r1: LatchFile
    r2: LatchFile

class ReadType(Enum):
    single = "single"
    paired = "paired"

class Strandedness(Enum):
    auto = "auto"

Replicate = Union[SingleEndReads, PairedEndReads]

@dataclass_json
@dataclass
class Sample:
    name: str
    strandedness: Strandedness
    replicates: List[Replicate]

@dataclass_json
@dataclass
class ATACseqInput:
    sample_name: str
    replicates: List[Replicate]
    run_name: str
    base_remote_output_dir: str
    genome_index: LatchDir
    genome_size: int
    clip_r1: Optional[int] = None
    clip_r2: Optional[int] = None
    three_prime_clip_r1: Optional[int] = None
    three_prime_clip_r2: Optional[int] = None

@dataclass_json
@dataclass
class ATACseqOutput:
    sample_name: str
    trimgalore_reports: List[LatchFile]
    fastqc_zip: List[LatchFile]
    fastqc_html: List[LatchFile]
    bam_file: LatchFile
    bed_file: LatchFile
    bdg_file: LatchFile

#
# Task specific exception classes
#

class TrimgaloreError(Exception):
    pass

class Bowtie2IndexError(Exception):
    pass

class Bowtie2AlignmentError(Exception):
    pass

class MalformedBowtie2Index(Exception):
    pass
        
class Macs3Error(Exception):
    pass

# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)

LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)

def slugify(value: str) -> str:
    return value.lower().replace(" ", "_")

def _remote_output_dir(custom_output_dir: Optional[LatchDir]) -> str:
    if custom_output_dir is None:
        return "/ATAC-Seq Outputs/"
    remote_path = custom_output_dir.remote_path
    assert remote_path is not None
    if remote_path[-1] != "/":
        remote_path += "/"
    return remote_path

def _get_96_spot_pod() -> Pod:
    """[ "c6i.24xlarge", "c5.24xlarge", "c5.metal", "c5d.24xlarge", "c5d.metal" ]"""

    primary_container = V1Container(name="primary")
    resources = V1ResourceRequirements(
        requests={"cpu": "90", "memory": "170Gi"},
        limits={"cpu": "96", "memory": "192Gi"},
    )
    primary_container.resources = resources

    return Pod(
        pod_spec=V1PodSpec(
            containers=[primary_container],
            tolerations=[
                V1Toleration(effect="NoSchedule", key="ng", value="cpu-96-spot")
            ],
        ),
        primary_container_name="primary",
    )

large_spot_task = task(task_config=_get_96_spot_pod(), retries=3)

def piped_process(command: List[str], stdin: _io.BufferedReader) -> subprocess.Popen:
    print(command)
    return subprocess.Popen(
        command,
        stdin=stdin,
        stdout=subprocess.PIPE,
    )

def piped_run(commands: List[List[str]]) -> Tuple[List[int], str]:
    processes = [piped_process(commands[0], stdin=None)]

    for cmd in commands[1:]:
        processes.append(
            piped_process(cmd, stdin=processes[-1].stdout)
        )

    (stdout, _stderr) = processes[-1].communicate()    
    return_code = processes[-1].returncode

    return (return_code, stdout.decode())

def _merge_replicates(
    replicates: List[Replicate], sample_name: str
) -> Union[Tuple[Path], Tuple[Path, Path]]:
    local_r1_path = f"{slugify(sample_name)}_r1_merged.fq"
    r1 = _concatenate_files((str(x.r1.path) for x in replicates), local_r1_path)

    if isinstance(replicates[0], SingleEndReads):
        return (r1,)

    if not all(isinstance(x, PairedEndReads) for x in replicates):
        raise RuntimeError("Not all technical replicates were paired end")

    local_r2_path = f"{slugify(sample_name)}_r2_merged.fq"
    r2 = _concatenate_files((str(x.r2.path) for x in replicates), local_r2_path)
    return (r1, r2)

def _concatenate_files(filepaths: Iterable[str], output_path: str) -> Path:
    path = Path(output_path).resolve()
    with path.open("w") as output_file:
        for p in filepaths:
            p = p.removesuffix(".gz")
            with open(p, "r") as f:
                shutil.copyfileobj(f, output_file)
    return path

@large_spot_task
def prepare_inputs(
    samples: List[Sample],
    run_name: str,
    latch_genome: LatchGenome,
    clip_r1: Optional[int] = None,
    clip_r2: Optional[int] = None,
    three_prime_clip_r1: Optional[int] = None,
    three_prime_clip_r2: Optional[int] = None,
    custom_output_dir: Optional[LatchDir] = None,
    custom_ref_genome: Optional[LatchFile] = None,
) -> List[ATACseqInput]:

    # Compute bowtie2 index
    (genome_index, gsize) = build_bowtie2_index(
        latch_genome,
        custom_ref_genome
    )

    # Bundle parameters before sending to main process
    return [
        ATACseqInput(
            sample_name=sample.name,
            replicates=sample.replicates,
            run_name=run_name,
            base_remote_output_dir=_remote_output_dir(custom_output_dir),
            genome_index=LatchDir(str(genome_index)),
            genome_size=gsize,
            clip_r1=clip_r1,
            clip_r2=clip_r2,
            three_prime_clip_r1=three_prime_clip_r1,
            three_prime_clip_r2=three_prime_clip_r2,
        )
        for sample in samples
    ]

def do_trimgalore(
    input: ATACseqInput,
    replicate_index: int,
    reads: Replicate,
) -> Tuple[Replicate, List[LatchFile], List[LatchFile], List[LatchFile]]:
    def _flag(name: str) -> List[str]:
        value = getattr(input, name)
        return [f"--{name}", value] if value is not None else []

    flags = [*_flag("clip_r1"), *_flag("three_prime_clip_r1")]
    read_paths = [reads.r1.local_path]
    if isinstance(reads, PairedEndReads):
        flags += ["--paired", *_flag("clip_r2"), *_flag("three_prime_clip_r2")]
        read_paths.append(reads.r2.local_path)

    local_output = f"{slugify(input.sample_name)}_replicate_{replicate_index}"

    return_code, stdout = piped_run(
        [[
            "trim_galore",
            "--cores",
            str(8),
            "--fastqc",
            "--dont_gzip",
            "--output_dir",
            f"./{local_output}",
            *flags,
            *read_paths,
        ]]
    )

    if return_code != 0:
        stdout = stdout.rstrip()
        stdout = stdout[stdout.rindex("\n") + 1 :]
        assert reads.r1.remote_path is not None
        path_name = reads.r1.remote_path.split("/")[-1]
        identifier = f"sample {input.sample_name}, replicate {path_name}"
        message(
            "error",
            {"title": f"Trimgalore error for {identifier}", "body": stdout},
        )
        raise TrimgaloreError(stdout)

    def remote(middle: str) -> str:
        base = f"{input.base_remote_output_dir}{input.run_name}"
        tail = f"{input.sample_name}/replicate_{replicate_index}/"
        return f"latch:///{base}/Quality Control/Trimming {middle}/{tail}"

    reads_directory = remote("Reads")
    if isinstance(reads, SingleEndReads):
        (r1,) = file_glob(f"{local_output}/*trimmed.fq*", reads_directory)
        trimmed_replicate = SingleEndReads(r1=r1)
    else:
        # File glob sorts files alphanumerically
        r1, r2 = file_glob(f"{local_output}/*val*.fq*", reads_directory)
        trimmed_replicate = PairedEndReads(r1=r1, r2=r2)

    # Delete unneeded files to free disk space
    os.remove(reads.r1.local_path)
    if isinstance(reads, PairedEndReads):
        os.remove(reads.r2.local_path)

    reports_directory = remote("Reports")
    reports = file_glob(f"{local_output}/*trimming_report.txt", reports_directory)

    fastqc_directory = remote("FastQC")
    fastqc_html = file_glob(f"{local_output}/*.html", fastqc_directory)
    fastqc_zip = file_glob(f"{local_output}/*.zip", fastqc_directory)

    return trimmed_replicate, reports, fastqc_zip, fastqc_html

def build_bowtie2_index(
    latch_genome: LatchGenome,
    custom_ref_genome: Optional[LatchFile],
) -> (Path, int):
    # 
    # Build or fetch bowtie2 index
    #

    if custom_ref_genome is None:
        # Fetch the genome if already in db
        gm = lgenome.GenomeManager(latch_genome.name)
        custom_ref_genome = gm.download_ref_genome()
    else:
        custom_ref_genome = custom_ref_genome.local_path()

    gsize = sum(len(line.strip()) for line in open(custom_ref_genome)
                if not line.startswith(">"))

    # # Run bowtie2-build
    local_index_dir = Path("bowtie2_index").resolve()
    
    print(f"Building bowtie2 index for {custom_ref_genome}")
    Path(local_index_dir).mkdir(exist_ok=True)

    try:
        (exit_code, stdout) = piped_run([
            ["bowtie2-build",
            custom_ref_genome,
            f"{local_index_dir}/ref",
            "--threads",
            "96"]
        ])
    except subprocess.CalledProcessError as e:
        exit_code = 1
        stdout = e

    if exit_code != 0:
        raise Bowtie2IndexError(stdout)

    return (local_index_dir, gsize)

@large_spot_task
def trim_align_callpeak(input: ATACseqInput) -> Optional[ATACseqOutput]:
    print(f"Processing {input.sample_name}")
    REMOTE_PATH = f"latch:///{input.base_remote_output_dir}{input.run_name}"
    BAM_REMOTE_PATH = f"{REMOTE_PATH}/Alignments (Bowtie2)/"
    MACS3_REMOTE_PATH = f"{REMOTE_PATH}/Peak-calling (MACS3)/"

    #
    # Read trimming
    # 
    try:
        outputs = [do_trimgalore(input, i, x) for i, x in enumerate(input.replicates)]
    except TrimgaloreError as e:
        print(f"Handling failure in trimming {input.sample_name} gracefully.")
        print(f"\tTrimming error ~ {e}")
        return

    trimmed_replicates = [x[0] for x in outputs]
    trimgalore_reports = [y for x in outputs for y in x[1]]
    fastqc_zip = [y for x in outputs for y in x[2]]
    fastqc_html = [y for x in outputs for y in x[3]]

    #
    # Merge replicates before alignment
    #
    print("Merging replicates")
    merged = _merge_replicates(trimmed_replicates, input.sample_name)
    flags = ["-U"] if len(merged) == 1 else ["-1", "-2"]
    reads = [x for (flag, path) in zip(flags, merged) for x in [flag, path]]

    # Free space.
    for rep in trimmed_replicates:
        os.remove(rep.r1.path)
        if isinstance(rep, PairedEndReads):
            os.remove(rep.r2.path)
    
    # 
    # Alignment with Bowtie2 piped with samtools to save space
    #
    bam_output = Path(f"bowtie2/{input.sample_name}.bam").resolve()
    bam_output.parent.mkdir(exist_ok=True)

    # Flag: -F 1804: remove reads where
    # - self or PE mate unmapped
    # - not primary alignment
    # - fails platform/vendor quality checks
    # - PCR/optical duplicate
    # Flag: -f 2: mapped in proper pair
    try:
        (exit_code, stdout) = piped_run([
            ["bowtie2",
             "--threads", "96",
             "-x", f"{input.genome_index.local_path}/ref"] + reads,
            ["samtools",
             "view", "-b", "-h",
             "-@", "96",
             "-q", "30",
             "-F", "1804",
             "-f", "2",
             "-o", bam_output]
        ])
    except subprocess.CalledProcessError as e:
        exit_code = 1
        stdout = e

    if exit_code > 0:
        raise Bowtie2AlignmentError(stdout)
    
    # 
    # Peak calling with MACS3
    #
    macs3_output = Path("MACS3").resolve()
    macs3_output.mkdir(exist_ok=True)
    
    return_code, stdout = piped_run([
        [
            "macs3", "callpeak",
            "--treatment", bam_output,
            "--name", input.sample_name,
            "--gsize", str(input.genome_size),
            "--bdg",
            "--qvalue", "0.01",
            "--outdir", "MACS3"
        ] + ["-f", "BAMPE"] * isinstance(input.replicates[0], PairedEndReads)
    ])

    if return_code != 0:
        raise Macs3Error(stdout)

    peak_file = Path(macs3_output, f"{input.sample_name}_summits.bed")
    bdg_file = Path(macs3_output, f"{input.sample_name}_treat_pileup.bdg")

    return ATACseqOutput(
        sample_name=input.sample_name,
        trimgalore_reports=trimgalore_reports,
        fastqc_html=fastqc_html,
        fastqc_zip=fastqc_zip,
        bam_file=LatchFile(str(bam_output), BAM_REMOTE_PATH + bam_output.name),
        bed_file=LatchFile(str(peak_file), MACS3_REMOTE_PATH + peak_file.name),
        bdg_file=LatchFile(str(bdg_file), MACS3_REMOTE_PATH + bdg_file.name),
    )

# @small_task
# def multiqc(
#     run_name: str,
#     fastqc_outputs: List[LatchFile],
#     output_directory: Optional[LatchDir],
# ) -> Optional[LatchFile]:

#     multiqc_report_file = None

#     REMOTE_PATH = f"latch:///{_remote_output_dir(output_directory)}{run_name}/"
#     """Remote path prefix for LatchFiles + LatchDirs"""

#     try:
#         subprocess.run(["multiqc", "."], check=True)
#         multiqc_report_file = LatchFile(
#             "/root/multiqc_report.html",
#             REMOTE_PATH + "multiqc_report.html",
#         )
#     except subprocess.CalledProcessError as e:
#         print(f"Error occurred while generating MultiQC report -> {e}")
#         message(
#             "error",
#             {
#                 "title": "Unable to generate MultiQC report",
#                 "body": "See logs for more information",
#             },
#         )

#     return multiqc_report_file

    
@workflow
def atacseq(
    samples: List[Sample],
    ref_selection_fork: str,
    output_location_fork: str,
    run_name: str,
    latch_genome: LatchGenome,
    custom_ref_genome: Optional[LatchFile] = None,
    custom_output_dir: Optional[LatchDir] = None,
):
    """Perform alignment and Peak calling on ATAC-Sequencing reads

    ATAC-Seq ()
    ----

    This workflow will produce [...]

    # Workflow Anatomy

    # Disclaimer

    This workflow assumes that your sequencing reads were derived from *short-read
    cDNA seqeuncing (as opposed to long-read cDNA/direct RNA sequencing). If in
    doubt, you can likely make the same assumption, as it is by far the most common
    form of "RNA-sequencing".

    # Brief Summary of ATAC-seq

    This workflow ingests short-read sequencing files (in FastQ format) that came
    from the following sequence of steps[^1]:

      - RNA extraction from sample
      - cDNA synthesis from extracted RNA
      - adaptor ligation / library prep
      - (likely) PCR amplification of library
      - sequencing of library

    You will likely end up with one or more FastQ files from this process that hold
    the sequencing reads in raw text form. This will be the starting point of our
    workflow.

    (If you have a `.bcl` file, this holds the raw output of a sequencing machine.
    There are there are [external
    tools](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
    that can convert these files to FastQ format, which you will need before you can
    proceed).

    # Quality Control

    As a pre-processing step, its important to check the quality of your sequencing
    files. FastQC is the industry staple for generating a report of useful summary
    statistics[^2] and is available if you double-click on a file on the [LatchBio
    platform](https://console.latch.bio).

    The following are the most useful of these statistics:

      - *Per base sequence quality* gives the per-site distribution over the length
    of the read
      - *Sequence duplication levels* reveals duplicated reads, indicating degraded
    RNA samples or aggressive PCR cycling[^1]

    For a full breakdown of the values and their interpretation, we refer the
    reader to this
    [tutorial](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html).

    # Trimming

    Short-read sequencing introduces adapters, small sequences attached to the 5'
    and 3' end of cDNA fragments, that are present as artifacts in our FastQ files
    and must be removed.

    We have yet to identify a comprehensive review of the various trimming tools, so
    we have selected [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
    trusted by researchers we work with out of UCSF and Stanford, until we are able
    to do so ourself.

    # Alignment

    Alignment is the process of assigning a sequencing read a location on a
    reference genome or transcriptome. It is the most computationally expensive step
    of the workflow, requiring a comparison against the entire reference sequence
    for each of millions of reads.

    Transcript alignment is conducted similarly to genomic alignment, using tools like
    [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to rigorously
    recover reference coordinates for each read.

    # Peak calling

    __metadata__:
        display_name: ATAC-seq
        wiki_url: https://www.latch.wiki/
        video_tutorial: 
        author:
            name: LatchBio
            email: help@latch.bio
            github: github.com/latchbio
        repository: github.com/latch-verified/atac-seq
        license:
            id: MIT
        flow:
        - section: Samples
          flow:
            - text: >-
                  Sample files can be provided and their read type can be
                  inferred from their name or this information can be specified manually.
                  Sample strandedness is inferred automatically (learn more).

            - params:
                - samples
        - section: Alignment
          flow:
            - text: >-
                  Alignment of trimmed reads on reference genome with Bowtie2
            - fork: ref_selection_fork
              flows:
                from_db:
                  display_name: Select from Latch Genome Database
                  flow:
                    - text: >-
                        We have curated a set of reference
                        genome data for ease and
                        reproducibility. More information about
                        these managed files can be found
                        [here](https://github.com/latchbio/latch-genomes).
                    - params:
                        - latch_genome
                custom:
                    display_name: Provide Custom Genome
                    flow:
                        - params:
                            - custom_ref_genome
        - section: Output Location
          flow:
          - params:
              - run_name
          - fork: output_location_fork
            flows:
                default:
                    display_name: Default
                    flow:
                    - text:
                        Output will be at default location in the data
                        viewer - ATAC-Seq Outputs/"Run Name"
                custom:
                    display_name: Specify Custom Path
                    _tmp_unwrap_optionals:
                        - custom_output_dir
                    flow:
                    - params:
                        - custom_output_dir
    Args:

        samples:
            Here you can organize your FastQ files by sample and add technical
            replicates for each sample.  Biological replicates should be
            organized as separate samples.

          __metadata__:
            display_name: Sample Sheet
            batch_table_column: true
            _tmp:
                custom_ingestion: auto

        latch_genome:
          Curated reference files for specific genome sources and builds.

          __metadata__:
            batch_table_column: true
            display_name: Genome Database Option

        ref_selection_fork:
          Select a reference genome from our curated database or provide your own.

          __metadata__:
            display_name: Reference Genome Source

        custom_ref_genome:
          The reference genome you want to align you samples to.

          __metadata__:
            display_name: Reference Genome File
            appearance:
                detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        run_name:
          A name for this analysis run, this will be used to name outputs from
          this run.

          __metadata__:
            batch_table_column: true
            display_name: Run Name

        output_location_fork:

        custom_output_dir:
          You can provide a custom location where this run's analysis outputs
          will be located.

          __metadata__:
            display_name: Custom Output Location
    """
    inputs = prepare_inputs(
        samples=samples,
        run_name=run_name,
        latch_genome=latch_genome,
        clip_r1=None,
        clip_r2=None,
        three_prime_clip_r1=None,
        three_prime_clip_r2=None,
        custom_output_dir=custom_output_dir,
        custom_ref_genome=custom_ref_genome,
    )

    map_task(trim_align_callpeak)(input=inputs)

    # multiqc_file = multiqc(
    #     run_name=run_name,
    #     fastqc_outputs=[f for output in outputs for f in output.fastqc
    #                     if output is not None],
    #     output_directory=custom_output_dir,
    # )
    # multiqc_report_file = multiqc(
    #     run_name=run_name,
    #     ts_outputs=outputs,
    #     output_directory=custom_output_dir,
    # )

LaunchPlan(
    atacseq,
    "Test Data - arp8delta vs wild-type",
    {
        "samples": [
            Sample(
                name="Control",
                strandedness=Strandedness.auto,
                replicates=[
                    PairedEndReads(
                        r1=LatchFile(
                            "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/WT_R1.fastq.gz",
                        ),
                        r2=LatchFile(
                            "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/WT_R2.fastq.gz",
                        )
                    ),
                ],
            ),
            Sample(
                name="Arp8delta",
                strandedness=Strandedness.auto,
                replicates=[
                    PairedEndReads(
                        r1=LatchFile(
                            "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/arp8delta_rep1_R1.fastq.gz",
                        ),
                        r2=LatchFile(
                            "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/arp8delta_rep1_R2.fastq.gz",
                        )
                    ),
                    PairedEndReads(
                        r1=LatchFile(
                            "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/arp8delta_rep2_R1.fastq.gz",
                        ),
                        r2=LatchFile(
                            "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/arp8delta_rep2_R2.fastq.gz",
                        )
                    ),
                ],
            ),
        ],
        "run_name": "ATAC-seq-test",
        "latch_genome": LatchGenome.RefSeq_R64,
    },
)

if __name__ == '__main__':
    
    samples = [
        # Sample(
        #     name="Control",
        #     strandedness=Strandedness.auto,
        #     replicates=[
        #         PairedEndReads(
        #             r1=LatchFile(
        #                 "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/WT_R1.fastq.gz",
        #             ),
        #             r2=LatchFile(
        #                 "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/WT_R2.fastq.gz",
        #             )
        #         ),
        #     ],
        # ),
        Sample(
            name="Arp8delta",
            strandedness=Strandedness.auto,
            replicates=[
                PairedEndReads(
                    r1=LatchFile(
                        "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/arp8delta_rep1_R1.fastq.gz",
                    ),
                    r2=LatchFile(
                        "s3://latch-public/test-data/7482/SCerevisae-wt-vs-arp8d/arp8delta_rep1_R2.fastq.gz",
                    )
                ),
            ],
        ),
    ]

    atacseq(samples=samples,
            ref_selection_fork="from_db",
            output_location_fork="default",
            run_name="ATAC-seq-test",
            latch_genome=LatchGenome.RefSeq_R64,
            custom_ref_genome=None,
            custom_output_dir=None)
