#!/usr/bin/env python
import argparse
import datetime
import os
import re
import subprocess
import time

import utilities.log_util as ut_log, get_logger, log_command
import utilities.s3_util as s3u
from utilities.alignment.run_star_and_htseq import reference_genomes, deprecated

import boto3
from boto3.s3.transfer import TransferConfig


CURR_MIN_VER = datetime.datetime(2018, 10, 1, tzinfo=datetime.timezone.utc)

# reference genome bucket name for different regions
S3_REFERENCE = {"east": "czbiohub-reference-east", "west": "czbiohub-reference", "krasnow" : "czbiohub-reference-krasnow"}

# valid and deprecated reference genomes
reference_genomes = {
    "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1": "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1",
}
deprecated = {
    "homo": "hg38-plus",
    "mus": "mm10-plus",
    "mus-premrna": "mm10-1.2.0-premrna",
}

# 10x barcodes for different versions
barcodes_10x = {
    "10x3v3": "3M-february-2018.txt",
    "10x3v2": "737K-august-2016.txt",
    "10x3v1": "737K-april-2014_rc.txt",
    "10x5v1": "737K-august-2016.txt",
    "10x5v2": "737K-august-2016.txt"
}
# other helpful constants
STAR = "STAR"
S3_RETRY = 5


def get_default_requirements():
    return argparse.Namespace(vcpus=64, memory=384000, storage=2000, image="starsolo-aligner")


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = argparse.ArgumentParser(
        prog="run_10x_velocyto.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run expression dynamics (RNA velocity) analysis on 10x fastq data using STAR solo",
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes.keys()),
        help="Reference genome for the velocyto run. Choose the same genome used in the alignment job.",
    )

    requiredNamed.add_argument(
        "--version_10x",
        required=True,
        choices=list(barcodes_10x.keys()),
        help="A barcode whitelist is the list of all known barcode sequences that have been included in the assay kit and are available during library preparation, corresponding to 10x technolgoy used.",
    )

    requiredNamed.add_argument(
        "--sample_prefix", required=False, help="Specify sample prefix"
    )

    requiredNamed.add_argument(
        "--s3_input_path",
        required=True,
        help="The folder with fastq.gz files to align",
    )

    requiredNamed.add_argument(
        "--s3_output_path", required=True, help="Location for output",
    )

    requiredNamed.add_argument(
        "--num_partitions",
        type=int,
        required=True,
        default=10,
        help="Number of velocyto jobs to launch on the STAR " "alignment outputs",
    )

    requiredNamed.add_argument(
        "--partition_id", type=int, required=True, help="Index of velocyto job group",
    )
    
    requiredNamed.add_argument('--by_folder', action='store_true')

    requiredNamed.add_argument( # what does this do?
        "--input_dirs",
        nargs="+",
        required=True,
        help="List of input folders to process",
    )

    # optional arguments
    parser.add_argument(
        "--force_redo",
        action="store_true",
        help="Process files even when results already exist",
    )

    parser.add_argument("--glacier", action="store_true")
    parser.add_argument("--root_dir", default="/mnt")
    return parser

def main(logger):
    """ Download reference genome, run velocyto jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()

    args = parser.parse_args()

    args.root_dir = pathlib.Path(args.root_dir)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        args.root_dir = args.root_dir / os.environ["AWS_BATCH_JOB_ID"]

    # local directories
    if args.s3_input_path.endswith("/"):
        args.s3_input_path = args.s3_input_path[:-1]

    sample_id = os.path.basename(args.s3_input_path)
    result_path = args.root_dir / "data" / sample_id
    if args.legacy:
        fastq_path = result_path / "fastqs"
    else:
        fastq_path = result_path
    fastq_path.mkdir(parents=True)

    genome_base_dir = args.root_dir / "genome" / "STAR-2.7.9a"
    genome_base_dir.mkdir(parents=True)

    barcode_base_dir = args.root_dir / "barcodes"
    barcode_base_dir.mkdir(parents=True)

    # check if the input genome and region are valid
    if args.taxon in reference_genomes:
        genome_name = reference_genomes[args.taxon]
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    # check if the input genome and region are valid
    if args.version_10x in barcodes_10x:

        barcode_name = barcodes_10x[args.version_10x]
    else:
        raise ValueError(f"unknown 10x version {args.version_10x}")

    genome_dir = genome_base_dir / genome_name
    ref_genome_10x_file = f"STAR-2.7.9a/{genome_name}.tgz"

    barcode_dir = barcode_base_dir / version_10x
    barcode_10x_file = f"STAR-2.7.9a/{barcode_name}"


    logger.info(
        f"""Run Info: partition {args.partition_id} out of {args.num_partitions}
                   genome_dir:\t{genome_dir}
         ref_genome_10x_file:\t{ref_genome_10x_file}
                        taxon:\t{args.taxon}
                        10x version:\t{args.version_10x}
                s3_input_path:\t{args.s3_input_path}"""
    )

    sys.exit(0)


    s3 = boto3.resource("s3")

    # download the reference genome data
    logger.info(f"Downloading and extracting genome data {genome_name}")

    s3_genome_object = s3.Object(S3_REFERENCE[args.region], ref_genome_10x_file)

    with tarfile.open(fileobj=s3_genome_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=genome_base_dir)

    sys.stdout.flush()

    # download the fastq files
    command = [
        "aws",
        "s3",
        "cp",
        "--no-progress",
        "--recursive",
        "--exclude",
        "'*'",
        "--include",
        f"'*{args.sample_prefix}*'",
        "--force-glacier-transfer" if args.glacier else "",
        args.s3_input_path,
        f"{fastq_path}",
    ]
    log_command(logger, command, shell=True)

    logger.info(f"Running partition {args.partition_id} of {args.num_partitions}")

    # Set variables
    # Read in R1 files
    R1_files = glob.glob(fastq_path + '/*R[1]*.fastq.gz')
    R1_files.sort()

    R1_list = ",".join(R1_files)

    # Read in R2 files, sort to ensure order
    R2_files = glob.glob(fastq_path + '/*R[2]*.fastq.gz')
    R2_files.sort()
    
    R2_list = ",".join(R2_files)

    # Run cellranger
    os.chdir(result_path)

    sys.exit(0)



    if args.by_folder:
        
        command = [
        CELLRANGER,
        "count",
        "--localmem=240", # By default, will use 90% of mem
        "--nosecondary",
        "--disable-ui",
        f"--expect-cells={args.cell_count}",
        f"--id={args.sample_prefix}",
        f"--fastqs={fastq_path}",
        f"--transcriptome={genome_dir}" # no sample_prefix: run all samples in folder
        ]

    else:

        command = [
            CELLRANGER,
            "count",
            "--localmem=240", # By default, will use 90% of mem
            "--nosecondary",
            "--disable-ui",
            f"--expect-cells={args.cell_count}",
            f"--id={args.sample_prefix}",
            f"--fastqs={fastq_path}",
            f"--transcriptome={genome_dir}",
            f"--sample={args.sample_prefix}",
        ]

    failed = log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    if failed:
        raise RuntimeError("cellranger count failed")

    # Move outs folder to S3
    command = [
        "aws",
        "s3",
        "sync",
        "--no-progress",
        os.path.join(result_path, args.sample_prefix, "outs"),
        posixpath.join(args.s3_output_path, args.sample_prefix),
    ]
    for i in range(S3_RETRY):
        if not log_command(logger, command, shell=True):
            break
        logger.info(f"retrying sync")
    else:
        raise RuntimeError(f"couldn't sync output")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)

    main(mainlogger)