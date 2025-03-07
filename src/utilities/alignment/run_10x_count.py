#!/usr/bin/env python

# Example: TAXON=homo CELL_COUNT=3000 S3_DIR=s3://biohub-spyros/data/10X_data/CK_Healthy/ ./10x_count.py
import argparse
import os
import pathlib
import sys
import subprocess
import tarfile
import posixpath

from utilities.log_util import get_logger, log_command


import boto3


# reference genome bucket name for different regions
S3_REFERENCE = {"east": "czbiohub-reference-east", "west": "czbiohub-reference", "krasnow" : "czbiohub-reference-krasnow"}

# valid and deprecated reference genomes
reference_genomes = {
    "homo": "HG38-PLUS",
    "hg38-plus": "HG38-PLUS",
    "homo.gencode.v30.annotation.ERCC92": "homo.gencode.v30.annotation.ERCC92",
    "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1": "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1",
    #"homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1_cellranger-3.0": "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1_cellranger-3.0",
    "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1.GCA_009937905.1" : "homo.gencode.v30.annotation.ERCC92_and_sars.cov2.wa1.GCA_009937905.1",
    "mus": "MM10-PLUS",
    "mm10-plus": "MM10-PLUS",
    "mm10-1.2.0": "mm10-1.2.0",
    "mus-premrna": "mm10-1.2.0-premrna",
    "mm10-1.2.0-premrna": "mm10-1.2.0-premrna",
    "hg19-mm10-3.0.0": "hg19-mm10-3.0.0",
    "microcebus": "MicMur3-PLUS",
    "gencode.vM19": "gencode.vM19",
    "GRCh38_premrna": "GRCh38_premrna",
    "zebrafish-plus": "danio_rerio_plus_STAR2.6.1d",
    "botryllus": "botryllus",
    "Mmur3-cellranger-7": "Mmur3-cellranger-7"
}
deprecated = {
    "homo": "hg38-plus",
    "mus": "mm10-plus",
    "mus-premrna": "mm10-1.2.0-premrna",
}

# other helpful constants
CELLRANGER = "cellranger"
S3_RETRY = 5


def get_default_requirements():
    return argparse.Namespace(
        vcpus=64, memory=384000, storage=2000#, image="thsuanwu/cellranger_7.0"
    )


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = argparse.ArgumentParser(
        prog="run_10x_count.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run alignment jobs using 10x",
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes.keys()),
        help="Reference genome for the alignment run",
    )

    requiredNamed.add_argument(
        "--sample_prefix", required=False, help="Specify sample prefix"
    )

    requiredNamed.add_argument(
        "--s3_input_path", required=True, help="The folder with fastq.gz files to align"
    )

    requiredNamed.add_argument(
        "--s3_output_path",
        required=True,
        help="The folder to store the alignment results",
    )

    requiredNamed.add_argument(
        "--num_partitions",
        type=int,
        required=True,
        default=10,
        help="Number of groups to divide samples "
        "into for the alignment run. Enter 10 as the default "
        "value here since we don't divide a single sample",
    )

    requiredNamed.add_argument(
        "--partition_id",
        type=int,
        required=True,
        help="Index of sample group. Enter 0 as "
        "the default value here since we only have one sample",
    )

    requiredNamed.add_argument('--by_folder', action='store_true')

    # optional arguments
    parser.add_argument("--cell_count", type=int, default=3000)

    parser.add_argument(
        "--legacy",
        action="store_true",
        help="Use if 10x run was not demuxed locally (pre November 2019)",
    )

    parser.add_argument(
        "--region",
        default="krasnow",
        choices=("east", "west", "krasnow"),
        help=(
            "Region you're running jobs in."
            " Should match the location of"
            " the fastq.gz files"
        ),
    )

    parser.add_argument("--glacier", action="store_true")
    parser.add_argument("--root_dir", default="/mnt")

    return parser


def main(logger):
    """ Download reference genome, run alignment jobs, and upload results to S3.

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

    genome_base_dir = args.root_dir / "genome" / "cellranger"
    genome_base_dir.mkdir(parents=True)

    # check if the input genome and region are valid
    if args.taxon in reference_genomes:
        if args.taxon in deprecated:
            logger.warn(
                f"The name '{args.taxon}' will be removed in the future,"
                f" start using '{deprecated[args.taxon]}'"
            )

        genome_name = reference_genomes[args.taxon]
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    genome_dir = genome_base_dir / genome_name
    ref_genome_10x_file = f"cellranger/{genome_name}.tgz"

    #if args.region != "west" and genome_name not in ("HG38-PLUS", "MM10-PLUS"):
    #    raise ValueError(f"you must use --region west for {genome_name}")

    if args.region == "east":
        ref_genome_10x_file = f"ref-genome/{ref_genome_10x_file}"

    logger.info(
        f"""Run Info: partition {args.partition_id} out of {args.num_partitions}
                   genome_dir:\t{genome_dir}
         ref_genome_10x_file:\t{ref_genome_10x_file}
                        taxon:\t{args.taxon}
                s3_input_path:\t{args.s3_input_path}"""
    )

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

    # check the input folder for existing runs
    #sample_name = {
    #    os.path.basename(fn).rsplit("_", 4)[0] for fn in fastq_path.glob("*fastq.gz")
    #}
    #assert len(sample_name) == 1, "Should only have one sample name to process"
    #sample_name = sample_name.pop()

    # Run cellranger
    os.chdir(result_path)

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
