#!/usr/bin/env python3

import io
import glob
import os

import setuptools


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ).read()


setuptools.setup(
    name="czb-util",
    version="0.4.0",
    license="MIT License",
    description="A collection of scripts for some common Biohub tasks",
    long_description=read("README.md"),
    author="James Webber",
    author_email="james.webber@czbiohub.org",
    url="https://github.com/thsuanwu/utilities",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("src/*.py")
    ],
    include_package_data=True,
    zip_safe=False,
    install_requires=["boto3 >= 1.7.41"],
    extras_require={
        "evros": ["aegea >= 3.6", "awscli >= 1.15.41", "awscli-cwlogs >= 1.4.4"],
        "h5ad": ["anndata"],
    },
    entry_points={
        "console_scripts": [
            "aws_star = utilities.scripts.aws_star:main [evros]",
            "aws_10x = utilities.scripts.aws_10x:main [evros]",
            "aws_velocyto = utilities.scripts.aws_velocyto:main [evros]",
            "batch_samplesheet = utilities.scripts.batch_samplesheet:main",
            "evros = utilities.scripts.evros:main [evros]",
            "frython = utilities.scripts.frython:main",
            "gene_cell_table = utilities.scripts.gene_cell_table:main [evros]",
            "starfails = utilities.scripts.starfails:main [evros]",
        ]
    },
    scripts=glob.glob("scripts/*"),
)
