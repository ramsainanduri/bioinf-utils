#!/usr/bin/env python
"""
Script to extract and convert specific features from a GTF file into a BED file.

Author: Ram Nanduri
Version: 1.0.0
License: MIT
Contact: sai.rama7412@gmail.com or Ram.Nanduri@skane.se
Date: 2025-01-20
GitHub: https://github.com/ramnanduri/gtf-to-bed

Description:
    This script provides a flexible way to convert GTF annotations to BED format
    with options to filter and format the output. Users can extract specific
    features (gene, transcript, exon) and apply advanced filters based on
    regex patterns, feature types, and custom fields.

Usage:
    python gtf2bed_converter.py --gtf <input_gtf_file> --output <output_file> --feature <feature_type>
                        [--filter_key <key>] [--include_patterns <pattern> ...]
                        [--exclude_patterns <pattern> ...] [--bed_columns <simple|detailed>]
                        [--name_key <key>]

Features:
    - Extract genes, transcripts, or exons with custom filters.
    - Supports flexible filtering using regex patterns.
    - Outputs BED format in simple (4 columns) or detailed (6 columns).
    - Allows customization of the name field in the BED output.

Dependencies:
    - Python 3.6 or higher

Example:
    python gtf2bed_converter.py --gtf gencode.v33.annotation.gtf --output output.bed --feature gene --filter_key gene_type --include_patterns protein_coding --bed_columns detailed --name_key hgnc_id
"""


import argparse
import re
from typing import List, Optional


def gtf_to_bed(
    input_gtf: str,
    output_file: str,
    feature_type: str,
    filter_key: Optional[str] = None,
    include_patterns: Optional[List[str]] = None,
    exclude_patterns: Optional[List[str]] = None,
    bed_columns: Optional[str] = "simple",
    name_key: Optional[str] = None,
) -> None:
    """
    Extracts information from a GTF (Gene Transfer Format) file and converts it into a BED (Browser Extensible Data) file format.

    This function reads a GTF file, filters the entries based on the specified feature type and optional filtering criteria,
    and writes the relevant information to a BED file. The BED file can be in a simple format with 4 columns or a detailed format with 6 columns.

    -----------
    input_gtf : str
        Path to the input GTF file.
    output_file : str
        Path to the output BED file.
    feature_type : str
        The type of feature to extract (e.g., gene, transcript, exon).
    filter_key : str, optional
        Key to filter on (e.g., gene_type, transcript_type). Default is gene_type.
    include_patterns : list or str, optional
        Patterns to include based on filter_key (supports regex). Default is None.
    exclude_patterns : list or str, optional
        Patterns to exclude based on filter_key (supports regex). Default is None.
    bed_columns : str, optional
        Format of the BED file, "simple" (4 columns) or "detailed" (6 columns). Default is "simple".
    name_key : str, optional
        The key to use for the name field in the BED file. Default is None.

    Notes:
    ------
    - The GTF file is expected to have at least 9 tab-separated columns.
    - The function skips lines starting with '#' and lines with fewer than 9 columns.
    - The start coordinate in the GTF file is 1-based, while the BED file uses 0-based start coordinates.
    - If `name_key` is not provided, the function uses default naming conventions based on the feature type.
    - Filtering is applied based on the `filter_key`, `include_patterns`, and `exclude_patterns` parameters.
    - The `bed_columns` parameter determines whether the output BED file includes strand and type information.

    """
    with open(input_gtf, "r") as gtf, open(output_file, "w") as bed:
        for line in gtf:
            # Skip comments or metadata lines
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = (
                fields
            )

            # Filter based on user-selected feature type
            if feature != feature_type:
                continue

            # Parse attributes into a dictionary
            attr_dict = {}
            for attr in attributes.split(";"):
                if attr.strip():
                    key, value = attr.strip().split(" ", 1)
                    attr_dict[key] = value.strip('"')

            # Apply filtering based on include_patterns and exclude_patterns
            if filter_key:
                value = attr_dict.get(filter_key)
                if include_patterns:
                    if isinstance(include_patterns, str):
                        include_patterns = [include_patterns]
                    if not any(
                        re.search(pattern, value) for pattern in include_patterns
                    ):
                        continue
                if exclude_patterns:
                    if isinstance(exclude_patterns, str):
                        exclude_patterns = [exclude_patterns]
                    if any(re.search(pattern, value) for pattern in exclude_patterns):
                        continue

            # Extract desired BED columns
            start = int(start) - 1  # Convert GTF 1-based start to 0-based BED start

            # Determine the name to use in the BED file
            if name_key:
                name = attr_dict.get(name_key, feature_type)
            elif feature_type == "gene":
                name = attr_dict.get("gene_name", "gene")
            elif feature_type == "transcript":
                name = attr_dict.get("transcript_name", "transcript")
            elif feature_type == "exon":
                _gene = attr_dict.get("gene_name", "gene")
                _exon = attr_dict.get("exon_number", "exon")
                name = f"{_gene};Exon{_exon}"
            else:
                name = feature_type

            # Optionally include strand and type in BED file
            if bed_columns == "detailed":
                type_info = attr_dict.get(filter_key, "NA")
                bed.write(f"{chrom}\t{start}\t{end}\t{strand}\t{name}\t{type_info}\n")
            else:
                bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract information from GTF file to BED file based on feature type and optional filtering."
    )
    parser.add_argument("--gtf", required=True, help="Path to the input GTF file.")
    parser.add_argument("--output", required=True, help="Path to the output BED file.")
    parser.add_argument(
        "--feature",
        required=True,
        choices=["gene", "transcript", "exon"],
        help="Feature type to extract (e.g., gene, transcript, exon).",
    )
    parser.add_argument(
        "--filter_key",
        help="Key to filter on (e.g., gene_type, transcript_type).",
        default="gene_type",
    )
    parser.add_argument(
        "--include_patterns",
        nargs="+",
        help="Patterns or regex to include based on the filter_key.",
    )
    parser.add_argument(
        "--exclude_patterns",
        nargs="+",
        help="Patterns or regex to exclude based on the filter_key.",
    )
    parser.add_argument(
        "--bed_columns",
        choices=["simple", "detailed"],
        default="simple",
        help="Format of the BED file: 'simple' (4 columns) or 'detailed' (6 columns). Default is 'simple'.",
    )
    parser.add_argument(
        "--name_key",
        help="Key to use for the name field in the BED file (e.g., gene_name, transcript_id).",
        default="gene_name",
    )

    args = parser.parse_args()

    gtf_to_bed(
        args.gtf,
        args.output,
        args.feature,
        filter_key=args.filter_key,
        include_patterns=args.include_patterns,
        exclude_patterns=args.exclude_patterns,
        bed_columns=args.bed_columns,
        name_key=args.name_key,
    )
