import argparse
from pathlib import Path
from altair import Type
import pandas as pd
import numpy as np
from multiprocessing import Pool

def main(args):

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cores = int(args.cores)

    # read domain file
    domain_df = pd.read_csv(args.domain, sep="\t", header=0)

    # read nucleotide bed file
    nucleotide_df = pd.read_csv(args.nucleotide_bed, sep="\t", header=0)

    # calculate the affected residues
    nucleotide_df["Residue_affected"] = nucleotide_df["Accumulated_CDS_bases"] // 3 + 1
    nucleotide_df["Residue_frame"] = nucleotide_df["Accumulated_CDS_bases"] % 3

    nucleotide_df[["domain_id", "domain_residues", "domain_type"]] = nucleotide_df.apply(
        lambda row: assign_protein_domain(row, domain_df), axis=1, result_type="expand"
    )

    nucleotide_df.to_csv(
        args.output,
        sep="\t",
        header=True,
        index=False
    )

    # multi-processing test

    # chunk_size = int(len(nucleotide_df)/cores)
    # chunks = np.array_split(nucleotide_df, len(nucleotide_df) // chunk_size + 1)

    # pool = Pool(cores)
    # chunked_nucleotide_df = pool.map(assign_domain_parallelly, chunks, domain_df)
    # pool.close()

    # all_nucleotide_df = pd.concat(chunked_nucleotide_df, ignore_index=True).sort_values(["Systematic ID", "domain_residues"])

    # all_nucleotide_df.to_csv(
    #     # args.output.parent / "all_nucleotide_df.tsv",
    #     args.output,
    #     sep="\t",
    #     header=True,
    #     index=False
    # )


def assign_protein_domain(row, domain):
    Gene = row["Systematic ID"]
    if Gene in domain["Systematic ID"].values:
        residue = int(row["Residue_affected"])
        domains = domain[domain["Systematic ID"] == Gene].copy()
        for idx, row in domains.iterrows():
            for iDomain in row.loc["domain_residues"].split(","):
                domain_start = int(iDomain.split("-")[0])
                domain_end = int(iDomain.split("-")[1])
                domain_id = row.loc["domain_id"]

                if domain_start <= residue <= domain_end:
                    return domain_id, row.loc["domain_residues"], row.loc["domain_type"]
        return np.nan, np.nan, np.nan
    else:
        return np.nan, np.nan, np.nan

def assign_domain_parallelly(chunk, domain_df):
    chunk[["domain_id", "domain_residues", "domain_type"]] = chunk.apply(
        lambda row: assign_protein_domain(row, domain_df), axis=1, result_type="expand"
    )
    return chunk


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Map DNA sequence with protein domain'
    )
    parser.add_argument(
        '-d',
        '--domain',
        required=True,
        type=Path,
        help='Protein domain tsv file'
    )
    parser.add_argument(
        '-n',
        '--nucleotide-bed',
        dest="nucleotide_bed",
        required=True,
        type=Path,
        help='Bed file for each nucleotide of the coding region in the genome'
    )
    parser.add_argument(
        '-c',
        '--cores',
        dest="cores",
        required=True,
        type=int,
        help='Number of cores'
    )
    parser.add_argument(
        '-o',
        '--output',
        dest="output",
        required=True,
        type=Path,
        help='Output file'
    )

    args = parser.parse_args()

    main(args)