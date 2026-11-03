import argparse
from pathlib import Path
import pandas as pd

def main():

    input = "/data/c/yangyusheng_optimized/AFdb_based_domain_extraction/results/UP000002485_284812_SCHPO_v4_merged_DNA_domain.tsv"

    merged_DNA_domains = pd.read_csv(input, sep="\t", header=0, usecols=list(range(10)), index_col=list(range(9)))
    
    merged_DNA_domains = merged_DNA_domains["Merged_Start_End"].str.split(
        ",", expand=True
        ).stack().droplevel(-1).str.split(
            "-", expand=True
        ).astype(int).rename(
            {
                0: "Start",
                1: "End"
            }, axis=1
        )

    merged_DNA_domains["pseudoScore"] = "."


    reordered_columns = [
        "#Chr",
        "Start_Region",
        "End_Region",
        "domain_residues",
        "score",
        "Strand",
        "Start",
        "End",
        "itemRgb"
    ]


    pass


if __name__ == '__main__':

    # parser = argparse.ArgumentParser(
    #     description="Create genome bed file for protein domains"
    # )

    # parser.add_argument(
    #     "-i",
    #     "--input",
    #     type=Path,
    #     help="Input file path",
    #     required=True
    # )

    # parser.add_argument(
    #     "-o",
    #     "--output",
    #     type=Path,
    #     help="Output file path",
    #     required=True
    # )

    # args = parser.parse_args()

    # main(args)
    main()