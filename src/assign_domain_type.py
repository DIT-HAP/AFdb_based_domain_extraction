import argparse
from pathlib import Path
import pandas as pd

def main(args):

    domains = pd.read_csv(args.input, header=0, sep="\t")

    domains[["location", "domain_type"]] = domains.apply(
        lambda row: assign_type_and_loc(row), axis=1, result_type="expand"
    )

    domains.to_csv(args.output, sep="\t", index=False)

def assign_type_and_loc(row):

    start_and_ends = [ 
        i.split("-") for i in row["domain_residues"].split(",")
        ]

    # use the min and max of the start and end positions to get the start and end of the domain
    start = min([int(i[0]) for i in start_and_ends])
    end = max([int(i[1]) for i in start_and_ends])

    domain_length = int(row["domain_length"])
    protein_length = int(row["protein_length"])

    location = "middle"
    disorder = "disordered"
    if domain_length == protein_length:
        location = "Whole"
        if row["mean_pLDDT"] > 70:
            disorder = "domain"
    elif start == 1 and end != protein_length:
        location = "N-terminal"
        if row["mean_pLDDT"] > 70:
            disorder = "domain"
    elif start != 1 and end == protein_length:
        location = "C-terminal"
        if row["mean_pLDDT"] > 70:
            disorder = "domain"
    elif start != 1 and end != protein_length:
        location = "middle"
        if row["mean_pLDDT"] > 70:
            disorder = "domain"
    elif start == 1 and end == protein_length:
        location = "middle"
        if row["mean_pLDDT"] > 70:
            disorder = "domain"

    return [location, disorder]

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        required=True,
        help="Input tsv file storing the PAE domains",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        required=True,
        help="Output tsv file for storing reformated PAE domains based on pLDDT",
    )

    args = parser.parse_args()
    main(args)