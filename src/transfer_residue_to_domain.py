from itertools import groupby
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def main(input_dir, output_file):
    domain_files = input_dir.glob("*.domain")

    protein_domains = pd.DataFrame(columns=["uniprot_id", "version", "domain_id", "domain_residues", "protein_length"])
    for domain_file in domain_files:
        uniprot_id = domain_file.name.split("-F1-")[0].lstrip("AF-")
        version = domain_file.name.rstrip(".domain").split("_")[-1]
        domains = transfer_residue_to_domain_per_file(domain_file)
        protein_length = max([int(domain.split("-")[-1]) for domain in domains])
        for domain_id, domain in enumerate(domains):
            ipd = pd.DataFrame.from_dict({
                "uniprot_id": [uniprot_id],
                "version": [version],
                "domain_id": [domain_id+1],
                "domain_residues": [domain],
                "protein_length": [protein_length]
            })
            protein_domains = pd.concat([protein_domains, ipd], ignore_index=True)

    protein_domains.sort_values(by=["uniprot_id", "version", "domain_id"], inplace=True)

    protein_domains.to_csv(output_file, index=False, header=True, sep="\t")


def transfer_residue_to_domain_per_file(file):
    with open(file, "r") as f:
        lines = f.readlines()

        domains = []
        for line in lines:
            residues = np.array(line.strip().split(","))
            domain_residues = residues[residues != ""].astype(int)

            # extract the sub domain
            sub_domains = []
            for k, v in groupby(enumerate(domain_residues), lambda x: x[1] - x[0]):
                # switch to 1-based index
                consecutive_residues = [y + 1 for x, y in list(v)]
                sub_domain = f"{consecutive_residues[0]}-{consecutive_residues[-1]}"
                sub_domains.append(sub_domain)

            domains.append(sub_domains)

    domains = list(map(lambda x: ",".join(x), domains))

    return domains


if __name__ == "__main__":
    # add arguments
    parser = argparse.ArgumentParser(description="Transfer residue to domain")
    parser.add_argument(
        "-i", "--input_dir", type=Path, required=True, help="input directory"
    )
    parser.add_argument(
        "-o", "--output-file", type=Path, required=True, help="output tsv file"
    )
    args = parser.parse_args()

    main(args.input_dir, args.output_file)
