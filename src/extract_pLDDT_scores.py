import sys
from pathlib import Path
from itertools import groupby
import argparse
import numpy as np
import pandas as pd
from util.extract_pLDDT_functions import *
from tqdm import tqdm


def main(pLDDT_dir, output_file):
    all_pdb_files = list((pLDDT_dir).glob("*.pdb.gz")) + list((pLDDT_dir).glob("*.cif"))
    tqdm_iterator = tqdm(all_pdb_files)
    pLDDTs = pd.DataFrame(columns=["uniprot_id", "protein_length","pLDDT"])

    for pdb_file in tqdm_iterator:
        uniprot_id = pdb_file.name.split("-F1-")[0].lstrip("AF-")
        # pLDDT = np.array(extract_pLDDT_pdb_gz(pdb_file))
        pLDDT = np.array(extract_pLDDT(pdb_file))
        length_protein = len(pLDDT)
        ipLDDT = pd.DataFrame.from_dict({
            "uniprot_id": [uniprot_id],
            "protein_length": [length_protein],
            "pLDDT": [",".join(pLDDT.astype(str))]
        })
        pLDDTs = pd.concat([pLDDTs, ipLDDT], ignore_index=True)


    pLDDTs.to_csv(output_file, index=False, sep="\t")


if __name__ == "__main__":
    # add arguments
    parser = argparse.ArgumentParser(
        description="Extract pLDDT scores from AF database"
    )

    parser.add_argument(
        "-p",
        "--pLDDT_dir",
        type=Path,
        required=True,
        help="Path to directory containing pLDDT files",
    )
    parser.add_argument(
        "-o", "--output_file", type=Path, required=True, help="Path to output tsv file"
    )

    args = parser.parse_args()

    main(args.pLDDT_dir, args.output_file)
