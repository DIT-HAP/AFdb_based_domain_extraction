"""
# Introduction: Methods for extracting confident domains from AFdb

### ["Real-time structure search and structure classification for AlphaFold protein models"](https://www.nature.com/articles/s42003-022-03261-8)

##### Method: Extraction of confident domain regions in AlphaFold2 models

- To extract a confident domain in an AlphaFold2 model, we first extracted all contiguous regions of more than 50 confident residues that have a pLDDT score greater than 70.0. 
- Then, confident regions separated by at most 5 non-confident residues were merged, along with the intervening residues regardless of confidence level. 
- AlphaFold2 models were discarded if they have no confident domains.

##### Results

In total, this procedure yielded 508,787 domains. 83,615 (22.9%) models out of 365,198 total AlphaFold2 models contain no confident domains. The statistics of model counts is provided in Table 1. In terms of total residues, the domain dataset in 3D-AF-Surfer contains 48.8% (78,133,986 residues) of residues among the residues in all the AlphaFold2 models (160,235,650 residues).
"""
# %%
import sys
from pathlib import Path
from itertools import groupby
import argparse
import numpy as np
import pandas as pd
from util.domain_related_functions import (
    calculate_mean_pLDDT,
    merge_intervals,
    check_overlap,
)

def main(args):

    pLDDT_scores = pd.read_csv(args.pLDDT_score, header=0, sep="\t")
    sys_id_to_uniprot = pd.read_csv(args.systematic_ID, sep="\t", header=0)

    confident_domains_disordered_domains = pd.DataFrame(columns=["uniprot_id", "domain_residues", "domain_type"])

    for _, row in pLDDT_scores.iterrows():

        uniprot_id = row["uniprot_id"]
        pLDDT = np.array(row["pLDDT"].split(",")).astype(float)
        protein_length = int(row["protein_length"])

        confident_domains = extract_confident_domains(pLDDT)
        merged_confident_domains = merge_intervals(confident_domains, distance=5+1)
        non_confident_domains = extract_non_confident_domains(protein_length, merged_confident_domains)

        # append each element in the list to a row in the dataframe, and then add the row to the dataframe
        for i in merged_confident_domains:
            new_row = pd.Series([uniprot_id, i, "confident"], index=confident_domains_disordered_domains.columns).to_frame().T
            confident_domains_disordered_domains = pd.concat([confident_domains_disordered_domains, new_row], ignore_index=True)
        for i in non_confident_domains:
            new_row = pd.Series([uniprot_id, i, np.nan], index=confident_domains_disordered_domains.columns).to_frame().T
            confident_domains_disordered_domains = pd.concat([confident_domains_disordered_domains, new_row], ignore_index=True)
    
    # sort the domain
    confident_domains_disordered_domains = (
        confident_domains_disordered_domains.groupby("uniprot_id")
        .apply(
            lambda sub_df: sub_df.sort_values(
                "domain_residues",
                key=lambda x: x.str.split("-", expand=True)[0].astype(int),
            ).reset_index(drop=False)[["domain_residues", "domain_type"]]
        )
        .rename_axis(["uniprot_id", "domain_id"])
        .reset_index()
    )
    confident_domains_disordered_domains["domain_id"] = confident_domains_disordered_domains["domain_id"] + 1

    # add the pLDDT column
    confident_domains_disordered_domains = confident_domains_disordered_domains.merge(pLDDT_scores, on="uniprot_id")

    # add the domain_length, mean_pLDDT, fraction_gt_70 columns
    confident_domains_disordered_domains[
        ["domain_length", "mean_pLDDT", "fraction_gt_70"]
    ] = confident_domains_disordered_domains.apply(
        lambda row: calculate_mean_pLDDT(row["domain_residues"], row["pLDDT"]),
        axis=1,
        result_type="expand",
    )

    confident_domains_disordered_domains_with_sysID = (
        sys_id_to_uniprot[["gene_systematic_id", "uniprot_id"]]
        .merge(confident_domains_disordered_domains, on="uniprot_id", how="right")
        .rename(columns={"gene_systematic_id": "Systematic ID"})
        .drop(columns="pLDDT")
    )
    
    # save file
    confident_domains_disordered_domains_with_sysID.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

def extract_confident_domains(pLDDT):
    pLDDT_ge_70 = np.where(pLDDT >= 70)[0]
    confident_domains = []
    for k, v in groupby(enumerate(pLDDT_ge_70), lambda x: x[0] - x[1]):
        # switch the 0-based index to 1 based coordinate
        residues_of_confident_domain = [y + 1 for x, y in list(v)]
        domain_start = residues_of_confident_domain[0]
        domain_end = residues_of_confident_domain[-1]
        domain_length = domain_end - domain_start + 1
        if domain_length >= 50:
            confident_domains.append(f"{domain_start}-{domain_end}")
    return confident_domains

def extract_non_confident_domains(protein_length, confident_domains):
    non_confident_domains = []
    start = 1
    for domain in confident_domains:
        domain_start, domain_end = map(int, domain.split("-"))
        # Add section before confident domain
        if start < domain_start:
            non_confident_domains.append(f"{start}-{domain_start-1}")
        # Update start to end of confident domain
        start = domain_end + 1
        # Add section after last confident domain  
    if start <= protein_length:
        non_confident_domains.append(f"{start}-{protein_length}")
    return non_confident_domains

if __name__ == "__main__":

    # add argument
    parser = argparse.ArgumentParser(description="Extract confident domains from AFdb")
    parser.add_argument(
        "-pLDDT",
        "--pLDDT-score",
        dest="pLDDT_score",
        type=Path,
        required=True,
        help="Path to the pLDDT score file",
    )
    parser.add_argument(
        "-sysID",
        "--systematic-ID",
        dest="systematic_ID",
        type=Path,
        required=True,
        help="TSV file storing the systematic ID and uniprot ID",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        required=True,
        help="Path to the output tsv file",
    )
    args = parser.parse_args()

    main(args)