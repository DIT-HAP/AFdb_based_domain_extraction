# %%
from ast import arg, main
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from functools import reduce
from util.domain_related_functions import (
    calculate_mean_pLDDT,
    merge_intervals,
    check_overlap,
)

def main(args):
    extracted_domains_based_on_PAE = pd.read_csv(args.PAE_domain, header=0, sep="\t")
    pLDDT_scores = pd.read_csv(args.pLDDT_score, header=0, sep="\t")
    if ".tsv" in args.systematic_ID.name:
        sys_id_to_uniprot = pd.read_csv(args.systematic_ID, sep="\t")
    elif ".csv" in args.systematic_ID.name:
        sys_id_to_uniprot = pd.read_csv(args.systematic_ID, sep=",")
    elif ".xlsx" in args.systematic_ID.name:
        sys_id_to_uniprot = pd.read_excel(args.systematic_ID)
    else:
        raise ValueError("Systematic ID file must be a tsv, csv or xlsx file")

    domains_from_PAE_with_pLDDT = extracted_domains_based_on_PAE.merge(pLDDT_scores)
    domains_from_PAE_with_pLDDT[
        ["domain_length", "mean_pLDDT", "fraction_gt_70"]
    ] = domains_from_PAE_with_pLDDT.apply(
        lambda row: calculate_mean_pLDDT(row["domain_residues"], row["pLDDT"]),
        axis=1,
        result_type="expand",
    )

    reformat_domains = (
        domains_from_PAE_with_pLDDT[["uniprot_id", "version", "protein_length", "domain_residues", "mean_pLDDT"]]
        .groupby(["uniprot_id", "version", "protein_length"])
        .apply(merge_domain)
    )

    reformat_domains = (
        reformat_domains.explode()
        .apply(lambda x: ",".join(x) if isinstance(x, list) else x)
        .to_frame("domain_residues")
    )
    reformat_domains = (
        reformat_domains.groupby(["uniprot_id", "version", "protein_length"])
        .apply(
            lambda sub_df: sub_df.sort_values(
                "domain_residues",
                key=lambda x: x.str.split("-", expand=True)[0].astype(int),
            ).reset_index(drop=False)["domain_residues"]
        )
        .rename_axis(["uniprot_id", "version", "protein_length", "domain_id"])
        .reset_index()
    )
    reformat_domains["domain_id"] = reformat_domains["domain_id"] + 1
    reformat_domains = reformat_domains.merge(pLDDT_scores, how="left")
    reformat_domains[
        ["domain_length", "mean_pLDDT", "fraction_gt_70"]
    ] = reformat_domains.apply(
        lambda row: calculate_mean_pLDDT(row["domain_residues"], row["pLDDT"]),
        axis=1,
        result_type="expand",
    )

    reformated_domains_with_sysID = (
        sys_id_to_uniprot[["gene_systematic_id", "uniprot_id"]]
        .merge(reformat_domains, on="uniprot_id", how="right")
        .rename(columns={"gene_systematic_id": "Systematic ID"})
        .drop(columns="pLDDT")
    )

    reformated_domains_with_sysID["Systematic ID"] = reformated_domains_with_sysID.apply(
        lambda row: row["Systematic ID"] if isinstance(row["Systematic ID"], str) else row["uniprot_id"], axis=1
    )


    # save file
    reformated_domains_with_sysID.to_csv(
        args.output,
        sep="\t",
        index=False,
    )


# # 3.Merge the PAE domains based on the criteria:
# - Merge the continuous disordered domains (mean_pLDDT < 70, including contiguous domains with mean_pLDDT less than 70) first
# - Merge the sub-domains between which the distance is less than 10 residues for each contiguous domain
# - Remove the extra domains that overlap with the merged domains
def merge_domain(sub_df):
    has_contiguous_domain_or_not = sub_df["domain_residues"].str.contains(",") & (
        sub_df[sub_df["domain_residues"].str.contains(",")]["mean_pLDDT"] >= 70
    )

    if sub_df.shape[0] == 1:
        return sub_df["domain_residues"].tolist()
    elif has_contiguous_domain_or_not.sum() == 0:
        return merge_only_non_contiguous_domains(sub_df)
    elif has_contiguous_domain_or_not.sum() == 1:
        return merge_only_one_contiguous_domain(sub_df, has_contiguous_domain_or_not)
    elif has_contiguous_domain_or_not.sum() >= 2:
        return merge_multi_contiguous_domains(sub_df, has_contiguous_domain_or_not)
    else:
        return np.nan


def merge_only_non_contiguous_domains(sub_df):
    disordered_region = sub_df[sub_df["mean_pLDDT"] < 70].copy()
    non_disordered_region = sub_df[sub_df["mean_pLDDT"] >= 70].copy()

    if disordered_region.shape[0] == 0:
        merged_disordered_region = []
    else:
        disordered_intervals = reduce(
            lambda x, y: x + y,
            disordered_region["domain_residues"].str.split(",").tolist(),
        )
        merged_disordered_region = merge_intervals(disordered_intervals)

    return merged_disordered_region + non_disordered_region["domain_residues"].tolist()


def merge_only_one_contiguous_domain(sub_df, has_contiguous_domain_or_not):
    disordered_region = sub_df[sub_df["mean_pLDDT"] < 70].copy()
    non_disordered_region = sub_df[
        (sub_df["mean_pLDDT"] >= 70) & (~has_contiguous_domain_or_not)
    ].copy()
    contiguous_domain = sub_df[has_contiguous_domain_or_not].copy()

    if disordered_region.shape[0] == 0:
        merged_disordered_region = []
    else:
        disordered_intervals = reduce(
            lambda x, y: x + y,
            disordered_region["domain_residues"].str.split(",").tolist(),
        )
        merged_disordered_region = merge_intervals(disordered_intervals)

    contiguous_domain_intervals = (
        contiguous_domain["domain_residues"].values[0].split(",")
    )
    # merge the sub domains in the contiguous domain that are separated by less than 11 residues
    merged_contiguous_domain = merge_intervals(contiguous_domain_intervals, 11)

    # drop the regions in merged_disordered_region and non_disordered_region that are overlapped with the merged_contiguous_domain
    filtered_merged_disordered_region = [
        x
        for x in merged_disordered_region
        if not check_overlap(x, merged_contiguous_domain)
    ]
    filtered_non_disordered_region = [
        x
        for x in non_disordered_region["domain_residues"].tolist()
        if not check_overlap(x, merged_contiguous_domain)
    ]

    if len(merged_contiguous_domain) == 1:
        return (
            filtered_merged_disordered_region
            + filtered_non_disordered_region
            + merged_contiguous_domain
        )
    else:
        return (
            filtered_merged_disordered_region
            + filtered_non_disordered_region
            + [merged_contiguous_domain]
        )


def merge_multi_contiguous_domains(sub_df, has_contiguous_domain_or_not):
    disordered_region = sub_df[sub_df["mean_pLDDT"] < 70].copy()
    non_disordered_region = sub_df[
        (sub_df["mean_pLDDT"] >= 70) & (~has_contiguous_domain_or_not)
    ].copy()
    contiguous_domain = sub_df[has_contiguous_domain_or_not].copy()

    if disordered_region.shape[0] == 0:
        merged_disordered_region = []
    else:
        disordered_intervals = reduce(
            lambda x, y: x + y,
            disordered_region["domain_residues"].str.split(",").tolist(),
        )
        merged_disordered_region = merge_intervals(disordered_intervals)

    # merge the sub domains in the contiguous domain that are separated by less than 11 residues
    merged_contiguous_domain = [
        merge_intervals(x.split(","), 11)
        for x in contiguous_domain["domain_residues"].tolist()
    ]
    merged_contiguous_domain_in_one_list = reduce(
        lambda x, y: x + y, merged_contiguous_domain
    )
    merged_contiguous_domain_in_one_list_modified = merged_contiguous_domain_in_one_list
    overlapping_contiguous_domains = []
    for x in merged_contiguous_domain:
        merged_contiguous_domain_in_one_list_modified = [
            m for m in merged_contiguous_domain_in_one_list_modified if m not in x
        ]
        for y in x:
            if check_overlap(y, merged_contiguous_domain_in_one_list_modified):
                overlapping_contiguous_domains.append(y)
    overlapping_contiguous_domains = list(set(overlapping_contiguous_domains))

    # drop the regions in merged_disordered_region and non_disordered_region that are overlapped with the merged_contiguous_domain
    filtered_merged_disordered_region = [
        x
        for x in merged_disordered_region
        if not check_overlap(x, merged_contiguous_domain_in_one_list)
    ]
    filtered_non_disordered_region = [
        x
        for x in non_disordered_region["domain_residues"].tolist()
        if not check_overlap(x, merged_contiguous_domain_in_one_list)
    ]

    if not overlapping_contiguous_domains:
        return (
            filtered_merged_disordered_region
            + filtered_non_disordered_region
            + merged_contiguous_domain
        )
    else:
        return (
            merged_disordered_region
            + non_disordered_region["domain_residues"].tolist()
            + contiguous_domain["domain_residues"].tolist()
        )


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-pae",
        "--PAE-domain",
        dest="PAE_domain",
        type=Path,
        required=True,
        help="CSV file storing the PAE domains",
    )
    parser.add_argument(
        "-pLDDT",
        "--pLDDT-score",
        dest="pLDDT_score",
        type=Path,
        required=True,
        help="CSV file storing the pLDDT scores",
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
        help="Output tsv file for storing reformated PAE domains based on pLDDT",
    )

    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    main(args)
