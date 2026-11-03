
import argparse
from pathlib import Path

import pandas as pd

from util.domain_related_functions import merge_intervals

def main(args):

    args.output.parent.mkdir(parents=True, exist_ok=True)

    DNA_nucleotides_mapped_with_domains = pd.read_csv(args.domain, sep="\t", header=0)

    DNA_nucleotides_mapped_with_domains["Start_End"] = DNA_nucleotides_mapped_with_domains.apply(
        lambda row: f"{row['Start']}-{row['End']}",
        axis=1
    )

    merged_nucleotides = DNA_nucleotides_mapped_with_domains.groupby(
                            [
                                "#Chr",
                                "Systematic ID",
                                "Strand",
                                "Transcript",
                                "Start_Region",
                                "End_Region",
                                "Length_Region",
                                "domain_id",
                                "domain_residues",
                                "domain_type"
                            ]
                        )["Start_End"].apply(list).to_frame()
    
    merged_nucleotides["Merged_Start_End"] = merged_nucleotides.apply(
        lambda row: ",".join(merge_intervals(row["Start_End"])),
        axis=1
    )

    merged_nucleotides["Merged_Length"] = merged_nucleotides.apply(
        lambda row: sum([int(i.split("-")[1]) - int(i.split("-")[0]) for i in row["Merged_Start_End"].split(",")]),
        axis=1
    )

    merged_nucleotides.drop(columns="Start_End", inplace=True)

    merged_nucleotides.reset_index(inplace=True)

    merged_nucleotides["domain_fraction"] = merged_nucleotides.apply(cal_region_ratio, axis=1)

    merged_nucleotides.to_csv(args.output, sep="\t", index=False, header=True)
    

def cal_region_ratio(row):
    # calculate the ratio of the region affected by the domain
    region_strand = row["Strand"]
    if region_strand == "+":
        region_start = int(row["Start_Region"])
        region_end = int(row["End_Region"])
    elif region_strand == "-":
        region_start = int(row["End_Region"])
        region_end = int(row["Start_Region"])

    domain_start_end = row["Merged_Start_End"].split(",")

    distance_ratio = []
    for start_end in domain_start_end:
        start = int(start_end.split("-")[0])
        end = int(start_end.split("-")[1])
        if region_strand == "+":
            start, end = start, end
        elif region_strand == "-":
            start, end = end, start
        
        start_distance_to_region_start = abs(start - region_start)
        end_distance_to_region_start = abs(end - region_start)

        frac_start = round(start_distance_to_region_start / row["Length_Region"],2)
        frac_end = round(end_distance_to_region_start / row["Length_Region"], 2)

        distance_ratio.append(f"{frac_start}-{frac_end}")
    
    return ",".join(distance_ratio)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Organize DNA domains'
    )
    parser.add_argument(
        '-d',
        '--domain',
        required=True,
        type=Path,
        help='DNA nucleotides mapped with protein domains'
    )
    parser.add_argument(
        '-o',
        '--output',
        dest="output",
        required=True,
        type=Path,
        help='Output tsv file'
    )

    args = parser.parse_args()

    main(args)