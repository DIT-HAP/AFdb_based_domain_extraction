def calculate_mean_pLDDT(domain, pLDDTs):
    pLDDT_each_residue = list(map(float, pLDDTs.split(",")))

    domain_length = 0
    domain_pLDDT = []

    for sub_domain in domain.split(","):
        start, end = sub_domain.split("-")
        # switch to 0-based index and left open right close
        start = int(start) - 1
        end = int(end)
        sub_domain_length = end - start
        domain_length += sub_domain_length
        domain_pLDDT += pLDDT_each_residue[start:end]

    mean_pLDDT = round(sum(domain_pLDDT) / domain_length, 2)
    frac_gt_70 = round(
        sum([1 for pLDDT in domain_pLDDT if pLDDT >= 70]) / domain_length, 2
    )

    return domain_length, mean_pLDDT, frac_gt_70


def merge_intervals(input, distance=1):
    # Convert the input string into a list of lists
    intervals = [list(map(int, x.split("-"))) for x in input]
    # Sort the list by the start of each interval
    intervals.sort(key=lambda x: x[0])
    # Initialize an empty list to store the merged intervals
    merged = []
    # Loop through the sorted list
    for interval in intervals:
        # If the merged list is empty or the current interval is separated by more than one unit from the last one in the merged list
        if not merged or interval[0] > merged[-1][1] + distance:
            # Append the current interval to the merged list
            merged.append(interval)
        else:
            # Otherwise, update the end of the last interval in the merged list to be the maximum of both ends
            merged[-1][1] = max(merged[-1][1], interval[1])
    # Return the merged list as a string
    return ["-".join(map(str, x)) for x in merged]


def check_overlap(interval1, interval_list):
    start1, end1 = interval1.split("-")
    start1 = int(start1)
    end1 = int(end1)
    for interval2 in interval_list:
        start2, end2 = interval2.split("-")
        start2 = int(start2)
        end2 = int(end2)
        if start1 <= start2 <= end2 <= end1:
            return True
        elif start2 <= start1 <= end1 <= end2:
            return True
        elif start1 <= start2 < end1 <= end2:
            return True
        elif start2 <= start1 < end2 <= end1:
            return True
    return False
