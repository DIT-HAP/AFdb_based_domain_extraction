#!/bin/bash

# Path to your script directory to make output relative path work as intended
# cd /data/c/yangyusheng_optimized/AFdb_based_domain_extraction/src/

POMBASE_FILE="/data/c/yangyusheng_optimized/resource/Pombase_FTP/release/20250501/gene_IDs_names_products.tsv"
TED_FILE="/data/c/yangyusheng_optimized/AFdb_based_domain_extraction/references/TED/ted_365m.domain_summary.cath.globularity.taxid.tsv"

# Adjust output path if running from a different PWD or use absolute paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
OUTPUT_DIR_AWK="${SCRIPT_DIR}/../results/" # Relative to AFdb_based_domain_extraction/src/
mkdir -p "$OUTPUT_DIR_AWK"
OUTPUT_FILE_AWK="${OUTPUT_DIR_AWK}pombe_ted_domains_detailed.tsv"
TEMP_UNIPROT_IDS="${SCRIPT_DIR}/uniprot_ids_to_keep.tmp" # Place temp file in script dir

echo "Extracting UniProt IDs for awk..."
# Assuming uniprot_id is the 6th column in the Pombase TSV file
awk -F'\t' 'NR > 1 && $6 != "" {print $6}' "$POMBASE_FILE" > "$TEMP_UNIPROT_IDS"

total_target_ids=$(wc -l < "$TEMP_UNIPROT_IDS" | tr -d ' ') # Get the count of IDs to find

if [[ "$total_target_ids" -eq 0 ]]; then
    echo "Warning: No UniProt IDs found in $POMBASE_FILE to filter by. Output file will be empty."
    touch "$OUTPUT_FILE_AWK" # Create empty output file
    rm "$TEMP_UNIPROT_IDS"
    exit 0
fi

echo "Found $total_target_ids UniProt IDs to filter by."
echo "Processing $TED_FILE with awk... Output to $OUTPUT_FILE_AWK"

awk -v total_n_ids="$total_target_ids" -F'\t' '
    BEGIN {
        # Load IDs to keep into an associative array
        while ((getline id < "'"$TEMP_UNIPROT_IDS"'") > 0) {
            ids_to_keep[id] = 1
        }
        close("'"$TEMP_UNIPROT_IDS"'");

        found_target_ids_count = 0;
        update_interval = 100000; # Print progress every 100,000 lines from TED_FILE
        
        # For formatting lines processed
        if (NR > 0) { # Check if NR is available, though in BEGIN it is 0
             lines_processed_display = sprintf("%.1fM", NR / 1000000.0);
        } else {
             lines_processed_display = "0.0M";
        }
        
        # Initial progress message
        percentage = 0;
        if (total_n_ids > 0) {
            percentage = (found_target_ids_count / total_n_ids) * 100;
        }
        printf("\rProgress: %d/%d unique IDs found (%.2f%%). Processed: %s lines", found_target_ids_count, total_n_ids, percentage, lines_processed_display) > "/dev/stderr";

    }
    {
        # $1 is the first field. Match AF-(UNIPROTID)-F pattern.
        split($1, ted_id_parts, "-");
        if (length(ted_id_parts) >= 3 && ted_id_parts[1] == "AF" && ted_id_parts[3] ~ /^F/) {
            ted_uniprot_id = ted_id_parts[2];
            if (ted_uniprot_id in ids_to_keep) {
                print $0; # Print the whole line to stdout (goes to OUTPUT_FILE_AWK)

                # If this target ID hasn_t been counted yet for progress
                if (!(ted_uniprot_id in already_counted_for_progress)) {
                    already_counted_for_progress[ted_uniprot_id] = 1;
                    found_target_ids_count++;
                }
            }
        }

        # Progress update
        if (NR % update_interval == 0) {
            percentage = 0;
            if (total_n_ids > 0) {
                percentage = (found_target_ids_count / total_n_ids) * 100;
            }
            lines_processed_display = sprintf("%.1fM", NR / 1000000.0);
            printf("\rProgress: %d/%d unique IDs found (%.2f%%). Processed: %s lines", found_target_ids_count, total_n_ids, percentage, lines_processed_display) > "/dev/stderr";
        }
    }
    END {
        percentage = 0;
        if (total_n_ids > 0) {
            percentage = (found_target_ids_count / total_n_ids) * 100;
        }
        lines_processed_display = sprintf("%.1fM", NR / 1000000.0);
        # Clear the line with \r and spaces then print final message
        printf("\r%-100s\r", " ") > "/dev/stderr"; # Clear the line
        printf("Final Progress: %d/%d unique IDs found (%.2f%%). Total Processed: %s lines from TED file.\n", found_target_ids_count, total_n_ids, percentage, lines_processed_display) > "/dev/stderr";
    }
' "$TED_FILE" > "$OUTPUT_FILE_AWK"

rm "$TEMP_UNIPROT_IDS"
echo "Awk processing complete. Output at $OUTPUT_FILE_AWK"