import csv
import re
import os

def extract_and_filter_ted_domains(pombase_file_path, ted_file_path, output_file_path):
    uniprot_ids_to_keep = set()
    
    print(f"Reading UniProt IDs from {pombase_file_path}...")
    try:
        with open(pombase_file_path, 'r', newline='') as p_file:
            reader = csv.reader(p_file, delimiter='\t')
            header = next(reader) # Skip header
            try:
                uniprot_id_col_index = header.index("uniprot_id")
            except ValueError:
                print(f"Error: 'uniprot_id' column not found in the header of the pombase file: {header}")
                return

            for row in reader:
                if len(row) > uniprot_id_col_index and row[uniprot_id_col_index]:
                    uniprot_ids_to_keep.add(row[uniprot_id_col_index])
        print(f"Found {len(uniprot_ids_to_keep)} unique UniProt IDs to filter by.")
    except FileNotFoundError:
        print(f"Error: Pombase file not found at {pombase_file_path}")
        return
    except Exception as e:
        print(f"An error occurred while reading {pombase_file_path}: {e}")
        return

    if not uniprot_ids_to_keep:
        print("No UniProt IDs found in the pombase file. Exiting.")
        return

    print(f"Processing {ted_file_path} and writing matches to {output_file_path}...")
    lines_written = 0
    
    # Compile the regex for extracting UniProt ID from the first field of TED file lines.
    # This pattern looks for "AF-" followed by one or more alphanumeric characters (the UniProt ID),
    # followed by "-F". This is more specific than searching the whole line with ".*".
    ted_uniprot_extractor = re.compile(r'AF-([A-Z0-9]+)-F')

    try:
        # Ensure output directory exists
        output_dir = os.path.dirname(output_file_path)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                print(f"Created output directory: {output_dir}")
            except OSError as e:
                print(f"Error creating output directory {output_dir}: {e}")
                return


        with open(ted_file_path, 'r', newline='') as ted_in_file, \
             open(output_file_path, 'w', newline='') as ted_out_file:
            
            # Iterate through all lines in the TED file.
            # If the TED file has a header row that you wish to preserve unconditionally,
            # you would read it before this loop and write it to ted_out_file.
            # This script assumes all lines are data to be filtered, or that a header
            # would not match the UniProt ID pattern anyway.

            for line_number, line in enumerate(ted_in_file, start=1):
                stripped_line = line.strip()
                if not stripped_line: # Skip empty or whitespace-only lines
                    continue

                try:
                    # Extract the first field (column) by splitting the line by tab.
                    # We only split once to get the first field, which is more efficient
                    # than splitting the entire line if it has many columns.
                    first_field = stripped_line.split('\t', 1)[0]
                    
                    match = ted_uniprot_extractor.search(first_field)
                    
                    if match:
                        ted_uniprot_id = match.group(1)
                        if ted_uniprot_id in uniprot_ids_to_keep:
                            ted_out_file.write(line) # Write the original line (with newline)
                            lines_written += 1
                
                except AttributeError: 
                    # This occurs if match is None, meaning the regex didn't find a pattern.
                    # This is expected for lines not matching the format, so we can just continue.
                    # print(f"Debug: Line {line_number} did not match UniProt pattern in first field: '{first_field[:50]}...'")
                    continue 
                except IndexError:
                    # This could happen if a line is unexpectedly empty after stripping, though covered by `if not stripped_line`.
                    # Or if split somehow results in an empty list (very unlikely).
                    print(f"Warning: Line {line_number} in {ted_file_path} might be malformed. Line: '{stripped_line[:100]}...'. Skipping.")
                    continue
                except Exception as e:
                    # Catch any other unexpected error during the processing of a single line.
                    print(f"An error occurred processing line {line_number} ('{stripped_line[:100]}...') in {ted_file_path}: {e}. Skipping.")
                    continue

        print(f"\nExtraction complete.")
        print(f"{lines_written} matching domain entries written to {output_file_path}")
        if lines_written == 0 and len(uniprot_ids_to_keep) > 0:
            print("No matching UniProt IDs were found in the TED file.")

    except FileNotFoundError:
        print(f"Error: TED file not found at {ted_file_path}")
    except Exception as e:
        print(f"An critical error occurred: {e}")

# --- Main execution ---
if __name__ == "__main__":
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Define file paths relative to the workspace root or using absolute paths
    # Assuming workspace root is /data/c/yangyusheng_optimized based on previous context
    workspace_root = "/data/c/yangyusheng_optimized" 
    
    pombase_file = os.path.join(workspace_root, "resource/Pombase_FTP/release/20250501/gene_IDs_names_products.tsv")
    ted_domains_file = os.path.join(workspace_root, "AFdb_based_domain_extraction/references/TED/ted_214m_per_chain_segmentation.tsv")
    
    # Output file path from your script: "../results/pombe_ted_domains.tsv"
    # This will be relative to script_dir: AFdb_based_domain_extraction/src/
    # So, output path becomes: AFdb_based_domain_extraction/results/pombe_ted_domains.tsv
    output_filtered_file = os.path.join(script_dir, "..", "results", "pombe_ted_domains.tsv")
    
    # Normalize the path to resolve ".."
    output_filtered_file = os.path.normpath(output_filtered_file)

    extract_and_filter_ted_domains(pombase_file, ted_domains_file, output_filtered_file)