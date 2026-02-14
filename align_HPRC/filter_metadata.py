import os
import sys
import csv
import argparse

def main(argv=None):
    parser = argparse.ArgumentParser(description='Filter metadata CSV based on coverage')
    parser.add_argument('--infile', '-i', required=True, help='Input metadata CSV file')
    parser.add_argument('--outfile', '-o', required=True, help='Output filtered metadata CSV file')
    parser.add_argument('--countfile', '-c', required=True, help='Output file to write count of filtered samples')
    parser.add_argument('--num_samples', type=int, default=None, help='Optional maximum number of samples to keep (default: no limit)')
    parser.add_argument('--chemistry', type=str, default=None, help='Optional sequencing chemistry to filter by (e.g. "R9" or "R10", default: no chemistry filter)')
    parser.add_argument('--min_cov', type=float, default=10.0, help='Minimum coverage threshold (default: 10.0)')
    parser.add_argument('--max_cov', type=float, default=1000.0, help='Maximum coverage threshold (default: 1000.0)')
    parser.add_argument('--keep_cols', nargs='*', default=['sample_ID', 'coverage', 'path', 'sequencing_chemistry'],
                        help='List of columns to keep. First three must be sample_ID, coverage, file path, or equivalent.')

    args = parser.parse_args(argv)
    # sequencing_chemistry column check if values starts with R9 or R10 and ignore additional values after that

    filter_metadata(infile=args.infile, outfile=args.outfile, countfile=args.countfile, min_cov=args.min_cov, max_cov=args.max_cov,
                    chemistry=args.chemistry, keep_cols=args.keep_cols, num_samples=args.num_samples)


def filter_metadata(infile='metadata.csv', outfile='filtered_metadata.csv', countfile='filtered_count.txt', min_cov=10.0, max_cov=1000.0,
                    chemistry=None, keep_cols=None, num_samples=None):
    """Filter metadata CSV based on coverage.

    Keeps the highest-coverage row per sample where coverage >= min_cov.
    Expects the CSV to contain at least 'sample_ID' and 'path' columns. If a
    'coverage' column exists it will be used for filtering and choosing the
    best sample; otherwise rows without coverage are skipped.
    """
    skip_samples = ["HG002", "HG005"]

    newheader = list(keep_cols) + ['readgroup']
    if len(newheader) < 3:
        raise SystemExit('keep_cols must include at least sample_ID, coverage and path column names')
    id_col = newheader[0]
    cov_col = newheader[1]
    path_col = newheader[2]

    if not os.path.exists(infile):
        raise SystemExit(f"input file not found: {infile}")

    with open(infile, newline='') as inf:
        reader = list(csv.reader(inf))
    if not reader:
        raise SystemExit(1)
    header = reader[0]
    # Build index map for requested keep_cols. Missing cols will be written
    # as empty values in the output. Require that id_col and path_col exist.
    idx_map = {col: (header.index(col) if col in header else None) for col in newheader}
    if idx_map[id_col] is None or idx_map[path_col] is None:
        raise SystemExit('metadata CSV must contain sample_ID and path')
    idx_id = idx_map[id_col]
    idx_path = idx_map[path_col]
    idx_cov = idx_map[cov_col] if cov_col in header else -1
    idx_chemistry = header.index('sequencing_chemistry') if 'sequencing_chemistry' in header else -1

    rows = []
    for row in reader[1:]:
        # allow rows shorter than header (treated as missing values)
        if idx_id is None or idx_path is None:
            continue
        sid = (row[idx_id].strip() if idx_id < len(row) else '').strip()
        if not sid:
            continue
        if sid in skip_samples:
            continue
        path = (row[idx_path].strip() if idx_path < len(row) else '').strip()
        cov = None
        if idx_cov != -1 and idx_cov is not None and idx_cov < len(row):
            try:
                cov = float(row[idx_cov])
            except Exception:
                cov = None
        if not path:
            continue
        if cov is None or cov < float(min_cov) or cov > float(max_cov):
            continue
        # Filter by chemistry if specified
        if chemistry is not None and idx_chemistry != -1:
            chem = (row[idx_chemistry].strip() if idx_chemistry < len(row) else '').strip()
            if not chem.startswith(chemistry):
                continue
        # store the full original row so we can output requested columns
        rows.append((sid, cov, path, row))

    # pick highest coverage per sample and remember the original row
    best = {}
    for sid, cov, path, row in rows:
        if sid not in best or cov > best[sid][0]:
            best[sid] = (cov, row)

    # write filtered metadata preserving requested keep_cols order
    with open(outfile, 'w', newline='') as outf:
        writer = csv.writer(outf)
        writer.writerow(newheader)
        best_items = list(best.items())
        if num_samples is not None:
            best_items = best_items[:num_samples]
        for sid, (cov, row) in best_items:
            cov_out = int(cov) if float(cov).is_integer() else cov
            out_vals = []
            for col in newheader:
                idx = idx_map.get(col)
                if idx is not None and idx < len(row):
                    out_vals.append(row[idx])
                else:
                    # if this is the coverage column but it was missing in input,
                    # use the cov we parsed/selected
                    if col == cov_col:
                        out_vals.append(str(cov_out))
                    else:
                        out_vals.append('')
            # derive read group string from metadata
            try:
                rg = derive_rg_from_metadata({header[i]: row[i] for i in range(len(row)) if i < len(header)})
            except Exception as e:
                row_data = {header[i]: row[i] for i in range(len(row)) if i < len(header)}
                raise SystemExit(f"Error deriving read group for sample {sid}: {e} \nrun_index: {row_data.get('run_index')} \nfilename: {row_data.get('filename')}")
            out_vals[-1] = rg
            writer.writerow(out_vals)

    with open(countfile, 'w') as cf:
        cf.write(str(len(best_items)) + "\n")

import re

def extract_run_index_from_filename(filename):
    """
    Extract run index from filenames like:
      - ..._1_...
      - ...part02_...
    Returns a normalized integer string (no leading zeros).
    """

    # Case 1: partXX
    m = re.search(r'(?:^|[_.])part(\d+)(?:[_.]|$)', filename, re.IGNORECASE)
    if m:
        return str(int(m.group(1)))

    # Case 2: underscore-delimited integer token (fallback)
    m = re.search(r'_(\d+)_', filename)
    if m:
        return str(int(m.group(1)))
    
    return "NA"  # default if no index found


def derive_rg_from_metadata(row):
    """
    Derive RG PU and ID from structured metadata.

    Required fields in row:
      - sample_ID
      - library_ID
      - platform
      - sequencing_chemistry
      - filename (only for run index fallback)
    """

    sample = row["sample_ID"]
    platform = row["platform"].replace("OXFORD_NANOPORE", "ONT")
    library_id = row["library_ID"]

    # Try to get run index from metadata; fallback to filename
    run_idx = row.get("run_index")

    if run_idx is None:
        run_idx = extract_run_index_from_filename(row["filename"])

    # Platform unit
    pu = f"{library_id}_{run_idx}"

    # Read group ID
    rg_id = f"{sample}_{platform}_{pu}"

    read_group = f"@RG\\tID:{rg_id}\\tLB:{library_id}\\tSM:{sample}\\tPL:{platform}\\tPU:{pu}"

    return read_group


if __name__ == '__main__':
    main()