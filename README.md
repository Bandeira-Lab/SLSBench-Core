# Construct a balanced TP/FP benchamrk for spectral library search

Command-line tool to balance true-positive (TP) and false-positive (FP) query
spectra for a spectral-library-search benchmark. It partitions queries into
TP/FP groups via connected components of a query↔library graph built on
charge + I/L-and-ambiguity-collapsed peptide sequences, and outputs the
remapped scan-id lists into TSV files to build the balanced benchmark.

---

## 1. Download the input data

Two TSV files are required. They are provided together in a single Google
Drive folder:

**Google Drive folder:** <https://drive.google.com/drive/folders/1XMOwPECxu8bfphlQBne6gwD2oibsw56R?usp=sharing>

The folder contains:

- `query_info_with_msgf_ids.tsv`
- `library_info.tsv`

Download both files and place them under `data/` in this directory so the
layout looks like:

```text
core_code/
├── balance_tp_fp.py
├── README.md
└── data/
    ├── query_info_with_msgf_ids.tsv   <-- download into here
    └── library_info.tsv               <-- download into here
```

Create the `data/` directory if it does not already exist:

```bash
mkdir -p data
# then move the downloaded files into ./data/
```

---

## 2. Required input columns

### `data/query_info_with_msgf_ids.tsv` (TSV with header)

One row per query spectrum. Required columns:

| Column           | Type | Description                                                                                                                              |
| ---------------- | ---- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `query_scan`     | int  | Unique scan id of the query spectrum.                                                                                                    |
| `query_charge`   | int  | Precursor charge state of the query.                                                                                                     |
| `query_seq`      | str  | Top-1 MSGF+ identified peptide sequence of the query spectrum.                                                                           |
| `query_seq_topk` | str  | Top-k MSGF+ identified peptide sequences of the query spectrum, joined by `\|`.                                                          |
| `in_recall`      | 0/1  | 1 if the library contains a spectrum with the same precursor (same charge, precursor m/z difference within 0.05 Da) for the top-1 seq.  |
| `in_recall_topk` | 0/1  | 1 if the library contains a spectrum with the same precursor for any of the top-k seqs. Rows with 0 are forced into FP.                 |

### `data/library_info.tsv` (TSV with header)

One row per library spectrum. Required columns:

| Column   | Type | Description                                  |
| -------- | ---- | -------------------------------------------- |
| `scan`   | int  | Unique scan id of the library spectrum.      |
| `seq`    | str  | Peptide sequence of the library entry.       |
| `charge` | int  | Precursor charge state of the library entry. |

Extra columns are ignored — but only these are needed.

---

## 3. Run the command

The project uses the `data` conda environment (pandas required):

```bash
python balance_tp_fp.py \
    --query_info_path   data/query_info_with_msgf_ids.tsv \
    --library_info_path data/library_info.tsv \
    --output_dir        output/
```

Run `python balance_tp_fp.py --help` for the full argument reference.

---

## 4. Output files (written to `--output_dir`)

The output directory is created if it does not exist.

| File                       | Purpose                                                                                                                                                 |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `isolated_subgraphs.pkl`   | Pickled `list[set[str]]`. Each set is the node set of one connected component of the query↔library graph. Nodes encoded as `q_<query_scan>` / `l_<library_scan>`. Saved for inspection / reproducibility of the partition. |
| `query_scan_to_use.tsv`    | Plain text, one **original** `query_scan` id per line, sorted ascending. Concatenation of the final TP and FP sets; defines the old→new scan-id remap (the line index, 1-based, is the new scan id of the original scan on that line). |
| `tp_new_scan.tsv`          | Plain text, one **remapped (new)** scan id per line, listing the queries assigned to the TP set.                                                        |
| `fp_new_scan.tsv`          | Plain text, one **remapped (new)** scan id per line, listing the queries assigned to the FP set.                                                        |
| `library_scan_to_use.tsv`  | Plain text, one library `scan` id per line, sorted ascending, listing library entries to retain (all library scans minus those dropped during FP assignment). |
