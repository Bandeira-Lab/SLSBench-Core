"""
balance_tp_fp.py
================

Command-line tool to balance true-positive (TP) and false-positive (FP) query spectra for a spectral-library-search benchmark by partitioning queries into TP/FP groups via connected components of a query<->library graph built on charge + I/L-and-ambiguity-collapsed peptide sequences.

The script reads a query-info table and a library-info table, builds a graph between query scans and library scans that share the same charge and a (possibly ambiguous) demodified sequence, finds isolated subgraphs, and greedily assigns each subgraph's queries to TP or FP so that the two sets stay roughly balanced. Library spectra that would otherwise turn an FP-assigned subgraph into a recoverable TP are recorded as "to drop".

------------------------------------------------------------------------------
Required input files and columns
------------------------------------------------------------------------------

1) --query_info_path : TSV (tab-separated, with header)
   One row per query spectrum. Required columns:
     * query_scan     (int)  : unique scan id of the query spectrum.
     * query_seq      (str)  : top 1 MSGF+ identified peptide sequence of the query spectrum
     * query_charge   (int)  : precursor charge state of the query.
     * query_seq_topk (str)  : top k MSGF+ identified peptide sequences of the query spectrum, joined by '|'.
     * in_recall      (0/1)  : 1 if library contains a spectrum with the same precursor for query with top 1 msgf+ sequence (same charge, precursor m/z difference within 0.05 Da)
     * in_recall_topk (0/1)  : 1 if library contains a spectrum with the same precursor for query with any of the top k msgf+ sequences

2) --library_info_path : TSV (tab-separated, with header)
   One row per library spectrum. Required columns:
     * scan   (int) : unique scan id of the library spectrum.
     * seq    (str) : peptide sequence of the library entry (same format
                      conventions as query_seq).
     * charge (int) : precursor charge state of the library entry.

------------------------------------------------------------------------------
Output directory contents (--output_dir)
------------------------------------------------------------------------------

The directory is created if it does not exist. The following files are written:

  * isolated_subgraphs.pkl
        Pickled list[set[str]]. Each set is the node set of one connected
        component of the query<->library graph. Query nodes are encoded as
        "q_<query_scan>" and library nodes as "l_<library_scan>". Saved for
        downstream inspection / reproducibility of the partition.

  * query_scan_to_use.tsv
        Plain-text file, one original query_scan id per line, sorted
        ascending. Concatenation of the final TP and FP sets; defines the
        old->new scan-id remapping (line index i (1-based) becomes the new
        scan id of the original scan on that line).

  * tp_new_scan.tsv
        Plain-text file, one *remapped* (new) scan id per line, listing the
        queries assigned to the TP set.

  * fp_new_scan.tsv
        Plain-text file, one *remapped* (new) scan id per line, listing the
        queries assigned to the FP set.

  * library_scan_to_use.tsv
        Plain-text file, one library `scan` id per line, sorted ascending,
        listing the library entries that should be retained (i.e. all
        library scans minus those marked for dropping during FP assignment).

------------------------------------------------------------------------------
Example
------------------------------------------------------------------------------

  python balance_tp_fp.py \\
      --query_info_path   data/query_info_with_msgf_ids.tsv \\
      --library_info_path data/library_info.tsv \\
      --output_dir        output/
"""

import argparse
import os
import pickle
import re
from collections import defaultdict

import pandas as pd


# =============== sequence / graph utilities ===============

def get_all_ambiguous_sequence_matches_demodified(seq):
    """Given a peptide sequence, return the set of all demodified sequences
    obtained by stripping non-uppercase characters, collapsing I->L, and
    enumerating the standard ambiguous substitutions (Q<->K, Q<->E, D<->N).
    """
    seq = seq.replace('I', 'L')
    seq = re.sub(r'[^A-Z]', '', seq.replace('I', 'L'))

    seq_mod_queue = [seq]
    for idx, s in enumerate(seq):
        if s == 'Q':
            seq_mod_queue.append(seq[:idx] + 'K' + seq[idx + 1:])
            seq_mod_queue.append(seq[:idx] + 'E' + seq[idx + 1:])
        elif s == 'K':
            seq_mod_queue.append(seq[:idx] + 'Q' + seq[idx + 1:])
        elif s == 'E':
            seq_mod_queue.append(seq[:idx] + 'Q' + seq[idx + 1:])
        elif s == 'D':
            seq_mod_queue.append(seq[:idx] + 'N' + seq[idx + 1:])
        elif s == 'N':
            seq_mod_queue.append(seq[:idx] + 'D' + seq[idx + 1:])

    return set(seq_mod_queue)


def generate_charge_to_demodified_seq_to_scan_set_for_query_topk(query_df, query_scan_to_topk_id_list=None):
    """Build, for the top-k candidates of each query, mappings:
        charge -> demodified_seq -> set(query_scan)
        query_scan -> set(demodified_seq)
    """
    if query_scan_to_topk_id_list is None:
        query_scan_to_topk_id_list = defaultdict(list)

    if len(query_scan_to_topk_id_list) == 0:
        for s, top_k in query_df[["query_scan", "query_seq_topk"]].values:
            query_scan_to_topk_id_list[s] = top_k.split('|')

    query_scan_to_demodified_seq_set_topk = defaultdict(set)
    query_charge_to_demodified_seq_to_scan_set_topk = defaultdict(lambda: defaultdict(set))
    for s, _, c in query_df[["query_scan", "query_seq", "query_charge"]].values:
        if s not in query_scan_to_topk_id_list:
            print(f"Query scan {s} not in topk list")
            continue
        seq_list = query_scan_to_topk_id_list[s]
        for q in seq_list:
            query_scan_to_demodified_seq_set_topk[s].update(get_all_ambiguous_sequence_matches_demodified(q))
        for dq in query_scan_to_demodified_seq_set_topk[s]:
            query_charge_to_demodified_seq_to_scan_set_topk[c][dq].add(s)
    return query_charge_to_demodified_seq_to_scan_set_topk, query_scan_to_demodified_seq_set_topk


def generate_charge_to_demodified_seq_to_scan_set(library_df, mode='library'):
    """Build mappings charge -> demodified_seq -> set(scan) and
    scan -> set(demodified_seq) for either the library table (columns
    scan/seq/charge) or the query table (columns query_scan/query_seq/
    query_charge).
    """
    if mode == 'library':
        tmp = library_df[["scan", "seq", "charge"]].values
    elif mode == 'query':
        tmp = library_df[["query_scan", "query_seq", "query_charge"]].values
    else:
        raise ValueError('mode should be either "library" or "query"')

    library_charge_to_demodified_seq_to_scan_set = defaultdict(lambda: defaultdict(set))
    library_scan_to_demodified_seq = defaultdict(set)
    for s, q, c in tmp:
        demodified_seq_list = get_all_ambiguous_sequence_matches_demodified(q)
        for dq in demodified_seq_list:
            library_charge_to_demodified_seq_to_scan_set[c][dq].add(s)
        library_scan_to_demodified_seq[s].update(demodified_seq_list)
    return library_charge_to_demodified_seq_to_scan_set, library_scan_to_demodified_seq


def find_isolated_subgraphs(graph):
    """Return all connected components of an undirected adjacency-list graph
    as a list of node sets, using iterative DFS.
    """
    visited = set()
    all_nodes = set(graph.keys())
    subgraphs = []

    while len(visited) < len(all_nodes):
        for node in all_nodes:
            if node not in visited:
                start_node = node
                break

        stack = [start_node]
        subgraph = set()

        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                subgraph.add(node)
                for neighbor in graph[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
            if len(visited) % 20000 == 0:
                print(f"Visited {len(visited)} nodes - current number of subgraphs: {len(subgraphs)}")

        subgraphs.append(subgraph)

    return subgraphs


# =============== main pipeline ===============

def balance_tp_fp(query_info_path, library_info_path, output_dir):
    """Main pipeline to balance TP/FP queries by partitioning the query<->library graph into connected components and greedily assigning each component to TP or FP based on current set sizes. Writes outputs to output_dir.
    """

    os.makedirs(output_dir, exist_ok=True)

    # ----- load -----
    query_df = pd.read_csv(query_info_path, sep='\t')
    print("{} queries, {} in recall, {} in recall topk".format(
        query_df.shape[0], query_df.in_recall.sum(), query_df.in_recall_topk.sum()))
    query_df.sort_values('query_scan', inplace=True, ignore_index=True)

    library_df = pd.read_csv(library_info_path, sep='\t')
    print("{} library spectra".format(library_df.shape[0]))
    library_df.sort_values('scan', inplace=True, ignore_index=True)

    # ----- TP/FP initialization -----
    tp_scan_set_final = set()
    fp_scan_set_final = set(query_df[query_df.in_recall_topk == 0].query_scan)
    library_to_drop = set()
    print(f"TP: {len(tp_scan_set_final)}, FP: {len(fp_scan_set_final)}")

    # ----- top-k query side -----
    query_charge_to_demodified_seq_to_scan_set_topk, query_scan_to_demodified_seq_set_topk = \
        generate_charge_to_demodified_seq_to_scan_set_for_query_topk(query_df[query_df.in_recall_topk == 1])
    query_df['num_demodified_seq_topk'] = [
        len(query_scan_to_demodified_seq_set_topk[s]) for s in query_df.query_scan
    ]

    # ----- library side -----
    library_charge_to_demodified_seq_to_scan_set, _ = \
        generate_charge_to_demodified_seq_to_scan_set(library_df, mode='library')

    # ----- build bipartite query<->library graph -----
    id_query_library_graph = defaultdict(set)  # "q_<scan>" <-> "l_<scan>"
    for c, v in query_charge_to_demodified_seq_to_scan_set_topk.items():
        for id_, s_set in v.items():
            for s in s_set:
                for ls in library_charge_to_demodified_seq_to_scan_set[c][id_]:
                    id_query_library_graph[f"q_{s}"].add(f"l_{ls}")
                    id_query_library_graph[f"l_{ls}"].add(f"q_{s}")
    print(f"Number of nodes: {len(id_query_library_graph)}")
    print(f"number of query nodes: {len([k for k in id_query_library_graph.keys() if k.startswith('q_')])}")
    print(f"number of library nodes: {len([k for k in id_query_library_graph.keys() if k.startswith('l_')])}")
    print(f"Number of edges: {int(sum([len(v) for v in id_query_library_graph.values()]) / 2)}")
    assert len([k for k in id_query_library_graph.keys() if k.startswith('q_')]) + len(fp_scan_set_final) == query_df.shape[0]

    # ----- isolated subgraphs -----
    isolated_subgraphs = find_isolated_subgraphs(id_query_library_graph)
    with open(os.path.join(output_dir, 'isolated_subgraphs.pkl'), 'wb') as f:
        pickle.dump(isolated_subgraphs, f)
    print(f"Number of isolated subgraphs: {len(isolated_subgraphs)}")

    # ----- greedy TP/FP assignment, balanced by current set sizes -----
    sorted_isolated_subgraphs = sorted(
        isolated_subgraphs,
        key=lambda x: len([n for n in x if n.startswith('q_')]),
        reverse=True,
    )
    for subgraph in sorted_isolated_subgraphs:
        query_nodes = set(int(n.split("_")[1]) for n in subgraph if n.startswith('q_'))
        library_nodes = set(int(n.split("_")[1]) for n in subgraph if n.startswith('l_'))

        if len(query_nodes) == 0:
            print("No query nodes in the subgraph")
        elif len(library_nodes) == 0:
            fp_scan_set_final.update(query_nodes)
        else:
            if len(tp_scan_set_final) <= len(fp_scan_set_final):
                tp_scan_set_final.update(query_nodes)
            else:
                fp_scan_set_final.update(query_nodes)
                library_to_drop.update(library_nodes)
    print(f"TP: {len(tp_scan_set_final)}, FP: {len(fp_scan_set_final)}, library_to_drop: {len(library_to_drop)}")
    assert len(tp_scan_set_final.intersection(fp_scan_set_final)) == 0
    assert len(tp_scan_set_final) + len(fp_scan_set_final) == query_df.shape[0]

    # ----- write outputs -----
    query_scan_to_use = list(tp_scan_set_final) + list(fp_scan_set_final)
    query_scan_to_use.sort()
    out_path = os.path.join(output_dir, 'query_scan_to_use.tsv')
    with open(out_path, 'w') as f:
        f.write('\n'.join([str(s) for s in query_scan_to_use]))
    print(f"number of query scan to use: {len(query_scan_to_use)}")

    old_scan_to_current = {j: i + 1 for i, j in enumerate(query_scan_to_use)}

    out_path = os.path.join(output_dir, 'tp_new_scan.tsv')
    with open(out_path, 'w') as f:
        f.write('\n'.join([str(old_scan_to_current[s]) for s in tp_scan_set_final]))
    print(f"number of TP scan to use: {len(tp_scan_set_final)}")

    out_path = os.path.join(output_dir, 'fp_new_scan.tsv')
    with open(out_path, 'w') as f:
        f.write('\n'.join([str(old_scan_to_current[s]) for s in fp_scan_set_final]))
    print(f"number of FP scan to use: {len(fp_scan_set_final)}")

    library_scan_to_use = list(set(library_df.scan).difference(library_to_drop))
    library_scan_to_use.sort()
    out_path = os.path.join(output_dir, 'library_scan_to_use.tsv')
    with open(out_path, 'w') as f:
        f.write('\n'.join([str(s) for s in library_scan_to_use]))
    print(f"number of library scan to use: {len(library_scan_to_use)}")


# =============== CLI ===============

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Balance TP/FP queries for spectral-library-search benchmarks by "
            "partitioning queries via connected components of a charge + "
            "demodified-sequence query<->library graph."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Required input columns:\n"
            "  query_info_path (TSV): query_scan, query_seq, query_charge,\n"
            "                         query_seq_topk ('|'-separated top-k seqs),\n"
            "                         in_recall (0/1), in_recall_topk (0/1)\n"
            "  library_info_path (TSV): scan, seq, charge\n\n"
            "Output files written into output_dir:\n"
            "  isolated_subgraphs.pkl   - pickled list[set[str]] of graph components\n"
            "                             (nodes 'q_<scan>' / 'l_<scan>')\n"
            "  query_scan_to_use.tsv    - sorted original query_scan ids kept\n"
            "                             (defines the old->new scan remap)\n"
            "  tp_new_scan.tsv          - remapped (new) scan ids of TP queries\n"
            "  fp_new_scan.tsv          - remapped (new) scan ids of FP queries\n"
            "  library_scan_to_use.tsv  - sorted library `scan` ids retained\n"
        ),
    )
    parser.add_argument(
        "--query_info_path", required=True,
        help=(
            "Path to the query-info TSV. Must contain columns: query_scan (int), "
            "query_seq (str peptide), query_charge (int), query_seq_topk "
            "(str, top-k candidate peptides joined by '|'), in_recall (0/1), "
            "in_recall_topk (0/1, queries with 0 are forced to FP)."
        ),
    )
    parser.add_argument(
        "--library_info_path", required=True,
        help=(
            "Path to the library-info TSV. Must contain columns: scan (int), "
            "seq (str peptide), charge (int)."
        ),
    )
    parser.add_argument(
        "--output_dir", required=True,
        help=(
            "Directory to write outputs into (created if missing). Files: "
            "isolated_subgraphs.pkl (component dump), query_scan_to_use.tsv "
            "(sorted kept original query scans defining old->new remap), "
            "tp_new_scan.tsv / fp_new_scan.tsv (remapped TP / FP scan ids), "
            "library_scan_to_use.tsv (sorted retained library scan ids)."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    balance_tp_fp(
        query_info_path=args.query_info_path,
        library_info_path=args.library_info_path,
        output_dir=args.output_dir,
    )
