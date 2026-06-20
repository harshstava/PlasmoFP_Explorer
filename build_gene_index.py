#!/usr/bin/env python3
"""
Build the Explorer's primary search data from the per-species
*_gene_dict_with_PFP.pkl files produced by run_puf_inference.py.

Reads:
  - with_PFP_predictions_complete_2/*_gene_dict_with_PFP.pkl
  - go_terms.json  (GO ID -> name, for resolving original-annotation names)

Writes (next to this script):
  - optimized_gene_index.json : per-gene species, original GO annotations,
                                and PFP predictions by ontology aspect + eFDR
  - search_index.json         : uppercased gene-ID prefixes (len 3-9) -> gene IDs

After running this, run preprocess_enhanced_data.py to rebuild the secondary
search indices (product/GO-term lookups) from optimized_gene_index.json.
"""

import glob
import json
import os
import pickle
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKL_DIR = HERE / "with_PFP_predictions_complete_2"
GO_TERMS_PATH = HERE / "go_terms.json"
GENE_INDEX_OUT = HERE / "optimized_gene_index.json"
SEARCH_INDEX_OUT = HERE / "search_index.json"

# (experimental slot, IEA slot, PFP slot) per ontology aspect
ASPECT_SLOTS = {
    "MF": ("GO Function",  "GO IEA Function",  "PFP MF"),
    "BP": ("GO Process",   "GO IEA Process",   "PFP BP"),
    "CC": ("GO Component", "GO IEA Component", "PFP CC"),
}

PREFIX_MIN = 3
PREFIX_MAX = 9


def species_from_filename(path):
    name = os.path.basename(path)
    return name.replace("PlasmoDB-68_", "").replace("_gene_dict_with_PFP.pkl", "")


def build_gene_index(pkl_dir, go_terms):
    gene_index = {}
    for fp in sorted(glob.glob(str(Path(pkl_dir) / "*_gene_dict_with_PFP.pkl"))):
        species = species_from_filename(fp)
        with open(fp, "rb") as f:
            records = pickle.load(f)

        for gene_id, rec in records.items():
            original = {}
            pfp = {}
            for aspect, (exp_slot, iea_slot, pfp_slot) in ASPECT_SLOTS.items():
                # Original annotations: experimental first, then IEA (no dedup).
                anns = list(rec.get(exp_slot, [])) + list(rec.get(iea_slot, []))
                if anns:
                    original[aspect] = [
                        {"id": go_id, "name": go_terms.get(go_id, go_id)}
                        for go_id in anns
                    ]

                # PFP predictions: per eFDR threshold, drop empties, sort by
                # score descending. Stored tuples are (go_id, score, name).
                fdr_map = {}
                for fdr, terms in rec.get(pfp_slot, {}).items():
                    if not terms:
                        continue
                    # Round scores to 6 dp: keeps all float32-meaningful
                    # precision (the app only displays round(score, 4)) while
                    # keeping optimized_gene_index.json under GitHub's 100 MB
                    # per-file limit.
                    fdr_map[str(fdr)] = sorted(
                        ({"id": go_id, "name": name, "score": round(float(score), 6)}
                         for go_id, score, name in terms),
                        key=lambda t: t["score"],
                        reverse=True,
                    )
                if fdr_map:
                    pfp[aspect] = fdr_map

            # Keep only genes with at least one annotation or prediction.
            if original or pfp:
                gene_index[gene_id] = {
                    "species": species,
                    "original_annotations": original,
                    "pfp_predictions": pfp,
                }

    return gene_index


def build_search_index(gene_index):
    search_index = defaultdict(list)
    for gene_id in gene_index:
        upper = gene_id.upper()
        for i in range(PREFIX_MIN, min(len(upper), PREFIX_MAX) + 1):
            search_index[upper[:i]].append(gene_id)
    return search_index


def main():
    with open(GO_TERMS_PATH) as f:
        go_terms = json.load(f)

    print(f"Reading pkls from {PKL_DIR} ...")
    gene_index = build_gene_index(PKL_DIR, go_terms)
    print(f"  {len(gene_index)} genes retained")

    search_index = build_search_index(gene_index)
    print(f"  {len(search_index)} search-index prefixes")

    # Compact separators (no whitespace) to minimise file size.
    with open(GENE_INDEX_OUT, "w") as f:
        json.dump(gene_index, f, separators=(",", ":"))
    with open(SEARCH_INDEX_OUT, "w") as f:
        json.dump(search_index, f, separators=(",", ":"))
    print(f"Wrote {GENE_INDEX_OUT.name} and {SEARCH_INDEX_OUT.name}")


if __name__ == "__main__":
    main()
