#!/usr/bin/env python3
"""
extract_binding_site.py  – receptor–peptide interface mapping  (May 2025)

For each PDB/mmCIF file:
• Treat the largest protein chain as the *receptor*.
• Treat every other protein chain as a *peptide*.
• Write every receptor residue whose atom is ≤ cut-off Å from any
  peptide atom.

Usage (example)
--------------------------------------
    python extract_binding_site.py *.pdb \
           --cutoff 4.0 --threads 6 --out peptide_interfaces.csv
"""
from __future__ import annotations
import argparse, csv, gzip, os, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path
from typing import List, Tuple

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, Selection   # pip install biopython
from tqdm import tqdm                                                   # pip install tqdm

WATERS = {"HOH", "WAT", "DOD", "D2O"}

# ──────────────── file helpers (PDB / CIF, plain or .gz) ─────────────── #
def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")

def load_structure(path: Path):
    txt  = open_text(path).read()
    hndl = StringIO(txt)
    if path.suffix.lower().endswith(".cif") or ".cif." in path.name:
        return MMCIFParser(QUIET=True).get_structure(path.stem, hndl)
    return PDBParser(QUIET=True).get_structure(path.stem, hndl)

# ───────────────────────── interface calculation ─────────────────────── #
def largest_chain(struct):
    """Return chain with most residues (ties → most atoms)."""
    chains = list(struct.get_chains())
    chains.sort(key=lambda c: (len(c), len(list(c.get_atoms()))), reverse=True)
    return chains[0]

def protein_atoms(chain):
    """All standard-residue atoms from a chain (no HETATM)."""
    return [a for a in chain.get_atoms() if a.get_parent().id[0] == " "]

def find_contacts(path: Path, cutoff: float) -> List[Tuple]:
    st     = load_structure(path)
    rec    = largest_chain(st)
    rec_atoms = protein_atoms(rec)

    pep_atoms = [a for c in st.get_chains() if c != rec
                   for a in protein_atoms(c)]
    if not pep_atoms:
        return []          # nothing but the receptor present

    neigh = NeighborSearch(rec_atoms)          # KD-tree on receptor
    hits  = set()
    for a in pep_atoms:
        for nb in neigh.search(a.coord, cutoff):
            rres  = nb.get_parent()
            pchn  = a.get_parent().get_parent().id
            rcid  = rec.id
            resid, icode = rres.id[1], rres.id[2].strip() or "-"
            hits.add((path.name, rcid, pchn, resid, icode, rres.resname))
    return sorted(hits, key=lambda x: (x[1], x[2], x[3]))

# ───────────────────────────────── main ───────────────────────────────── #
def parse_cli(argv=None):
    ap = argparse.ArgumentParser(description="Map receptor–peptide interfaces")
    ap.add_argument("pdb", nargs="*", default=["**/*.pdb", "**/*.cif", "**/*.pdb.gz", "**/*.cif.gz"],
                    help="files / globs / directories (defaults recurse for *.pdb/*cif)")
    ap.add_argument("--cutoff", type=float, default=4.0,
                    help="distance cut-off in Å  [4.0]")
    ap.add_argument("--threads", type=int, default=os.cpu_count(),
                    help="parallel workers [CPU cores]")
    ap.add_argument("--out", default="peptide_interfaces.csv",
                    help="output CSV filename [peptide_interfaces.csv]")
    return ap.parse_args(argv)

def gather_paths(items):
    pats = []
    for it in items:
        p = Path(it)
        pats.append(str(p / "**/*") if p.is_dir() else it)
    seen = set()
    for pat in pats:
        for f in Path().glob(pat):
            if f.is_file():
                seen.add(f)
    return sorted(seen)

def main(argv=None):
    args      = parse_cli(argv)
    pdb_files = gather_paths(args.pdb)
    if not pdb_files:
        sys.exit("✗ No structure files found – check paths/globs.")

    with ThreadPoolExecutor(max_workers=args.threads) as pool, \
         open(args.out, "w", newline="") as fh:
        writer  = csv.writer(fh)
        writer.writerow(["pdb", "receptor_chain",
                         "peptide_chain", "resi", "icode", "resn"])
        futs = {pool.submit(find_contacts, p, args.cutoff): p for p in pdb_files}
        for fut in tqdm(as_completed(futs), total=len(futs),
                        unit="file", desc="Processing"):
            writer.writerows(fut.result())

    print(f"✔ {len(pdb_files)} file(s) → {args.out}")

if __name__ == "__main__":
    main()
