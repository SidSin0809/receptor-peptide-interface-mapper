# receptor-peptide-interface-mapper
High-throughput identification of receptor–peptide contact residues from PDB/mmCIF structures (Python 3, Biopython, multithreaded)


**extract_binding_site.py** is a command-line utility that locates all receptor
residues in close contact with peptide chains for any number of PDB/mmCIF
structures (optionally .gz-compressed).  
It assumes the *largest* protein chain is the receptor and treats all other
protein chains as peptides.  Results are written to a single tidy CSV file.

| Key points | |
|-----------------------|---------------------------------------------------------------------------------------------------------------------------------|
| **Fast & parallel**   | Uses Python `ThreadPoolExecutor`; scales to hundreds of structures in minutes                                                     |
| **Flexible input**    | Accepts individual files, unix‐style globs, or directories (recurses automatically)                                             |
| **Robust parsing**    | Handles `.pdb`, `.cif`, and their `.gz` versions via Biopython                                                                  |
| **Clean output**      | Produces a 6-column CSV ready for downstream analysis or visualisation                                                          |
| **Pure Python 3**     | Only external dependencies are `biopython`, `numpy` (pulled by Biopython), and `tqdm` for a neat progress bar                   |

---

## Installation

```bash
# 1. Clone the repo
git clone https://github.com/SidSin0809/receptor-peptide-interface-mapper.git

# 2. Install requirements
pip install -r requirements.txt

Quick start
# Example: process all *.pdb files in current & nested folders,
#          using 8 threads and a 4 Å cut-off
python extract_binding_site.py "**/*.pdb" \
       --cutoff 4.0 \
       --threads 8 \
       --out peptide_interfaces.csv

CLI options
| Flag        | Default | Description                                  |
| ----------- | ------- | -------------------------------------------- |
| `pdb ...`   | (none)  | File paths, globs, or directories to analyse |
| `--cutoff`  |   4.0   | Atom–atom distance threshold (Å)             |
| `--threads` |   CPU   | Number of worker threads (≤ CPU cores)       |
| `--out`     |   csv   | Output filename                              |


Citation
If this tool contributes to academic work, please cite:

Singh S. (2025) Receptor–Peptide Interface Mapper. GitHub repository.
https://github.com/SidSin0809/receptor-peptide-interface-mapper

