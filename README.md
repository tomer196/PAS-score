# Motif-based Synthetic Accessibility Score (PASscore)
PASscore is a synthetic accessibility metric for PASs (polycyclic aromatic systems) based on motif frequency analysis. The method analyzes all PASs in the PubChem dataset and counts the frequency of each ring-based motif (the fundamental building blocks). Each motif receives a score based on its relative frequency. To score a new molecule, all its motifs are extracted, each motif's score is retrieved, and the scores are combined as described above to yield the molecule's PASscore.

## Detailed Formulation

### Motif Score Calculation
Each motif is assigned a synthetic accessibility (SA) score based on its frequency in the PubChem PAS dataset. The process is as follows:
1. **Count Motif Occurrences:** For each motif (aromatic ring-based substructure), count its total appearances across all PAS molecules (see `sum_motifs_from_backbones` in `extract_PASs_motifs.py`).
2. **Score Assignment:** The motif score is calculated as the logarithm of its count:  
	$$\text{SA score per motif} = \log_{10}(\text{motif count})$$
	(see `calc_pas_score_per_motif` in `extract_PASs_motifs.py`).
3. **Normalization:** Optionally, scores can be normalized. Because the finale score is normalized to be between 0 to 10 it doesn't has any effect.

### PAS Score Calculation
To score a new molecule:
1. **Extract Motifs:** All motifs are extracted from the molecule (see `get_all_motifs` in `extract_motifs.py`).
2. **Score Aggregation:** For each motif, retrieve its precomputed score. Motifs not found in the database are assigned a default penalty score.
3. **Combine Scores:**
   - First, group all motifs by their number of rings.
   - For each group (same ring count), compute the average motif score.
   - The final PAS score for the molecule is the average of these per-ring-count averages (see `PASScore` class in `pas_score.py`).
   - This two-step averaging ensures that motifs with different ring sizes contribute equally, regardless of their abundance.
4. **Normalization:** The score is normalized to a [0, 10] scale (higher = less accessible):
    $$\text{PASscore} = \frac{\text{max} - \text{score}}{\text{max} - \text{min}} \times 10$$
    where `min` and `max` are the observed motif score range (hardcoded).

---
# File-by-file Explanation
The typical workflow is:
1. **Extract motif scores** using `extract_PASs_motifs.py` (process PubChem, extract motifs, assign scores).
2. **Score molecules** using `pas_score.py` (extract motifs from new molecules, aggregate motif scores).
## extract_PASs_motifs.py
This script orchestrates the full pipeline for motif-based scoring:
- **filter_pas_from_pubchem:** Filters PAS molecules from a PubChem SMILES file, saving PAS-only entries.
- **extract_backbones:** Extracts backbone structures (fused aromatic ring systems) from PAS molecules.
- **count_backbones:** Counts unique backbone structures.
- **filter_catacondensed_backbones:** Filters for catacondensed backbones (no atom in more than 2 rings).
- **sum_motifs_from_backbones:** Extracts and counts all motifs from catacondensed backbones.
- **calc_pas_score_per_motif:** Assigns a score to each motif based on its frequency (log10 of count).  

Each function processes large datasets in parallel for efficiency. The main script runs all steps sequentially to produce the motif score database.

## pas_score.py
Defines the `PASScore` class for scoring molecules:
- **PASScore:** Loads the motif score database and provides a callable interface to score molecules (SMILES or InChI). Extracts motifs, retrieves their scores, averages per ring count, and normalizes the result.
- **n_rings:** Utility to count rings in a motif.   

This is the main interface for applying the motif-based PAS score to new molecules.

## extract_motifs.py
Contains helper functions for motif and backbone extraction:
- **extract_backbone_inchis:** Extracts all backbone InChIs (fused aromatic clusters) from a molecule.
- **get_all_motifs:** Extracts all possible aromatic motifs (connected ring subgraphs) from a molecule or backbone.
- **extract_induced_submol:** Builds submolecules for motifs, using placeholder atoms (e.g., F, Se, Cl) to cap external bonds or represent fused rings. Placeholders are replaced with hydrogens for visualization.

This file also handles molecule sanitization and ring system graph construction.

### About Atom Placeholders in Motif Extraction

During motif extraction (see `extract_motifs.py`), special placeholder atoms are used to represent certain chemical features when building submolecules for motifs:

- **Why are placeholders used?**
  - When extracting a motif (a fused ring system) from a larger molecule, the motif may have bonds to atoms outside the motif. To preserve the motif's chemical structure and valency, these external connections are replaced with placeholder atoms.
  - This allows the motif to be represented as a chemically valid molecule, suitable for conversion to InChI/SMILES and for further analysis.

- **Which placeholders are used and when?**
  - **F (Fluorine):** Used as a generic placeholder for most external single bonds (substituents) attached to the motif.
  - **Se (Selenium):** Used as a placeholder for external double bonds (e.g., exocyclic double bonds to the motif).
  - **Cl (Chlorine):** Used to represent fused ring connections when extracting motifs in 'ring mode', marking where rings are fused together.

- **How are placeholders handled?**
  - Placeholders are only used internally for motif extraction and representation. When visualizing motifs, placeholders may be replaced with hydrogens for clarity.
  - The use of uncommon atoms (like Se or Cl) as placeholders ensures they do not conflict with typical organic atoms in PAS motifs, making motif identification and deduplication robust.

## analyze_motif_count.py
Provides visualization and analysis tools:
- **plot_top_motifs:** Plots the most common motifs as molecule grids.
- **plot_motifs_frequancy:** Plots the frequency distribution of motif counts (log scale).
- **analyze:** Prints statistics on motif count distributions (e.g., how many motifs appear above certain thresholds).
These tools help interpret motif diversity and prevalence in the dataset.

## analyze_pas_scores.py
Performs analysis and visualization of PAS scores:
- **calculate_all_backbone_pas_scores:** Computes PAS scores for all backbones in the dataset.
- **plot_pas_score_distribution:** Plots the distribution of PAS scores.
- **sample_molecules_by_score:** Samples molecules across the PAS score range for expert evaluation, saving images and scores.
- **analyze_score_correlation:** Compares PAS scores to other metrics (e.g., SAscore), computes correlations, and plots scatter/rank plots.
This script is useful for benchmarking and expert review.

## other_metrics/add_other_scores.py
Calculates and compares additional synthetic accessibility metrics from the literature:
- **get_sascore, get_rascore, get_rscore, get_scscore, get_fscore:** Compute various published SA metrics for a list of molecules.
- **add_scores2csv:** Adds all scores to a CSV file for comparison.
- **plot_scores_correlation:** Plots and computes correlation between each metric and expert scores.
This enables comprehensive evaluation of the motif-based PAS score against established methods.

## Cloning and Setup Instructions

To clone this repository and set up the environment:

1. **Clone the repository including submodules:**
   ```bash
   git clone --recurse-submodules git@github.com:tomer196/PAS-score.git
   ```
   This will clone the main repository and all submodules (including `other_metrics/fsscore`).

2. **Create a Python virtual environment:**
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Cloning with Git LFS

Some files in this repository are tracked with Git Large File Storage (LFS).
To ensure you get all large files:

1. **Install Git LFS** (if not already installed):
   https://git-lfs.github.com

2. **Clone the repository (recommended):**
   ```bash
   git clone --recurse-submodules <repo-url>
   ```

3. **Fetch LFS files after cloning (if needed):**
   ```bash
   git lfs pull
   ```

> All collaborators must have Git LFS installed to access large files.
