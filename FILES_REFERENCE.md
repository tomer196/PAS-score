# PAS Motif-Based Score: File-by-File Reference

This document describes the purpose and main functions of each script/module in the PASscore project.

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
