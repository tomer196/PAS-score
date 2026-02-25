
# Motif-based Synthetic Accessibility Score (PASscore)

PASscore is a synthetic accessibility metric for polycyclic aromatic systems (PASs) based on motif frequency analysis. It quantifies the synthetic accessibility of a molecule by decomposing it into ring-based motifs and scoring them according to their frequency in a large chemical database.

---

## Setup

1. **Clone the repository:**
   ```bash
   git clone git@github.com:tomer196/PAS-score.git
   cd PAS-score/motif_based
   ```

2. **Create a Python virtual environment:**
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Get the motif score data:**
   - Download the provided zip file (e.g., `CATACONDENSED-MOTIFS-INCHI-PAS-SCORES.zip`).
   - Unzip it into the `data/` directory:
     ```bash
     mkdir -p data
     unzip CATACONDENSED-MOTIFS-INCHI-PAS-SCORES.zip -d data/
     ```

---

## Inference: Scoring Molecules

To score a molecule using the PASscore metric:
```bash
python pas_score.py
```

---

## Training: Generating Motif Scores from Raw Data

If you want to regenerate the motif scores from raw PubChem data:

1. Download the raw PubChem data file `CID-SMILES.gz` from [here](https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz). Place the downloaded file in the `data/` directory.

2. Run the full pipeline:
   ```bash
   python extract_PASs_motifs.py
   ```
   This will process the raw data, extract motifs, and generate the motif score database.

---

## More Information

For detailed descriptions of each script and module, see [FILES_REFERENCE.md](FILES_REFERENCE.md).

