import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import os

data_folder = "data_v3/"
out_folder = "output_plots"
os.makedirs(out_folder, exist_ok=True)

def plot_top_motifs(top_n=10, n_rings=None):
    file = data_folder + "CATACONDENSED-MOTIFS-INCHI-COUNTS.txt"
    inchis = []
    counts = []
    motif_mols = []
    with open(file, "r") as f:
        for line in f:
            count, inchi = line.strip().split()
            motif_mol = Chem.MolFromInchi(inchi)
            if n_rings is not None:
                ri = motif_mol.GetRingInfo().NumRings()
                if ri < n_rings:
                    continue
            inchis.append(inchi)
            counts.append(int(count))
            motif_mols.append(motif_mol)
            if len(inchis) >= top_n:
                break

    # Plot all motifs in a grid
    img = Draw.MolsToGridImage(motif_mols, molsPerRow=4, subImgSize=(300, 300), legends=list(map(str, counts)))
    plt.figure(figsize=(12, 3 * ((len(motif_mols) - 1) // 4 + 1)))
    plt.imshow(img)
    plt.title(f"Top {top_n} Most Common Motifs")
    plt.axis('off')
    plt.tight_layout()
    if n_rings is not None:
        plt.savefig(f'{out_folder}/top_motifs_min_rings_{n_rings}.jpg')
    else:
        plt.savefig(f'{out_folder}/top_motifs.jpg')

def plot_motifs_frequancy():
    file = data_folder + "CATACONDENSED-MOTIFS-INCHI-COUNTS.txt"
    counts = []
    with open(file, "r") as f:
        for line in f:
            count, inchi = line.strip().split()
            counts.append(int(count))

    counts = np.array(counts)
    sorted_counts = np.sort(counts)[::-1]
    plt.figure(figsize=(10, 6))
    plt.plot(sorted_counts, marker='o', markersize=3)
    plt.yscale('log')
    plt.xlabel('Motif Rank')
    plt.ylabel('Motif Count (log scale)')
    plt.title('Motif Frequencies in Dataset')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{out_folder}/motif_frequencies_log.jpg')

    score_ratio = 0.8
    sum_motifs = np.sum(sorted_counts)
    cumsum = np.cumsum(sorted_counts)
    threshold = sum_motifs * score_ratio
    num_fragments_80_percent = np.searchsorted(cumsum, threshold) + 1
    
    print(f"Total motifs: {len(sorted_counts)}")
    print(f"Number of fragment types forming 80% of all fragments: {num_fragments_80_percent}")
    
    # Calculate PAS score: log10(actual_count / num_fragments_80_percent)
    pas_score_per_motif = np.log10(sorted_counts / num_fragments_80_percent)
    # Plot PAS score per motif
    plt.figure(figsize=(10, 6))
    plt.plot(pas_score_per_motif, marker='o', markersize=3)
    plt.xlabel('Motif Rank')
    plt.ylabel('PAS Score')
    plt.title('PAS Score per Motif')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{out_folder}/pas_score_per_motif.jpg')

def analyze():
    file = data_folder + "CATACONDENSED-MOTIFS-INCHI-COUNTS.txt"
    counts = []
    with open(file, "r") as f:
        for line in f:
            count, inchi = line.strip().split()
            counts.append(int(count))

    total_motifs = len(counts)
    for limit in [1, 10, 100, 1000, 10000, 100000]:
        num_above_limit = sum(1 for c in counts if c > limit)
        fraction = num_above_limit / total_motifs
        print(f"Motifs with count > {limit:<7}: {num_above_limit:<5} ({fraction:.2%} of total)")

if __name__ == "__main__":
    plot_top_motifs()
    plot_top_motifs(n_rings=2)
    plot_top_motifs(n_rings=3)
    plot_top_motifs(n_rings=4)
    plot_motifs_frequancy()
    analyze()