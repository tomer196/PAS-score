import pickle
import shutil, sys
from typing import Union
from rdkit import Chem
from rdkit.Chem import Draw, RDConfig
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
from scipy.stats import pearsonr, spearmanr

from extract_motifs import get_all_motifs
from pas_score import PASScore
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


def calculate_all_backbone_pas_scores(output_file):
    """Calculate PAS scores for all molecules in CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"""
    pas = PASScore()
    file = pas.data_folder + "CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"
    failed_inchis_file = pas.data_folder + "BACKBONE-PAS-SCORE-FAILED-INCHIS.txt"
    
    inchis = []
    scores = []
    failed = []
    
    # First pass: count total lines
    print("Loading molecules...")
    with open(file, "r") as f:
        total_lines = sum(1 for _ in f)
    
    # Second pass: calculate scores with progress bar
    with open(file, "r") as f:
        for line in tqdm(f, total=total_lines, desc="Calculating PAS scores"):
            count, inchi = line.strip().split()
            inchis.append(inchi)
            
            try:
                score = pas(inchi)
                if score is None:
                    score = np.nan
                    failed.append(inchi)
            except Exception as e:
                score = np.nan
                failed.append(inchi)
                print(f"Failed to calculate PAS score for InChI: {inchi}, Error: {e}")
            scores.append(score)
    
    inchis_array = np.array(inchis)
    scores_array = np.array(scores)
    
    # Save to file
    np.savez(output_file, inchis=inchis_array, scores=scores_array)
    print(f"\nSaved PASscore results to {output_file}")
    if len(failed) > 0:
        with open(failed_inchis_file, "w") as f:
            for inchi in failed:
                f.write(inchi + "\n")
        print(f"Saved {len(failed)} failed InChIs to {failed_inchis_file}")
    
    return inchis_array, scores_array

def calculate_all_backbone_sa_scores(output_file):
    """Calculate PAS scores for all molecules in CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"""
    file = pas.data_folder + "CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"
    
    inchis = []
    scores = []
    
    # First pass: count total lines
    print("Loading molecules...")
    with open(file, "r") as f:
        total_lines = sum(1 for _ in f)
    
    # Second pass: calculate scores with progress bar
    with open(file, "r") as f:
        for line in tqdm(f, total=total_lines, desc="Calculating SA scores"):
            count, inchi = line.strip().split()
            inchis.append(inchi)
            
            try:
                score = sascorer.calculateScore(Chem.MolFromInchi(inchi))
                if score is None:
                    score = np.nan
            except Exception as e:
                score = np.nan
            
            scores.append(score)
    
    inchis_array = np.array(inchis)
    scores_array = np.array(scores)
    
    # Save to file
    np.savez(output_file, inchis=inchis_array, scores=scores_array)
    print(f"\nSaved SAscore results to {output_file}")
    
    return inchis_array, scores_array


def plot_pas_score_distribution(scores=None, data_file=None, out_folder='output_plots'):
    """Plot PAS score distribution similar to plot_motifs_frequancy"""
    # Load from file if provided
    if data_file is not None:
        data = np.load(data_file)
        scores = data['scores']
    
    # Filter out NaN values
    valid_scores = scores[~np.isnan(scores)]
    
    if len(valid_scores) == 0:
        print("No valid PAS scores to plot")
        return
    
    # Sort scores in descending order
    sorted_scores = np.sort(valid_scores)[::-1]
    
    plt.figure(figsize=(10, 6))
    plt.plot(sorted_scores, marker='o', markersize=3)
    plt.xlabel('Molecule Rank')
    plt.ylabel('PAS Score')
    plt.title(f'PAS Scores Distribution for Catacondensed Backbone Molecules (n={len(valid_scores)})')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{out_folder}/backbone_pas_scores.jpg')
    print(f"Saved plot to {out_folder}/backbone_pas_scores.jpg")
    plt.close()
    
    # Print statistics
    print(f"\nPAS Score Statistics:")
    print(f"  Valid scores: {len(valid_scores)}")
    print(f"  Invalid (NaN) scores: {len(scores) - len(valid_scores)}")
    print(f"  Min score: {np.min(valid_scores):.4f}")
    print(f"  Max score: {np.max(valid_scores):.4f}")
    print(f"  Mean score: {np.mean(valid_scores):.4f}")
    print(f"  Median score: {np.median(valid_scores):.4f}")


def sample_molecules_by_score(inchis=None, scores=None, data_file=None, n_samples=20, min_rings=None, max_rings=None, out_folder='sampled_molecules'):
    """Sample molecules with equal spacing across the score range"""
    # Load from file if provided
    if data_file is not None:
        data = np.load(data_file)
        inchis = data['inchis']
        scores = data['scores']
    
    # Filter out NaN values
    valid_mask = ~np.isnan(scores)
    valid_inchis = inchis[valid_mask]
    valid_scores = scores[valid_mask]
    
    if len(valid_scores) == 0:
        print("No valid scores to sample from")
        return
    
    # Filter by ring count if specified
    if min_rings is not None or max_rings is not None:
        print(f"Filtering molecules by ring count (min: {min_rings}, max: {max_rings})...")
        filtered_inchis = []
        filtered_scores = []
        
        for inchi, score in tqdm(zip(valid_inchis, valid_scores), total=len(valid_inchis), desc="Filtering by rings"):
            try:
                mol = Chem.MolFromInchi(inchi)
                if mol is not None:
                    num_rings = mol.GetRingInfo().NumRings()
                    if (min_rings is None or num_rings >= min_rings) and \
                       (max_rings is None or num_rings <= max_rings):
                        filtered_inchis.append(inchi)
                        filtered_scores.append(score)
            except:
                pass
        
        valid_inchis = np.array(filtered_inchis)
        valid_scores = np.array(filtered_scores)
        print(f"After filtering: {len(valid_scores)} molecules")
        
        if len(valid_scores) == 0:
            print("No molecules match the ring count criteria")
            return
    
    # Sort by score
    sorted_indices = np.argsort(valid_scores)
    sorted_inchis = valid_inchis[sorted_indices]
    sorted_scores = valid_scores[sorted_indices]
    
    # Sample with equal spacing
    indices = np.linspace(0, len(sorted_scores) - 1, n_samples+2, dtype=int)
    indices = indices[1:-1]  # Exclude first and last to avoid extremes
    sampled_inchis = sorted_inchis[indices]
    sampled_scores = sorted_scores[indices]
    
    print(f"\n{'='*80}")
    print(f"Sampled {n_samples} molecules with equal spacing across score range:")
    print(f"{'='*80}")
    
    # Create directory for sampled molecules if it doesn't exist
    shutil.rmtree(out_folder, ignore_errors=True)
    os.makedirs(out_folder, exist_ok=True)
    
    csv_data = []
    for i, (inchi, score) in enumerate(zip(sampled_inchis, sampled_scores)):
        print(f"{i+1:2d}. Score: {score:7.4f} | InChI: {inchi}")
        
        # Convert InChI to molecule and draw it
        try:
            mol = Chem.MolFromInchi(inchi, treatWarningAsError=False)
            if mol is not None:
                # replace any SUBSTITUENT_PLACEHOLDER with H
                rw = Chem.RWMol(mol)
                for atom in rw.GetAtoms():
                    if atom.GetSymbol() == 'F':
                        atom.SetAtomicNum(1)  # change to H
                mol = rw.GetMol()
                mol = Chem.RemoveHs(mol)
                img = Draw.MolToImage(mol, size=(600, 600))
                sa_score = sascorer.calculateScore(mol)
                
                # Create figure with title
                fig, ax = plt.subplots(figsize=(8, 8))
                ax.imshow(img)
                ax.set_title(inchi, fontsize=8, wrap=True)
                ax.axis('off')
                plt.tight_layout()
                
                # Save figure
                filename = f'{out_folder}/molecule_{i+1:02d}_PASscore_{score:.4f}_SAscore_{sa_score:.4f}.png'
                plt.savefig(filename, dpi=150, bbox_inches='tight')
                plt.close()

                csv_data.append({
                    'molecule_number': i,
                    'inchi': inchi,
                    'pas_score': score,
                    'sa_score': sa_score
                })
            else:
                print(f"    Warning: Could not convert InChI to molecule")
        except Exception as e:
            print(f"    Error drawing molecule: {e}")

    df = pd.DataFrame(csv_data)
    csv_file = f'{out_folder}/molecule_scores.csv'
    df.to_csv(csv_file, index=False)
    print(f"Saved CSV with molecule scores to {csv_file}")
    
    print(f"{'='*80}")
    print(f"Saved {n_samples} molecule images to {out_folder} directory")
    
    return sampled_inchis, sampled_scores


def analyze_score_correlation(pas_file, sa_file, out_folder='output_plots'):
    # Load
    pas = np.load(pas_file, allow_pickle=True)
    sa = np.load(sa_file, allow_pickle=True)

    pas_dict = dict(zip(pas["inchis"], pas["scores"]))
    sa_dict = dict(zip(sa["inchis"], sa["scores"]))

    # Align by common inchis
    common_inchis = sorted(set(pas_dict) & set(sa_dict))
    pas_scores = np.array([pas_dict[i] for i in common_inchis])
    sa_scores = np.array([sa_dict[i] for i in common_inchis])
    finite_both = np.isfinite(pas_scores) & np.isfinite(sa_scores)
    pas_scores = pas_scores[finite_both]
    sa_scores = sa_scores[finite_both]
    common_inchis = [inchi for i, inchi in enumerate(common_inchis) if finite_both[i]]

    print(f"Common molecules: {len(common_inchis)}")

    # Correlations
    pearson = pearsonr(pas_scores, sa_scores)
    spearman = spearmanr(pas_scores, sa_scores)

    print(f"Pearson r:   {pearson.statistic:.4f} (p={pearson.pvalue:.2e})")
    print(f"Spearman ρ:  {spearman.statistic:.4f} (p={spearman.pvalue:.2e})")

    # Scatter plot
    plt.figure()
    plt.scatter(pas_scores, sa_scores, alpha=0.5)
    plt.xlabel("PAS score")
    plt.ylabel("SA score")
    plt.title("PAS vs SA scores")
    plt.savefig(f"{out_folder}/pas_vs_sa_scatter.jpg")

    # Rank comparison
    pas_rank = np.argsort(np.argsort(pas_scores))
    sa_rank = np.argsort(np.argsort(sa_scores))

    plt.figure()
    plt.scatter(pas_rank, sa_rank, alpha=0.5)
    plt.xlabel("PAS rank")
    plt.ylabel("SA rank")
    plt.title("Rank correlation")
    plt.savefig(f"{out_folder}/rank_correlation_scatter.jpg")


if __name__ == '__main__':
    pas = PASScore()
    
    pas_file = f"{pas.data_folder}/BACKBONE-PAS-SCORES.npz"
    if os.path.exists(pas_file):
        print(f"Loading existing PAS scores from {pas_file}")
        print("To recalculate, delete the file and run again.\n")
    else:
        print("Calculating PAS scores for all molecules...")
        inchis, scores = calculate_all_backbone_pas_scores(pas_file)

    sa_file = f"{pas.data_folder}/BACKBONE-SA-SCORES.npz"
    if os.path.exists(sa_file):
        print(f"Loading existing SA scores from {sa_file}")
        print("To recalculate, delete the file and run again.\n")
    else:
        print("Calculating SA scores for all molecules...")
        # Calculate SA scores for all molecules
        inchis, scores = calculate_all_backbone_sa_scores(sa_file)
    
    # # Plot the distribution
    plot_pas_score_distribution(data_file=pas_file)
    
    # # Sample 20 molecules with equal spacing
    # sample_molecules_by_score(data_file=pas_file, n_samples=20, out_folder='sampled_molecules')
    # sample_molecules_by_score(data_file=pas_file, n_samples=20, min_rings=5, max_rings=7, out_folder='sampled_molecules2')

    # Analyze correlation between PAS and SA scores
    # analyze_score_correlation(pas_file, sa_file)