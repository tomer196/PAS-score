
import os, sys
from typing import Union
from rdkit import Chem
from rdkit import RDLogger
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
from collections import Counter
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from extract_motifs import extract_induced_submol


data_folder = "data_v3/"
MAX_RINGS_ATOMS = 8
ring_str2mol = {}

def canonical_bracelet(chars):
    s = list(chars)
    n = len(s)

    def min_rotation(seq):
        # lexicographically smallest rotation (O(n^2), fine for small rings)
        return min(seq[i:] + seq[:i] for i in range(n))

    a = min_rotation(s)
    b = min_rotation(list(reversed(s)))
    canon = min(a, b)
    return "".join(canon)

# print(canonical_bracelet(['C','C','C','N','C','C']))          # same as any rotation/mirror
# print(canonical_bracelet(['C','N','C','C','C','C']))          # same as any rotation/mirror
# print(canonical_bracelet(['C','O','C','N','C','C']))  
# print(canonical_bracelet(['C','C','C','N','C','O']))  

def count_string(chars):
    c = Counter(chars)
    return "".join(f"{ch}{c[ch]}" for ch in sorted(c))

def get_all_rings(mol: Union[Chem.Mol, str]):
    # 1. Get molecule from InChI
    if isinstance(mol, str):
        inchi = mol
        mol = Chem.MolFromInchi(inchi, sanitize=False)
    else:
        inchi = 'unknown'
    if mol is None:
        print(f"Failed to parse InChI: {mol}")
        return []
    unsanitize_mol = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(mol)
    except:
        print(f"Failed to sanitize molecule from InChI: {mol}")
        return []
    
    # 2. Find all aromatic rings and build graph
    ring_info = mol.GetRingInfo()
    aromatic_rings = [set(ring) for ring in ring_info.AtomRings() if len(ring) <= MAX_RINGS_ATOMS] # in this stage all acount as romatic
    
    aromatic_rings_atoms = [[mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring] for ring in aromatic_rings]
    ring_strs = [canonical_bracelet(atoms) for atoms in aromatic_rings_atoms]
    for ring_str, atoms in zip(ring_strs, aromatic_rings):
        if ring_str not in ring_str2mol:
            mol = extract_induced_submol(unsanitize_mol, atoms)
            # remove any 'F' or 'Se' from the mol
            rw = Chem.RWMol(mol)
            for atom in rw.GetAtoms():
                if atom.GetSymbol() in ['F', 'Se']:
                    atom.SetAtomicNum(1)  # change to H
            submol = rw.GetMol()
            submol = Chem.RemoveHs(submol)
            ring_str2mol[ring_str] = submol
    return ring_strs
    

def process_line_count_rings(line):
    
    parts = line.strip().split("\t")
    if len(parts) != 2:
        return []
    
    count_str, inchi = parts
    try:
        count = int(count_str)
    except ValueError:
        return []
    
    motifs = get_all_rings(inchi)
    return [(motif, count) for motif in motifs]

def count_rings(plot=True):
    source_file = data_folder + "CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"
    target_file = data_folder + "RING-COUNTS.txt"
    
    # Get total lines for progress bar
    with open(source_file) as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines to process: {total_lines}")

    rings_count = Counter()

    # with Pool(cpu_count()) as pool:
    #     ring_lists = list(tqdm(pool.imap(process_line_count_rings, open(source_file)), total=total_lines, desc="Extracting motifs"))
    #     for ring_list in ring_lists:
    #         for ring_str, count in ring_list:
    #             rings_count[ring_str] += count
    
    for line in tqdm(open(source_file), total=total_lines, desc="Extracting motifs"):
        ring_list = process_line_count_rings(line)
        for ring_str, count in ring_list:
            rings_count[ring_str] += count
    
    print(f"Found {len(rings_count)} unique ring structures")
    
    # Write to file
    with open(target_file, "w") as out_f:
        for ring_str, count in rings_count.items():
            out_f.write(f"{count}\t{ring_str}\n")
    
    print(f"Saved motif counts to {target_file}")
    if plot:
        import matplotlib.pyplot as plt
        from rdkit.Chem import Draw
        rings_mols = []
        ring_strs = list(rings_count.keys())
        for ring_str in ring_strs:
            mol = ring_str2mol[ring_str]
            rings_mols.append(mol)

        img = Draw.MolsToGridImage(rings_mols, molsPerRow=4, subImgSize=(300, 300), legends=[f"{rings_count[ring_str]:,}" for ring_str in ring_strs])
        plt.figure(figsize=(12, 3 * ((len(rings_mols) - 1) // 4 + 1)))
        plt.imshow(img)
        plt.title(f"All Motifs ({len(rings_mols)} total)")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(f'output_plots/ring_counts.jpg')
        plt.close()


if __name__ == '__main__':
    count_rings()   