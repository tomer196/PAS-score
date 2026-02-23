import os
from rdkit import Chem
from rdkit import RDLogger
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
from collections import Counter
import numpy as np
import pickle

from extract_motifs import get_all_motifs, extract_backbone_inchis

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')
BACKBONE_ATOM_TYPES = {'C', 'H', 'O', 'S', 'N', 'P', 'Si', 'B'}
SUBSTITUENT_ATOM_TYPES = {'C', 'H', 'O', 'S', 'N', 'P', 'Si', 'B', 'F', 'Cl', 'Br', 'I'}
MAX_RINGS_ATOMS = 8
data_folder = "data/"
os.makedirs(data_folder, exist_ok=True)

def is_pas_molecule(smiles: str, min_backbone_rings=2) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    # Ensure proper sanitization
    try:
        Chem.SanitizeMol(mol)
    except:
        return False
    if mol is None:
        return False

    # 1. Must contain carbon
    if not any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
        return False

    ring_info = mol.GetRingInfo()
    aromatic_rings = [
        set(ring) for ring in ring_info.AtomRings()
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        and len(ring) <= MAX_RINGS_ATOMS
    ]

    # 2. Backbone: at least 2 fused aromatic rings
    fused_clusters = []
    visited = set()

    for i, r1 in enumerate(aromatic_rings):
        if i in visited:
            continue
        cluster = set(r1)
        stack = [i]
        visited.add(i)

        while stack:
            k = stack.pop()
            for j, r2 in enumerate(aromatic_rings):
                if j not in visited and len(aromatic_rings[k] & r2) >= 2:
                    visited.add(j)
                    cluster |= r2
                    stack.append(j)

        fused_clusters.append(cluster)

    backbone_atoms = max(fused_clusters, key=len, default=None)
    if backbone_atoms is None or len(fused_clusters) == 0:
        return False

    # Require at least 2 rings in the backbone
    backbone_ring_count = sum(
        1 for r in aromatic_rings if r.issubset(backbone_atoms)
    )
    if backbone_ring_count < min_backbone_rings:
        return False

    # 3. Backbone atom types constraint
    for idx in backbone_atoms:
        if mol.GetAtomWithIdx(idx).GetSymbol() not in BACKBONE_ATOM_TYPES:
            return False

    # 4. Substituent atom types constraint
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in backbone_atoms:
            if atom.GetSymbol() not in SUBSTITUENT_ATOM_TYPES:
                return False

    return True

def process_line_pas(line):
    cid, smiles = line.strip().split("\t")
    if is_pas_molecule(smiles):
        return line
    return None

def filter_pas_from_pubchem():
    source_file = data_folder + "CID-SMILES"
    target_file = data_folder + "PAS-CID-SMILES.txt"
    # Get total lines for progress bar (fast single pass)
    with open(source_file) as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines to process: {total_lines}")

    # Stream processing - write results directly to file
    with open(source_file) as f, \
         open(target_file, "w") as out_f, \
         Pool(cpu_count()) as pool:
        
        results = pool.imap_unordered(process_line_pas, f, chunksize=1000)
        pas_count = 0
        
        for result in tqdm(results, total=total_lines, desc="Processing molecules"):
            if result:
                out_f.write(result)
                pas_count += 1

    print(f"Found {pas_count} PAS compounds")

def process_line_backbone(line, min_backbone_rings=2):
    cid, smiles = line.strip().split("\t")
    try:
        backbones = extract_backbone_inchis(smiles, min_backbone_rings)
    except:
        print(f"ERROR processing CID {cid} with SMILES {smiles}")
        # raise KeyError
        return None
    if not backbones:
        return None
    return cid + "\t" + "\t".join(backbones) + "\n"

def extract_backbones(parallel=True):
    source_file = data_folder + "PAS-CID-SMILES.txt"
    target_file = data_folder + "BACKBONE-INCHI.txt"
    # Get total lines for progress bar (fast single pass)
    with open(source_file) as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines to process: {total_lines}")

    if not parallel:
        # Non-parallel processing
        with open(source_file) as f, \
             open(target_file, "w") as out_f:
            
            for line in tqdm(f, total=total_lines, desc="Extracting backbones"):
                result = process_line_backbone(line, min_backbone_rings=2)
                if result:
                    out_f.write(result)
    else:
        # Parallel processing
        with open(source_file) as f, \
             open(target_file, "w") as out_f, \
             Pool(cpu_count()) as pool:

            process_func = partial(process_line_backbone, min_backbone_rings=2)
            results = pool.imap_unordered(
                process_func,
                f,
                chunksize=1000
            )

            for r in tqdm(results, total=total_lines, desc="Extracting backbones"):
                if r:
                    out_f.write(r)


def process_line_count_backbone(line):
    parts = line.strip().split("\t")
    if len(parts) < 2:
        return set()
    
    cid = parts[0]
    inchis_list = parts[1:]
    
    inchis = set(inchis_list)
    return inchis

def count_backbones():
    from collections import Counter
    
    source_file = data_folder + "BACKBONE-INCHI.txt"
    target_file = data_folder + "BACKBONE-INCHI-COUNTS.txt"
    
    # Get total lines for progress bar
    with open(source_file) as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines to process: {total_lines}")

    # Process all lines in parallel and count on the go
    inchi_counts = Counter()
    with open(source_file) as f, Pool(cpu_count()) as pool:
        results = pool.imap_unordered(
            process_line_count_backbone,
            f,
            chunksize=1000
        )
        
        for inchi_set in tqdm(results, total=total_lines, desc="Summing backbones"):
            inchi_counts.update(inchi_set)
    
    print(f"Found {len(inchi_counts)} unique backbone structures")
    
    # Write to file
    with open(target_file, "w") as out_f:
        for inchi, count in inchi_counts.most_common():
            out_f.write(f"{count}\t{inchi}\n")
    
    print(f"Saved counts to {target_file}")


def is_catacondensed(mol) -> bool:
    """Check if a molecule is catacondensed (no atom in more than 2 rings)"""
    if mol is None:
        return False
    
    try:
        Chem.SanitizeMol(mol)
    except:
        return False
    
    ring_info = mol.GetRingInfo()
    
    # Check each atom - if any atom is in 3+ rings, it's pericondensed
    for atom in mol.GetAtoms():
        if ring_info.NumAtomRings(atom.GetIdx()) >= 3:
            return False
    
    return True

def process_line_filter_catacondensed(line):
    parts = line.strip().split("\t")
    if len(parts) != 2:
        return None
    
    count, inchi = parts
    
    # Convert InChI to molecule and check if catacondensed
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is not None:
            if is_catacondensed(mol):
                return line
    except:
        pass
    
    return None

def filter_catacondensed_backbones():
    source_file = data_folder + "BACKBONE-INCHI-COUNTS.txt"
    target_file = data_folder + "CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"
    
    # Get total lines for progress bar
    with open(source_file) as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines to process: {total_lines}")

    # Process all lines in parallel
    catacondensed_count = 0
    with open(source_file) as f, \
         open(target_file, "w") as out_f, \
         Pool(cpu_count()) as pool:
        
        results = pool.imap_unordered(
            process_line_filter_catacondensed,
            f,
            chunksize=1000
        )
        
        for result in tqdm(results, total=total_lines, desc="Filtering catacondensed"):
            if result:
                out_f.write(result)
                catacondensed_count += 1
    
    print(f"Found {catacondensed_count} catacondensed backbones")
    print(f"Saved to {target_file}")

def process_line_extract_motifs(line):
    """Process a single line: extract motifs and return them with their weight"""
    import signal
    
    def timeout_handler(signum, frame):
        raise TimeoutError("Processing took too long")
    
    parts = line.strip().split("\t")
    if len(parts) != 2:
        return []
    
    count_str, inchi = parts
    try:
        count = int(count_str)
    except ValueError:
        return []
    
    # Set a timeout of 60 seconds per molecule
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(60)
    
    try:
        # Get all motifs for this backbone
        motifs = get_all_motifs(inchi)
        signal.alarm(0)  # Cancel the alarm
        
        # Return list of (motif_inchi, count) tuples
        return [(motif, count) for motif in motifs]
    except TimeoutError:
        signal.alarm(0)
        print(f"TIMEOUT: Skipping molecule with count={count}, InChI={inchi}")
        return []
    except Exception as e:
        signal.alarm(0)
        print(f"ERROR processing count={count}: {e}")
        return []

def sum_motifs_from_backbones():
    source_file = data_folder + "CATACONDENSED-BACKBONE-INCHI-COUNTS.txt"
    target_file = data_folder + "CATACONDENSED-MOTIFS-INCHI-COUNTS.txt"
    
    # Get total lines for progress bar
    with open(source_file) as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines to process: {total_lines}")

    # Process all lines in parallel and accumulate motif counts
    motif_counts = Counter()
    with open(source_file) as f, Pool(cpu_count()) as pool:
        results = pool.imap_unordered(
            process_line_extract_motifs,
            f,
            chunksize=100
        )
        
        for motif_list in tqdm(results, total=total_lines, desc="Extracting motifs"):
            for motif_inchi, count in motif_list:
                motif_counts[motif_inchi] += count
    
    print(f"Found {len(motif_counts)} unique motif structures")
    
    # Write to file
    with open(target_file, "w") as out_f:
        for inchi, count in motif_counts.most_common():
            out_f.write(f"{count}\t{inchi}\n")
    
    print(f"Saved motif counts to {target_file}")

def calc_pas_score_per_motif():
    file = data_folder + "CATACONDENSED-MOTIFS-INCHI-COUNTS.txt"
    target_file = data_folder + "CATACONDENSED-MOTIFS-INCHI-PAS-SCORES.pkl"
    counts = []
    inchis = []
    with open(file, "r") as f:
        for line in f:
            count, inchi = line.strip().split()
            counts.append(int(count))
            inchis.append(inchi)

    counts = np.array(counts)
    inchis = np.array(inchis)
    sorted_inchis = inchis[np.argsort(-counts)]
    sorted_counts = np.sort(counts)[::-1]
    # score_ratio = 0.8
    # sum_motifs = np.sum(sorted_counts)
    # cumsum = np.cumsum(sorted_counts)
    # threshold = sum_motifs * score_ratio
    # num_fragments_80_percent = np.searchsorted(cumsum, threshold) + 1
    
    print(f"Total motifs: {len(sorted_counts)}")
    # print(f"Number of fragment types forming 80% of all fragments: {num_fragments_80_percent}")
    
    # sa_score_per_motif = np.log10(counts / num_fragments_80_percent)
    sa_score_per_motif = np.log10(sorted_counts)
    sa_score = {inchi: score for inchi, score in zip(inchis, sa_score_per_motif)}
    with open(target_file, "wb") as f:
        pickle.dump(sa_score, f)
    print(f"Saved SA scores to {target_file}")

    #plot scores
    scores = [sa_score[inchi] for inchi in sorted_inchis]
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(scores, marker='o', markersize=3)
    plt.xlabel('Motif Rank')
    plt.ylabel('SA Score')
    plt.title('SA Score per Motif')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('sa_score_per_motif.jpg')

if __name__ == '__main__':
    filter_pas_from_pubchem()             # filter out only PAS save in PAS-CID-SMILES.txt
    extract_backbones()                   # extract backbones from PAS save in BACKBONE-INCHI.txt
    count_backbones()                     # count unique backbones save in BACKBONE-INCHI-COUNTS.txt
    filter_catacondensed_backbones()      # filter catacondensed backbones save in CATACONDENSED-BACKBONE-INCHI-COUNTS.txt
    sum_motifs_from_backbones()           # sum motifs from backbones save in CATACONDENSED-MOTIFS-INCHI-COUNTS.txt
    calc_pas_score_per_motif()               # calculate PAS score per motif and save in CATACONDENSED-MOTIFS-INCHI-PAS-SCORES.pkl