from collections import Counter
import os
from typing import List, Union
from rdkit import Chem
from rdkit import RDLogger
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
from itertools import combinations
import networkx as nx
from time import time

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')
BACKBONE_ATOM_TYPES = {'C', 'H', 'O', 'S', 'N', 'P', 'Si', 'B'}
SUBSTITUENT_ATOM_TYPES = {'C', 'H', 'O', 'S', 'N', 'P', 'Si', 'B', 'F', 'Cl', 'Br', 'I'}
SUBSTITUENT_DOUBLE_BOND_PLACEHOLDER = 'Se'
SUBSTITUENT_PLACEHOLDER = 'F'
FUSED_RING_PLACEHOLDER = 'Cl'

MAX_RINGS_ATOMS = 8


def extract_induced_submol(mol, atom_indices, ring_mode=False):
    atom_indices_set = set(atom_indices)
    em = Chem.EditableMol(Chem.Mol())
    # if ring mode also keep placeholders that are connected to one of the atoms in atom_indices
    if ring_mode:
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in atom_indices_set and atom.GetSymbol() in [SUBSTITUENT_PLACEHOLDER, SUBSTITUENT_DOUBLE_BOND_PLACEHOLDER]:
                # check if connected to any atom in atom_indices
                neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                if any(nbr in atom_indices_set for nbr in neighbors):
                    atom_indices_set.add(idx)
    
    # Add atoms from atom_indices
    old_to_new_idx = {}
    for old_idx in atom_indices_set:
        atom = mol.GetAtomWithIdx(old_idx)
        new_idx = em.AddAtom(atom)
        old_to_new_idx[old_idx] = new_idx
    
    # Add bonds between atoms in atom_indices and cap external bonds with F
    def get_placeholder_atom(bond):
        if ring_mode:
            return FUSED_RING_PLACEHOLDER, Chem.BondType.SINGLE
        else:
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                return SUBSTITUENT_DOUBLE_BOND_PLACEHOLDER, Chem.BondType.DOUBLE
            else:
                return SUBSTITUENT_PLACEHOLDER, Chem.BondType.SINGLE
            
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx in atom_indices_set and end_idx in atom_indices_set:
            # Both atoms in set, add the bond
            em.AddBond(old_to_new_idx[begin_idx], old_to_new_idx[end_idx], bond.GetBondType())
        elif begin_idx in atom_indices_set and end_idx not in atom_indices_set:
            placeholder, bond_type = get_placeholder_atom(bond) 
            # if placeholder == SUBSTITUENT_PLACEHOLDER: continue  # skip 
            f_idx = em.AddAtom(Chem.Atom(placeholder))
            em.AddBond(old_to_new_idx[begin_idx], f_idx, bond_type)
        elif end_idx in atom_indices_set and begin_idx not in atom_indices_set:
            placeholder, bond_type = get_placeholder_atom(bond)
            # if placeholder == SUBSTITUENT_PLACEHOLDER: continue  # skip 
            f_idx = em.AddAtom(Chem.Atom(placeholder))
            em.AddBond(old_to_new_idx[end_idx], f_idx, bond_type)
    if ring_mode:
        # replace any SUBSTITUENT_PLACEHOLDER with H
        rw = Chem.RWMol(em.GetMol())
        for atom in rw.GetAtoms():
            if atom.GetSymbol() == SUBSTITUENT_PLACEHOLDER:
                atom.SetAtomicNum(1)  # change to H
        submol = rw.GetMol()
    else:
        submol = em.GetMol()
    
    return submol

def extract_backbone_inchis(mol: Union[Chem.Mol, str], min_backbone_rings=2):    
    if isinstance(mol, str):
        string_rep = mol
        if mol.startswith("InChI="):
            mol = Chem.MolFromInchi(string_rep, sanitize=False)
        else:
            mol = Chem.MolFromSmiles(string_rep, sanitize=False)
    else:
        string_rep = 'unknown'
    if mol is None:
        print(f"Failed to parse InChI: {mol}")
        return []
    unsanitize_mol = Chem.Mol(mol)
    if mol is None:
        return []

    try:
        Chem.SanitizeMol(mol)
    except:
        return []

    ring_info = mol.GetRingInfo()
    aromatic_rings = [
        set(r) for r in ring_info.AtomRings()
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
        and len(r) <= MAX_RINGS_ATOMS
    ]
    if not aromatic_rings:
        return []

    # find all fused aromatic clusters
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

    backbones = []
    for cluster_atoms in fused_clusters:
        ring_count = sum(1 for r in aromatic_rings if r.issubset(cluster_atoms))
        if ring_count >= min_backbone_rings:
            submol = extract_induced_submol(unsanitize_mol, cluster_atoms, ring_mode=False)

            inchi = Chem.MolToInchi(submol)
            backbones.append(inchi)

    return backbones

def get_all_motifs(mol: Union[Chem.Mol, str], plot: bool = False, plot_prefix: str = "") -> List[str]:
    """
    Extract all possible aromatic motifs from a molecule.
    
    Args:
        inchi: InChI string of the molecule
        plot: If True, plot the original molecule and all submolecules
    
    Returns:
        List of InChI strings for all motifs
    """
    
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
    aromatic_rings = [set(ring) for ring in ring_info.AtomRings()] # in this stage all acount as romatic
    
    if len(aromatic_rings) == 0:
        return []
    
    # Build adjacency list for ring graph (rings are fused if they share >=2 atoms)
    n_rings = len(aromatic_rings)
    G = nx.Graph()
    G.add_nodes_from(range(n_rings))
    
    for i in range(n_rings):
        for j in range(i + 1, n_rings):
            if len(aromatic_rings[i] & aromatic_rings[j]) >= 2:
                G.add_edge(i, j)
    
    # 3. Find all connected subgraphs with at least 2 nodes using NetworkX
    def connected_subgraphs(G):
        nodes = sorted(G.nodes())

        def expand(subgraph_nodes, frontier, root):
            yield subgraph_nodes.copy()

            for v in list(frontier):
                new_subgraph = subgraph_nodes | {v}
                new_frontier = frontier | set(G.neighbors(v))
                new_frontier -= new_subgraph

                # ordering constraint avoids duplicates
                if v > root:
                    yield from expand(new_subgraph, new_frontier, root)

                frontier.remove(v)

        for root in nodes:
            subgraph = {root}
            frontier = set(G.neighbors(root))
            yield from expand(subgraph, frontier, root)

    connected_components = [g for g in connected_subgraphs(G) if len(g) >=1]
    
    motif_inchis = []
    for i, ring_subset in enumerate(connected_components):
        atoms_in_motif = set()
        for ring_idx in ring_subset:
            atoms_in_motif |= aromatic_rings[ring_idx]
        
        submol = extract_induced_submol(unsanitize_mol, atoms_in_motif, ring_mode=True)
        
        try:
            motif_inchis.append(Chem.MolToInchi(submol))
        except Exception as e:
            print(f"Motif {i+1}: Direct InChI failed - {e}")
    
    # Remove duplicates to get unique substructures
    unique_motifs = list(set(motif_inchis))
    if len(connected_components) != len(motif_inchis):
        print(f"Warning: Expected {len(connected_components)} motifs, but got {len(motif_inchis)} after InChI conversion for molecule {inchi}")

    if plot and unique_motifs:
        from rdkit.Chem import Draw
        import matplotlib.pyplot as plt
        out_folder = "output_plots"
        os.makedirs(out_folder, exist_ok=True)
        unique_motifs = sorted(list(unique_motifs), key=lambda x: len(x))  # sort by length of InChI for better visualization
        motif_count = Counter(motif_inchis)
        # Convert motif InChIs back to molecules for plotting
        motif_mols = []
        for motif_inchi in unique_motifs:
            motif_mol = Chem.MolFromInchi(motif_inchi)
            if motif_mol:
                motif_mols.append(motif_mol)
        
        # Plot original molecule
        img = Draw.MolToImage(mol)
        plt.figure(figsize=(4, 4))
        plt.imshow(img)
        plt.title("Original Molecule")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(f'{out_folder}/{plot_prefix}orig.jpg')
        plt.close()
        
        # Plot all motifs in a grid
        if motif_mols:
            img = Draw.MolsToGridImage(motif_mols, molsPerRow=4, subImgSize=(300, 300), legends=[str(motif_count[inchi]) for inchi in unique_motifs])
            plt.figure(figsize=(12, 3 * ((len(motif_mols) - 1) // 4 + 1)))
            plt.imshow(img)
            plt.title(f"All Motifs ({len(motif_mols)} total)")
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(f'{out_folder}/{plot_prefix}motifs.jpg')
            plt.close()
    
    return motif_inchis

def plot_mol(mol: Chem.Mol, title: str="", filename: str="output.jpg"):
    from rdkit.Chem import Draw
    import matplotlib.pyplot as plt

    img = Draw.MolToImage(mol)
    plt.figure(figsize=(4, 4))
    plt.imshow(img)
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

if __name__ == '__main__':    
    # s = time()
    # get_all_motifs('InChI=1S/C15H9N3/c1-2-6-11-10(5-1)9-14-15(17-18-16-14)13-8-4-3-7-12(11)13/h1-9H', plot=True)

    # # get_all_motifs('InChI=1S/C84H32/c1-2-6-34-33(5-1)37-9-41-42(10-38(34)37)46-14-50-49(13-45(41)46)53-17-57-58(18-54(50)53)62-22-66-65(21-61(57)62)69-25-73-74(26-70(66)69)78-30-82-81(29-77(73)78)83-31-79-75-27-71-67-23-63-59-19-55-51-15-47-43-11-39-35-7-3-4-8-36(35)40(39)12-44(43)48(47)16-52(51)56(55)20-60(59)64(63)24-68(67)72(71)28-76(75)80(79)32-84(82)83/h1-32H', plot=True)
    # print(f"Time taken: {time() - s} seconds")

    smiles = 'CC1=CC(=O)OC2=C1C=C3C=C(OC3=C2C)C'
    smiles = 'C[C@@]12[C@](C[C@@H](O1)N3C4C=CC=CC4=C5C3=C6N2C7=CC=CCC7=C6C8=C5C(=O)NC8)(C(=O)OC)O'
    # smiles = 'C12=C(C(=C(C(=C1Cl)Cl)Cl)Cl)C3=NC4=NC(=NC5=NC(=NC6=C7C(=C(N6)NC2=N3)C(=C(C(=C7Cl)Cl)Cl)Cl)C8=C5C(=C(C(=C8Cl)Cl)Cl)Cl)C9=C4C(=C(C(=C9Cl)Cl)Cl)Cl'
    # plot_mol(Chem.MolFromSmiles(smiles), filename="output_plots/test_mol.jpg")
    smiles = "InChI=1S/C22H8F2N4Se2/c23-13-15-17(27-19(25-15)9-5-1-3-7-11(9)21(27)29)14(24)18-16(13)26-20-10-6-2-4-8-12(10)22(30)28(18)20/h1-8H"
    inchis = extract_backbone_inchis(smiles)
    for i, inchi in enumerate(inchis):
        get_all_motifs(inchi, plot=True)

    # #read first line from PAS-CID-SMILES.txt
    # with open("data_v3/PAS-CID-SMILES.txt", "r") as f:
    #     line = f.readline()
    #     cid, smiles = line.strip().split("\t")
    # print(f"Processing CID: {cid}, SMILES: {smiles}")
    # inchis = extract_backbone_inchis(smiles)
    # print(f"Extracted {len(inchis)} backbone InChIs:")
    # plot_mol(Chem.MolFromSmiles(smiles), filename="output_plots/sample_mol.jpg")
    # for i, inchi in enumerate(inchis):
    #     plot_mol(Chem.MolFromInchi(inchi), filename=f"output_plots/backbone_{i+1}.jpg")
    #     get_all_motifs(inchi, plot=True)

    # smiles = [        
    #     "C(C1CCCC2CCCCC12)C1CCCC2CCCCC12", # test two fused 2 benzene rings
    #     "C(C1CCCC2COCCC12)C1CCCC2COCCC12",  # test two fused benzene with O
    #     "C(C1CCCC2CCCCC12)C1CCCC2COCCC12"  # test fused benzene and benzene with O

    #     ]
    # for backbone_smiles in smiles:
    #     mol = Chem.MolFromSmiles(backbone_smiles)
    #     plot_mol(mol, filename="output_plots/test_mol.jpg")
    #     inchi = Chem.MolToInchi(mol)
    #     inchis = extract_backbone_inchis(inchi)
    #     inchis = get_all_motifs(inchi, plot=True)