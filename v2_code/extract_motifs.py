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
MAX_RINGS_ATOMS = 8


def sanitize_molecule(mol):
    """
    Comprehensively sanitize a molecule to remove radicals and fix bonding patterns.
    Makes consecutive double bonds alternate and saturates radicals with hydrogens.
    """
    if mol is None:
        return None
    
    try:
        # Step 1: Clear aromatic flags to work with explicit bonds
        rwmol = Chem.RWMol(mol)
        for atom in rwmol.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetNoImplicit(False)
        for bond in rwmol.GetBonds():
            bond.SetIsAromatic(False)
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.SetBondType(Chem.BondType.SINGLE)
        
        # Step 2: Convert triple bonds to double bonds
        for bond in rwmol.GetBonds():
            if bond.GetBondType() == Chem.BondType.TRIPLE:
                bond.SetBondType(Chem.BondType.DOUBLE)
        
        mol = rwmol.GetMol()
        
        # Step 3: Try to sanitize with explicit valence checking disabled
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except:
            print("1. faild SanitizeMol")
            pass
        
        # Step 4: Set implicit hydrogens to satisfy valences (this fixes radicals)
        rwmol = Chem.RWMol(mol)
        for atom in rwmol.GetAtoms():
            atom.SetNoImplicit(False)
            # Calculate how many hydrogens are needed
            try:
                atom.UpdatePropertyCache()
            except:
                print("2. faild UpdatePropertyCache")
                pass
        
        mol = rwmol.GetMol()
        
        # Step 5: Add explicit hydrogens to saturate everything
        try:
            mol = Chem.AddHs(mol, addCoords=False)
        except:
            print("3. faild AddHs")
            pass
        
        # Step 6: Final sanitization
        try:
            Chem.SanitizeMol(mol, catchErrors=True)
        except:
            print("4. faild SanitizeMol")
            pass
        
        # Step 7: Remove and re-add Hs to clean up
        try:
            mol = Chem.RemoveHs(mol, sanitize=False)
            mol = Chem.AddHs(mol)
        except:
            print("5. faild AddHs")
            pass
        
        return mol
        
    except Exception as e:
        print(f"Sanitization error: {e}")
        return mol

def extract_induced_submol(mol, atom_indices):
    em = Chem.EditableMol(Chem.Mol(mol))
    to_remove = sorted(
        [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in atom_indices],
        reverse=True
    )
    for idx in to_remove:
        em.RemoveAtom(idx)
    
    submol = em.GetMol()
    
    # Comprehensively sanitize the molecule
    # This will: remove radicals, fix consecutive double bonds, convert triple bonds, saturate with H
    submol = sanitize_molecule(submol)

    return submol

def get_all_motifs(mol: Union[Chem.Mol, str], plot: bool = False) -> List[str]:
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
        mol = Chem.MolFromInchi(mol)
    if mol is None:
        print(f"Failed to parse InChI: {mol}")
        return []
    
    try:
        Chem.SanitizeMol(mol)
    except:
        print(f"Failed to sanitize molecule from InChI: {mol}")
        return []
    
    # 2. Find all aromatic rings and build graph
    ring_info = mol.GetRingInfo()
    aromatic_rings = [set(ring) for ring in ring_info.AtomRings()] # in this stage all acount as romatic
    
    if len(aromatic_rings) < 2:
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

    connected_components = [g for g in connected_subgraphs(G) if len(g) >=2]
    # print(f"Found {len(connected_components)} connected aromatic subgraphs")
    
    motif_inchis = []
    for i, ring_subset in enumerate(connected_components):
        # 4. Extract atoms from these rings
        atoms_in_motif = set()
        for ring_idx in ring_subset:
            atoms_in_motif |= aromatic_rings[ring_idx]
        
        # Get submolecule
        submol = extract_induced_submol(mol, atoms_in_motif)
        
        # Try multiple strategies to convert to InChI
        motif_inchi = None
        
        # Strategy 1: Direct conversion
        try:
            motif_inchi = Chem.MolToInchi(submol)
        except Exception as e:
            print(f"Motif {i+1}: Direct InChI failed - {e}")
        
        # Strategy 2: Convert via SMILES first
        if not motif_inchi:
            try:
                smiles = Chem.MolToSmiles(submol)
                submol_fixed = Chem.MolFromSmiles(smiles)
                if submol_fixed:
                    motif_inchi = Chem.MolToInchi(submol_fixed)
                    print(f"Motif {i+1}: Success via SMILES conversion")
            except Exception as e:
                print(f"Motif {i+1}: SMILES route failed - {e}")
        
        # Strategy 3: Try with sanitize=False in SMILES
        if not motif_inchi:
            try:
                # Get SMILES without sanitization
                Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
                submol_copy = Chem.Mol(submol)
                smiles = Chem.MolToSmiles(submol_copy, kekuleSmiles=True)
                submol_fixed = Chem.MolFromSmiles(smiles, sanitize=True)
                if submol_fixed:
                    motif_inchi = Chem.MolToInchi(submol_fixed)
                    print(f"Motif {i+1}: Success via Kekule SMILES")
            except Exception as e:
                print(f"Motif {i+1}: Kekule SMILES route failed - {e}")
        
        if motif_inchi:
            motif_inchis.append(motif_inchi)
        else:
            print(f"Motif {i+1}: All conversion strategies failed")
    
    # Remove duplicates to get unique substructures
    unique_motifs = list(set(motif_inchis))
    if len(connected_components) != len(motif_inchis):
        print(f"Warning: Expected {len(connected_components)} motifs, but got {len(motif_inchis)} after InChI conversion for molecule {mol}")
    # if len(connected_components) > 100:
    #     print(f"Warning: Large number of motifs ({len(connected_components)}) for molecule {inchi}")
   
    
    # Plot if requested
    if plot and unique_motifs:
        from rdkit.Chem import Draw
        import matplotlib.pyplot as plt
        out_folder = "output_plots"
        os.makedirs(out_folder, exist_ok=True)
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
        plt.savefig(f'{out_folder}/orig.jpg')
        
        # Plot all motifs in a grid
        if motif_mols:
            img = Draw.MolsToGridImage(motif_mols, molsPerRow=4, subImgSize=(300, 300), legends=[str(motif_count[inchi]) for inchi in unique_motifs])
            plt.figure(figsize=(12, 3 * ((len(motif_mols) - 1) // 4 + 1)))
            plt.imshow(img)
            plt.title(f"All Motifs ({len(motif_mols)} total)")
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(f'{out_folder}/motifs.jpg')
    
    return motif_inchis

if __name__ == '__main__':
    # get_all_motifs('InChI=1S/C15H9N3/c1-2-6-11-10(5-1)9-14-15(17-18-16-14)13-8-4-3-7-12(11)13/h1-9H', plot=True)
    s = time()
    get_all_motifs('InChI=1S/C84H32/c1-2-6-34-33(5-1)37-9-41-42(10-38(34)37)46-14-50-49(13-45(41)46)53-17-57-58(18-54(50)53)62-22-66-65(21-61(57)62)69-25-73-74(26-70(66)69)78-30-82-81(29-77(73)78)83-31-79-75-27-71-67-23-63-59-19-55-51-15-47-43-11-39-35-7-3-4-8-36(35)40(39)12-44(43)48(47)16-52(51)56(55)20-60(59)64(63)24-68(67)72(71)28-76(75)80(79)32-84(82)83/h1-32H', plot=True)
    print(f"Time taken: {time() - s} seconds")
    # test('InChI=1S/C84H32/c1-2-6-34-33(5-1)37-9-41-42(10-38(34)37)46-14-50-49(13-45(41)46)53-17-57-58(18-54(50)53)62-22-66-65(21-61(57)62)69-25-73-74(26-70(66)69)78-30-82-81(29-77(73)78)83-31-79-75-27-71-67-23-63-59-19-55-51-15-47-43-11-39-35-7-3-4-8-36(35)40(39)12-44(43)48(47)16-52(51)56(55)20-60(59)64(63)24-68(67)72(71)28-76(75)80(79)32-84(82)83/h1-32H')