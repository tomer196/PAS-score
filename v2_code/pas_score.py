
import pickle
from typing import Union
from rdkit import Chem

from motif_based.v2_code.extract_motifs import get_all_motifs

class PASScore:
    data_folder = "data_v2/"
    def __init__(self):
        target_file = self.data_folder + "CATACONDENSED-MOTIFS-INCHI-PAS-SCORES.pkl"
        with open(target_file, "rb") as f:
            self.pas_score = pickle.load(f)
        self.default_missing_score = min(self.pas_score.values()) - 0.5 

    def __call__(self, mol: Union[Chem.Mol, str]) -> float:
        if isinstance(mol, str):
            if mol.startswith("InChI="):
                mol = Chem.MolFromInchi(mol)
            else:
                mol = Chem.MolFromSmiles(mol)
        if mol is None:
            raise ValueError("Invalid molecule input.")
        motifs = get_all_motifs(mol)
        if len(motifs) == 0:
            return None  
        score = 0
        for motif_inchi in motifs:
            score += self.pas_score.get(motif_inchi, self.default_missing_score)

        return score / len(motifs) 

if __name__ == '__main__':
    pas = PASScore()
    test_mols = [
        "InChI=1S/C84H32/c1-2-6-34-33(5-1)37-9-41-42(10-38(34)37)46-14-50-49(13-45(41)46)53-17-57-58(18-54(50)53)62-22-66-65(21-61(57)62)69-25-73-74(26-70(66)69)78-30-82-81(29-77(73)78)83-31-79-75-27-71-67-23-63-59-19-55-51-15-47-43-11-39-35-7-3-4-8-36(35)40(39)12-44(43)48(47)16-52(51)56(55)20-60(59)64(63)24-68(67)72(71)28-76(75)80(79)32-84(82)83/h1-32H",
        "InChI=1S/C18H30N2/c1-3-7-15-11(5-1)13-9-14-12-6-2-4-8-16(12)20-18(14)10-17(13)19-15/h11-20H,1-10H2",
        "InChI=1S/C20H32O/c1-3-7-15-13(5-1)10-12-18-19(15)17-11-9-14-6-2-4-8-16(14)20(17)21-18/h13-20H,1-12H2",
        "InChI=1S/C22H36/c1-3-7-18-15(5-1)11-14-21-20(18)13-12-17-10-9-16-6-2-4-8-19(16)22(17)21/h15-22H,1-14H2",
        "InChI=1S/C22H38N2/c1-3-7-17-11-21-19(9-15(17)5-1)13-23-22-12-18-8-4-2-6-16(18)10-20(22)14-24-21/h15-24H,1-14H2",
        "CC(=O)OC1=CC=CC=C1C(=O)O", # Aspirin
        "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCSC(=O)CC(=O)O)O",
        "C(C1CCCC2CCCCC12)C1CCCC2CCCCC12", # test two fused 2 benzene rings
        "C(C1CCCC2COCCC12)C1CCCC2COCCC12",  # test two fused benzene with O
        "C(C1CCCC2CCCCC12)C1CCCC2COCCC12"  # test fused benzene and benzene with O

    ]
    for mol in test_mols:
        score = pas(mol)
        if score is None:
            print(f"Mol: {mol}, PAS Score: N/A (no motifs found)")
        else:
            print(f"Mol: {mol}, PAS Score: {score:.4f}")    


    