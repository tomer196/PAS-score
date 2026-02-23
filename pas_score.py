
from collections import defaultdict
import pickle
from typing import Union
from rdkit import Chem

from extract_motifs import get_all_motifs, extract_backbone_inchis, plot_mol

def n_rings(inchi: str) -> int:
    mol = Chem.MolFromInchi(inchi)
    num_rings = mol.GetRingInfo().NumRings()
    return num_rings

class PASScore:
    min_value = -0.5
    max_value = 3
    average_per_n_rings = True

    def __init__(self, plot=False):
        target_file = "CATACONDENSED-MOTIFS-INCHI-PAS-SCORES.pkl"
        with open(target_file, "rb") as f:
            self.pas_score = pickle.load(f)
        self.default_missing_score = min(self.pas_score.values()) - 0.5 
        self.plot = plot

    def normalize_score(self, score: float):
        # convert from range [min_value, max_value] to [0, 10] in opposite direction
        norm_score = (self.max_value - score) / (self.max_value - self.min_value) * 10.0
        return norm_score

    def __call__(self, mol: Union[Chem.Mol, str]) -> float:
        if isinstance(mol, str):
            if mol.startswith("InChI="):
                mol = Chem.MolFromInchi(mol, sanitize=False)
            else:
                mol = Chem.MolFromSmiles(mol, sanitize=False)
        if mol is None:
            raise ValueError("Invalid molecule input.")
        backbone_inchis = extract_backbone_inchis(mol, min_backbone_rings=1)
        motifs = []
        for i, backbone_inchi in enumerate(backbone_inchis):   
            motifs.extend(get_all_motifs(backbone_inchi, plot=self.plot, plot_prefix=str(i)))
        if len(motifs) == 0:
            return None  
        # score = 0
        scores_per_n_rings = defaultdict(list)
        for motif_inchi in motifs:
            motif_score = self.pas_score.get(motif_inchi, self.default_missing_score)
            # score += motif_score
            nr = n_rings(motif_inchi)
            scores_per_n_rings[nr].append(motif_score)
        # scores_per_n_rings_filtered = {k: v for k, v in scores_per_n_rings.items() if k >= 5}
        scores_per_n_rings_filtered = scores_per_n_rings

        if self.average_per_n_rings:
            for nr, scores in scores_per_n_rings_filtered.items():
                scores_per_n_rings_filtered[nr] = sum(scores) / len(scores)
            score = sum(scores_per_n_rings_filtered.values()) / len(scores_per_n_rings_filtered)
        else:
            score = sum([sum(scores) for scores in scores_per_n_rings_filtered.values()])
            score /= len(motifs)  # normalize by number of motifs to avoid bias towards larger molecules with more motifs

        # return score
        return self.normalize_score(score) 

if __name__ == '__main__':
    pas = PASScore(plot=False)
    test_mols = [
        "InChI=1S/C84H32/c1-2-6-34-33(5-1)37-9-41-42(10-38(34)37)46-14-50-49(13-45(41)46)53-17-57-58(18-54(50)53)62-22-66-65(21-61(57)62)69-25-73-74(26-70(66)69)78-30-82-81(29-77(73)78)83-31-79-75-27-71-67-23-63-59-19-55-51-15-47-43-11-39-35-7-3-4-8-36(35)40(39)12-44(43)48(47)16-52(51)56(55)20-60(59)64(63)24-68(67)72(71)28-76(75)80(79)32-84(82)83/h1-32H",
        "InChI=1S/C12H6F3N/c13-7-1-3-11-9(5-7)10-6-8(14)2-4-12(10)16(11)15/h1-6H",
        "InChI=1S/C12H6F3N/c13-9-5-8-7-3-1-2-4-11(7)16(15)12(8)6-10(9)14/h1-6H",
        "InChI=1S/C14F10/c15-5-1-2(8(18)12(22)11(21)7(1)17)6(16)4-3(5)9(19)13(23)14(24)10(4)20",
        "InChI=1S/C16H9FO/c17-13-6-3-7-14-15(13)12-9-8-10-4-1-2-5-11(10)16(12)18-14/h1-9H",
        "CC(=O)OC1=CC=CC=C1C(=O)O", # Aspirin
        "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCSC(=O)CC(=O)O)O",
    ]
    test_mols = [
"C1=CC=CC2=C1[N](C3=C2C=CC7=C3C4=C(C5=C(C=C4)[N](C6=CC=CC=C56)[H])S7)[H]",
    ]
    # for mol in test_mols:
    #     score = pas(mol)
    #     if score is None:
    #         print(f"Mol: {mol}, PAS Score: N/A (no motifs found)")
    #     else:
    #         print(f"Mol: {mol}, PAS Score: {score:.4f}")    


    test_mols = {
        # 1: "InChI=1S/C18H8F3NO/c19-12-6-3-5-10-15-16-11(8-13(20)18(15)23-17(10)12)9-4-1-2-7-14(9)22(16)21/h1-8H",
        # 2: "InChI=1S/C20H11FO/c21-19-15-8-4-3-7-14(15)10-17-16-9-12-5-1-2-6-13(12)11-18(16)22-20(17)19/h1-11H",
        # 3: "InChI=1S/C18H9FOS/c19-13-9-12-10-5-1-3-7-14(10)20-17(12)18-16(13)11-6-2-4-8-15(11)21-18/h1-9H",
        # 4: "InChI=1S/C20H6F6N2/c21-7-1-3-9-10-4-2-8(22)6-12(10)18-17(11(9)5-7)27-19-15(25)13(23)14(24)16(26)20(19)28-18/h1-6H",
        # 5: "InChI=1S/C19H10F2N2/c20-15-10-16(21)23-18(15)17-12-5-2-1-4-11(12)7-8-13(17)14-6-3-9-22-19(14)23/h1-10H",
        # 6: "InChI=1S/C24H10F5N3/c25-13-3-1-5-15-19(13)11-7-9-17-21(23(11)31(15)28)22-18(30(17)27)10-8-12-20-14(26)4-2-6-16(20)32(29)24(12)22/h1-10H",
        # 7: "InChI=1S/C19H9F2N3Se2/c20-23-15-7-8-22-9-14(15)11-5-6-13-16(17(11)23)10-3-1-2-4-12(10)18(25)24(21)19(13)26/h1-9H",
        # 8: "InChI=1S/C18H9FN4O/c19-23-13-5-2-1-4-12(13)22-16-11(21-18(22)23)8-7-10-15-14(24-17(10)16)6-3-9-20-15/h1-9H",
        # 9: "InChI=1S/C21H11F2N3/c22-12-9-10-17-16(11-12)13-5-1-2-6-14(13)20-24-21-19(25(17)20)15-7-3-4-8-18(15)26(21)23/h1-11H",
        # 10: "InChI=1S/C20H8F2O2Se2/c21-9-2-5-14-13(7-9)19(25)12-4-6-15-17(18(12)24-14)20(26)11-3-1-10(22)8-16(11)23-15/h1-8H", 
        # 11: "InChI=1S/C22H11FS2/c23-15-9-5-11-17-20(15)18-12-6-1-2-7-13(12)21-19(22(18)25-17)14-8-3-4-10-16(14)24-21/h1-11H",
        # 12: "InChI=1S/C22H12FNO/c23-12-9-10-15-17(11-12)24-21-19(15)13-5-1-2-6-14(13)20-16-7-3-4-8-18(16)25-22(20)21/h1-11,24H",
        # 13: "InChI=1S/C21H12FN/c22-19-11-18-20-14(10-9-13-5-1-2-6-15(13)20)12-23-21(18)17-8-4-3-7-16(17)19/h1-12H",
        # 14: "InChI=1S/C24H12F3N/c25-15-6-8-22-19(11-15)20-12-18-14(10-23(20)28(22)27)5-7-17-21(26)9-13-3-1-2-4-16(13)24(17)18/h1-12H",
        # 15: "InChI=1S/C24H14N2S/c1-3-7-17-13(5-1)14-10-12-20-22(23(14)26-17)16-9-11-19-21(24(16)27-20)15-6-2-4-8-18(15)25-19/h1-12,25-26H",
        # 16: "InChI=1S/C23H10F3N3S/c24-11-1-4-16-14(9-11)20-18(28(16)25)5-2-13-21-19(30-23(13)20)6-3-12-15-10-27-8-7-17(15)29(26)22(12)21/h1-10H",
        # 17: "InChI=1S/C22H12N2S/c1-3-7-15-13(5-1)9-11-17-18-12-10-14-6-2-4-8-16(14)20(18)22-21(19(15)17)23-25-24-22/h1-12H",
        # 18: "InChI=1S/C25H11F3N4/c26-16-9-10-29-18-11-17(27)22-20(19(16)18)24(28)32-23-15-8-4-2-6-13(15)12-5-1-3-7-14(12)21(23)30-25(32)31-22/h1-11H",
        19: "InChI=1S/C22H8F2N4Se2/c23-13-15-17(27-19(25-15)9-5-1-3-7-11(9)21(27)29)14(24)18-16(13)26-20-10-6-2-4-8-12(10)22(30)28(18)20/h1-8H", 
        20: "InChI=1S/C20H11FN2Se/c21-19-18-15(14-7-3-4-8-16(14)22-19)11-17-13-6-2-1-5-12(13)9-10-23(17)20(18)24/h1-11H"
    }
    scores = []
    for key, mol in test_mols.items():
        score = pas(mol)
        if score is None:
            print(f"Mol: {key}, PAS Score: N/A (no motifs found)")
        else:
            print(f"Mol: {key}, PAS Score: {score:.4f}")    
        scores.append(score.item())
    print(",".join(map(str, scores)))

