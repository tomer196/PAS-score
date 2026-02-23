from statistics import correlation
import os, sys
import pandas as pd
import rdkit.Chem as Chem

sys.path.insert(0, '/Users/tomer/private/chem/motif_based')
from other_metrics.RAscore.RAscore_NN import RAScorerNN
from other_metrics.RSPred.predictorRS import RSPredictor
from other_metrics.scscore.standalone_model_numpy import SCScorer
from pas_score import PASScore

sys.path.insert(0, '/Users/tomer/private/chem/motif_based/other_metrics/fsscore/src')
from fsscore.score import Scorer
from fsscore.models.ranknet import LitRankNet

sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, 'SA_Score'))
import sascorer

def get_passcore(inchis):
    pas_scorer = PASScore()
    scores = []
    for inchi in inchis:
        score = pas_scorer(inchi)
        scores.append(score)
    print(f"PAS scores calculated for {len(scores)} molecules")
    return scores

def get_sascore(smiles):
    scores = []
    for smi in smiles:
        score = sascorer.calculateScore(Chem.MolFromSmiles(smi))
        scores.append(score)
    print(f"SA scores calculated for {len(scores)} molecules")
    return scores

def get_rascore(smiles):
    model = RAScorerNN()
    scores = []
    for smi in smiles:
        scores.append(model.predict(smi))
    print(f"RA scores calculated for {len(scores)} molecules")
    return scores

def get_rscore(smiles):
    model = RSPredictor()
    scores = []
    for smi in smiles:
        scores.append(model.predict(smi))
    print(f"RScorePred scores calculated for {len(scores)} molecules")
    return scores

def get_scscore(smiles):
    model = SCScorer()
    model.restore('/Users/tomer/private/chem/motif_based/other_metrics/scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz')
    scores = []
    for smi in smiles:
        score = model.get_score_from_smi(smi)
        scores.append(score[1])
    print(f"SC scores calculated for {len(scores)} molecules")
    return scores

def get_fscore(smiles):
    model = LitRankNet.load_from_checkpoint('/Users/tomer/private/chem/motif_based/other_metrics/fsscore/models/pretrain_graph_GGLGGL_ep242_best_valloss.ckpt', weights_only=False)
    scorer = Scorer(model=model)
    scores = scorer.score(smiles)
    print(f"FS scores calculated for {len(scores)} molecules")
    return scores

def add_scores2csv(csv_path, output_csv_path=None):
    df = pd.read_csv(csv_path)
    inchis = df['inchi'].tolist()
    smiles = []
    for inchi in inchis:
        mol = Chem.MolFromInchi(inchi, treatWarningAsError=False)
        if mol is None:
            smiles.append(None)
            continue
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True)
        smiles.append(smi)
    df['smiles'] = smiles
    df['PASscore'] = get_passcore(inchis)
    df['SAscore'] = get_sascore(smiles)
    df['RAscore'] = get_rascore(smiles)
    df['RScorePred'] = get_rscore(smiles)
    df['SCScore'] = get_scscore(smiles)
    df['FSScore'] = get_fscore(smiles)
    if output_csv_path is None:
        output_csv_path = csv_path
    df.to_csv(output_csv_path, index=False)

def plot_scores_correlation(csv_path, out_path='/Users/tomer/private/chem/motif_based/output_plots'):
    # For each metric plot correlation with the ExpertScore - plot scatter and best fit line, and print correlation coefficient (r^2)
    df = pd.read_csv(csv_path)
    import matplotlib.pyplot as plt
    import seaborn as sns
    metrics = ['PASscore', 'SAscore', 'RAscore', 'RScorePred', 'SCScore', 'FSScore']
    for metric in metrics:
        correlation = df[metric].corr(df['ExpertScore'])
        print(f"Correlation between {metric} and ExpertScore: {correlation:.4f}")
        print(df[metric].to_list())
        plt.figure(figsize=(6, 6))
        sns.regplot(x=metric, y='ExpertScore', data=df, scatter_kws={'s': 10}, line_kws={'color': 'red'}, ci=False)
        plt.title(f"{metric} vs ExpertScore (r^2={correlation**2:.3f})")
        plt.xlabel(metric)
        plt.ylabel('ExpertScore')
        plt.tight_layout()
        plt.savefig(os.path.join(out_path, f"{metric}_vs_ExpertScore.png"))

if __name__ == '__main__':
    # add_scores2csv('/Users/tomer/private/chem/motif_based/sampled_molecules/molecule_scores.csv')
    add_scores2csv('/Users/tomer/private/chem/motif_based/sampled_molecules2/expert_test.csv')
    # plot_scores_correlation('/Users/tomer/private/chem/motif_based/sampled_molecules2/expert_test.csv')
