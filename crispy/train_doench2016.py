"""
Train a Doench 2016 Rule Set 2 model from scratch.

Reads the training data (FC_plus_RES_withPredictions.csv), extracts 627
features from each 30mer, trains an AdaBoostRegressor (100 trees, lr=0.1),
validates against CRISPOR, and exports the model to JSON.

Output: doench2016_model.json (same format as before, no pickle dependency)

Usage:
    python train_doench2016.py [--validate] [--export-onnx]
"""

import argparse
import json
import itertools
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
import warnings
warnings.filterwarnings('ignore')

# ── Feature extraction (same as doench2016.py) ──────────────────────────────

ALPHABET = ['A', 'T', 'C', 'G']

def get_alphabet(order):
    return ["".join(i) for i in itertools.product(ALPHABET, repeat=order)]

def nucleotide_features(seq, order):
    alphabet = get_alphabet(order)
    n_a = len(alphabet)
    n_pos = len(seq) - order + 1
    pd = np.zeros(n_a * n_pos)
    pi = np.zeros(n_a)
    for pos in range(n_pos):
        kmer = seq[pos:pos + order]
        if kmer in alphabet:
            idx = alphabet.index(kmer)
            pd[idx + pos * n_a] = 1.0
            pi[idx] += 1.0
    return pd, pi

def featurize_30mer(s):
    """Extract 627 features from a 30mer string."""
    pd1, pi1 = nucleotide_features(s, 1)
    pd2, pi2 = nucleotide_features(s, 2)

    # GC features on guide (positions 5-24)
    guide = s[5:25]
    gc = sum(1 for c in guide if c in 'GC')
    gc_feat = np.array([1 if gc > 10 else 0, 1 if gc < 10 else 0, gc])

    # NGGX interaction
    alphabet2 = get_alphabet(2)
    nx = s[24] + s[27]
    nggx = np.zeros(len(alphabet2))
    if nx in alphabet2:
        nggx[alphabet2.index(nx)] = 1.0

    # Melting temperature
    try:
        from Bio.SeqUtils import MeltingTemp as Tm
        tm = np.array([
            Tm.Tm_NN(s),
            Tm.Tm_NN(s[19:24]),
            Tm.Tm_NN(s[11:19]),
            Tm.Tm_NN(s[6:11]),
        ])
    except ImportError:
        tm = np.zeros(4)

    return np.concatenate([pd1, pi1, pd2, pi2, gc_feat, nggx, tm])

# ── Model export ─────────────────────────────────────────────────────────────

def export_model_to_json(model, path):
    """Export a fitted GradientBoostingRegressor to JSON."""
    data = {
        'init': float(model.init_.predict(np.zeros((1, 627)))[0]),
        'learning_rate': float(model.learning_rate),
        'n_estimators': int(model.n_estimators),
        'max_depth': int(model.max_depth),
        'estimators': [],
    }
    for i in range(model.n_estimators):
        tree = model.estimators_[i, 0].tree_
        data['estimators'].append({
            'n_nodes': int(tree.node_count),
            'children_left': tree.children_left.astype(int).tolist(),
            'children_right': tree.children_right.astype(int).tolist(),
            'feature': tree.feature.astype(int).tolist(),
            'threshold': tree.threshold.astype(float).tolist(),
            'value': tree.value.flatten().astype(float).tolist(),
        })
    with open(path, 'w') as f:
        json.dump(data, f)
    print("Exported model to %s (%d trees, %d estimators)" % (path, len(data['estimators']), model.n_estimators))

# ── CRISPOR validation ──────────────────────────────────────────────────────

CRISPOR_GUIDES = {
    "TATGAACTTACTGTGATAGC": 51,
    "ATGAACTTACTGTGATAGCA": 70,
    "AGCGGGGAGCAGTTTACATC": 51,
    "AAGCGGGGAGCAGTTTACAT": 57,
}
INPUT_SEQ = "GTTTCATAGGGCCATGCATTTTTACATCTAGGTGTATGAACTTACTGTGATAGCAGGGTGTAGGATATTATACCACAGACAGGAAGCGGGGAGCAGTTTACATCGGGTCGTCTGCGGAAGCTTGGCGGCTCTTTGACGGCGAGATCTCCGTAG"

def validate_against_crispor(model):
    print("\nValidation against CRISPOR website:")
    print("-" * 60)
    total_error = 0
    for guide, crispor_score in CRISPOR_GUIDES.items():
        pos = INPUT_SEQ.find(guide)
        if pos == -1:
            print("  %s: NOT FOUND" % guide)
            continue
        s30 = INPUT_SEQ[pos-5:pos+25]
        X = np.array([featurize_30mer(s30)])
        pred = np.clip(model.predict(X)[0], 0.0, 1.0)
        our_score = round(pred * 100)
        diff = our_score - crispor_score
        total_error += abs(diff)
        print("  %s  Ours: %2d  CRISPOR: %2d  %s" % (
            guide, our_score, crispor_score,
            "OK" if diff == 0 else "err=%+d" % diff))
    print("-" * 60)
    print("Mean absolute error: %.1f" % (total_error / len(CRISPOR_GUIDES)))

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Train Doench 2016 model")
    parser.add_argument("--validate", action="store_true", help="Run CRISPOR validation")
    parser.add_argument("--export-onnx", action="store_true", help="Also export ONNX format")
    parser.add_argument("--output", default="doench2016_model.json", help="Output model path")
    parser.add_argument("--data", default="../../src/crisporWebsite/bin/Azimuth-2.0/azimuth/data/FC_plus_RES_withPredictions.csv")
    args = parser.parse_args()

    # Load training data
    print("Loading training data from %s ..." % args.data)
    df = pd.read_csv(args.data)
    print("  %d guides, %d genes" % (len(df), df['Target gene'].nunique()))

    # Extract features
    print("Extracting features ...")
    X = np.array([featurize_30mer(s) for s in df['30mer']])
    y = df['score_drug_gene_rank'].values
    print("  Feature matrix: %s" % str(X.shape))

    # Check for NaN/inf
    nan_count = np.isnan(X).sum()
    inf_count = np.isinf(X).sum()
    if nan_count > 0:
        print("  WARNING: %d NaN values in features, replacing with 0" % nan_count)
        X = np.nan_to_num(X)
    if inf_count > 0:
        print("  WARNING: %d inf values in features, replacing with 0" % inf_count)
        X = np.nan_to_num(X)

    # Train model (same hyperparameters as original Doench 2016)
    print("Training GradientBoostingRegressor (100 trees, lr=0.1, max_depth=3) ...")
    model = GradientBoostingRegressor(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=3,
        subsample=1.0,
        min_samples_split=2,
        min_samples_leaf=1,
        random_state=42,
    )
    model.fit(X, y)

    # Cross-validation
    print("Running 5-fold cross-validation ...")
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2')
    print("  R2 scores: %s" % ["%.3f" % s for s in cv_scores])
    print("  Mean R2: %.3f +/- %.3f" % (cv_scores.mean(), cv_scores.std()))

    spearman_cv = cross_val_score(model, X, y, cv=5, scoring='r2')
    print("  CV R2 mean: %.3f" % spearman_cv.mean())

    # Export
    export_model_to_json(model, args.output)

    # Also save as joblib for reference
    import pickle
    joblib_path = args.output.replace('.json', '.pkl')
    with open(joblib_path, 'wb') as f:
        pickle.dump(model, f)
    print("Exported sklearn model to %s" % joblib_path)

    # Optional ONNX export
    if args.export_onnx:
        try:
            from skl2onnx import convert_sklearn
            from skl2onnx.common.data_types import FloatTensorType
            initial_type = [('input', FloatTensorType([None, 627]))]
            onnx_model = convert_sklearn(model, initial_types=initial_type)
            onnx_path = args.output.replace('.json', '.onnx')
            with open(onnx_path, 'wb') as f:
                f.write(onnx_model.SerializeToString())
            print("Exported ONNX model to %s" % onnx_path)
        except ImportError:
            print("WARNING: skl2onnx not installed, skipping ONNX export")

    # Validate
    if args.validate:
        validate_against_crispor(model)

    print("\nDone. Model saved to %s" % args.output)

if __name__ == '__main__':
    main()
"""
Train a Doench 2016 Rule Set 2 model from scratch.

Reads the training data (FC_plus_RES_withPredictions.csv), extracts 629
features from each 30mer (627 base + 2 position), trains a
GradientBoostingRegressor (100 trees, lr=0.1), validates against
CRISPOR, and exports the model to JSON.

Output: doench2016_model_v2.json (no pickle dependency, no sklearn needed at runtime)

Usage:
    python train_doench2016.py [--validate] [--output doench2016_model_v2.json]
"""

import argparse
import json
import itertools
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
import warnings
warnings.filterwarnings('ignore')

# ── Feature extraction (same as doench2016.py) ──────────────────────────────

ALPHABET = ['A', 'T', 'C', 'G']

def get_alphabet(order):
    return ["".join(i) for i in itertools.product(ALPHABET, repeat=order)]

def nucleotide_features(seq, order):
    alphabet = get_alphabet(order)
    n_a = len(alphabet)
    n_pos = len(seq) - order + 1
    pd = np.zeros(n_a * n_pos)
    pi = np.zeros(n_a)
    for pos in range(n_pos):
        kmer = seq[pos:pos + order]
        if kmer in alphabet:
            idx = alphabet.index(kmer)
            pd[idx + pos * n_a] = 1.0
            pi[idx] += 1.0
    return pd, pi

def featurize_30mer(s, aa_cut=None, pct_peptide=None):
    """Extract features from a 30mer string.

    Base features (627): nucleotide, GC, NGGX, Tm
    + Position features (2): Amino Acid Cut position, Percent Peptide
    Total: 629
    """
    pd1, pi1 = nucleotide_features(s, 1)
    pd2, pi2 = nucleotide_features(s, 2)

    # GC features on guide (positions 5-24)
    guide = s[5:25]
    gc = sum(1 for c in guide if c in 'GC')
    gc_feat = np.array([1 if gc > 10 else 0, 1 if gc < 10 else 0, gc])

    # NGGX interaction
    alphabet2 = get_alphabet(2)
    nx = s[24] + s[27]
    nggx = np.zeros(len(alphabet2))
    if nx in alphabet2:
        nggx[alphabet2.index(nx)] = 1.0

    # Melting temperature
    try:
        from Bio.SeqUtils import MeltingTemp as Tm
        tm = np.array([
            Tm.Tm_NN(s),
            Tm.Tm_NN(s[19:24]),
            Tm.Tm_NN(s[11:19]),
            Tm.Tm_NN(s[6:11]),
        ])
    except ImportError:
        tm = np.zeros(4)

    # Position features: Amino Acid Cut position and Percent Peptide
    pos_feats = np.array([
        float(aa_cut) if aa_cut is not None and not np.isnan(float(aa_cut)) else 0.0,
        float(pct_peptide) if pct_peptide is not None and not np.isnan(float(pct_peptide)) else 101.0,
    ])

    return np.concatenate([pd1, pi1, pd2, pi2, gc_feat, nggx, tm, pos_feats])

# ── Model export ─────────────────────────────────────────────────────────────

def export_model_to_json(model, path):
    """Export a fitted GradientBoostingRegressor to JSON."""
    data = {
        'init': float(model.init_.predict(np.zeros((1, 629)))[0]),
        'learning_rate': float(model.learning_rate),
        'n_estimators': int(model.n_estimators),
        'max_depth': int(model.max_depth),
        'estimators': [],
    }
    for i in range(model.n_estimators):
        tree = model.estimators_[i, 0].tree_
        data['estimators'].append({
            'n_nodes': int(tree.node_count),
            'children_left': tree.children_left.astype(int).tolist(),
            'children_right': tree.children_right.astype(int).tolist(),
            'feature': tree.feature.astype(int).tolist(),
            'threshold': tree.threshold.astype(float).tolist(),
            'value': tree.value.flatten().astype(float).tolist(),
        })
    with open(path, 'w') as f:
        json.dump(data, f)
    print("Exported model to %s (%d trees)" % (path, len(data['estimators'])))

# ── CRISPOR validation ──────────────────────────────────────────────────────

CRISPOR_GUIDES = {
    "TATGAACTTACTGTGATAGC": 51,
    "ATGAACTTACTGTGATAGCA": 70,
    "AGCGGGGAGCAGTTTACATC": 51,
    "AAGCGGGGAGCAGTTTACAT": 57,
}
INPUT_SEQ = "GTTTCATAGGGCCATGCATTTTTACATCTAGGTGTATGAACTTACTGTGATAGCAGGGTGTAGGATATTATACCACAGACAGGAAGCGGGGAGCAGTTTACATCGGGTCGTCTGCGGAAGCTTGGCGGCTCTTTGACGGCGAGATCTCCGTAG"

def validate_against_crispor(model):
    print("\nValidation against CRISPOR website:")
    print("-" * 60)
    print("(Position features set to defaults: aa_cut=NaN, pct_peptide=NaN)")
    total_error = 0
    for guide, crispor_score in CRISPOR_GUIDES.items():
        pos = INPUT_SEQ.find(guide)
        if pos == -1:
            print("  %s: NOT FOUND" % guide)
            continue
        s30 = INPUT_SEQ[pos-5:pos+25]
        X = np.array([featurize_30mer(s30, aa_cut=None, pct_peptide=None)])
        pred = np.clip(model.predict(X)[0], 0.0, 1.0)
        our_score = round(pred * 100)
        diff = our_score - crispor_score
        total_error += abs(diff)
        print("  %s  Ours: %2d  CRISPOR: %2d  %s" % (
            guide, our_score, crispor_score,
            "OK" if diff == 0 else "err=%+d" % diff))
    print("-" * 60)
    print("Mean absolute error: %.1f" % (total_error / len(CRISPOR_GUIDES)))

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Train Doench 2016 model")
    parser.add_argument("--validate", action="store_true", help="Run CRISPOR validation")
    parser.add_argument("--output", default="doench2016_model_v2.json", help="Output model path")
    parser.add_argument("--data", default="../../src/crisporWebsite/bin/Azimuth-2.0/azimuth/data/FC_plus_RES_withPredictions.csv")
    args = parser.parse_args()

    # Load training data
    print("Loading training data from %s ..." % args.data)
    df = pd.read_csv(args.data)
    print("  %d guides, %d genes" % (len(df), df['Target gene'].nunique()))

    # Extract features
    print("Extracting features ...")
    X = np.array([
        featurize_30mer(row['30mer'], row['Amino Acid Cut position'], row['Percent Peptide'])
        for _, row in df.iterrows()
    ])
    y = df['score_drug_gene_rank'].values
    print("  Feature matrix: %s  (627 base + 2 position = 629)" % str(X.shape))

    # Check for NaN/inf
    nan_count = np.isnan(X).sum()
    if nan_count > 0:
        print("  WARNING: %d NaN values in features, replacing with 0" % nan_count)
        X = np.nan_to_num(X)

    # Train model (same hyperparameters as original Doench 2016)
    print("Training GradientBoostingRegressor (100 trees, lr=0.1, max_depth=3) ...")
    model = GradientBoostingRegressor(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=3,
        subsample=1.0,
        min_samples_split=2,
        min_samples_leaf=1,
        random_state=42,
    )
    model.fit(X, y)

    # Cross-validation
    print("Running 5-fold cross-validation ...")
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2')
    print("  R2 scores: %s" % ["%.3f" % s for s in cv_scores])
    print("  Mean R2: %.3f +/- %.3f" % (cv_scores.mean(), cv_scores.std()))

    # Export
    export_model_to_json(model, args.output)

    # Validate
    if args.validate:
        validate_against_crispor(model)

    print("\nDone. Model saved to %s" % args.output)
    print("NOTE: Do NOT update doench2016_model.json or pipeline.py until")
    print("this model has been validated end-to-end in the pipeline.")

if __name__ == '__main__':
    main()
