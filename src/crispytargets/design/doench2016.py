"""
Doench 2016 Rule Set 2 on-target scoring.

Pure numpy implementation using the pre-trained GradientBoostingRegressor
exported as JSON. No sklearn, no onnxruntime — works on any Python 3.x
with just numpy.

Requirements: numpy

Usage:
    from doench2016 import score_guides
    scores = score_guides(["ATCGATCGATCGATCGATCG", ...])
"""

import os
import json
import numpy as np
import itertools

_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(_DIR, "data", "doench2016_model.json")
RAW_ALPHABET = ['A', 'T', 'C', 'G']

# ── feature extraction ───────────────────────────────────────────────────────

def _get_alphabet(order):
    return ["".join(i) for i in itertools.product(RAW_ALPHABET, repeat=order)]


def _nucleotide_features(seq, order):
    """Position-dependent and position-independent nucleotide features."""
    alphabet = _get_alphabet(order)
    n = len(alphabet)
    n_pos = len(seq) - order + 1
    pd = np.zeros(n * n_pos)
    pi = np.zeros(n)
    for pos in range(n_pos):
        kmer = seq[pos:pos + order]
        if kmer in alphabet:
            idx = alphabet.index(kmer)
            pd[idx + pos * n] = 1.0
            pi[idx] += 1.0
    return pd, pi


def _featurize_guide(guide_20bp, flank_5prime=None, flank_3prime=None, pam_seq=None):
    """Extract 627 features from a 20bp guide sequence.

    Constructs a 30mer: 5bp 5' flank + 20bp guide + PAM + 2bp 3' flank.
    Uses actual PAM and flanking sequences when provided.
    The model was trained on 30mers with the guide at positions 5-24.
    """
    guide = guide_20bp.upper().replace('U', 'T')
    pam = pam_seq.upper() if pam_seq else 'NGG'
    left = flank_5prime[-5:] if flank_5prime and len(flank_5prime) >= 5 else 'N' * 5
    right = flank_3prime[:2] if flank_3prime and len(flank_3prime) >= 2 else 'N' * 2
    s = left + guide + pam + right
    assert len(s) == 30

    pd1, pi1 = _nucleotide_features(s, 1)
    pd2, pi2 = _nucleotide_features(s, 2)

    guide = s[5:25]
    gc = sum(1 for c in guide if c in 'GC')
    gc_feat = np.array([1 if gc > 10 else 0, 1 if gc < 10 else 0, gc])

    alphabet = _get_alphabet(2)
    nx = s[24] + s[27]
    nggx = np.zeros(len(alphabet))
    if nx in alphabet:
        nggx[alphabet.index(nx)] = 1.0

    # Melting temperature features (must be computed, not zeros)
    try:
        import Bio.SeqUtils.MeltingTemp as Tm
        tm_segments = [(19, 24), (11, 19), (6, 11)]
        tm_feats = np.zeros(4)
        tm_feats[0] = Tm.Tm_NN(s)
        tm_feats[1] = Tm.Tm_NN(s[tm_segments[0][0]:tm_segments[0][1]])
        tm_feats[2] = Tm.Tm_NN(s[tm_segments[1][0]:tm_segments[1][1]])
        tm_feats[3] = Tm.Tm_NN(s[tm_segments[2][0]:tm_segments[2][1]])
    except ImportError:
        tm_feats = np.zeros(4)

    return np.concatenate([pd1, pi1, pd2, pi2, gc_feat, nggx, tm_feats])


# ── model loading and prediction ─────────────────────────────────────────────

_model = None


def _load_model():
    global _model
    if _model is None:
        with open(MODEL_PATH) as f:
            _model = json.load(f)
    return _model


def _predict_tree(tree, X):
    """Traverse a single decision tree for a batch of samples (vectorized)."""
    n = X.shape[0]
    node = np.zeros(n, dtype=int)
    children_left = np.array(tree['children_left'])
    children_right = np.array(tree['children_right'])
    features = np.array(tree['feature'])
    thresholds = np.array(tree['threshold'])
    values = np.array(tree['value'])

    for _ in range(tree['n_nodes']):
        is_leaf = children_left[node] == -1
        go_left = (~is_leaf) & (X[np.arange(n), features[node]] <= thresholds[node])
        go_right = (~is_leaf) & (~go_left)
        node[go_left] = children_left[node[go_left]]
        node[go_right] = children_right[node[go_right]]

    return values[node]


# ── public API ───────────────────────────────────────────────────────────────

def score_guides(guides, flanks=None, pams=None):
    """Score a list of 20bp guide sequences using Doench 2016 Rule Set 2.

    Args:
        guides: list of 20bp strings (no PAM)
        flanks: optional list of (left_flank, right_flank) tuples
        pams: optional list of 3bp PAM strings

    Returns:
        numpy array of scores (0-1, higher = better)
    """
    model = _load_model()
    if flanks:
        X = np.array([_featurize_guide(g, l, r, p)
                       for g, (l, r), p in zip(guides, flanks, pams)])
    elif pams:
        X = np.array([_featurize_guide(g, pam_seq=p) for g, p in zip(guides, pams)])
    else:
        X = np.array([_featurize_guide(g) for g in guides])

    pred = np.full(X.shape[0], model['init'])
    lr = model['learning_rate']
    for est in model['estimators']:
        pred += lr * _predict_tree(est, X)

    return np.clip(pred, 0.0, 1.0)


def score_guide(guide):
    """Score a single 20bp guide sequence."""
    return score_guides([guide])[0]
