"""
Rule Set 3 on-target scoring.

Uses the GPP/rs3 model (Doench et al.) via native LightGBM JSON.
No sklearn dependency for prediction.

Requirements: numpy, pandas, lightgbm, sglearn

Usage:
    from rs3score import score_guides
    scores = score_guides(["AGGTGTATGAACTTACTGTGATAGCAGGGTG", ...])
"""

import os
import numpy as np
import pandas as pd

_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(_DIR, "data", "RuleSet3.json")

_model = None


def _load_model():
    global _model
    if _model is None:
        import lightgbm as lgb
        _model = lgb.Booster(model_file=MODEL_PATH)
    return _model


def score_guides(context_sequences, sequence_tracr="Hsu2013", n_jobs=1):
    """Score 30mer context sequences using Rule Set 3.

    Args:
        context_sequences: list of 30bp strings
            Format: 5bp 5' flank + 20bp guide + 3bp PAM + 2bp 3' flank
        sequence_tracr: tracrRNA version ("Hsu2013" or "Chen2013")
        n_jobs: parallel jobs for featurization

    Returns:
        numpy array of raw scores (not 0-1 bounded)
    """
    import sglearn

    for s in context_sequences:
        if len(s) != 30:
            raise ValueError(
                "All context sequences must be 30 nucleotides, "
                "got %d: %s" % (len(s), s)
            )

    model = _load_model()

    featurized = sglearn.featurize_guides(
        context_sequences, n_jobs=n_jobs
    )

    ref_tracrs = ["Hsu2013", "Chen2013"]
    for tracr in ref_tracrs:
        if isinstance(sequence_tracr, str):
            featurized[tracr + " tracr"] = int(sequence_tracr == tracr)
        else:
            featurized[tracr + " tracr"] = (
                pd.Series(sequence_tracr) == tracr
            ).astype(int).tolist()

    return model.predict(featurized)
