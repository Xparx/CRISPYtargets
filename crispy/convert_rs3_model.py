"""
convert_rs3_model.py

One-time conversion: load RuleSet3.pkl and dump as native LightGBM JSON.

Run: conda activate rs3test
     python convert_rs3_model.py --rs3-dir ~/projects/src/rs3
"""

import os
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Convert RuleSet3.pkl to LightGBM JSON"
    )
    parser.add_argument(
        "--rs3-dir", required=True,
        help="Path to the rs3 repo (e.g., ~/projects/src/rs3)"
    )
    parser.add_argument(
        "--output", default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "RuleSet3.json"
        ),
        help="Output path for JSON model"
    )
    args = parser.parse_args()

    import joblib
    import lightgbm as lgb

    pkl_path = os.path.join(args.rs3_dir, "rs3", "RuleSet3.pkl")
    print("Loading: %s" % pkl_path)
    model = joblib.load(pkl_path)
    print("Type: %s" % type(model).__name__)

    # Get the LightGBM booster
    if hasattr(model, "booster_"):
        booster = model.booster_
    elif hasattr(model, "dump_model"):
        booster = model
    else:
        print("ERROR: Unknown model type %s" % type(model).__name__)
        sys.exit(1)

    booster.save_model(args.output)
    print("Saved: %s" % args.output)

    # Verify by reloading the saved model
    lgb_model = lgb.Booster(model_file=args.output)
    print("Verification: loaded JSON model, %d trees" % lgb_model.num_trees())


if __name__ == "__main__":
    main()
