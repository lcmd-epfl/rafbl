import sys

import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import BayesianRidge
from sklearn.model_selection import LeaveOneOut

from modsel.selector import BLFeatureSelector
from modsel.utils import get_combinations
from modsel.featurizer import (
    bl_featurize,
    e_feats,
    sa_feats,
    sp_feats,
    ta_feats,
    tp_feats,
)

lig_feats = pd.read_csv("ligs/lit_pool.csv")
lig_resps = pd.read_csv("ligs/responses.csv")

i = int(sys.argv[1])
tasks = ["oa", "cp", "cc", "da_f"]

X, Y = bl_featurize(
    lig_feats, lig_resps, uid="c_smiles", task=tasks[i], target="ddg", full=False
)

x = X.reset_index(drop=True)
y = Y.reset_index(drop=True)

candidate_features = get_combinations(
    ((e_feats, sa_feats, tp_feats), (e_feats, sp_feats, ta_feats))
)


pipe = Pipeline(
    [
        (
            "selector",
            BLFeatureSelector(
                x.columns,
                candidate_features,
                model_selection=LeaveOneOut(),
                save_experiment=f"fm_{i+1}",
            ),
        ),
        ("scaler", StandardScaler()),
        ("estimator", BayesianRidge()),
    ]
)

sel = pipe.named_steps["selector"]

print(pipe.fit(x, y).score(x, y))
print(x.columns[sel.mask_])
print(sel.mask_)
