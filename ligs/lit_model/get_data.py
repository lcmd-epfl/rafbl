import sys

sys.path.append("../..")

import pandas as pd
from modsel.featurizer import bl_featurize, bl_pool_featurize

lig_feats = pd.read_csv("../lit_pool.csv")
lig_resps = pd.read_csv("../responses.csv")

tasks = ["oa", "cp", "cc", "da_f"]

for task in tasks:
    x_full, feats = bl_featurize(
        lig_feats, lig_resps, uid="c_smiles", task=task, target="ddg", full=True
    )

    x_full["name"] = x_full["name"].map(lambda x: str(x)[:-3] + "xyz")
    x_full.to_csv(f"{task}_full.csv", index=False)

x_full, feats = bl_pool_featurize(lig_feats, symm=True, full=True)
x_full.to_csv(f"lit_full.csv", index=False)
