import sys
import pandas as pd
from modsel.fs_selector import fs_selector


task = int(sys.argv[1])
paths = [
    "ligs/lit_model/oa_full.csv",
    "ligs/lit_model/cp_full.csv",
    "ligs/lit_model/cc_full.csv",
    "ligs/lit_model/da_f_full.csv",
]

data = pd.read_csv(paths[task])
x_full = data.drop(columns=["set", "smiles", "c_smiles", "ee", "temp", "order", "name"])
target = "ddg"

models = fs_selector(x_full, target)

s_models = models.sort_values(by=["adj_score"], ascending=False)
s_models.to_csv(f"results/results_{task}.csv", index=False)
