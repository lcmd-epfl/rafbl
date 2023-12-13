import sys
import pandas as pd
import numpy as np
from ast import literal_eval
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import BayesianRidge
from sklearn.pipeline import make_pipeline
from modsel.featurizer import bl_pool_featurize
from modsel.utils import ei

task = int(sys.argv[1])

### Fit best model
models = pd.read_csv("results/results_0.csv")
data = pd.read_csv("ligs/lit_model/oa_full.csv")
target = "ddg"
best_feats = list(literal_eval(models.iloc[0]["Model"]))

reg = make_pipeline(StandardScaler(), BayesianRidge())

x, y = data[best_feats], data[target]
reg.fit(x, y)

### Evaluate candidate pool
if task == 1:
    pool = pd.read_csv("ligs/csd_pool.csv")
elif task == 2:
    pool = pd.read_csv("ligs/lit_pool.csv")
else:
    print(
        f"Options: 1 -> csd ligands, 2 -> literature ligands. '{task}' is not an option."
    )
    exit(0)

plfs = bl_pool_featurize(pool)

y_p, y_u = reg.predict(plfs[best_feats], return_std=True)

print(np.array([pool["name"][x] for x in np.argsort(ei(y_p, y_u, y_p.max()))[::-1]]))

si = np.argsort(ei(y_p, y_u, y_p.max()))[::-1]
print(y_p[si], y_u[si], si)
