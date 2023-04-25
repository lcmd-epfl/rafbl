import sys

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import BayesianRidge

from modsel.utils import ei
from modsel.featurizer import bl_featurize, bl_pool_featurize

# indices corresponding to: 'nbo1', 'f2ka', 'onhh1'
feat = [0, 127, 203]

lig_feats = pd.read_csv("ligs/lit_pool.csv")
lig_resps = pd.read_csv("ligs/responses.csv")

X, Y = bl_featurize(lig_feats, lig_resps, uid="c_smiles", task="oa", target="ddg")

used_in_tr = X.index.to_list()

X = X.reset_index(drop=True)
Y = Y.reset_index(drop=True)

x = X.to_numpy()
y = Y.to_numpy()

x = x[:, feat]
y = y

brr = BayesianRidge()
sts = StandardScaler()

x = sts.fit_transform(x)

print(brr.fit(x, y).score(x, y))
print(brr.coef_, brr.intercept_)
print(np.asarray(X.columns)[feat])

if sys.argv[1] == "1":
    pool = pd.read_csv("ligs/csd_pool.csv")
elif sys.argv[1] == "2":
    pool = lig_feats

plfs = bl_pool_featurize(pool, symm=False)
flfs = bl_pool_featurize(lig_feats)

x_p = plfs.to_numpy()[:, feat]
x_p = sts.transform(x_p)

y_p, y_s = brr.predict(x_p, return_std=True)

if sys.argv[1] == "1":
    print(
        np.array([pool["name"][x] for x in np.argsort(ei(y_p, y_s, y_p.max()))[::-1]])
    )

elif sys.argv[1] == "2":
    # Gives file number of ligand in lit_pool.csv
    print(
        np.array(
            [
                pool["name"][x]
                for x in np.argsort(ei(y_p, y_s, y_p.max()))[::-1]
                if x not in used_in_tr
            ]
        )
    )
