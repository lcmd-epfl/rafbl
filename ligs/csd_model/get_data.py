import sys
sys.path.append("../..")

import pandas as pd
from modsel.featurizer import bl_pool_featurize

csd_feats = pd.read_csv("../csd_pool.csv")

x_full, feats = bl_pool_featurize(csd_feats,symm=True, full=True)
x_full.to_csv(f"csd_full.csv", index=False)