import sys
import numpy as np
from sklearn.metrics import r2_score
from sklearn.model_selection import LeavePOut
from sklearn.linear_model import BayesianRidge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

try:
    from ForwardStepCandidates_2 import ForwardStep_py
except ModuleNotFoundError:
    print("'ForwardStepCandidates_2.py' not installed. Download it and 'loo_q2.py' from\nhttps://github.com/SigmanGroup/enzyme-MLR-GluER/tree/main/aMD-and-IFD-models (last accessed 13.12.2023).")
    exit(1)

sys.path.append("../")


def adjr2(y, y_pred, n, k):
    if n - 2 < k:
        return np.inf
    else:
        return 1 - (((1 - r2_score(y, y_pred)) * (n - 1)) / (n - k - 1))


def fs_selector(x_full, target):
    res1, _, _, _, _ = ForwardStep_py(
        x_full, target, n_steps=3, n_candidates=50, reg=BayesianRidge()
    )

    cv = LeavePOut(1)
    estimator = BayesianRidge()
    scoring = adjr2
    scaler = StandardScaler()

    reg = make_pipeline(scaler, estimator)

    res1_sub = res1.copy()

    adj_scores = []  # adjusted score of CV left-out
    fin_scores = []  # final scores, all training data

    for mod in res1_sub["Model"]:
        x_cv = x_full[list(mod)].to_numpy()
        y_cv = x_full[target].to_numpy()

        preds = np.array([])
        tests = np.array([])

        for i_train, i_test in cv.split(x_cv, y=y_cv):
            x_train = x_cv[i_train]
            x_test = x_cv[i_test]

            y_train = y_cv[i_train]
            y_test = y_cv[i_test]

            reg.fit(x_train, y_train)
            preds = np.concatenate((preds, reg.predict(x_test)))

            tests = np.concatenate((tests, y_test))

        adj_scores.append(scoring(tests, preds, x_cv.shape[0], x_cv.shape[1]))

        reg.fit(x_cv, y_cv)
        fin_sc = scoring(y_cv, reg.predict(x_cv), x_cv.shape[0], x_cv.shape[1])
        fin_scores.append(fin_sc)

    res1_sub[f"adj_score"] = adj_scores
    res1_sub[f"fin_score"] = fin_scores

    return res1_sub
