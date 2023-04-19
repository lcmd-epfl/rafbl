import numpy as np

from sklearn.base import TransformerMixin
from sklearn.linear_model import LinearRegression, BayesianRidge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_absolute_error, mean_squared_error

class BLFeatureSelector(TransformerMixin):
    def __init__(
        self,
        full_features,
        candidate_features,
        model_selection=LeaveOneOut(),
        save_experiment=None,
    ):
        self.feats = full_features

        self.candidate_features = candidate_features
        self.cfs = np.asarray(list(candidate_features))

        self.model_selection = model_selection

        self.save_experiment = save_experiment

    def fit(self, X, y=None):

        mlr = LinearRegression()
        brr = BayesianRidge(compute_score=True)
        sts = StandardScaler()
        scores = []

        y = np.asarray(y)
        X = np.asarray(X)
        ss_tot = np.power(y - np.mean(y), 2).mean()

        for comb in self.cfs:
            ind = np.argwhere(np.isin(self.feats, comb)).reshape(-1)
            x = sts.fit_transform(X[:, ind])
            score = mlr.fit(x, y).score(x, y)
            scores += (score,)

        self.mask_ = None
        self.ols_scores_ = np.asarray(scores)
        self.ols_tops_ = self.ols_scores_ > self.ols_scores_.max() * 0.8


        self.maes_ = []
        self.rmses_ = []
        self.r2s_ = []

        self.sd_maes_ = []
        self.sd_rmses_ = []

        self.masks_ = []

        self.coefs_ = []
        self.intercepts_ = []
        self.lmls_ = []
        self.f_r2s_ = []
        self.f_maes_ = []
        self.f_rmses_ = []

        for comb in self.cfs[self.ols_tops_]:
            ind = np.argwhere(np.isin(self.feats, comb)).reshape(-1)
            self.masks_ += (ind,)
            x = X[:, ind]

            errors = []

            for i_train, i_test in self.model_selection.split(x):

                x_train = sts.fit_transform(x[i_train]) 
                x_test = sts.transform(x[i_test])

                brr.fit(x_train, y[i_train])
                yp = brr.predict(x_test)

                errors += yp - y[i_test],
    
            x = sts.fit_transform(x)
            
            brr.fit(x, y)
            y_p = brr.predict(x)
            self.coefs_ += brr.coef_,
            self.intercepts_ += brr.intercept_,
            self.lmls_ += brr.scores_[-1],
            self.f_r2s_ += brr.score(x,y),
            self.f_maes_ += mean_absolute_error(y,y_p),
            self.f_rmses_ += np.sqrt(mean_squared_error(y,y_p)),

            ae = np.abs(errors).mean(axis=1)
            se = np.power(errors, 2).mean(axis=1)

            ss_res = se.mean()
            self.maes_ += (ae.mean(),)
            self.rmses_ += (np.sqrt(ss_res),)
            self.r2s_ += (1 - ss_res / ss_tot,)

            self.sd_maes_ += (ae.std(),)
            self.sd_rmses_ += (se.std(),)

        self.mask_ = np.argwhere(
            np.isin(self.feats, self.cfs[self.ols_tops_][np.argmax(self.r2s_)])
        ).reshape(-1)

        """
        self.mask_ = np.argwhere(
            np.isin(self.feats, self.cfs[self.ols_tops_][np.argmin(self.rmses_)])
        ).reshape(-1)
        """

        if self.save_experiment is not None:
            from pandas import DataFrame

            DataFrame(
                data={
                    "combs": [str(list(x)) for x in self.cfs[self.ols_tops_]],
                    "mask": [str(list(x)) for x in self.masks_],
                    "coef": [str(list(np.around(x,6))) for x in self.coefs_],
                    "inters": self.intercepts_,
                    "f_r2": self.f_r2s_,
                    "f_mae": self.f_maes_,
                    "f_rmse": self.f_rmses_,
                    "lml": self.lmls_,
                    "mae": self.maes_,
                    "sd_mae": self.sd_maes_,
                    "rmse": self.rmses_,
                    "sd_rmse": self.sd_rmses_,
                    "r2": self.r2s_,
                }
            ).to_csv(self.save_experiment + ".csv", index=False)

        return self

    def transform(self, X):
        # check_is_fitted(self)
        X = np.asarray(X)
        if self.mask_ is not None:
            X = X[:, self.mask_]

        return X

    def adjr2(r2, n, k):
        return 1 - (((1 - r2) * (n - 1)) / (n - k - 1))