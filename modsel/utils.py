from itertools import combinations, chain, product

from scipy.stats import norm

def get_combinations(features, n=None):

    if type(features) == list:
        if type(n) != int:
            raise ValueError(
                f"A signle list of features requires 'n' to be an integer. 'n' is currently '{n}'."
            )
        combs = combinations(features, n)
    elif type(features) == tuple:
        combs = chain( *(product(*feats) for feats in features) )
        # combs = product(*features)
    else:
        raise NotImplementedError()
    
    return combs


def ei(mu,si,b_mu):
    u = (mu - b_mu) / si
    pdf = norm.pdf(u)
    cdf = norm.cdf(u)
    
    res = si * (u * cdf + pdf)
    res[si == 0] = 0

    return res

def adjr2(r2, n, k):
    return 1 - (((1 - r2) * (n - 1)) / (n - k - 1))