import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_classif, chi2, VarianceThreshold
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from skrebate import ReliefF, SURF, MultiSURF, TuRF
from sklearn_relief import Relief
from scipy.sparse import lil_matrix, diags

### ========== 通用 Filter 方法 ========== ###

def calculate_information_gain(X, y):
    mi = mutual_info_classif(X, y)
    return sorted(zip(X.columns, mi), key=lambda x: x[1], reverse=True)

def sorted_chi_square_scores(X, y):
    X_scaled = MinMaxScaler().fit_transform(X)
    scores, _ = chi2(X_scaled, y)
    return sorted(zip(X.columns, scores), key=lambda x: x[1], reverse=True)

def calculate_mutual_information(X, y):
    return calculate_information_gain(X, y)

def calculate_permutation_importance(X, y, model=None, n_repeats=10):
    if model is None:
        model = RandomForestClassifier(random_state=42)
    model.fit(X, y)
    result = permutation_importance(model, X, y, n_repeats=n_repeats)
    return sorted(zip(X.columns, result.importances_mean), key=lambda x: x[1], reverse=True)

def low_variance_filter(X, y=None, threshold=0.01):
    selector = VarianceThreshold(threshold=threshold)
    selector.fit(X)
    return sorted(zip(X.columns, selector.variances_), key=lambda x: x[1], reverse=True)

def mean_absolute_difference(X, y=None):
    mad = X.mad()
    return sorted(zip(X.columns, mad), key=lambda x: x[1], reverse=True)

def calculate_missing_ratios(X, y=None):
    missing_ratios = X.isnull().mean()
    return sorted(zip(X.columns, missing_ratios), key=lambda x: x[1])

def calculate_highest_correlation_per_feature(X, y=None):
    corr = X.corr().abs()
    highest = corr.apply(lambda x: sorted(x[x < 1.0])[-1] if len(x[x < 1.0]) else 0, axis=1)
    return sorted(zip(X.columns, highest), key=lambda x: x[1], reverse=True)

def dispersion_ratio(X, y):
    overall_var = np.var(X, axis=0)
    within_var = np.zeros(X.shape[1])
    for cls in np.unique(y):
        samples = X[y == cls]
        cls_var = np.var(samples, axis=0)
        proportion = len(samples) / len(X)
        within_var += proportion * cls_var
    overall_var[overall_var == 0] = np.finfo(float).eps
    ratios = within_var / overall_var
    return sorted(zip(X.columns, ratios), key=lambda x: x[1], reverse=True)

### ========== FCBF & Entropy ========== ###

def entropy(x):
    vals, counts = np.unique(x, return_counts=True)
    p = counts / counts.sum()
    return -np.sum(p * np.log2(p + 1e-12))

def symmetricalUncertain(x, y):
    n = float(len(y))
    vals_y = np.unique(y)
    Hx = entropy(x)
    Hy = entropy(y)
    Hxy = sum(entropy(x[y == val]) * (y == val).sum() / n for val in vals_y)
    IG = Hx - Hxy
    return 2 * IG / (Hx + Hy + 1e-12)

class FCBF:
    def __init__(self, threshold=0.01):
        self.threshold = threshold
        self.selected_features = []
        self.importances = []

    def fit(self, X, y):
        X_vals = X.values
        y_vals = y.values
        su = np.array([symmetricalUncertain(X_vals[:, i], y_vals) for i in range(X.shape[1])])
        self.importances = su
        idx_sorted = np.argsort(su)[::-1]
        selected = idx_sorted[su[idx_sorted] >= self.threshold]

        while len(selected) > 0:
            current = selected[0]
            self.selected_features.append(current)
            selected = selected[1:]
            selected = [i for i in selected if symmetricalUncertain(X_vals[:, current], X_vals[:, i]) < su[current]]

    def fit_transform(self, X, y):
        self.fit(X, y)
        feature_names = X.columns
        return sorted(zip(feature_names, self.importances), key=lambda x: x[1], reverse=True)

### ========== Fisher Score ========== ###

def construct_W(X, y):
    n = X.shape[0]
    labels = np.unique(y)
    W = lil_matrix((n, n), dtype=np.float64)
    for lbl in labels:
        idx = np.where(y == lbl)[0]
        for i in idx:
            for j in idx:
                W[i, j] = 1.0 / len(idx) if i != j else 1
    return W.tocsc()

def fisher_score(X, y):
    X_vals = X.values
    W = construct_W(X_vals, y)
    D = np.array(W.sum(axis=0)).flatten()
    D_mat = diags(D)

    Xt = X_vals.T
    wm = (Xt @ D_mat @ np.ones(X.shape[0]) / D.sum()).reshape(-1, 1)
    fr_hat = Xt - wm

    num = np.sum((fr_hat * D) * fr_hat, axis=1)
    L = D_mat - W
    den = np.sum((fr_hat * L.diagonal()) * fr_hat, axis=1)

    scores = num / np.maximum(den, 1e-12)
    return sorted(zip(X.columns, scores), key=lambda x: x[1], reverse=True)

### ========== Relief 系列 ========== ###

def calculate_relief_feature_importance(X, y, method="ReliefF"):
    fs_map = {
        'ReliefF': ReliefF(),
        'SURF': SURF(),
        'MultiSURF': MultiSURF()
    }
    fs = fs_map[method]
    fs.fit(X.values, y.values)
    return sorted(zip(X.columns, fs.feature_importances_), key=lambda x: x[1], reverse=True)

def calculate_turf_relieff_feature_importance(X, y):
    from skrebate import TuRF
    X_numeric = X.select_dtypes(include=[np.number])
    X_numeric = X_numeric.reindex(sorted(X_numeric.columns), axis=1)

    X_np = X_numeric.to_numpy()
    y_np = y.to_numpy()
    headers = list(X_numeric.columns)

    print(f"[TuRF] X shape: {X_np.shape}")
    print(f"[TuRF] Headers sample: {headers[:5]}")

    fs = TuRF(core_algorithm='ReliefF')
    fs.fit(X_np, y_np, headers=headers)

    return sorted(zip(headers, fs.feature_importances_), key=lambda x: x[1], reverse=True)