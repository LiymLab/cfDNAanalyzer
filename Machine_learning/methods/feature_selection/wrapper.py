import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
from mlxtend.feature_selection import ExhaustiveFeatureSelector as EFS
from boruta import BorutaPy


def boruta_feature_selection(X, y, percentage=0.2):
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    selector = BorutaPy(rf, n_estimators='auto', random_state=42, verbose=0)
    selector.fit(X.values, y.values.ravel())

    rankings = selector.ranking_
    # top_n = max(1, int(len(rankings) * percentage))
    if (percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if (percentage >= 1):
            top_n = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            top_n = max(1, int(len(rankings) * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")
              
    top_features_idx = np.argsort(rankings)[:top_n]
    return list(X.columns[top_features_idx])


def forward_feature_selection(X, y, percentage=0.2):
    # n_features = max(1, int(X.shape[1] * percentage))
    if (percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if (percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")
      
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    n_splits = min(5, max(2, min(np.bincount(y))))

    sfs = SFS(rf, k_features=n_features, forward=True, floating=False,
              scoring='accuracy', cv=StratifiedKFold(n_splits=n_splits), n_jobs=-1)
    sfs.fit(X.values, y.values)
    return list(X.columns[list(sfs.k_feature_idx_)])


def backward_feature_selection(X, y, percentage=0.2):
    # n_features = max(1, int(X.shape[1] * percentage))
    if (percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if (percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")
      
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    n_splits = min(5, max(2, min(np.bincount(y))))

    sfs = SFS(rf, k_features=n_features, forward=False, floating=False,
              scoring='accuracy', cv=StratifiedKFold(n_splits=n_splits), n_jobs=-1)
    sfs.fit(X.values, y.values)
    return list(X.columns[list(sfs.k_feature_idx_)])


def exhaustive_feature_selection(X, y, percentage=0.2):
    # n_features = max(1, int(X.shape[1] * percentage))
    if (percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if (percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")
      
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    n_splits = min(5, max(2, min(np.bincount(y))))

    efs = EFS(rf, min_features=n_features, max_features=n_features,
              scoring='accuracy', cv=StratifiedKFold(n_splits=n_splits), n_jobs=-1)
    efs.fit(X.values, y.values)
    return list(X.columns[list(efs.best_idx_)])


def recursive_feature_elimination(X, y, percentage=0.2):
    # n_features = max(1, int(X.shape[1] * percentage))
    if (percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if (percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")
      
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    n_splits = min(5, max(2, min(np.bincount(y))))

    rfecv = RFECV(estimator=rf, step=1, cv=StratifiedKFold(n_splits=n_splits),
                  scoring='accuracy', n_jobs=-1)
    rfecv.fit(X.values, y.values)

    rankings = rfecv.ranking_
    top_n = max(1, int(len(rankings) * percentage))
    top_features_idx = np.argsort(rankings)[:top_n]
    return list(X.columns[top_features_idx])
