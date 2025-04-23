import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import linear_kernel, polynomial_kernel, rbf_kernel, sigmoid_kernel
from sklearn.preprocessing import StandardScaler
from methods.machine_learning.classifiers import get_classifier
from run_feature_selection import run_feature_selection, infer_fs_type
from snf import compute
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def transform_features(X_train_fs, X_test_fs, method):
    if method == "pca":
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train_fs)
        X_test_scaled = scaler.transform(X_test_fs)
        # pca = PCA(n_components=min(10, X_train_scaled.shape[1]))
        pca = PCA(n_components=min(10, X_train_scaled.shape[0], X_train_scaled.shape[1]))
        return pca.fit_transform(X_train_scaled), pca.transform(X_test_scaled)
    elif method == "linear":
        return linear_kernel(X_train_fs), linear_kernel(X_test_fs, X_train_fs)
    elif method == "polynomial":
        return polynomial_kernel(X_train_fs), polynomial_kernel(X_test_fs, X_train_fs)
    elif method == "rbf":
        return rbf_kernel(X_train_fs), rbf_kernel(X_test_fs, X_train_fs)
    elif method == "sigmoid":
        return sigmoid_kernel(X_train_fs), sigmoid_kernel(X_test_fs, X_train_fs)
    elif method == "snf":
        train_affinity = compute.make_affinity(X_train_fs)
        full_affinity = compute.make_affinity(np.vstack([X_train_fs, X_test_fs]))
        test_affinity = full_affinity[-len(X_test_fs):, :-len(X_test_fs)]
        return train_affinity, test_affinity
    else:
        raise ValueError("Unsupported trans method")

def process_fold_trans(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, method, args):
    y_train = y.iloc[train_idx]
    y_test = y.iloc[test_idx]

    trans_train_list, trans_test_list = [], []

    for i, X in enumerate(Xs):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        fs_type = infer_fs_type(args)
        idx_list = run_feature_selection(X_train, y_train, args, fs_type,
                                         verbose=False, fs_tag=f"{clf_name}_fold{fold_idx}_mod{i}")
        X_train_fs = X_train.iloc[:, idx_list].values
        X_test_fs = X_test.iloc[:, idx_list].values

        X_train_trans, X_test_trans = transform_features(X_train_fs, X_test_fs, method)
        trans_train_list.append(X_train_trans)
        trans_test_list.append(X_test_trans)

    X_train_all = np.concatenate(trans_train_list, axis=1)
    X_test_all = np.concatenate(trans_test_list, axis=1)

    clf = get_classifier(clf_name)
    clf.fit(X_train_all, y_train)

    if hasattr(clf, "predict_proba"):
        y_probs = clf.predict_proba(X_test_all)
    else:
        pred_labels = clf.predict(X_test_all)
        y_probs = []
        for label in pred_labels:
            prob = [0.0] * args.classNum
            prob[int(label)] = 1.0
            y_probs.append(prob)

    rows = []
    for pos, idx in enumerate(test_idx):
        sample_test_id = sample_ids[idx]
        true_label = y.iloc[idx]
        row = {'SampleID': sample_test_id, 'TrueLabel': true_label, 'Classifier': clf_name}
        for cls in range(args.classNum):
            row[f'Prob_Class{cls}'] = y_probs[pos][cls]
        rows.append(row)
    return rows

def run_trans_fusion(Xs, clf_name, cv, y, sample_ids, args, fs_label, method, use_thread=False):
    fold_results = []
    if use_thread:
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_fold_trans, fold_idx, train_idx, test_idx,
                                Xs, y, sample_ids, clf_name, method, args)
                for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y))
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing folds"):
                fold_results.extend(future.result())
    else:
        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y)):
            result = process_fold_trans(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, method, args)
            fold_results.extend(result)
    return pd.DataFrame(fold_results)
