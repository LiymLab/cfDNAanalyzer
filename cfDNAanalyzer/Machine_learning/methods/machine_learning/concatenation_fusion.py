import pandas as pd
import numpy as np
from sklearn.base import clone
from methods.machine_learning.classifiers import get_classifier
from run_feature_selection import run_feature_selection, infer_fs_type
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def process_fold_concat(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args):
    y_train = y.iloc[train_idx]
    y_test = y.iloc[test_idx]

    X_train_parts, X_test_parts = [], []
    for i, X in enumerate(Xs):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        fs_type = infer_fs_type(args)
        idx_list = run_feature_selection(X_train, y_train, args, fs_type,
                                         verbose=False, fs_tag=f"{clf_name}_fold{fold_idx}_mod{i}")
        X_train_fs = X_train.iloc[:, idx_list]
        X_test_fs = X_test.iloc[:, idx_list]
        X_train_parts.append(X_train_fs.values)
        X_test_parts.append(X_test_fs.values)

    X_train_concat = np.concatenate(X_train_parts, axis=1)
    X_test_concat = np.concatenate(X_test_parts, axis=1)

    clf = get_classifier(clf_name)
    clf.fit(X_train_concat, y_train)

    if hasattr(clf, "predict_proba"):
        y_probs = clf.predict_proba(X_test_concat)
    else:
        pred_labels = clf.predict(X_test_concat)
        y_probs = []
        for label in pred_labels:
            prob = [0.0] * args.classNum
            prob[int(label)] = 1.0
            y_probs.append(prob)

    rows = []
    for i, idx in enumerate(test_idx):
        sample_id = sample_ids[idx]
        true_label = y.iloc[idx]
        prob = y_probs[i]
        row = {'SampleID': sample_id, 'TrueLabel': true_label, 'Classifier': clf_name}
        for cls in range(args.classNum):
            row[f'Prob_Class{cls}'] = prob[cls]
        rows.append(row)

    return rows

def run_concatenation_fusion(Xs, clf_name, cv, y, sample_ids, args, fs_label, use_thread=False):
    fold_results = []
    if use_thread:
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_fold_concat, fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args)
                for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y))
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing folds"):
                #fold_results.append(future.result())
                fold_results.extend(future.result())
    else:
        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y)):
            result = process_fold_concat(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args)
            fold_results.append(result)

    return pd.DataFrame(fold_results)
