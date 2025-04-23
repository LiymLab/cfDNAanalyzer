import numpy as np
import pandas as pd
from collections import Counter
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from methods.machine_learning.classifiers import get_classifier
from run_feature_selection import run_feature_selection, infer_fs_type
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def process_fold_voting(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args, method):
    n_samples = len(test_idx)
    modality_probs = []   
    modality_preds = []   
    modality_aucs = []  

    y_train = y.iloc[train_idx]
    y_test = y.iloc[test_idx]

    for i, X in enumerate(Xs):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        fs_type = infer_fs_type(args)
        idx_list = run_feature_selection(X_train, y_train, args, fs_type,
                                         verbose=False, fs_tag=f"{clf_name}_fold{fold_idx}_mod{i}")
        X_train_fs = X_train.iloc[:, idx_list]
        X_test_fs = X_test.iloc[:, idx_list]

        clf = get_classifier(clf_name)
        clf.fit(X_train_fs, y_train)

        if hasattr(clf, "predict_proba"):
            prob = clf.predict_proba(X_test_fs)  # (n_samples, classNum)
            prob_train = clf.predict_proba(X_train_fs)
        else:
            pred_labels = clf.predict(X_test_fs)  # (n_samples,)
            prob = np.zeros((len(pred_labels), args.classNum))
            for j, label in enumerate(pred_labels):
                prob[j, int(label)] = 1.0
            pred_train = clf.predict(X_train_fs)
            prob_train = np.eye(args.classNum)[pred_train]

        modality_probs.append(prob)
        modality_preds.append(np.argmax(prob, axis=1))
        try:
            if args.classNum == 2:
                auc = roc_auc_score(y_train, prob_train[:, 1])
            else:
                auc = 1.0
        except Exception:
            auc = 1.0
        modality_aucs.append(auc)

    n_modalities = len(Xs)
    fused_rows = []
    for pos, sample_idx in enumerate(test_idx):
        sample_id = sample_ids[sample_idx]
        true_label = y.iloc[sample_idx]
        if method == "average":
            probs = np.array([modality_probs[m][pos] for m in range(n_modalities)])
            fused_prob = np.mean(probs, axis=0)
        elif method == "weighted":
            probs = np.array([modality_probs[m][pos] for m in range(n_modalities)])
            weights = np.array(modality_aucs)
            weights = weights / np.sum(weights)
            fused_prob = np.average(probs, axis=0, weights=weights)
        elif method == "majority":
            preds = [modality_preds[m][pos] for m in range(n_modalities)]
            majority_label = Counter(preds).most_common(1)[0][0]
            fused_prob = np.zeros(args.classNum)
            fused_prob[int(majority_label)] = 1.0
        else:
            raise ValueError("Unsupported voting_method")

        row = {'SampleID': sample_id, 'TrueLabel': true_label, 'Classifier': clf_name}
        for cls in range(args.classNum):
            row[f'Prob_Class{cls}'] = fused_prob[cls]
        fused_rows.append(row)
    return fused_rows

def process_fold_stacking(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args):
    y_train = y.iloc[train_idx]
    y_test = y.iloc[test_idx]

    base_train_preds = []  
    base_test_preds = []   

    for i, X in enumerate(Xs):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        fs_type = infer_fs_type(args)
        idx_list = run_feature_selection(X_train, y_train, args, fs_type,
                                         verbose=False, fs_tag=f"{clf_name}_fold{fold_idx}_mod{i}")
        X_train_fs = X_train.iloc[:, idx_list]
        X_test_fs = X_test.iloc[:, idx_list]

        clf = get_classifier(clf_name)
        clf.fit(X_train_fs, y_train)

        base_train_preds.append(clf.predict_proba(X_train_fs))  # shape: (n_train, classNum)
        base_test_preds.append(clf.predict_proba(X_test_fs))    # shape: (n_test, classNum)

    X_meta_train = np.hstack(base_train_preds)
    X_meta_test = np.hstack(base_test_preds)

    meta_clf = LogisticRegression()
    meta_clf.fit(X_meta_train, y_train)
    y_meta_probs = meta_clf.predict_proba(X_meta_test)  # shape: (n_test, classNum)

    rows = []
    for pos, sample_idx in enumerate(test_idx):
        sample_id = sample_ids[sample_idx]
        true_label = y.iloc[sample_idx]
        fused_prob = y_meta_probs[pos]
        row = {'SampleID': sample_id, 'TrueLabel': true_label, 'Classifier': clf_name}
        for cls in range(args.classNum):
            row[f'Prob_Class{cls}'] = fused_prob[cls]
        rows.append(row)
    return rows

def run_voting_fusion(Xs, clf_name, cv, y, sample_ids, args, fs_label, method, use_thread=False):
    fold_results = []
    if use_thread:
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_fold_voting, fold_idx, train_idx, test_idx,
                                Xs, y, sample_ids, clf_name, args, method)
                for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y))
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing folds"):
                fold_results.extend(future.result())
    else:
        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y)):
            result = process_fold_voting(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args, method)
            fold_results.extend(result)
    return pd.DataFrame(fold_results)

def run_stacking_fusion(Xs, clf_name, cv, y, sample_ids, args, fs_label, use_thread=False):
    fold_results = []
    if use_thread:
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_fold_stacking, fold_idx, train_idx, test_idx,
                                Xs, y, sample_ids, clf_name, args)
                for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y))
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing folds"):
                fold_results.extend(future.result())
    else:
        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(Xs[0], y)):
            result = process_fold_stacking(fold_idx, train_idx, test_idx, Xs, y, sample_ids, clf_name, args)
            fold_results.extend(result)
    return pd.DataFrame(fold_results)
