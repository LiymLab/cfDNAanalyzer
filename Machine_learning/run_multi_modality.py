"""
This script performs multi-modality fusion analysis for classification.

1. Parse arguments (input_dir, fusion type, FS method, classifier list, CV type, etc.)
2. Load all modality files from input_dir (each .csv is a modality)
3. Build feature selection combinations (filter, wrapper, embedded, hybrid)
4. For each fusion method (concat / model / trans):
    a. For each FS combination:
        b. For each classifier:
            i. Record runtime and memory usage
            ii. For each fold:
                - Apply feature selection on each modality using training data
                - Use same selected features on test data
                - Fuse feature representations or predictions
                - Train classifier and predict probabilities
            iii. Collect prediction results
            iv. Evaluate metrics (AUC, F1, etc.)
            v. Save confusion matrix and ROC (for binary classification)
5. Output results and metrics to fusion-specific folders under output_dir

Fusion methods:
- concat: feature concatenation
- model: average / weighted / majority voting, stacking
- trans: PCA, kernels (linear, rbf, poly, sigmoid), SNF

Usage example:
1. generate npy cache to increase operation speed
python generate_npy_cache.py --input_dir /home/zjp/projects/202402_cfDNAIntegratedTool/250402_downstream_result/two-class/01preprocessing/0317_format_var_std/CNA --cache_dir /home/zky/projects/202312_cfDNA_integrated_tools/2504DA-npy/testdata

2. run the main script
nohup python run_multi_modality.py \
  --modality multi \
  --input_dir /home/zky/projects/202312_cfDNA_integrated_tools/2504DA-figure/multi-modality/multi-class/input \
  --fusion_type concat model trans \
  --model_method average weighted stack \
  --trans_method pca rbf snf \
  --classifierMulti KNN SVM \
  --cvMulti LOO \
  --nsplitMulti 5 \
  --filterMethod IG CHI FS FCBF CC LVF MAD DR MI RLF SURF MSURF \
  --filterFrac 0.023618328 \
  --wrapperMethod BOR \
  --wrapperFrac 0.023618328 \
  --embeddedMethod LASSO RIDGE ELASTICNET RF \
  --embeddedFrac 0.023618328 \
  --hybridType FE \
  --hybridMethod1 IG \
  --hybridFrac1 0.023618328 \
  --hybridMethod2 LASSO \
  --hybridFrac2 0.023618328\
  --DA_output_dir /home/zky/projects/202312_cfDNA_integrated_tools/2504DA-figure/multi-modality/multi-classs/results > multi_modality-multiclass.out 2>&1 &
"""
#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import time
import tracemalloc
import psutil
from tqdm import tqdm
from config_parser import get_args
from methods.machine_learning.train_test_split import get_cv
from methods.machine_learning.generate_metrics import evaluate_classification
from methods.machine_learning.visualization import plot_confusion_matrix, plot_roc_curve
from run_feature_selection import build_fs_combinations_from_args
from methods.machine_learning.model_based_fusion import run_voting_fusion, run_stacking_fusion
from methods.machine_learning.concatenation_fusion import run_concatenation_fusion
from methods.machine_learning.transformation_based_fusion import run_trans_fusion

import os
os.environ["OMP_NUM_THREADS"] = "1"  
os.environ["OPENBLAS_NUM_THREADS"] = "1" 

def load_modalities_smart(input_dir):
    """Load multiple modalities from CSV or pre-saved .npz files (preferred)."""
    npz_files = sorted(glob.glob(os.path.join(input_dir, '*.npz')))
    csv_files = sorted(glob.glob(os.path.join(input_dir, '*.csv')))

    modality_names = []
    modality_data = []

    if npz_files:
        print(f"Loading .npz files from {input_dir}")
        for npz_file in npz_files:
            name = os.path.splitext(os.path.basename(npz_file))[0]
            with np.load(npz_file, allow_pickle=True) as data:
                X = pd.DataFrame(data["data"])
                y = pd.Series(data["label"])
                sample_ids = data['sample_ids']
            modality_names.append(name)
            modality_data.append({"X": X, "y": y, "sample_ids": sample_ids})
    else:
        print(f"Loading .csv files from {input_dir}")
        for csv_file in csv_files:
            name = os.path.splitext(os.path.basename(csv_file))[0]
            df = pd.read_csv(csv_file)
            df.columns = df.columns.str.lower()

            y = df["label"]
            if "sample" in df.columns:
                sample_ids = df["sample"].tolist()
            else:
                sample_ids = df.index.tolist()

            X = df.drop(columns=["label", "sample"], errors="ignore")

            modality_names.append(name)
            modality_data.append({"X": X, "y": y, "sample_ids": sample_ids})

    return modality_names, modality_data


def merge_classifier_outputs(method_dir, classifier_list):
    all_preds, all_metrics = [], []
    for clf_name in classifier_list:
        pred_path = os.path.join(method_dir, f"{clf_name}_predictions.csv")
        metrics_path = os.path.join(method_dir, f"{clf_name}_metrics.csv")
        if os.path.exists(pred_path):
            all_preds.append(pd.read_csv(pred_path))
        if os.path.exists(metrics_path):
            all_metrics.append(pd.read_csv(metrics_path))
    if all_preds:
        pd.concat(all_preds).to_csv(os.path.join(method_dir, "all_classifiers_predictions.csv"), index=False)
    if all_metrics:
        pd.concat(all_metrics).to_csv(os.path.join(method_dir, "all_classifiers_metrics.csv"), index=False)


def merge_all_methods_outputs(base_dir, method_list, classifier_list):
    all_preds, all_metrics = [], []
    for method in method_list:
        method_dir = os.path.join(base_dir, method)
        pred_file = os.path.join(method_dir, "all_classifiers_predictions.csv")
        metric_file = os.path.join(method_dir, "all_classifiers_metrics.csv")
        if os.path.exists(pred_file):
            all_preds.append(pd.read_csv(pred_file))
        if os.path.exists(metric_file):
            all_metrics.append(pd.read_csv(metric_file))
    # if all_preds:
    #     all_preds = all_preds.drop_duplicates()
    #     pd.concat(all_preds).to_csv(os.path.join(base_dir, "all_methods_all_classifiers_predictions.csv"), index=False)
    # if all_metrics:
    #     all_metrics = all_metrics.drop_duplicates()
    #     pd.concat(all_metrics).to_csv(os.path.join(base_dir, "all_methods_all_classifiers_metrics.csv"), index=False)
    if all_preds:
        combined_preds = pd.concat(all_preds, ignore_index=True)
        combined_preds = combined_preds.drop_duplicates()
        combined_preds.to_csv(os.path.join(base_dir, "all_methods_all_classifiers_predictions.csv"), index=False)
    if all_metrics:
        combined_metrics = pd.concat(all_metrics, ignore_index=True)
        combined_metrics = combined_metrics.drop_duplicates()
        combined_metrics.to_csv(os.path.join(base_dir, "all_methods_all_classifiers_metrics.csv"), index=False)


def save_outputs(
    prediction_df, clf_name, fs_label, method_dir, args,
    fusion_type=None, fusion_method=None,
    elapsed_time=None, peak_memory=None
):
    process = psutil.Process()
    mem_used = round(process.memory_info().rss / 1024 / 1024, 4)
    total_time = elapsed_time if elapsed_time is not None else None
    peak_mem   = peak_memory if peak_memory is not None else None

    prediction_df["Classifier"] = clf_name
    prediction_df["FS_Combination"] = fs_label
    prediction_df["FusionType"] = fusion_type
    prediction_df["FusionMethod"] = fusion_method

    prob_cols = [f'Prob_Class{i}' for i in range(args.classNum)]
    for col in prob_cols:
        if col in prediction_df.columns:
            prediction_df[col] = prediction_df[col].round(6)
    # pred_order = ['SampleID', 'TrueLabel', 'Classifier', 'FS_Combination', 'FusionType', 'FusionMethod'] + prob_cols
    pred_order = ['SampleID', 'TrueLabel', 'FusionType', 'FusionMethod', 'FS_Combination', 'Classifier'] + prob_cols
  
    available_cols = [c for c in pred_order if c in prediction_df.columns]
    other_cols = [c for c in prediction_df.columns if c not in available_cols]
    prediction_df = prediction_df[available_cols + other_cols]

    metrics = evaluate_classification(prediction_df, class_num=args.classNum)
    metrics["Classifier"] = clf_name
    metrics["FS_Combination"] = fs_label
    metrics["FusionType"] = fusion_type
    metrics["FusionMethod"] = fusion_method
    metrics["TotalTime_sec"]  = total_time
    metrics["FinalMemory_MB"] = mem_used
    metrics["PeakMemory_MB"]  = peak_mem

    metrics_df = pd.DataFrame([metrics])
    # metric_order = ['Classifier', 'FS_Combination', 'FusionType', 'FusionMethod']
    metric_order = ['FusionType', 'FusionMethod', 'FS_Combination', 'Classifier']
   
    other_metrics = [c for c in metrics_df.columns if c not in metric_order]
    metrics_df = metrics_df[metric_order + other_metrics]
    
    prediction_df.to_csv(os.path.join(method_dir, f"{clf_name}_predictions.csv"), index=False)
    metrics_df.to_csv(os.path.join(method_dir, f"{clf_name}_metrics.csv"), index=False)

    # plot_confusion_matrix(
    #     prediction_df, class_num=args.classNum,
    #     title=f"{clf_name}_{fs_label}",
    #     save_path=os.path.join(method_dir, f"{clf_name}_confusion_matrix.pdf")
    # )
    # if args.classNum == 2:
    #     plot_roc_curve(
    #         prediction_df, title=f"{clf_name}_{fs_label}",
    #         save_path=os.path.join(method_dir, f"{clf_name}_roc.pdf")
    #     )

    print(f"Saved results to: {method_dir}")
    print(f"Time = {total_time}s | Memory = {mem_used}MB | Peak = {peak_mem}MB")


def run_multi_modality_analysis(modality_names, modality_data, args):
    fs_combinations = build_fs_combinations_from_args(args)
    if not fs_combinations:
        raise ValueError("No valid FS parameters provided.")

    for fusion_type in args.fusion_type:
        # output_base_dir = os.path.join(args.DA_output_dir, "multi_modality", fusion_type)
        output_base_dir = os.path.join(args.DA_output_dir, fusion_type)
        os.makedirs(output_base_dir, exist_ok=True)
        print(f"\nFusing modalities using: {fusion_type}")

        for fs_combo in fs_combinations:
            for key in fs_combo:
                setattr(args, key, fs_combo[key])
            fs_label = fs_combo.get("fs_label", "default")
            method_names = []

            for clf_name in args.classifierMulti:
                Xs, ys, sample_ids = [], None, None
                for mdata in modality_data:
                    Xs.append(mdata["X"])
                    if ys is None:
                        ys = mdata["y"]
                        sample_ids = mdata["sample_ids"]
                cv = get_cv(args.cvMulti, args.nsplitMulti, ys)

                # -------- Concat --------
                if fusion_type == "concat":
                    method_name = "concat"
                    method_dir = os.path.join(output_base_dir, fs_label, method_name)
                    os.makedirs(method_dir, exist_ok=True)
                    method_names.append(method_name)

                    print(f"\n>>> [Concat] Classifier: {clf_name}, FS: {fs_label}")
                    start_time = time.time()
                    tracemalloc.start()

                    prediction_df = run_concatenation_fusion(
                        Xs, clf_name, cv, ys, sample_ids,
                        args, fs_label, use_thread=True
                    )
                    
                    end_time = time.time()
                    current, peak = tracemalloc.get_traced_memory()
                    tracemalloc.stop()

                    elapsed_time = round(end_time - start_time, 4)
                    peak_mem    = round(peak / 1024 / 1024, 4)

                    save_outputs(
                        prediction_df, clf_name, fs_label, method_dir, args,
                        fusion_type="concat", fusion_method="concat",
                        elapsed_time=elapsed_time, peak_memory=peak_mem
                    )

                # -------- Model-based --------
                elif fusion_type == "model":
                    for model_method in args.model_method:
                        method_dir = os.path.join(output_base_dir, fs_label, model_method)
                        os.makedirs(method_dir, exist_ok=True)
                        if model_method not in method_names:
                            method_names.append(model_method)

                        print(f"\n>>> [Model-{model_method}] Classifier: {clf_name}, FS: {fs_label}")
                        start_time = time.time()
                        tracemalloc.start()

                        if model_method in ["average", "weighted", "majority"]:
                            prediction_df = run_voting_fusion(
                                Xs, clf_name, cv, ys, sample_ids,
                                args, fs_label, method=model_method, use_thread=True
                            )
                        elif model_method == "stack":
                            prediction_df = run_stacking_fusion(
                                Xs, clf_name, cv, ys, sample_ids,
                                args, fs_label, use_thread=True
                            )
                        else:
                            raise ValueError("Unsupported model fusion method.")

                        end_time = time.time()
                        current, peak = tracemalloc.get_traced_memory()
                        tracemalloc.stop()

                        elapsed_time = round(end_time - start_time, 4)
                        peak_mem    = round(peak / 1024 / 1024, 4)

                        save_outputs(
                            prediction_df, clf_name, fs_label, method_dir, args,
                            fusion_type="model", fusion_method=model_method,
                            elapsed_time=elapsed_time, peak_memory=peak_mem
                        )

                # -------- Transformation-based --------
                elif fusion_type == "trans":
                    for trans_method in args.trans_method:
                        method_dir = os.path.join(output_base_dir, fs_label, trans_method)
                        os.makedirs(method_dir, exist_ok=True)
                        if trans_method not in method_names:
                            method_names.append(trans_method)

                        print(f"\n>>> [Trans-{trans_method}] Classifier: {clf_name}, FS: {fs_label}")
                        start_time = time.time()
                        tracemalloc.start()

                        prediction_df = run_trans_fusion(
                            Xs, clf_name, cv, ys, sample_ids,
                            args, fs_label, method=trans_method, use_thread=True
                        )

                        end_time = time.time()
                        current, peak = tracemalloc.get_traced_memory()
                        tracemalloc.stop()

                        elapsed_time = round(end_time - start_time, 4)
                        peak_mem    = round(peak / 1024 / 1024, 4)

                        save_outputs(
                            prediction_df, clf_name, fs_label, method_dir, args,
                            fusion_type="trans", fusion_method=trans_method,
                            elapsed_time=elapsed_time, peak_memory=peak_mem
                        )

                else:
                    raise ValueError(f"Unsupported fusion type: {fusion_type}")

            for method in method_names:
                method_dir = os.path.join(output_base_dir, fs_label, method)
                merge_classifier_outputs(method_dir, args.classifierMulti)

            merge_all_methods_outputs(
                os.path.join(output_base_dir, fs_label), method_names, args.classifierMulti
            )


def cleanup_intermediate_csvs(base_dir):
    for f in glob.glob(os.path.join(base_dir, "**/*_predictions.csv"), recursive=True):
        if not f.endswith("all_classifiers_predictions.csv"):
            os.remove(f)
    for f in glob.glob(os.path.join(base_dir, "**/*_metrics.csv"), recursive=True):
        if not f.endswith("all_classifiers_metrics.csv"):
            os.remove(f)


def main():
    args = get_args()
    os.makedirs(args.DA_output_dir, exist_ok=True)
    modality_names, modality_data = load_modalities_smart(args.input_dir)
    print(f"Loaded modalities: {modality_names}")
    run_multi_modality_analysis(modality_names, modality_data, args)
    # cleanup_intermediate_csvs(args.DA_output_dir)


if __name__ == '__main__':
    main()
