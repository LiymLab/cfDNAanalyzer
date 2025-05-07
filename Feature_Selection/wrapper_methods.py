import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
# from sklearn.feature_selection import SequentialFeatureSelector as SFS
from mlxtend.feature_selection import ExhaustiveFeatureSelector as EFS
from boruta import BorutaPy
import argparse
from load_data import load_data
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import factorial
import time
import psutil
import datetime
import json
from contextlib import contextmanager

def comb(n, r):
    if n < r:
        return 0
    return factorial(n) // (factorial(r) * factorial(n - r))

@contextmanager
def timer_and_memory_monitor():
    """Monitor execution time and memory usage of code block"""
    process = psutil.Process()
    start_time = time.time()
    start_memory = process.memory_info().rss / 1024 / 1024  # MB
    peak_memory = start_memory
    
    try:
        yield lambda: (process.memory_info().rss / 1024 / 1024, peak_memory)
        
    finally:
        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024
        execution_time = end_time - start_time
        peak_memory = max(peak_memory, end_memory)
        
        return {
            'execution_time': execution_time,
            'start_memory': start_memory,
            'end_memory': end_memory,
            'peak_memory': peak_memory
        }

def feature_importances(classifier, X, y):
    classifier.fit(X, y)
    importances = classifier.feature_importances_
    feature_importance = pd.Series(importances, index=X.columns).sort_values(ascending=False)
    return feature_importance

# def print_features_and_importances(features, importances, method):
    # print(f"************ {method} Feature Selection Results ************", flush=True)
    # print("Selected features:", ", ".join(features.tolist()), flush=True)
    # print(f"{method} Feature Importance:", flush=True)
    # for idx in features:
    #     print(f"{idx}: {importances.loc[idx]:.4f}", flush=True)

def boruta_feature_selection(X, y, classifier, percentage, callback=None):
    boruta_selector = BorutaPy(
        classifier,
        n_estimators='auto',
        random_state=42,
        verbose=0
    )
    boruta_selector.fit(X.values, y.values.ravel())

    # Extract feature rankings
    rankings = boruta_selector.ranking_

    # Calculate top percentage features count
    n_features = len(rankings)
    # top_percent_count = max(1, int(percentage * n_features))  # Ensure at least 1 feature is selected
    if ( percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if ( percentage >= 1):
            top_percent_count = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            top_percent_count = max(1, int(percentage * n_features))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")    

    # Get top ranked features
    top_features = np.argsort(rankings)[:top_percent_count]  # Return indices
    selected_features = X.columns[top_features]  # Select feature names by indices

    if callback:
        callback()
    return selected_features
  
def forward_feature_selection(X, y, classifier, percentage, callback=None):
    # n_features = max(1, int(X.shape[1] * percentage))
    if ( percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if ( percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")  

    n_splits = max(2, min(5, min(np.bincount(y))))
    skf = StratifiedKFold(n_splits=n_splits, random_state=42, shuffle=True)
    sfs = SFS(
        classifier, 
        k_features=n_features, 
        forward=True, 
        floating=False, 
        scoring='accuracy', 
        cv=skf,
        n_jobs=-1,
        verbose=0
    )
    sfs.fit(X, y)
    if callback:
        callback()
    return X.columns[list(sfs.k_feature_idx_)]

def backward_feature_selection(X, y, classifier, percentage, callback=None):
    # n_features = max(1, int(X.shape[1] * percentage))
    if ( percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if ( percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1") 
      
    n_splits = max(2, min(5, min(np.bincount(y))))
    skf = StratifiedKFold(n_splits=n_splits, random_state=42, shuffle=True)
    sfs = SFS(
        classifier,
        k_features=n_features,
        forward=False,
        floating=False,
        scoring='accuracy',
        cv=skf,
        n_jobs=-1,
        verbose=0
    )
    sfs.fit(X, y)
    if callback:
        callback()
    return X.columns[list(sfs.k_feature_idx_)]

def exhaustive_feature_selection(X, y, classifier, percentage, callback=None):
    # n_features = max(1, int(X.shape[1] * percentage))
    if ( percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if ( percentage >= 1):
            n_features = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            n_features = max(1, int(X.shape[1] * percentage))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1") 
      
    n_splits = max(2, min(5, min(np.bincount(y))))
    skf = StratifiedKFold(n_splits=n_splits, random_state=42, shuffle=True)
    efs = EFS(
        classifier, 
        min_features=n_features, 
        max_features=n_features, 
        scoring='accuracy', 
        cv=skf,
        n_jobs=-1
    )
    efs.fit(X, y)
    if callback:
        callback()
    return X.columns[list(efs.best_idx_)]

def recursive_feature_elimination(X, y, classifier, percentage, callback=None):
  
    n_splits = max(2, min(5, min(np.bincount(y))))
    skf = StratifiedKFold(n_splits=n_splits, random_state=42, shuffle=True)
    
    class CustomRFECV(RFECV):
        def _fit(self, X, y, step, **fit_params):
            result = super()._fit(X, y, step, **fit_params)
            if callback:
                callback()
            return result
    
    rfecv = CustomRFECV(
        estimator=classifier, 
        step=1, 
        cv=skf, 
        scoring='accuracy',
        n_jobs=-1
    )
    rfecv.fit(X, y)

    rankings = rfecv.ranking_

    # Calculate top percentage features count
    n_features = len(rankings)
    # top_percent_count = max(1, int(percentage * n_features))  # Ensure at least 1 feature is selected
    if ( percentage >= 1) or (isinstance(percentage, (int, float)) and 0 < percentage < 1):
        if ( percentage >= 1):
            top_percent_count = int(percentage)
        elif (isinstance(percentage, (int, float)) and 0 < percentage < 1) :
            top_percent_count = max(1, int(percentage * n_features))
    else:
        raise ValueError(" wrapperNum must be positive integer or float from 0 to 1")  
      
    # Get top ranked features
    top_features = np.argsort(rankings)[:top_percent_count]  # Return indices
    selected_features = X.columns[top_features]  # Select feature names by indices

    if callback:
        callback()
    return selected_features

def load_and_validate_data(file_path, transpose=False):
    """Load and validate data integrity"""
    data, X, y, feature_names = load_data(file_path, transpose=transpose)
    
    # Validate data
    if X is None or X.empty:
        raise ValueError("Feature matrix X is empty")
    if y is None or len(y) == 0:
        raise ValueError("Label vector y is empty")
    if X.shape[0] != len(y):
        raise ValueError(f"Sample count mismatch: X has {X.shape[0]} rows, y has {len(y)} labels")
        
    # Check for sufficient samples and features
    if X.shape[0] < 2:
        raise ValueError(f"Insufficient samples: only {X.shape[0]} samples")
    if X.shape[1] < 2:
        raise ValueError(f"Insufficient features: only {X.shape[1]} features")
        
    # Check data types
    if not np.issubdtype(X.dtypes.iloc[0], np.number):
        raise ValueError("Feature matrix contains non-numeric data")
        
    # Check for missing values
    if X.isnull().any().any():
        raise ValueError("Feature matrix contains missing values")
        
    # Check labels
    unique_labels = np.unique(y)
    if len(unique_labels) < 2:
        raise ValueError(f"Insufficient labels: only {len(unique_labels)} classes")
        
    return data, X, y, feature_names

def process_single_file(method_func, file_path, classifier, method_output_dir, method_name, percentage):
    filename = os.path.basename(file_path)
    print(f"\nProcessing file: {filename} using method: {method_name}", flush=True)
    
    performance_stats = {
        'filename': filename,
        'method': method_name,
        'start_time': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'status': 'failed',  # Default to failed, updated on success
        'error': None
    }
    
    with timer_and_memory_monitor() as monitor:
        try:
            # Load data
            data, X, y, feature_names = load_and_validate_data(file_path, transpose=False)
            print(f"Successfully loaded data: {X.shape[0]} samples, {X.shape[1]} features", flush=True)
            
            if X.shape[1] > 1000:
                print(f"Warning: High number of features ({X.shape[1]}), may require longer processing time", flush=True)
            
            position = 0
            leave = True
            
            # Feature selection process
            if method_name == 'BOR':
                with tqdm(total=100, desc=f"Boruta - {filename}", position=position, leave=leave) as pbar:
                    selected_features = method_func(X, y, classifier, percentage)
                    pbar.update(100)
            elif method_name in ['FFS', 'BFS']:
                total_features = X.shape[1]
                with tqdm(total=total_features, desc=f"{method_name} - {filename}", position=position, leave=leave) as pbar:
                    def update_progress(*args):
                        pbar.update(1)
                    selected_features = method_func(X, y, classifier, percentage, callback=update_progress)
            elif method_name == 'EFS':
                total_combinations = sum(comb(X.shape[1], r) for r in range(1, 6))
                with tqdm(total=total_combinations, desc=f"EFS - {filename}", position=position, leave=leave) as pbar:
                    def update_progress(*args):
                        pbar.update(1)
                    selected_features = method_func(X, y, classifier, percentage, callback=update_progress)
            else:  # RFE
                with tqdm(total=X.shape[1], desc=f"RFE - {filename}", position=position, leave=leave) as pbar:
                    def update_progress(*args):
                        pbar.update(1)
                    selected_features = method_func(X, y, classifier, percentage, callback=update_progress)

            # Validate selected features
            if len(selected_features) == 0:
                raise ValueError("No features were selected")

            selected_X = X[selected_features]
            importances = feature_importances(classifier, selected_X, y)
            
            final_data = pd.concat([
                data.index.to_series(name="Sample"),
                pd.DataFrame(y, columns=["label"]),
                selected_X
            ], axis=1)

            input_filename = os.path.splitext(filename)[0]
            output_file = os.path.join(method_output_dir, f"{input_filename}_{method_name}_selected.csv")
            final_data.to_csv(output_file, index=False)
            
            # Update performance stats
            performance_stats.update({
                'status': 'success',
                'n_samples': X.shape[0],
                'n_features': X.shape[1],
                'n_selected_features': len(selected_features),
                'selected_features': selected_features.tolist()
            })
            
            # print(f"\nFeature selection results for file {filename}:", flush=True)
            # print_features_and_importances(selected_features, importances, method_name)
            
        except Exception as e:
            performance_stats['error'] = str(e)
            raise e
                
    # Get monitoring results
    monitoring_results = monitor()
    performance_stats.update({
        'end_time': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'execution_time_seconds': monitoring_results['execution_time'],
        'start_memory_mb': monitoring_results['start_memory'],
        'end_memory_mb': monitoring_results['end_memory'],
        'peak_memory_mb': monitoring_results['peak_memory']
    })
        
    return performance_stats

def process_method(method, input_dir, output_dir, max_workers, percentage):
    methods_dict = {
        'BOR': boruta_feature_selection,
        'FFS': forward_feature_selection,
        'BFS': backward_feature_selection,
        'EFS': exhaustive_feature_selection,
        'RFS': recursive_feature_elimination
    }

    method_upper = method.upper()
    if method_upper not in methods_dict:
        print(f"Unsupported method {method}", flush=True)
        return

    method_func = methods_dict[method_upper]
    method_output_dir = output_dir
    
    os.makedirs(method_output_dir, exist_ok=True)
    
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
    
    method_stats = []
    # with tqdm(total=len(files), desc=f"Overall Progress - {method_upper}", position=0, leave=True) as main_pbar:
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for filename in files:
            file_path = os.path.join(input_dir, filename)
            futures.append(
                executor.submit(
                    process_single_file,
                    method_func,
                    file_path,
                    RandomForestClassifier(n_estimators=100, random_state=42),
                    method_output_dir,
                    method_upper,
                    percentage
                )
            )

        for future in as_completed(futures):
            stats = future.result()
            method_stats.append(stats)
            main_pbar.update(1)
    
    # Save method-level statistics
    method_summary = {
        'method': method_upper,
        'total_files': len(files),
        'successful_files': sum(1 for s in method_stats if s['status'] == 'success'),
        'failed_files': sum(1 for s in method_stats if s['status'] == 'failed'),
        'total_execution_time': sum(s.get('execution_time_seconds', 0) for s in method_stats),
        'average_execution_time': sum(s.get('execution_time_seconds', 0) for s in method_stats) / len(files),
        'max_peak_memory': max(s.get('peak_memory_mb', 0) for s in method_stats),
        'file_stats': method_stats
    }
    
    summary_file = os.path.join(output_dir, f"{method_upper}_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(method_summary, f, indent=4)

def main(input_dir, output_dir, methods, max_workers, percentage):
    os.makedirs(output_dir, exist_ok=True)
    
    # Use thread pool to process multiple methods concurrently
    with ThreadPoolExecutor(max_workers=len(methods)) as executor:
        futures = []
        for method in methods:
            futures.append(
                executor.submit(
                    process_method,
                    method,
                    input_dir,
                    output_dir,
                    max_workers,
                    percentage
                )
            )
        
        # for future in as_completed(futures):
        #     try:
        #         future.result()
        #     except Exception as e:
        #         print(f"Method execution failed: {e}", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform feature selection on multiple files using embedded methods.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help="Path to input data directory")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="Path to output data directory")
    parser.add_argument('-m', '--methods', type=str, nargs='+', required=True, help="List of methods to use (e.g. 'BOR')")
    parser.add_argument('-w', '--workers', type=int, default=1, help="Number of threads for parallel processing")
    parser.add_argument('-p', '--percentage', type=float, default=0.2, help="Percentage of features to select")

    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.methods, args.workers, args.percentage)
