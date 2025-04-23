import os
import sys
import pandas as pd
import numpy as np
import argparse
from sklearn.feature_selection import mutual_info_classif, chi2
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from skrebate import ReliefF, SURF, MultiSURF, TuRF
from sklearn_relief import Relief
from scipy.sparse import csc_matrix, lil_matrix, diags
from load_data import load_data
from multiprocessing import Pool, cpu_count

# Missing Ratio
# def calculate_missing_ratios(X):
#     if isinstance(X, pd.DataFrame):
#         feature_names = X.columns
#         missing_ratios = X.isnull().mean()
#         sorted_missing_ratios = sorted(zip(feature_names, missing_ratios), key=lambda x: x[1])
#         return sorted_missing_ratios

# Correlation Coefficient
def calculate_highest_correlation_per_feature(X):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
        corr_matrix = X.corr().abs()
        highest_corr = corr_matrix.apply(lambda x: sorted(x[x < 1.0])[-1] if len(x[x < 1.0]) > 0 else 0, axis=1)
        sorted_highest_corr = sorted(zip(feature_names, highest_corr), key=lambda x: x[1], reverse=True)
        return sorted_highest_corr

# Information Gain
def calculate_information_gain(X, y):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
    mi_scores = mutual_info_classif(X, y)
    sorted_mi_scores = sorted(zip(feature_names, mi_scores), key=lambda x: x[1], reverse=True)
    return sorted_mi_scores

# Chi-square Test
def sorted_chi_square_scores(X, y):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
        X = X.values
    X_scaled = MinMaxScaler().fit_transform(X)
    chi_scores, p_values = chi2(X_scaled, y)
    sorted_scores = sorted(zip(feature_names, chi_scores), key=lambda x: x[1], reverse=True)
    return sorted_scores

# Permutation Feature Importance
def calculate_permutation_importance(model, X, y, n_repeats=10):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
    if isinstance(y, pd.Series):
        y = y.values
    result = permutation_importance(model, X, y, n_repeats=n_repeats)
    feature_importances = result.importances_mean
    sorted_importances = sorted(zip(feature_names, feature_importances), key=lambda x: x[1], reverse=True)
    return sorted_importances

# Low Variance Filter
def low_variance_filter(X, threshold=0.01):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
    selector = VarianceThreshold(threshold=threshold)
    selector.fit(X)
    variances = selector.variances_
    sorted_variances = sorted(zip(feature_names, variances), key=lambda x: x[1], reverse=True)
    return sorted_variances

# Mutual Information
def calculate_mutual_information(X, y):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
    mi_scores = mutual_info_classif(X, y)
    sorted_mi_scores = sorted(zip(feature_names, mi_scores), key=lambda x: x[1], reverse=True)
    return sorted_mi_scores

# Mean Absolute Difference
def mean_absolute_difference(X):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
    mad_scores = X.mad()
    sorted_mad_scores = sorted(zip(feature_names, mad_scores), key=lambda x: x[1], reverse=True)
    return sorted_mad_scores

# Dispersion Ratio
def dispersion_ratio(X, y):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
    overall_variance = np.var(X, axis=0)
    classes = np.unique(y)
    within_class_variance = np.zeros(X.shape[1])
    for cls in classes:
        class_samples = X[y == cls]
        class_variance = np.var(class_samples, axis=0)
        proportion = len(class_samples) / len(X)
        within_class_variance += proportion * class_variance

    overall_variance[overall_variance == 0] = np.finfo(float).eps
    ratios = within_class_variance / overall_variance
    sorted_ratios = sorted(zip(feature_names, ratios), key=lambda x: x[1], reverse=True)
    return sorted_ratios

# FCBF
def symmetricalUncertain(x, y):
    n = float(y.size)
    vals_y = np.unique(y)
    Hx = entropy(x)
    Hy = entropy(y)
    Hxy = sum(entropy(x[y == val]) * (y == val).sum() / n for val in vals_y)
    IG = Hx - Hxy
    return 2 * IG / (Hx + Hy)

def entropy(x):
    vals, counts = np.unique(x, return_counts=True)
    probabilities = counts / counts.sum()
    return -np.sum(probabilities * np.log2(probabilities))

class FCBF:
    def __init__(self, threshold=0.01):
        self.threshold = threshold
        self.selected_features = []
        self.feature_importances = []

    def fit(self, X, y):
        if isinstance(X, pd.DataFrame):
            X = X.values
        if isinstance(y, pd.Series):
            y = y.values

        n_features = X.shape[1]
        su = np.array([symmetricalUncertain(X[:, i], y) for i in range(n_features)])
        self.feature_importances = su  # Store the feature importances
        sorted_indices = np.argsort(su)[::-1]
        selected = sorted_indices[su[sorted_indices] >= self.threshold]

        while len(selected) > 0:
            current_feature = selected[0]
            self.selected_features.append(current_feature)
            selected = selected[1:]
            selected = [idx for idx in selected if symmetricalUncertain(X[:, current_feature], X[:, idx]) < su[current_feature]]

    def fit_transform(self, X, y, feature_names):
        self.fit(X, y)
        feature_importances_with_names = sorted(zip(feature_names, self.feature_importances), key=lambda x: x[1], reverse=True)
        return feature_importances_with_names

# Fisher's score
def construct_W(X, y):
    n_samples = X.shape[0]
    labels = np.unique(y)
    W = lil_matrix((n_samples, n_samples), dtype=np.float64) 
    for label in labels:
        indices = np.where(y == label)[0]
        class_size = len(indices)
        for i in indices:
            for j in indices:
                W[i, j] = 1.0 / class_size if i != j else 1
    return W.tocsc() 

def fisher_score(X, y):
    if isinstance(X, pd.DataFrame):
        feature_names = X.columns
        X = X.values

    W = construct_W(X, y)
    D = np.array(W.sum(axis=0)).flatten()
    D_matrix = diags(D)

    Xt = X.T
    weighted_mean = (Xt @ D_matrix @ np.ones(X.shape[0]) / D.sum()).reshape(-1, 1)
    fr_hat = Xt - weighted_mean
    
    fr_hat_D = fr_hat * D
    num = np.sum(fr_hat_D * fr_hat, axis=1)
    
    L = D_matrix - W
    fr_hat_L = fr_hat * L.diagonal()
    den = np.sum(fr_hat_L * fr_hat, axis=1)
    
    scores = num / np.maximum(den, 1e-12)
    sorted_scores = sorted(zip(feature_names, scores), key=lambda x: x[1], reverse=True)

    return sorted_scores

# Relief methods
def calculate_feature_importance(features, labels, algorithm, feature_names):
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

    if isinstance(X_train, pd.DataFrame):
        X_train = X_train.values
    if isinstance(y_train, pd.Series):
        y_train = y_train.values

    if algorithm == 'ReliefF':
        fs = ReliefF()
    elif algorithm == 'SURF':
        fs = SURF()
    elif algorithm == 'MultiSURF':
        fs = MultiSURF()
    else:
        raise ValueError("Invalid algorithm!")

    fs.fit(X_train, y_train)
    feature_importances = fs.feature_importances_
    sorted_features = sorted(zip(feature_names, feature_importances), key=lambda x: x[1], reverse=True)

    return sorted_features

# TuRF-ReliefF
def calculate_turf_relieff_feature_importance(features, labels, feature_names):
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

    if isinstance(X_train, pd.DataFrame):
        X_train = X_train.values
    if isinstance(y_train, pd.Series):
        y_train = y_train.values

    fs = TuRF(core_algorithm='ReliefF')
    fs.fit(X_train, y_train, headers=feature_names) 

    feature_importances_with_names = sorted(zip(feature_names, fs.feature_importances_), key=lambda x: x[1], reverse=True)
    return feature_importances_with_names

# Main function
def process_single_file(args):
    input_file, output_dir, methods, percentages = args
    try:
        # Get the filename (without the path)
        filename = os.path.basename(input_file)
        # Create output file prefix (remove the extension)
        output_file_prefix = os.path.splitext(filename)[0]
        
        # Load data
        data, X, y, feature_names = load_data(input_file, transpose=False)  # Load data

        percentages_local = []

        if len(percentages) == len(methods):
            percentages_local = percentages
        elif len(percentages) == 1:
            percentages_local = percentages * len(methods)
        elif len(percentages) == 0:
            percentages_local = [1.0] * len(methods) 
        else:
            raise ValueError("Number of percentages must match number of methods, be 1, or be empty.")

        all_results = {}

        for method, percentage in zip(methods, percentages_local):
            if method not in ['CC', 'IG', 'CHI', 'PI', 'LVF', 'MI', 'DR', 'MAD', 'FCBF', 'FS', 'RLF', 'SURF', 'MSURF', 'TRF']:
                raise ValueError(f"Method {method} is not supported.")
            
            # if method == 'MR':
            #     results = calculate_missing_ratios(X)
            elif method == 'CC':
                results = calculate_highest_correlation_per_feature(X)
            elif method == 'IG':
                results = calculate_information_gain(X, y)
            elif method == 'CHI':
                results = sorted_chi_square_scores(X, y)
            elif method == 'PI':
                model = RandomForestClassifier(random_state=42)
                model.fit(X, y)
                results = calculate_permutation_importance(model, X, y)
            elif method == 'LVF':
                results = low_variance_filter(X)
            elif method == 'MI':
                results = calculate_mutual_information(X, y)
            elif method == 'DR':
                results = dispersion_ratio(X, y)
            elif method == 'MAD':
                results = mean_absolute_difference(X)
            elif method == 'FCBF':
                fcbf = FCBF(threshold=0.01)
                results = fcbf.fit_transform(X, y, feature_names)
            elif method == 'FS':
                results = fisher_score(X, y)
            elif method == 'RLF':
                results = calculate_feature_importance(X, y, 'ReliefF', feature_names)
            elif method == 'SURF':
                results = calculate_feature_importance(X, y, 'SURF', feature_names)
            elif method == 'MSURF':
                results = calculate_feature_importance(X, y, 'MultiSURF', feature_names)
            elif method == 'TRF':
                results = calculate_turf_relieff_feature_importance(X, y, feature_names)
            else:
                raise NotImplementedError(f"Method {method} is not implemented yet.")

            # Determine the number of top features to select
            total_features = len(X.columns)
            threshold = max(1, round(total_features * percentage))  # Ensure at least one feature is selected
            selected_features = [feature for feature, score in results[:threshold]]

            # Create the final dataset with selected features
            final_X = X[selected_features]
            final_data = pd.concat([data.index.to_series(name="Sample"), pd.DataFrame(y, columns=["label"]), final_X], axis=1)

            # Save the final dataset with method suffix
            output_file = os.path.join(output_dir, f"{output_file_prefix}_{method}_selected.csv")
            final_data.to_csv(output_file, index=False)
            print(f"[{filename}] Top features selected by {method} saved to {output_file}")

            # Store feature importance ranking results for all methods
            all_results[method] = results

        # Save feature importance rankings for all methods to a file, adding file prefix
        # output_file_importance = os.path.join(output_dir, f"{output_file_prefix}_importance_scores.csv")
        # with open(output_file_importance, 'w', encoding='utf-8') as f:
        #     for method, results in all_results.items():
        #         f.write(f"Method: {method}\n")
        #         if isinstance(results, pd.Series):
        #             for feature, score in results.items():
        #                 f.write(f"{feature},{score}\n")
        #         else:
        #             for feature, score in results:
        #                 f.write(f"{feature},{score}\n")
        #         f.write("\n")
        # print(f"[{filename}] All methods' scores saved to {output_file_importance}")
    
    except Exception as e:
        print(f"Error processing file {input_file}: {e}")

def main(input_dir, output_dir, methods, percentages):
    # Check if the input directory exists
    if not os.path.isdir(input_dir):
        print(f"Input directory {input_dir} does not exist. Please check the path.")
        sys.exit(1)
    
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all files in the input directory (assuming all files need to be processed)
    input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
    
    if not input_files:
        print(f"No files to process in input directory {input_dir}.")
        sys.exit(1)
    
    # Prepare arguments for multiprocessing
    pool_args = [(file, output_dir, methods, percentages) for file in input_files]
    
    # Use multiprocessing to speed up processing
    num_processes = min(cpu_count(), 8)  # Adjust the maximum number of processes as needed
    with Pool(processes=num_processes) as pool:
        pool.map(process_single_file, pool_args)
    
    print(f"All files processed, results saved in the {output_dir} directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter methods for feature selection on multiple files.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help="Path to the input data directory.")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="Path to the output data directory.")
    parser.add_argument('-m', '--methods', type=str, nargs='+', required=True, help="List of methods to use (e.g., 'IG CHI').")
    parser.add_argument('-p', '--percentages', type=float, nargs='*', default=[], help="List of percentages of top features to select for each method (0-1). Default: 1.0")

    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.methods, args.percentages)
