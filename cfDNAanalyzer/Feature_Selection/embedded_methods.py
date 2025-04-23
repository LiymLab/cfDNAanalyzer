import os
import pandas as pd
from sklearn.linear_model import Lasso, Ridge, ElasticNet
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from load_data import load_data
import argparse
from multiprocessing import Pool, cpu_count

# Function to apply regularization methods (LASSO, RIDGE, ELASTICNET)
def apply_regularization(X, y, alpha=1.0, l1_ratio=0.5, method='LASSO'):
    """
    Apply regularization method to the data.

    Parameters:
    - X: Feature matrix (pandas DataFrame).
    - y: Target vector.
    - alpha: Regularization strength.
    - l1_ratio: The ElasticNet mixing parameter, with 0 <= l1_ratio <= 1.
    - method: Type of regularization ('LASSO', 'RIDGE', 'ELASTICNET').

    Returns:
    - Sorted list of tuples containing feature names and their coefficients, sorted by absolute value in descending order.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    if method == 'LASSO':
        model = Lasso(alpha=alpha)
    elif method == 'RIDGE':
        model = Ridge(alpha=alpha)
    elif method == 'ELASTICNET':
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
    else:
        raise ValueError("Method must be 'LASSO', 'RIDGE', or 'ELASTICNET'")
    
    model.fit(X_scaled, y)
    feature_names = X.columns
    coef_scores = model.coef_
    sorted_coef_scores = sorted(zip(feature_names, coef_scores), key=lambda x: abs(x[1]), reverse=True)

    return sorted_coef_scores

# Function to calculate Random Forest feature importances
def random_forest_importance(X, y, n_estimators=100, random_state=None):
    """
    Calculate feature importances using Random Forest.

    Parameters:
    - X: Feature matrix (pandas DataFrame).
    - y: Target vector.
    - n_estimators: Number of trees in the forest.
    - random_state: Random state for reproducibility.

    Returns:
    - Sorted list of tuples containing feature names and their importances, sorted in descending order.
    """
    model = RandomForestClassifier(n_estimators=n_estimators, random_state=random_state)
    model.fit(X, y)
    feature_names = X.columns
    importances = model.feature_importances_
    sorted_importances = sorted(zip(feature_names, importances), key=lambda x: x[1], reverse=True)
    
    return sorted_importances

# Function to process a single file with a specific method
def process_single_method(args):
    """
    Process a single method across all files.

    Parameters:
    - args: Tuple containing (method, file_path, output_dir, percentage, classifier_params).

    Returns:
    - None
    """
    method, file_path, output_dir, percentage, classifier_params = args
    method_upper = method.upper()
    
    try:
        # Load data
        data, X, y, feature_names = load_data(file_path, transpose=False)

        # Apply the feature selection method
        if method_upper in ['LASSO', 'RIDGE', 'ELASTICNET']:
            if method_upper == 'ELASTICNET':
                results = apply_regularization(X, y, alpha=0.5, l1_ratio=0.5, method=method_upper)
            else:
                results = apply_regularization(X, y, alpha=1.0, method=method_upper)
        elif method_upper == 'RF':
            results = random_forest_importance(X, y, **classifier_params)
        else:
            raise ValueError(f"Method {method} is not supported.")
        
        # Determine the number of top features to select based on the percentage
        total_features = len(X.columns)
        threshold = max(1, round(total_features * percentage))  # Ensure at least one feature is selected
        selected_features = [feature for feature, score in results[:threshold]]

        # Create the final dataset with selected features
        final_X = X[selected_features]
        final_data = pd.concat([data.index.to_series(name="Sample"), pd.DataFrame(y, columns=["label"]), final_X], axis=1)

        # Generate the output file name
        input_filename = os.path.splitext(os.path.basename(file_path))[0]
        output_file = os.path.join(output_dir, f"{input_filename}_{method_upper}_selected.csv")

        # Save the selected features to a CSV file
        final_data.to_csv(output_file, index=False)
        print(f"[{method_upper}] Top {percentage*100}% features selected and saved to {output_file}", flush=True)

        # Save feature importance rankings for the method to a separate file
        # output_file_importance = os.path.join(output_dir, f"{input_filename}_{method_upper}_importance_scores.csv")
        # with open(output_file_importance, 'w', encoding='utf-8') as f:
        #     f.write(f"Method: {method_upper}\n")
        #     for feature, score in results:
        #         f.write(f"{feature},{score}\n")
        # print(f"[{method_upper}] Feature importance scores saved to {output_file_importance}", flush=True)

    except Exception as e:
        print(f"[{method_upper}] Failed to process {file_path}: {e}", flush=True)

# Function to process all methods across all files
def process_method(method, input_dir, output_dir, percentages, classifier_params):
    """
    Apply a specific feature selection method to all files in the input directory.

    Parameters:
    - method: Feature selection method name.
    - input_dir: Directory containing input files.
    - output_dir: Base directory to save output files.
    - percentages: List of percentages of top features to select for each method.
    - classifier_params: Dictionary of parameters for the classifier.

    Returns:
    - None
    """
    method_upper = method.upper()
    supported_methods = ['LASSO', 'RIDGE', 'ELASTICNET', 'RF']
    
    if method_upper not in supported_methods:
        print(f"Method {method} is not supported.", flush=True)
        return

    # Create a separate output directory for the current method
    # method_output_dir = os.path.join(output_dir, method_upper)
    method_output_dir = output_dir
    
    os.makedirs(method_output_dir, exist_ok=True)

    print(f"=== Starting feature selection using {method_upper} on all files ===", flush=True)

    # Get all files in the input directory
    input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

    if not input_files:
        print(f"No files found in the input directory {input_dir}.", flush=True)
        return

    # Prepare arguments for multiprocessing
    pool_args = []
    for file_path in input_files:
        # Determine the percentage for the current method
        if len(percentages) == 1:
            percentage = percentages[0]
        else:
            # If multiple percentages are provided, use the first one for simplicity
            # You can implement more complex logic if needed
            percentage = percentages[0]
        
        pool_args.append((method, file_path, method_output_dir, percentage, classifier_params))

    # Determine the number of processes to use
    num_processes = min(cpu_count(), 8)  # Adjust based on your system's capabilities

    # Use multiprocessing Pool to process files in parallel
    with Pool(processes=num_processes) as pool:
        pool.map(process_single_method, pool_args)

    print(f"=== Finished feature selection using {method_upper} on all files ===\n", flush=True)

# Main function to orchestrate the feature selection process
def main(input_dir, output_dir, methods, percentages, n_estimators, random_state):
    """
    Main function to apply feature selection methods to all files.

    Parameters:
    - input_dir: Directory containing input files.
    - output_dir: Base directory to save output files.
    - methods: List of feature selection methods to apply.
    - percentages: List of percentages indicating the proportion of top features to select for each method.
    - n_estimators: Number of trees for Random Forest.
    - random_state: Random state for reproducibility.

    Returns:
    - None
    """
    # Ensure the input directory exists
    if not os.path.isdir(input_dir):
        print(f"Input directory {input_dir} does not exist. Please check the path.", flush=True)
        sys.exit(1)

    # Create the base output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Validate percentages
    if len(percentages) not in [0, 1, len(methods)]:
        print("Number of percentages must match number of methods, be 1, or be empty.", flush=True)
        sys.exit(1)
    
    # If no percentages provided, default to 1.0 for all methods
    if len(percentages) == 0:
        percentages = [1.0] * len(methods)
    elif len(percentages) == 1:
        percentages = percentages * len(methods)
    
    # Prepare classifier parameters
    classifier_params = {
        'n_estimators': n_estimators,
        'random_state': random_state
    }

    # Iterate over each specified feature selection method
    for method, percentage in zip(methods, percentages):
        process_method(method, input_dir, output_dir, [percentage], classifier_params)

    # print(f"All methods processed. Results are saved in the {output_dir} directory.", flush=True)

# Entry point of the script
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Embedded methods for feature selection on multiple files.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help="Path to the input data directory.")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="Path to the output data directory.")
    parser.add_argument('-m', '--methods', type=str, nargs='+', required=True, help="List of methods to use (e.g., 'LASSO RF').")
    parser.add_argument('-p', '--percentages', type=float, nargs='*', default=[], help="List of percentages of top features to select for each method (0-1, Default: 1.0)")
    parser.add_argument('--n_estimators', type=int, default=100, help="Number of trees for Random Forest (default: 100)")
    parser.add_argument('--random_state', type=int, default=42, help="Random state for reproducibility (default: 42)")
    
    args = parser.parse_args()

    # Execute the main function with parsed arguments
    main(args.input_dir, args.output_dir, args.methods, args.percentages, args.n_estimators, args.random_state)
