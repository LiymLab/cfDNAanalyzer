import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder

def load_data(file_path, transpose=False):
    """
    Load a CSV file and return data, feature matrix X, labels y, and feature names.

    Parameters:
    - file_path: Path to the CSV file.
    - transpose: Whether to transpose the data.

    Returns:
    - data: Original data (including index and labels).
    - X: Feature matrix (pandas DataFrame).
    - y: Label vector (pandas Series).
    - feature_names: List of feature names.
    """
    try:
        # Read the CSV file, assuming the first column is the index
        data = pd.read_csv(file_path, index_col=0)
        
        if transpose:
            data = data.T
        
        # Check if the 'label' column exists
        if 'label' not in data.columns:
            raise KeyError("Missing 'label' column.")
        
        y = data['label']
        X = data.drop(columns=['label'])
        feature_names = X.columns.tolist()

        # Convert all features to numeric, set non-convertible values to NaN
        X = X.apply(pd.to_numeric, errors='coerce')

        # Apply label encoding to any remaining non-numeric features
        for col in X.columns:
            if X[col].dtype == 'object':
                le = LabelEncoder()
                X[col] = le.fit_transform(X[col].astype(str))

        # print(f"Data types after conversion:\n{X.dtypes}", flush=True)

        # Drop columns that contain any missing values
        initial_feature_count = X.shape[1]
        X.dropna(axis=1, how='any', inplace=True)
        final_feature_count = X.shape[1]
        dropped_features = initial_feature_count - final_feature_count

        if dropped_features > 0:
            print(f"Dropped columns containing missing values: {initial_feature_count} -> {final_feature_count} (dropped {dropped_features} columns)", flush=True)

        # Update feature_names
        feature_names = X.columns.tolist()

        # Ensure all features are numeric
        if not np.all([np.issubdtype(dtype, np.number) for dtype in X.dtypes]):
            raise ValueError("There are non-numeric features in the dataset.")

        return data, X, y, feature_names
    except KeyError as e:
        raise KeyError(f"Missing expected column: {e}")
    except ValueError as ve:
        raise ValueError(f"Invalid data types in dataset: {ve}")
    except Exception as e:
        raise Exception(f"Error loading data: {e}")