import pandas as pd

# Load data
def load_data(file_path, transpose=False):
    # Default: column---feature, row---sample
    # Tranpose: row---feature, column---sample
    data = pd.read_csv(file_path, index_col=0)

    if transpose:
        data = data.T 

    if 'label' in data.columns:
        y = data['label']
        X = data.drop('label', axis=1)
    else:
        # If there is no "label" column, we need some other way to determine the label column,
        # assuming here that the last column is a label.
        y = data.iloc[:, -1]
        X = data.iloc[:, :-1]

    feature_names = X.columns.tolist()
    return data, X, y, feature_names

# print feature scores
def print_feature_scores(feature_scores, title):
    print(f"\n************ {title} ************")
    for feature, score in sorted(feature_scores, key=lambda item: item[1], reverse=True):
        print(f"{feature}: {title.split()[0]}={score:.4f}")