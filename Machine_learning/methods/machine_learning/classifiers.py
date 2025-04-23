from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from xgboost import XGBClassifier


def get_classifier(name):
    # name = name.upper()
    classifiers = {
        'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
        'KNN': KNeighborsClassifier(n_neighbors=5),
        'GaussianNB': GaussianNB(),
        'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
        'SVM': svm.SVC(probability=True, kernel='rbf', random_state=42),
        'XGB': XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
    }
    if name not in classifiers:
        raise ValueError(f"Unsupported classifier: {name}")
    return classifiers[name]
