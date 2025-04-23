import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Lasso, Ridge, ElasticNet
from sklearn.ensemble import RandomForestClassifier

def apply_regularization(X, y, alpha=1.0, l1_ratio=0.5, method='LASSO'):
    """
    适用于 LASSO / RIDGE / ELASTICNET 三种正则回归嵌入式方法
    返回特征重要性排序 [(feature, score)]
    """
    X_scaled = StandardScaler().fit_transform(X)

    if method.upper() == 'LASSO':
        model = Lasso(alpha=alpha)
    elif method.upper() == 'RIDGE':
        model = Ridge(alpha=alpha)
    elif method.upper() == 'ELASTICNET':
        model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
    else:
        raise ValueError(f"Unsupported method: {method}")

    model.fit(X_scaled, y)
    coef_scores = model.coef_
    return sorted(zip(X.columns, coef_scores), key=lambda x: abs(x[1]), reverse=True)


def random_forest_importance(X, y, n_estimators=100, random_state=42):
    """
    嵌入式方法：基于随机森林的特征重要性
    返回特征排序 [(feature, importance_score)]
    """
    model = RandomForestClassifier(n_estimators=n_estimators, random_state=random_state)
    model.fit(X, y)
    importances = model.feature_importances_
    return sorted(zip(X.columns, importances), key=lambda x: x[1], reverse=True)