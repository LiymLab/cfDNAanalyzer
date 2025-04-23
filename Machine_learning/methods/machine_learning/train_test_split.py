from sklearn.model_selection import LeaveOneOut, StratifiedKFold

def get_cv(method, nsplit, y):
    """
    :param method: 'LOO' or 'KFold'
    :param nsplit: KFold's fold
    :param y: 标签序列（用于Stratify）
    """
    if method.upper() == 'LOO':
        return LeaveOneOut()
    elif method.upper() == 'KFOLD':
        return StratifiedKFold(n_splits=nsplit, shuffle=True, random_state=42)
    else:
        raise ValueError("Only 'LOO' and 'KFold' are supported.")
