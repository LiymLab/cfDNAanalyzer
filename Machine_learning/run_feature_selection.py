"""
This script controls feature selection logic.

1. Parse arguments (FS method, data path, output path)
2. Load data via load_data (external)
3. Based on fs_type, delegate to corresponding module:
   - filter -> filter.py
   - embedded -> embedded.py
   - wrapper -> wrapper.py
   - hybrid -> FW or FE
4. Output:
   - Return ranked feature-score list or
   - Return selected top-N feature index list for use in modeling
"""

from methods.feature_selection import filter, embedded, wrapper
import pandas as pd

# --- Method maps ---
FILTER_METHODS = {
    "IG": filter.calculate_information_gain,
    "CHI": filter.sorted_chi_square_scores,
    "FS": filter.fisher_score,
    "FCBF": lambda X, y: filter.FCBF(threshold=0.01).fit_transform(X, y),
    "PI": filter.calculate_permutation_importance,
    "CC": filter.calculate_highest_correlation_per_feature,
    "LVF": filter.low_variance_filter,
    "MAD": filter.mean_absolute_difference,
    "DR": filter.dispersion_ratio,
    "MI": filter.calculate_mutual_information,
    "RLF": lambda X, y: filter.calculate_relief_feature_importance(X, y, method="ReliefF"),
    "SURF": lambda X, y: filter.calculate_relief_feature_importance(X, y, method="SURF"),
    "MSURF": lambda X, y: filter.calculate_relief_feature_importance(X, y, method="MultiSURF"),
    "TRF": lambda X, y: filter.calculate_turf_relieff_feature_importance(X, y)
}

EMBEDDED_METHODS = {
    'LASSO': lambda X, y: embedded.apply_regularization(X, y, alpha=1.0, method='LASSO'),
    'RIDGE': lambda X, y: embedded.apply_regularization(X, y, alpha=1.0, method='RIDGE'),
    'ELASTICNET': lambda X, y: embedded.apply_regularization(X, y, alpha=0.5, l1_ratio=0.5, method='ELASTICNET'),
    'RF': lambda X, y: embedded.random_forest_importance(X, y, n_estimators=100, random_state=42)
}

WRAPPER_METHODS = {
    'BOR': wrapper.boruta_feature_selection,
    'FFS': wrapper.forward_feature_selection,
    'BFS': wrapper.backward_feature_selection,
    'EFS': wrapper.exhaustive_feature_selection,
    'RFS': wrapper.recursive_feature_elimination
}

def infer_fs_type(args):
    """Auto-detect feature selection type from args"""
    if args.hybridType in ['FE', 'FW'] and args.hybridMethod1 and args.hybridMethod2:
        return args.hybridType
    elif args.wrapperMethod and len(args.wrapperMethod) > 0:
        return 'wrapper'
    elif args.embeddedMethod and len(args.embeddedMethod) > 0:
        return 'embedded'
    elif args.filterMethod and len(args.filterMethod) > 0:
        return 'filter'
    else:
        raise ValueError("No valid feature selection method provided. Please specify at least one: filterMethod, embeddedMethod, wrapperMethod, or hybridType.")

def build_fs_combinations_from_args(args):
    """
    Build list of feature selection combinations based on args
    Supports auto-matching methods and fractions
    """
    fs_combinations = []

    def parse_fs_combo(methods, fracs, method_type):
        if not methods:
            return []
        if isinstance(fracs, (float, int)):
            fracs = [fracs]
        if not isinstance(methods, list):
            methods = [methods]
        if len(fracs) == 1:
            return [(m, fracs[0]) for m in methods]
        elif len(fracs) == len(methods):
            return zip(methods, fracs)
        else:
            raise ValueError(f"For {method_type}, number of methods ({len(methods)}) and fractions ({len(fracs)}) must match, or only one fraction can be provided.")

    # Filter FS combinations
    for filt, frac in parse_fs_combo(args.filterMethod, args.filterFrac, "filter"):
        fs_combinations.append({
            "filterMethod": [filt], "filterFrac": frac,
            "wrapperMethod": [], "wrapperFrac": 0.0,
            "embeddedMethod": [], "embeddedFrac": 0.0,
            "hybridType": "", "hybridMethod1": "", "hybridMethod2": "",
            "hybridFrac1": 0.0, "hybridFrac2": 0.0,
            "fs_label": f"filter_{filt}_{frac}"
        })

    # Wrapper FS combinations
    for wrap, frac in parse_fs_combo(args.wrapperMethod, args.wrapperFrac, "wrapper"):
        fs_combinations.append({
            "filterMethod": [], "filterFrac": 0.0,
            "wrapperMethod": [wrap], "wrapperFrac": frac,
            "embeddedMethod": [], "embeddedFrac": 0.0,
            "hybridType": "", "hybridMethod1": "", "hybridMethod2": "",
            "hybridFrac1": 0.0, "hybridFrac2": 0.0,
            "fs_label": f"wrapper_{wrap}_{frac}"
        })

    # Embedded FS combinations
    for emb, frac in parse_fs_combo(args.embeddedMethod, args.embeddedFrac, "embedded"):
        fs_combinations.append({
            "filterMethod": [], "filterFrac": 0.0,
            "wrapperMethod": [], "wrapperFrac": 0.0,
            "embeddedMethod": [emb], "embeddedFrac": frac,
            "hybridType": "", "hybridMethod1": "", "hybridMethod2": "",
            "hybridFrac1": 0.0, "hybridFrac2": 0.0,
            "fs_label": f"embedded_{emb}_{frac}"
        })

    # Hybrid FS combination
    if args.hybridType and args.hybridMethod1 and args.hybridMethod2:
        fs_combinations.append({
            "filterMethod": [], "filterFrac": 0.0,
            "wrapperMethod": [], "wrapperFrac": 0.0,
            "embeddedMethod": [], "embeddedFrac": 0.0,
            "hybridType": args.hybridType,
            "hybridMethod1": args.hybridMethod1,
            "hybridMethod2": args.hybridMethod2,
            "hybridFrac1": args.hybridFrac1,
            "hybridFrac2": args.hybridFrac2,
            "fs_label": f"hybrid_{args.hybridType}_{args.hybridMethod1}_{args.hybridMethod2}_{args.hybridFrac1}_{args.hybridFrac2}"
        })

    return fs_combinations

def run_feature_selection(X, y, args, fs_type=None, verbose=True, fs_tag="modality"):
    """Core function for running selected feature selection method"""
    if fs_type is None:
        fs_type = infer_fs_type(args)

    selected_features = []
    scored_features = []

    if fs_type == 'filter':
        for method in args.filterMethod:
            if method not in FILTER_METHODS:
                raise ValueError(f"[filter] Unsupported method: {method}")
            ranked = FILTER_METHODS[method](X, y)
            # k = int(len(ranked) * args.filterFrac)
            if ( args.filterFrac >= 1) or (isinstance(args.filterFrac, (int, float)) and 0 < args.filterFrac < 1):
                if ( args.filterFrac >= 1):
                    k = int(args.filterFrac)
                elif (isinstance(args.filterFrac, (int, float)) and 0 < args.filterFrac < 1) :
                    k = int(len(ranked) * args.filterFrac)
            else:
                raise ValueError(" filterNum must be positive integer or float from 0 to 1")        
                
            top_features = [f for f, _ in ranked[:k]]
            # print(top_features)
            selected_features.extend(top_features)
            scored_features.extend(ranked)

    elif fs_type == 'embedded':
        for method in args.embeddedMethod:
            if method not in EMBEDDED_METHODS:
                raise ValueError(f"[filter] Unsupported method: {method}")
            ranked = EMBEDDED_METHODS[method](X, y)
            # k = int(len(ranked) * args.embeddedFrac)
            if ( args.embeddedFrac >= 1) or (isinstance(args.embeddedFrac, (int, float)) and 0 < args.embeddedFrac < 1):
                if (args.embeddedFrac >= 1):
                    k = int(args.embeddedFrac)
                elif (isinstance(args.embeddedFrac, (int, float)) and 0 < args.embeddedFrac < 1) :
                    k = int(len(ranked) * args.embeddedFrac)
            else:
                raise ValueError(" embeddedNum must be positive integer or float from 0 to 1")
            
            top_features = [f for f, _ in ranked[:k]]
            # print(top_features)
            selected_features.extend(top_features)
            scored_features.extend(ranked)

    elif fs_type == 'wrapper':
        for method in args.wrapperMethod:
            if method not in WRAPPER_METHODS:
                raise ValueError(f"[wrapper] Unsupported method: {method}")
            top_features = WRAPPER_METHODS[method](X, y, percentage=args.wrapperFrac)
            # print(top_features)
            selected_features.extend(top_features)

    elif fs_type in ['FE', 'FW']:
        method1 = args.hybridMethod1
        method2 = args.hybridMethod2
        # k1 = int(len(X.columns) * args.hybridFrac1)
        # k2 = int(len(X.columns) * args.hybridFrac2)
        if ( args.hybridFrac1 >= 1) or (isinstance(args.hybridFrac1, (int, float)) and 0 < args.hybridFrac1 < 1):
            if (args.hybridFrac1 >= 1):
                k1 = int(args.hybridFrac1)
            elif (isinstance(args.hybridFrac1, (int, float)) and 0 < args.hybridFrac1 < 1) :
                k1 = int(len(X.columns) * args.hybridFrac1)
        else:
            raise ValueError(" hybridNum1 must be positive integer or float from 0 to 1")
        
        if ( args.hybridFrac2 >= 1) or (isinstance(args.hybridFrac2, (int, float)) and 0 < args.hybridFrac2 < 1):
            if (args.hybridFrac2 >= 1):
                k2 = int(args.hybridFrac2)
            elif (isinstance(args.hybridFrac2, (int, float)) and 0 < args.hybridFrac2 < 1) :
                k2 = int(len(X.columns) * args.hybridFrac2)
        else:
            raise ValueError(" hybridNum2 must be positive integer or float from 0 to 1")
        

        ranked1 = FILTER_METHODS[method1](X, y)
        top_feats = [f for f, _ in ranked1[:k1]]

        if fs_type == 'FE':
            ranked2 = EMBEDDED_METHODS[method2](X[top_feats], y)
            selected_features = [f for f, _ in ranked2[:k2]]
        else:
            selected_features = WRAPPER_METHODS[method2](X[top_feats], y)

    else:
        raise ValueError("Invalid feature selection type.")

    # Remove duplicates while preserving order
    selected_features = list(dict.fromkeys(selected_features))
    feature_index_list = [X.columns.get_loc(f) for f in selected_features if f in X.columns]

    if verbose:
        print(f"\n [{fs_tag}] Selected Features ({len(selected_features)}):")
        for f in selected_features:
            if f in X.columns:
                print(f" - {f} (Index: {X.columns.get_loc(f)})")

    return feature_index_list
