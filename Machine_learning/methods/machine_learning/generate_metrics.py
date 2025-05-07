from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, classification_report
import numpy as np

def evaluate_classification(prediction_df, class_num=2, verbose=True):
    y_true = prediction_df['TrueLabel'].values

    if class_num == 2:
        y_prob = prediction_df['Prob_Class1'].values
        y_pred = (y_prob >= 0.5).astype(int)

        acc = accuracy_score(y_true, y_pred)
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)

        if len(np.unique(y_true)) < 2:
            auc = None
            if verbose:
                print("Warning: Only one class present in y_true. AUC is undefined.")
        else:
            auc = roc_auc_score(y_true, y_prob)

        acc = round(acc, 6)
        precision = round(precision, 6)
        recall = round(recall, 6)
        f1 = round(f1, 6)
        auc = round(auc, 6) if class_num == 2 and len(set(y_true)) > 1 else None

        if verbose:
            print(f"Accuracy: {acc:.6f}")
            print(f"Precision: {precision:.6f}")
            print(f"Recall: {recall:.6f}")
            print(f"F1 Score: {f1:.6f}")
            if auc is not None:
                print(f"AUC: {auc:.6f}")

        return {
            'accuracy': acc,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'auc': auc if auc is not None else "N/A"
        }

    else:
        probs = prediction_df[[col for col in prediction_df.columns if col.startswith("Prob_Class")]].values
        y_pred = np.argmax(probs, axis=1)
        acc = accuracy_score(y_true, y_pred)
        f1_macro = f1_score(y_true, y_pred, average='macro', zero_division=0)
        precision_macro = precision_score(y_true, y_pred, average='macro', zero_division=0)
        recall_macro = recall_score(y_true, y_pred, average='macro', zero_division=0)

        acc = round(acc, 6)
        f1_macro = round(f1_macro, 6)
        precision_macro = round(precision_macro, 6)
        recall_macro = round(recall_macro, 6)

        if verbose:
            print(f"Accuracy: {acc:.6f}")
            print(f"Precision (macro): {precision_macro:.6f}")
            print(f"Recall (macro): {recall_macro:.6f}")
            print(f"F1 Score (macro): {f1_macro:.6f}")

        return {
            'accuracy': acc,
            'precision_macro': precision_macro,
            'recall_macro': recall_macro,
            'f1_macro': f1_macro
        }
