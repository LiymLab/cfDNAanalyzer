# visualization.py

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, roc_curve, auc
import numpy as np
import os


def plot_confusion_matrix(prediction_df, class_num, title="", save_path=None):
    y_true = prediction_df['TrueLabel'].values
    if class_num == 2:
        y_pred = (prediction_df['Prob_Class1'].values >= 0.5).astype(int)
    else:
        probs = prediction_df[[col for col in prediction_df.columns if col.startswith("Prob_Class")]].values
        y_pred = np.argmax(probs, axis=1)

    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title(f'Confusion Matrix - {title}')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    if save_path:
        plt.savefig(save_path)

def plot_roc_curve(prediction_df, title="", save_path=None):
    y_true = prediction_df['TrueLabel'].values
    y_score = prediction_df['Prob_Class1'].values
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve - {title}')
    plt.legend(loc="lower right")
    if save_path:
        plt.savefig(save_path)
