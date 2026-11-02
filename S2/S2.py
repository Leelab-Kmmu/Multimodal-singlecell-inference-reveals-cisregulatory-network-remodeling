# -*- coding: utf-8 -*-
##### S2 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = 'S2'
OUT = OUTPUT_ROOT / FIGURE_NAME
ERRORS = []
PANELS = []

def missing(panel, reason):
    ERRORS.append({'panel': panel, 'reason': reason})
    print(f'[{panel}] {reason}')

def panel_file(label, path):
    path = Path(path)
    if path.exists():
        PANELS.append((label, path))
    else:
        missing(label, '未生成面板文件：' + str(path))

def finish(expected):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / '运行报告.json').write_text(json.dumps({'figure': FIGURE_NAME, 'expected_panels': expected, 'generated': [x[0] for x in PANELS], 'issues': ERRORS}, ensure_ascii=False, indent=2), encoding='utf-8')
    if PANELS:
        import matplotlib.pyplot as plt
        from PIL import Image
        import math
        n = len(PANELS)
        cols = 2 if n > 1 else 1
        rows = math.ceil(n / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(8 * cols, 6 * rows), squeeze=False)
        for ax in axes.flat:
            ax.axis('off')
        for ax, (label, path) in zip(axes.flat, PANELS):
            with Image.open(path) as im:
                ax.imshow(im.convert('RGB'))
            ax.set_title(label, loc='left', fontweight='bold', fontsize=18)
        fig.tight_layout()
        stem = FIGURE_NAME if not ERRORS and len(PANELS) == expected else FIGURE_NAME + '_partial'
        fig.savefig(OUT / (stem + '.png'), dpi=200, bbox_inches='tight')
        fig.savefig(OUT / (stem + '.pdf'), bbox_inches='tight')
        plt.close(fig)
    print('输出目录：', OUT)
    if ERRORS:
        print('部分面板未完成，详见运行报告；未用占位图替代。')
from pathlib import Path as validation_Path
import warnings as validation_warnings
import math as validation_math
import numpy as validation_np
import pandas as validation_pd
import matplotlib.pyplot as validation_plt
validation_warnings.filterwarnings('ignore')
validation_BASE_DIR = validation_Path(str(DATA_DIR))
validation_ATACMOTIF_DIR = validation_BASE_DIR / '4.0result' / 'ATACmotif'
validation_OUT_DIR = validation_ATACMOTIF_DIR / '53_model_rigour_validation_no_baseline'
validation_FIG_DIR = validation_OUT_DIR / 'figures'
validation_TABLE_DIR = validation_OUT_DIR / 'tables'
validation_PREDICTION_CANDIDATES = [validation_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'predictions_all_folds_all_label_edges.csv', validation_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'predictions_all_folds.csv', validation_ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'predictions_all_folds_all_label_edges.csv']
validation_SUMMARY_CANDIDATES = [validation_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'model_comparison_summary.csv', validation_ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'model_comparison_summary.csv']
validation_RANDOM_SEED = 42
validation_N_SHUFFLE = 1000
validation_N_BOOTSTRAP = 1000
validation_N_BINS = 30
validation_FIXED_THRESHOLD = None

def validation_first_existing(paths):
    for p in paths:
        if validation_Path(p).exists():
            return validation_Path(p)
    return None

def validation_savefig(fig, name):
    png = validation_FIG_DIR / f'{name}.png'
    pdf = validation_FIG_DIR / f'{name}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    validation_plt.close(fig)
    print('Saved:', png)

def validation_infer_col(df, candidates, required=True, name='column'):
    lower = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    for c in df.columns:
        cl = str(c).lower()
        if any((cand.lower() in cl for cand in candidates)):
            return c
    if required:
        raise ValueError(f'Cannot infer {name}. Available columns: {list(df.columns)}')
    return None

def validation_safe_auc(y, s):
    from sklearn.metrics import roc_auc_score as roc_auc_score
    y = validation_np.asarray(y).astype(int)
    if len(validation_np.unique(y)) < 2:
        return validation_np.nan
    return float(roc_auc_score(y, s))

def validation_safe_aupr(y, s):
    from sklearn.metrics import average_precision_score as average_precision_score
    y = validation_np.asarray(y).astype(int)
    if len(validation_np.unique(y)) < 2:
        return validation_np.nan
    return float(average_precision_score(y, s))

def validation_get_threshold(y, s):
    from sklearn.metrics import roc_curve as roc_curve
    if validation_FIXED_THRESHOLD is not None:
        return float(validation_FIXED_THRESHOLD)
    if len(validation_np.unique(y)) < 2:
        return 0.5
    fpr, tpr, thr = roc_curve(y, s)
    idx = int(validation_np.nanargmax(tpr - fpr))
    val = float(thr[idx])
    return val if validation_np.isfinite(val) and 0 <= val <= 1 else 0.5

def validation_binary_metrics(y, s, thr):
    from sklearn.metrics import accuracy_score as accuracy_score, precision_score as precision_score, recall_score as recall_score, f1_score as f1_score, matthews_corrcoef as matthews_corrcoef
    pred = (validation_np.asarray(s) >= thr).astype(int)
    y = validation_np.asarray(y).astype(int)
    return {'threshold': float(thr), 'accuracy': float(accuracy_score(y, pred)), 'precision': float(precision_score(y, pred, zero_division=0)), 'recall': float(recall_score(y, pred, zero_division=0)), 'f1': float(f1_score(y, pred, zero_division=0)), 'mcc': float(matthews_corrcoef(y, pred)) if len(validation_np.unique(pred)) > 1 and len(validation_np.unique(y)) > 1 else validation_np.nan}

def validation_separation(y, s):
    from scipy.stats import ks_2samp as ks_2samp, mannwhitneyu as mannwhitneyu
    y = validation_np.asarray(y).astype(int)
    s = validation_np.asarray(s).astype(float)
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return {}
    mean_pos, mean_neg = (float(validation_np.mean(pos)), float(validation_np.mean(neg)))
    pooled = validation_math.sqrt(((len(pos) - 1) * validation_np.var(pos, ddof=1) + (len(neg) - 1) * validation_np.var(neg, ddof=1)) / max(1, len(pos) + len(neg) - 2))
    cohend = (mean_pos - mean_neg) / pooled if pooled > 0 else validation_np.nan
    ks = ks_2samp(pos, neg)
    try:
        mwp = float(mannwhitneyu(pos, neg, alternative='greater').pvalue)
    except Exception:
        mwp = validation_np.nan
    hp, bins = validation_np.histogram(pos, bins=validation_N_BINS, range=(0, 1), density=True)
    hn, _ = validation_np.histogram(neg, bins=bins, density=True)
    overlap = float(validation_np.minimum(hp, hn).sum() * (bins[1] - bins[0]))
    return {'n_pos': int(len(pos)), 'n_neg': int(len(neg)), 'mean_pos': mean_pos, 'mean_neg': mean_neg, 'mean_separation': mean_pos - mean_neg, 'median_pos': float(validation_np.median(pos)), 'median_neg': float(validation_np.median(neg)), 'cohen_d': float(cohend), 'ks_statistic': float(ks.statistic), 'ks_pvalue': float(ks.pvalue), 'mannwhitney_pvalue_pos_greater': mwp, 'hist_overlap_coefficient': overlap}

def validation_load_predictions():
    p = validation_first_existing(validation_PREDICTION_CANDIDATES)
    if p is None:
        raise FileNotFoundError('Prediction file not found. Edit PREDICTION_CANDIDATES.')
    df = validation_pd.read_csv(p)
    label_col = validation_infer_col(df, ['label', 'y_true', 'true_label', 'target', 'is_edge'], name='label column')
    score_col = validation_infer_col(df, ['pred_prob', 'prob', 'score', 'y_prob', 'prediction', 'mean_score'], name='score column')
    fold_col = validation_infer_col(df, ['fold', 'cv_fold', 'split'], required=False, name='fold column')
    out = df.copy()
    out['label_binary'] = validation_pd.to_numeric(out[label_col], errors='coerce').fillna(0).astype(int)
    out['pred_score'] = validation_pd.to_numeric(out[score_col], errors='coerce')
    out = out.dropna(subset=['pred_score'])
    out['pred_score'] = out['pred_score'].clip(0, 1)
    print('Using:', p)
    print('Rows:', out.shape, 'label:', label_col, 'score:', score_col, 'fold:', fold_col)
    return (out, p, label_col, score_col, fold_col)

def validation_per_fold(df, fold_col):
    if fold_col is None or fold_col not in df.columns:
        return validation_pd.DataFrame()
    rows = []
    for fold, sub in df.groupby(fold_col):
        y = sub['label_binary'].values
        s = sub['pred_score'].values
        if len(validation_np.unique(y)) < 2:
            continue
        thr = validation_get_threshold(y, s)
        row = {'fold': fold, 'n': len(sub), 'n_pos': int((y == 1).sum()), 'n_neg': int((y == 0).sum()), 'auroc': validation_safe_auc(y, s), 'aupr': validation_safe_aupr(y, s)}
        row.update(validation_binary_metrics(y, s, thr))
        row.update(validation_separation(y, s))
        rows.append(row)
    res = validation_pd.DataFrame(rows)
    if not res.empty:
        res.to_csv(validation_TABLE_DIR / 'per_fold_metric_summary.csv', index=False, encoding='utf-8-sig')
    return res

##### S2A #####
def validation_plot_per_fold(per_fold_df):
    if per_fold_df.empty:
        return
    metrics = [m for m in ['auroc', 'aupr', 'accuracy', 'f1'] if m in per_fold_df.columns]
    if not metrics:
        return
    fig, ax = validation_plt.subplots(figsize=(7.2, 4.8))
    data = [per_fold_df[m].dropna().values for m in metrics]
    ax.boxplot(data, positions=validation_np.arange(len(metrics)), widths=0.55, patch_artist=True)
    rng = validation_np.random.default_rng(validation_RANDOM_SEED)
    for i, vals in enumerate(data):
        ax.scatter(rng.normal(i, 0.035, len(vals)), vals, s=25, alpha=0.75)
    ax.set_xticks(validation_np.arange(len(metrics)))
    ax.set_xticklabels([m.upper() for m in metrics])
    ax.set_ylabel('Score')
    ax.set_title('Cross-fold metric stability', fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    validation_savefig(fig, 'Fig3_3E_cross_fold_stability')
validation_OUT_DIR = OUT / 'validation' / 'OUT_DIR'
validation_OUT_DIR.mkdir(parents=True, exist_ok=True)
validation_TABLE_DIR = OUT / 'validation' / 'TABLE_DIR'
validation_TABLE_DIR.mkdir(parents=True, exist_ok=True)
validation_FIG_DIR = OUT / 'validation' / 'FIG_DIR'
validation_FIG_DIR.mkdir(parents=True, exist_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, _, _, _, fold_col = validation_load_predictions()
    if fold_col is None:
        raise ValueError('预测表缺fold列，无法计算S2A')
    fold = validation_per_fold(df, fold_col)
    validation_plot_per_fold(fold)
    panel_file('A', validation_FIG_DIR / 'Fig3_3E_cross_fold_stability.png')
    missing('B', '未定位论文共识支持箱图的原始绘图函数；不能以重新设计的统计代码冒充提取源码。')
if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        missing('运行', f'{type(exc).__name__}: {exc}')
        traceback.print_exc()
    finally:
        finish(2)

if __name__ == "__main__" and ERRORS:
    raise SystemExit(1)
