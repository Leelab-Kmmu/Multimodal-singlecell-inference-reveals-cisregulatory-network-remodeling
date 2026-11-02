# -*- coding: utf-8 -*-
##### 图3 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = '图3'
OUT = OUTPUT_ROOT / FIGURE_NAME
ERRORS = []
PANELS = []

def missing(panel, reason):
    ERRORS.append({'panel': panel, 'reason': reason})
    print(f'[{panel}] {reason}')

def attempt(panel, action):
    try:
        return action()
    except Exception as exc:
        missing(panel, f'{type(exc).__name__}: {exc}')
        traceback.print_exc()
        return None

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
from pathlib import Path as performance_Path
import math as performance_math
import warnings as performance_warnings
import numpy as performance_np
import pandas as performance_pd
import matplotlib.pyplot as performance_plt
from matplotlib.patches import Circle as performance_Circle
from sklearn.metrics import roc_curve as performance_roc_curve, precision_recall_curve as performance_precision_recall_curve, roc_auc_score as performance_roc_auc_score, average_precision_score as performance_average_precision_score
try:
    from scipy.stats import ks_2samp as ks_2samp
except Exception:
    ks_2samp = None
performance_MODEL_DIR = performance_Path(str(DATA_DIR / '4.0result' / 'ATACmotif' / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile'))
performance_SUMMARY_PATH = performance_MODEL_DIR / 'model_comparison_summary.csv'
performance_TEST_PRED_PATH = performance_MODEL_DIR / 'predictions_all_folds_test_edges.csv'
performance_ALL_PRED_PATH = performance_MODEL_DIR / 'predictions_all_folds_all_label_edges.csv'
performance_OUT_DIR = performance_MODEL_DIR / 'paper_like_performance_figure'
performance_OUT_PNG = performance_OUT_DIR / 'scMultiomeGRN_paper_like_performance.png'
performance_OUT_PDF = performance_OUT_DIR / 'scMultiomeGRN_paper_like_performance.pdf'
performance_MODE_ORDER = ['rna', 'atac', 'multi']
performance_MODE_LABEL = {'rna': 'RNA', 'atac': 'ATAC', 'multi': 'Multi'}
performance_MODE_COLORS = {'rna': '#4C72B0', 'atac': '#DD8452', 'multi': '#55A868'}
performance_BAR_AUC_COLOR = '#4C72B0'
performance_BAR_AUPR_COLOR = '#DD8452'
performance_POS_COLOR = '#D62728'
performance_NEG_COLOR = '#1F77B4'
performance_HIGH_CONF_MIN_COUNT = 6
performance_FALLBACK_PRED_THRESHOLD = 'fold_median'

def performance_require(path, label):
    if not performance_Path(path).exists():
        raise FileNotFoundError(f'{label} not found: {path}')

def performance_normalize_mode(s):
    return str(s).strip().lower()

def performance_clean_predictions(df):
    df = df.copy()
    if 'mode' not in df.columns:
        raise ValueError('Prediction file must contain column: mode')
    if 'label' not in df.columns:
        raise ValueError('Prediction file must contain column: label')
    if 'y_prob' not in df.columns:
        raise ValueError('Prediction file must contain column: y_prob')
    df['mode'] = df['mode'].map(performance_normalize_mode)
    df['label'] = df['label'].astype(int)
    df['y_prob'] = performance_pd.to_numeric(df['y_prob'], errors='coerce')
    df = df.dropna(subset=['y_prob'])
    return df

def performance_parse_summary(path):
    performance_require(path, 'model_comparison_summary.csv')
    try:
        s = performance_pd.read_csv(path, header=[0, 1], index_col=0)
        if 'mode' in s.index:
            s = s.drop(index='mode')
        s.index = [performance_normalize_mode(x) for x in s.index]
        out = {}
        for mode in s.index:
            out[mode] = {}
            for metric in ['auroc', 'auprc', 'accuracy', 'f1']:
                if (metric, 'mean') in s.columns:
                    out[mode][metric] = float(s.loc[mode, (metric, 'mean')])
                if (metric, 'std') in s.columns:
                    out[mode][metric + '_std'] = float(s.loc[mode, (metric, 'std')])
        return performance_pd.DataFrame(out).T
    except Exception:
        s = performance_pd.read_csv(path)
        first = s.columns[0]
        if first.lower().startswith('unnamed'):
            s = s.rename(columns={first: 'mode'})
        s['mode'] = s['mode'].map(performance_normalize_mode)
        s = s.set_index('mode')
        return s

def performance_mode_subset(df, mode):
    return df[df['mode'] == mode].copy()

def performance_smooth_hist_density(values, bins=performance_np.linspace(0, 1, 101), sigma=2.0):
    values = performance_np.asarray(values, dtype=float)
    values = values[performance_np.isfinite(values)]
    if values.size == 0:
        centers = (bins[:-1] + bins[1:]) / 2
        return (centers, performance_np.zeros_like(centers))
    hist, edges = performance_np.histogram(values, bins=bins, density=False)
    hist = hist.astype(float)
    bin_width = edges[1] - edges[0]
    area = hist.sum() * bin_width
    if area > 0:
        hist = hist / area
    radius = int(4 * sigma + 1)
    x = performance_np.arange(-radius, radius + 1)
    kernel = performance_np.exp(-0.5 * (x / sigma) ** 2)
    kernel = kernel / kernel.sum()
    smooth = performance_np.convolve(hist, kernel, mode='same')
    centers = (edges[:-1] + edges[1:]) / 2
    return (centers, smooth)

def performance_kl_divergence_from_probs(pos, neg, bins=performance_np.linspace(0, 1, 101), eps=1e-12):
    p, _ = performance_np.histogram(pos, bins=bins, density=False)
    q, _ = performance_np.histogram(neg, bins=bins, density=False)
    p = p.astype(float) + eps
    q = q.astype(float) + eps
    p = p / p.sum()
    q = q / q.sum()
    return float(performance_np.sum(p * performance_np.log(p / q)))

def performance_cohens_d(pos, neg):
    pos = performance_np.asarray(pos, dtype=float)
    neg = performance_np.asarray(neg, dtype=float)
    pos = pos[performance_np.isfinite(pos)]
    neg = neg[performance_np.isfinite(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return performance_np.nan
    n1, n0 = (len(pos), len(neg))
    s1 = performance_np.var(pos, ddof=1)
    s0 = performance_np.var(neg, ddof=1)
    pooled = performance_math.sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / (n1 + n0 - 2))
    if pooled == 0:
        return performance_np.nan
    return float((performance_np.mean(pos) - performance_np.mean(neg)) / pooled)

def performance_ks_test(pos, neg):
    if ks_2samp is None:
        return performance_np.nan
    if len(pos) == 0 or len(neg) == 0:
        return performance_np.nan
    return float(ks_2samp(pos, neg).pvalue)

def performance_sci_p(p):
    if performance_pd.isna(p):
        return 'NA'
    if p == 0:
        return '0.0e+00'
    return f'{p:.1e}'

def performance_add_edge_key(df):
    df = df.copy()
    if {'source_id', 'target_id'}.issubset(df.columns):
        df['edge_key'] = df['source_id'].astype(str) + '->' + df['target_id'].astype(str)
    elif {'regulator_tf', 'target_tf'}.issubset(df.columns):
        df['edge_key'] = df['regulator_tf'].astype(str) + '->' + df['target_tf'].astype(str)
    else:
        raise ValueError('Need source_id/target_id or regulator_tf/target_tf to define edge_key')
    return df

def performance_ensure_pred_label(df):
    df = df.copy()
    if 'paper_pred_label' in df.columns:
        df['__pred_label__'] = df['paper_pred_label'].astype(int)
        return df
    if 'pred_label' in df.columns:
        df['__pred_label__'] = df['pred_label'].astype(int)
        return df
    df['__pred_label__'] = 0
    if performance_FALLBACK_PRED_THRESHOLD == 'fixed_0.5':
        df['__pred_label__'] = (df['y_prob'] >= 0.5).astype(int)
        return df
    if 'fold' in df.columns:
        for (mode, fold), sub in df.groupby(['mode', 'fold']):
            cutoff = sub['y_prob'].median()
            df.loc[sub.index, '__pred_label__'] = (sub['y_prob'] >= cutoff).astype(int)
    else:
        for mode, sub in df.groupby('mode'):
            cutoff = sub['y_prob'].median()
            df.loc[sub.index, '__pred_label__'] = (sub['y_prob'] >= cutoff).astype(int)
    return df

def performance_high_confidence_sets(all_pred):
    all_pred = performance_clean_predictions(all_pred)
    all_pred = performance_add_edge_key(all_pred)
    all_pred = performance_ensure_pred_label(all_pred)
    sets = {}
    details = []
    for mode in performance_MODE_ORDER:
        sub = all_pred[all_pred['mode'] == mode].copy()
        if sub.empty:
            sets[mode] = set()
            continue
        n_runs = sub['fold'].nunique() if 'fold' in sub.columns else 1
        threshold = min(performance_HIGH_CONF_MIN_COUNT, n_runs)
        agg = sub.groupby('edge_key').agg(prediction_count=('__pred_label__', 'sum'), mean_score=('y_prob', 'mean'), n_runs=('__pred_label__', 'count')).reset_index()
        high = agg[agg['prediction_count'] >= threshold].copy()
        sets[mode] = set(high['edge_key'].tolist())
        high['mode'] = mode
        high['effective_min_count'] = threshold
        details.append(high)
    if details:
        detail_df = performance_pd.concat(details, axis=0, ignore_index=True)
        detail_df.to_csv(performance_OUT_DIR / 'high_confidence_edge_sets_used_for_venn.csv', index=False)
    return sets

##### 图3A #####
def performance_panel_a_bar(ax, summary):
    x_labels = ['RNA', 'ATAC', 'Multi']
    x = performance_np.arange(len(performance_MODE_ORDER))
    width = 0.36
    auc = []
    aupr = []
    auc_std = []
    aupr_std = []
    for mode in performance_MODE_ORDER:
        if mode in summary.index:
            auc.append(float(summary.loc[mode, 'auroc']))
            aupr.append(float(summary.loc[mode, 'auprc']))
            auc_std.append(float(summary.loc[mode, 'auroc_std']) if 'auroc_std' in summary.columns else 0)
            aupr_std.append(float(summary.loc[mode, 'auprc_std']) if 'auprc_std' in summary.columns else 0)
        else:
            auc.append(performance_np.nan)
            aupr.append(performance_np.nan)
            auc_std.append(0)
            aupr_std.append(0)
    ax.bar(x - width / 2, auc, width, yerr=auc_std, capsize=3, label='AUC', color=performance_BAR_AUC_COLOR, alpha=0.9)
    ax.bar(x + width / 2, aupr, width, yerr=aupr_std, capsize=3, label='AUPR', color=performance_BAR_AUPR_COLOR, alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.set_ylabel('Score')
    ymin = max(0.45, performance_np.nanmin(auc + aupr) - 0.08)
    ax.set_ylim(ymin, 1.02)
    ax.legend(title='Metric', frameon=True)
    ax.grid(axis='y', alpha=0.35)

##### 图3B #####
def performance_panel_b_roc(ax, test_pred):
    test_pred = performance_clean_predictions(test_pred)
    for mode in performance_MODE_ORDER:
        sub = performance_mode_subset(test_pred, mode)
        if sub.empty or sub['label'].nunique() < 2:
            continue
        y = sub['label'].values
        p = sub['y_prob'].values
        fpr, tpr, _ = performance_roc_curve(y, p)
        auc = performance_roc_auc_score(y, p)
        ax.plot(fpr, tpr, label=f'{performance_MODE_LABEL[mode]} (AUC={auc:.3f})', color=performance_MODE_COLORS[mode], lw=2)
    ax.plot([0, 1], [0, 1], linestyle='--', color='#C44E52', lw=1.8)
    ax.set_title('(A) ROC Curve', fontsize=11)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(alpha=0.35)

def performance_panel_b_pr(ax, test_pred):
    test_pred = performance_clean_predictions(test_pred)
    for mode in performance_MODE_ORDER:
        sub = performance_mode_subset(test_pred, mode)
        if sub.empty or sub['label'].nunique() < 2:
            continue
        y = sub['label'].values
        p = sub['y_prob'].values
        precision, recall, _ = performance_precision_recall_curve(y, p)
        aupr = performance_average_precision_score(y, p)
        ax.plot(recall, precision, label=f'{performance_MODE_LABEL[mode]} (AUPR={aupr:.3f})', color=performance_MODE_COLORS[mode], lw=2)
    ax.set_title('(B) Precision-Recall Curve', fontsize=11)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(alpha=0.35)

##### 图3C #####
def performance_panel_c_density(ax, test_pred, mode, star=False):
    test_pred = performance_clean_predictions(test_pred)
    sub = performance_mode_subset(test_pred, mode)
    pos = sub.loc[sub['label'] == 1, 'y_prob'].values
    neg = sub.loc[sub['label'] == 0, 'y_prob'].values
    bins = performance_np.linspace(0, 1, 101)
    cx_pos, dy_pos = performance_smooth_hist_density(pos, bins=bins, sigma=2.0)
    cx_neg, dy_neg = performance_smooth_hist_density(neg, bins=bins, sigma=2.0)
    ax.plot(cx_pos, dy_pos, color=performance_POS_COLOR, lw=2, label='Positive')
    ax.fill_between(cx_pos, 0, dy_pos, color=performance_POS_COLOR, alpha=0.25)
    ax.plot(cx_neg, dy_neg, color=performance_NEG_COLOR, lw=2, label='Negative')
    ax.fill_between(cx_neg, 0, dy_neg, color=performance_NEG_COLOR, alpha=0.25)
    kl = performance_kl_divergence_from_probs(pos, neg, bins=bins)
    ks = performance_ks_test(pos, neg)
    d = performance_cohens_d(pos, neg)
    title_prefix = '★ ' if star else ''
    ax.set_title(f'{title_prefix}{performance_MODE_LABEL[mode]}\nKL={kl:.2f} | KS p={performance_sci_p(ks)}\nd={d:.2f}', fontsize=10)
    ax.set_xlabel('Predicted probability')
    ax.set_ylabel('Density')
    ax.set_xlim(0, 1)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.2)

##### 图3D #####
def performance_draw_venn2(ax, set_a, set_b, label_a, label_b, color_a='#C7B9E8', color_b='#B7C4C7'):
    a_only = len(set_a - set_b)
    b_only = len(set_b - set_a)
    inter = len(set_a & set_b)
    ax.add_patch(performance_Circle((0.42, 0.52), 0.3, color=color_a, alpha=0.55, ec='black', lw=1.0))
    ax.add_patch(performance_Circle((0.62, 0.52), 0.3, color=color_b, alpha=0.55, ec='black', lw=1.0))
    ax.text(0.29, 0.52, str(a_only), ha='center', va='center', fontsize=10)
    ax.text(0.52, 0.52, str(inter), ha='center', va='center', fontsize=10)
    ax.text(0.75, 0.52, str(b_only), ha='center', va='center', fontsize=10)
    ax.annotate(label_a, xy=(0.24, 0.82), xytext=(0.02, 0.95), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=12, ha='left')
    ax.annotate(label_b, xy=(0.78, 0.24), xytext=(0.98, 0.06), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=12, ha='right')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

def performance_draw_venn3(ax, A, B, C, labels, colors=('#C7B9E8', '#B7C4C7', '#B7CDE8')):
    only_A = len(A - B - C)
    only_B = len(B - A - C)
    only_C = len(C - A - B)
    AB = len((A & B) - C)
    AC = len((A & C) - B)
    BC = len((B & C) - A)
    ABC = len(A & B & C)
    ax.add_patch(performance_Circle((0.42, 0.6), 0.27, color=colors[0], alpha=0.55, ec='black', lw=1.0))
    ax.add_patch(performance_Circle((0.62, 0.6), 0.27, color=colors[1], alpha=0.55, ec='black', lw=1.0))
    ax.add_patch(performance_Circle((0.52, 0.4), 0.27, color=colors[2], alpha=0.55, ec='black', lw=1.0))
    ax.text(0.3, 0.63, str(only_A), ha='center', va='center', fontsize=9)
    ax.text(0.74, 0.63, str(only_B), ha='center', va='center', fontsize=9)
    ax.text(0.52, 0.25, str(only_C), ha='center', va='center', fontsize=9)
    ax.text(0.52, 0.7, str(AB), ha='center', va='center', fontsize=9)
    ax.text(0.42, 0.45, str(AC), ha='center', va='center', fontsize=9)
    ax.text(0.62, 0.45, str(BC), ha='center', va='center', fontsize=9)
    ax.text(0.52, 0.52, str(ABC), ha='center', va='center', fontsize=9)
    ax.annotate(labels[0], xy=(0.28, 0.84), xytext=(0.02, 0.96), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=11, ha='left')
    ax.annotate(labels[1], xy=(0.78, 0.83), xytext=(0.98, 0.95), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=11, ha='right')
    ax.annotate(labels[2], xy=(0.67, 0.2), xytext=(0.98, 0.05), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=11, ha='right')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
performance_OUT_DIR = OUT / 'performance' / 'OUT_DIR'
performance_OUT_DIR.mkdir(parents=True, exist_ok=True)
performance_OUT_PDF = OUT / 'performance.pdf'
performance_OUT_PNG = OUT / 'performance.png'
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

def validation_kl_to_uniform(scores, n_bins=validation_N_BINS, eps=1e-12):
    scores = validation_np.asarray(scores, dtype=float)
    scores = scores[validation_np.isfinite(scores)]
    hist, edges = validation_np.histogram(scores, bins=n_bins, range=(0, 1), density=False)
    p = hist.astype(float) + eps
    p = p / p.sum()
    q = validation_np.ones_like(p) / len(p)
    kl = float(validation_np.sum(p * validation_np.log(p / q)))
    entropy = float(-(p[p > 0] * validation_np.log2(p[p > 0])).sum())
    max_entropy = validation_math.log2(n_bins)
    return {'kl': kl, 'entropy_bits': entropy, 'info_gain_bits': max_entropy - entropy, 'hist': hist, 'edges': edges, 'prob': p}

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

def validation_original_summary(df):
    y = df['label_binary'].values
    s = df['pred_score'].values
    thr = validation_get_threshold(y, s)
    out = {'auroc': validation_safe_auc(y, s), 'aupr': validation_safe_aupr(y, s)}
    out.update(validation_binary_metrics(y, s, thr))
    out.update(validation_separation(y, s))
    k_all = validation_kl_to_uniform(s)
    out.update({'kl_score_to_uniform': k_all['kl'], 'entropy_score_bits': k_all['entropy_bits'], 'information_gain_bits': k_all['info_gain_bits'], 'fraction_score_ge_0.8': float(validation_np.mean(s >= 0.8)), 'fraction_score_le_0.2': float(validation_np.mean(s <= 0.2)), 'fraction_uncertain_0.4_0.6': float(validation_np.mean((s >= 0.4) & (s <= 0.6)))})
    out['kl_pos_to_uniform'] = validation_kl_to_uniform(s[y == 1])['kl'] if validation_np.any(y == 1) else validation_np.nan
    out['kl_neg_to_uniform'] = validation_kl_to_uniform(s[y == 0])['kl'] if validation_np.any(y == 0) else validation_np.nan
    validation_pd.DataFrame([out]).to_csv(validation_TABLE_DIR / 'original_prediction_metric_summary.csv', index=False, encoding='utf-8-sig')
    return out

def validation_shuffle_test(df):
    rng = validation_np.random.default_rng(validation_RANDOM_SEED)
    y = df['label_binary'].values.astype(int)
    s = df['pred_score'].values.astype(float)
    thr = validation_get_threshold(y, s)
    rows = []
    for i in range(validation_N_SHUFFLE):
        ys = rng.permutation(y)
        row = {'iter': i + 1, 'auroc': validation_safe_auc(ys, s), 'aupr': validation_safe_aupr(ys, s)}
        row.update(validation_binary_metrics(ys, s, thr))
        row.update(validation_separation(ys, s))
        rows.append(row)
    res = validation_pd.DataFrame(rows)
    res.to_csv(validation_TABLE_DIR / 'shuffle_label_null_metrics.csv', index=False, encoding='utf-8-sig')
    orig = validation_original_summary(df)
    rows = []
    for m in ['auroc', 'aupr', 'accuracy', 'f1', 'mean_separation', 'cohen_d', 'ks_statistic']:
        if m in res.columns and m in orig:
            vals = res[m].dropna().values
            if len(vals) == 0:
                continue
            rows.append({'metric': m, 'observed': orig[m], 'null_mean': float(validation_np.mean(vals)), 'null_sd': float(validation_np.std(vals, ddof=1)), 'null_q025': float(validation_np.quantile(vals, 0.025)), 'null_q975': float(validation_np.quantile(vals, 0.975)), 'empirical_p_null_ge_observed': float((validation_np.sum(vals >= orig[m]) + 1) / (len(vals) + 1)), 'z_score_vs_null': float((orig[m] - validation_np.mean(vals)) / validation_np.std(vals, ddof=1)) if validation_np.std(vals, ddof=1) > 0 else validation_np.nan})
    summ = validation_pd.DataFrame(rows)
    summ.to_csv(validation_TABLE_DIR / 'shuffle_label_test_summary.csv', index=False, encoding='utf-8-sig')
    return (res, summ)

def validation_bootstrap(df):
    rng = validation_np.random.default_rng(validation_RANDOM_SEED)
    y = df['label_binary'].values.astype(int)
    s = df['pred_score'].values.astype(float)
    n = len(y)
    thr = validation_get_threshold(y, s)
    rows = []
    for i in range(validation_N_BOOTSTRAP):
        idx = rng.integers(0, n, size=n)
        yb, sb = (y[idx], s[idx])
        if len(validation_np.unique(yb)) < 2:
            continue
        row = {'iter': i + 1, 'auroc': validation_safe_auc(yb, sb), 'aupr': validation_safe_aupr(yb, sb)}
        row.update(validation_binary_metrics(yb, sb, thr))
        row.update(validation_separation(yb, sb))
        rows.append(row)
    res = validation_pd.DataFrame(rows)
    res.to_csv(validation_TABLE_DIR / 'bootstrap_metrics.csv', index=False, encoding='utf-8-sig')
    rows = []
    for c in res.columns:
        if c == 'iter':
            continue
        vals = validation_pd.to_numeric(res[c], errors='coerce').dropna()
        if vals.empty:
            continue
        rows.append({'metric': c, 'mean': float(vals.mean()), 'sd': float(vals.std(ddof=1)), 'ci_2.5': float(vals.quantile(0.025)), 'ci_50': float(vals.quantile(0.5)), 'ci_97.5': float(vals.quantile(0.975))})
    summ = validation_pd.DataFrame(rows)
    summ.to_csv(validation_TABLE_DIR / 'bootstrap_metric_summary.csv', index=False, encoding='utf-8-sig')
    return (res, summ)

##### 图3E #####
def validation_plot_score_kl(df):
    y = df['label_binary'].values
    s = df['pred_score'].values
    pos, neg = (s[y == 1], s[y == 0])
    k = validation_kl_to_uniform(s)
    kp = validation_kl_to_uniform(pos) if len(pos) else {'kl': validation_np.nan}
    kn = validation_kl_to_uniform(neg) if len(neg) else {'kl': validation_np.nan}
    sep = validation_separation(y, s)
    fig, axes = validation_plt.subplots(1, 2, figsize=(12, 4.6))
    bins = validation_np.linspace(0, 1, validation_N_BINS + 1)
    axes[0].hist(neg, bins=bins, density=True, alpha=0.55, label='Negative / unlabeled')
    axes[0].hist(pos, bins=bins, density=True, alpha=0.55, label='Positive')
    axes[0].set_xlabel('Predicted probability')
    axes[0].set_ylabel('Density')
    axes[0].set_title('Positive and negative score distributions', fontweight='bold')
    axes[0].legend(frameon=False)
    axes[0].text(0.02, 0.98, f"AUROC={validation_safe_auc(y, s):.3f}\nAUPR={validation_safe_aupr(y, s):.3f}\nSeparation={sep.get('mean_separation', validation_np.nan):.3f}\nKS={sep.get('ks_statistic', validation_np.nan):.3f}", transform=axes[0].transAxes, va='top', fontsize=9, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC'))
    centers = (k['edges'][:-1] + k['edges'][1:]) / 2
    uniform = validation_np.ones_like(k['prob']) / len(k['prob'])
    axes[1].plot(centers, k['prob'], marker='o', linewidth=2, label='Observed')
    axes[1].plot(centers, uniform, linestyle='--', linewidth=1.5, label='Uniform baseline')
    axes[1].set_xlabel('Predicted probability bin')
    axes[1].set_ylabel('Probability mass')
    axes[1].set_title('Information shift from uniform baseline', fontweight='bold')
    axes[1].legend(frameon=False)
    axes[1].text(0.02, 0.98, f"KL(all)={k['kl']:.3f}\nKL(pos)={kp['kl']:.3f}\nKL(neg)={kn['kl']:.3f}\nInfo gain={k['info_gain_bits']:.3f} bits", transform=axes[1].transAxes, va='top', fontsize=9, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC'))
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle('Information-theoretic and distributional validation', fontsize=14, fontweight='bold')
    validation_savefig(fig, 'Fig3_3A_score_distribution_and_KL')

##### 图3F #####
def validation_plot_roc_pr_shuffle(df):
    from sklearn.metrics import roc_curve as roc_curve, precision_recall_curve as precision_recall_curve
    rng = validation_np.random.default_rng(validation_RANDOM_SEED)
    y = df['label_binary'].values.astype(int)
    s = df['pred_score'].values.astype(float)
    ys = rng.permutation(y)
    fig, axes = validation_plt.subplots(1, 2, figsize=(11.2, 4.6))
    for yy, lab, lw in [(y, 'Original labels', 2.5), (ys, 'One shuffled label', 1.8)]:
        fpr, tpr, _ = roc_curve(yy, s)
        axes[0].plot(fpr, tpr, linewidth=lw, label=f'{lab} AUC={validation_safe_auc(yy, s):.3f}')
        pr, rc, _ = precision_recall_curve(yy, s)
        axes[1].plot(rc, pr, linewidth=lw, label=f'{lab} AUPR={validation_safe_aupr(yy, s):.3f}')
    axes[0].plot([0, 1], [0, 1], '--', color='gray', linewidth=1)
    axes[0].set_xlabel('False positive rate')
    axes[0].set_ylabel('True positive rate')
    axes[0].set_title('ROC curve', fontweight='bold')
    axes[1].set_xlabel('Recall')
    axes[1].set_ylabel('Precision')
    axes[1].set_title('Precision-Recall curve', fontweight='bold')
    for ax in axes:
        ax.legend(frameon=False)
        ax.grid(alpha=0.2)
    fig.suptitle('Original-label prediction versus shuffled-label null example', fontsize=14, fontweight='bold')
    validation_savefig(fig, 'Fig3_3B_ROC_PR_original_vs_shuffled')

##### 图3G #####
def validation_plot_bootstrap_ci(boot_summ):
    keep = ['auroc', 'aupr', 'accuracy', 'f1', 'mean_separation', 'ks_statistic']
    df = boot_summ[boot_summ['metric'].isin(keep)].copy()
    if df.empty:
        return
    df['metric'] = validation_pd.Categorical(df['metric'], categories=[m for m in keep if m in set(df['metric'])], ordered=True)
    df = df.sort_values('metric')
    y = validation_np.arange(df.shape[0])
    x = df['mean'].values
    xerr = validation_np.vstack([x - df['ci_2.5'].values, df['ci_97.5'].values - x])
    fig, ax = validation_plt.subplots(figsize=(6.8, 4.8))
    ax.errorbar(x, y, xerr=xerr, fmt='o', capsize=4, linewidth=1.8)
    ax.set_yticks(y)
    ax.set_yticklabels(df['metric'])
    ax.set_xlabel('Bootstrap estimate with 95% CI')
    ax.set_title('Bootstrap robustness of model metrics', fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    validation_savefig(fig, 'Fig3_3D_bootstrap_confidence_intervals')
validation_OUT_DIR = OUT / 'validation' / 'OUT_DIR'
validation_OUT_DIR.mkdir(parents=True, exist_ok=True)
validation_TABLE_DIR = OUT / 'validation' / 'TABLE_DIR'
validation_TABLE_DIR.mkdir(parents=True, exist_ok=True)
validation_FIG_DIR = OUT / 'validation' / 'FIG_DIR'
validation_FIG_DIR.mkdir(parents=True, exist_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    import matplotlib.pyplot as plt
    summary = performance_parse_summary(performance_SUMMARY_PATH)
    test = pd.read_csv(performance_TEST_PRED_PATH)
    all_pred = pd.read_csv(performance_ALL_PRED_PATH)
    sets = performance_high_confidence_sets(all_pred)

    def draw_a():
        fig, ax = plt.subplots(figsize=(7, 5))
        performance_panel_a_bar(ax, summary)
        save_panel(fig, 'A')

    def save_panel(fig, label):
        path = OUT / (label + '.png')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
        plt.close(fig)
        panel_file(label, path)

    def draw_b():
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        performance_panel_b_roc(axes[0], test)
        performance_panel_b_pr(axes[1], test)
        save_panel(fig, 'B')

    def draw_c():
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        for ax, mode in zip(axes, ['rna', 'atac', 'multi']):
            performance_panel_c_density(ax, test, mode, star=mode == 'multi')
        save_panel(fig, 'C')

    def draw_d():
        fig, axes = plt.subplots(3, 1, figsize=(7, 15))
        performance_draw_venn2(axes[0], sets.get('multi', set()), sets.get('atac', set()), 'ATAC+scRNA feature', 'ATAC feature')
        performance_draw_venn2(axes[1], sets.get('multi', set()), sets.get('rna', set()), 'ATAC+scRNA feature', 'scRNA feature')
        performance_draw_venn3(axes[2], sets.get('multi', set()), sets.get('atac', set()), sets.get('rna', set()), labels=('ATAC+scRNA feature', 'ATAC feature', 'scRNA feature'))
        save_panel(fig, 'D')
    for label, fn in [('A', draw_a), ('B', draw_b), ('C', draw_c), ('D', draw_d)]:
        attempt(label, fn)

    def validation():
        df, *_ = validation_load_predictions()
        validation_original_summary(df)
        validation_shuffle_test(df)
        validation_plot_score_kl(df)
        panel_file('E', validation_FIG_DIR / 'Fig3_3A_score_distribution_and_KL.png')
        validation_plot_roc_pr_shuffle(df)
        panel_file('F', validation_FIG_DIR / 'Fig3_3B_ROC_PR_original_vs_shuffled.png')
        _, boot = validation_bootstrap(df)
        validation_plot_bootstrap_ci(boot)
        panel_file('G', validation_FIG_DIR / 'Fig3_3D_bootstrap_confidence_intervals.png')
    attempt('E–G', validation)
if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        missing('运行', f'{type(exc).__name__}: {exc}')
        traceback.print_exc()
    finally:
        finish(7)

if __name__ == "__main__" and ERRORS:
    raise SystemExit(1)
