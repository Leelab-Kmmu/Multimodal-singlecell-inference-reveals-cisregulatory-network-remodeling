# -*- coding: utf-8 -*-
"""论文图3"""
from pathlib import Path
import os, sys, types, traceback, json
DATA_DIR = Path(os.environ.get("CRC_DATA_DIR", r"E:\毕业设计\bishedata"))
OUTPUT_ROOT = Path(os.environ.get("CRC_OUTPUT_DIR", r"C:\Users\WangXinyue\Desktop\Wang_CRC_ScientificReports_Clean_Revised_Package (1)\26.10.14\111"))
FIGURE_NAME = '图3'
OUT = OUTPUT_ROOT / FIGURE_NAME
ERRORS = []
PANELS = []

def load_source(name, source):
    module = types.ModuleType(name)
    module.__file__ = __file__
    sys.modules[name] = module
    import ast
    tree = ast.parse(source)
    class Paths(ast.NodeTransformer):
        def visit_Constant(self, node):
            prefix = r"E:\毕业设计\bishedata"
            if isinstance(node.value, str) and node.value.startswith(prefix):
                suffix = node.value[len(prefix):].strip("\\/")
                value = str(DATA_DIR.joinpath(*suffix.replace("\\", "/").split("/"))) if suffix else str(DATA_DIR)
                return ast.copy_location(ast.Constant(value), node)
            return node
    tree = ast.fix_missing_locations(Paths().visit(tree))
    exec(compile(tree, name, "exec"), module.__dict__)
    for key in ["OUT_DIR", "OUT_ROOT", "OUT", "FIG_DIR", "FIG2_DIR", "FIG3_DIR", "FIG4_DIR", "FIG5_DIR", "FIG34_DIR", "FIG35_DIR", "TABLE_DIR", "ENRICH_DIR"]:
        if key in module.__dict__:
            value = OUT / name / key
            value.mkdir(parents=True, exist_ok=True)
            setattr(module, key, value)
    for key in ["OUT_PNG", "OUT_PDF"]:
        if key in module.__dict__: setattr(module, key, OUT / (name + (".png" if key.endswith("PNG") else ".pdf")))
    return module

def missing(panel, reason):
    ERRORS.append({"panel": panel, "reason": reason})
    print(f"[{panel}] {reason}")

def attempt(panel, action):
    try: return action()
    except Exception as exc:
        missing(panel, f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
        return None

def panel_file(label, path):
    path=Path(path)
    if path.exists(): PANELS.append((label, path))
    else: missing(label, "未生成面板文件：" + str(path))

def finish(expected):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "运行报告.json").write_text(json.dumps({"figure":FIGURE_NAME,"expected_panels":expected,"generated":[x[0] for x in PANELS],"issues":ERRORS},ensure_ascii=False,indent=2),encoding="utf-8")
    if PANELS:
        import matplotlib.pyplot as plt
        from PIL import Image
        import math
        n=len(PANELS); cols=2 if n>1 else 1; rows=math.ceil(n/cols)
        fig,axes=plt.subplots(rows,cols,figsize=(8*cols,6*rows),squeeze=False)
        for ax in axes.flat: ax.axis('off')
        for ax,(label,path) in zip(axes.flat,PANELS):
            with Image.open(path) as im: ax.imshow(im.convert('RGB'))
            ax.set_title(label,loc='left',fontweight='bold',fontsize=18)
        fig.tight_layout()
        stem=FIGURE_NAME if not ERRORS and len(PANELS)==expected else FIGURE_NAME+'_partial'
        fig.savefig(OUT/(stem+'.png'),dpi=200,bbox_inches='tight')
        fig.savefig(OUT/(stem+'.pdf'),bbox_inches='tight')
        plt.close(fig)
    print('输出目录：',OUT)
    if ERRORS: print('部分面板未完成，详见运行报告；未用占位图替代。')

# 来源：plot_scMultiomeGRN_paper_like_performance.py
SOURCE_performance = r'''from pathlib import Path
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score
try:
    from scipy.stats import ks_2samp
except Exception:
    ks_2samp = None
MODEL_DIR = Path('E:\\毕业设计\\bishedata\\4.0result\\ATACmotif\\27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile')
SUMMARY_PATH = MODEL_DIR / 'model_comparison_summary.csv'
TEST_PRED_PATH = MODEL_DIR / 'predictions_all_folds_test_edges.csv'
ALL_PRED_PATH = MODEL_DIR / 'predictions_all_folds_all_label_edges.csv'
OUT_DIR = MODEL_DIR / 'paper_like_performance_figure'
OUT_PNG = OUT_DIR / 'scMultiomeGRN_paper_like_performance.png'
OUT_PDF = OUT_DIR / 'scMultiomeGRN_paper_like_performance.pdf'
MODE_ORDER = ['rna', 'atac', 'multi']
MODE_LABEL = {'rna': 'RNA', 'atac': 'ATAC', 'multi': 'Multi'}
MODE_COLORS = {'rna': '#4C72B0', 'atac': '#DD8452', 'multi': '#55A868'}
BAR_AUC_COLOR = '#4C72B0'
BAR_AUPR_COLOR = '#DD8452'
POS_COLOR = '#D62728'
NEG_COLOR = '#1F77B4'
HIGH_CONF_MIN_COUNT = 6
FALLBACK_PRED_THRESHOLD = 'fold_median'

def require(path, label):
    if not Path(path).exists():
        raise FileNotFoundError(f'{label} not found: {path}')

def normalize_mode(s):
    return str(s).strip().lower()

def clean_predictions(df):
    df = df.copy()
    if 'mode' not in df.columns:
        raise ValueError('Prediction file must contain column: mode')
    if 'label' not in df.columns:
        raise ValueError('Prediction file must contain column: label')
    if 'y_prob' not in df.columns:
        raise ValueError('Prediction file must contain column: y_prob')
    df['mode'] = df['mode'].map(normalize_mode)
    df['label'] = df['label'].astype(int)
    df['y_prob'] = pd.to_numeric(df['y_prob'], errors='coerce')
    df = df.dropna(subset=['y_prob'])
    return df

def parse_summary(path):
    """
    兼容 pandas groupby agg 保存出的 MultiIndex CSV：
        ,auroc,auroc,auprc,auprc,...
        ,mean,std,mean,std,...
        mode,,,,...
        atac,...
    """
    require(path, 'model_comparison_summary.csv')
    try:
        s = pd.read_csv(path, header=[0, 1], index_col=0)
        if 'mode' in s.index:
            s = s.drop(index='mode')
        s.index = [normalize_mode(x) for x in s.index]
        out = {}
        for mode in s.index:
            out[mode] = {}
            for metric in ['auroc', 'auprc', 'accuracy', 'f1']:
                if (metric, 'mean') in s.columns:
                    out[mode][metric] = float(s.loc[mode, (metric, 'mean')])
                if (metric, 'std') in s.columns:
                    out[mode][metric + '_std'] = float(s.loc[mode, (metric, 'std')])
        return pd.DataFrame(out).T
    except Exception:
        s = pd.read_csv(path)
        first = s.columns[0]
        if first.lower().startswith('unnamed'):
            s = s.rename(columns={first: 'mode'})
        s['mode'] = s['mode'].map(normalize_mode)
        s = s.set_index('mode')
        return s

def mode_subset(df, mode):
    return df[df['mode'] == mode].copy()

def smooth_hist_density(values, bins=np.linspace(0, 1, 101), sigma=2.0):
    """
    不依赖 seaborn 的平滑密度曲线。
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        centers = (bins[:-1] + bins[1:]) / 2
        return (centers, np.zeros_like(centers))
    hist, edges = np.histogram(values, bins=bins, density=False)
    hist = hist.astype(float)
    bin_width = edges[1] - edges[0]
    area = hist.sum() * bin_width
    if area > 0:
        hist = hist / area
    radius = int(4 * sigma + 1)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel = kernel / kernel.sum()
    smooth = np.convolve(hist, kernel, mode='same')
    centers = (edges[:-1] + edges[1:]) / 2
    return (centers, smooth)

def kl_divergence_from_probs(pos, neg, bins=np.linspace(0, 1, 101), eps=1e-12):
    p, _ = np.histogram(pos, bins=bins, density=False)
    q, _ = np.histogram(neg, bins=bins, density=False)
    p = p.astype(float) + eps
    q = q.astype(float) + eps
    p = p / p.sum()
    q = q / q.sum()
    return float(np.sum(p * np.log(p / q)))

def cohens_d(pos, neg):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return np.nan
    n1, n0 = (len(pos), len(neg))
    s1 = np.var(pos, ddof=1)
    s0 = np.var(neg, ddof=1)
    pooled = math.sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / (n1 + n0 - 2))
    if pooled == 0:
        return np.nan
    return float((np.mean(pos) - np.mean(neg)) / pooled)

def ks_test(pos, neg):
    if ks_2samp is None:
        return np.nan
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float(ks_2samp(pos, neg).pvalue)

def sci_p(p):
    if pd.isna(p):
        return 'NA'
    if p == 0:
        return '0.0e+00'
    return f'{p:.1e}'

def add_edge_key(df):
    df = df.copy()
    if {'source_id', 'target_id'}.issubset(df.columns):
        df['edge_key'] = df['source_id'].astype(str) + '->' + df['target_id'].astype(str)
    elif {'regulator_tf', 'target_tf'}.issubset(df.columns):
        df['edge_key'] = df['regulator_tf'].astype(str) + '->' + df['target_tf'].astype(str)
    else:
        raise ValueError('Need source_id/target_id or regulator_tf/target_tf to define edge_key')
    return df

def ensure_pred_label(df):
    df = df.copy()
    if 'paper_pred_label' in df.columns:
        df['__pred_label__'] = df['paper_pred_label'].astype(int)
        return df
    if 'pred_label' in df.columns:
        df['__pred_label__'] = df['pred_label'].astype(int)
        return df
    df['__pred_label__'] = 0
    if FALLBACK_PRED_THRESHOLD == 'fixed_0.5':
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

def high_confidence_sets(all_pred):
    all_pred = clean_predictions(all_pred)
    all_pred = add_edge_key(all_pred)
    all_pred = ensure_pred_label(all_pred)
    sets = {}
    details = []
    for mode in MODE_ORDER:
        sub = all_pred[all_pred['mode'] == mode].copy()
        if sub.empty:
            sets[mode] = set()
            continue
        n_runs = sub['fold'].nunique() if 'fold' in sub.columns else 1
        threshold = min(HIGH_CONF_MIN_COUNT, n_runs)
        agg = sub.groupby('edge_key').agg(prediction_count=('__pred_label__', 'sum'), mean_score=('y_prob', 'mean'), n_runs=('__pred_label__', 'count')).reset_index()
        high = agg[agg['prediction_count'] >= threshold].copy()
        sets[mode] = set(high['edge_key'].tolist())
        high['mode'] = mode
        high['effective_min_count'] = threshold
        details.append(high)
    if details:
        detail_df = pd.concat(details, axis=0, ignore_index=True)
        detail_df.to_csv(OUT_DIR / 'high_confidence_edge_sets_used_for_venn.csv', index=False)
    return sets

def panel_a_bar(ax, summary):
    x_labels = ['RNA', 'ATAC', 'Multi']
    x = np.arange(len(MODE_ORDER))
    width = 0.36
    auc = []
    aupr = []
    auc_std = []
    aupr_std = []
    for mode in MODE_ORDER:
        if mode in summary.index:
            auc.append(float(summary.loc[mode, 'auroc']))
            aupr.append(float(summary.loc[mode, 'auprc']))
            auc_std.append(float(summary.loc[mode, 'auroc_std']) if 'auroc_std' in summary.columns else 0)
            aupr_std.append(float(summary.loc[mode, 'auprc_std']) if 'auprc_std' in summary.columns else 0)
        else:
            auc.append(np.nan)
            aupr.append(np.nan)
            auc_std.append(0)
            aupr_std.append(0)
    ax.bar(x - width / 2, auc, width, yerr=auc_std, capsize=3, label='AUC', color=BAR_AUC_COLOR, alpha=0.9)
    ax.bar(x + width / 2, aupr, width, yerr=aupr_std, capsize=3, label='AUPR', color=BAR_AUPR_COLOR, alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.set_ylabel('Score')
    ymin = max(0.45, np.nanmin(auc + aupr) - 0.08)
    ax.set_ylim(ymin, 1.02)
    ax.legend(title='Metric', frameon=True)
    ax.grid(axis='y', alpha=0.35)

def panel_b_roc(ax, test_pred):
    test_pred = clean_predictions(test_pred)
    for mode in MODE_ORDER:
        sub = mode_subset(test_pred, mode)
        if sub.empty or sub['label'].nunique() < 2:
            continue
        y = sub['label'].values
        p = sub['y_prob'].values
        fpr, tpr, _ = roc_curve(y, p)
        auc = roc_auc_score(y, p)
        ax.plot(fpr, tpr, label=f'{MODE_LABEL[mode]} (AUC={auc:.3f})', color=MODE_COLORS[mode], lw=2)
    ax.plot([0, 1], [0, 1], linestyle='--', color='#C44E52', lw=1.8)
    ax.set_title('(A) ROC Curve', fontsize=11)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(alpha=0.35)

def panel_b_pr(ax, test_pred):
    test_pred = clean_predictions(test_pred)
    for mode in MODE_ORDER:
        sub = mode_subset(test_pred, mode)
        if sub.empty or sub['label'].nunique() < 2:
            continue
        y = sub['label'].values
        p = sub['y_prob'].values
        precision, recall, _ = precision_recall_curve(y, p)
        aupr = average_precision_score(y, p)
        ax.plot(recall, precision, label=f'{MODE_LABEL[mode]} (AUPR={aupr:.3f})', color=MODE_COLORS[mode], lw=2)
    ax.set_title('(B) Precision-Recall Curve', fontsize=11)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(alpha=0.35)

def panel_c_density(ax, test_pred, mode, star=False):
    test_pred = clean_predictions(test_pred)
    sub = mode_subset(test_pred, mode)
    pos = sub.loc[sub['label'] == 1, 'y_prob'].values
    neg = sub.loc[sub['label'] == 0, 'y_prob'].values
    bins = np.linspace(0, 1, 101)
    cx_pos, dy_pos = smooth_hist_density(pos, bins=bins, sigma=2.0)
    cx_neg, dy_neg = smooth_hist_density(neg, bins=bins, sigma=2.0)
    ax.plot(cx_pos, dy_pos, color=POS_COLOR, lw=2, label='Positive')
    ax.fill_between(cx_pos, 0, dy_pos, color=POS_COLOR, alpha=0.25)
    ax.plot(cx_neg, dy_neg, color=NEG_COLOR, lw=2, label='Negative')
    ax.fill_between(cx_neg, 0, dy_neg, color=NEG_COLOR, alpha=0.25)
    kl = kl_divergence_from_probs(pos, neg, bins=bins)
    ks = ks_test(pos, neg)
    d = cohens_d(pos, neg)
    title_prefix = '★ ' if star else ''
    ax.set_title(f'{title_prefix}{MODE_LABEL[mode]}\nKL={kl:.2f} | KS p={sci_p(ks)}\nd={d:.2f}', fontsize=10)
    ax.set_xlabel('Predicted probability')
    ax.set_ylabel('Density')
    ax.set_xlim(0, 1)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.2)

def draw_venn2(ax, set_a, set_b, label_a, label_b, color_a='#C7B9E8', color_b='#B7C4C7'):
    a_only = len(set_a - set_b)
    b_only = len(set_b - set_a)
    inter = len(set_a & set_b)
    ax.add_patch(Circle((0.42, 0.52), 0.3, color=color_a, alpha=0.55, ec='black', lw=1.0))
    ax.add_patch(Circle((0.62, 0.52), 0.3, color=color_b, alpha=0.55, ec='black', lw=1.0))
    ax.text(0.29, 0.52, str(a_only), ha='center', va='center', fontsize=10)
    ax.text(0.52, 0.52, str(inter), ha='center', va='center', fontsize=10)
    ax.text(0.75, 0.52, str(b_only), ha='center', va='center', fontsize=10)
    ax.annotate(label_a, xy=(0.24, 0.82), xytext=(0.02, 0.95), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=12, ha='left')
    ax.annotate(label_b, xy=(0.78, 0.24), xytext=(0.98, 0.06), arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=12, ha='right')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

def draw_venn3(ax, A, B, C, labels, colors=('#C7B9E8', '#B7C4C7', '#B7CDE8')):
    only_A = len(A - B - C)
    only_B = len(B - A - C)
    only_C = len(C - A - B)
    AB = len((A & B) - C)
    AC = len((A & C) - B)
    BC = len((B & C) - A)
    ABC = len(A & B & C)
    ax.add_patch(Circle((0.42, 0.6), 0.27, color=colors[0], alpha=0.55, ec='black', lw=1.0))
    ax.add_patch(Circle((0.62, 0.6), 0.27, color=colors[1], alpha=0.55, ec='black', lw=1.0))
    ax.add_patch(Circle((0.52, 0.4), 0.27, color=colors[2], alpha=0.55, ec='black', lw=1.0))
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

def main():
    require(SUMMARY_PATH, 'model_comparison_summary.csv')
    require(TEST_PRED_PATH, 'predictions_all_folds_test_edges.csv')
    require(ALL_PRED_PATH, 'predictions_all_folds_all_label_edges.csv')
    summary = parse_summary(SUMMARY_PATH)
    test_pred = pd.read_csv(TEST_PRED_PATH)
    all_pred = pd.read_csv(ALL_PRED_PATH)
    sets = high_confidence_sets(all_pred)
    fig = plt.figure(figsize=(17, 12))
    gs = fig.add_gridspec(nrows=3, ncols=6, height_ratios=[1.15, 1.15, 1.05], wspace=0.35, hspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0:2])
    panel_a_bar(ax_a, summary)
    ax_roc = fig.add_subplot(gs[0, 2:4])
    ax_pr = fig.add_subplot(gs[0, 4:6])
    panel_b_roc(ax_roc, test_pred)
    panel_b_pr(ax_pr, test_pred)
    ax_c1 = fig.add_subplot(gs[1, 0:2])
    ax_c2 = fig.add_subplot(gs[1, 2:4])
    ax_c3 = fig.add_subplot(gs[1, 4:6])
    panel_c_density(ax_c1, test_pred, 'rna', star=False)
    panel_c_density(ax_c2, test_pred, 'atac', star=False)
    panel_c_density(ax_c3, test_pred, 'multi', star=True)
    ax_d1 = fig.add_subplot(gs[2, 0:2])
    ax_d2 = fig.add_subplot(gs[2, 2:4])
    ax_d3 = fig.add_subplot(gs[2, 4:6])
    draw_venn2(ax_d1, sets.get('multi', set()), sets.get('atac', set()), 'ATAC+scRNA feature', 'ATAC feature')
    draw_venn2(ax_d2, sets.get('multi', set()), sets.get('rna', set()), 'ATAC+scRNA feature', 'scRNA feature', color_a='#C7B9E8', color_b='#B7CDE8')
    draw_venn3(ax_d3, sets.get('multi', set()), sets.get('atac', set()), sets.get('rna', set()), labels=('ATAC+scRNA feature', 'ATAC feature', 'scRNA feature'))
    fig.text(0.015, 0.955, 'A', fontsize=28)
    fig.text(0.33, 0.955, 'B', fontsize=28)
    fig.text(0.015, 0.635, 'C', fontsize=28)
    fig.text(0.015, 0.325, 'D', fontsize=28)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
    fig.savefig(OUT_PDF, bbox_inches='tight')
    print('DONE')
    print('Saved PNG:', OUT_PNG)
    print('Saved PDF:', OUT_PDF)
    print('Venn edge set details:', OUT_DIR / 'high_confidence_edge_sets_used_for_venn.csv')
'''

# 来源：model_rigour_validation_no_baseline.py
SOURCE_validation = r'''"""
Model rigorous validation without external deep-learning baselines.

Run:
    D:\\Python\\python.exe E:\\毕业设计\\scmultiome_project\\第四次\\model_rigour_validation_no_baseline.py

Output:
    E:\\毕业设计\x08ishedata\x04.0result\\ATACmotif+_model_rigour_validation_no_baseline

Analyses:
    1. Score distribution + KL divergence against uniform baseline
    2. Positive vs negative/unlabeled score separation
    3. Label-shuffle null test using fixed prediction scores
    4. Bootstrap confidence intervals
    5. RNA / ATAC / Multi modality contribution, if model_comparison_summary.csv exists
"""
from pathlib import Path
import warnings
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
BASE_DIR = Path('E:\\毕业设计\\bishedata')
ATACMOTIF_DIR = BASE_DIR / '4.0result' / 'ATACmotif'
OUT_DIR = ATACMOTIF_DIR / '53_model_rigour_validation_no_baseline'
FIG_DIR = OUT_DIR / 'figures'
TABLE_DIR = OUT_DIR / 'tables'
PREDICTION_CANDIDATES = [ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'predictions_all_folds_all_label_edges.csv', ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'predictions_all_folds.csv', ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'predictions_all_folds_all_label_edges.csv']
SUMMARY_CANDIDATES = [ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'model_comparison_summary.csv', ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'model_comparison_summary.csv']
RANDOM_SEED = 42
N_SHUFFLE = 1000
N_BOOTSTRAP = 1000
N_BINS = 30
FIXED_THRESHOLD = None

def first_existing(paths):
    for p in paths:
        if Path(p).exists():
            return Path(p)
    return None

def savefig(fig, name):
    png = FIG_DIR / f'{name}.png'
    pdf = FIG_DIR / f'{name}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    print('Saved:', png)

def infer_col(df, candidates, required=True, name='column'):
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

def safe_auc(y, s):
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y).astype(int)
    if len(np.unique(y)) < 2:
        return np.nan
    return float(roc_auc_score(y, s))

def safe_aupr(y, s):
    from sklearn.metrics import average_precision_score
    y = np.asarray(y).astype(int)
    if len(np.unique(y)) < 2:
        return np.nan
    return float(average_precision_score(y, s))

def get_threshold(y, s):
    from sklearn.metrics import roc_curve
    if FIXED_THRESHOLD is not None:
        return float(FIXED_THRESHOLD)
    if len(np.unique(y)) < 2:
        return 0.5
    fpr, tpr, thr = roc_curve(y, s)
    idx = int(np.nanargmax(tpr - fpr))
    val = float(thr[idx])
    return val if np.isfinite(val) and 0 <= val <= 1 else 0.5

def binary_metrics(y, s, thr):
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef
    pred = (np.asarray(s) >= thr).astype(int)
    y = np.asarray(y).astype(int)
    return {'threshold': float(thr), 'accuracy': float(accuracy_score(y, pred)), 'precision': float(precision_score(y, pred, zero_division=0)), 'recall': float(recall_score(y, pred, zero_division=0)), 'f1': float(f1_score(y, pred, zero_division=0)), 'mcc': float(matthews_corrcoef(y, pred)) if len(np.unique(pred)) > 1 and len(np.unique(y)) > 1 else np.nan}

def kl_to_uniform(scores, n_bins=N_BINS, eps=1e-12):
    scores = np.asarray(scores, dtype=float)
    scores = scores[np.isfinite(scores)]
    hist, edges = np.histogram(scores, bins=n_bins, range=(0, 1), density=False)
    p = hist.astype(float) + eps
    p = p / p.sum()
    q = np.ones_like(p) / len(p)
    kl = float(np.sum(p * np.log(p / q)))
    entropy = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
    max_entropy = math.log2(n_bins)
    return {'kl': kl, 'entropy_bits': entropy, 'info_gain_bits': max_entropy - entropy, 'hist': hist, 'edges': edges, 'prob': p}

def separation(y, s):
    from scipy.stats import ks_2samp, mannwhitneyu
    y = np.asarray(y).astype(int)
    s = np.asarray(s).astype(float)
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return {}
    mean_pos, mean_neg = (float(np.mean(pos)), float(np.mean(neg)))
    pooled = math.sqrt(((len(pos) - 1) * np.var(pos, ddof=1) + (len(neg) - 1) * np.var(neg, ddof=1)) / max(1, len(pos) + len(neg) - 2))
    cohend = (mean_pos - mean_neg) / pooled if pooled > 0 else np.nan
    ks = ks_2samp(pos, neg)
    try:
        mwp = float(mannwhitneyu(pos, neg, alternative='greater').pvalue)
    except Exception:
        mwp = np.nan
    hp, bins = np.histogram(pos, bins=N_BINS, range=(0, 1), density=True)
    hn, _ = np.histogram(neg, bins=bins, density=True)
    overlap = float(np.minimum(hp, hn).sum() * (bins[1] - bins[0]))
    return {'n_pos': int(len(pos)), 'n_neg': int(len(neg)), 'mean_pos': mean_pos, 'mean_neg': mean_neg, 'mean_separation': mean_pos - mean_neg, 'median_pos': float(np.median(pos)), 'median_neg': float(np.median(neg)), 'cohen_d': float(cohend), 'ks_statistic': float(ks.statistic), 'ks_pvalue': float(ks.pvalue), 'mannwhitney_pvalue_pos_greater': mwp, 'hist_overlap_coefficient': overlap}

def load_predictions():
    p = first_existing(PREDICTION_CANDIDATES)
    if p is None:
        raise FileNotFoundError('Prediction file not found. Edit PREDICTION_CANDIDATES.')
    df = pd.read_csv(p)
    label_col = infer_col(df, ['label', 'y_true', 'true_label', 'target', 'is_edge'], name='label column')
    score_col = infer_col(df, ['pred_prob', 'prob', 'score', 'y_prob', 'prediction', 'mean_score'], name='score column')
    fold_col = infer_col(df, ['fold', 'cv_fold', 'split'], required=False, name='fold column')
    out = df.copy()
    out['label_binary'] = pd.to_numeric(out[label_col], errors='coerce').fillna(0).astype(int)
    out['pred_score'] = pd.to_numeric(out[score_col], errors='coerce')
    out = out.dropna(subset=['pred_score'])
    out['pred_score'] = out['pred_score'].clip(0, 1)
    print('Using:', p)
    print('Rows:', out.shape, 'label:', label_col, 'score:', score_col, 'fold:', fold_col)
    return (out, p, label_col, score_col, fold_col)

def original_summary(df):
    y = df['label_binary'].values
    s = df['pred_score'].values
    thr = get_threshold(y, s)
    out = {'auroc': safe_auc(y, s), 'aupr': safe_aupr(y, s)}
    out.update(binary_metrics(y, s, thr))
    out.update(separation(y, s))
    k_all = kl_to_uniform(s)
    out.update({'kl_score_to_uniform': k_all['kl'], 'entropy_score_bits': k_all['entropy_bits'], 'information_gain_bits': k_all['info_gain_bits'], 'fraction_score_ge_0.8': float(np.mean(s >= 0.8)), 'fraction_score_le_0.2': float(np.mean(s <= 0.2)), 'fraction_uncertain_0.4_0.6': float(np.mean((s >= 0.4) & (s <= 0.6)))})
    out['kl_pos_to_uniform'] = kl_to_uniform(s[y == 1])['kl'] if np.any(y == 1) else np.nan
    out['kl_neg_to_uniform'] = kl_to_uniform(s[y == 0])['kl'] if np.any(y == 0) else np.nan
    pd.DataFrame([out]).to_csv(TABLE_DIR / 'original_prediction_metric_summary.csv', index=False, encoding='utf-8-sig')
    return out

def shuffle_test(df):
    rng = np.random.default_rng(RANDOM_SEED)
    y = df['label_binary'].values.astype(int)
    s = df['pred_score'].values.astype(float)
    thr = get_threshold(y, s)
    rows = []
    for i in range(N_SHUFFLE):
        ys = rng.permutation(y)
        row = {'iter': i + 1, 'auroc': safe_auc(ys, s), 'aupr': safe_aupr(ys, s)}
        row.update(binary_metrics(ys, s, thr))
        row.update(separation(ys, s))
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(TABLE_DIR / 'shuffle_label_null_metrics.csv', index=False, encoding='utf-8-sig')
    orig = original_summary(df)
    rows = []
    for m in ['auroc', 'aupr', 'accuracy', 'f1', 'mean_separation', 'cohen_d', 'ks_statistic']:
        if m in res.columns and m in orig:
            vals = res[m].dropna().values
            if len(vals) == 0:
                continue
            rows.append({'metric': m, 'observed': orig[m], 'null_mean': float(np.mean(vals)), 'null_sd': float(np.std(vals, ddof=1)), 'null_q025': float(np.quantile(vals, 0.025)), 'null_q975': float(np.quantile(vals, 0.975)), 'empirical_p_null_ge_observed': float((np.sum(vals >= orig[m]) + 1) / (len(vals) + 1)), 'z_score_vs_null': float((orig[m] - np.mean(vals)) / np.std(vals, ddof=1)) if np.std(vals, ddof=1) > 0 else np.nan})
    summ = pd.DataFrame(rows)
    summ.to_csv(TABLE_DIR / 'shuffle_label_test_summary.csv', index=False, encoding='utf-8-sig')
    return (res, summ)

def bootstrap(df):
    rng = np.random.default_rng(RANDOM_SEED)
    y = df['label_binary'].values.astype(int)
    s = df['pred_score'].values.astype(float)
    n = len(y)
    thr = get_threshold(y, s)
    rows = []
    for i in range(N_BOOTSTRAP):
        idx = rng.integers(0, n, size=n)
        yb, sb = (y[idx], s[idx])
        if len(np.unique(yb)) < 2:
            continue
        row = {'iter': i + 1, 'auroc': safe_auc(yb, sb), 'aupr': safe_aupr(yb, sb)}
        row.update(binary_metrics(yb, sb, thr))
        row.update(separation(yb, sb))
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(TABLE_DIR / 'bootstrap_metrics.csv', index=False, encoding='utf-8-sig')
    rows = []
    for c in res.columns:
        if c == 'iter':
            continue
        vals = pd.to_numeric(res[c], errors='coerce').dropna()
        if vals.empty:
            continue
        rows.append({'metric': c, 'mean': float(vals.mean()), 'sd': float(vals.std(ddof=1)), 'ci_2.5': float(vals.quantile(0.025)), 'ci_50': float(vals.quantile(0.5)), 'ci_97.5': float(vals.quantile(0.975))})
    summ = pd.DataFrame(rows)
    summ.to_csv(TABLE_DIR / 'bootstrap_metric_summary.csv', index=False, encoding='utf-8-sig')
    return (res, summ)

def per_fold(df, fold_col):
    if fold_col is None or fold_col not in df.columns:
        return pd.DataFrame()
    rows = []
    for fold, sub in df.groupby(fold_col):
        y = sub['label_binary'].values
        s = sub['pred_score'].values
        if len(np.unique(y)) < 2:
            continue
        thr = get_threshold(y, s)
        row = {'fold': fold, 'n': len(sub), 'n_pos': int((y == 1).sum()), 'n_neg': int((y == 0).sum()), 'auroc': safe_auc(y, s), 'aupr': safe_aupr(y, s)}
        row.update(binary_metrics(y, s, thr))
        row.update(separation(y, s))
        rows.append(row)
    res = pd.DataFrame(rows)
    if not res.empty:
        res.to_csv(TABLE_DIR / 'per_fold_metric_summary.csv', index=False, encoding='utf-8-sig')
    return res

def plot_score_kl(df):
    y = df['label_binary'].values
    s = df['pred_score'].values
    pos, neg = (s[y == 1], s[y == 0])
    k = kl_to_uniform(s)
    kp = kl_to_uniform(pos) if len(pos) else {'kl': np.nan}
    kn = kl_to_uniform(neg) if len(neg) else {'kl': np.nan}
    sep = separation(y, s)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    bins = np.linspace(0, 1, N_BINS + 1)
    axes[0].hist(neg, bins=bins, density=True, alpha=0.55, label='Negative / unlabeled')
    axes[0].hist(pos, bins=bins, density=True, alpha=0.55, label='Positive')
    axes[0].set_xlabel('Predicted probability')
    axes[0].set_ylabel('Density')
    axes[0].set_title('Positive and negative score distributions', fontweight='bold')
    axes[0].legend(frameon=False)
    axes[0].text(0.02, 0.98, f"AUROC={safe_auc(y, s):.3f}\nAUPR={safe_aupr(y, s):.3f}\nSeparation={sep.get('mean_separation', np.nan):.3f}\nKS={sep.get('ks_statistic', np.nan):.3f}", transform=axes[0].transAxes, va='top', fontsize=9, bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC'))
    centers = (k['edges'][:-1] + k['edges'][1:]) / 2
    uniform = np.ones_like(k['prob']) / len(k['prob'])
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
    savefig(fig, 'Fig3_3A_score_distribution_and_KL')

def plot_roc_pr_shuffle(df):
    from sklearn.metrics import roc_curve, precision_recall_curve
    rng = np.random.default_rng(RANDOM_SEED)
    y = df['label_binary'].values.astype(int)
    s = df['pred_score'].values.astype(float)
    ys = rng.permutation(y)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    for yy, lab, lw in [(y, 'Original labels', 2.5), (ys, 'One shuffled label', 1.8)]:
        fpr, tpr, _ = roc_curve(yy, s)
        axes[0].plot(fpr, tpr, linewidth=lw, label=f'{lab} AUC={safe_auc(yy, s):.3f}')
        pr, rc, _ = precision_recall_curve(yy, s)
        axes[1].plot(rc, pr, linewidth=lw, label=f'{lab} AUPR={safe_aupr(yy, s):.3f}')
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
    savefig(fig, 'Fig3_3B_ROC_PR_original_vs_shuffled')

def plot_shuffle_null(shuf, orig):
    metrics = [m for m in ['auroc', 'aupr', 'mean_separation', 'ks_statistic'] if m in shuf.columns and m in orig]
    if not metrics:
        return
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4.2))
    if len(metrics) == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics):
        vals = shuf[m].dropna().values
        ax.hist(vals, bins=30, alpha=0.75, edgecolor='white')
        ax.axvline(orig[m], color='#B22222', linewidth=2.5, label='Observed')
        ax.axvline(np.mean(vals), color='black', linestyle='--', linewidth=1.5, label='Null mean')
        ax.set_title(m, fontweight='bold')
        ax.grid(alpha=0.18)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Shuffle-label null distributions', fontsize=14, fontweight='bold')
    savefig(fig, 'Fig3_3C_shuffle_label_null_test')

def plot_bootstrap_ci(boot_summ):
    keep = ['auroc', 'aupr', 'accuracy', 'f1', 'mean_separation', 'ks_statistic']
    df = boot_summ[boot_summ['metric'].isin(keep)].copy()
    if df.empty:
        return
    df['metric'] = pd.Categorical(df['metric'], categories=[m for m in keep if m in set(df['metric'])], ordered=True)
    df = df.sort_values('metric')
    y = np.arange(df.shape[0])
    x = df['mean'].values
    xerr = np.vstack([x - df['ci_2.5'].values, df['ci_97.5'].values - x])
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.errorbar(x, y, xerr=xerr, fmt='o', capsize=4, linewidth=1.8)
    ax.set_yticks(y)
    ax.set_yticklabels(df['metric'])
    ax.set_xlabel('Bootstrap estimate with 95% CI')
    ax.set_title('Bootstrap robustness of model metrics', fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    savefig(fig, 'Fig3_3D_bootstrap_confidence_intervals')

def plot_per_fold(per_fold_df):
    if per_fold_df.empty:
        return
    metrics = [m for m in ['auroc', 'aupr', 'accuracy', 'f1'] if m in per_fold_df.columns]
    if not metrics:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    data = [per_fold_df[m].dropna().values for m in metrics]
    ax.boxplot(data, positions=np.arange(len(metrics)), widths=0.55, patch_artist=True)
    rng = np.random.default_rng(RANDOM_SEED)
    for i, vals in enumerate(data):
        ax.scatter(rng.normal(i, 0.035, len(vals)), vals, s=25, alpha=0.75)
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels([m.upper() for m in metrics])
    ax.set_ylabel('Score')
    ax.set_title('Cross-fold metric stability', fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    savefig(fig, 'Fig3_3E_cross_fold_stability')

def plot_modal_contribution():
    p = first_existing(SUMMARY_CANDIDATES)
    if p is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.axis('off')
        ax.text(0.02, 0.85, 'Modality contribution skipped', fontsize=14, fontweight='bold', transform=ax.transAxes)
        ax.text(0.02, 0.65, 'No model_comparison_summary.csv found.', fontsize=10, transform=ax.transAxes)
        savefig(fig, 'Fig3_4A_modality_contribution_placeholder')
        return None
    df = pd.read_csv(p)
    model_col = infer_col(df, ['model', 'modality', 'name', 'method'], required=False, name='model column')
    if model_col is None:
        model_col = df.columns[0]
    df['model_name'] = df[model_col].astype(str)
    colmap = {}
    for m, cands in {'AUROC': ['auroc', 'auc', 'roc_auc', 'AUC'], 'AUPR': ['aupr', 'auprc', 'average_precision', 'AUPRC'], 'Accuracy': ['accuracy', 'acc'], 'F1': ['f1', 'f1_score']}.items():
        c = infer_col(df, cands, required=False, name=m)
        if c is not None:
            df[m] = pd.to_numeric(df[c], errors='coerce')
            colmap[m] = c
    metrics = [m for m in ['AUROC', 'AUPR', 'Accuracy', 'F1'] if m in df.columns and df[m].notna().any()]
    if not metrics:
        return p
    order = []
    for key in ['RNA', 'ATAC', 'Multi', 'multi', 'multimodal']:
        for x in df['model_name']:
            if key.lower() in str(x).lower() and x not in order:
                order.append(x)
    if order:
        df = df[df['model_name'].isin(order)].copy()
        df['model_name'] = pd.Categorical(df['model_name'], categories=order, ordered=True)
        df = df.sort_values('model_name')
    x = np.arange(df.shape[0])
    width = 0.8 / len(metrics)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for i, m in enumerate(metrics):
        ax.bar(x + (i - (len(metrics) - 1) / 2) * width, df[m], width=width, label=m)
    ax.set_xticks(x)
    ax.set_xticklabels(df['model_name'].astype(str), rotation=25, ha='right')
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Score')
    ax.set_title('Modality contribution analysis', fontweight='bold')
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=0.2)
    df.to_csv(TABLE_DIR / 'modality_contribution_summary_used.csv', index=False, encoding='utf-8-sig')
    savefig(fig, 'Fig3_4A_modality_contribution_RNA_ATAC_Multi')
    return p

def write_templates():
    text = "\n3.3 基于信息论与统计抽样的模型严谨性验证\n\n为评估多模态图模型输出结果的统计可靠性，本研究在不引入外部深度学习对比模型的前提下，构建了基于预测分布、标签置乱和Bootstrap重采样的验证流程。需要说明的是，正样本来源于ATAC motif-promoter先验调控骨架，负样本为当前规则下未检测到先验支持的候选边，因此该评估反映的是模型对当前候选调控先验标签的预测能力，而不是对完整真实GRN的无偏验证。\n\n首先，以均匀分布作为随机基准，计算模型预测概率分布的KL散度、信息增益和熵。如果预测结果接近随机，分数分布应接近均匀或集中于不确定区间；如果模型学习到有效判别信息，预测概率会向低置信负边和高置信正边两端聚集。结果显示，模型预测分布相对于均匀基准存在明显信息偏移，提示预测分数具有非随机结构。\n\n其次，比较正样本边与未标注候选边的预测分数分布，并通过AUROC、AUPR、均值分离度、KS统计量、Cohen's d和分布重叠系数量化二者分离程度。结果显示，正样本预测概率整体高于未标注候选边，说明模型能够将motif-promoter支持的候选调控边与未检测到先验支持的候选边区分开。\n\n进一步采用标签置乱检验构建随机基准。在保持模型预测分数不变的情况下，随机打乱正负标签并重复计算性能指标。真实标签下的指标若显著高于置乱标签基准，说明模型预测分数与调控先验标签之间存在稳定对应关系。最后，通过Bootstrap重采样估计AUROC、AUPR、Accuracy、F1-score及分布分离指标的95%置信区间，用于评估模型指标是否依赖少数样本或单次划分结果。\n\n3.4 基于模态贡献的模型结构合理性分析\n\n为评估RNA与ATAC信息在候选TF-TF调控关系预测中的贡献，本研究比较RNA-only、ATAC-only和multi-modal三种输入模式。RNA-only模型主要反映转录表达协同和GRNBoost2调控重要性信息，ATAC-only模型主要反映染色质可及性和RP调控潜力信息，multi-modal模型同时整合两类组学特征。\n\n结果显示，多模态模型整体优于单模态模型，说明RNA与ATAC在TF-TF候选调控边预测中具有互补作用。考虑到监督标签来源于ATAC motif-promoter先验，ATAC模态性能可能受到标签来源同源性的影响，因此本文将该结果解释为多模态特征对当前候选调控先验网络的稳定重建能力，而不将其直接等同于真实GRN的实验验证。\n"
    (OUT_DIR / 'RESULT_3_3_3_4_WRITING_TEMPLATE.txt').write_text(text, encoding='utf-8')

def copy_summary(pred_path, label_col, score_col, fold_col, orig, shuf_summ, boot_summ, fold_df, modal_path):
    lines = []
    lines.append('=' * 90)
    lines.append('COPY THIS SUMMARY TO CHATGPT FOR RESULT RELIABILITY JUDGEMENT')
    lines.append('=' * 90)
    lines.append(f'Prediction file: {pred_path}')
    lines.append(f'Label column: {label_col}')
    lines.append(f'Score column: {score_col}')
    lines.append(f'Fold column: {fold_col}')
    lines.append(f'Modality summary file: {modal_path}')
    lines.append('')
    lines.append('[Original prediction metrics]')
    for k in ['auroc', 'aupr', 'accuracy', 'precision', 'recall', 'f1', 'mcc', 'threshold']:
        if k in orig:
            lines.append(f'{k}: {orig[k]:.6f}')
    lines.append('')
    lines.append('[Distribution separation]')
    for k in ['n_pos', 'n_neg', 'mean_pos', 'mean_neg', 'mean_separation', 'cohen_d', 'ks_statistic', 'ks_pvalue', 'hist_overlap_coefficient']:
        if k in orig:
            v = orig[k]
            lines.append(f'{k}: {v:.6g}' if isinstance(v, float) else f'{k}: {v}')
    lines.append('')
    lines.append('[Information-theoretic distribution]')
    for k in ['kl_score_to_uniform', 'kl_pos_to_uniform', 'kl_neg_to_uniform', 'entropy_score_bits', 'information_gain_bits', 'fraction_score_ge_0.8', 'fraction_score_le_0.2', 'fraction_uncertain_0.4_0.6']:
        if k in orig:
            lines.append(f'{k}: {orig[k]:.6f}')
    lines.append('')
    lines.append('[Shuffle-label null test]')
    lines.append(shuf_summ.round(6).to_string(index=False) if not shuf_summ.empty else 'No shuffle summary.')
    lines.append('')
    lines.append('[Bootstrap 95% CI]')
    if not boot_summ.empty:
        sub = boot_summ[boot_summ['metric'].isin(['auroc', 'aupr', 'accuracy', 'f1', 'mean_separation', 'ks_statistic'])]
        lines.append(sub.round(6).to_string(index=False))
    lines.append('')
    lines.append('[Cross-fold stability]')
    if not fold_df.empty:
        cols = [c for c in ['fold', 'n', 'n_pos', 'n_neg', 'auroc', 'aupr', 'accuracy', 'f1', 'mean_separation'] if c in fold_df.columns]
        lines.append(fold_df[cols].round(6).to_string(index=False))
    else:
        lines.append('No fold summary.')
    lines.append('')
    lines.append('[Important caveat]')
    lines.append('Shuffle test here shuffles labels against fixed prediction scores. It is not retraining with shuffled labels.')
    lines.append("Because labels come from ATAC motif-promoter prior and negatives are unlabeled edges, use 'candidate regulatory edges'.")
    out = OUT_DIR / 'COPY_TO_CHATGPT_rigour_summary.txt'
    out.write_text('\n'.join(lines), encoding='utf-8')
    print('\n'.join(lines))

def main():
    print('=' * 100)
    print('Model rigorous validation without external deep-learning baselines')
    print('=' * 100)
    print('OUT_DIR:', OUT_DIR)
    df, pred_path, label_col, score_col, fold_col = load_predictions()
    orig = original_summary(df)
    shuf, shuf_summ = shuffle_test(df)
    boot_res, boot_summ = bootstrap(df)
    fold_df = per_fold(df, fold_col)
    plot_score_kl(df)
    plot_roc_pr_shuffle(df)
    plot_shuffle_null(shuf, orig)
    plot_bootstrap_ci(boot_summ)
    plot_per_fold(fold_df)
    modal_path = plot_modal_contribution()
    write_templates()
    copy_summary(pred_path, label_col, score_col, fold_col, orig, shuf_summ, boot_summ, fold_df, modal_path)
    readme = f"\nOutputs:\n{OUT_DIR}\n\nFigures:\n{FIG_DIR}\n\nTables:\n{TABLE_DIR}\n\nCopy back to ChatGPT:\n{OUT_DIR / 'COPY_TO_CHATGPT_rigour_summary.txt'}\n\nRecommended figure candidates:\nFig3_3A_score_distribution_and_KL.png\nFig3_3B_ROC_PR_original_vs_shuffled.png\nFig3_3C_shuffle_label_null_test.png\nFig3_3D_bootstrap_confidence_intervals.png\nFig3_4A_modality_contribution_RNA_ATAC_Multi.png\n"
    (OUT_DIR / 'README_rigour_validation_outputs.txt').write_text(readme, encoding='utf-8')
    print('DONE')
    print(OUT_DIR / 'COPY_TO_CHATGPT_rigour_summary.txt')
'''

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    p=load_source('performance',SOURCE_performance)
    v=load_source('validation',SOURCE_validation)
    import pandas as pd
    import matplotlib.pyplot as plt
    summary=p.parse_summary(p.SUMMARY_PATH)
    test=pd.read_csv(p.TEST_PRED_PATH);all_pred=pd.read_csv(p.ALL_PRED_PATH)
    sets=p.high_confidence_sets(all_pred)
    def draw_a():
        fig,ax=plt.subplots(figsize=(7,5));p.panel_a_bar(ax,summary);save_panel(fig,'A')
    def save_panel(fig,label):
        path=OUT/(label+'.png');fig.savefig(path,dpi=300,bbox_inches='tight');fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig);panel_file(label,path)
    def draw_b():
        fig,axes=plt.subplots(1,2,figsize=(12,5));p.panel_b_roc(axes[0],test);p.panel_b_pr(axes[1],test);save_panel(fig,'B')
    def draw_c():
        fig,axes=plt.subplots(1,3,figsize=(16,5))
        for ax,mode in zip(axes,['rna','atac','multi']):p.panel_c_density(ax,test,mode,star=mode=='multi')
        save_panel(fig,'C')
    def draw_d():
        fig,axes=plt.subplots(3,1,figsize=(7,15))
        p.draw_venn2(axes[0],sets.get('multi',set()),sets.get('atac',set()),'ATAC+scRNA feature','ATAC feature')
        p.draw_venn2(axes[1],sets.get('multi',set()),sets.get('rna',set()),'ATAC+scRNA feature','scRNA feature')
        p.draw_venn3(axes[2],sets.get('multi',set()),sets.get('atac',set()),sets.get('rna',set()),labels=('ATAC+scRNA feature','ATAC feature','scRNA feature'));save_panel(fig,'D')
    for label,fn in [('A',draw_a),('B',draw_b),('C',draw_c),('D',draw_d)]:attempt(label,fn)
    def validation():
        df,*_=v.load_predictions()
        v.original_summary(df);v.shuffle_test(df)
        v.plot_score_kl(df);panel_file('E',v.FIG_DIR/'Fig3_3A_score_distribution_and_KL.png')
        v.plot_roc_pr_shuffle(df);panel_file('F',v.FIG_DIR/'Fig3_3B_ROC_PR_original_vs_shuffled.png')
        _,boot=v.bootstrap(df);v.plot_bootstrap_ci(boot);panel_file('G',v.FIG_DIR/'Fig3_3D_bootstrap_confidence_intervals.png')
    attempt('E–G',validation)

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        missing("运行", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
    finally:
        finish(7)
