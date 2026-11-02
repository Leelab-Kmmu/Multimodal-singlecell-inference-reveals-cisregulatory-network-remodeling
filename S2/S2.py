# -*- coding: utf-8 -*-
"""论文S2：同图子图代码合并。独立文件，不读取原脚本。缺失处不生成替代图。"""
from pathlib import Path
import os, sys, types, traceback, json
# 只需修改此处，或设置环境变量 CRC_DATA_DIR / CRC_OUTPUT_DIR。
DATA_DIR = Path(os.environ.get("CRC_DATA_DIR", r"E:\毕业设计\bishedata"))
OUTPUT_ROOT = Path(os.environ.get("CRC_OUTPUT_DIR", r"C:\Users\WangXinyue\Desktop\Wang_CRC_ScientificReports_Clean_Revised_Package (1)\26.10.14\111"))
FIGURE_NAME = 'S2'
OUT = OUTPUT_ROOT / FIGURE_NAME
ERRORS = []
PANELS = []

def load_source(name, source):
    # 各来源脚本使用独立命名空间，防止同名函数、路径变量互相覆盖。
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
    # 将输出重定向到当前论文图号，候选输入路径保留原有位置。
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
        # 子图保留原始矢量PDF；汇总版为PNG栅格拼版，不声称与Word手工排版完全相同。
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
    v=load_source('validation',SOURCE_validation)
    df,_,_,_,fold_col=v.load_predictions()
    if fold_col is None: raise ValueError('预测表缺fold列，无法计算S2A')
    fold=v.per_fold(df,fold_col);v.plot_per_fold(fold);panel_file('A',v.FIG_DIR/'Fig3_3E_cross_fold_stability.png')
    missing('B','未定位论文共识支持箱图的原始绘图函数；不能以重新设计的统计代码冒充提取源码。')

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        missing("运行", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
    finally:
        finish(2)
