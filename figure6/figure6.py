# -*- coding: utf-8 -*-
##### 图6 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = '图6'
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
from pathlib import Path as modules_Path
import warnings as modules_warnings
import math as modules_math
import numpy as modules_np
import pandas as modules_pd
import matplotlib.pyplot as modules_plt
import networkx as modules_nx
from matplotlib.patches import FancyArrowPatch as modules_FancyArrowPatch
from matplotlib.lines import Line2D as modules_Line2D
modules_warnings.filterwarnings('ignore')
modules_BASE_DIR = modules_Path(str(DATA_DIR))
modules_ATACMOTIF_DIR = modules_BASE_DIR / '4.0result' / 'ATACmotif'
modules_OUT_ROOT = modules_ATACMOTIF_DIR / '51_results_3_4_3_5_PATZ1_novel_edges'
modules_FIG34_DIR = modules_OUT_ROOT / '3_4_PATZ1_representative_TF_modules'
modules_FIG35_DIR = modules_OUT_ROOT / '3_5_novel_edges_functional_enrichment'
modules_TABLE_DIR = modules_OUT_ROOT / 'tables'
modules_ENRICH_DIR = modules_FIG35_DIR / 'gene_sets_for_enrichment'
modules_HIGH_CONF_EDGE_CANDIDATES = [modules_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', modules_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '07_high_confidence_edges_with_stage_weights.csv', modules_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '01_high_confidence_edges_used.csv', modules_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'high_confidence_edges_multimodal.csv', modules_ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'edgeCNN_result_analysis_figures' / 'edgeCNN_high_conf_edges_used.csv']
modules_STAGE_EDGE_CANDIDATES = [modules_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', modules_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '07_high_confidence_edges_with_stage_weights.csv']
modules_CORE_SCORE_CANDIDATES = [modules_ATACMOTIF_DIR / '47_final_results_no_celltype' / 'tables' / 'TF_centrality_core_scores.csv', modules_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '05_stage_specific_TF_role_scores.csv', modules_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '04_coreTF_scores_recomputed.csv', modules_ATACMOTIF_DIR / '30_coreTF_backbone_refined' / '03_recomputed_coreTF_scores.csv']
modules_INITIAL_ADJ_CANDIDATES = [modules_ATACMOTIF_DIR / '10_TF_TF_adjacency_safe' / 'TF_TF_initial_adjacency_matrix.csv']
modules_OPTIONAL_ENRICH_CANDIDATES = [modules_FIG35_DIR / 'enrichment_results_combined.csv', modules_ENRICH_DIR / 'enrichment_results_combined.csv']
modules_STAGES = ['Normal', 'Polyp', 'CRC']
modules_FOCUS_TFS = ['PATZ1', 'KLF5', 'SP5', 'KLF4']
modules_EXTRA_TFS = ['ZNF148', 'ZNF281', 'EGR1', 'MAZ', 'KLF16', 'ZFX']
modules_ALL_INTEREST_TFS = modules_FOCUS_TFS + modules_EXTRA_TFS
modules_TOP_CORE_N = 25
modules_EGO_MAX_TARGETS = 28
modules_NOVEL_TOP_N = 35
modules_DPI = 300
modules_COLORS = {'PATZ1': '#D1495B', 'KLF5': '#F58518', 'SP5': '#6A4C93', 'KLF4': '#2A9D8F', 'ZNF148': '#4C78A8', 'ZNF281': '#72B7B2', 'other': '#4C78A8', 'novel': '#D1495B', 'initial': '#72B7B2'}

def modules_first_existing(paths):
    for p in paths:
        p = modules_Path(p)
        if p.exists():
            return p
    return None

def modules_norm_gene(x):
    return str(x).strip().upper()

def modules_savefig(fig, path):
    path = modules_Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=modules_DPI, bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    modules_plt.close(fig)
    print('Saved:', path)

def modules_minmax01(x):
    x = modules_pd.Series(x, dtype=float).replace([modules_np.inf, -modules_np.inf], modules_np.nan)
    if x.dropna().empty:
        return modules_pd.Series(modules_np.zeros(len(x)), index=x.index)
    x = x.fillna(x.median())
    if float(x.max()) == float(x.min()):
        return modules_pd.Series(modules_np.ones(len(x)) * 0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())

def modules_clean_edges(df):
    df = df.copy()
    if 'regulator_tf' not in df.columns:
        for c in ['source_tf', 'source', 'regulator', 'TF', 'from', 'from_tf', 'tf1', 'u']:
            if c in df.columns:
                df['regulator_tf'] = df[c]
                break
    if 'target_tf' not in df.columns:
        for c in ['target_tf', 'target', 'target_gene', 'gene', 'to', 'to_tf', 'tf2', 'v']:
            if c in df.columns:
                df['target_tf'] = df[c]
                break
    if 'regulator_tf' not in df.columns or 'target_tf' not in df.columns:
        raise ValueError('Edge table must include regulator/source and target columns.')
    df['regulator_tf'] = df['regulator_tf'].map(modules_norm_gene)
    df['target_tf'] = df['target_tf'].map(modules_norm_gene)
    if 'mean_score' not in df.columns:
        for c in ['score', 'prob', 'pred_prob', 'y_prob', 'prediction_score', 'weight', 'confidence', 'consensus_score']:
            if c in df.columns:
                df['mean_score'] = modules_pd.to_numeric(df[c], errors='coerce')
                break
    if 'mean_score' not in df.columns:
        df['mean_score'] = 1.0
    df['mean_score'] = modules_pd.to_numeric(df['mean_score'], errors='coerce').fillna(1.0)
    df = df[(df['regulator_tf'] != '') & (df['target_tf'] != '')]
    df = df[df['regulator_tf'] != df['target_tf']].copy()
    df = df.drop_duplicates(['regulator_tf', 'target_tf']).reset_index(drop=True)
    return df

def modules_infer_stage_activity_cols(df):
    if df is None:
        return {}
    cols = list(df.columns)
    mapping = {}
    for stage in modules_STAGES:
        st = stage.lower()
        exact = [f'edge_activity_{stage}', f'{stage}_edge_activity', f'activity_{stage}', f'{stage}_activity', f'edge_score_{stage}', f'{stage}_score', f'mean_score_{stage}', f'{stage}_mean_score', stage]
        for e in exact:
            if e in cols:
                mapping[stage] = e
                break
        if stage in mapping:
            continue
        candidates = []
        for c in cols:
            cl = c.lower()
            if st in cl and 'delta' not in cl and any((k in cl for k in ['activity', 'score', 'weight', 'prob'])):
                candidates.append(c)
        if candidates:
            candidates = sorted(candidates, key=lambda x: ('activity' not in x.lower(), len(x)))
            mapping[stage] = candidates[0]
    return mapping

def modules_load_high_conf_edges():
    p = modules_first_existing(modules_HIGH_CONF_EDGE_CANDIDATES)
    if p is None:
        raise FileNotFoundError('Cannot find high-confidence edge file.')
    edges = modules_clean_edges(modules_pd.read_csv(p))
    print('Using high-confidence edge file:', p)
    print('Edges:', edges.shape)
    return (edges, p)

def modules_load_stage_edges():
    p = modules_first_existing(modules_STAGE_EDGE_CANDIDATES)
    if p is None:
        print('No stage edge table found.')
        return (None, None)
    edges = modules_clean_edges(modules_pd.read_csv(p))
    print('Using stage edge file:', p)
    return (edges, p)

def modules_build_graph(edges):
    G = modules_nx.DiGraph()
    for r in edges.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, weight=float(r.mean_score))
    return G

def modules_load_core_score_table():
    p = modules_first_existing(modules_CORE_SCORE_CANDIDATES)
    if p is None:
        return (None, None)
    df = modules_pd.read_csv(p, index_col=0)
    df.index = df.index.astype(str).str.upper().str.strip()
    for c in df.columns:
        df[c] = modules_pd.to_numeric(df[c], errors='coerce')
    print('Using core score table:', p)
    return (df, p)

def modules_compute_centrality(edges, core_table=None):
    G = modules_build_graph(edges)
    nodes = sorted(G.nodes())
    out_degree = dict(G.out_degree())
    in_degree = dict(G.in_degree())
    weighted_out = {n: 0.0 for n in nodes}
    weighted_in = {n: 0.0 for n in nodes}
    for u, v, d in G.edges(data=True):
        w = float(d.get('weight', 1.0))
        weighted_out[u] += w
        weighted_in[v] += w
    try:
        pagerank = modules_nx.pagerank(G, weight='weight')
    except Exception:
        pagerank = {n: 0.0 for n in nodes}
    cent = modules_pd.DataFrame(index=nodes)
    cent['out_degree'] = cent.index.map(out_degree).fillna(0).astype(float)
    cent['in_degree'] = cent.index.map(in_degree).fillna(0).astype(float)
    cent['weighted_out_degree'] = cent.index.map(weighted_out).fillna(0).astype(float)
    cent['weighted_in_degree'] = cent.index.map(weighted_in).fillna(0).astype(float)
    cent['pagerank'] = cent.index.map(pagerank).fillna(0).astype(float)
    if core_table is not None:
        for c in core_table.columns:
            if c not in cent.columns:
                cent[c] = cent.index.map(core_table[c]).fillna(modules_np.nan)
    fallback = 0.35 * modules_minmax01(cent['out_degree']) + 0.35 * modules_minmax01(cent['weighted_out_degree']) + 0.2 * modules_minmax01(cent['pagerank']) + 0.1 * modules_minmax01(cent['in_degree'])
    if 'core_score' not in cent.columns or cent['core_score'].isna().all():
        cent['core_score'] = fallback
    else:
        cent['core_score'] = modules_pd.to_numeric(cent['core_score'], errors='coerce')
        cent['core_score'] = cent['core_score'].fillna(fallback)
    cent = cent.sort_values('core_score', ascending=False)
    cent.to_csv(modules_TABLE_DIR / 'TF_centrality_core_scores.csv', encoding='utf-8-sig')
    return cent

def modules_compute_tf_stage_activity(stage_edges):
    if stage_edges is None:
        return modules_pd.DataFrame()
    stage_cols = modules_infer_stage_activity_cols(stage_edges)
    if len(stage_cols) != 3:
        print('Cannot infer three stage activity columns. Stage-trajectory panels may be skipped.')
        return modules_pd.DataFrame()
    df = stage_edges.copy()
    for st, c in stage_cols.items():
        df[c] = modules_pd.to_numeric(df[c], errors='coerce').fillna(0)
    mat = df.groupby('regulator_tf')[[stage_cols[s] for s in modules_STAGES]].sum()
    mat.columns = modules_STAGES
    mat = mat.reindex(sorted(set(df['regulator_tf']).union(df['target_tf']))).fillna(0)
    mat.to_csv(modules_TABLE_DIR / 'TF_stage_out_activity_matrix.csv', encoding='utf-8-sig')
    return mat

def modules_get_module_edges(edges, tf, max_targets=modules_EGO_MAX_TARGETS):
    sub = edges[edges['regulator_tf'] == tf].sort_values('mean_score', ascending=False).head(max_targets).copy()
    sub['module_tf'] = tf
    return sub

def modules_module_overlap_matrix(edges, tfs):
    targets = {tf: set(edges.loc[edges['regulator_tf'] == tf, 'target_tf']) for tf in tfs}
    mat = modules_pd.DataFrame(index=tfs, columns=tfs, dtype=float)
    jac = modules_pd.DataFrame(index=tfs, columns=tfs, dtype=float)
    for a in tfs:
        for b in tfs:
            inter = len(targets[a] & targets[b])
            union = len(targets[a] | targets[b])
            mat.loc[a, b] = inter
            jac.loc[a, b] = inter / union if union else 0
    mat.to_csv(modules_TABLE_DIR / 'module_target_overlap_counts.csv', encoding='utf-8-sig')
    jac.to_csv(modules_TABLE_DIR / 'module_target_overlap_jaccard.csv', encoding='utf-8-sig')
    return (mat, jac)

##### 图6B #####
def modules_plot_tf_stage_trajectories(tf_stage, out_path, tfs=modules_FOCUS_TFS):
    fig, ax = modules_plt.subplots(figsize=(6.5, 4.8))
    plotted = False
    x = modules_np.arange(len(modules_STAGES))
    for tf in tfs:
        if tf not in tf_stage.index:
            continue
        vals = tf_stage.loc[tf, modules_STAGES].astype(float)
        vals = modules_minmax01(vals).values
        ax.plot(x, vals, marker='o', linewidth=2.5, label=tf, color=modules_COLORS.get(tf, None))
        plotted = True
    if not plotted:
        ax.text(0.5, 0.5, 'No stage activity available', ha='center', va='center')
    ax.set_xticks(x)
    ax.set_xticklabels(modules_STAGES)
    ax.set_ylabel('Normalized out-activity')
    ax.set_title('Stage trajectories of PATZ1 and representative TFs', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=9)
    modules_savefig(fig, out_path)

def modules_draw_ego(ax, edges, cent, tf, max_targets=22):
    sub = modules_get_module_edges(edges, tf, max_targets=max_targets)
    if sub.empty:
        ax.text(0.5, 0.5, f'No outgoing module for {tf}', ha='center', va='center')
        ax.axis('off')
        return
    G = modules_nx.DiGraph()
    for r in sub.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, weight=float(r.mean_score))
    targets = [r.target_tf for r in sub.itertuples(index=False)]
    pos = {tf: modules_np.array([0.0, 0.0])}
    angles = modules_np.linspace(0, 2 * modules_np.pi, len(targets), endpoint=False)
    for t, a in zip(targets, angles):
        pos[t] = modules_np.array([modules_np.cos(a), modules_np.sin(a)])
    weights = sub['mean_score'].astype(float).values
    widths = 0.9 + 3.3 * modules_minmax01(modules_pd.Series(weights)).values
    for idx, r in enumerate(sub.itertuples(index=False)):
        patch = modules_FancyArrowPatch(pos[r.regulator_tf], pos[r.target_tf], arrowstyle='-|>', mutation_scale=11, linewidth=float(widths[idx]), color=modules_COLORS.get(tf, '#4C78A8'), alpha=0.48, shrinkA=10, shrinkB=7, connectionstyle='arc3,rad=0.05')
        ax.add_patch(patch)
    sizes = []
    colors = []
    for n in G.nodes():
        if n == tf:
            sizes.append(800)
            colors.append(modules_COLORS.get(tf, '#D1495B'))
        else:
            cs = float(cent.loc[n, 'core_score']) if n in cent.index else 0.1
            sizes.append(90 + 420 * cs)
            colors.append('#7AA6C2')
    ax.scatter([pos[n][0] for n in G.nodes()], [pos[n][1] for n in G.nodes()], s=sizes, color=colors, edgecolor='white', linewidth=0.8, zorder=3)
    for n in G.nodes():
        ax.text(pos[n][0], pos[n][1], n, fontsize=7, ha='center', va='center', fontweight='bold' if n == tf else 'normal', zorder=4)
    ax.set_title(f'{tf}-centered module', fontsize=11, fontweight='bold')
    ax.axis('off')
    ax.set_aspect('equal')

##### 图6C #####
def modules_plot_representative_ego_grid(edges, cent, out_path):
    tfs = modules_FOCUS_TFS
    fig, axes = modules_plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.ravel()
    for ax, tf in zip(axes, tfs):
        modules_draw_ego(ax, edges, cent, tf, max_targets=18)
    fig.suptitle('PATZ1 and representative TF-centered regulatory modules', fontsize=15, fontweight='bold')
    modules_savefig(fig, out_path)

##### 图6A #####
def modules_plot_module_overlap(jaccard, out_path):
    if jaccard.empty:
        return
    fig, ax = modules_plt.subplots(figsize=(5.6, 4.8))
    im = ax.imshow(jaccard.values.astype(float), cmap='YlGnBu', vmin=0, vmax=max(0.01, float(jaccard.values.max())))
    ax.set_xticks(range(jaccard.shape[1]))
    ax.set_xticklabels(jaccard.columns, rotation=30, ha='right')
    ax.set_yticks(range(jaccard.shape[0]))
    ax.set_yticklabels(jaccard.index)
    for i in range(jaccard.shape[0]):
        for j in range(jaccard.shape[1]):
            val = float(jaccard.iloc[i, j])
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', color='white' if val > 0.35 else 'black', fontsize=9)
    ax.set_title('Target-overlap between TF modules', fontsize=13, fontweight='bold')
    cbar = modules_plt.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label('Jaccard index', fontsize=9)
    modules_savefig(fig, out_path)
modules_TABLE_DIR = OUT / 'modules' / 'TABLE_DIR'
modules_TABLE_DIR.mkdir(parents=True, exist_ok=True)
modules_ENRICH_DIR = OUT / 'modules' / 'ENRICH_DIR'
modules_ENRICH_DIR.mkdir(parents=True, exist_ok=True)
modules_FIG34_DIR = OUT / 'modules' / 'FIG34_DIR'
modules_FIG34_DIR.mkdir(parents=True, exist_ok=True)
modules_FIG35_DIR = OUT / 'modules' / 'FIG35_DIR'
modules_FIG35_DIR.mkdir(parents=True, exist_ok=True)
modules_OUT_ROOT = OUT / 'modules' / 'OUT_ROOT'
modules_OUT_ROOT.mkdir(parents=True, exist_ok=True)
from pathlib import Path as network_Path
import warnings as network_warnings
import math as network_math
import numpy as network_np
import pandas as network_pd
import matplotlib.pyplot as network_plt
import networkx as network_nx
from matplotlib.patches import FancyArrowPatch as network_FancyArrowPatch
from matplotlib.lines import Line2D as network_Line2D
network_warnings.filterwarnings('ignore')
network_BASE_DIR = network_Path(str(DATA_DIR))
network_ATACMOTIF_DIR = network_BASE_DIR / '4.0result' / 'ATACmotif'
network_OUT_ROOT = network_ATACMOTIF_DIR / '47_final_results_no_celltype'
network_FIG2_DIR = network_OUT_ROOT / 'Figure2_high_conf_network_coreTF'
network_FIG3_DIR = network_OUT_ROOT / 'Figure3_stage_network_remodeling'
network_FIG4_DIR = network_OUT_ROOT / 'Figure4_PATZ1_deep_analysis'
network_FIG5_DIR = network_OUT_ROOT / 'Figure5_representative_TF_enrichment'
network_TABLE_DIR = network_OUT_ROOT / 'tables'
network_ENRICH_DIR = network_FIG5_DIR / 'gene_sets_for_enrichment'
network_HIGH_CONF_EDGE_CANDIDATES = [network_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', network_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '07_high_confidence_edges_with_stage_weights.csv', network_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '01_high_confidence_edges_used.csv', network_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'high_confidence_edges_multimodal.csv', network_ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'edgeCNN_result_analysis_figures' / 'edgeCNN_high_conf_edges_used.csv']
network_STAGE_EDGE_CANDIDATES = [network_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', network_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '07_high_confidence_edges_with_stage_weights.csv']
network_CORE_SCORE_CANDIDATES = [network_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '05_stage_specific_TF_role_scores.csv', network_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '04_coreTF_scores_recomputed.csv', network_ATACMOTIF_DIR / '30_coreTF_backbone_refined' / '03_recomputed_coreTF_scores.csv']
network_INITIAL_ADJ_CANDIDATES = [network_ATACMOTIF_DIR / '10_TF_TF_adjacency_safe' / 'TF_TF_initial_adjacency_matrix.csv']
network_PREDICTION_CANDIDATES = [network_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'predictions_all_folds_all_label_edges.csv', network_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'all_predictions_multimodal.csv']
network_OPTIONAL_ACTIVITY_CANDIDATES = [network_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / 'TF_stage_multimodal_activity.csv', network_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '03_TF_stage_activity_scores.csv', network_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '03_TF_stage_activity_scores.csv', network_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / 'TF_stage_multimodal_activity.csv']
network_STAGES = ['Normal', 'Polyp', 'CRC']
network_FOCUS_TFS = ['PATZ1', 'KLF5', 'SP5', 'KLF4', 'ZNF148', 'ZNF281', 'EGR1', 'MAZ']
network_REPRESENTATIVE_TFS = ['KLF5', 'SP5', 'KLF4']
network_DPI = 300
network_COLOR_STAGE = {'Normal': '#4C78A8', 'Polyp': '#F58518', 'CRC': '#E45756'}
network_COLOR_ROLE = {'stable_backbone_score': '#6A4C93', 'early_transition_score': '#F2C14E', 'late_transition_score': '#E45756', 'crc_enhanced_score': '#D1495B', 'crc_suppressed_score': '#2A9D8F'}
network_ROLE_LABELS = {'stable_backbone_score': 'Stable\nbackbone', 'early_transition_score': 'Early\ntransition', 'late_transition_score': 'Late\ntransition', 'crc_enhanced_score': 'CRC\nenhanced', 'crc_suppressed_score': 'CRC\nsuppressed'}

def network_first_existing(paths):
    for p in paths:
        if network_Path(p).exists():
            return network_Path(p)
    return None

def network_savefig(fig, path):
    path = network_Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=network_DPI, bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    network_plt.close(fig)
    print('Saved:', path)

##### 图6D #####
def network_fig5D_enrichment_dotplot_or_placeholder(out_path):
    candidates = [network_FIG5_DIR / 'enrichment_results_combined.csv', network_ENRICH_DIR / 'enrichment_results_combined.csv']
    p = network_first_existing(candidates)
    if p is None:
        summary_path = network_ENRICH_DIR / 'gene_sets_summary.csv'
        if summary_path.exists():
            summary = network_pd.read_csv(summary_path)
        else:
            summary = network_pd.DataFrame({'gene_set': [], 'n_genes': []})
        fig, ax = network_plt.subplots(figsize=(7.2, 4.8))
        ax.axis('off')
        ax.text(0.02, 0.88, 'GO/KEGG enrichment input prepared', fontsize=14, fontweight='bold', transform=ax.transAxes)
        ax.text(0.02, 0.75, 'Run the provided R script to generate the final enrichment dotplot:', fontsize=10, transform=ax.transAxes)
        ax.text(0.02, 0.65, 'run_enrichment_clusterProfiler.R', fontsize=11, fontweight='bold', color='#D1495B', transform=ax.transAxes)
        if not summary.empty:
            y = 0.52
            ax.text(0.02, y, 'Exported gene sets:', fontsize=10, fontweight='bold', transform=ax.transAxes)
            y -= 0.07
            for _, r in summary.head(8).iterrows():
                ax.text(0.05, y, f"{r['gene_set']}: {int(r['n_genes'])} genes", fontsize=9, transform=ax.transAxes)
                y -= 0.055
        network_savefig(fig, out_path)
        return
    enrich = network_pd.read_csv(p)
    desc_col = 'Description' if 'Description' in enrich.columns else 'term'
    padj_col = 'p.adjust' if 'p.adjust' in enrich.columns else 'padj'
    count_col = 'Count' if 'Count' in enrich.columns else 'count'
    set_col = 'gene_set' if 'gene_set' in enrich.columns else 'GeneSet'
    source_col = 'source' if 'source' in enrich.columns else 'ontology'
    enrich[padj_col] = network_pd.to_numeric(enrich[padj_col], errors='coerce')
    enrich[count_col] = network_pd.to_numeric(enrich[count_col], errors='coerce').fillna(1)
    enrich = enrich.dropna(subset=[padj_col])
    enrich = enrich[enrich[padj_col] <= 0.1].copy()
    if enrich.empty:
        print('No enrichment terms with adjusted p <= 0.10.')
        return
    enrich = enrich.sort_values([set_col, padj_col]).groupby(set_col).head(5).copy()
    enrich['label'] = enrich[set_col].astype(str) + ' | ' + enrich[desc_col].astype(str)
    enrich = enrich.sort_values(padj_col, ascending=False).tail(30)
    enrich['neglog10_padj'] = -network_np.log10(enrich[padj_col].clip(lower=1e-300))
    fig, ax = network_plt.subplots(figsize=(9.5, max(5.5, 0.24 * enrich.shape[0])))
    y = network_np.arange(enrich.shape[0])
    sc = ax.scatter(enrich['neglog10_padj'], y, s=30 + 18 * enrich[count_col], c=enrich['neglog10_padj'], cmap='YlOrRd', edgecolor='white', linewidth=0.5, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(enrich['label'], fontsize=7)
    ax.set_xlabel('-log10(adjusted P value)')
    ax.set_title('Functional enrichment of model-derived TF sets', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    cbar = network_plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('-log10(adj. P)', fontsize=9)
    network_savefig(fig, out_path)
network_TABLE_DIR = OUT / 'network' / 'TABLE_DIR'
network_TABLE_DIR.mkdir(parents=True, exist_ok=True)
network_ENRICH_DIR = OUT / 'network' / 'ENRICH_DIR'
network_ENRICH_DIR.mkdir(parents=True, exist_ok=True)
network_FIG3_DIR = OUT / 'network' / 'FIG3_DIR'
network_FIG3_DIR.mkdir(parents=True, exist_ok=True)
network_FIG2_DIR = OUT / 'network' / 'FIG2_DIR'
network_FIG2_DIR.mkdir(parents=True, exist_ok=True)
network_OUT_ROOT = OUT / 'network' / 'OUT_ROOT'
network_OUT_ROOT.mkdir(parents=True, exist_ok=True)
network_FIG4_DIR = OUT / 'network' / 'FIG4_DIR'
network_FIG4_DIR.mkdir(parents=True, exist_ok=True)
network_FIG5_DIR = OUT / 'network' / 'FIG5_DIR'
network_FIG5_DIR.mkdir(parents=True, exist_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    edges, _ = modules_load_high_conf_edges()
    stage, _ = modules_load_stage_edges()
    core, _ = modules_load_core_score_table()
    cent = modules_compute_centrality(edges, core)
    activity = modules_compute_tf_stage_activity(stage)
    _, overlap = modules_module_overlap_matrix(edges, modules_FOCUS_TFS)
    for label, fn in [('A', lambda path: modules_plot_module_overlap(overlap, path)), ('B', lambda path: modules_plot_tf_stage_trajectories(activity, path)), ('C', lambda path: modules_plot_representative_ego_grid(edges, cent, path))]:
        path = OUT / (label + '.png')
        attempt(label, lambda label=label, fn=fn, path=path: (fn(path), panel_file(label, path)))

    def draw_d():
        candidates = [OUT / 'enrichment_results_combined.csv', DATA_DIR / '4.0result' / 'ATACmotif' / '47_final_results_no_celltype' / 'Figure5_representative_TF_enrichment' / 'enrichment_results_combined.csv', DATA_DIR / '4.0result' / 'ATACmotif' / '47_final_results_no_celltype' / 'Figure5_representative_TF_enrichment' / 'gene_sets_for_enrichment' / 'enrichment_results_combined.csv']
        import shutil
        found = next((x for x in candidates if x.exists()), None)
        if found is None:
            raise FileNotFoundError('缺真实enrichment_results_combined.csv及原R脚本，无法生成D；可将已有CSV放入本图输出目录')
        shutil.copy2(found, network_FIG5_DIR / 'enrichment_results_combined.csv')
        path = OUT / 'D.png'
        network_fig5D_enrichment_dotplot_or_placeholder(path)
        panel_file('D', path)
    attempt('D', draw_d)
if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        missing('运行', f'{type(exc).__name__}: {exc}')
        traceback.print_exc()
    finally:
        finish(4)

if __name__ == "__main__" and ERRORS:
    raise SystemExit(1)
