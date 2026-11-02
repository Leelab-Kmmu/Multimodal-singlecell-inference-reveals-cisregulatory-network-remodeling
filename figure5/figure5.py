# -*- coding: utf-8 -*-
##### 图5 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = '图5'
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

def network_norm_gene(x):
    return str(x).strip().upper()

def network_savefig(fig, path):
    path = network_Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=network_DPI, bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    network_plt.close(fig)
    print('Saved:', path)

def network_minmax01(x):
    x = network_pd.Series(x, dtype=float).replace([network_np.inf, -network_np.inf], network_np.nan)
    if x.dropna().empty:
        return network_pd.Series(network_np.zeros(len(x)), index=x.index)
    x = x.fillna(x.median())
    if float(x.max()) == float(x.min()):
        return network_pd.Series(network_np.ones(len(x)) * 0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())

def network_zscore_rows(df):
    out = df.astype(float).copy()
    for idx in out.index:
        x = out.loc[idx].astype(float)
        sd = x.std()
        out.loc[idx] = 0 if sd == 0 or network_pd.isna(sd) else (x - x.mean()) / sd
    return out.replace([network_np.inf, -network_np.inf], 0).fillna(0)

def network_clean_edges(df):
    df = df.copy()
    if 'regulator_tf' not in df.columns:
        for c in ['source_tf', 'source', 'regulator', 'TF', 'from', 'from_tf', 'tf1']:
            if c in df.columns:
                df['regulator_tf'] = df[c]
                break
    if 'target_tf' not in df.columns:
        for c in ['target', 'target_gene', 'gene', 'to', 'to_tf', 'tf2']:
            if c in df.columns:
                df['target_tf'] = df[c]
                break
    if 'regulator_tf' not in df.columns or 'target_tf' not in df.columns:
        raise ValueError('Edge table must contain regulator_tf/source and target_tf/target columns.')
    df['regulator_tf'] = df['regulator_tf'].map(network_norm_gene)
    df['target_tf'] = df['target_tf'].map(network_norm_gene)
    if 'mean_score' not in df.columns:
        for c in ['score', 'prob', 'y_prob', 'pred_prob', 'prediction_score', 'weight', 'consensus_score']:
            if c in df.columns:
                df['mean_score'] = network_pd.to_numeric(df[c], errors='coerce')
                break
    if 'mean_score' not in df.columns:
        df['mean_score'] = 1.0
    df['mean_score'] = network_pd.to_numeric(df['mean_score'], errors='coerce').fillna(1.0)
    df = df[df['regulator_tf'] != df['target_tf']].copy()
    df = df.drop_duplicates(['regulator_tf', 'target_tf']).reset_index(drop=True)
    return df

def network_load_high_conf_edges():
    p = network_first_existing(network_HIGH_CONF_EDGE_CANDIDATES)
    if p is None:
        raise FileNotFoundError('Cannot find high confidence edge file.')
    edges = network_clean_edges(network_pd.read_csv(p))
    edges.to_csv(network_TABLE_DIR / 'high_confidence_edges_used.csv', index=False, encoding='utf-8-sig')
    print('Using high-confidence edge file:', p)
    print('High-confidence edges:', edges.shape)
    return (edges, p)

def network_load_stage_edges():
    p = network_first_existing(network_STAGE_EDGE_CANDIDATES)
    if p is None:
        print('No stage edge activity file found. Stage-specific plots may be limited.')
        return (None, None)
    edges = network_clean_edges(network_pd.read_csv(p))
    print('Using stage edge file:', p)
    return (edges, p)

def network_build_graph(edges):
    G = network_nx.DiGraph()
    for r in edges.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, weight=float(r.mean_score))
    return G

def network_compute_centrality(edges, core_score_table=None):
    G = network_build_graph(edges)
    nodes = sorted(G.nodes())
    out_deg = dict(G.out_degree())
    in_deg = dict(G.in_degree())
    weighted_out = {n: 0.0 for n in nodes}
    weighted_in = {n: 0.0 for n in nodes}
    for u, v, d in G.edges(data=True):
        w = float(d.get('weight', 1.0))
        weighted_out[u] += w
        weighted_in[v] += w
    try:
        pagerank = network_nx.pagerank(G, weight='weight')
    except Exception:
        pagerank = {n: 0.0 for n in nodes}
    cent = network_pd.DataFrame(index=nodes)
    cent['out_degree'] = cent.index.map(out_deg).fillna(0).astype(float)
    cent['in_degree'] = cent.index.map(in_deg).fillna(0).astype(float)
    cent['weighted_out_degree'] = cent.index.map(weighted_out).fillna(0).astype(float)
    cent['weighted_in_degree'] = cent.index.map(weighted_in).fillna(0).astype(float)
    cent['pagerank'] = cent.index.map(pagerank).fillna(0).astype(float)
    if core_score_table is not None:
        for c in core_score_table.columns:
            if c not in cent.columns:
                cent[c] = cent.index.map(core_score_table[c]).fillna(network_np.nan)
    if 'core_score' not in cent.columns or cent['core_score'].isna().all():
        cent['core_score'] = 0.35 * network_minmax01(cent['out_degree']) + 0.35 * network_minmax01(cent['weighted_out_degree']) + 0.2 * network_minmax01(cent['pagerank']) + 0.1 * network_minmax01(cent['in_degree'])
    else:
        cent['core_score'] = network_pd.to_numeric(cent['core_score'], errors='coerce')
        missing = cent['core_score'].isna()
        if missing.any():
            fallback = 0.35 * network_minmax01(cent['out_degree']) + 0.35 * network_minmax01(cent['weighted_out_degree']) + 0.2 * network_minmax01(cent['pagerank']) + 0.1 * network_minmax01(cent['in_degree'])
            cent.loc[missing, 'core_score'] = fallback.loc[missing]
    cent = cent.sort_values('core_score', ascending=False)
    cent.to_csv(network_TABLE_DIR / 'TF_centrality_core_scores.csv', encoding='utf-8-sig')
    return cent

def network_load_core_score_table():
    p = network_first_existing(network_CORE_SCORE_CANDIDATES)
    if p is None:
        print('Core score file not found. Core scores will be computed from network.')
        return (None, None)
    df = network_pd.read_csv(p, index_col=0)
    df.index = df.index.astype(str).str.upper().str.strip()
    for c in df.columns:
        df[c] = network_pd.to_numeric(df[c], errors='ignore')
    print('Using core score file:', p)
    return (df, p)

def network_infer_stage_activity_cols(df):
    if df is None:
        return {}
    cols = list(df.columns)
    mapping = {}
    for stage in network_STAGES:
        st = stage.lower()
        exacts = [f'edge_activity_{stage}', f'{stage}_edge_activity', f'activity_{stage}', f'{stage}_activity', f'edge_score_{stage}', f'{stage}_score', stage]
        for e in exacts:
            if e in cols:
                mapping[stage] = e
                break
        if stage in mapping:
            continue
        candidates = []
        for c in cols:
            cl = c.lower()
            if st in cl and 'delta' not in cl and any((k in cl for k in ['activity', 'weight', 'score', 'out', 'reg'])):
                candidates.append(c)
        if candidates:
            candidates = sorted(candidates, key=lambda x: ('edge' not in x.lower(), 'activity' not in x.lower(), len(x)))
            mapping[stage] = candidates[0]
    return mapping

def network_compute_stage_tf_activity(stage_edges):
    if stage_edges is None:
        return network_pd.DataFrame()
    stage_cols = network_infer_stage_activity_cols(stage_edges)
    if len(stage_cols) == 3:
        tmp = stage_edges.copy()
        for st, col in stage_cols.items():
            tmp[col] = network_pd.to_numeric(tmp[col], errors='coerce').fillna(0)
        mat = tmp.groupby('regulator_tf')[[stage_cols[s] for s in network_STAGES]].sum()
        mat.columns = network_STAGES
        mat = mat.reindex(sorted(set(stage_edges['regulator_tf']).union(stage_edges['target_tf']))).fillna(0)
        mat.to_csv(network_TABLE_DIR / 'TF_stage_out_activity_matrix.csv', encoding='utf-8-sig')
        return mat
    print('Could not infer explicit stage activity columns from stage edge table.')
    return network_pd.DataFrame()

##### 图5C #####
def network_fig3C_early_late_scatter(tf_stage, core_scores, cent, out_path):
    if tf_stage.empty:
        print('Skip Fig3C: no TF stage activity matrix.')
        return
    df = tf_stage.copy()
    for s in network_STAGES:
        df[s] = network_pd.to_numeric(df[s], errors='coerce').fillna(0)
    df['early_delta_Polyp_vs_Normal'] = df['Polyp'] - df['Normal']
    df['late_delta_CRC_vs_Polyp'] = df['CRC'] - df['Polyp']
    df['core_score'] = df.index.map(cent['core_score']).fillna(0)
    df.to_csv(network_TABLE_DIR / 'TF_early_late_transition_scores.csv', encoding='utf-8-sig')
    fig, ax = network_plt.subplots(figsize=(7.2, 6.2))
    sizes = 25 + 350 * network_minmax01(df['core_score']).values
    colors = df['late_delta_CRC_vs_Polyp'].values
    sc = ax.scatter(df['early_delta_Polyp_vs_Normal'], df['late_delta_CRC_vs_Polyp'], s=sizes, c=colors, cmap='coolwarm', alpha=0.75, edgecolor='white', linewidth=0.4)
    ax.axhline(0, color='#555555', lw=0.8, ls='--')
    ax.axvline(0, color='#555555', lw=0.8, ls='--')
    ax.set_xlabel('Early transition: Polyp - Normal')
    ax.set_ylabel('Late transition: CRC - Polyp')
    ax.set_title('Early vs late TF regulatory activity transition', fontsize=14, fontweight='bold')
    label_tfs = list(dict.fromkeys(network_FOCUS_TFS + df.sort_values('core_score', ascending=False).head(10).index.tolist()))
    for tf in label_tfs:
        if tf in df.index:
            ax.text(df.loc[tf, 'early_delta_Polyp_vs_Normal'], df.loc[tf, 'late_delta_CRC_vs_Polyp'], f' {tf}', fontsize=8)
    cbar = network_plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label('Late transition delta', fontsize=10)
    ax.grid(alpha=0.2)
    network_savefig(fig, out_path)

##### 图5B #####
def network_fig3D_dynamic_edge_heatmap(stage_edges, out_path, top_n=35):
    if stage_edges is None:
        print('Skip Fig3D: no stage edge table.')
        return
    stage_cols = network_infer_stage_activity_cols(stage_edges)
    if len(stage_cols) != 3:
        print('Skip Fig3D: cannot infer stage edge activity columns.')
        return
    df = stage_edges.copy()
    for s, c in stage_cols.items():
        df[c] = network_pd.to_numeric(df[c], errors='coerce').fillna(0)
    vals = df[[stage_cols[s] for s in network_STAGES]].copy()
    vals.columns = network_STAGES
    df['dynamic_range'] = vals.max(axis=1) - vals.min(axis=1)
    sub = df.sort_values('dynamic_range', ascending=False).head(top_n).copy()
    mat = sub[[stage_cols[s] for s in network_STAGES]].copy()
    mat.columns = network_STAGES
    mat.index = sub['regulator_tf'] + '→' + sub['target_tf']
    z = network_zscore_rows(mat)
    sub.to_csv(network_TABLE_DIR / 'dynamic_edges_used_for_Fig3D.csv', index=False, encoding='utf-8-sig')
    fig, ax = network_plt.subplots(figsize=(5.8, max(6, 0.22 * top_n)))
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(network_STAGES)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index, fontsize=6.5)
    ax.set_title('Dynamic TF-TF edge activity across stages', fontsize=14, fontweight='bold')
    cbar = network_plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label('Row z-score', fontsize=10)
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
from pathlib import Path as roles_Path
import numpy as roles_np
import pandas as roles_pd
import matplotlib.pyplot as roles_plt
roles_BASE_DIR = roles_Path(str(DATA_DIR))
roles_ATACMOTIF_DIR = roles_BASE_DIR / '4.0result' / 'ATACmotif'
roles_IN_DIR = roles_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks'
roles_OUT_DIR = roles_Path('C:\\Users\\WangXinyue\\Desktop\\面试个人资料\\毕设图\\4.0\\斜体')
roles_NODE_SCORE_PATH = roles_IN_DIR / '05_stage_specific_TF_role_scores.csv'
roles_SELECTED_TF_PATH = roles_IN_DIR / '07_data_driven_selected_TFs_for_plots.csv'
roles_TOP_PER_ROLE = 6
roles_MAX_TFS = 24
roles_X_GAP = 1.75
roles_Y_GAP = 1.28
roles_MIN_SIZE = 100
roles_MAX_SIZE = 1200
roles_FIG_WIDTH = 13.5
roles_BASE_HEIGHT = 5.0
roles_HEIGHT_PER_TF = 0.42
roles_STAR_SIZE = 110
roles_LABEL_FONTSIZE = 14
roles_TITLE_FONTSIZE = 18
roles_ROLE_COLS = ['stable_backbone_score', 'early_transition_score', 'late_transition_score', 'crc_enhanced_score', 'crc_suppressed_score']
roles_ROLE_LABELS = ['Stable\nbackbone', 'Early\ntransition', 'Late\ntransition', 'CRC\nenhanced', 'CRC\nsuppressed']
roles_ROLE_NAME_MAP = {'stable_backbone_score': 'Stable_backbone', 'early_transition_score': 'Early_transition', 'late_transition_score': 'Late_transition', 'crc_enhanced_score': 'CRC_enhanced', 'crc_suppressed_score': 'CRC_suppressed'}
roles_ROLE_COLORS = {'Stable_backbone': '#6A4C93', 'Early_transition': '#F2C14E', 'Late_transition': '#E45756', 'CRC_enhanced': '#D1495B', 'CRC_suppressed': '#2A9D8F', 'Other': '#B8B8B8'}

def roles_select_tfs(node):
    tfs = []
    for col in roles_ROLE_COLS:
        if col in node.columns:
            tfs += node.sort_values(col, ascending=False).head(roles_TOP_PER_ROLE).index.tolist()
    if roles_SELECTED_TF_PATH.exists():
        selected = roles_pd.read_csv(roles_SELECTED_TF_PATH).iloc[:, 0].astype(str).str.upper().str.strip().tolist()
        tfs += selected[:12]
    tfs = list(dict.fromkeys([tf for tf in tfs if tf in node.index]))
    return tfs[:roles_MAX_TFS]

##### 图5A #####
def roles_main():
    if not roles_NODE_SCORE_PATH.exists():
        raise FileNotFoundError(f'Cannot find: {roles_NODE_SCORE_PATH}')
    node = roles_pd.read_csv(roles_NODE_SCORE_PATH, index_col=0)
    node.index = node.index.astype(str).str.upper().str.strip()
    if 'dominant_role' not in node.columns:
        node['dominant_role'] = node[roles_ROLE_COLS].idxmax(axis=1).map(roles_ROLE_NAME_MAP)
    tfs = roles_select_tfs(node)
    if len(tfs) == 0:
        raise ValueError('No TFs selected for plotting.')
    x_positions = roles_np.arange(len(roles_ROLE_COLS)) * roles_X_GAP
    y_positions = roles_np.arange(len(tfs)) * roles_Y_GAP
    fig_h = max(7.0, roles_BASE_HEIGHT + roles_HEIGHT_PER_TF * len(tfs))
    fig, ax = roles_plt.subplots(figsize=(roles_FIG_WIDTH, fig_h))
    for i, tf in enumerate(tfs):
        for j, col in enumerate(roles_ROLE_COLS):
            role_name = roles_ROLE_NAME_MAP[col]
            score = float(node.loc[tf, col]) if col in node.columns else 0.0
            score = max(0.0, min(1.0, score))
            size = roles_MIN_SIZE + (roles_MAX_SIZE - roles_MIN_SIZE) * score
            ax.scatter(x_positions[j], y_positions[i], s=size, color=roles_ROLE_COLORS[role_name], alpha=0.72, edgecolor='white', linewidth=0.8, zorder=3)
    for i, tf in enumerate(tfs):
        role = node.loc[tf, 'dominant_role']
        if role not in roles_ROLE_COLORS:
            continue
        role_to_col = {v: k for k, v in roles_ROLE_NAME_MAP.items()}
        if role not in role_to_col:
            continue
        j = roles_ROLE_COLS.index(role_to_col[role])
        ax.scatter(x_positions[j], y_positions[i], s=roles_STAR_SIZE, color='black', marker='*', zorder=6)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(roles_ROLE_LABELS, fontsize=12)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(tfs, fontsize=roles_LABEL_FONTSIZE)
    ax.invert_yaxis()
    ax.set_xlim(x_positions[0] - 0.85, x_positions[-1] + 1.45)
    ax.set_ylim(y_positions[-1] + 1.05, -1.05)
    ax.set_title('Data-driven TF role scores\nbubble size = role score; star = dominant role', fontsize=roles_TITLE_FONTSIZE, fontweight='bold', pad=14)
    ax.grid(alpha=0.18, linewidth=0.8)
    ax.set_axisbelow(True)
    legend_x = x_positions[-1] + 0.95
    legend_y0 = y_positions[int(len(tfs) * 0.7)] if len(tfs) > 4 else y_positions[-1]
    for k, score in enumerate([0.25, 0.5, 0.75]):
        size = roles_MIN_SIZE + (roles_MAX_SIZE - roles_MIN_SIZE) * score
        yy = legend_y0 + k * roles_Y_GAP * 1.1
        ax.scatter(legend_x, yy, s=size, color='#7F7F7F', alpha=0.35, edgecolor='white', linewidth=0.6, zorder=3)
        ax.text(legend_x + 0.45, yy, f'{score:.2f}', va='center', fontsize=10, color='black')
    ax.text(legend_x - 0.05, legend_y0 - roles_Y_GAP * 0.95, 'Role score', fontsize=10, ha='left', va='center')
    out_png = roles_OUT_DIR / 'Fig1_TF_role_score_bubble_heatmap_spaced.png'
    out_pdf = roles_OUT_DIR / 'Fig1_TF_role_score_bubble_heatmap_spaced.pdf'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    roles_plt.close(fig)
    print('Saved:')
    print(out_png)
    print(out_pdf)
    print('\nIf still crowded, increase X_GAP/Y_GAP or decrease MAX_SIZE at the top of this script.')
roles_OUT_DIR = OUT / 'roles' / 'OUT_DIR'
roles_OUT_DIR.mkdir(parents=True, exist_ok=True)
from pathlib import Path as circular_Path
import numpy as circular_np
import pandas as circular_pd
import matplotlib.pyplot as circular_plt
from matplotlib.patches import FancyArrowPatch as circular_FancyArrowPatch
import networkx as circular_nx
circular_BASE_DIR = circular_Path(str(DATA_DIR))
circular_ATACMOTIF_DIR = circular_BASE_DIR / '4.0result' / 'ATACmotif'
circular_IN_DIR = circular_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks'
circular_OUT_DIR = circular_ATACMOTIF_DIR / '39_stage_network_nonbar_figures'
circular_EDGE_PATH = circular_IN_DIR / '04_stage_specific_edge_activity.csv'
circular_NODE_PATH = circular_IN_DIR / '05_stage_specific_TF_role_scores.csv'
circular_TITLE_FONT_SIZE = 22
circular_NODE_LABEL_FONT_SIZE = 10
circular_LEGEND_FONT_SIZE = 12
circular_LEGEND_TITLE_SIZE = 12
circular_NODE_SIZE_BASE = 120
circular_NODE_SIZE_CORE_SCALE = 520
circular_EDGE_ALPHA = 0.62
circular_EDGE_BASE_WIDTH = 1.2
circular_EDGE_WIDTH_SCALE = 3.2
circular_ARROW_SIZE = 11
circular_FIG_SIZE = (11, 11)
circular_RADIUS_LABEL = 1.13
circular_RADIUS_NODE = 1.0
circular_TOP_N_PER_EDGE_TYPE = 12
circular_ARC_RAD = 0.2
circular_UNIFORM_EDGE_WIDTH = False
circular_DIFF_COLORS = {'Early gain\nPolyp>Normal': '#F2C14E', 'Late gain\nCRC>Polyp': '#E45756', 'CRC gain\nCRC>Normal': '#D1495B', 'CRC loss\nCRC<Normal': '#2A9D8F'}
circular_ROLE_COLORS = {'Stable_backbone': '#6A4C93', 'Early_transition': '#F2C14E', 'Late_transition': '#E45756', 'CRC_enhanced': '#D1495B', 'CRC_suppressed': '#2A9D8F', 'Other': '#B8B8B8'}

def circular_load_inputs():
    if not circular_EDGE_PATH.exists():
        raise FileNotFoundError(f'Cannot find edge file: {circular_EDGE_PATH}')
    if not circular_NODE_PATH.exists():
        raise FileNotFoundError(f'Cannot find node role file: {circular_NODE_PATH}')
    edges = circular_pd.read_csv(circular_EDGE_PATH)
    node = circular_pd.read_csv(circular_NODE_PATH, index_col=0)
    node.index = node.index.astype(str).str.upper().str.strip()
    edges['regulator_tf'] = edges['regulator_tf'].astype(str).str.upper().str.strip()
    edges['target_tf'] = edges['target_tf'].astype(str).str.upper().str.strip()
    if 'mean_score' not in edges.columns:
        edges['mean_score'] = 1.0
    edges['mean_score'] = circular_pd.to_numeric(edges['mean_score'], errors='coerce').fillna(1.0)
    return (edges, node)

def circular_select_differential_edges(edges):
    early = edges.sort_values('edge_delta_Polyp_vs_Normal', ascending=False).head(circular_TOP_N_PER_EDGE_TYPE).copy()
    early['diff_type'] = 'Early gain\nPolyp>Normal'
    early['diff_value'] = early['edge_delta_Polyp_vs_Normal']
    late = edges.sort_values('edge_delta_CRC_vs_Polyp', ascending=False).head(circular_TOP_N_PER_EDGE_TYPE).copy()
    late['diff_type'] = 'Late gain\nCRC>Polyp'
    late['diff_value'] = late['edge_delta_CRC_vs_Polyp']
    crc_gain = edges.sort_values('edge_delta_CRC_vs_Normal', ascending=False).head(circular_TOP_N_PER_EDGE_TYPE).copy()
    crc_gain['diff_type'] = 'CRC gain\nCRC>Normal'
    crc_gain['diff_value'] = crc_gain['edge_delta_CRC_vs_Normal']
    crc_loss = edges.sort_values('edge_delta_CRC_vs_Normal', ascending=True).head(circular_TOP_N_PER_EDGE_TYPE).copy()
    crc_loss['diff_type'] = 'CRC loss\nCRC<Normal'
    crc_loss['diff_value'] = -crc_loss['edge_delta_CRC_vs_Normal']
    sub = circular_pd.concat([early, late, crc_gain, crc_loss], axis=0)
    sub = sub.drop_duplicates(['regulator_tf', 'target_tf', 'diff_type']).copy()
    sub.to_csv(circular_OUT_DIR / 'Fig5_edges_used_circular_differential_map_solid.csv', index=False)
    return sub

def circular_build_graph(sub):
    G = circular_nx.DiGraph()
    for r in sub.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, diff_type=r.diff_type, mean_score=float(r.mean_score), diff_value=float(r.diff_value))
    return G

def circular_circular_positions(nodes):
    nodes = list(nodes)
    nodes = sorted(nodes)
    n = len(nodes)
    theta = circular_np.linspace(0, 2 * circular_np.pi, n, endpoint=False)
    pos = {node_name: circular_np.array([circular_RADIUS_NODE * circular_np.cos(t), circular_RADIUS_NODE * circular_np.sin(t)]) for node_name, t in zip(nodes, theta)}
    label_pos = {node_name: circular_np.array([circular_RADIUS_LABEL * circular_np.cos(t), circular_RADIUS_LABEL * circular_np.sin(t)]) for node_name, t in zip(nodes, theta)}
    return (pos, label_pos)

def circular_draw_solid_circular_map(edges, node):
    sub = circular_select_differential_edges(edges)
    G = circular_build_graph(sub)
    if G.number_of_edges() == 0:
        raise ValueError('No edges selected for plotting.')
    pos, label_pos = circular_circular_positions(G.nodes())
    fig, ax = circular_plt.subplots(figsize=circular_FIG_SIZE)
    ax.axis('off')
    ax.set_aspect('equal')
    diff_values = circular_np.array([G[u][v]['diff_value'] for u, v in G.edges()], dtype=float)
    if diff_values.max() == diff_values.min():
        edge_widths = circular_np.ones_like(diff_values) * (circular_EDGE_BASE_WIDTH + circular_EDGE_WIDTH_SCALE * 0.5)
    else:
        edge_widths = circular_EDGE_BASE_WIDTH + circular_EDGE_WIDTH_SCALE * (diff_values - diff_values.min()) / (diff_values.max() - diff_values.min())
    if circular_UNIFORM_EDGE_WIDTH:
        edge_widths = circular_np.ones_like(edge_widths) * (circular_EDGE_BASE_WIDTH + 1.2)
    for idx, (u, v, data) in enumerate(G.edges(data=True)):
        color = circular_DIFF_COLORS.get(data['diff_type'], '#888888')
        rad = circular_ARC_RAD
        if hash(u + v) % 2 == 0:
            rad = -circular_ARC_RAD
        arr = circular_FancyArrowPatch(pos[u], pos[v], connectionstyle=f'arc3,rad={rad}', arrowstyle='-|>', mutation_scale=circular_ARROW_SIZE, linewidth=float(edge_widths[idx]), linestyle='solid', color=color, alpha=circular_EDGE_ALPHA, shrinkA=8, shrinkB=8, zorder=2)
        ax.add_patch(arr)
    for nd in G.nodes():
        if nd in node.index:
            role = node.loc[nd, 'dominant_role'] if 'dominant_role' in node.columns else 'Other'
            core = float(node.loc[nd, 'core_score']) if 'core_score' in node.columns else 0.3
        else:
            role = 'Other'
            core = 0.3
        node_color = circular_ROLE_COLORS.get(role, circular_ROLE_COLORS['Other'])
        node_size = circular_NODE_SIZE_BASE + circular_NODE_SIZE_CORE_SCALE * core
        ax.scatter(pos[nd][0], pos[nd][1], s=node_size, color=node_color, edgecolor='white', linewidth=0.9, zorder=5)
    for nd in G.nodes():
        x, y = label_pos[nd]
        ha = 'left' if x >= 0 else 'right'
        va = 'center'
        ax.text(x, y, nd, fontsize=circular_NODE_LABEL_FONT_SIZE, ha=ha, va=va, color='#111111', zorder=6)
    edge_handles = [circular_plt.Line2D([0], [0], color=color, lw=4, linestyle='solid', label=label.replace('\n', ' ')) for label, color in circular_DIFF_COLORS.items()]
    role_handles = [circular_plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markeredgecolor='white', markersize=9, label=role.replace('_', ' ')) for role, color in circular_ROLE_COLORS.items() if role != 'Other']
    leg1 = ax.legend(handles=edge_handles, loc='lower center', bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False, fontsize=circular_LEGEND_FONT_SIZE, title='Stage-differential edge type', title_fontsize=circular_LEGEND_TITLE_SIZE)
    ax.add_artist(leg1)
    ax.legend(handles=role_handles, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False, fontsize=circular_LEGEND_FONT_SIZE, title='TF role', title_fontsize=circular_LEGEND_TITLE_SIZE)
    ax.set_title('Circular map of stage-differential TF→TF relationships', fontsize=circular_TITLE_FONT_SIZE, fontweight='bold', pad=18)
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    out_png = circular_OUT_DIR / 'Fig5_circular_stage_differential_edge_map_solid.png'
    out_pdf = circular_OUT_DIR / 'Fig5_circular_stage_differential_edge_map_solid.pdf'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    circular_plt.close(fig)
    print('Saved:')
    print(out_png)
    print(out_pdf)

##### 图5D #####
def circular_main():
    edges, node = circular_load_inputs()
    circular_draw_solid_circular_map(edges, node)
circular_OUT_DIR = OUT / 'circular' / 'OUT_DIR'
circular_OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    edges, _ = network_load_high_conf_edges()
    core, _ = network_load_core_score_table()
    cent = network_compute_centrality(edges, core)
    stage, _ = network_load_stage_edges()
    activity = network_compute_stage_tf_activity(stage)

    def draw_a():
        roles_main()
        panel_file('A', roles_OUT_DIR / 'Fig1_TF_role_score_bubble_heatmap_spaced.png')
    attempt('A', draw_a)
    for label, fn in [('B', lambda path: network_fig3D_dynamic_edge_heatmap(stage, path)), ('C', lambda path: network_fig3C_early_late_scatter(activity, core, cent, path))]:
        path = OUT / (label + '.png')
        attempt(label, lambda label=label, fn=fn, path=path: (fn(path), panel_file(label, path)))

    def draw_d():
        circular_main()
        panel_file('D', circular_OUT_DIR / 'Fig5_circular_stage_differential_edge_map_solid.png')
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
