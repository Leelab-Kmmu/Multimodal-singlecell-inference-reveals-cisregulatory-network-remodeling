# -*- coding: utf-8 -*-
##### 图4 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = '图4'
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

def network_get_top_core_tfs(cent, n=20):
    return cent.sort_values('core_score', ascending=False).head(n).index.tolist()

##### 图4A #####
def network_fig2A_network(edges, cent, out_path):
    top_core = network_get_top_core_tfs(cent, 25)
    e1 = edges.sort_values('mean_score', ascending=False).head(550)
    e2 = edges[edges['regulator_tf'].isin(top_core) | edges['target_tf'].isin(top_core)].copy()
    sub = network_pd.concat([e1, e2], axis=0).drop_duplicates(['regulator_tf', 'target_tf'])
    if sub.shape[0] > 1000:
        sub = sub.sort_values('mean_score', ascending=False).head(1000)
    G = network_build_graph(sub)
    if G.number_of_edges() == 0:
        print('No edges for Fig2A.')
        return
    pos = network_nx.spring_layout(G, seed=42, k=0.55, iterations=160, weight='weight')
    fig, ax = network_plt.subplots(figsize=(10.5, 8.2))
    ax.axis('off')
    weights = network_np.array([G[u][v].get('weight', 1.0) for u, v in G.edges()], dtype=float)
    widths = 0.25 + 2.2 * network_minmax01(network_pd.Series(weights)).values
    network_nx.draw_networkx_edges(G, pos, ax=ax, width=widths, edge_color='#A9A9A9', alpha=0.28, arrows=True, arrowsize=5, connectionstyle='arc3,rad=0.05')
    node_sizes = []
    node_colors = []
    for n in G.nodes():
        cs = float(cent.loc[n, 'core_score']) if n in cent.index else 0
        node_sizes.append(30 + 700 * cs)
        node_colors.append(cs)
    nodes = network_nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors, cmap='viridis', edgecolors='white', linewidths=0.7)
    label_tfs = [tf for tf in top_core[:18] if tf in G.nodes()]
    labels = {tf: tf for tf in label_tfs}
    network_nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold', ax=ax)
    cbar = network_plt.colorbar(nodes, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('Core score', fontsize=10)
    ax.set_title('Final high-confidence TF-TF regulatory network', fontsize=15, fontweight='bold')
    network_savefig(fig, out_path)

##### 图4C #####
def network_fig2B_core_ranking(cent, out_path, n=25):
    df = cent.sort_values('core_score', ascending=False).head(n).iloc[::-1].copy()
    fig, ax = network_plt.subplots(figsize=(7.2, max(6, 0.28 * n)))
    y = network_np.arange(df.shape[0])
    ax.hlines(y, 0, df['core_score'], color='#D0D0D0', linewidth=2)
    colors = ['#D1495B' if tf in network_FOCUS_TFS else '#4C78A8' for tf in df.index]
    ax.scatter(df['core_score'], y, s=55, color=colors, edgecolor='white', linewidth=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index, fontsize=9)
    ax.set_xlabel('Core score')
    ax.set_title('Core TF ranking in the high-confidence network', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    ax.set_xlim(0, max(1.0, df['core_score'].max() * 1.05))
    network_savefig(fig, out_path)

##### 图4B #####
def network_fig2C_centrality_heatmap(cent, out_path, n=25):
    cols = ['out_degree', 'in_degree', 'weighted_out_degree', 'weighted_in_degree', 'pagerank', 'core_score']
    df = cent.sort_values('core_score', ascending=False).head(n)[cols].copy()
    scaled = df.copy()
    for c in cols:
        scaled[c] = network_minmax01(scaled[c]).values
    fig, ax = network_plt.subplots(figsize=(7.6, max(6, 0.3 * n)))
    im = ax.imshow(scaled.values, aspect='auto', cmap='YlGnBu', vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(['Out-degree', 'In-degree', 'Weighted\nout-degree', 'Weighted\nin-degree', 'PageRank', 'Core\nscore'], rotation=30, ha='right')
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index, fontsize=8)
    ax.set_title('Centrality profiles of top core TFs', fontsize=14, fontweight='bold')
    for i in range(df.shape[0]):
        for j in range(len(cols)):
            if scaled.iloc[i, j] > 0.55:
                color = 'white'
            else:
                color = 'black'
            val = df.iloc[i, j]
            text = f'{val:.2f}' if cols[j] in ['pagerank', 'core_score'] else f'{val:.0f}'
            ax.text(j, i, text, ha='center', va='center', fontsize=6.5, color=color)
    cbar = network_plt.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label('Column-scaled value', fontsize=10)
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
    edges, _ = network_load_high_conf_edges()
    core, _ = network_load_core_score_table()
    cent = network_compute_centrality(edges, core)
    for label, fn in [('A', network_fig2A_network), ('B', network_fig2C_centrality_heatmap), ('C', network_fig2B_core_ranking)]:
        path = OUT / (label + '.png')
        attempt(label, lambda label=label, fn=fn, path=path: (fn(edges, cent, path) if label == 'A' else fn(cent, path), panel_file(label, path)))
if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        missing('运行', f'{type(exc).__name__}: {exc}')
        traceback.print_exc()
    finally:
        finish(3)

if __name__ == "__main__" and ERRORS:
    raise SystemExit(1)
