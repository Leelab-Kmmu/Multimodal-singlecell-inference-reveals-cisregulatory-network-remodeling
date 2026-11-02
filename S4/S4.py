# -*- coding: utf-8 -*-
##### S4 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = 'S4'
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

##### S4A #####
def network_fig4D_patz1_target_heatmap(stage_edges, tf_stage, out_path, top_targets=30):
    pat = stage_edges[stage_edges['regulator_tf'] == 'PATZ1'].sort_values('mean_score', ascending=False).head(top_targets).copy() if stage_edges is not None else network_pd.DataFrame()
    if pat.empty:
        print('Skip Fig4D: no PATZ1 outgoing edges.')
        return
    stage_cols = network_infer_stage_activity_cols(pat)
    if len(stage_cols) == 3:
        mat = pat[[stage_cols[s] for s in network_STAGES]].copy()
        mat.columns = network_STAGES
        mat.index = pat['target_tf'].values
    elif not tf_stage.empty:
        targets = [t for t in pat['target_tf'] if t in tf_stage.index]
        mat = tf_stage.loc[targets, network_STAGES].copy()
    else:
        print('Skip Fig4D: no stage activities.')
        return
    mat = mat.groupby(mat.index).mean()
    z = network_zscore_rows(mat)
    z.to_csv(network_TABLE_DIR / 'PATZ1_target_stage_activity_matrix.csv', encoding='utf-8-sig')
    fig, ax = network_plt.subplots(figsize=(5.6, max(5, 0.25 * z.shape[0])))
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(network_STAGES)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index, fontsize=7)
    ax.set_title('Stage activity of PATZ1 target TFs', fontsize=13, fontweight='bold')
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

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    edges, _ = network_load_high_conf_edges()
    core, _ = network_load_core_score_table()
    cent = network_compute_centrality(edges, core)
    stage, _ = network_load_stage_edges()
    activity = network_compute_stage_tf_activity(stage)
    path = OUT / 'A.png'
    attempt('A', lambda: (network_fig4D_patz1_target_heatmap(stage, activity, path), panel_file('A', path)))
    missing('B', '相关分析有Enrichr点图函数，但论文为横向柱图，最终版原始绘图函数未定位。')
    missing('C', '未找到论文生物学主题汇总柱图的原始绘图函数。')
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
