# -*- coding: utf-8 -*-
"""论文S4：同图子图代码合并。独立文件，不读取原脚本。缺失处不生成替代图。"""
from pathlib import Path
import os, sys, types, traceback, json
# 只需修改此处，或设置环境变量 CRC_DATA_DIR / CRC_OUTPUT_DIR。
DATA_DIR = Path(os.environ.get("CRC_DATA_DIR", r"E:\毕业设计\bishedata"))
OUTPUT_ROOT = Path(os.environ.get("CRC_OUTPUT_DIR", r"C:\Users\WangXinyue\Desktop\Wang_CRC_ScientificReports_Clean_Revised_Package (1)\26.10.14\111"))
FIGURE_NAME = 'S4'
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

# 来源：plot_final_results_no_celltype.py
SOURCE_network = r'''"""
Generate final result figures WITHOUT cell-type annotation.

Figures:
Figure 2: High-confidence TF-TF network and core TFs
    A. final high-confidence TF-TF network
    B. core TF ranking
    C. centrality heatmap
    D. core TF ego networks

Figure 3: Normal-Polyp-CRC stage-specific network remodeling
    A. stage-specific network similarity heatmap
    B. TF role score profile plot
    C. early vs late transition scatter
    D. dynamic edge heatmap
    E. circular differential edge map

Figure 4: PATZ1 deep analysis
    A. PATZ1 core score ranking
    B. PATZ1 activity trajectory
    C. PATZ1-centered ego network
    D. PATZ1 target TF stage activity heatmap
    E. PATZ1 novel candidate edges

Figure 5: Representative TFs and functional interpretation
    A. KLF5 / SP5 / KLF4 activity trajectories
    B. representative TF ego networks
    C. novel edge summary
    D. enrichment dotplot if R enrichment output exists; otherwise exports gene sets for R.

Run:
    D:\\Python\\python.exe E:\\毕业设计\\scmultiome_project\\第四次\\plot_final_results_no_celltype.py

Dependencies:
    D:\\Python\\python.exe -m pip install pandas numpy matplotlib networkx scipy scikit-learn

Output:
    E:\\毕业设计\x08ishedata\x04.0result\\ATACmotif'_final_results_no_celltype
"""
from pathlib import Path
import warnings
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')
BASE_DIR = Path('E:\\毕业设计\\bishedata')
ATACMOTIF_DIR = BASE_DIR / '4.0result' / 'ATACmotif'
OUT_ROOT = ATACMOTIF_DIR / '47_final_results_no_celltype'
FIG2_DIR = OUT_ROOT / 'Figure2_high_conf_network_coreTF'
FIG3_DIR = OUT_ROOT / 'Figure3_stage_network_remodeling'
FIG4_DIR = OUT_ROOT / 'Figure4_PATZ1_deep_analysis'
FIG5_DIR = OUT_ROOT / 'Figure5_representative_TF_enrichment'
TABLE_DIR = OUT_ROOT / 'tables'
ENRICH_DIR = FIG5_DIR / 'gene_sets_for_enrichment'
HIGH_CONF_EDGE_CANDIDATES = [ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '07_high_confidence_edges_with_stage_weights.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '01_high_confidence_edges_used.csv', ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'high_confidence_edges_multimodal.csv', ATACMOTIF_DIR / '28_scMultiomeGRN_edgeCNN_optimized' / 'edgeCNN_result_analysis_figures' / 'edgeCNN_high_conf_edges_used.csv']
STAGE_EDGE_CANDIDATES = [ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '07_high_confidence_edges_with_stage_weights.csv']
CORE_SCORE_CANDIDATES = [ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '05_stage_specific_TF_role_scores.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '04_coreTF_scores_recomputed.csv', ATACMOTIF_DIR / '30_coreTF_backbone_refined' / '03_recomputed_coreTF_scores.csv']
INITIAL_ADJ_CANDIDATES = [ATACMOTIF_DIR / '10_TF_TF_adjacency_safe' / 'TF_TF_initial_adjacency_matrix.csv']
PREDICTION_CANDIDATES = [ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'predictions_all_folds_all_label_edges.csv', ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'all_predictions_multimodal.csv']
OPTIONAL_ACTIVITY_CANDIDATES = [ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / 'TF_stage_multimodal_activity.csv', ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '03_TF_stage_activity_scores.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '03_TF_stage_activity_scores.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / 'TF_stage_multimodal_activity.csv']
STAGES = ['Normal', 'Polyp', 'CRC']
FOCUS_TFS = ['PATZ1', 'KLF5', 'SP5', 'KLF4', 'ZNF148', 'ZNF281', 'EGR1', 'MAZ']
REPRESENTATIVE_TFS = ['KLF5', 'SP5', 'KLF4']
DPI = 300
COLOR_STAGE = {'Normal': '#4C78A8', 'Polyp': '#F58518', 'CRC': '#E45756'}
COLOR_ROLE = {'stable_backbone_score': '#6A4C93', 'early_transition_score': '#F2C14E', 'late_transition_score': '#E45756', 'crc_enhanced_score': '#D1495B', 'crc_suppressed_score': '#2A9D8F'}
ROLE_LABELS = {'stable_backbone_score': 'Stable\nbackbone', 'early_transition_score': 'Early\ntransition', 'late_transition_score': 'Late\ntransition', 'crc_enhanced_score': 'CRC\nenhanced', 'crc_suppressed_score': 'CRC\nsuppressed'}

def first_existing(paths):
    for p in paths:
        if Path(p).exists():
            return Path(p)
    return None

def norm_gene(x):
    return str(x).strip().upper()

def savefig(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    plt.close(fig)
    print('Saved:', path)

def minmax01(x):
    x = pd.Series(x, dtype=float).replace([np.inf, -np.inf], np.nan)
    if x.dropna().empty:
        return pd.Series(np.zeros(len(x)), index=x.index)
    x = x.fillna(x.median())
    if float(x.max()) == float(x.min()):
        return pd.Series(np.ones(len(x)) * 0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())

def zscore_rows(df):
    out = df.astype(float).copy()
    for idx in out.index:
        x = out.loc[idx].astype(float)
        sd = x.std()
        out.loc[idx] = 0 if sd == 0 or pd.isna(sd) else (x - x.mean()) / sd
    return out.replace([np.inf, -np.inf], 0).fillna(0)

def clean_edges(df):
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
    df['regulator_tf'] = df['regulator_tf'].map(norm_gene)
    df['target_tf'] = df['target_tf'].map(norm_gene)
    if 'mean_score' not in df.columns:
        for c in ['score', 'prob', 'y_prob', 'pred_prob', 'prediction_score', 'weight', 'consensus_score']:
            if c in df.columns:
                df['mean_score'] = pd.to_numeric(df[c], errors='coerce')
                break
    if 'mean_score' not in df.columns:
        df['mean_score'] = 1.0
    df['mean_score'] = pd.to_numeric(df['mean_score'], errors='coerce').fillna(1.0)
    df = df[df['regulator_tf'] != df['target_tf']].copy()
    df = df.drop_duplicates(['regulator_tf', 'target_tf']).reset_index(drop=True)
    return df

def load_high_conf_edges():
    p = first_existing(HIGH_CONF_EDGE_CANDIDATES)
    if p is None:
        raise FileNotFoundError('Cannot find high confidence edge file.')
    edges = clean_edges(pd.read_csv(p))
    edges.to_csv(TABLE_DIR / 'high_confidence_edges_used.csv', index=False, encoding='utf-8-sig')
    print('Using high-confidence edge file:', p)
    print('High-confidence edges:', edges.shape)
    return (edges, p)

def load_stage_edges():
    p = first_existing(STAGE_EDGE_CANDIDATES)
    if p is None:
        print('No stage edge activity file found. Stage-specific plots may be limited.')
        return (None, None)
    edges = clean_edges(pd.read_csv(p))
    print('Using stage edge file:', p)
    return (edges, p)

def load_initial_edge_set():
    p = first_existing(INITIAL_ADJ_CANDIDATES)
    if p is None:
        print('Initial adjacency matrix not found; all high-confidence edges will be treated as unknown initial status.')
        return (set(), None)
    adj = pd.read_csv(p, index_col=0)
    adj.index = adj.index.astype(str).str.upper().str.strip()
    adj.columns = adj.columns.astype(str).str.upper().str.strip()
    edge_set = set()
    arr = adj.values
    rows = list(adj.index)
    cols = list(adj.columns)
    for i, u in enumerate(rows):
        nz = np.where(pd.to_numeric(pd.Series(arr[i]), errors='coerce').fillna(0).values > 0)[0]
        for j in nz:
            edge_set.add((u, cols[j]))
    print('Using initial adjacency:', p)
    print('Initial positive edges:', len(edge_set))
    return (edge_set, p)

def add_initial_status(edges, initial_set):
    edges = edges.copy()
    if len(initial_set) == 0:
        edges['initial_supported'] = np.nan
        edges['is_novel'] = np.nan
    else:
        edges['initial_supported'] = [(u, v) in initial_set for u, v in zip(edges['regulator_tf'], edges['target_tf'])]
        edges['is_novel'] = ~edges['initial_supported']
    return edges

def build_graph(edges):
    G = nx.DiGraph()
    for r in edges.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, weight=float(r.mean_score))
    return G

def compute_centrality(edges, core_score_table=None):
    G = build_graph(edges)
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
        pagerank = nx.pagerank(G, weight='weight')
    except Exception:
        pagerank = {n: 0.0 for n in nodes}
    cent = pd.DataFrame(index=nodes)
    cent['out_degree'] = cent.index.map(out_deg).fillna(0).astype(float)
    cent['in_degree'] = cent.index.map(in_deg).fillna(0).astype(float)
    cent['weighted_out_degree'] = cent.index.map(weighted_out).fillna(0).astype(float)
    cent['weighted_in_degree'] = cent.index.map(weighted_in).fillna(0).astype(float)
    cent['pagerank'] = cent.index.map(pagerank).fillna(0).astype(float)
    if core_score_table is not None:
        for c in core_score_table.columns:
            if c not in cent.columns:
                cent[c] = cent.index.map(core_score_table[c]).fillna(np.nan)
    if 'core_score' not in cent.columns or cent['core_score'].isna().all():
        cent['core_score'] = 0.35 * minmax01(cent['out_degree']) + 0.35 * minmax01(cent['weighted_out_degree']) + 0.2 * minmax01(cent['pagerank']) + 0.1 * minmax01(cent['in_degree'])
    else:
        cent['core_score'] = pd.to_numeric(cent['core_score'], errors='coerce')
        missing = cent['core_score'].isna()
        if missing.any():
            fallback = 0.35 * minmax01(cent['out_degree']) + 0.35 * minmax01(cent['weighted_out_degree']) + 0.2 * minmax01(cent['pagerank']) + 0.1 * minmax01(cent['in_degree'])
            cent.loc[missing, 'core_score'] = fallback.loc[missing]
    cent = cent.sort_values('core_score', ascending=False)
    cent.to_csv(TABLE_DIR / 'TF_centrality_core_scores.csv', encoding='utf-8-sig')
    return cent

def load_core_score_table():
    p = first_existing(CORE_SCORE_CANDIDATES)
    if p is None:
        print('Core score file not found. Core scores will be computed from network.')
        return (None, None)
    df = pd.read_csv(p, index_col=0)
    df.index = df.index.astype(str).str.upper().str.strip()
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='ignore')
    print('Using core score file:', p)
    return (df, p)

def infer_stage_activity_cols(df):
    """Infer Normal/Polyp/CRC stage activity columns in edge or TF table."""
    if df is None:
        return {}
    cols = list(df.columns)
    mapping = {}
    for stage in STAGES:
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

def infer_delta_cols(df):
    if df is None:
        return {}
    cols = list(df.columns)
    candidates = {'Polyp_vs_Normal': ['edge_delta_Polyp_vs_Normal', 'delta_Polyp_vs_Normal', 'Polyp_vs_Normal_delta', 'delta_polyp_vs_normal', 'Polyp_vs_Normal'], 'CRC_vs_Polyp': ['edge_delta_CRC_vs_Polyp', 'delta_CRC_vs_Polyp', 'CRC_vs_Polyp_delta', 'delta_crc_vs_polyp', 'CRC_vs_Polyp'], 'CRC_vs_Normal': ['edge_delta_CRC_vs_Normal', 'delta_CRC_vs_Normal', 'CRC_vs_Normal_delta', 'delta_crc_vs_normal', 'CRC_vs_Normal']}
    out = {}
    for key, names in candidates.items():
        for n in names:
            if n in cols:
                out[key] = n
                break
        if key not in out:
            a, b = key.lower().split('_vs_')
            for c in cols:
                cl = c.lower()
                if 'delta' in cl and a in cl and (b in cl):
                    out[key] = c
                    break
    return out

def compute_stage_tf_activity(stage_edges):
    """Return TF x stage composite out-activity matrix."""
    if stage_edges is None:
        return pd.DataFrame()
    stage_cols = infer_stage_activity_cols(stage_edges)
    if len(stage_cols) == 3:
        tmp = stage_edges.copy()
        for st, col in stage_cols.items():
            tmp[col] = pd.to_numeric(tmp[col], errors='coerce').fillna(0)
        mat = tmp.groupby('regulator_tf')[[stage_cols[s] for s in STAGES]].sum()
        mat.columns = STAGES
        mat = mat.reindex(sorted(set(stage_edges['regulator_tf']).union(stage_edges['target_tf']))).fillna(0)
        mat.to_csv(TABLE_DIR / 'TF_stage_out_activity_matrix.csv', encoding='utf-8-sig')
        return mat
    print('Could not infer explicit stage activity columns from stage edge table.')
    return pd.DataFrame()

def load_optional_multimodal_activity():
    """Try to load RNA/RP/composite stage activity table if the user has already generated it."""
    p = first_existing(OPTIONAL_ACTIVITY_CANDIDATES)
    if p is None:
        return (None, None)
    df = pd.read_csv(p)
    cols = [c.lower() for c in df.columns]
    if 'tf' not in cols:
        if df.columns[0].lower() in ['gene', 'regulator_tf', 'transcription_factor']:
            df = df.rename(columns={df.columns[0]: 'tf'})
        else:
            return (None, None)
    for c in df.columns:
        if c.lower() == 'tf':
            df = df.rename(columns={c: 'tf'})
    df['tf'] = df['tf'].map(norm_gene)
    print('Using optional multimodal activity table:', p)
    return (df, p)

def get_top_core_tfs(cent, n=20):
    return cent.sort_values('core_score', ascending=False).head(n).index.tolist()

def fig2A_network(edges, cent, out_path):
    top_core = get_top_core_tfs(cent, 25)
    e1 = edges.sort_values('mean_score', ascending=False).head(550)
    e2 = edges[edges['regulator_tf'].isin(top_core) | edges['target_tf'].isin(top_core)].copy()
    sub = pd.concat([e1, e2], axis=0).drop_duplicates(['regulator_tf', 'target_tf'])
    if sub.shape[0] > 1000:
        sub = sub.sort_values('mean_score', ascending=False).head(1000)
    G = build_graph(sub)
    if G.number_of_edges() == 0:
        print('No edges for Fig2A.')
        return
    pos = nx.spring_layout(G, seed=42, k=0.55, iterations=160, weight='weight')
    fig, ax = plt.subplots(figsize=(10.5, 8.2))
    ax.axis('off')
    weights = np.array([G[u][v].get('weight', 1.0) for u, v in G.edges()], dtype=float)
    widths = 0.25 + 2.2 * minmax01(pd.Series(weights)).values
    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, edge_color='#A9A9A9', alpha=0.28, arrows=True, arrowsize=5, connectionstyle='arc3,rad=0.05')
    node_sizes = []
    node_colors = []
    for n in G.nodes():
        cs = float(cent.loc[n, 'core_score']) if n in cent.index else 0
        node_sizes.append(30 + 700 * cs)
        node_colors.append(cs)
    nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors, cmap='viridis', edgecolors='white', linewidths=0.7)
    label_tfs = [tf for tf in top_core[:18] if tf in G.nodes()]
    labels = {tf: tf for tf in label_tfs}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold', ax=ax)
    cbar = plt.colorbar(nodes, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('Core score', fontsize=10)
    ax.set_title('Final high-confidence TF-TF regulatory network', fontsize=15, fontweight='bold')
    savefig(fig, out_path)

def fig2B_core_ranking(cent, out_path, n=25):
    df = cent.sort_values('core_score', ascending=False).head(n).iloc[::-1].copy()
    fig, ax = plt.subplots(figsize=(7.2, max(6, 0.28 * n)))
    y = np.arange(df.shape[0])
    ax.hlines(y, 0, df['core_score'], color='#D0D0D0', linewidth=2)
    colors = ['#D1495B' if tf in FOCUS_TFS else '#4C78A8' for tf in df.index]
    ax.scatter(df['core_score'], y, s=55, color=colors, edgecolor='white', linewidth=0.6, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index, fontsize=9)
    ax.set_xlabel('Core score')
    ax.set_title('Core TF ranking in the high-confidence network', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    ax.set_xlim(0, max(1.0, df['core_score'].max() * 1.05))
    savefig(fig, out_path)

def fig2C_centrality_heatmap(cent, out_path, n=25):
    cols = ['out_degree', 'in_degree', 'weighted_out_degree', 'weighted_in_degree', 'pagerank', 'core_score']
    df = cent.sort_values('core_score', ascending=False).head(n)[cols].copy()
    scaled = df.copy()
    for c in cols:
        scaled[c] = minmax01(scaled[c]).values
    fig, ax = plt.subplots(figsize=(7.6, max(6, 0.3 * n)))
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
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label('Column-scaled value', fontsize=10)
    savefig(fig, out_path)

def plot_ego_network_for_tf(ax, edges, cent, tf, max_targets=22, title_suffix=''):
    sub = edges[edges['regulator_tf'] == tf].sort_values('mean_score', ascending=False).head(max_targets)
    if sub.empty:
        ax.text(0.5, 0.5, f'No outgoing edges for {tf}', ha='center', va='center')
        ax.axis('off')
        return
    G = nx.DiGraph()
    for r in sub.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, weight=float(r.mean_score))
    targets = [v for u, v in G.edges() if u == tf]
    pos = {tf: np.array([0.0, 0.0])}
    angles = np.linspace(0, 2 * np.pi, len(targets), endpoint=False)
    for t, a in zip(targets, angles):
        pos[t] = np.array([np.cos(a), np.sin(a)])
    weights = np.array([G[u][v]['weight'] for u, v in G.edges()], dtype=float)
    widths = 0.8 + 3.0 * minmax01(pd.Series(weights)).values
    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, edge_color='#999999', alpha=0.55, arrows=True, arrowsize=9)
    sizes = []
    colors = []
    for n in G.nodes():
        if n == tf:
            sizes.append(700)
            colors.append('#D1495B')
        else:
            cs = float(cent.loc[n, 'core_score']) if n in cent.index else 0.1
            sizes.append(80 + 450 * cs)
            colors.append('#4C78A8')
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=colors, edgecolors='white', linewidths=0.8)
    labels = {n: n for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, font_weight='bold', ax=ax)
    ax.set_title(f'{tf}-centered module{title_suffix}', fontsize=11, fontweight='bold')
    ax.axis('off')

def fig2D_core_ego_networks(edges, cent, out_path):
    top = [tf for tf in FOCUS_TFS if tf in cent.index]
    more = [tf for tf in get_top_core_tfs(cent, 12) if tf not in top]
    tfs = (top + more)[:4]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.ravel()
    for ax, tf in zip(axes, tfs):
        plot_ego_network_for_tf(ax, edges, cent, tf, max_targets=20)
    for ax in axes[len(tfs):]:
        ax.axis('off')
    fig.suptitle('Representative core TF-centered regulatory modules', fontsize=15, fontweight='bold')
    savefig(fig, out_path)

def make_figure2(edges, cent):
    print('\n=== Figure 2 ===')
    fig2A_network(edges, cent, FIG2_DIR / 'Fig2A_final_high_confidence_TF_TF_network.png')
    fig2B_core_ranking(cent, FIG2_DIR / 'Fig2B_core_TF_ranking_lollipop.png')
    fig2C_centrality_heatmap(cent, FIG2_DIR / 'Fig2C_centrality_heatmap.png')
    fig2D_core_ego_networks(edges, cent, FIG2_DIR / 'Fig2D_core_TF_ego_networks.png')

def fig3A_stage_similarity(tf_stage, out_path):
    if tf_stage.empty or tf_stage.shape[1] < 2:
        print('Skip Fig3A: no TF stage activity matrix.')
        return
    corr = tf_stage[STAGES].corr(method='pearson')
    corr.to_csv(TABLE_DIR / 'stage_network_similarity_pearson.csv', encoding='utf-8-sig')
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    im = ax.imshow(corr.values, cmap='magma', vmin=0, vmax=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels(STAGES)
    ax.set_yticks(range(3))
    ax.set_yticklabels(STAGES)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f'{corr.iloc[i, j]:.2f}', ha='center', va='center', color='white' if corr.iloc[i, j] > 0.55 else 'black', fontsize=12, fontweight='bold')
    ax.set_title('Stage-specific network similarity', fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label('Pearson correlation', fontsize=10)
    savefig(fig, out_path)

def fig3B_role_profile(core_scores, cent, out_path, n=14):
    role_cols = [c for c in COLOR_ROLE if core_scores is not None and c in core_scores.columns]
    if len(role_cols) == 0:
        print('Skip Fig3B: role score columns not found.')
        return
    selected = []
    for c in role_cols:
        selected += core_scores.sort_values(c, ascending=False).head(4).index.tolist()
    selected += get_top_core_tfs(cent, 10)
    selected = list(dict.fromkeys([x for x in selected if x in core_scores.index]))[:n]
    df = core_scores.loc[selected, role_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    df.to_csv(TABLE_DIR / 'TF_role_scores_used_for_Fig3B.csv', encoding='utf-8-sig')
    fig, ax = plt.subplots(figsize=(9.5, max(5.5, 0.35 * len(selected))))
    x = np.arange(len(role_cols))
    for i, tf in enumerate(selected):
        y = df.loc[tf].values.astype(float)
        ax.plot(x, y, color='#B8B8B8', alpha=0.35, linewidth=1.2)
        dom = int(np.argmax(y))
        ax.scatter(x, y, s=22, color='#B8B8B8', alpha=0.7, zorder=2)
        ax.scatter(x[dom], y[dom], s=70, color=COLOR_ROLE[role_cols[dom]], edgecolor='white', linewidth=0.6, zorder=3)
        ax.text(len(role_cols) - 0.05, y[-1], f' {tf}', va='center', fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([ROLE_LABELS[c] for c in role_cols], fontsize=10)
    ax.set_ylabel('Role score')
    ax.set_title('Data-driven TF role score profiles', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    ax.set_xlim(-0.25, len(role_cols) + 1.2)
    savefig(fig, out_path)

def fig3C_early_late_scatter(tf_stage, core_scores, cent, out_path):
    if tf_stage.empty:
        print('Skip Fig3C: no TF stage activity matrix.')
        return
    df = tf_stage.copy()
    for s in STAGES:
        df[s] = pd.to_numeric(df[s], errors='coerce').fillna(0)
    df['early_delta_Polyp_vs_Normal'] = df['Polyp'] - df['Normal']
    df['late_delta_CRC_vs_Polyp'] = df['CRC'] - df['Polyp']
    df['core_score'] = df.index.map(cent['core_score']).fillna(0)
    df.to_csv(TABLE_DIR / 'TF_early_late_transition_scores.csv', encoding='utf-8-sig')
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    sizes = 25 + 350 * minmax01(df['core_score']).values
    colors = df['late_delta_CRC_vs_Polyp'].values
    sc = ax.scatter(df['early_delta_Polyp_vs_Normal'], df['late_delta_CRC_vs_Polyp'], s=sizes, c=colors, cmap='coolwarm', alpha=0.75, edgecolor='white', linewidth=0.4)
    ax.axhline(0, color='#555555', lw=0.8, ls='--')
    ax.axvline(0, color='#555555', lw=0.8, ls='--')
    ax.set_xlabel('Early transition: Polyp - Normal')
    ax.set_ylabel('Late transition: CRC - Polyp')
    ax.set_title('Early vs late TF regulatory activity transition', fontsize=14, fontweight='bold')
    label_tfs = list(dict.fromkeys(FOCUS_TFS + df.sort_values('core_score', ascending=False).head(10).index.tolist()))
    for tf in label_tfs:
        if tf in df.index:
            ax.text(df.loc[tf, 'early_delta_Polyp_vs_Normal'], df.loc[tf, 'late_delta_CRC_vs_Polyp'], f' {tf}', fontsize=8)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label('Late transition delta', fontsize=10)
    ax.grid(alpha=0.2)
    savefig(fig, out_path)

def fig3D_dynamic_edge_heatmap(stage_edges, out_path, top_n=35):
    if stage_edges is None:
        print('Skip Fig3D: no stage edge table.')
        return
    stage_cols = infer_stage_activity_cols(stage_edges)
    if len(stage_cols) != 3:
        print('Skip Fig3D: cannot infer stage edge activity columns.')
        return
    df = stage_edges.copy()
    for s, c in stage_cols.items():
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)
    vals = df[[stage_cols[s] for s in STAGES]].copy()
    vals.columns = STAGES
    df['dynamic_range'] = vals.max(axis=1) - vals.min(axis=1)
    sub = df.sort_values('dynamic_range', ascending=False).head(top_n).copy()
    mat = sub[[stage_cols[s] for s in STAGES]].copy()
    mat.columns = STAGES
    mat.index = sub['regulator_tf'] + '→' + sub['target_tf']
    z = zscore_rows(mat)
    sub.to_csv(TABLE_DIR / 'dynamic_edges_used_for_Fig3D.csv', index=False, encoding='utf-8-sig')
    fig, ax = plt.subplots(figsize=(5.8, max(6, 0.22 * top_n)))
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(STAGES)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index, fontsize=6.5)
    ax.set_title('Dynamic TF-TF edge activity across stages', fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label('Row z-score', fontsize=10)
    savefig(fig, out_path)

def fig3E_circular_differential_edges(stage_edges, out_path, top_n_per_type=10):
    if stage_edges is None:
        print('Skip Fig3E: no stage edge table.')
        return
    delta_cols = infer_delta_cols(stage_edges)
    df = stage_edges.copy()
    stage_cols = infer_stage_activity_cols(df)
    if len(delta_cols) < 3 and len(stage_cols) == 3:
        df['delta_Polyp_vs_Normal_auto'] = pd.to_numeric(df[stage_cols['Polyp']], errors='coerce').fillna(0) - pd.to_numeric(df[stage_cols['Normal']], errors='coerce').fillna(0)
        df['delta_CRC_vs_Polyp_auto'] = pd.to_numeric(df[stage_cols['CRC']], errors='coerce').fillna(0) - pd.to_numeric(df[stage_cols['Polyp']], errors='coerce').fillna(0)
        df['delta_CRC_vs_Normal_auto'] = pd.to_numeric(df[stage_cols['CRC']], errors='coerce').fillna(0) - pd.to_numeric(df[stage_cols['Normal']], errors='coerce').fillna(0)
        delta_cols = {'Polyp_vs_Normal': 'delta_Polyp_vs_Normal_auto', 'CRC_vs_Polyp': 'delta_CRC_vs_Polyp_auto', 'CRC_vs_Normal': 'delta_CRC_vs_Normal_auto'}
    if len(delta_cols) < 3:
        print('Skip Fig3E: cannot infer delta columns.')
        return
    sets = []
    early = df.sort_values(delta_cols['Polyp_vs_Normal'], ascending=False).head(top_n_per_type).copy()
    early['diff_type'] = 'Early gain\nPolyp>Normal'
    early['diff_value'] = early[delta_cols['Polyp_vs_Normal']]
    sets.append(early)
    late = df.sort_values(delta_cols['CRC_vs_Polyp'], ascending=False).head(top_n_per_type).copy()
    late['diff_type'] = 'Late gain\nCRC>Polyp'
    late['diff_value'] = late[delta_cols['CRC_vs_Polyp']]
    sets.append(late)
    gain = df.sort_values(delta_cols['CRC_vs_Normal'], ascending=False).head(top_n_per_type).copy()
    gain['diff_type'] = 'CRC gain\nCRC>Normal'
    gain['diff_value'] = gain[delta_cols['CRC_vs_Normal']]
    sets.append(gain)
    loss = df.sort_values(delta_cols['CRC_vs_Normal'], ascending=True).head(top_n_per_type).copy()
    loss['diff_type'] = 'CRC loss\nCRC<Normal'
    loss['diff_value'] = -loss[delta_cols['CRC_vs_Normal']]
    sets.append(loss)
    sub = pd.concat(sets, axis=0).drop_duplicates(['regulator_tf', 'target_tf', 'diff_type'])
    sub.to_csv(TABLE_DIR / 'differential_edges_used_for_Fig3E.csv', index=False, encoding='utf-8-sig')
    G = nx.DiGraph()
    for r in sub.itertuples(index=False):
        G.add_edge(r.regulator_tf, r.target_tf, diff_type=r.diff_type, diff_value=float(r.diff_value), mean_score=float(r.mean_score))
    if G.number_of_edges() == 0:
        print('Skip Fig3E: selected graph empty.')
        return
    nodes = sorted(G.nodes())
    theta = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
    pos = {n: np.array([np.cos(t), np.sin(t)]) for n, t in zip(nodes, theta)}
    labpos = {n: np.array([1.14 * np.cos(t), 1.14 * np.sin(t)]) for n, t in zip(nodes, theta)}
    colors = {'Early gain\nPolyp>Normal': '#F2C14E', 'Late gain\nCRC>Polyp': '#E45756', 'CRC gain\nCRC>Normal': '#D1495B', 'CRC loss\nCRC<Normal': '#2A9D8F'}
    vals = np.array([G[u][v]['diff_value'] for u, v in G.edges()], dtype=float)
    widths = 0.9 + 3.5 * minmax01(pd.Series(vals)).values
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis('off')
    ax.set_aspect('equal')
    for idx, (u, v, d) in enumerate(G.edges(data=True)):
        rad = 0.18 if hash(u + v) % 2 else -0.18
        patch = FancyArrowPatch(pos[u], pos[v], connectionstyle=f'arc3,rad={rad}', arrowstyle='-|>', mutation_scale=10, linewidth=float(widths[idx]), linestyle='solid', color=colors.get(d['diff_type'], '#888888'), alpha=0.62, shrinkA=8, shrinkB=8, zorder=2)
        ax.add_patch(patch)
    ax.scatter([pos[n][0] for n in nodes], [pos[n][1] for n in nodes], s=90, color='#4C78A8', edgecolor='white', linewidth=0.8, zorder=3)
    for n in nodes:
        x, y = labpos[n]
        ax.text(x, y, n, ha='left' if x >= 0 else 'right', va='center', fontsize=8)
    handles = [Line2D([0], [0], color=c, lw=4, label=k.replace('\n', ' ')) for k, c in colors.items()]
    ax.legend(handles=handles, frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=9)
    ax.set_title('Stage-differential TF-TF relationships', fontsize=15, fontweight='bold')
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    savefig(fig, out_path)

def make_figure3(stage_edges, core_scores, cent, tf_stage):
    print('\n=== Figure 3 ===')
    fig3A_stage_similarity(tf_stage, FIG3_DIR / 'Fig3A_stage_network_similarity_heatmap.png')
    fig3B_role_profile(core_scores, cent, FIG3_DIR / 'Fig3B_TF_role_score_profile.png')
    fig3C_early_late_scatter(tf_stage, core_scores, cent, FIG3_DIR / 'Fig3C_early_vs_late_transition_scatter.png')
    fig3D_dynamic_edge_heatmap(stage_edges, FIG3_DIR / 'Fig3D_dynamic_edge_heatmap.png')
    fig3E_circular_differential_edges(stage_edges, FIG3_DIR / 'Fig3E_circular_differential_edge_map.png')

def fig4A_patz1_ranking(cent, out_path, n=30):
    df = cent.sort_values('core_score', ascending=False).head(n).iloc[::-1].copy()
    fig, ax = plt.subplots(figsize=(7.0, max(6, 0.25 * n)))
    y = np.arange(df.shape[0])
    colors = ['#D1495B' if tf == 'PATZ1' else '#4C78A8' for tf in df.index]
    ax.hlines(y, 0, df['core_score'], color='#D0D0D0', linewidth=1.8)
    ax.scatter(df['core_score'], y, s=[90 if tf == 'PATZ1' else 45 for tf in df.index], color=colors, edgecolor='white', linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index, fontsize=8)
    ax.set_xlabel('Core score')
    ax.set_title('PATZ1 position in core TF ranking', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    if 'PATZ1' in df.index:
        rank = cent.sort_values('core_score', ascending=False).index.get_loc('PATZ1') + 1
        ax.text(0.98, 0.05, f'PATZ1 rank = {rank}', transform=ax.transAxes, ha='right', fontsize=11, fontweight='bold', color='#D1495B')
    savefig(fig, out_path)

def plot_tf_trajectory(tf, tf_stage, optional_activity, out_path, title=None):
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    x = np.arange(len(STAGES))
    plotted = False
    if optional_activity is not None and 'tf' in optional_activity.columns:
        sub = optional_activity[optional_activity['tf'] == tf].copy()
        if 'stage' in [c.lower() for c in sub.columns]:
            for c in sub.columns:
                if c.lower() == 'stage':
                    sub = sub.rename(columns={c: 'stage'})
            metric_cols = [c for c in sub.columns if c not in ['tf', 'stage'] and pd.api.types.is_numeric_dtype(pd.to_numeric(sub[c], errors='coerce'))]
            for metric in metric_cols[:4]:
                vals = []
                for st in STAGES:
                    ss = sub[sub['stage'].astype(str).str.lower() == st.lower()]
                    vals.append(pd.to_numeric(ss[metric], errors='coerce').mean() if not ss.empty else np.nan)
                if not np.all(pd.isna(vals)):
                    ax.plot(x, vals, marker='o', linewidth=2, label=metric)
                    plotted = True
        else:
            for metric in ['RNA', 'RP', 'ATAC', 'composite', 'Composite']:
                cols = []
                for st in STAGES:
                    candidates = [c for c in sub.columns if metric.lower() in c.lower() and st.lower() in c.lower()]
                    cols.append(candidates[0] if candidates else None)
                if all((c is not None for c in cols)):
                    vals = [pd.to_numeric(sub[c], errors='coerce').mean() for c in cols]
                    ax.plot(x, vals, marker='o', linewidth=2, label=metric)
                    plotted = True
    if not tf_stage.empty and tf in tf_stage.index:
        vals = tf_stage.loc[tf, STAGES].astype(float).values
        vals_norm = minmax01(pd.Series(vals, index=STAGES)).values
        ax.plot(x, vals_norm, marker='o', linewidth=2.5, color='#D1495B', label='Composite out-activity')
        plotted = True
    if not plotted:
        ax.text(0.5, 0.5, f'No stage activity available for {tf}', ha='center', va='center')
    ax.set_xticks(x)
    ax.set_xticklabels(STAGES)
    ax.set_ylabel('Normalized activity')
    ax.set_title(title or f'{tf} activity trajectory', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)
    savefig(fig, out_path)

def fig4C_patz1_ego(edges, cent, out_path):
    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    plot_ego_network_for_tf(ax, edges, cent, 'PATZ1', max_targets=28)
    savefig(fig, out_path)

def fig4D_patz1_target_heatmap(stage_edges, tf_stage, out_path, top_targets=30):
    pat = stage_edges[stage_edges['regulator_tf'] == 'PATZ1'].sort_values('mean_score', ascending=False).head(top_targets).copy() if stage_edges is not None else pd.DataFrame()
    if pat.empty:
        print('Skip Fig4D: no PATZ1 outgoing edges.')
        return
    stage_cols = infer_stage_activity_cols(pat)
    if len(stage_cols) == 3:
        mat = pat[[stage_cols[s] for s in STAGES]].copy()
        mat.columns = STAGES
        mat.index = pat['target_tf'].values
    elif not tf_stage.empty:
        targets = [t for t in pat['target_tf'] if t in tf_stage.index]
        mat = tf_stage.loc[targets, STAGES].copy()
    else:
        print('Skip Fig4D: no stage activities.')
        return
    mat = mat.groupby(mat.index).mean()
    z = zscore_rows(mat)
    z.to_csv(TABLE_DIR / 'PATZ1_target_stage_activity_matrix.csv', encoding='utf-8-sig')
    fig, ax = plt.subplots(figsize=(5.6, max(5, 0.25 * z.shape[0])))
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(STAGES)
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index, fontsize=7)
    ax.set_title('Stage activity of PATZ1 target TFs', fontsize=13, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label('Row z-score', fontsize=10)
    savefig(fig, out_path)

def fig4E_patz1_novel_edges(edges, out_path, top_n=18):
    if 'is_novel' not in edges.columns or edges['is_novel'].isna().all():
        print('Skip Fig4E: novel status unavailable.')
        return
    sub = edges[(edges['regulator_tf'] == 'PATZ1') & (edges['is_novel'] == True)].sort_values('mean_score', ascending=False).head(top_n).copy()
    if sub.empty:
        print('No novel PATZ1 outgoing edges found; using top PATZ1 edges instead.')
        sub = edges[edges['regulator_tf'] == 'PATZ1'].sort_values('mean_score', ascending=False).head(top_n).copy()
        sub['is_novel'] = False
    sub.to_csv(TABLE_DIR / 'PATZ1_novel_candidate_edges.csv', index=False, encoding='utf-8-sig')
    df = sub.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.0, max(4.5, 0.28 * df.shape[0])))
    y = np.arange(df.shape[0])
    labels = df['regulator_tf'] + '→' + df['target_tf']
    colors = ['#D1495B' if x else '#4C78A8' for x in df['is_novel']]
    ax.hlines(y, 0, df['mean_score'], color='#D0D0D0', linewidth=1.7)
    ax.scatter(df['mean_score'], y, s=55, color=colors, edgecolor='white', linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('Mean model confidence')
    ax.set_title('PATZ1 candidate novel regulatory edges', fontsize=13, fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor='#D1495B', label='Novel vs initial motif network', markersize=7), Line2D([0], [0], marker='o', color='w', markerfacecolor='#4C78A8', label='Initial-supported / unknown', markersize=7)]
    ax.legend(handles=legend, frameon=False, fontsize=8, loc='lower right')
    savefig(fig, out_path)

def make_figure4(edges, stage_edges, cent, tf_stage, optional_activity):
    print('\n=== Figure 4 ===')
    fig4A_patz1_ranking(cent, FIG4_DIR / 'Fig4A_PATZ1_core_score_ranking.png')
    plot_tf_trajectory('PATZ1', tf_stage, optional_activity, FIG4_DIR / 'Fig4B_PATZ1_activity_trajectory.png', 'PATZ1 multi-stage activity trajectory')
    fig4C_patz1_ego(edges, cent, FIG4_DIR / 'Fig4C_PATZ1_centered_ego_network.png')
    fig4D_patz1_target_heatmap(stage_edges, tf_stage, FIG4_DIR / 'Fig4D_PATZ1_target_stage_activity_heatmap.png')
    fig4E_patz1_novel_edges(edges, FIG4_DIR / 'Fig4E_PATZ1_novel_candidate_edges.png')

def fig5A_representative_trajectories(tf_stage, optional_activity, out_path):
    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    x = np.arange(len(STAGES))
    colors = {'KLF5': '#D1495B', 'SP5': '#F58518', 'KLF4': '#2A9D8F'}
    for tf in REPRESENTATIVE_TFS:
        if not tf_stage.empty and tf in tf_stage.index:
            vals = tf_stage.loc[tf, STAGES].astype(float)
            vals = minmax01(vals).values
            ax.plot(x, vals, marker='o', linewidth=2.3, label=tf, color=colors.get(tf, None))
        else:
            print(f'Missing stage activity for {tf}')
    ax.set_xticks(x)
    ax.set_xticklabels(STAGES)
    ax.set_ylabel('Normalized composite out-activity')
    ax.set_title('Representative TF stage trajectories', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    savefig(fig, out_path)

def fig5B_representative_ego(edges, cent, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, tf in zip(axes, REPRESENTATIVE_TFS):
        plot_ego_network_for_tf(ax, edges, cent, tf, max_targets=18)
    fig.suptitle('Representative TF-centered regulatory modules', fontsize=15, fontweight='bold')
    savefig(fig, out_path)

def fig5C_novel_edge_summary(edges, cent, out_path):
    if 'is_novel' not in edges.columns or edges['is_novel'].isna().all():
        print('Skip Fig5C: novel status unavailable.')
        return
    top_core = set(get_top_core_tfs(cent, 20))
    edges = edges.copy()
    edges['core_related'] = edges['regulator_tf'].isin(top_core) | edges['target_tf'].isin(top_core)
    rows = []
    rows.append({'category': 'All high-confidence edges', 'n_edges': edges.shape[0]})
    rows.append({'category': 'Initial-supported edges', 'n_edges': int((edges['is_novel'] == False).sum())})
    rows.append({'category': 'Novel predicted edges', 'n_edges': int((edges['is_novel'] == True).sum())})
    rows.append({'category': 'Novel core-related edges', 'n_edges': int(((edges['is_novel'] == True) & edges['core_related']).sum())})
    rows.append({'category': 'PATZ1 novel outgoing edges', 'n_edges': int(((edges['is_novel'] == True) & (edges['regulator_tf'] == 'PATZ1')).sum())})
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLE_DIR / 'novel_edge_summary.csv', index=False, encoding='utf-8-sig')
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    x = np.arange(summary.shape[0])
    colors = ['#4C78A8', '#72B7B2', '#D1495B', '#F58518', '#6A4C93']
    ax.bar(x, summary['n_edges'], color=colors, edgecolor='white', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(summary['category'], rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Number of edges')
    ax.set_title('Summary of model-predicted novel regulatory edges', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    for i, v in enumerate(summary['n_edges']):
        ax.text(i, v, str(v), ha='center', va='bottom', fontsize=9)
    savefig(fig, out_path)

def export_gene_sets_for_enrichment(edges, cent, tf_stage):
    """Export target TF sets for external GO/KEGG enrichment with R/clusterProfiler."""
    ENRICH_DIR.mkdir(parents=True, exist_ok=True)
    gene_sets = {}
    top_core = get_top_core_tfs(cent, 20)
    core_targets = sorted(set(edges.loc[edges['regulator_tf'].isin(top_core), 'target_tf']))
    gene_sets['coreTF_downstream_targets'] = core_targets
    gene_sets['PATZ1_downstream_targets'] = sorted(set(edges.loc[edges['regulator_tf'] == 'PATZ1', 'target_tf']))
    for tf in REPRESENTATIVE_TFS:
        gene_sets[f'{tf}_downstream_targets'] = sorted(set(edges.loc[edges['regulator_tf'] == tf, 'target_tf']))
    if 'is_novel' in edges.columns and (not edges['is_novel'].isna().all()):
        gene_sets['novel_edge_targets'] = sorted(set(edges.loc[edges['is_novel'] == True, 'target_tf']))
        gene_sets['novel_core_related_targets'] = sorted(set(edges.loc[(edges['is_novel'] == True) & edges['regulator_tf'].isin(top_core), 'target_tf']))
    if not tf_stage.empty:
        tmp = tf_stage.copy()
        tmp['CRC_vs_Normal'] = tmp['CRC'] - tmp['Normal']
        crc_enhanced = tmp.sort_values('CRC_vs_Normal', ascending=False).head(50).index.tolist()
        gene_sets['CRC_enhanced_TFs'] = crc_enhanced
    rows = []
    for name, genes in gene_sets.items():
        genes = [g for g in genes if isinstance(g, str) and g not in ['', 'NAN', 'NONE']]
        (ENRICH_DIR / f'{name}.txt').write_text('\n'.join(genes), encoding='utf-8')
        rows.append({'gene_set': name, 'n_genes': len(genes), 'file': str(ENRICH_DIR / f'{name}.txt')})
    pd.DataFrame(rows).to_csv(ENRICH_DIR / 'gene_sets_summary.csv', index=False, encoding='utf-8-sig')
    print('Exported enrichment gene sets to:', ENRICH_DIR)
    return pd.DataFrame(rows)

def fig5D_enrichment_dotplot_or_placeholder(out_path):
    candidates = [FIG5_DIR / 'enrichment_results_combined.csv', ENRICH_DIR / 'enrichment_results_combined.csv']
    p = first_existing(candidates)
    if p is None:
        summary_path = ENRICH_DIR / 'gene_sets_summary.csv'
        if summary_path.exists():
            summary = pd.read_csv(summary_path)
        else:
            summary = pd.DataFrame({'gene_set': [], 'n_genes': []})
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
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
        savefig(fig, out_path)
        return
    enrich = pd.read_csv(p)
    desc_col = 'Description' if 'Description' in enrich.columns else 'term'
    padj_col = 'p.adjust' if 'p.adjust' in enrich.columns else 'padj'
    count_col = 'Count' if 'Count' in enrich.columns else 'count'
    set_col = 'gene_set' if 'gene_set' in enrich.columns else 'GeneSet'
    source_col = 'source' if 'source' in enrich.columns else 'ontology'
    enrich[padj_col] = pd.to_numeric(enrich[padj_col], errors='coerce')
    enrich[count_col] = pd.to_numeric(enrich[count_col], errors='coerce').fillna(1)
    enrich = enrich.dropna(subset=[padj_col])
    enrich = enrich[enrich[padj_col] <= 0.1].copy()
    if enrich.empty:
        print('No enrichment terms with adjusted p <= 0.10.')
        return
    enrich = enrich.sort_values([set_col, padj_col]).groupby(set_col).head(5).copy()
    enrich['label'] = enrich[set_col].astype(str) + ' | ' + enrich[desc_col].astype(str)
    enrich = enrich.sort_values(padj_col, ascending=False).tail(30)
    enrich['neglog10_padj'] = -np.log10(enrich[padj_col].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(9.5, max(5.5, 0.24 * enrich.shape[0])))
    y = np.arange(enrich.shape[0])
    sc = ax.scatter(enrich['neglog10_padj'], y, s=30 + 18 * enrich[count_col], c=enrich['neglog10_padj'], cmap='YlOrRd', edgecolor='white', linewidth=0.5, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(enrich['label'], fontsize=7)
    ax.set_xlabel('-log10(adjusted P value)')
    ax.set_title('Functional enrichment of model-derived TF sets', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.2)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('-log10(adj. P)', fontsize=9)
    savefig(fig, out_path)

def make_figure5(edges, cent, tf_stage, optional_activity):
    print('\n=== Figure 5 ===')
    fig5A_representative_trajectories(tf_stage, optional_activity, FIG5_DIR / 'Fig5A_KLF5_SP5_KLF4_stage_trajectory.png')
    fig5B_representative_ego(edges, cent, FIG5_DIR / 'Fig5B_representative_TF_ego_networks.png')
    fig5C_novel_edge_summary(edges, cent, FIG5_DIR / 'Fig5C_novel_edge_summary.png')
    export_gene_sets_for_enrichment(edges, cent, tf_stage)
    fig5D_enrichment_dotplot_or_placeholder(FIG5_DIR / 'Fig5D_GO_KEGG_enrichment_dotplot.png')

def main():
    print('=' * 100)
    print('Generate final result figures without cell-type annotation')
    print('=' * 100)
    print('OUT_ROOT:', OUT_ROOT)
    high_edges, edge_path = load_high_conf_edges()
    initial_set, initial_path = load_initial_edge_set()
    high_edges = add_initial_status(high_edges, initial_set)
    core_scores, core_path = load_core_score_table()
    cent = compute_centrality(high_edges, core_scores)
    stage_edges, stage_path = load_stage_edges()
    if stage_edges is not None:
        stage_edges = add_initial_status(stage_edges, initial_set)
    tf_stage = compute_stage_tf_activity(stage_edges)
    optional_activity, optional_path = load_optional_multimodal_activity()
    high_edges.to_csv(TABLE_DIR / 'high_confidence_edges_with_initial_status.csv', index=False, encoding='utf-8-sig')
    cent.to_csv(TABLE_DIR / 'TF_centrality_core_scores.csv', encoding='utf-8-sig')
    if not tf_stage.empty:
        tf_stage.to_csv(TABLE_DIR / 'TF_stage_out_activity_matrix.csv', encoding='utf-8-sig')
    make_figure2(high_edges, cent)
    make_figure3(stage_edges, core_scores, cent, tf_stage)
    make_figure4(high_edges, stage_edges, cent, tf_stage, optional_activity)
    make_figure5(high_edges, cent, tf_stage, optional_activity)
    readme = f'\nFinal result figures without cell-type annotation\n\nInputs:\n    high-confidence edges: {edge_path}\n    initial adjacency: {initial_path}\n    core scores: {core_path}\n    stage edge activity: {stage_path}\n    optional RNA/RP/composite activity table: {optional_path}\n\nMain output:\n    {OUT_ROOT}\n\nFigure 2:\n    {FIG2_DIR}\n\nFigure 3:\n    {FIG3_DIR}\n\nFigure 4:\n    {FIG4_DIR}\n\nFigure 5:\n    {FIG5_DIR}\n\nNotes:\n    - Cell-type annotation is not used.\n    - Stage-specific analyses are based on Normal/Polyp/CRC activity in the model-derived network.\n    - Enrichment gene sets are exported under:\n      {ENRICH_DIR}\n    - Run run_enrichment_clusterProfiler.R to generate GO/KEGG enrichment dotplot.\n'
    (OUT_ROOT / 'README_final_figures_no_celltype.txt').write_text(readme, encoding='utf-8')
    print('\nDONE.')
    print('Output root:', OUT_ROOT)
'''

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    d=load_source('network', SOURCE_network)
    edges,_=d.load_high_conf_edges()
    core,_=d.load_core_score_table()
    cent=d.compute_centrality(edges,core)
    stage,_=d.load_stage_edges();activity=d.compute_stage_tf_activity(stage)
    path=OUT/'A.png';attempt('A',lambda:(d.fig4D_patz1_target_heatmap(stage,activity,path),panel_file('A',path)))
    missing('B','相关分析有Enrichr点图函数，但论文为横向柱图，最终版原始绘图函数未定位。')
    missing('C','未找到论文生物学主题汇总柱图的原始绘图函数。')

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        missing("运行", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
    finally:
        finish(3)
