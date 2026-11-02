# -*- coding: utf-8 -*-
##### S1 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = 'S1'
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
from pathlib import Path as celltype_final_marker_utils_Path
import json as celltype_final_marker_utils_json
import re as celltype_final_marker_utils_re
import numpy as celltype_final_marker_utils_np
import pandas as celltype_final_marker_utils_pd
celltype_final_marker_utils_FINAL_MARKERS_JSON = None
celltype_final_marker_utils_MARKER_PATH_CANDIDATES = [celltype_final_marker_utils_Path(__file__).resolve().parent / 'final_markers.json', celltype_final_marker_utils_Path('E:\\毕业设计\\scmultiome_project\\第四次\\final_markers.json'), celltype_final_marker_utils_Path(str(DATA_DIR / 'reference' / 'final_markers.json')), celltype_final_marker_utils_Path(str(DATA_DIR / '4.0result' / 'ATACmotif' / '33_celltype_annotation' / 'final_markers.json'))]
celltype_final_marker_utils_GENE_ALIASES = {'CD45': 'PTPRC', 'CD31': 'PECAM1', 'CD144': 'CDH5', 'CD20': 'MS4A1', 'CD56': 'NCAM1', 'CD25': 'IL2RA', 'CD103': 'ITGAE', 'CD39': 'ENTPD1', 'CD73': 'NT5E', 'CD206': 'MRC1', 'CD169': 'SIGLEC1', 'UPAR': 'PLAUR', 'KI67': 'MKI67', 'KI-67': 'MKI67', 'ASMA': 'ACTA2', 'VIMENTIN': 'VIM', 'E-CADHERIN': 'CDH1', 'CLAUDIN-1': 'CLDN1', 'CK20': 'KRT20', 'PANCK': 'EPCAM', 'ESA': 'EPCAM', 'PD1': 'PDCD1', 'CHROMOGRANIN': 'CHGA', 'XBP1': 'XBP1', 'XPB1': 'XBP1', 'OCT3/4': 'POU5F1', 'OCT4': 'POU5F1', '4-OCT': 'POU5F1', 'NESTIN': 'NES', 'ALDEFLUOR': 'ALDH1A1', 'ALDH1': 'ALDH1A1', 'GPR49': 'LGR5', 'GLI-1': 'GLI1', 'MUSASHI-1': 'MSI1', 'MHC-II': 'HLA-DRA', 'HLA-DR': 'HLA-DRA', 'CD3DELTA': 'CD3D', 'CD3GAMA': 'CD3G'}
celltype_final_marker_utils_DROP_MARKERS = {'', 'NA', 'NONE', 'NAN', 'TOMATO', 'CLASS', 'FAMILY', 'II', 'MHC', 'PAN', 'MARKER'}

def celltype_final_marker_utils_norm_stage(x):
    s = str(x).strip().lower()
    if 'normal' in s or 'b001' in s or 'b004' in s:
        return 'Normal'
    if 'polyp' in s or 'a018' in s or 'a022' in s:
        return 'Polyp'
    if 'crc' in s or 'tumor' in s or 'cancer' in s or ('8810' in s) or ('15564' in s) or ('11773' in s):
        return 'CRC'
    if s in ['0', 'n']:
        return 'Normal'
    if s in ['1', 'p']:
        return 'Polyp'
    if s in ['2', 'c']:
        return 'CRC'
    return celltype_final_marker_utils_np.nan

def celltype_final_marker_utils_sparse_to_dense(X):
    try:
        from scipy import sparse as sparse
        if sparse.issparse(X):
            return X.toarray()
    except Exception:
        pass
    return celltype_final_marker_utils_np.asarray(X)

def celltype_final_marker_utils_minmax01(x):
    x = celltype_final_marker_utils_pd.Series(x, dtype=float).replace([celltype_final_marker_utils_np.inf, -celltype_final_marker_utils_np.inf], celltype_final_marker_utils_np.nan)
    if x.dropna().empty:
        return celltype_final_marker_utils_pd.Series(celltype_final_marker_utils_np.zeros(len(x)), index=x.index)
    x = x.fillna(x.median())
    if x.max() == x.min():
        return celltype_final_marker_utils_pd.Series(celltype_final_marker_utils_np.ones(len(x)) * 0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())
from pathlib import Path as celltype_Path
import warnings as celltype_warnings
import numpy as celltype_np
import pandas as celltype_pd
import matplotlib.pyplot as celltype_plt
import networkx as celltype_nx
celltype_warnings.filterwarnings('ignore')
celltype_BASE_DIR = celltype_Path(str(DATA_DIR))
celltype_ATACMOTIF_DIR = celltype_BASE_DIR / '4.0result' / 'ATACmotif'
celltype_OUT_DIR = celltype_ATACMOTIF_DIR / '40_coreTF_stage_celltype_like_gkf138_finalmarkers'
celltype_META_CANDIDATES = [celltype_ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'cell_metadata_with_final_marker_cell_type.csv', celltype_ATACMOTIF_DIR / '33_celltype_annotation' / 'cell_metadata_with_cell_type.csv']
celltype_H5AD_CANDIDATES = [celltype_ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'RNA_annotated_celltypes_final_markers.h5ad', celltype_ATACMOTIF_DIR / '33_celltype_annotation' / 'RNA_annotated_celltypes.h5ad']
celltype_EDGE_CANDIDATES = [celltype_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', celltype_ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'high_confidence_edges_multimodal.csv']
celltype_CORE_SCORE_CANDIDATES = [celltype_ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '05_stage_specific_TF_role_scores.csv', celltype_ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '04_coreTF_scores_recomputed.csv']
celltype_STAGES = ['Normal', 'Polyp', 'CRC']
celltype_MANUAL_FOCUS_TFS = ['PATZ1', 'KLF5', 'SP5', 'KLF4', 'ZNF148', 'ZNF281', 'EGR1', 'MAZ']
celltype_TOP_CORE_N = 15
celltype_MIN_CELLS_PER_GROUP = 30
celltype_EXPR_THRESHOLD = 0.0
celltype_MAX_TARGETS_PER_EGO = 24
celltype_MAX_FOCUS_TF_NETWORKS = 5
celltype_STAGE_COLORS = {'Normal': '#4C78A8', 'Polyp': '#F58518', 'CRC': '#E45756'}
celltype_CELLTYPE_COLORS = {'Epithelial': '#4C78A8', 'T_NK': '#54A24B', 'B_Plasma': '#72B7B2', 'Myeloid': '#F58518', 'Stromal': '#9D755D', 'Endothelial': '#B279A2', 'Mast_cell': '#ECA82C', 'Enteroendocrine': '#8CD17D', 'Unknown': '#B8B8B8'}
celltype_DEFAULT_COLORS = ['#4C78A8', '#F58518', '#54A24B', '#E45756', '#72B7B2', '#B279A2', '#9D755D', '#BAB0AC']

def celltype_first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None

def celltype_norm_gene(x):
    return str(x).strip().upper()

def celltype_savefig(fig, name):
    fig.savefig(celltype_OUT_DIR / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(celltype_OUT_DIR / f'{name}.pdf', bbox_inches='tight')
    celltype_plt.close(fig)

def celltype_load_meta_h5ad():
    import anndata as ad
    meta_path = celltype_first_existing(celltype_META_CANDIDATES)
    h5ad_path = celltype_first_existing(celltype_H5AD_CANDIDATES)
    if meta_path is None or h5ad_path is None:
        raise FileNotFoundError('Run annotate_celltypes_with_final_markers.py first.')
    meta = celltype_pd.read_csv(meta_path, index_col=0)
    meta.index = meta.index.astype(str)
    if 'cell_type_final' in meta.columns:
        meta['cell_type'] = meta['cell_type_final'].astype(str)
    if 'cell_type' not in meta.columns:
        raise ValueError('No cell_type / cell_type_final in metadata.')
    if 'stage' not in meta.columns:
        if 'stage_label' in meta.columns:
            meta['stage'] = meta['stage_label'].map(celltype_final_marker_utils_norm_stage)
        elif 'sample' in meta.columns:
            meta['stage'] = meta['sample'].map(celltype_final_marker_utils_norm_stage)
        else:
            meta['stage'] = [celltype_final_marker_utils_norm_stage(x) for x in meta.index]
    else:
        meta['stage'] = meta['stage'].map(celltype_final_marker_utils_norm_stage)
    meta = meta[meta['stage'].isin(celltype_STAGES)].copy()
    adata = ad.read_h5ad(h5ad_path)
    adata.obs.index = adata.obs.index.astype(str)
    adata.var_names = adata.var_names.astype(str)
    common = adata.obs.index.intersection(meta.index)
    adata = adata[common].copy()
    meta = meta.loc[common].copy()
    return (meta, adata)

def celltype_load_edges():
    p = celltype_first_existing(celltype_EDGE_CANDIDATES)
    if p is None:
        raise FileNotFoundError('No high-confidence edge file found.')
    edges = celltype_pd.read_csv(p)
    if 'regulator_tf' not in edges.columns:
        for c in ['source_tf', 'source', 'TF', 'regulator']:
            if c in edges.columns:
                edges['regulator_tf'] = edges[c]
                break
    if 'target_tf' not in edges.columns:
        for c in ['target', 'target_gene', 'gene']:
            if c in edges.columns:
                edges['target_tf'] = edges[c]
                break
    if 'mean_score' not in edges.columns:
        for c in ['score', 'y_prob', 'prob', 'weight']:
            if c in edges.columns:
                edges['mean_score'] = celltype_pd.to_numeric(edges[c], errors='coerce')
                break
    if 'mean_score' not in edges.columns:
        edges['mean_score'] = 1.0
    edges['regulator_tf'] = edges['regulator_tf'].map(celltype_norm_gene)
    edges['target_tf'] = edges['target_tf'].map(celltype_norm_gene)
    edges['mean_score'] = celltype_pd.to_numeric(edges['mean_score'], errors='coerce').fillna(1.0)
    edges = edges[edges['regulator_tf'] != edges['target_tf']].drop_duplicates(['regulator_tf', 'target_tf'])
    edges.to_csv(celltype_OUT_DIR / 'edges_used.csv', index=False)
    return edges

def celltype_load_core_tfs(edges):
    p = celltype_first_existing(celltype_CORE_SCORE_CANDIDATES)
    if p is not None:
        core = celltype_pd.read_csv(p, index_col=0)
        core.index = core.index.astype(str).str.upper().str.strip()
    else:
        core = celltype_pd.DataFrame(index=sorted(set(edges['regulator_tf']).union(edges['target_tf'])))
    if 'core_score' not in core.columns:
        out_degree = edges.groupby('regulator_tf').size()
        weighted_out = edges.groupby('regulator_tf')['mean_score'].sum()
        core['out_degree'] = core.index.map(out_degree).fillna(0)
        core['weighted_out_degree'] = core.index.map(weighted_out).fillna(0)
        core['core_score'] = 0.5 * celltype_final_marker_utils_minmax01(core['out_degree']) + 0.5 * celltype_final_marker_utils_minmax01(core['weighted_out_degree'])
    top = core.sort_values('core_score', ascending=False).head(celltype_TOP_CORE_N).index.tolist()
    tf_list = list(dict.fromkeys([celltype_norm_gene(x) for x in celltype_MANUAL_FOCUS_TFS + top]))
    return (tf_list, core)

def celltype_extract_expr_for_genes(adata, genes):
    var_upper = celltype_pd.Index([str(x).strip().upper() for x in adata.var_names])
    gene_to_orig = {}
    for orig, up in zip(adata.var_names, var_upper):
        if up not in gene_to_orig:
            gene_to_orig[up] = orig
    present = [g for g in genes if g in gene_to_orig]
    orig = [gene_to_orig[g] for g in present]
    X = celltype_final_marker_utils_sparse_to_dense(adata[:, orig].X).astype(float)
    if X.size and celltype_np.nanmax(X) > 30:
        X = celltype_np.log1p(X)
    return celltype_pd.DataFrame(X.T, index=present, columns=adata.obs.index)

def celltype_compute_stats(expr, meta, tf_list):
    rows = []
    for (stage, ct), submeta in meta.groupby(['stage', 'cell_type']):
        cells = submeta.index.intersection(expr.columns)
        if len(cells) < celltype_MIN_CELLS_PER_GROUP:
            continue
        for tf in tf_list:
            if tf not in expr.index:
                continue
            vals = expr.loc[tf, cells].astype(float).values
            rows.append({'tf': tf, 'stage': stage, 'cell_type': ct, 'group': f'{stage}|{ct}', 'n_cells': len(cells), 'mean_expr': float(celltype_np.mean(vals)), 'frac_expr': float(celltype_np.mean(vals > celltype_EXPR_THRESHOLD))})
    stats = celltype_pd.DataFrame(rows)
    stats.to_csv(celltype_OUT_DIR / 'coreTF_stage_celltype_expression_stats.csv', index=False)
    dominant = stats.sort_values(['tf', 'stage', 'frac_expr', 'mean_expr'], ascending=[True, True, False, False]).groupby(['tf', 'stage']).head(1)
    dominant.to_csv(celltype_OUT_DIR / 'coreTF_dominant_celltype_per_stage.csv', index=False)
    return (stats, dominant)

def celltype_edge_activity(edges, stats):
    frac = {(r.tf, r.stage, r.cell_type): r.frac_expr for r in stats.itertuples(index=False)}
    groups = stats[['stage', 'cell_type']].drop_duplicates()
    rows = []
    for g in groups.itertuples(index=False):
        for e in edges.itertuples(index=False):
            rf = frac.get((e.regulator_tf, g.stage, g.cell_type), 0.0)
            tf = frac.get((e.target_tf, g.stage, g.cell_type), 0.0)
            val = float(e.mean_score) * (0.7 * rf + 0.3 * tf)
            if val > 0:
                rows.append({'stage': g.stage, 'cell_type': g.cell_type, 'group': f'{g.stage}|{g.cell_type}', 'regulator_tf': e.regulator_tf, 'target_tf': e.target_tf, 'mean_score': e.mean_score, 'reg_frac_expr': rf, 'target_frac_expr': tf, 'stage_celltype_edge_activity': val})
    out = celltype_pd.DataFrame(rows)
    out.to_csv(celltype_OUT_DIR / 'stage_celltype_edge_activity.csv', index=False)
    deg = out.groupby(['group', 'regulator_tf'])['stage_celltype_edge_activity'].sum().reset_index().pivot(index='regulator_tf', columns='group', values='stage_celltype_edge_activity').fillna(0)
    deg.to_csv(celltype_OUT_DIR / 'stage_celltype_TF_out_activity_matrix.csv')
    if deg.shape[1] >= 2:
        deg.corr().to_csv(celltype_OUT_DIR / 'stage_celltype_GRN_similarity_correlation.csv')
    return (out, deg)

##### S1A #####
def celltype_plot_dotplot(stats, tf_list):
    group_info = stats[['group', 'stage', 'cell_type']].drop_duplicates()
    group_info['stage_order'] = group_info['stage'].map({s: i for i, s in enumerate(celltype_STAGES)})
    group_info = group_info.sort_values(['stage_order', 'cell_type'])
    groups = group_info['group'].tolist()
    tf_list = [x for x in tf_list if x in set(stats['tf'])]
    mean_mat = stats.pivot(index='tf', columns='group', values='mean_expr').reindex(tf_list)[groups].fillna(0)
    frac_mat = stats.pivot(index='tf', columns='group', values='frac_expr').reindex(tf_list)[groups].fillna(0)
    scaled = mean_mat.copy()
    for tf in scaled.index:
        scaled.loc[tf] = celltype_final_marker_utils_minmax01(scaled.loc[tf]).values
    fig, ax = celltype_plt.subplots(figsize=(max(11, 0.45 * len(groups)), max(6, 0.35 * len(tf_list))))
    for i, tf in enumerate(tf_list):
        for j, group in enumerate(groups):
            ax.scatter(j, i, s=20 + 650 * frac_mat.loc[tf, group], c=[[scaled.loc[tf, group], 0.3, 1 - scaled.loc[tf, group]]], alpha=0.78, edgecolor='white', linewidth=0.3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=50, ha='right', fontsize=8)
    ax.set_yticks(range(len(tf_list)))
    ax.set_yticklabels(tf_list, fontsize=8)
    ax.invert_yaxis()
    ax.set_title('Core TF expression by stage × final-marker cell type', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.12)
    celltype_savefig(fig, 'Fig_coreTF_stage_celltype_dotplot_finalmarkers')

##### S1B #####
def celltype_plot_similarity(deg):
    if deg.shape[1] < 2:
        return
    corr = deg.corr()
    fig, ax = celltype_plt.subplots(figsize=(max(7, 0.42 * corr.shape[1]), max(6, 0.42 * corr.shape[0])))
    im = ax.imshow(corr.values, cmap='magma', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(corr.shape[1]))
    ax.set_xticklabels(corr.columns, rotation=50, ha='right', fontsize=8)
    ax.set_yticks(range(corr.shape[0]))
    ax.set_yticklabels(corr.index, fontsize=8)
    ax.set_title('Stage × cell-type regulatory profile similarity', fontsize=12, fontweight='bold')
    fig.colorbar(im, ax=ax, label='Pearson correlation')
    celltype_savefig(fig, 'Fig_stage_celltype_GRN_similarity_finalmarkers')
celltype_OUT_DIR = OUT / 'celltype' / 'OUT_DIR'
celltype_OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta, adata = celltype_load_meta_h5ad()
    edges = celltype_load_edges()
    tf_list, _ = celltype_load_core_tfs(edges)
    all_tfs = sorted(set(tf_list).union(edges['regulator_tf']).union(edges['target_tf']))
    expr = celltype_extract_expr_for_genes(adata, all_tfs)
    tf_list = [x for x in tf_list if x in expr.index][:celltype_TOP_CORE_N]
    stats, _ = celltype_compute_stats(expr, meta, all_tfs)
    _, deg = celltype_edge_activity(edges, stats)
    attempt('A', lambda: (celltype_plot_dotplot(stats[stats['tf'].isin(tf_list)].copy(), tf_list), panel_file('A', celltype_OUT_DIR / 'Fig_coreTF_stage_celltype_dotplot_finalmarkers.png')))
    attempt('B', lambda: (celltype_plot_similarity(deg), panel_file('B', celltype_OUT_DIR / 'Fig_stage_celltype_GRN_similarity_finalmarkers.png')))
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
