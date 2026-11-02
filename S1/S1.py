# -*- coding: utf-8 -*-
"""论文S1：同图子图代码合并。独立文件，不读取原脚本。缺失处不生成替代图。"""
from pathlib import Path
import os, sys, types, traceback, json
# 只需修改此处，或设置环境变量 CRC_DATA_DIR / CRC_OUTPUT_DIR。
DATA_DIR = Path(os.environ.get("CRC_DATA_DIR", r"E:\毕业设计\bishedata"))
OUTPUT_ROOT = Path(os.environ.get("CRC_OUTPUT_DIR", r"C:\Users\WangXinyue\Desktop\Wang_CRC_ScientificReports_Clean_Revised_Package (1)\26.10.14\111"))
FIGURE_NAME = 'S1'
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

# 来源：celltype_final_marker_utils.py
SOURCE_celltype_final_marker_utils = r'''"""
Shared utilities for cell-type annotation using final_markers.json.

IMPORTANT:
When level="major", all fine cell-type labels mapped to the same major class
are merged, and all marker genes under those fine labels are used after:
    1. cleaning obvious encoding artifacts,
    2. converting common protein/antibody aliases to official gene symbols where possible,
    3. removing only clearly non-gene tokens.
It does NOT discard valid genes such as CD4, CD8A, CD14, CD34, CD44, etc.
"""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
FINAL_MARKERS_JSON = None
MARKER_PATH_CANDIDATES = [Path(__file__).resolve().parent / 'final_markers.json', Path('E:\\毕业设计\\scmultiome_project\\第四次\\final_markers.json'), Path('E:\\毕业设计\\bishedata\\reference\\final_markers.json'), Path('E:\\毕业设计\\bishedata\\4.0result\\ATACmotif\\33_celltype_annotation\\final_markers.json')]
GENE_ALIASES = {'CD45': 'PTPRC', 'CD31': 'PECAM1', 'CD144': 'CDH5', 'CD20': 'MS4A1', 'CD56': 'NCAM1', 'CD25': 'IL2RA', 'CD103': 'ITGAE', 'CD39': 'ENTPD1', 'CD73': 'NT5E', 'CD206': 'MRC1', 'CD169': 'SIGLEC1', 'UPAR': 'PLAUR', 'KI67': 'MKI67', 'KI-67': 'MKI67', 'ASMA': 'ACTA2', 'VIMENTIN': 'VIM', 'E-CADHERIN': 'CDH1', 'CLAUDIN-1': 'CLDN1', 'CK20': 'KRT20', 'PANCK': 'EPCAM', 'ESA': 'EPCAM', 'PD1': 'PDCD1', 'CHROMOGRANIN': 'CHGA', 'XBP1': 'XBP1', 'XPB1': 'XBP1', 'OCT3/4': 'POU5F1', 'OCT4': 'POU5F1', '4-OCT': 'POU5F1', 'NESTIN': 'NES', 'ALDEFLUOR': 'ALDH1A1', 'ALDH1': 'ALDH1A1', 'GPR49': 'LGR5', 'GLI-1': 'GLI1', 'MUSASHI-1': 'MSI1', 'MHC-II': 'HLA-DRA', 'HLA-DR': 'HLA-DRA', 'CD3DELTA': 'CD3D', 'CD3GAMA': 'CD3G'}
DROP_MARKERS = {'', 'NA', 'NONE', 'NAN', 'TOMATO', 'CLASS', 'FAMILY', 'II', 'MHC', 'PAN', 'MARKER'}

def find_marker_path(path=None):
    if path is not None:
        p = Path(path)
        if p.exists():
            return p
    if FINAL_MARKERS_JSON is not None:
        p = Path(FINAL_MARKERS_JSON)
        if p.exists():
            return p
    for p in MARKER_PATH_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError('Cannot find final_markers.json. Put it in the same directory as the scripts or set FINAL_MARKERS_JSON / MARKER_PATH_CANDIDATES.')

def clean_gene_symbol(g):
    s = str(g).strip().upper()
    s = s.replace('ÃŽÂ', '')
    s = s.replace('Â', '')
    s = s.replace('ÃŽÂ', '')
    s = s.replace('GAMA', 'G')
    s = s.replace('DELTA', 'D')
    s = s.strip()
    s = GENE_ALIASES.get(s, s)
    s = re.sub('[^A-Z0-9\\-\\/]', '', s)
    s = GENE_ALIASES.get(s, s)
    if s in DROP_MARKERS:
        return ''
    return s

def fine_to_major_celltype(label):
    x = str(label).lower()
    if 'mast cell' in x or 'eosinophil' in x:
        return 'Mast_cell'
    if any((k in x for k in ['endothelial', 'vein', 'artery', 'capillary', 'vascular', 'lymphatic vessel', 'tip cell', 'tumor endothelial'])):
        return 'Endothelial'
    if any((k in x for k in ['fibroblast', 'stromal', 'mesenchymal', 'myofibroblast', 'pericyte', 'smooth muscle', 'muscle', 'adventitial', 'cajal', 'mesothelial'])):
        return 'Stromal'
    if any((k in x for k in ['macrophage', 'monocyte', 'myeloid', 'dendritic', 'neutrophil', 'granulocyte', 'myelocyte', 'phagocyte'])):
        return 'Myeloid'
    if any((k in x for k in ['b cell', 'plasma', 'plasmablast', 'follicular', 'germinal center', 'precursor b', 'b-1 b'])):
        return 'B_Plasma'
    if any((k in x for k in ['t cell', 'natural killer', 'nk cell', ' nk ', 'mait', 'lymphoid', 'alpha-beta', 'gamma-delta', 'cytotoxic', 'helper', 'treg', 'tfh', 'th17', 'th2', 'th0', 'trm', 'nkt', 'innate lymphoid'])):
        return 'T_NK'
    if any((k in x for k in ['enteroendocrine', 'endocrine', 'chromaffin', 'type a entero', 'type d entero', 'type i entero', 'type l entero', 'type m entero', 'type n entero'])):
        return 'Enteroendocrine'
    if any((k in x for k in ['neuron', 'glial'])):
        return 'Neural_glial'
    if any((k in x for k in ['epithelial', 'enterocyte', 'goblet', 'coloncyte', 'colonocyte', 'paneth', 'tuft', 'brush cell', 'stem cell', 'progenitor', 'transit', 'cancer cell', 'colorectal stem', 'intestinal crypt', 'columnar', 'simple columnar'])):
        return 'Epithelial'
    if 'immune cell' in x:
        return 'Immune'
    return str(label).replace(' ', '_')

def load_marker_dict(marker_path=None, level='major', min_markers=1):
    """
    level:
        "major" = merge all fine labels into broad, stable cell types.
        "fine"  = keep original labels from final_markers.json.

    Returns:
        marker_dict, marker_path
    """
    p = find_marker_path(marker_path)
    with open(p, 'r', encoding='utf-8') as f:
        raw = json.load(f)
    out = {}
    mapping_rows = []
    for label, genes in raw.items():
        major = fine_to_major_celltype(label)
        ct = major if level == 'major' else str(label).replace(' ', '_')
        cleaned = []
        removed = []
        for g in genes:
            cg = clean_gene_symbol(g)
            if cg:
                cleaned.append(cg)
            else:
                removed.append(str(g))
        if not cleaned:
            mapping_rows.append({'source_label': label, 'major_label': major, 'used_label': ct, 'n_raw_markers': len(genes), 'n_cleaned_markers': 0, 'n_removed_tokens': len(removed), 'removed_tokens': ','.join(removed)})
            continue
        out.setdefault(ct, [])
        out[ct].extend(cleaned)
        mapping_rows.append({'source_label': label, 'major_label': major, 'used_label': ct, 'n_raw_markers': len(genes), 'n_cleaned_markers': len(cleaned), 'n_removed_tokens': len(removed), 'removed_tokens': ','.join(removed)})
    final = {}
    for ct, genes in out.items():
        uniq = list(dict.fromkeys(genes))
        if len(uniq) >= min_markers:
            final[ct] = uniq
    try:
        audit_path = p.parent / f'final_markers_{level}_merge_audit.csv'
        pd.DataFrame(mapping_rows).to_csv(audit_path, index=False, encoding='utf-8-sig')
    except Exception:
        pass
    return (final, p)

def norm_stage(x):
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
    return np.nan

def sparse_to_dense(X):
    try:
        from scipy import sparse
        if sparse.issparse(X):
            return X.toarray()
    except Exception:
        pass
    return np.asarray(X)

def minmax01(x):
    x = pd.Series(x, dtype=float).replace([np.inf, -np.inf], np.nan)
    if x.dropna().empty:
        return pd.Series(np.zeros(len(x)), index=x.index)
    x = x.fillna(x.median())
    if x.max() == x.min():
        return pd.Series(np.ones(len(x)) * 0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())
'''

# 来源：coreTF_stage_celltype_like_gkf138_finalmarkers.py
SOURCE_celltype = r'''"""
Core TF stage × cell-type analysis using final_markers.json-based annotation.

Run after annotate_celltypes_with_final_markers.py:
    D:\\Python\\python.exe E:\\毕业设计\\scmultiome_project\\第四次\\coreTF_stage_celltype_like_gkf138_finalmarkers.py
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
warnings.filterwarnings('ignore')
from celltype_final_marker_utils import norm_stage, sparse_to_dense, minmax01
BASE_DIR = Path('E:\\毕业设计\\bishedata')
ATACMOTIF_DIR = BASE_DIR / '4.0result' / 'ATACmotif'
OUT_DIR = ATACMOTIF_DIR / '40_coreTF_stage_celltype_like_gkf138_finalmarkers'
META_CANDIDATES = [ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'cell_metadata_with_final_marker_cell_type.csv', ATACMOTIF_DIR / '33_celltype_annotation' / 'cell_metadata_with_cell_type.csv']
H5AD_CANDIDATES = [ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'RNA_annotated_celltypes_final_markers.h5ad', ATACMOTIF_DIR / '33_celltype_annotation' / 'RNA_annotated_celltypes.h5ad']
EDGE_CANDIDATES = [ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '04_stage_specific_edge_activity.csv', ATACMOTIF_DIR / '27_scMultiomeGRN_paper_faithful_optimized_valauc_quantile' / 'high_confidence_edges_multimodal.csv']
CORE_SCORE_CANDIDATES = [ATACMOTIF_DIR / '38_stage_specific_TF_TF_networks' / '05_stage_specific_TF_role_scores.csv', ATACMOTIF_DIR / '31_sections_3_3_to_3_6_downstream_results' / '04_coreTF_scores_recomputed.csv']
STAGES = ['Normal', 'Polyp', 'CRC']
MANUAL_FOCUS_TFS = ['PATZ1', 'KLF5', 'SP5', 'KLF4', 'ZNF148', 'ZNF281', 'EGR1', 'MAZ']
TOP_CORE_N = 15
MIN_CELLS_PER_GROUP = 30
EXPR_THRESHOLD = 0.0
MAX_TARGETS_PER_EGO = 24
MAX_FOCUS_TF_NETWORKS = 5
STAGE_COLORS = {'Normal': '#4C78A8', 'Polyp': '#F58518', 'CRC': '#E45756'}
CELLTYPE_COLORS = {'Epithelial': '#4C78A8', 'T_NK': '#54A24B', 'B_Plasma': '#72B7B2', 'Myeloid': '#F58518', 'Stromal': '#9D755D', 'Endothelial': '#B279A2', 'Mast_cell': '#ECA82C', 'Enteroendocrine': '#8CD17D', 'Unknown': '#B8B8B8'}
DEFAULT_COLORS = ['#4C78A8', '#F58518', '#54A24B', '#E45756', '#72B7B2', '#B279A2', '#9D755D', '#BAB0AC']

def first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None

def norm_gene(x):
    return str(x).strip().upper()

def get_colors(labels):
    labels = sorted(pd.Series(labels).dropna().astype(str).unique())
    out, j = ({}, 0)
    for lab in labels:
        if lab in CELLTYPE_COLORS:
            out[lab] = CELLTYPE_COLORS[lab]
        else:
            out[lab] = DEFAULT_COLORS[j % len(DEFAULT_COLORS)]
            j += 1
    return out

def savefig(fig, name):
    fig.savefig(OUT_DIR / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)

def load_meta_h5ad():
    import anndata as ad
    meta_path = first_existing(META_CANDIDATES)
    h5ad_path = first_existing(H5AD_CANDIDATES)
    if meta_path is None or h5ad_path is None:
        raise FileNotFoundError('Run annotate_celltypes_with_final_markers.py first.')
    meta = pd.read_csv(meta_path, index_col=0)
    meta.index = meta.index.astype(str)
    if 'cell_type_final' in meta.columns:
        meta['cell_type'] = meta['cell_type_final'].astype(str)
    if 'cell_type' not in meta.columns:
        raise ValueError('No cell_type / cell_type_final in metadata.')
    if 'stage' not in meta.columns:
        if 'stage_label' in meta.columns:
            meta['stage'] = meta['stage_label'].map(norm_stage)
        elif 'sample' in meta.columns:
            meta['stage'] = meta['sample'].map(norm_stage)
        else:
            meta['stage'] = [norm_stage(x) for x in meta.index]
    else:
        meta['stage'] = meta['stage'].map(norm_stage)
    meta = meta[meta['stage'].isin(STAGES)].copy()
    adata = ad.read_h5ad(h5ad_path)
    adata.obs.index = adata.obs.index.astype(str)
    adata.var_names = adata.var_names.astype(str)
    common = adata.obs.index.intersection(meta.index)
    adata = adata[common].copy()
    meta = meta.loc[common].copy()
    return (meta, adata)

def load_edges():
    p = first_existing(EDGE_CANDIDATES)
    if p is None:
        raise FileNotFoundError('No high-confidence edge file found.')
    edges = pd.read_csv(p)
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
                edges['mean_score'] = pd.to_numeric(edges[c], errors='coerce')
                break
    if 'mean_score' not in edges.columns:
        edges['mean_score'] = 1.0
    edges['regulator_tf'] = edges['regulator_tf'].map(norm_gene)
    edges['target_tf'] = edges['target_tf'].map(norm_gene)
    edges['mean_score'] = pd.to_numeric(edges['mean_score'], errors='coerce').fillna(1.0)
    edges = edges[edges['regulator_tf'] != edges['target_tf']].drop_duplicates(['regulator_tf', 'target_tf'])
    edges.to_csv(OUT_DIR / 'edges_used.csv', index=False)
    return edges

def load_core_tfs(edges):
    p = first_existing(CORE_SCORE_CANDIDATES)
    if p is not None:
        core = pd.read_csv(p, index_col=0)
        core.index = core.index.astype(str).str.upper().str.strip()
    else:
        core = pd.DataFrame(index=sorted(set(edges['regulator_tf']).union(edges['target_tf'])))
    if 'core_score' not in core.columns:
        out_degree = edges.groupby('regulator_tf').size()
        weighted_out = edges.groupby('regulator_tf')['mean_score'].sum()
        core['out_degree'] = core.index.map(out_degree).fillna(0)
        core['weighted_out_degree'] = core.index.map(weighted_out).fillna(0)
        core['core_score'] = 0.5 * minmax01(core['out_degree']) + 0.5 * minmax01(core['weighted_out_degree'])
    top = core.sort_values('core_score', ascending=False).head(TOP_CORE_N).index.tolist()
    tf_list = list(dict.fromkeys([norm_gene(x) for x in MANUAL_FOCUS_TFS + top]))
    return (tf_list, core)

def extract_expr_for_genes(adata, genes):
    var_upper = pd.Index([str(x).strip().upper() for x in adata.var_names])
    gene_to_orig = {}
    for orig, up in zip(adata.var_names, var_upper):
        if up not in gene_to_orig:
            gene_to_orig[up] = orig
    present = [g for g in genes if g in gene_to_orig]
    orig = [gene_to_orig[g] for g in present]
    X = sparse_to_dense(adata[:, orig].X).astype(float)
    if X.size and np.nanmax(X) > 30:
        X = np.log1p(X)
    return pd.DataFrame(X.T, index=present, columns=adata.obs.index)

def compute_stats(expr, meta, tf_list):
    rows = []
    for (stage, ct), submeta in meta.groupby(['stage', 'cell_type']):
        cells = submeta.index.intersection(expr.columns)
        if len(cells) < MIN_CELLS_PER_GROUP:
            continue
        for tf in tf_list:
            if tf not in expr.index:
                continue
            vals = expr.loc[tf, cells].astype(float).values
            rows.append({'tf': tf, 'stage': stage, 'cell_type': ct, 'group': f'{stage}|{ct}', 'n_cells': len(cells), 'mean_expr': float(np.mean(vals)), 'frac_expr': float(np.mean(vals > EXPR_THRESHOLD))})
    stats = pd.DataFrame(rows)
    stats.to_csv(OUT_DIR / 'coreTF_stage_celltype_expression_stats.csv', index=False)
    dominant = stats.sort_values(['tf', 'stage', 'frac_expr', 'mean_expr'], ascending=[True, True, False, False]).groupby(['tf', 'stage']).head(1)
    dominant.to_csv(OUT_DIR / 'coreTF_dominant_celltype_per_stage.csv', index=False)
    return (stats, dominant)

def edge_activity(edges, stats):
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
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / 'stage_celltype_edge_activity.csv', index=False)
    deg = out.groupby(['group', 'regulator_tf'])['stage_celltype_edge_activity'].sum().reset_index().pivot(index='regulator_tf', columns='group', values='stage_celltype_edge_activity').fillna(0)
    deg.to_csv(OUT_DIR / 'stage_celltype_TF_out_activity_matrix.csv')
    if deg.shape[1] >= 2:
        deg.corr().to_csv(OUT_DIR / 'stage_celltype_GRN_similarity_correlation.csv')
    return (out, deg)

def plot_dotplot(stats, tf_list):
    group_info = stats[['group', 'stage', 'cell_type']].drop_duplicates()
    group_info['stage_order'] = group_info['stage'].map({s: i for i, s in enumerate(STAGES)})
    group_info = group_info.sort_values(['stage_order', 'cell_type'])
    groups = group_info['group'].tolist()
    tf_list = [x for x in tf_list if x in set(stats['tf'])]
    mean_mat = stats.pivot(index='tf', columns='group', values='mean_expr').reindex(tf_list)[groups].fillna(0)
    frac_mat = stats.pivot(index='tf', columns='group', values='frac_expr').reindex(tf_list)[groups].fillna(0)
    scaled = mean_mat.copy()
    for tf in scaled.index:
        scaled.loc[tf] = minmax01(scaled.loc[tf]).values
    fig, ax = plt.subplots(figsize=(max(11, 0.45 * len(groups)), max(6, 0.35 * len(tf_list))))
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
    savefig(fig, 'Fig_coreTF_stage_celltype_dotplot_finalmarkers')

def plot_dominant(dominant, tf_list):
    dom = dominant[dominant['tf'].isin(tf_list)].copy()
    tfs = [x for x in tf_list if x in set(dom['tf'])]
    mat_ct = dom.pivot(index='tf', columns='stage', values='cell_type').reindex(tfs)[STAGES]
    mat_frac = dom.pivot(index='tf', columns='stage', values='frac_expr').reindex(tfs)[STAGES]
    all_ct = sorted(pd.unique(mat_ct.values.ravel()))
    all_ct = [x for x in all_ct if pd.notna(x)]
    colors = get_colors(all_ct)
    import matplotlib.colors as mcolors
    img = np.zeros((len(tfs), len(STAGES), 3))
    for i, tf in enumerate(tfs):
        for j, st in enumerate(STAGES):
            img[i, j] = mcolors.to_rgb(colors.get(mat_ct.loc[tf, st], '#B8B8B8'))
    fig, ax = plt.subplots(figsize=(6.5, max(5, 0.35 * len(tfs))))
    ax.imshow(img, aspect='auto')
    ax.set_xticks(range(len(STAGES)))
    ax.set_xticklabels(STAGES)
    ax.set_yticks(range(len(tfs)))
    ax.set_yticklabels(tfs, fontsize=8)
    for i, tf in enumerate(tfs):
        for j, st in enumerate(STAGES):
            ct = mat_ct.loc[tf, st]
            fr = mat_frac.loc[tf, st]
            ax.text(j, i, f'{ct}\n{fr:.0%}' if pd.notna(ct) else 'NA', ha='center', va='center', fontsize=7)
    ax.set_title('Dominant cell type of core TFs in each stage', fontsize=13, fontweight='bold')
    handles = [plt.Line2D([0], [0], marker='s', color='w', label=ct, markerfacecolor=colors[ct], markersize=8) for ct in all_ct]
    ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    savefig(fig, 'Fig_coreTF_dominant_celltype_per_stage_finalmarkers')

def plot_similarity(deg):
    if deg.shape[1] < 2:
        return
    corr = deg.corr()
    fig, ax = plt.subplots(figsize=(max(7, 0.42 * corr.shape[1]), max(6, 0.42 * corr.shape[0])))
    im = ax.imshow(corr.values, cmap='magma', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(corr.shape[1]))
    ax.set_xticklabels(corr.columns, rotation=50, ha='right', fontsize=8)
    ax.set_yticks(range(corr.shape[0]))
    ax.set_yticklabels(corr.index, fontsize=8)
    ax.set_title('Stage × cell-type regulatory profile similarity', fontsize=12, fontweight='bold')
    fig.colorbar(im, ax=ax, label='Pearson correlation')
    savefig(fig, 'Fig_stage_celltype_GRN_similarity_finalmarkers')

def main():
    meta, adata = load_meta_h5ad()
    edges = load_edges()
    tf_list, core = load_core_tfs(edges)
    all_tfs = sorted(set(tf_list).union(edges['regulator_tf']).union(edges['target_tf']))
    expr = extract_expr_for_genes(adata, all_tfs)
    tf_list = [x for x in tf_list if x in expr.index][:TOP_CORE_N]
    stats, dominant = compute_stats(expr, meta, all_tfs)
    core_stats = stats[stats['tf'].isin(tf_list)].copy()
    core_stats.to_csv(OUT_DIR / 'coreTF_only_stage_celltype_expression_stats.csv', index=False)
    edges_ct, deg = edge_activity(edges, stats)
    plot_dotplot(core_stats, tf_list)
    plot_dominant(dominant, tf_list)
    plot_similarity(deg)
    print('DONE')
    print(OUT_DIR)
'''

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    load_source('celltype_final_marker_utils',SOURCE_celltype_final_marker_utils)
    c=load_source('celltype',SOURCE_celltype)
    meta,adata=c.load_meta_h5ad();edges=c.load_edges();tf_list,_=c.load_core_tfs(edges)
    all_tfs=sorted(set(tf_list).union(edges['regulator_tf']).union(edges['target_tf']))
    expr=c.extract_expr_for_genes(adata,all_tfs);tf_list=[x for x in tf_list if x in expr.index][:c.TOP_CORE_N]
    stats,_=c.compute_stats(expr,meta,all_tfs);_,deg=c.edge_activity(edges,stats)
    attempt('A',lambda:(c.plot_dotplot(stats[stats['tf'].isin(tf_list)].copy(),tf_list),panel_file('A',c.OUT_DIR/'Fig_coreTF_stage_celltype_dotplot_finalmarkers.png')))
    attempt('B',lambda:(c.plot_similarity(deg),panel_file('B',c.OUT_DIR/'Fig_stage_celltype_GRN_similarity_finalmarkers.png')))

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        missing("运行", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
    finally:
        finish(2)
