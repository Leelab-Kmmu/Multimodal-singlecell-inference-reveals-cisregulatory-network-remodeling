# -*- coding: utf-8 -*-
"""论文图2"""
from pathlib import Path
import os, sys, types, traceback, json
DATA_DIR = Path(os.environ.get("CRC_DATA_DIR", r"E:\毕业设计\bishedata"))
OUTPUT_ROOT = Path(os.environ.get("CRC_OUTPUT_DIR", r"C:\Users\WangXinyue\Desktop\Wang_CRC_ScientificReports_Clean_Revised_Package (1)\26.10.14\111"))
FIGURE_NAME = '图2'
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
    if ERRORS: print('部分面板未完成')

# 生成 Figure 1：
    A. UMAP by stage
    B. UMAP by broad lineage / major cell type, with Polyp-enriched state highlighted
    C. Fraction of Polyp-enriched transition state across Normal / Polyp / CRC
    D. Literature-inspired transition-state signature score dotplot
    E. Literature-derived key-gene heatmap

运行：
    D:\\Python\\python.exe E:\\毕业设计\\scmultiome_project\\第四次\\plot_result1_overall_polyp_transition.py

依赖：
    D:\\Python\\python.exe -m pip install pandas numpy matplotlib scipy anndata scanpy

输出目录：
    E:\\毕业设计\x08ishedata\x04.0result\\ATACmotif(_result1_overall_polyp_transition
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
BASE_DIR = Path('E:\\毕业设计\\bishedata')
ATACMOTIF_DIR = BASE_DIR / '4.0result' / 'ATACmotif'
OUT_DIR = ATACMOTIF_DIR / '50_result1_overall_polyp_transition'
H5AD_CANDIDATES = [ATACMOTIF_DIR / '49_polyp_transition_signature_becker2022' / 'RNA_with_polyp_transition_state_becker2022.h5ad', ATACMOTIF_DIR / '48_polyp_enriched_unclassified_state' / 'RNA_with_polyp_enriched_unclassified_state.h5ad', ATACMOTIF_DIR / '46_layered_cell_annotation' / '01_major_lineage_all_data' / 'RNA_major_lineage_annotated.h5ad', ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'RNA_annotated_celltypes_final_markers.h5ad', ATACMOTIF_DIR / '33_celltype_annotation' / 'RNA_annotated_celltypes.h5ad', BASE_DIR / '4.0版本数据' / 'GRNBoost2_input_FINAL' / 'GRNBoost2_ready.h5ad']
META_CANDIDATES = [ATACMOTIF_DIR / '49_polyp_transition_signature_becker2022' / 'cell_metadata_with_polyp_transition_state.csv', ATACMOTIF_DIR / '48_polyp_enriched_unclassified_state' / 'cell_metadata_with_polyp_enriched_state.csv', ATACMOTIF_DIR / '46_layered_cell_annotation' / '01_major_lineage_all_data' / 'cell_metadata_major_lineage.csv', ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'cell_metadata_with_final_marker_cell_type.csv', ATACMOTIF_DIR / '43_celltype_unknown_diagnosis_reannotation' / 'cell_metadata_with_cell_type_reannotated.csv', ATACMOTIF_DIR / '33_celltype_annotation' / 'cell_metadata_with_cell_type.csv']
STAGES = ['Normal', 'Polyp', 'CRC']
UNKNOWN_LABELS = {'unknown', 'unclassified', 'unassigned', 'other', 'nan', 'none', 'low_confidence', 'low_confidence_unknown'}
MIN_POLYP_UNKNOWN_CELLS = 20
POINT_SIZE = 4
RANDOM_SEED = 42
COMPUTE_UMAP_IF_MISSING = True
SIGNATURES = {'Stem-like epithelial': ['LGR5', 'SMOC2', 'RGMB', 'PTPRO', 'EPHB2', 'LRIG1', 'ASCL2', 'OLFM4', 'SOX9'], 'Malignancy continuum': ['OLFM4', 'GPX2', 'ASCL2', 'HNF4A', 'SOX4', 'KRT18', 'NORAD'], 'Mature enterocyte': ['CA1', 'CA2', 'AQP8', 'SLC26A3', 'GUCA2A', 'GUCA2B', 'RAB6B', 'BEST4', 'CA7', 'OTOP2', 'OTOP3', 'SLC4A4', 'NR5A2'], 'Goblet / immature goblet': ['MUC2', 'TFF1', 'TFF3', 'FCGBP', 'TBX10', 'KLK1', 'ITLN1', 'WFDC2', 'CLCA1', 'AGR2'], 'Proliferative TA': ['MKI67', 'TOP2A', 'PCNA', 'STMN1', 'UBE2C', 'CENPF', 'HMGB2', 'MCM2', 'MCM5'], 'preCAF-like': ['FAP', 'VCAN', 'COL1A2', 'COL11A1', 'THBS2', 'TWIST1', 'WNT2', 'RUNX1', 'SERPINE1', 'INHBA', 'PLXDC2', 'BICC1', 'FBN1'], 'Treg-like': ['FOXP3', 'CTLA4', 'IL2RA', 'TNFRSF4', 'LAIR2', 'IKZF2'], 'CRC-like epithelial': ['GPX2', 'HNF4A', 'CEACAM5', 'CEACAM6', 'CLDN1', 'MMP7', 'KRT19', 'TACSTD2'], 'General epithelial': ['EPCAM', 'TACSTD2', 'KRT8', 'KRT18', 'KRT19', 'KRT20', 'CDH1'], 'Stromal': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'ACTA2', 'TAGLN', 'VIM', 'PDGFRA'], 'Myeloid': ['LYZ', 'LST1', 'FCGR3A', 'S100A8', 'S100A9', 'C1QA', 'C1QB', 'CD68'], 'T/NK': ['PTPRC', 'CD3D', 'CD3E', 'TRAC', 'IL7R', 'NKG7', 'GNLY', 'PRF1']}
KEY_GENES_FOR_HEATMAP = ['LGR5', 'SMOC2', 'ASCL2', 'OLFM4', 'SOX9', 'GPX2', 'HNF4A', 'CA1', 'CA2', 'AQP8', 'SLC26A3', 'BEST4', 'MUC2', 'TFF3', 'FCGBP', 'CLCA1', 'MKI67', 'TOP2A', 'PCNA', 'FAP', 'VCAN', 'COL1A2', 'COL11A1', 'WNT2', 'RUNX1', 'FOXP3', 'CTLA4', 'IL2RA', 'EPCAM', 'KRT8', 'KRT18', 'KRT19']
STAGE_COLORS = {'Normal': '#4C78A8', 'Polyp': '#F58518', 'CRC': '#E45756'}
LINEAGE_COLORS = {'Epithelial': '#4C78A8', 'T_NK': '#54A24B', 'T/NK': '#54A24B', 'B_Plasma': '#72B7B2', 'B/Plasma': '#72B7B2', 'Myeloid': '#F58518', 'Stromal': '#9D755D', 'Endothelial': '#B279A2', 'Mast_cell': '#ECA82C', 'Enteroendocrine': '#8CD17D', 'Unknown': '#B8B8B8', 'Unclassified': '#B8B8B8', 'Polyp-enriched transition': '#D1495B'}
DEFAULT_COLORS = ['#4C78A8', '#F58518', '#54A24B', '#E45756', '#72B7B2', '#B279A2', '#9D755D', '#BAB0AC', '#6A4C93', '#2A9D8F', '#ECA82C', '#8CD17D']

def first_existing(paths):
    for p in paths:
        p = Path(p)
        if p.exists():
            return p
    return None

def norm_gene(x):
    return str(x).strip().upper()

def norm_stage(x):
    s = str(x).strip().lower()
    if 'normal' in s or 'b001' in s or 'b004' in s:
        return 'Normal'
    if 'polyp' in s or 'adenoma' in s or 'a018' in s or ('a022' in s):
        return 'Polyp'
    if 'crc' in s or 'carcinoma' in s or 'tumor' in s or ('cancer' in s) or ('8810' in s) or ('15564' in s) or ('11773' in s):
        return 'CRC'
    if s in ['0', 'n']:
        return 'Normal'
    if s in ['1', 'p', 'a']:
        return 'Polyp'
    if s in ['2', 'c']:
        return 'CRC'
    return np.nan

def is_unknown_label(x):
    s = str(x).strip().lower()
    return s in UNKNOWN_LABELS or 'unknown' in s or 'unclassified' in s or ('unassigned' in s)

def sparse_to_dense(X):
    try:
        from scipy import sparse
        if sparse.issparse(X):
            return X.toarray()
    except Exception:
        pass
    return np.asarray(X)

def safe_X(adata):
    try:
        from scipy import sparse
        if sparse.issparse(adata.X):
            X = adata.X.copy()
            X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
            return X
    except Exception:
        pass
    return np.nan_to_num(np.asarray(adata.X, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)

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

def get_palette(labels, base_colors=None):
    labels = list(pd.Series(labels).dropna().astype(str).unique())
    colors = {}
    j = 0
    for lab in labels:
        if base_colors and lab in base_colors:
            colors[lab] = base_colors[lab]
        elif lab in LINEAGE_COLORS:
            colors[lab] = LINEAGE_COLORS[lab]
        else:
            colors[lab] = DEFAULT_COLORS[j % len(DEFAULT_COLORS)]
            j += 1
    return colors

def savefig(fig, name):
    png = OUT_DIR / f'{name}.png'
    pdf = OUT_DIR / f'{name}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    print('Saved:', png)

def load_data():
    import anndata as ad
    h5ad_path = first_existing(H5AD_CANDIDATES)
    if h5ad_path is None:
        raise FileNotFoundError('Cannot find h5ad. Please update H5AD_CANDIDATES.')
    print('Using h5ad:', h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    adata.obs.index = adata.obs.index.astype(str)
    adata.var_names = adata.var_names.astype(str)
    meta_path = first_existing(META_CANDIDATES)
    if meta_path is not None:
        print('Using metadata:', meta_path)
        meta = pd.read_csv(meta_path, index_col=0)
        meta.index = meta.index.astype(str)
        common = adata.obs.index.intersection(meta.index)
        print('Metadata overlap:', len(common))
        if len(common) > 0:
            adata = adata[common].copy()
            meta = meta.loc[common].copy()
            for col in meta.columns:
                adata.obs[col] = meta[col].astype(str).values
    obs = adata.obs.copy()
    if 'stage' in obs.columns:
        obs['stage'] = obs['stage'].map(norm_stage)
    elif 'stage_label' in obs.columns:
        obs['stage'] = obs['stage_label'].map(norm_stage)
    elif 'sample' in obs.columns:
        obs['stage'] = obs['sample'].map(norm_stage)
    elif 'sample_id' in obs.columns:
        obs['stage'] = obs['sample_id'].map(norm_stage)
    else:
        obs['stage'] = [norm_stage(x) for x in obs.index]
    lineage_col = None
    for c in ['cell_type_major', 'major_lineage', 'cell_type_final', 'cell_type_reannotated', 'cell_type', 'annotation', 'celltype', 'major_cell_type']:
        if c in obs.columns:
            lineage_col = c
            break
    if lineage_col is None:
        raise ValueError('Cannot find cell type / lineage column. Need cell_type_major or cell_type.')
    obs['lineage_label_raw'] = obs[lineage_col].astype(str)
    obs['lineage_label_plot'] = obs['lineage_label_raw'].replace({'T_NK': 'T/NK', 'B_Plasma': 'B/Plasma', 'Mast_cell': 'Mast cell', 'Enteroendocrine': 'Enteroendocrine'})
    if 'polyp_transition_state' in obs.columns:
        obs['polyp_transition_state'] = obs['polyp_transition_state'].astype(str)
    elif 'polyp_state_label' in obs.columns:
        obs['polyp_transition_state'] = obs['polyp_state_label'].astype(str)
    else:
        obs['is_unknown'] = obs['lineage_label_raw'].map(is_unknown_label)
        obs['polyp_transition_state'] = 'Other_cells'
        obs.loc[obs['is_unknown'] & (obs['stage'] != 'Polyp'), 'polyp_transition_state'] = 'Unknown_nonPolyp'
        obs.loc[obs['is_unknown'] & (obs['stage'] == 'Polyp'), 'polyp_transition_state'] = 'Polyp_enriched_unclassified_state'
    obs['is_polyp_transition'] = obs['polyp_transition_state'].astype(str).eq('Polyp_enriched_unclassified_state')
    obs['lineage_plus_polyp_state'] = obs['lineage_label_plot'].astype(str)
    obs.loc[obs['is_polyp_transition'], 'lineage_plus_polyp_state'] = 'Polyp-enriched transition'
    obs = obs[obs['stage'].isin(STAGES)].copy()
    adata = adata[obs.index].copy()
    adata.obs = obs
    n_target = int(obs['is_polyp_transition'].sum())
    print('Cells:', adata.n_obs)
    print('Stage counts:\n', obs['stage'].value_counts())
    print('Lineage column:', lineage_col)
    print('Polyp-enriched transition cells:', n_target)
    if n_target < MIN_POLYP_UNKNOWN_CELLS:
        print(f'Warning: only {n_target} Polyp-enriched transition cells found.')
    obs.to_csv(OUT_DIR / 'cell_metadata_used_for_result1.csv', encoding='utf-8-sig')
    return adata

def ensure_umap(adata):
    if 'X_umap' in adata.obsm:
        return True
    if not COMPUTE_UMAP_IF_MISSING:
        return False
    try:
        import scanpy as sc
        print('No UMAP found. Computing UMAP safely.')
        tmp = adata.copy()
        X = safe_X(tmp)
        try:
            from scipy import sparse
            if sparse.issparse(X):
                xmin = float(X.data.min()) if X.data.size else 0.0
                xmax = float(X.data.max()) if X.data.size else 0.0
            else:
                xmin = float(np.nanmin(X))
                xmax = float(np.nanmax(X))
        except Exception:
            xmin, xmax = (0.0, 0.0)
        if xmin >= 0 and xmax > 30:
            sc.pp.normalize_total(tmp, target_sum=10000.0)
            sc.pp.log1p(tmp)
        n_pcs = min(30, tmp.n_obs - 1, tmp.n_vars - 1)
        n_pcs = max(2, n_pcs)
        sc.pp.pca(tmp, n_comps=n_pcs)
        sc.pp.neighbors(tmp, n_neighbors=min(15, tmp.n_obs - 1), n_pcs=n_pcs)
        sc.tl.umap(tmp, random_state=RANDOM_SEED)
        adata.obsm['X_umap'] = tmp.obsm['X_umap']
        return True
    except Exception as e:
        print('Failed to compute UMAP:', e)
        return False

def group_definitions(obs):
    return {'Normal': np.where(obs['stage'].values == 'Normal')[0], 'Polyp other': np.where((obs['stage'].values == 'Polyp') & ~obs['is_polyp_transition'].values)[0], 'Polyp transition': np.where(obs['is_polyp_transition'].values)[0], 'CRC': np.where(obs['stage'].values == 'CRC')[0]}

def marker_presence(adata):
    genes = set([norm_gene(g) for g in adata.var_names])
    rows = []
    for sig, markers in SIGNATURES.items():
        clean = list(dict.fromkeys([norm_gene(g) for g in markers]))
        found = [g for g in clean if g in genes]
        rows.append({'signature': sig, 'n_total': len(clean), 'n_found': len(found), 'found_fraction': len(found) / max(1, len(clean)), 'markers_found': ','.join(found), 'markers_missing': ','.join([g for g in clean if g not in genes])})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / 'signature_marker_presence.csv', index=False, encoding='utf-8-sig')
    return df

def compute_signature_summary(adata):
    from scipy import sparse
    marker_presence(adata)
    X = safe_X(adata)
    genes = pd.Index([norm_gene(g) for g in adata.var_names])
    gene_to_idx = {}
    for i, g in enumerate(genes):
        if g not in gene_to_idx:
            gene_to_idx[g] = i
    if sparse.issparse(X):
        gene_mean = np.asarray(X.mean(axis=0)).ravel()
        gene_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        gene_std = np.sqrt(np.maximum(gene_sq - gene_mean ** 2, 1e-08))
    else:
        gene_mean = np.nanmean(X, axis=0)
        gene_std = np.nanstd(X, axis=0)
        gene_std[gene_std == 0] = 1.0
    groups = group_definitions(adata.obs)
    score_per_cell = pd.DataFrame(index=adata.obs.index)
    rows = []
    for sig, markers in SIGNATURES.items():
        present = [norm_gene(g) for g in markers if norm_gene(g) in gene_to_idx]
        if len(present) == 0:
            continue
        idx = [gene_to_idx[g] for g in present]
        sub = sparse_to_dense(X[:, idx]).astype(float)
        sub = np.nan_to_num(sub, nan=0.0, posinf=0.0, neginf=0.0)
        z = (sub - gene_mean[idx]) / gene_std[idx]
        score = np.nanmean(z, axis=1)
        score_per_cell[sig] = score
        for group, cell_idx in groups.items():
            if len(cell_idx) == 0:
                continue
            vals = score[cell_idx]
            rows.append({'signature': sig, 'group': group, 'n_cells': len(cell_idx), 'n_markers_present': len(present), 'markers_present': ','.join(present), 'mean_score': float(np.nanmean(vals)), 'median_score': float(np.nanmedian(vals))})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / 'signature_score_summary_by_group.csv', index=False, encoding='utf-8-sig')
    score_per_cell.to_csv(OUT_DIR / 'signature_scores_per_cell.csv.gz', compression='gzip')
    target = summary[summary['group'] == 'Polyp transition'].sort_values('mean_score', ascending=False)
    target.to_csv(OUT_DIR / 'signature_rank_in_polyp_transition_state.csv', index=False, encoding='utf-8-sig')
    return (summary, score_per_cell)

def compute_key_gene_group_matrix(adata):
    X = safe_X(adata)
    genes = pd.Index([norm_gene(g) for g in adata.var_names])
    gene_to_idx = {}
    for i, g in enumerate(genes):
        if g not in gene_to_idx:
            gene_to_idx[g] = i
    present = [g for g in KEY_GENES_FOR_HEATMAP if g in gene_to_idx]
    groups = group_definitions(adata.obs)
    group_order = ['Normal', 'Polyp other', 'Polyp transition', 'CRC']
    rows = []
    for g in present:
        idx = gene_to_idx[g]
        vals = []
        for group in group_order:
            cell_idx = groups.get(group, [])
            if len(cell_idx) == 0:
                vals.append(np.nan)
            else:
                arr = sparse_to_dense(X[cell_idx, idx])
                vals.append(float(np.nanmean(arr)))
        rows.append(vals)
    mat = pd.DataFrame(rows, index=present, columns=group_order)
    z = zscore_rows(mat).clip(-2, 2)
    mat.to_csv(OUT_DIR / 'literature_key_gene_group_mean_expression.csv', encoding='utf-8-sig')
    z.to_csv(OUT_DIR / 'literature_key_gene_group_mean_expression_zscore.csv', encoding='utf-8-sig')
    return (mat, z)

def panel_umap_stage(ax, adata):
    coords = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index, columns=['UMAP1', 'UMAP2'])
    obs = adata.obs.copy()
    for st in STAGES:
        idx = obs['stage'] == st
        if idx.sum() == 0:
            continue
        ax.scatter(coords.loc[idx, 'UMAP1'], coords.loc[idx, 'UMAP2'], s=POINT_SIZE, color=STAGE_COLORS[st], alpha=0.75, linewidths=0, label=st)
    ax.set_title('Stage distribution', fontsize=12, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.legend(frameon=False, fontsize=8, markerscale=3, loc='best')

def panel_umap_lineage(ax, adata):
    coords = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index, columns=['UMAP1', 'UMAP2'])
    obs = adata.obs.copy()
    labels = obs['lineage_plus_polyp_state'].astype(str)
    ordered = [x for x in sorted(labels.unique()) if x != 'Polyp-enriched transition']
    if 'Polyp-enriched transition' in set(labels):
        ordered.append('Polyp-enriched transition')
    colors = get_palette(ordered, LINEAGE_COLORS)
    for lab in ordered:
        idx = labels == lab
        if idx.sum() == 0:
            continue
        size = POINT_SIZE * 1.9 if lab == 'Polyp-enriched transition' else POINT_SIZE
        alpha = 0.9 if lab == 'Polyp-enriched transition' else 0.65
        ax.scatter(coords.loc[idx, 'UMAP1'], coords.loc[idx, 'UMAP2'], s=size, color=colors.get(lab, '#999999'), alpha=alpha, linewidths=0, label=lab)
    ax.set_title('Major lineage and Polyp-enriched state', fontsize=12, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.legend(frameon=False, fontsize=7, markerscale=3, bbox_to_anchor=(1.02, 1), loc='upper left')

def panel_polyp_fraction(ax, adata):
    obs = adata.obs.copy()
    count = pd.crosstab(obs['stage'], obs['polyp_transition_state']).reindex(STAGES).fillna(0)
    frac = count.div(count.sum(axis=1), axis=0).fillna(0)
    frac.to_csv(OUT_DIR / 'polyp_transition_fraction_by_stage.csv', encoding='utf-8-sig')
    count.to_csv(OUT_DIR / 'polyp_transition_count_by_stage.csv', encoding='utf-8-sig')
    categories = ['Other_cells', 'Unknown_nonPolyp', 'Polyp_enriched_unclassified_state']
    labels = {'Other_cells': 'Other cells', 'Unknown_nonPolyp': 'Unknown non-Polyp', 'Polyp_enriched_unclassified_state': 'Polyp-enriched\ntransition'}
    colors = {'Other_cells': '#D0D0D0', 'Unknown_nonPolyp': '#8F8F8F', 'Polyp_enriched_unclassified_state': '#D1495B'}
    x = np.arange(len(STAGES))
    bottom = np.zeros(len(STAGES))
    for c in categories:
        if c not in frac.columns:
            continue
        y = frac[c].values
        ax.bar(x, y, bottom=bottom, color=colors[c], edgecolor='white', linewidth=0.6, label=labels[c])
        bottom += y
    ax.set_xticks(x)
    ax.set_xticklabels(STAGES)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Fraction of cells')
    ax.set_title('Stage enrichment of transition state', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.2)
    ax.legend(frameon=False, fontsize=7, loc='upper left', bbox_to_anchor=(1.02, 1))

def panel_signature_dotplot(ax, summary):
    """
    Signature score dotplot with legends.

    Dot size and dot color both represent the row-scaled mean signature score.
    The number printed in the Polyp transition column is the original mean signature score.
    """
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    group_order = ['Normal', 'Polyp other', 'Polyp transition', 'CRC']
    sig_order = ['General epithelial', 'Stem-like epithelial', 'Proliferative TA', 'Malignancy continuum', 'CRC-like epithelial', 'Mature enterocyte', 'Goblet / immature goblet', 'preCAF-like', 'Treg-like', 'Stromal', 'Myeloid', 'T/NK']
    group_order = [g for g in group_order if g in set(summary['group'])]
    sig_order = [s for s in sig_order if s in set(summary['signature'])]
    mat = summary.pivot(index='signature', columns='group', values='mean_score').reindex(sig_order)[group_order].fillna(0)
    scaled = mat.copy()
    for sig in scaled.index:
        scaled.loc[sig] = minmax01(scaled.loc[sig]).values
    cmap = plt.cm.RdYlBu_r
    norm = Normalize(vmin=0, vmax=1)
    for i, sig in enumerate(sig_order):
        for j, group in enumerate(group_order):
            val = float(scaled.loc[sig, group])
            raw = float(mat.loc[sig, group])
            size = 45 + 440 * max(val, 0)
            ax.scatter(j, i, s=size, color=cmap(norm(val)), edgecolor='white', linewidth=0.4, alpha=0.88)
            if group == 'Polyp transition':
                ax.text(j, i, f'{raw:.2f}', ha='center', va='center', fontsize=5.5)
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels([g.replace(' ', '\n') for g in group_order], fontsize=7)
    ax.set_yticks(range(len(sig_order)))
    ax.set_yticklabels(sig_order, fontsize=7)
    ax.invert_yaxis()
    ax.set_title('Transition-state signature scores', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.12)
    size_values = [0.25, 0.5, 0.75, 1.0]
    size_handles = [ax.scatter([], [], s=45 + 440 * v, color=cmap(norm(v)), edgecolor='white', linewidth=0.4, alpha=0.88) for v in size_values]
    leg = ax.legend(size_handles, [f'{v:.2f}' for v in size_values], title='Row-scaled\nscore', frameon=False, fontsize=7, title_fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0, handletextpad=1.0, labelspacing=0.9)
    ax.add_artist(leg)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = ax.figure.colorbar(sm, ax=ax, fraction=0.04, pad=0.13)
    cbar.set_label('Row-scaled score', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

def panel_key_gene_heatmap(ax, zmat):
    if zmat.empty:
        ax.text(0.5, 0.5, 'No key genes found', ha='center', va='center')
        ax.axis('off')
        return
    im = ax.imshow(zmat.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    ax.set_xticks(range(zmat.shape[1]))
    ax.set_xticklabels([c.replace(' ', '\n') for c in zmat.columns], fontsize=7)
    ax.set_yticks(range(zmat.shape[0]))
    ax.set_yticklabels(zmat.index, fontsize=6)
    ax.set_title('Literature-derived key genes', fontsize=12, fontweight='bold')
    return im

def plot_individual_panels(adata, summary, zmat):
    if ensure_umap(adata):
        fig, ax = plt.subplots(figsize=(6.2, 5.5))
        panel_umap_stage(ax, adata)
        savefig(fig, 'Fig1A_UMAP_by_stage')
        fig, ax = plt.subplots(figsize=(7.6, 5.5))
        panel_umap_lineage(ax, adata)
        savefig(fig, 'Fig1B_UMAP_by_major_lineage_polyp_state')
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    panel_polyp_fraction(ax, adata)
    savefig(fig, 'Fig1C_polyp_transition_fraction_by_stage')
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    panel_signature_dotplot(ax, summary)
    savefig(fig, 'Fig1D_transition_signature_score_dotplot')
    fig, ax = plt.subplots(figsize=(6.2, max(6.0, 0.22 * max(1, zmat.shape[0]))))
    im = panel_key_gene_heatmap(ax, zmat)
    if im is not None:
        cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label('Row z-score', fontsize=8)
    savefig(fig, 'Fig1E_literature_key_gene_heatmap')

def plot_composite_figure1(adata, summary, zmat):
    if not ensure_umap(adata):
        print('Cannot make composite Figure 1 because UMAP is unavailable.')
        return
    fig = plt.figure(figsize=(17.5, 14.5))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.12], height_ratios=[1.0, 0.82, 1.2], wspace=0.38, hspace=0.42)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, :])
    panel_umap_stage(axA, adata)
    panel_umap_lineage(axB, adata)
    panel_polyp_fraction(axC, adata)
    panel_signature_dotplot(axD, summary)
    im = panel_key_gene_heatmap(axE, zmat)
    if im is not None:
        cbar = fig.colorbar(im, ax=axE, fraction=0.02, pad=0.01)
        cbar.set_label('Row z-score', fontsize=8)
    for ax, letter in zip([axA, axB, axC, axD, axE], list('ABCDE')):
        ax.text(-0.08, 1.08, letter, transform=ax.transAxes, fontsize=22, fontweight='bold', va='top', ha='right')
    fig.suptitle('Overall data distribution and Polyp-enriched transition state', fontsize=18, fontweight='bold', y=0.995)
    savefig(fig, 'Figure1_overall_data_distribution_polyp_transition_state')

def write_result_text_template(summary):
    target = summary[summary['group'] == 'Polyp transition'].sort_values('mean_score', ascending=False)
    top_sigs = target.head(5)['signature'].tolist()
    top_sig_text = '、'.join(top_sigs)
    text = f'\n3.1 Normal–Polyp–CRC 数据整体分布与 Polyp 阶段过渡状态识别\n\n为从整体层面评估 Normal、Polyp 和 CRC 三个阶段样本的细胞状态差异，本研究首先基于单细胞 RNA 数据进行低维可视化，并结合保守谱系 marker 对主要细胞群进行大类注释。UMAP 可视化结果显示，不同病理阶段的细胞在低维空间中呈现一定程度的分离，提示 Normal–Polyp–CRC 进展过程中伴随转录状态变化。主要细胞谱系注释显示，数据中包含上皮、免疫、基质和内皮等主要细胞成分，为后续阶段特异性调控网络分析提供了基础。\n\n值得注意的是，未能被稳定注释的细胞主要集中于 Polyp 阶段，而非随机分布于 Normal、Polyp 和 CRC 三个阶段。考虑到 Polyp 是 Normal 向 CRC 转化过程中的中间病理状态，且既往研究提示息肉阶段可能存在 stem-like epithelial、preCAF-like 和免疫抑制相关过渡状态，本研究未将该群体强行归入某一传统细胞类型，而是将其定义为 Polyp-enriched unclassified state，也即 Polyp 阶段富集的未分类过渡状态。\n\n进一步地，本研究基于文献报道的状态相关基因集计算 signature score，包括 stem-like epithelial、proliferative TA、mature enterocyte、goblet/immature goblet、preCAF-like、Treg-like 和 CRC-like epithelial 等状态特征。Polyp-enriched transition state 中得分较高的 signature 包括：{top_sig_text}。该结果提示，Polyp 阶段未分类细胞并非单纯的注释噪声，而可能反映息肉阶段特异的转录状态重塑。\n\n上述结果说明，Polyp 阶段可能已经出现介于正常组织与 CRC 之间的过渡性状态变化。因此，后续分析进一步从 TF–TF 调控网络层面解析 Normal–Polyp–CRC 进展过程中潜在的阶段特异性调控重塑，并识别可能参与该过程的核心转录因子。\n'
    (OUT_DIR / 'RESULT1_WRITING_TEMPLATE.txt').write_text(text, encoding='utf-8')

def main():
    print('=' * 100)
    print('Result 1: Overall data distribution and Polyp-enriched transition state')
    print('=' * 100)
    print('OUT_DIR:', OUT_DIR)
    adata = load_data()
    summary, scores = compute_signature_summary(adata)
    mat, zmat = compute_key_gene_group_matrix(adata)
    plot_individual_panels(adata, summary, zmat)
    plot_composite_figure1(adata, summary, zmat)
    write_result_text_template(summary)
    adata.obs.to_csv(OUT_DIR / 'cell_metadata_result1_final.csv', encoding='utf-8-sig')
    print('\nDONE')
    print('Main figure:')
    print(OUT_DIR / 'Figure1_overall_data_distribution_polyp_transition_state.png')
    print('Writing template:')
    print(OUT_DIR / 'RESULT1_WRITING_TEMPLATE.txt')
'''

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    a=load_source('landscape',SOURCE_landscape)
    adata=a.load_data()
    summary,_=a.compute_signature_summary(adata)
    _,zmat=a.compute_key_gene_group_matrix(adata)
    import matplotlib.pyplot as plt
    for label, fn in [('A',lambda ax:a.panel_umap_stage(ax,adata)),('B',lambda ax:a.panel_umap_lineage(ax,adata)),('C',lambda ax:a.panel_signature_dotplot(ax,summary)),('D',lambda ax:a.panel_key_gene_heatmap(ax,zmat))]:
        def draw(label=label, fn=fn):
            fig,ax=plt.subplots(figsize=(12,8));fn(ax)
            path=OUT/(label+'.png');fig.savefig(path,dpi=300,bbox_inches='tight');fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig);panel_file(label,path)
        attempt(label,draw)

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        missing("运行", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
    finally:
        finish(4)
