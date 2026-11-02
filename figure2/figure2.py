# -*- coding: utf-8 -*-
##### 图2 #####
from pathlib import Path
import os, traceback, json
DATA_DIR = Path(os.environ.get('CRC_DATA_DIR', r'E:\毕业设计\bishedata'))
OUTPUT_ROOT = Path(os.environ.get('CRC_OUTPUT_DIR', str(Path(__file__).resolve().parent / '绘图输出')))
FIGURE_NAME = '图2'
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
from pathlib import Path as landscape_Path
import warnings as landscape_warnings
import numpy as landscape_np
import pandas as landscape_pd
import matplotlib.pyplot as landscape_plt
landscape_warnings.filterwarnings('ignore')
landscape_BASE_DIR = landscape_Path(str(DATA_DIR))
landscape_ATACMOTIF_DIR = landscape_BASE_DIR / '4.0result' / 'ATACmotif'
landscape_OUT_DIR = landscape_ATACMOTIF_DIR / '50_result1_overall_polyp_transition'
landscape_H5AD_CANDIDATES = [landscape_ATACMOTIF_DIR / '49_polyp_transition_signature_becker2022' / 'RNA_with_polyp_transition_state_becker2022.h5ad', landscape_ATACMOTIF_DIR / '48_polyp_enriched_unclassified_state' / 'RNA_with_polyp_enriched_unclassified_state.h5ad', landscape_ATACMOTIF_DIR / '46_layered_cell_annotation' / '01_major_lineage_all_data' / 'RNA_major_lineage_annotated.h5ad', landscape_ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'RNA_annotated_celltypes_final_markers.h5ad', landscape_ATACMOTIF_DIR / '33_celltype_annotation' / 'RNA_annotated_celltypes.h5ad', landscape_BASE_DIR / '4.0版本数据' / 'GRNBoost2_input_FINAL' / 'GRNBoost2_ready.h5ad']
landscape_META_CANDIDATES = [landscape_ATACMOTIF_DIR / '49_polyp_transition_signature_becker2022' / 'cell_metadata_with_polyp_transition_state.csv', landscape_ATACMOTIF_DIR / '48_polyp_enriched_unclassified_state' / 'cell_metadata_with_polyp_enriched_state.csv', landscape_ATACMOTIF_DIR / '46_layered_cell_annotation' / '01_major_lineage_all_data' / 'cell_metadata_major_lineage.csv', landscape_ATACMOTIF_DIR / '44_celltype_final_marker_annotation' / 'cell_metadata_with_final_marker_cell_type.csv', landscape_ATACMOTIF_DIR / '43_celltype_unknown_diagnosis_reannotation' / 'cell_metadata_with_cell_type_reannotated.csv', landscape_ATACMOTIF_DIR / '33_celltype_annotation' / 'cell_metadata_with_cell_type.csv']
landscape_STAGES = ['Normal', 'Polyp', 'CRC']
landscape_UNKNOWN_LABELS = {'unknown', 'unclassified', 'unassigned', 'other', 'nan', 'none', 'low_confidence', 'low_confidence_unknown'}
landscape_MIN_POLYP_UNKNOWN_CELLS = 20
landscape_POINT_SIZE = 4
landscape_RANDOM_SEED = 42
landscape_COMPUTE_UMAP_IF_MISSING = True
landscape_SIGNATURES = {'Stem-like epithelial': ['LGR5', 'SMOC2', 'RGMB', 'PTPRO', 'EPHB2', 'LRIG1', 'ASCL2', 'OLFM4', 'SOX9'], 'Malignancy continuum': ['OLFM4', 'GPX2', 'ASCL2', 'HNF4A', 'SOX4', 'KRT18', 'NORAD'], 'Mature enterocyte': ['CA1', 'CA2', 'AQP8', 'SLC26A3', 'GUCA2A', 'GUCA2B', 'RAB6B', 'BEST4', 'CA7', 'OTOP2', 'OTOP3', 'SLC4A4', 'NR5A2'], 'Goblet / immature goblet': ['MUC2', 'TFF1', 'TFF3', 'FCGBP', 'TBX10', 'KLK1', 'ITLN1', 'WFDC2', 'CLCA1', 'AGR2'], 'Proliferative TA': ['MKI67', 'TOP2A', 'PCNA', 'STMN1', 'UBE2C', 'CENPF', 'HMGB2', 'MCM2', 'MCM5'], 'preCAF-like': ['FAP', 'VCAN', 'COL1A2', 'COL11A1', 'THBS2', 'TWIST1', 'WNT2', 'RUNX1', 'SERPINE1', 'INHBA', 'PLXDC2', 'BICC1', 'FBN1'], 'Treg-like': ['FOXP3', 'CTLA4', 'IL2RA', 'TNFRSF4', 'LAIR2', 'IKZF2'], 'CRC-like epithelial': ['GPX2', 'HNF4A', 'CEACAM5', 'CEACAM6', 'CLDN1', 'MMP7', 'KRT19', 'TACSTD2'], 'General epithelial': ['EPCAM', 'TACSTD2', 'KRT8', 'KRT18', 'KRT19', 'KRT20', 'CDH1'], 'Stromal': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'ACTA2', 'TAGLN', 'VIM', 'PDGFRA'], 'Myeloid': ['LYZ', 'LST1', 'FCGR3A', 'S100A8', 'S100A9', 'C1QA', 'C1QB', 'CD68'], 'T/NK': ['PTPRC', 'CD3D', 'CD3E', 'TRAC', 'IL7R', 'NKG7', 'GNLY', 'PRF1']}
landscape_KEY_GENES_FOR_HEATMAP = ['LGR5', 'SMOC2', 'ASCL2', 'OLFM4', 'SOX9', 'GPX2', 'HNF4A', 'CA1', 'CA2', 'AQP8', 'SLC26A3', 'BEST4', 'MUC2', 'TFF3', 'FCGBP', 'CLCA1', 'MKI67', 'TOP2A', 'PCNA', 'FAP', 'VCAN', 'COL1A2', 'COL11A1', 'WNT2', 'RUNX1', 'FOXP3', 'CTLA4', 'IL2RA', 'EPCAM', 'KRT8', 'KRT18', 'KRT19']
landscape_STAGE_COLORS = {'Normal': '#4C78A8', 'Polyp': '#F58518', 'CRC': '#E45756'}
landscape_LINEAGE_COLORS = {'Epithelial': '#4C78A8', 'T_NK': '#54A24B', 'T/NK': '#54A24B', 'B_Plasma': '#72B7B2', 'B/Plasma': '#72B7B2', 'Myeloid': '#F58518', 'Stromal': '#9D755D', 'Endothelial': '#B279A2', 'Mast_cell': '#ECA82C', 'Enteroendocrine': '#8CD17D', 'Unknown': '#B8B8B8', 'Unclassified': '#B8B8B8', 'Polyp-enriched transition': '#D1495B'}
landscape_DEFAULT_COLORS = ['#4C78A8', '#F58518', '#54A24B', '#E45756', '#72B7B2', '#B279A2', '#9D755D', '#BAB0AC', '#6A4C93', '#2A9D8F', '#ECA82C', '#8CD17D']

def landscape_first_existing(paths):
    for p in paths:
        p = landscape_Path(p)
        if p.exists():
            return p
    return None

def landscape_norm_gene(x):
    return str(x).strip().upper()

def landscape_norm_stage(x):
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
    return landscape_np.nan

def landscape_is_unknown_label(x):
    s = str(x).strip().lower()
    return s in landscape_UNKNOWN_LABELS or 'unknown' in s or 'unclassified' in s or ('unassigned' in s)

def landscape_sparse_to_dense(X):
    try:
        from scipy import sparse as sparse
        if sparse.issparse(X):
            return X.toarray()
    except Exception:
        pass
    return landscape_np.asarray(X)

def landscape_safe_X(adata):
    try:
        from scipy import sparse as sparse
        if sparse.issparse(adata.X):
            X = adata.X.copy()
            X.data = landscape_np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
            return X
    except Exception:
        pass
    return landscape_np.nan_to_num(landscape_np.asarray(adata.X, dtype=landscape_np.float32), nan=0.0, posinf=0.0, neginf=0.0)

def landscape_minmax01(x):
    x = landscape_pd.Series(x, dtype=float).replace([landscape_np.inf, -landscape_np.inf], landscape_np.nan)
    if x.dropna().empty:
        return landscape_pd.Series(landscape_np.zeros(len(x)), index=x.index)
    x = x.fillna(x.median())
    if float(x.max()) == float(x.min()):
        return landscape_pd.Series(landscape_np.ones(len(x)) * 0.5, index=x.index)
    return (x - x.min()) / (x.max() - x.min())

def landscape_zscore_rows(df):
    out = df.astype(float).copy()
    for idx in out.index:
        x = out.loc[idx].astype(float)
        sd = x.std()
        out.loc[idx] = 0 if sd == 0 or landscape_pd.isna(sd) else (x - x.mean()) / sd
    return out.replace([landscape_np.inf, -landscape_np.inf], 0).fillna(0)

def landscape_get_palette(labels, base_colors=None):
    labels = list(landscape_pd.Series(labels).dropna().astype(str).unique())
    colors = {}
    j = 0
    for lab in labels:
        if base_colors and lab in base_colors:
            colors[lab] = base_colors[lab]
        elif lab in landscape_LINEAGE_COLORS:
            colors[lab] = landscape_LINEAGE_COLORS[lab]
        else:
            colors[lab] = landscape_DEFAULT_COLORS[j % len(landscape_DEFAULT_COLORS)]
            j += 1
    return colors

def landscape_load_data():
    import anndata as ad
    h5ad_path = landscape_first_existing(landscape_H5AD_CANDIDATES)
    if h5ad_path is None:
        raise FileNotFoundError('Cannot find h5ad. Please update H5AD_CANDIDATES.')
    print('Using h5ad:', h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    adata.obs.index = adata.obs.index.astype(str)
    adata.var_names = adata.var_names.astype(str)
    meta_path = landscape_first_existing(landscape_META_CANDIDATES)
    if meta_path is not None:
        print('Using metadata:', meta_path)
        meta = landscape_pd.read_csv(meta_path, index_col=0)
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
        obs['stage'] = obs['stage'].map(landscape_norm_stage)
    elif 'stage_label' in obs.columns:
        obs['stage'] = obs['stage_label'].map(landscape_norm_stage)
    elif 'sample' in obs.columns:
        obs['stage'] = obs['sample'].map(landscape_norm_stage)
    elif 'sample_id' in obs.columns:
        obs['stage'] = obs['sample_id'].map(landscape_norm_stage)
    else:
        obs['stage'] = [landscape_norm_stage(x) for x in obs.index]
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
        obs['is_unknown'] = obs['lineage_label_raw'].map(landscape_is_unknown_label)
        obs['polyp_transition_state'] = 'Other_cells'
        obs.loc[obs['is_unknown'] & (obs['stage'] != 'Polyp'), 'polyp_transition_state'] = 'Unknown_nonPolyp'
        obs.loc[obs['is_unknown'] & (obs['stage'] == 'Polyp'), 'polyp_transition_state'] = 'Polyp_enriched_unclassified_state'
    obs['is_polyp_transition'] = obs['polyp_transition_state'].astype(str).eq('Polyp_enriched_unclassified_state')
    obs['lineage_plus_polyp_state'] = obs['lineage_label_plot'].astype(str)
    obs.loc[obs['is_polyp_transition'], 'lineage_plus_polyp_state'] = 'Polyp-enriched transition'
    obs = obs[obs['stage'].isin(landscape_STAGES)].copy()
    adata = adata[obs.index].copy()
    adata.obs = obs
    n_target = int(obs['is_polyp_transition'].sum())
    print('Cells:', adata.n_obs)
    print('Stage counts:\n', obs['stage'].value_counts())
    print('Lineage column:', lineage_col)
    print('Polyp-enriched transition cells:', n_target)
    if n_target < landscape_MIN_POLYP_UNKNOWN_CELLS:
        print(f'Warning: only {n_target} Polyp-enriched transition cells found.')
    obs.to_csv(landscape_OUT_DIR / 'cell_metadata_used_for_result1.csv', encoding='utf-8-sig')
    return adata

def landscape_group_definitions(obs):
    return {'Normal': landscape_np.where(obs['stage'].values == 'Normal')[0], 'Polyp other': landscape_np.where((obs['stage'].values == 'Polyp') & ~obs['is_polyp_transition'].values)[0], 'Polyp transition': landscape_np.where(obs['is_polyp_transition'].values)[0], 'CRC': landscape_np.where(obs['stage'].values == 'CRC')[0]}

def landscape_marker_presence(adata):
    genes = set([landscape_norm_gene(g) for g in adata.var_names])
    rows = []
    for sig, markers in landscape_SIGNATURES.items():
        clean = list(dict.fromkeys([landscape_norm_gene(g) for g in markers]))
        found = [g for g in clean if g in genes]
        rows.append({'signature': sig, 'n_total': len(clean), 'n_found': len(found), 'found_fraction': len(found) / max(1, len(clean)), 'markers_found': ','.join(found), 'markers_missing': ','.join([g for g in clean if g not in genes])})
    df = landscape_pd.DataFrame(rows)
    df.to_csv(landscape_OUT_DIR / 'signature_marker_presence.csv', index=False, encoding='utf-8-sig')
    return df

def landscape_compute_signature_summary(adata):
    from scipy import sparse as sparse
    landscape_marker_presence(adata)
    X = landscape_safe_X(adata)
    genes = landscape_pd.Index([landscape_norm_gene(g) for g in adata.var_names])
    gene_to_idx = {}
    for i, g in enumerate(genes):
        if g not in gene_to_idx:
            gene_to_idx[g] = i
    if sparse.issparse(X):
        gene_mean = landscape_np.asarray(X.mean(axis=0)).ravel()
        gene_sq = landscape_np.asarray(X.multiply(X).mean(axis=0)).ravel()
        gene_std = landscape_np.sqrt(landscape_np.maximum(gene_sq - gene_mean ** 2, 1e-08))
    else:
        gene_mean = landscape_np.nanmean(X, axis=0)
        gene_std = landscape_np.nanstd(X, axis=0)
        gene_std[gene_std == 0] = 1.0
    groups = landscape_group_definitions(adata.obs)
    score_per_cell = landscape_pd.DataFrame(index=adata.obs.index)
    rows = []
    for sig, markers in landscape_SIGNATURES.items():
        present = [landscape_norm_gene(g) for g in markers if landscape_norm_gene(g) in gene_to_idx]
        if len(present) == 0:
            continue
        idx = [gene_to_idx[g] for g in present]
        sub = landscape_sparse_to_dense(X[:, idx]).astype(float)
        sub = landscape_np.nan_to_num(sub, nan=0.0, posinf=0.0, neginf=0.0)
        z = (sub - gene_mean[idx]) / gene_std[idx]
        score = landscape_np.nanmean(z, axis=1)
        score_per_cell[sig] = score
        for group, cell_idx in groups.items():
            if len(cell_idx) == 0:
                continue
            vals = score[cell_idx]
            rows.append({'signature': sig, 'group': group, 'n_cells': len(cell_idx), 'n_markers_present': len(present), 'markers_present': ','.join(present), 'mean_score': float(landscape_np.nanmean(vals)), 'median_score': float(landscape_np.nanmedian(vals))})
    summary = landscape_pd.DataFrame(rows)
    summary.to_csv(landscape_OUT_DIR / 'signature_score_summary_by_group.csv', index=False, encoding='utf-8-sig')
    score_per_cell.to_csv(landscape_OUT_DIR / 'signature_scores_per_cell.csv.gz', compression='gzip')
    target = summary[summary['group'] == 'Polyp transition'].sort_values('mean_score', ascending=False)
    target.to_csv(landscape_OUT_DIR / 'signature_rank_in_polyp_transition_state.csv', index=False, encoding='utf-8-sig')
    return (summary, score_per_cell)

def landscape_compute_key_gene_group_matrix(adata):
    X = landscape_safe_X(adata)
    genes = landscape_pd.Index([landscape_norm_gene(g) for g in adata.var_names])
    gene_to_idx = {}
    for i, g in enumerate(genes):
        if g not in gene_to_idx:
            gene_to_idx[g] = i
    present = [g for g in landscape_KEY_GENES_FOR_HEATMAP if g in gene_to_idx]
    groups = landscape_group_definitions(adata.obs)
    group_order = ['Normal', 'Polyp other', 'Polyp transition', 'CRC']
    rows = []
    for g in present:
        idx = gene_to_idx[g]
        vals = []
        for group in group_order:
            cell_idx = groups.get(group, [])
            if len(cell_idx) == 0:
                vals.append(landscape_np.nan)
            else:
                arr = landscape_sparse_to_dense(X[cell_idx, idx])
                vals.append(float(landscape_np.nanmean(arr)))
        rows.append(vals)
    mat = landscape_pd.DataFrame(rows, index=present, columns=group_order)
    z = landscape_zscore_rows(mat).clip(-2, 2)
    mat.to_csv(landscape_OUT_DIR / 'literature_key_gene_group_mean_expression.csv', encoding='utf-8-sig')
    z.to_csv(landscape_OUT_DIR / 'literature_key_gene_group_mean_expression_zscore.csv', encoding='utf-8-sig')
    return (mat, z)

##### 图2A #####
def landscape_panel_umap_stage(ax, adata):
    coords = landscape_pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index, columns=['UMAP1', 'UMAP2'])
    obs = adata.obs.copy()
    for st in landscape_STAGES:
        idx = obs['stage'] == st
        if idx.sum() == 0:
            continue
        ax.scatter(coords.loc[idx, 'UMAP1'], coords.loc[idx, 'UMAP2'], s=landscape_POINT_SIZE, color=landscape_STAGE_COLORS[st], alpha=0.75, linewidths=0, label=st)
    ax.set_title('Stage distribution', fontsize=12, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.legend(frameon=False, fontsize=8, markerscale=3, loc='best')

##### 图2B #####
def landscape_panel_umap_lineage(ax, adata):
    coords = landscape_pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index, columns=['UMAP1', 'UMAP2'])
    obs = adata.obs.copy()
    labels = obs['lineage_plus_polyp_state'].astype(str)
    ordered = [x for x in sorted(labels.unique()) if x != 'Polyp-enriched transition']
    if 'Polyp-enriched transition' in set(labels):
        ordered.append('Polyp-enriched transition')
    colors = landscape_get_palette(ordered, landscape_LINEAGE_COLORS)
    for lab in ordered:
        idx = labels == lab
        if idx.sum() == 0:
            continue
        size = landscape_POINT_SIZE * 1.9 if lab == 'Polyp-enriched transition' else landscape_POINT_SIZE
        alpha = 0.9 if lab == 'Polyp-enriched transition' else 0.65
        ax.scatter(coords.loc[idx, 'UMAP1'], coords.loc[idx, 'UMAP2'], s=size, color=colors.get(lab, '#999999'), alpha=alpha, linewidths=0, label=lab)
    ax.set_title('Major lineage and Polyp-enriched state', fontsize=12, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.legend(frameon=False, fontsize=7, markerscale=3, bbox_to_anchor=(1.02, 1), loc='upper left')

##### 图2C #####
def landscape_panel_signature_dotplot(ax, summary):
    from matplotlib.cm import ScalarMappable as ScalarMappable
    from matplotlib.colors import Normalize as Normalize
    group_order = ['Normal', 'Polyp other', 'Polyp transition', 'CRC']
    sig_order = ['General epithelial', 'Stem-like epithelial', 'Proliferative TA', 'Malignancy continuum', 'CRC-like epithelial', 'Mature enterocyte', 'Goblet / immature goblet', 'preCAF-like', 'Treg-like', 'Stromal', 'Myeloid', 'T/NK']
    group_order = [g for g in group_order if g in set(summary['group'])]
    sig_order = [s for s in sig_order if s in set(summary['signature'])]
    mat = summary.pivot(index='signature', columns='group', values='mean_score').reindex(sig_order)[group_order].fillna(0)
    scaled = mat.copy()
    for sig in scaled.index:
        scaled.loc[sig] = landscape_minmax01(scaled.loc[sig]).values
    cmap = landscape_plt.cm.RdYlBu_r
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

##### 图2D #####
def landscape_panel_key_gene_heatmap(ax, zmat):
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
landscape_OUT_DIR = OUT / 'landscape' / 'OUT_DIR'
landscape_OUT_DIR.mkdir(parents=True, exist_ok=True)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    adata = landscape_load_data()
    summary, _ = landscape_compute_signature_summary(adata)
    _, zmat = landscape_compute_key_gene_group_matrix(adata)
    import matplotlib.pyplot as plt
    for label, fn in [('A', lambda ax: landscape_panel_umap_stage(ax, adata)), ('B', lambda ax: landscape_panel_umap_lineage(ax, adata)), ('C', lambda ax: landscape_panel_signature_dotplot(ax, summary)), ('D', lambda ax: landscape_panel_key_gene_heatmap(ax, zmat))]:

        def draw(label=label, fn=fn):
            fig, ax = plt.subplots(figsize=(12, 8))
            fn(ax)
            path = OUT / (label + '.png')
            fig.savefig(path, dpi=300, bbox_inches='tight')
            fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
            plt.close(fig)
            panel_file(label, path)
        attempt(label, draw)
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
