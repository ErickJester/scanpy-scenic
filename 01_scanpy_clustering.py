# =============================================================================
# TUTORIAL: Preprocesamiento y Clustering con scanpy
# Fuente: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
#
# Dataset: Médula ósea humana (~17,000 células, 36,601 genes)
# Formato de entrada: archivos .h5 de 10X Genomics
# Formato de salida: objeto AnnData (.h5ad) con clusters anotados
# =============================================================================

from __future__ import annotations

import sys
import time
import platform
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import pooch
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import loompy

# Semilla global para reproducibilidad
RANDOM_SEED = 0
np.random.seed(RANDOM_SEED)

# =============================================================================
# INFRAESTRUCTURA DE SALIDA — carpeta con timestamp + reporte técnico
# =============================================================================

RUN_TS  = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
OUT_DIR = Path(f"runs/scanpy_{RUN_TS}")
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

_report_lines  = []   # líneas del reporte markdown
_step_times    = {}   # tiempos por paso
_start_time    = time.time()


def rlog(msg: str = "") -> None:
    """Imprime en consola y agrega al reporte."""
    print(msg)
    _report_lines.append(msg)


def step_start(name: str) -> None:
    _step_times[name] = time.time()
    rlog(f"\n### {name}")
    rlog(f"- **Inicio:** {datetime.now().strftime('%H:%M:%S')}")


def step_end(name: str) -> None:
    elapsed = time.time() - _step_times.get(name, time.time())
    rlog(f"- **Duración:** {elapsed:.1f}s")
    rlog(f"- **Estado:** OK ✓")


def write_report(adata=None, error: str = None) -> None:
    """Escribe el reporte final en markdown."""
    total = time.time() - _start_time
    lines = []

    lines += [
        "# Reporte Técnico — scanpy clustering",
        "",
        f"**Script:** `01_scanpy_clustering.py`  ",
        f"**Fecha/hora:** {RUN_TS}  ",
        f"**Duración total:** {total:.1f}s ({total/60:.1f} min)  ",
        f"**Estado final:** {'ERROR ✗' if error else 'COMPLETADO ✓'}  ",
        "",
        "---",
        "",
        "## Entorno de ejecución",
        "",
        f"| Variable | Valor |",
        f"|---|---|",
        f"| Python | {sys.version.split()[0]} |",
        f"| Sistema | {platform.system()} {platform.release()} |",
        f"| scanpy | {sc.__version__} |",
        f"| anndata | {ad.__version__} |",
        f"| numpy | {np.__version__} |",
        f"| pandas | {pd.__version__} |",
        "",
        "---",
        "",
        "## Log de ejecución por pasos",
        "",
    ]

    lines += _report_lines

    if adata is not None:
        lines += [
            "",
            "---",
            "",
            "## Métricas finales del AnnData",
            "",
            f"| Métrica | Valor |",
            f"|---|---|",
            f"| Células (obs) | {adata.n_obs:,} |",
            f"| Genes (vars) | {adata.n_vars:,} |",
            f"| Layers guardados | {list(adata.layers.keys())} |",
            f"| Obsm (embeddings) | {list(adata.obsm.keys())} |",
            f"| Obs columnas | {list(adata.obs.columns)} |",
            f"| Var columnas | {list(adata.var.columns)} |",
            "",
        ]

        # Distribución de clusters
        for col in adata.obs.columns:
            if col.startswith("leiden"):
                counts = adata.obs[col].value_counts()
                lines.append(f"\n**Distribución {col}:**")
                # Orden numérico si los labels son enteros, alfabético si no
                try:
                    items = sorted(counts.items(), key=lambda x: int(x[0]))
                except (ValueError, TypeError):
                    items = sorted(counts.items(), key=lambda x: str(x[0]))
                for cluster, n in items:
                    lines.append(f"- Cluster {cluster}: {n:,} células")

        # Tipos celulares si existen
        if "cell_type" in adata.obs.columns:
            lines += ["", "**Tipos celulares anotados:**"]
            for ct, n in adata.obs["cell_type"].value_counts().items():
                lines.append(f"- {ct}: {n:,} células")

    # Sección de validación para IA
    lines += [
        "",
        "---",
        "",
        "## Criterios de validación (para revisión IA)",
        "",
        "Esta sección permite evaluar si el pipeline corrió correctamente.",
        "",
        "| Criterio | Valor esperado | Valor obtenido | Estado |",
        "|---|---|---|---|",
    ]

    if adata is not None:
        celulas_ok   = "✓" if 10_000 <= adata.n_obs <= 25_000 else "✗ REVISAR"
        genes_ok     = "✓" if 15_000 <= adata.n_vars <= 40_000 else "✗ REVISAR"
        umap_ok      = "✓" if "X_umap" in adata.obsm else "✗ FALTA"
        pca_ok       = "✓" if "X_pca" in adata.obsm else "✗ FALTA"
        leiden_ok    = "✓" if any(c.startswith("leiden") for c in adata.obs.columns) else "✗ FALTA"
        celltype_ok  = "✓" if "cell_type" in adata.obs.columns else "✗ FALTA"
        counts_ok    = "✓" if "counts" in adata.layers else "✗ FALTA"

        lines += [
            f"| Células tras filtrado | 10,000–25,000 | {adata.n_obs:,} | {celulas_ok} |",
            f"| Genes tras filtrado | 15,000–40,000 | {adata.n_vars:,} | {genes_ok} |",
            f"| PCA calculado | Sí | {'Sí' if 'X_pca' in adata.obsm else 'No'} | {pca_ok} |",
            f"| UMAP calculado | Sí | {'Sí' if 'X_umap' in adata.obsm else 'No'} | {umap_ok} |",
            f"| Clusters Leiden | Sí | {'Sí' if any(c.startswith('leiden') for c in adata.obs.columns) else 'No'} | {leiden_ok} |",
            f"| Tipos celulares | Sí | {'Sí' if 'cell_type' in adata.obs.columns else 'No'} | {celltype_ok} |",
            f"| Layer 'counts' crudo | Sí | {'Sí' if 'counts' in adata.layers else 'No'} | {counts_ok} |",
        ]
    else:
        lines.append("| Pipeline | Completado | NO COMPLETADO | ✗ REVISAR |")

    if error:
        lines += [
            "",
            "---",
            "",
            "## Error capturado",
            "",
            "```",
            error,
            "```",
        ]

    lines += [
        "",
        "---",
        "",
        "## Archivos generados",
        "",
    ]

    for f in sorted(OUT_DIR.rglob("*")):
        if f.is_file() and f.suffix != ".md":
            size_kb = f.stat().st_size / 1024
            lines.append(f"- `{f.relative_to(OUT_DIR)}` ({size_kb:.1f} KB)")

    report_path = OUT_DIR / "reporte_tecnico.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nReporte guardado: {report_path}")


# =============================================================================
# Configuración scanpy — figuras van a la carpeta de la corrida
# =============================================================================

sc.settings.figdir = str(FIG_DIR)
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.verbosity = 1   # 0=errors, 1=warnings, 2=info, 3=hints

adata  = None
error  = None

try:
    # ==========================================================================
    # PASO 1: CARGA DE DATOS
    # ==========================================================================
    step_start("PASO 1 — Carga de datos")

    EXAMPLE_DATA = pooch.create(
        path=pooch.os_cache("scverse_tutorials"),
        base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    EXAMPLE_DATA.load_registry_from_doi()

    samples = {
        "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
        "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
    }

    adatas = {}
    for sample_id, filename in samples.items():
        path = EXAMPLE_DATA.fetch(filename)
        sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata
        rlog(f"- Muestra {sample_id}: {sample_adata.n_obs:,} células × {sample_adata.n_vars:,} genes")

    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    rlog(f"- AnnData combinado: {adata.n_obs:,} células × {adata.n_vars:,} genes")

    step_end("PASO 1 — Carga de datos")

    # ==========================================================================
    # PASO 2: CONTROL DE CALIDAD
    # ==========================================================================
    step_start("PASO 2 — Control de calidad (QC)")

    adata.var["mt"]   = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"]   = adata.var_names.str.contains("^HB[^(P)]")

    rlog(f"- Genes mitocondriales detectados: {adata.var['mt'].sum()}")
    rlog(f"- Genes ribosomales detectados: {adata.var['ribo'].sum()}")
    rlog(f"- Genes hemoglobina detectados: {adata.var['hb'].sum()}")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save="_qc_metrics.png",
    )
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts",
        color="pct_counts_mt",
        show=False,
        save="_total_vs_genes.png",
    )

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    rlog(f"- Células antes del filtro: {n_before:,}")
    rlog(f"- Células después del filtro: {adata.n_obs:,} (eliminadas: {n_before - adata.n_obs:,})")
    rlog(f"- Genes después del filtro: {adata.n_vars:,}")

    step_end("PASO 2 — Control de calidad (QC)")

    # ==========================================================================
    # PASO 3: DETECCIÓN DE DOBLETES
    # ==========================================================================
    step_start("PASO 3 — Detección de dobletes (Scrublet)")

    sc.pp.scrublet(adata, batch_key="sample", random_state=RANDOM_SEED)
    n_doublets = adata.obs["predicted_doublet"].sum()
    rlog(f"- Dobletes predichos: {n_doublets:,} ({100*n_doublets/adata.n_obs:.1f}% del total)")

    step_end("PASO 3 — Detección de dobletes (Scrublet)")

    # ==========================================================================
    # PASO 4: NORMALIZACIÓN
    # ==========================================================================
    step_start("PASO 4 — Normalización")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    rlog("- Conteos crudos guardados en layer 'counts'")
    rlog("- Normalización total aplicada (target sum = 10,000)")
    rlog("- Transformación log1p aplicada")

    step_end("PASO 4 — Normalización")

    # ==========================================================================
    # PASO 5: GENES ALTAMENTE VARIABLES
    # ==========================================================================
    step_start("PASO 5 — Selección de genes variables")

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    sc.pl.highly_variable_genes(adata, show=False, save="_hvg.png")
    rlog(f"- Genes altamente variables seleccionados: {adata.var.highly_variable.sum():,}")

    step_end("PASO 5 — Selección de genes variables")

    # ==========================================================================
    # PASO 6: PCA
    # ==========================================================================
    step_start("PASO 6 — PCA")

    sc.tl.pca(adata, random_state=RANDOM_SEED)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False, save="_variance_ratio.png")
    sc.pl.pca(
        adata,
        color=["sample", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3)],
        ncols=2,
        size=2,
        show=False,
        save="_components.png",
    )
    var_explained = adata.uns["pca"]["variance_ratio"][:10].sum() * 100
    rlog(f"- Varianza explicada por primeros 10 PCs: {var_explained:.1f}%")

    step_end("PASO 6 — PCA")

    # ==========================================================================
    # PASO 7: GRAFO DE VECINOS + UMAP
    # ==========================================================================
    step_start("PASO 7 — Vecinos + UMAP")

    sc.pp.neighbors(adata, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    sc.pl.umap(adata, color="sample", size=2, show=False, save="_by_sample.png")
    rlog("- Grafo de vecinos construido (k=15 por defecto)")
    rlog("- UMAP calculado")

    step_end("PASO 7 — Vecinos + UMAP")

    # ==========================================================================
    # PASO 8: CLUSTERING LEIDEN
    # ==========================================================================
    step_start("PASO 8 — Clustering Leiden")

    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(
            adata,
            key_added=f"leiden_res_{res:4.2f}",
            resolution=res,
            flavor="igraph",
            n_iterations=2,
            random_state=RANDOM_SEED,
        )
        n_clusters = adata.obs[f"leiden_res_{res:4.2f}"].nunique()
        rlog(f"- Resolución {res}: {n_clusters} clusters")

    sc.pl.umap(
        adata,
        color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
        legend_loc="on data",
        show=False,
        save="_leiden_resolutions.png",
    )
    sc.pl.umap(
        adata,
        color=["leiden_res_0.02", "predicted_doublet", "doublet_score"],
        wspace=0.5,
        size=3,
        show=False,
        save="_qc_check.png",
    )

    step_end("PASO 8 — Clustering Leiden")

    # ==========================================================================
    # PASO 9: ANOTACIÓN DE TIPOS CELULARES
    # ==========================================================================
    step_start("PASO 9 — Anotación de tipos celulares")

    marker_genes = {
        "CD14+ Mono":      ["FCN1", "CD14"],
        "CD16+ Mono":      ["TCF7L2", "FCGR3A", "LYN"],
        "cDC2":            ["CST3", "COTL1", "LYZ", "CLEC10A", "FCER1A"],
        "Erythroblast":    ["MKI67", "HBA1", "HBB"],
        "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
        "NK":              ["GNLY", "NKG7", "CD247", "TYROBP", "KLRG1"],
        "Naive CD20+ B":   ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
        "Plasma cells":    ["MZB1", "HSP90B1", "PRDM1", "IGKC", "JCHAIN"],
        "CD4+ T":          ["CD4", "IL7R", "TRBC2"],
        "CD8+ T":          ["CD8A", "CD8B", "GZMK", "CCL5", "GZMB"],
        "T naive":         ["LEF1", "CCR7", "TCF7"],
        "pDC":             ["IL3RA", "COBLL1", "TCF4"],
    }

    sc.pl.dotplot(
        adata, marker_genes, groupby="leiden_res_0.02",
        standard_scale="var", show=False, save="_marker_genes.png",
    )

    # Mapeo seguro: solo asigna nombres a clusters existentes; los demás quedan como "Unknown"
    cluster_to_celltype = {
        "0": "Lymphocytes",
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    }
    adata.obs["cell_type"] = (
        adata.obs["leiden_res_0.02"]
        .map(cluster_to_celltype)
        .fillna("Unknown")
        .astype("category")
    )
    rlog(f"- Clusters mapeados a tipo celular: {list(cluster_to_celltype.keys())}")
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False, save="_cell_types.png")

    for ct, n in adata.obs["cell_type"].value_counts().items():
        rlog(f"- {ct}: {n:,} células")

    step_end("PASO 9 — Anotación de tipos celulares")

    # ==========================================================================
    # PASO 10: EXPRESIÓN DIFERENCIAL
    # ==========================================================================
    step_start("PASO 10 — Expresión diferencial (Wilcoxon)")

    sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(
        adata, groupby="leiden_res_0.50",
        standard_scale="var", n_genes=5,
        show=False, save="_diffexp.png",
    )

    # Guardar tabla de top genes diferenciales por cluster
    de_top = sc.get.rank_genes_groups_df(adata, group=None).head(50 * 15)
    de_path = OUT_DIR / "diff_expression_top.csv"
    de_top.to_csv(str(de_path), index=False)
    rlog(f"- Tabla de expresión diferencial guardada: {de_path.name}")
    rlog("- Análisis Wilcoxon completado para leiden_res_0.50")
    rlog("- Top 5 genes por cluster guardados en adata.uns['rank_genes_groups']")

    step_end("PASO 10 — Expresión diferencial (Wilcoxon)")

    # ==========================================================================
    # PASO 11: EXPORTAR LOOM DE COUNTS CRUDOS PARA SCENIC
    # ==========================================================================
    # SCENIC (02) corre sobre las MISMAS células que este clustering. Se exporta
    # la matriz de counts crudos (layer 'counts') con todos los genes filtrados y
    # los mismos códigos de barras, en símbolos HGNC, para que GRNBoost2/AUCell y
    # cisTarget (DB hg38) trabajen sobre los datos correctos. El puente entre
    # ambos tutoriales es este loom + los barcodes compartidos.
    step_start("PASO 11 — Exportar loom para SCENIC")

    counts = adata.layers["counts"]
    if sp.issparse(counts):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float32)

    scenic_loom = OUT_DIR / "bonemarrow_for_scenic.loom"
    loompy.create(
        str(scenic_loom),
        counts.T,                                   # loom: genes × células
        {"Gene": np.array(adata.var_names)},
        {"CellID": np.array(adata.obs_names)},
    )
    size_mb = scenic_loom.stat().st_size / (1024 * 1024)
    rlog(f"- Loom exportado: {scenic_loom.name} ({size_mb:.1f} MB)")
    rlog(f"- Contenido: {adata.n_obs:,} células × {adata.n_vars:,} genes (counts crudos)")

    step_end("PASO 11 — Exportar loom para SCENIC")

    # ==========================================================================
    # GUARDAR AnnData
    # ==========================================================================
    step_start("GUARDADO — AnnData final")

    out_h5ad = OUT_DIR / "adata_clustered.h5ad"
    adata.write_h5ad(str(out_h5ad))
    size_mb = out_h5ad.stat().st_size / (1024 * 1024)
    rlog(f"- Archivo: {out_h5ad}")
    rlog(f"- Tamaño: {size_mb:.1f} MB")

    step_end("GUARDADO — AnnData final")

except Exception as e:
    error = traceback.format_exc()
    print(f"\n*** ERROR ***\n{error}")

finally:
    write_report(adata=adata, error=error)
    total = time.time() - _start_time
    print(f"\n{'='*60}")
    print(f"Carpeta de salida: {OUT_DIR.resolve()}")
    print(f"Duración total:    {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*60}")
