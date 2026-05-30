# =============================================================================
# TUTORIAL 1: Preprocesamiento y Clustering con scanpy — PBMC3k (versión 2021)
# Fuente: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
#         (workflow clásico "Preprocessing and clustering 3k PBMCs")
#
# Dataset: 2,700 PBMCs de 10X Genomics (sc.datasets.pbmc3k)
# Salida:  - adata_clustered.h5ad       (objeto procesado con clusters)
#          - pbmc3k_for_scenic.loom      (counts crudos -> entrada del tutorial 2)
#
# Se ejecuta en el entorno conda `scenic-2021` (ver environment.yml), el mismo
# que usa el tutorial de SCENIC, para que ambos corran sobre las MISMAS células.
# =============================================================================

from __future__ import annotations

import os
import sys
import time
import platform
import traceback
from datetime import datetime
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import anndata as ad
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

_report_lines  = []
_step_times    = {}
_start_time    = time.time()


def rlog(msg: str = "") -> None:
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


def write_report(adata=None, scenic_loom=None, error: str = None) -> None:
    total = time.time() - _start_time
    lines = []

    lines += [
        "# Reporte Técnico — scanpy clustering (PBMC3k, tutorial 2021)",
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
            f"| Genes (vars, tras HVG) | {adata.n_vars:,} |",
            f"| Genes totales (raw) | {adata.raw.n_vars if adata.raw is not None else adata.n_vars:,} |",
            f"| Obsm (embeddings) | {list(adata.obsm.keys())} |",
            f"| Obs columnas | {list(adata.obs.columns)} |",
            "",
        ]

        for col in adata.obs.columns:
            if col == "leiden":
                counts = adata.obs[col].value_counts()
                lines.append(f"\n**Distribución {col}:**")
                try:
                    items = sorted(counts.items(), key=lambda x: int(x[0]))
                except (ValueError, TypeError):
                    items = sorted(counts.items(), key=lambda x: str(x[0]))
                for cluster, n in items:
                    lines.append(f"- Cluster {cluster}: {n:,} células")

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
        celulas_ok   = "✓" if 2_000 <= adata.n_obs <= 3_000 else "✗ REVISAR"
        umap_ok      = "✓" if "X_umap" in adata.obsm else "✗ FALTA"
        pca_ok       = "✓" if "X_pca" in adata.obsm else "✗ FALTA"
        leiden_ok    = "✓" if "leiden" in adata.obs.columns else "✗ FALTA"
        celltype_ok  = "✓" if "cell_type" in adata.obs.columns else "✗ FALTA"
        loom_ok      = "✓" if scenic_loom and Path(scenic_loom).exists() else "✗ FALTA"

        lines += [
            f"| Células tras filtrado | 2,000–3,000 | {adata.n_obs:,} | {celulas_ok} |",
            f"| PCA calculado | Sí | {'Sí' if 'X_pca' in adata.obsm else 'No'} | {pca_ok} |",
            f"| UMAP calculado | Sí | {'Sí' if 'X_umap' in adata.obsm else 'No'} | {umap_ok} |",
            f"| Clusters Leiden | Sí | {'Sí' if 'leiden' in adata.obs.columns else 'No'} | {leiden_ok} |",
            f"| Tipos celulares | Sí | {'Sí' if 'cell_type' in adata.obs.columns else 'No'} | {celltype_ok} |",
            f"| Loom para SCENIC | Sí | {'Sí' if scenic_loom and Path(scenic_loom).exists() else 'No'} | {loom_ok} |",
        ]
    else:
        lines.append("| Pipeline | Completado | NO COMPLETADO | ✗ REVISAR |")

    if error:
        lines += ["", "---", "", "## Error capturado", "", "```", error, "```"]

    lines += ["", "---", "", "## Archivos generados", ""]
    for f in sorted(OUT_DIR.rglob("*")):
        if f.is_file() and f.suffix != ".md":
            size_kb = f.stat().st_size / 1024
            lines.append(f"- `{f.relative_to(OUT_DIR)}` ({size_kb:.1f} KB)")

    report_path = OUT_DIR / "reporte_tecnico.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nReporte guardado: {report_path}")


# =============================================================================
# Configuración scanpy
# =============================================================================
sc.settings.figdir = str(FIG_DIR)
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.verbosity = 1

adata        = None
scenic_loom  = None
error        = None

try:
    # ==========================================================================
    # PASO 1: CARGA DE DATOS (PBMC3k)
    # ==========================================================================
    step_start("PASO 1 — Carga de datos (PBMC3k)")

    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    rlog(f"- Dataset PBMC3k: {adata.n_obs:,} células × {adata.n_vars:,} genes")

    step_end("PASO 1 — Carga de datos (PBMC3k)")

    # ==========================================================================
    # PASO 2: FILTRADO BÁSICO Y CONTROL DE CALIDAD
    # ==========================================================================
    step_start("PASO 2 — Filtrado y control de calidad")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    rlog(f"- Tras filtro básico: {adata.n_obs:,} células × {adata.n_vars:,} genes")

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, multi_panel=True, show=False, save="_qc_metrics.png",
    )

    n_before = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
    rlog(f"- Filtro QC (n_genes<2500, pct_mt<5): {n_before:,} → {adata.n_obs:,} células")

    step_end("PASO 2 — Filtrado y control de calidad")

    # ==========================================================================
    # PASO 3: EXPORTAR LOOM DE COUNTS CRUDOS PARA SCENIC
    # ==========================================================================
    # SCENIC (tutorial 2) necesita la matriz de counts crudos con TODOS los genes
    # filtrados y los nombres de gen en símbolos HGNC. Se exporta ANTES de
    # normalizar para que GRNBoost2/AUCell trabajen sobre los datos correctos.
    step_start("PASO 3 — Exportar loom para SCENIC")

    adata.layers["counts"] = adata.X.copy()

    counts = adata.X
    if sp.issparse(counts):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float32)

    scenic_loom = OUT_DIR / "pbmc3k_for_scenic.loom"
    loompy.create(
        str(scenic_loom),
        counts.T,                                   # loom: genes × células
        {"Gene": np.array(adata.var_names)},
        {"CellID": np.array(adata.obs_names)},
    )
    size_mb = scenic_loom.stat().st_size / (1024 * 1024)
    rlog(f"- Loom exportado: {scenic_loom.name} ({size_mb:.1f} MB)")
    rlog(f"- Contenido: {adata.n_obs:,} células × {adata.n_vars:,} genes (counts crudos)")

    step_end("PASO 3 — Exportar loom para SCENIC")

    # ==========================================================================
    # PASO 4: NORMALIZACIÓN
    # ==========================================================================
    step_start("PASO 4 — Normalización")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    rlog("- Normalización total (target sum = 10,000) + log1p aplicadas")

    step_end("PASO 4 — Normalización")

    # ==========================================================================
    # PASO 5: GENES ALTAMENTE VARIABLES
    # ==========================================================================
    step_start("PASO 5 — Selección de genes variables")

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, show=False, save="_hvg.png")
    rlog(f"- Genes altamente variables: {int(adata.var.highly_variable.sum()):,}")

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    rlog(f"- AnnData reducido a HVG: {adata.n_obs:,} células × {adata.n_vars:,} genes")
    rlog("- Matriz completa preservada en adata.raw")

    step_end("PASO 5 — Selección de genes variables")

    # ==========================================================================
    # PASO 6: REGRESIÓN, ESCALADO Y PCA
    # ==========================================================================
    step_start("PASO 6 — Regresión, escalado y PCA")

    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False, save="_variance_ratio.png")
    var_explained = adata.uns["pca"]["variance_ratio"][:10].sum() * 100
    rlog(f"- Varianza explicada por primeros 10 PCs: {var_explained:.1f}%")

    step_end("PASO 6 — Regresión, escalado y PCA")

    # ==========================================================================
    # PASO 7: GRAFO DE VECINOS + UMAP
    # ==========================================================================
    step_start("PASO 7 — Vecinos + UMAP")

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    rlog("- Grafo de vecinos (n_neighbors=10, n_pcs=40) + UMAP calculados")

    step_end("PASO 7 — Vecinos + UMAP")

    # ==========================================================================
    # PASO 8: CLUSTERING LEIDEN
    # ==========================================================================
    step_start("PASO 8 — Clustering Leiden")

    sc.tl.leiden(adata, random_state=RANDOM_SEED)
    n_clusters = adata.obs["leiden"].nunique()
    rlog(f"- Clusters Leiden encontrados: {n_clusters}")

    sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False, save="_leiden.png")

    step_end("PASO 8 — Clustering Leiden")

    # ==========================================================================
    # PASO 9: EXPRESIÓN DIFERENCIAL Y ANOTACIÓN
    # ==========================================================================
    step_start("PASO 9 — Expresión diferencial y anotación")

    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False, save="_ranked_genes.png")

    de_top = sc.get.rank_genes_groups_df(adata, group=None).head(25 * n_clusters)
    de_top.to_csv(str(OUT_DIR / "diff_expression_top.csv"), index=False)

    # Marcadores canónicos del tutorial PBMC3k
    marker_genes = [
        "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY",
        "NKG7", "FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP",
    ]
    marker_present = [g for g in marker_genes if g in adata.raw.var_names]
    sc.pl.dotplot(adata, marker_present, groupby="leiden", show=False, save="_marker_genes.png")

    # Anotación canónica del tutorial (8 clusters). Mapeo seguro: los clusters no
    # cubiertos quedan como "Unknown" para que el script no falle si Leiden
    # produce un número distinto de clusters según versión/semilla.
    new_cluster_names = {
        "0": "CD4 T",
        "1": "CD14 Monocytes",
        "2": "B",
        "3": "CD8 T",
        "4": "NK",
        "5": "FCGR3A Monocytes",
        "6": "Dendritic",
        "7": "Megakaryocytes",
    }
    adata.obs["cell_type"] = (
        adata.obs["leiden"].astype(str).map(new_cluster_names).fillna("Unknown").astype("category")
    )
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False, save="_cell_types.png")

    for ct, n in adata.obs["cell_type"].value_counts().items():
        rlog(f"- {ct}: {n:,} células")

    step_end("PASO 9 — Expresión diferencial y anotación")

    # ==========================================================================
    # GUARDAR AnnData
    # ==========================================================================
    step_start("GUARDADO — AnnData final")

    out_h5ad = OUT_DIR / "adata_clustered.h5ad"
    adata.write_h5ad(str(out_h5ad))
    size_mb = out_h5ad.stat().st_size / (1024 * 1024)
    rlog(f"- Archivo: {out_h5ad.name} ({size_mb:.1f} MB)")

    step_end("GUARDADO — AnnData final")

except Exception as e:
    error = traceback.format_exc()
    print(f"\n*** ERROR ***\n{error}")

finally:
    write_report(adata=adata, scenic_loom=scenic_loom, error=error)
    total = time.time() - _start_time
    print(f"\n{'='*60}")
    print(f"Carpeta de salida: {OUT_DIR.resolve()}")
    print(f"Duración total:    {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*60}")
