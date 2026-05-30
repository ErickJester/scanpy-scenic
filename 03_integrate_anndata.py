# =============================================================================
# Integración scanpy + SCENIC en un único AnnData / loom
# Fuente: https://github.com/aertslab/SCENICprotocol
#
# Toma los resultados de 01_scanpy_clustering.py y 02_scenic_pipeline.py,
# combina la matriz AUC de regulones con el AnnData de clustering,
# y exporta un loom integrado listo para SCope.
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
import numpy as np
import pandas as pd
import scanpy as sc
import loompy

RANDOM_SEED = 0
np.random.seed(RANDOM_SEED)

# =============================================================================
# INFRAESTRUCTURA DE SALIDA
# =============================================================================

RUN_TS  = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
OUT_DIR = Path(f"runs/integrate_{RUN_TS}")
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

_report_lines = []
_step_times   = {}
_start_time   = time.time()
error         = None


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


def _find_latest_run(prefix: str) -> Path | None:
    """Busca la carpeta runs/<prefix>_* más reciente."""
    runs = sorted(Path("runs").glob(f"{prefix}_*"), reverse=True)
    return runs[0] if runs else None


def write_report(
    adata=None,
    auc_matrix=None,
    scanpy_dir=None,
    scenic_dir=None,
    n_matched=0,
    error: str = None,
) -> None:
    total = time.time() - _start_time
    lines = []

    lines += [
        "# Reporte Técnico — Integración scanpy + SCENIC",
        "",
        f"**Script:** `03_integrate_anndata.py`  ",
        f"**Fecha/hora:** {RUN_TS}  ",
        f"**Duración total:** {total:.1f}s ({total/60:.1f} min)  ",
        f"**Estado final:** {'ERROR ✗' if error else 'COMPLETADO ✓'}  ",
        "",
        "---",
        "",
        "## Entorno de ejecución",
        "",
        "| Variable | Valor |",
        "|---|---|",
        f"| Python | {sys.version.split()[0]} |",
        f"| Sistema | {platform.system()} {platform.release()} |",
        f"| scanpy | {sc.__version__} |",
        f"| anndata | {ad.__version__} |",
        f"| numpy | {np.__version__} |",
        f"| pandas | {pd.__version__} |",
        "",
        "## Fuentes de datos",
        "",
        f"| Fuente | Carpeta |",
        f"|---|---|",
        f"| scanpy | `{scanpy_dir}` |",
        f"| SCENIC | `{scenic_dir}` |",
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
            "## Métricas del AnnData integrado",
            "",
            "| Métrica | Valor |",
            "|---|---|",
            f"| Células (obs) | {adata.n_obs:,} |",
            f"| Genes (vars) | {adata.n_vars:,} |",
            f"| Layers | {list(adata.layers.keys())} |",
            f"| Obsm | {list(adata.obsm.keys())} |",
        ]

        if auc_matrix is not None:
            lines += [
                f"| Regulones integrados | {auc_matrix.shape[1]:,} |",
                f"| Células con scores AUC | {n_matched:,} |",
            ]

        if "cell_type" in adata.obs.columns:
            lines += ["", "**Tipos celulares:**"]
            for ct, n in adata.obs["cell_type"].value_counts().items():
                lines.append(f"- {ct}: {n:,} células")

    # Validación
    lines += [
        "",
        "---",
        "",
        "## Criterios de validación (para revisión IA)",
        "",
        "| Criterio | Valor esperado | Valor obtenido | Estado |",
        "|---|---|---|---|",
    ]

    if adata is not None:
        has_scenic = "X_scenic" in adata.obsm
        has_umap   = "X_umap" in adata.obsm
        has_loom   = (OUT_DIR / "integrated_output.loom").exists()
        has_h5ad   = (OUT_DIR / "adata_integrated.h5ad").exists()

        lines += [
            f"| AnnData cargado | Sí | Sí | ✓ |",
            f"| X_scenic en obsm | Sí | {'Sí' if has_scenic else 'No'} | {'✓' if has_scenic else '✗'} |",
            f"| UMAP presente | Sí | {'Sí' if has_umap else 'No'} | {'✓' if has_umap else '✗'} |",
            f"| Loom integrado | Sí | {'Sí' if has_loom else 'No'} | {'✓' if has_loom else '✗'} |",
            f"| H5AD integrado | Sí | {'Sí' if has_h5ad else 'No'} | {'✓' if has_h5ad else '✗'} |",
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

adata      = None
auc_matrix = None
scanpy_dir = None
scenic_dir = None
n_matched  = 0

try:
    # ==========================================================================
    # PASO 1: LOCALIZAR CORRIDAS MÁS RECIENTES
    # ==========================================================================
    step_start("PASO 1 — Localizar corridas")

    scanpy_dir = _find_latest_run("scanpy")
    scenic_dir = _find_latest_run("scenic")

    if scanpy_dir is None:
        raise FileNotFoundError("No se encontró ninguna carpeta runs/scanpy_*. Ejecuta primero 01_scanpy_clustering.py")
    if scenic_dir is None:
        raise FileNotFoundError("No se encontró ninguna carpeta runs/scenic_*. Ejecuta primero 02_scenic_pipeline.py")

    scanpy_h5ad = scanpy_dir / "adata_clustered.h5ad"
    scenic_auc  = scenic_dir / "auc_matrix.csv"

    if not scanpy_h5ad.exists():
        raise FileNotFoundError(f"No se encontró {scanpy_h5ad}")
    if not scenic_auc.exists():
        raise FileNotFoundError(f"No se encontró {scenic_auc}")

    rlog(f"- Scanpy: `{scanpy_dir.name}`")
    rlog(f"- SCENIC: `{scenic_dir.name}`")

    step_end("PASO 1 — Localizar corridas")

    # ==========================================================================
    # PASO 2: CARGAR ANNDATA DE SCANPY
    # ==========================================================================
    step_start("PASO 2 — Cargar AnnData de scanpy")

    adata = sc.read_h5ad(str(scanpy_h5ad))
    rlog(f"- Células: {adata.n_obs:,}")
    rlog(f"- Genes: {adata.n_vars:,}")
    rlog(f"- Obsm: {list(adata.obsm.keys())}")
    rlog(f"- Clusters Leiden: {[c for c in adata.obs.columns if c.startswith('leiden')]}")

    if "cell_type" in adata.obs.columns:
        rlog(f"- Tipos celulares: {adata.obs['cell_type'].nunique()}")

    step_end("PASO 2 — Cargar AnnData de scanpy")

    # ==========================================================================
    # PASO 3: CARGAR MATRIZ AUC DE SCENIC
    # ==========================================================================
    step_start("PASO 3 — Cargar matriz AUC de SCENIC")

    auc_matrix = pd.read_csv(str(scenic_auc), index_col=0)
    rlog(f"- Matriz AUC: {auc_matrix.shape[0]:,} células × {auc_matrix.shape[1]:,} regulones")
    rlog(f"- Regulones: {list(auc_matrix.columns)}")
    rlog(f"- AUC media global: {auc_matrix.values.mean():.4f}")

    step_end("PASO 3 — Cargar matriz AUC de SCENIC")

    # ==========================================================================
    # PASO 4: INTEGRAR AUC EN ANNDATA
    # ==========================================================================
    step_start("PASO 4 — Integrar AUC en AnnData")

    # Intersectar células entre scanpy y SCENIC
    common_cells = adata.obs_names.intersection(auc_matrix.index)
    n_matched = len(common_cells)

    rlog(f"- Células en scanpy: {adata.n_obs:,}")
    rlog(f"- Células en SCENIC: {auc_matrix.shape[0]:,}")
    rlog(f"- Células en común: {n_matched:,}")

    if n_matched == 0:
        rlog("- ⚠️ No hay células en común — los pipelines usaron datasets distintos")
        rlog("  → Insertando scores AUC como NaN para todas las células de scanpy")
        rlog("  → (Cuando ambos pipelines usen los mismos datos, la integración será completa)")
        auc_aligned = pd.DataFrame(
            np.nan,
            index=adata.obs_names,
            columns=auc_matrix.columns,
        )
    else:
        rlog(f"- Cobertura: {100 * n_matched / adata.n_obs:.1f}% de las células de scanpy")
        auc_aligned = auc_matrix.reindex(adata.obs_names)

    adata.obsm["X_scenic"] = auc_aligned.values
    adata.uns["scenic_regulon_names"] = list(auc_aligned.columns)
    rlog(f"- Agregado adata.obsm['X_scenic'] con shape {auc_aligned.shape}")

    # Guardar scores por regulón individualmente en obs (para UMAP)
    for col in auc_aligned.columns:
        safe_name = col.replace("(", "_").replace(")", "").replace("+", "plus")
        adata.obs[f"scenic_{safe_name}"] = auc_aligned[col].values

    rlog(f"- Agregados {len(auc_aligned.columns)} scores individuales en adata.obs")

    step_end("PASO 4 — Integrar AUC en AnnData")

    # ==========================================================================
    # PASO 5: VISUALIZACIÓN
    # ==========================================================================
    step_start("PASO 5 — Visualización integrada")

    if n_matched > 0 and "X_umap" in adata.obsm:
        # Top 4 regulones por actividad media
        top_regs = auc_matrix[common_cells].mean().sort_values(ascending=False).head(4)
        top_obs_names = [
            f"scenic_{r.replace('(', '_').replace(')', '').replace('+', 'plus')}"
            for r in top_regs.index
        ]

        sc.pl.umap(
            adata,
            color=top_obs_names,
            ncols=2,
            title=[f"Regulón: {r}" for r in top_regs.index],
            show=False,
            save="_regulon_activity.png",
        )
        rlog(f"- UMAP de top regulones: {list(top_regs.index)}")

        # Heatmap: actividad media de regulones por tipo celular
        if "cell_type" in adata.obs.columns:
            scenic_cols = [c for c in adata.obs.columns if c.startswith("scenic_")]
            mean_by_ct = adata.obs.groupby("cell_type")[scenic_cols].mean()
            mean_by_ct.columns = [c.replace("scenic_", "") for c in mean_by_ct.columns]

            import matplotlib.pyplot as plt
            import seaborn as sns

            fig, ax = plt.subplots(figsize=(max(8, len(mean_by_ct.columns) * 0.6), 5))
            sns.heatmap(
                mean_by_ct,
                cmap="viridis",
                annot=True,
                fmt=".3f",
                ax=ax,
                linewidths=0.5,
            )
            ax.set_title("Actividad media de regulones por tipo celular")
            ax.set_ylabel("Tipo celular")
            ax.set_xlabel("Regulón")
            plt.tight_layout()
            heatmap_path = FIG_DIR / "heatmap_regulons_by_celltype.png"
            fig.savefig(str(heatmap_path), dpi=80, bbox_inches="tight")
            plt.close(fig)
            rlog(f"- Heatmap regulones × tipo celular guardado")

        # UMAP combinado: tipo celular + top regulón
        if "cell_type" in adata.obs.columns:
            sc.pl.umap(
                adata,
                color=["cell_type", top_obs_names[0]],
                ncols=2,
                title=["Tipo celular", f"Top regulón: {top_regs.index[0]}"],
                show=False,
                save="_celltype_vs_regulon.png",
            )
            rlog("- UMAP comparativo: tipo celular vs top regulón")

    else:
        rlog("- ⚠️ Saltando visualizaciones (no hay células en común o falta UMAP)")

    step_end("PASO 5 — Visualización integrada")

    # ==========================================================================
    # PASO 6: EXPORTAR LOOM INTEGRADO
    # ==========================================================================
    step_start("PASO 6 — Exportar loom integrado")

    loom_path = OUT_DIR / "integrated_output.loom"

    # Preparar datos para loom: expresión + metadatos
    row_attrs = {"Gene": np.array(adata.var_names)}

    col_attrs = {
        "CellID":   np.array(adata.obs_names),
        "nGene":    np.array(adata.obs["n_genes_by_counts"].values, dtype=np.int32)
                    if "n_genes_by_counts" in adata.obs.columns
                    else np.zeros(adata.n_obs, dtype=np.int32),
        "nUMI":     np.array(adata.obs["total_counts"].values, dtype=np.float32)
                    if "total_counts" in adata.obs.columns
                    else np.zeros(adata.n_obs, dtype=np.float32),
    }

    # Agregar clusters y tipo celular
    for col in adata.obs.columns:
        if col.startswith("leiden"):
            col_attrs[col] = np.array(adata.obs[col].astype(str))
    if "cell_type" in adata.obs.columns:
        col_attrs["cell_type"] = np.array(adata.obs["cell_type"].astype(str))

    # Agregar embeddings
    if "X_umap" in adata.obsm:
        col_attrs["UMAP_1"] = np.array(adata.obsm["X_umap"][:, 0], dtype=np.float32)
        col_attrs["UMAP_2"] = np.array(adata.obsm["X_umap"][:, 1], dtype=np.float32)

    # Agregar scores de regulones
    for col in adata.obs.columns:
        if col.startswith("scenic_"):
            vals = adata.obs[col].values
            col_attrs[col] = np.array(
                np.nan_to_num(vals, nan=0.0), dtype=np.float32
            )

    # Matriz de expresión (usar counts si existe, sino X)
    if "counts" in adata.layers:
        import scipy.sparse as sp
        X = adata.layers["counts"]
        if sp.issparse(X):
            X = X.toarray()
        expr = np.array(X.T, dtype=np.float32)
    else:
        import scipy.sparse as sp
        X = adata.X
        if sp.issparse(X):
            X = X.toarray()
        expr = np.array(X.T, dtype=np.float32)

    loompy.create(str(loom_path), expr, row_attrs, col_attrs)
    size_mb = loom_path.stat().st_size / (1024 * 1024)
    rlog(f"- Loom integrado: {loom_path.name} ({size_mb:.1f} MB)")
    rlog(f"- Atributos de fila: {list(row_attrs.keys())}")
    rlog(f"- Atributos de columna: {list(col_attrs.keys())}")

    step_end("PASO 6 — Exportar loom integrado")

    # ==========================================================================
    # PASO 7: GUARDAR ANNDATA INTEGRADO
    # ==========================================================================
    step_start("PASO 7 — Guardar AnnData integrado")

    h5ad_path = OUT_DIR / "adata_integrated.h5ad"
    adata.write_h5ad(str(h5ad_path))
    size_mb = h5ad_path.stat().st_size / (1024 * 1024)
    rlog(f"- Archivo: {h5ad_path.name} ({size_mb:.1f} MB)")
    rlog(f"- obsm keys: {list(adata.obsm.keys())}")

    step_end("PASO 7 — Guardar AnnData integrado")

except Exception as e:
    error = traceback.format_exc()
    print(f"\n*** ERROR ***\n{error}")

finally:
    write_report(
        adata=adata,
        auc_matrix=auc_matrix,
        scanpy_dir=scanpy_dir,
        scenic_dir=scenic_dir,
        n_matched=n_matched,
        error=error,
    )
    total = time.time() - _start_time
    print(f"\n{'='*60}")
    print(f"Carpeta de salida: {OUT_DIR.resolve()}")
    print(f"Duración total:    {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*60}")
