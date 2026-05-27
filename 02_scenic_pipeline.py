# =============================================================================
# TUTORIAL: Inferencia de Redes Regulatorias con pySCENIC
# Fuente: https://github.com/aertslab/SCENICprotocol
#
# Dataset de prueba: expr_mat_tiny.loom (~78 MB con bases de datos incluidas)
# Flujo: GRNBoost2 -> cisTarget -> AUCell -> visualización
# =============================================================================

import sys
import time
import platform
import traceback
import urllib.request
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
import scanpy as sc
import loompy

from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.export import export2loom

import pyscenic

# Semilla global para reproducibilidad
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# =============================================================================
# INFRAESTRUCTURA DE SALIDA — carpeta con timestamp + reporte técnico
# =============================================================================

RUN_TS  = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
OUT_DIR = Path(f"runs/scenic_{RUN_TS}")
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


def write_report(
    ex_matrix=None,
    tf_names=None,
    modules=None,
    regulons=None,
    auc_matrix=None,
    adata_scenic=None,
    error: str = None,
) -> None:
    total = time.time() - _start_time
    lines = []

    lines += [
        "# Reporte Técnico — SCENIC pipeline",
        "",
        f"**Script:** `02_scenic_pipeline.py`  ",
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
        f"| pyscenic | {pyscenic.__version__} |",
        f"| pandas | {pd.__version__} |",
        f"| numpy | {np.__version__} |",
        "",
        "---",
        "",
        "## Log de ejecución por pasos",
        "",
    ]

    lines += _report_lines

    # Métricas del pipeline
    lines += ["", "---", "", "## Métricas del pipeline SCENIC", ""]

    if ex_matrix is not None:
        lines += [
            "| Etapa | Métrica | Valor |",
            "|---|---|---|",
            f"| Datos de entrada | Células | {ex_matrix.shape[0]:,} |",
            f"| Datos de entrada | Genes | {ex_matrix.shape[1]:,} |",
        ]
        if tf_names:
            lines.append(f"| GRNBoost2 | TFs evaluados | {len(tf_names):,} |")
        if modules:
            lines.append(f"| GRNBoost2 | Módulos generados | {len(modules):,} |")
        if regulons:
            lines.append(f"| cisTarget | Regulones finales | {len(regulons):,} |")
            avg_targets = np.mean([len(r.genes) for r in regulons])
            lines.append(f"| cisTarget | Genes target promedio por regulón | {avg_targets:.1f} |")
        if auc_matrix is not None:
            lines.append(f"| AUCell | Células × Regulones | {auc_matrix.shape[0]:,} × {auc_matrix.shape[1]:,} |")
        if adata_scenic is not None:
            lines.append(f"| Clustering | Clusters Leiden | {adata_scenic.obs['leiden'].nunique()} |")

    # Top regulones por actividad media
    if auc_matrix is not None and len(auc_matrix) > 0:
        top_regulons = auc_matrix.mean().sort_values(ascending=False).head(10)
        lines += [
            "",
            "**Top 10 regulones por actividad media:**",
            "",
            "| Regulón | AUC medio |",
            "|---|---|",
        ]
        for reg, val in top_regulons.items():
            lines.append(f"| {reg} | {val:.4f} |")

    # Sección de validación para IA
    lines += [
        "",
        "---",
        "",
        "## Criterios de validación (para revisión IA)",
        "",
        "Esta sección permite evaluar si el pipeline SCENIC corrió correctamente.",
        "",
        "| Criterio | Valor esperado | Valor obtenido | Estado |",
        "|---|---|---|---|",
    ]

    if ex_matrix is not None:
        celulas_ok  = "✓" if ex_matrix.shape[0] > 0 else "✗ VACÍO"
        genes_ok    = "✓" if ex_matrix.shape[1] > 100 else "✗ MUY POCOS"
        modulos_ok  = "✓" if modules and len(modules) > 0 else "✗ FALTA"
        regulons_ok = "✓" if regulons and len(regulons) > 0 else "✗ FALTA — posible problema con cisTarget o DB"
        aucell_ok   = "✓" if auc_matrix is not None and not auc_matrix.empty else "✗ FALTA"
        umap_ok     = "✓" if adata_scenic is not None and "X_umap" in adata_scenic.obsm else "✗ FALTA"
        loom_ok     = "✓" if (OUT_DIR / "scenic_output.loom").exists() else "✗ NO GENERADO"

        lines += [
            f"| Matriz de expresión cargada | > 0 células | {ex_matrix.shape[0]:,} | {celulas_ok} |",
            f"| Genes en la matriz | > 100 | {ex_matrix.shape[1]:,} | {genes_ok} |",
            f"| Módulos GRN generados | > 0 | {len(modules) if modules else 0:,} | {modulos_ok} |",
            f"| Regulones cisTarget | > 0 | {len(regulons) if regulons else 0:,} | {regulons_ok} |",
            f"| Matriz AUCell calculada | Sí | {'Sí' if auc_matrix is not None else 'No'} | {aucell_ok} |",
            f"| UMAP sobre regulones | Sí | {'Sí' if adata_scenic is not None and 'X_umap' in adata_scenic.obsm else 'No'} | {umap_ok} |",
            f"| Loom exportado | Sí | {'Sí' if (OUT_DIR / 'scenic_output.loom').exists() else 'No'} | {loom_ok} |",
        ]

        if regulons and len(regulons) == 0:
            lines += [
                "",
                "> ⚠️ **ALERTA IA:** Si regulones = 0, verificar:",
                "> 1. Que el archivo `rankings.feather` corresponde al genoma correcto (hg38 para datos humanos)",
                "> 2. Que `motifs.tbl` es compatible con la versión de pySCENIC",
                "> 3. Que los nombres de genes en la matriz coinciden con los de la DB",
            ]
    else:
        lines.append("| Pipeline SCENIC | Completado | NO COMPLETADO | ✗ REVISAR |")

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
            "",
            "> ⚠️ **Para IA:** Analizar el traceback anterior e identificar en qué paso falló.",
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
# Configuración scanpy
# =============================================================================
sc.settings.figdir = str(FIG_DIR)
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.verbosity = 1

# Variables del pipeline
ex_matrix   = None
tf_names    = None
modules     = None
regulons    = None
auc_matrix  = None
adata_scenic = None

try:
    # ==========================================================================
    # PASO 0: DESCARGAR DATOS DE PRUEBA
    # ==========================================================================
    step_start("PASO 0 — Descarga de datos de prueba")

    DATA_DIR = Path("scenic_data")
    DATA_DIR.mkdir(exist_ok=True)

    URLS = {
        "expr_mat_tiny.loom": (
            "https://resources.aertslab.org/cistarget/tmp/pyscenic_tutorial/"
            "expr_mat_tiny.loom"
        ),
        "test_TFs_tiny.txt": (
            "https://resources.aertslab.org/cistarget/tmp/pyscenic_tutorial/"
            "test_TFs_tiny.txt"
        ),
        "motifs.tbl": (
            "https://resources.aertslab.org/cistarget/motif2tf/"
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        ),
        "rankings.feather": (
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/"
            "hg38/refseq_r80/mc_v10_clust/gene_based/"
            "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        ),
    }

    for filename, url in URLS.items():
        dest = DATA_DIR / filename
        if not dest.exists():
            rlog(f"- Descargando {filename}...")
            urllib.request.urlretrieve(url, str(dest))
            size_mb = dest.stat().st_size / (1024 * 1024)
            rlog(f"  → Guardado ({size_mb:.1f} MB)")
        else:
            size_mb = dest.stat().st_size / (1024 * 1024)
            rlog(f"- Ya existe: {filename} ({size_mb:.1f} MB)")

    LOOM_PATH    = DATA_DIR / "expr_mat_tiny.loom"
    TF_PATH      = DATA_DIR / "test_TFs_tiny.txt"
    MOTIFS_PATH  = DATA_DIR / "motifs.tbl"
    RANKING_PATH = DATA_DIR / "rankings.feather"

    step_end("PASO 0 — Descarga de datos de prueba")

    # ==========================================================================
    # PASO 1: CARGAR MATRIZ DE EXPRESIÓN
    # ==========================================================================
    step_start("PASO 1 — Carga de la matriz de expresión")

    with loompy.connect(str(LOOM_PATH)) as ds:
        ex_matrix = pd.DataFrame(
            data=ds[:, :].T,
            index=ds.ca["CellID"],
            columns=ds.ra["Gene"],
        )

    rlog(f"- Matriz cargada: {ex_matrix.shape[0]:,} células × {ex_matrix.shape[1]:,} genes")
    rlog(f"- Rango de expresión: min={ex_matrix.values.min():.2f}, max={ex_matrix.values.max():.2f}")
    rlog(f"- % de ceros: {(ex_matrix.values == 0).mean() * 100:.1f}%")

    tf_names = pd.read_csv(str(TF_PATH), header=None).squeeze().tolist()
    rlog(f"- Factores de transcripción cargados: {len(tf_names):,}")

    step_end("PASO 1 — Carga de la matriz de expresión")

    # ==========================================================================
    # PASO 2: GRNBoost2
    # ==========================================================================
    step_start("PASO 2 — Inferencia de red (GRNBoost2)")

    adjacencies = grnboost2(
        expression_data=ex_matrix,
        tf_names=tf_names,
        verbose=True,
        seed=RANDOM_SEED,
    )

    rlog(f"- Relaciones TF–gen encontradas: {len(adjacencies):,}")
    rlog(f"- Importancia máxima: {adjacencies['importance'].max():.4f}")
    rlog(f"- Importancia media: {adjacencies['importance'].mean():.4f}")
    rlog(f"- Top 5 TFs más conectados:")
    top_tfs = adjacencies.groupby("TF").size().sort_values(ascending=False).head(5)
    for tf, n in top_tfs.items():
        rlog(f"  • {tf}: {n} genes target")

    # Guardar adjacencias
    adj_path = OUT_DIR / "adjacencies.csv"
    adjacencies.to_csv(str(adj_path), index=False)
    rlog(f"- Adjacencias guardadas: {adj_path}")

    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    rlog(f"- Módulos de co-regulación: {len(modules):,}")

    step_end("PASO 2 — Inferencia de red (GRNBoost2)")

    # ==========================================================================
    # PASO 3: cisTarget
    # ==========================================================================
    step_start("PASO 3 — Pruning por motivos (cisTarget)")

    db = RankingDatabase(fname=str(RANKING_PATH), name="hg38_rankings")

    df_regulons = prune2df(
        rnkdbs=[db],
        modules=modules,
        motif_annotations_fname=str(MOTIFS_PATH),
        num_workers=4,
    )

    regulons = df2regulons(df_regulons)
    rlog(f"- Regulones identificados: {len(regulons):,}")

    if regulons:
        sizes = [len(r.genes) for r in regulons]
        rlog(f"- Genes target por regulón: min={min(sizes)}, max={max(sizes)}, promedio={np.mean(sizes):.1f}")
        rlog("- Top 5 regulones por número de genes target:")
        for r in sorted(regulons, key=lambda x: len(x.genes), reverse=True)[:5]:
            rlog(f"  • {r.transcription_factor}: {len(r.genes)} genes")
    else:
        rlog("- ⚠️ No se identificaron regulones — revisar compatibilidad de DB de motivos con la matriz")

    # Guardar tabla de regulones
    if not df_regulons.empty:
        reg_path = OUT_DIR / "regulons_table.csv"
        df_regulons.to_csv(str(reg_path), index=False)
        rlog(f"- Tabla de regulones guardada: {reg_path.name}")

    step_end("PASO 3 — Pruning por motivos (cisTarget)")

    # Si no hay regulones, no tiene sentido seguir con AUCell ni exportar
    if not regulons:
        raise RuntimeError(
            "No se identificaron regulones tras cisTarget. "
            "Verificar: (1) genoma de la DB coincide con la matriz, "
            "(2) versión de motifs.tbl compatible con pySCENIC, "
            "(3) gene symbols en la matriz coinciden con los de la DB."
        )

    # ==========================================================================
    # PASO 4: AUCell
    # ==========================================================================
    step_start("PASO 4 — Scoring de actividad (AUCell)")

    auc_matrix = aucell(
        expression_matrix=ex_matrix,
        signatures=regulons,
        num_workers=4,
    )

    rlog(f"- Matriz AUCell: {auc_matrix.shape[0]:,} células × {auc_matrix.shape[1]:,} regulones")
    rlog(f"- Score AUC: min={auc_matrix.values.min():.4f}, max={auc_matrix.values.max():.4f}, media={auc_matrix.values.mean():.4f}")

    # Guardar matriz AUCell
    auc_path = OUT_DIR / "auc_matrix.csv"
    auc_matrix.to_csv(str(auc_path))
    rlog(f"- Matriz guardada: {auc_path}")

    step_end("PASO 4 — Scoring de actividad (AUCell)")

    # ==========================================================================
    # PASO 5: VISUALIZACIÓN CON SCANPY
    # ==========================================================================
    step_start("PASO 5 — Visualización (scanpy)")

    # Construcción limpia del AnnData (sin warnings de obs_names asignados después)
    adata_scenic = sc.AnnData(
        X=auc_matrix.values,
        obs=pd.DataFrame(index=auc_matrix.index.astype(str)),
        var=pd.DataFrame(index=auc_matrix.columns.astype(str)),
    )

    sc.pp.neighbors(adata_scenic, random_state=RANDOM_SEED)
    sc.tl.umap(adata_scenic, random_state=RANDOM_SEED)
    sc.tl.leiden(
        adata_scenic,
        flavor="igraph",
        n_iterations=2,
        random_state=RANDOM_SEED,
    )

    n_clusters = adata_scenic.obs["leiden"].nunique()
    rlog(f"- Clusters Leiden sobre regulones: {n_clusters}")

    sc.pl.umap(
        adata_scenic, color="leiden",
        title="Clusters por actividad de regulones",
        show=False, save="_clusters.png",
    )

    if auc_matrix.shape[1] >= 4:
        # Top regulones por actividad media (más informativo que los primeros 4)
        top_regulon_names = auc_matrix.mean().sort_values(ascending=False).head(4).index.tolist()
        sc.pl.umap(
            adata_scenic,
            color=top_regulon_names,
            ncols=2,
            title=[f"Regulón: {r}" for r in top_regulon_names],
            show=False,
            save="_top_regulons.png",
        )
        rlog(f"- UMAP de top 4 regulones más activos guardado: {top_regulon_names}")

    step_end("PASO 5 — Visualización (scanpy)")

    # ==========================================================================
    # PASO 6: EXPORTAR A LOOM
    # ==========================================================================
    step_start("PASO 6 — Exportar a loom (SCope)")

    OUTPUT_LOOM = str(OUT_DIR / "scenic_output.loom")

    export2loom(
        ex_mtx=ex_matrix,
        regulons=regulons,
        cell_annotations=adata_scenic.obs["leiden"].to_dict(),
        out_fname=OUTPUT_LOOM,
        auc_mtx=auc_matrix,
        embeddings={
            "UMAP": pd.DataFrame(
                adata_scenic.obsm["X_umap"],
                index=adata_scenic.obs_names,
                columns=["UMAP_1", "UMAP_2"],
            )
        },
    )

    size_mb = Path(OUTPUT_LOOM).stat().st_size / (1024 * 1024)
    rlog(f"- Loom guardado: {OUTPUT_LOOM}")
    rlog(f"- Tamaño: {size_mb:.1f} MB")
    rlog("- Listo para visualizar en https://scope.aertslab.org")

    step_end("PASO 6 — Exportar a loom (SCope)")

except Exception as e:
    error = traceback.format_exc()
    print(f"\n*** ERROR ***\n{error}")

finally:
    write_report(
        ex_matrix=ex_matrix,
        tf_names=tf_names,
        modules=modules,
        regulons=regulons,
        auc_matrix=auc_matrix,
        adata_scenic=adata_scenic,
        error=error,
    )
    total = time.time() - _start_time
    print(f"\n{'='*60}")
    print(f"Carpeta de salida: {OUT_DIR.resolve()}")
    print(f"Duración total:    {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*60}")
