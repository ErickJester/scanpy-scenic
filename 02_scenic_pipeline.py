# =============================================================================
# TUTORIAL 2: Inferencia de Redes Regulatorias con pySCENIC — protocolo real
# Fuente: https://github.com/aertslab/SCENICprotocol
#
# Flujo OFICIAL tal cual, con los 3 comandos CLI de pySCENIC:
#     pyscenic grn     -> GRNBoost2 real (arboreto + dask)
#     pyscenic ctx     -> cisTarget (pruning por motivos)
#     pyscenic aucell  -> scoring de actividad de regulones
#
# Entrada: el loom de counts crudos que genera 01_scanpy_clustering.py
#          (médula ósea, las MISMAS células) -> runs/scanpy_*/bonemarrow_for_scenic.loom
# Salida:  auc_matrix.csv + pyscenic_output.loom
#
# Requiere el entorno conda `scenic-medula` (ver environment.yml). Con ese
# entorno pySCENIC corre nativo: no hace falta parchear numpy/pandas.
# =============================================================================

import os
import sys
import time
import platform
import traceback
import subprocess
import urllib.request
from datetime import datetime
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import pandas as pd
import numpy as np
import scanpy as sc
import loompy

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Número de workers para GRNBoost2/cisTarget (configurable por entorno).
N_WORKERS = int(os.environ.get("SCENIC_WORKERS", str(min(4, os.cpu_count() or 1))))

# =============================================================================
# INFRAESTRUCTURA DE SALIDA
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


def run_pyscenic(cmd: list[str], step_name: str) -> subprocess.CompletedProcess:
    """Ejecuta un comando pyscenic CLI y maneja errores."""
    rlog(f"- Comando: `{' '.join(cmd)}`")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=10800)
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "(sin output)"
        raise RuntimeError(f"{step_name} falló (exit code {result.returncode}):\n{detail}")
    if result.stderr.strip():
        for line in result.stderr.strip().split("\n")[-5:]:
            rlog(f"  [log] {line}")
    return result


def _get_pyscenic_version() -> str:
    try:
        r = subprocess.run(["pyscenic", "--version"], capture_output=True, text=True, timeout=10)
        v = (r.stdout.strip() or r.stderr.strip()).replace("pyscenic", "").strip()
        return v or "unknown"
    except Exception:
        return "unknown"


PYSCENIC_VERSION = _get_pyscenic_version()


def _find_latest_scanpy_loom() -> Path | None:
    """Busca el loom de SCENIC en la corrida runs/scanpy_* más reciente."""
    for run in sorted(Path("runs").glob("scanpy_*"), reverse=True):
        loom = run / "bonemarrow_for_scenic.loom"
        if loom.exists():
            return loom
    return None


def write_report(
    ex_matrix=None, tf_names=None, adjacencies=None,
    n_regulons=0, auc_matrix=None, adata_scenic=None, error: str = None,
) -> None:
    total = time.time() - _start_time
    lines = []

    lines += [
        "# Reporte Técnico — SCENIC pipeline (médula ósea, protocolo real)",
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
        f"| pyscenic (CLI) | {PYSCENIC_VERSION} |",
        f"| pandas | {pd.__version__} |",
        f"| numpy | {np.__version__} |",
        f"| workers GRNBoost2 | {N_WORKERS} |",
        "",
        "---",
        "",
        "## Log de ejecución por pasos",
        "",
    ]

    lines += _report_lines

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
        if adjacencies is not None:
            n_tfs = adjacencies["TF"].nunique() if "TF" in adjacencies.columns else 0
            lines.append(f"| GRNBoost2 | Relaciones TF-gen | {len(adjacencies):,} |")
            lines.append(f"| GRNBoost2 | TFs con targets | {n_tfs:,} |")
        if n_regulons > 0:
            lines.append(f"| cisTarget | Regulones finales | {n_regulons:,} |")
        if auc_matrix is not None:
            lines.append(f"| AUCell | Células × Regulones | {auc_matrix.shape[0]:,} × {auc_matrix.shape[1]:,} |")
        if adata_scenic is not None and "leiden" in adata_scenic.obs.columns:
            lines.append(f"| Clustering | Clusters Leiden | {adata_scenic.obs['leiden'].nunique()} |")

    if auc_matrix is not None and len(auc_matrix.columns) > 0:
        top_regulons = auc_matrix.mean().sort_values(ascending=False).head(10)
        lines += ["", "**Top 10 regulones por actividad media:**", "", "| Regulón | AUC medio |", "|---|---|"]
        for reg, val in top_regulons.items():
            lines.append(f"| {reg} | {val:.4f} |")

    lines += [
        "", "---", "", "## Criterios de validación (para revisión IA)", "",
        "| Criterio | Valor esperado | Valor obtenido | Estado |",
        "|---|---|---|---|",
    ]
    if ex_matrix is not None:
        celulas_ok  = "✓" if ex_matrix.shape[0] > 0 else "✗ VACÍO"
        genes_ok    = "✓" if ex_matrix.shape[1] > 100 else "✗ MUY POCOS"
        adj_ok      = "✓" if adjacencies is not None and len(adjacencies) > 0 else "✗ FALTA"
        regulons_ok = "✓" if n_regulons > 0 else "✗ FALTA — revisar cisTarget/DB"
        aucell_ok   = "✓" if auc_matrix is not None and not auc_matrix.empty else "✗ FALTA"
        umap_ok     = "✓" if adata_scenic is not None and "X_umap" in adata_scenic.obsm else "✗ FALTA"
        loom_ok     = "✓" if (OUT_DIR / "pyscenic_output.loom").exists() else "✗ NO GENERADO"

        lines += [
            f"| Matriz de expresión cargada | > 0 células | {ex_matrix.shape[0]:,} | {celulas_ok} |",
            f"| Genes en la matriz | > 100 | {ex_matrix.shape[1]:,} | {genes_ok} |",
            f"| Adjacencias GRN | > 0 | {len(adjacencies) if adjacencies is not None else 0:,} | {adj_ok} |",
            f"| Regulones cisTarget | > 0 | {n_regulons:,} | {regulons_ok} |",
            f"| Matriz AUCell calculada | Sí | {'Sí' if auc_matrix is not None else 'No'} | {aucell_ok} |",
            f"| UMAP sobre regulones | Sí | {'Sí' if adata_scenic is not None and 'X_umap' in adata_scenic.obsm else 'No'} | {umap_ok} |",
            f"| Loom pySCENIC exportado | Sí | {'Sí' if (OUT_DIR / 'pyscenic_output.loom').exists() else 'No'} | {loom_ok} |",
        ]
    else:
        lines.append("| Pipeline SCENIC | Completado | NO COMPLETADO | ✗ REVISAR |")

    if error:
        lines += [
            "", "---", "", "## Error capturado", "", "```", error, "```", "",
            "> ⚠️ **Para IA:** Analizar el traceback e identificar en qué paso falló.",
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

ex_matrix    = None
tf_names     = None
adjacencies  = None
n_regulons   = 0
auc_matrix   = None
adata_scenic = None

try:
    # ==========================================================================
    # PASO 0: LOCALIZAR LOOM DE 01 + DESCARGAR RECURSOS cisTarget
    # ==========================================================================
    step_start("PASO 0 — Entrada y recursos")

    LOOM_PATH = _find_latest_scanpy_loom()
    if LOOM_PATH is None:
        raise FileNotFoundError(
            "No se encontró runs/scanpy_*/bonemarrow_for_scenic.loom. "
            "Ejecuta primero 01_scanpy_clustering.py."
        )
    rlog(f"- Loom de entrada: `{LOOM_PATH}`")

    DATA_DIR = Path("scenic_data")
    DATA_DIR.mkdir(exist_ok=True)

    REMOTE_FILES = {
        "allTFs_hg38.txt": (
            "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
        ),
        "motifs.tbl": (
            "https://resources.aertslab.org/cistarget/motif2tf/"
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        ),
        "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather": (
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/"
            "hg38/refseq_r80/mc_v10_clust/gene_based/"
            "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        ),
    }

    def _remote_size(url: str) -> int:
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=30) as resp:
                return int(resp.headers.get("Content-Length", -1))
        except Exception:
            return -1

    for filename, url in REMOTE_FILES.items():
        dest = DATA_DIR / filename
        expected = _remote_size(url)
        if dest.exists():
            local = dest.stat().st_size
            if (expected > 0 and local != expected) or local == 0:
                rlog(f"- {filename} incompleto, re-descargando...")
                dest.unlink()
        if not dest.exists():
            rlog(f"- Descargando {filename}...")
            try:
                urllib.request.urlretrieve(url, str(dest))
            except Exception as dl_err:
                if dest.exists():
                    dest.unlink()
                raise RuntimeError(f"Falló la descarga de {filename}: {dl_err}") from dl_err
            actual = dest.stat().st_size
            if expected > 0 and actual != expected:
                raise RuntimeError(f"Tamaño inesperado para {filename}: {actual} B (esperado {expected} B).")
            rlog(f"  → Guardado ({actual/(1024*1024):.2f} MB)")
        else:
            rlog(f"- Ya existe: {filename} ({dest.stat().st_size/(1024*1024):.2f} MB)")

    TF_PATH      = DATA_DIR / "allTFs_hg38.txt"
    MOTIFS_PATH  = DATA_DIR / "motifs.tbl"
    RANKING_PATH = DATA_DIR / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

    step_end("PASO 0 — Entrada y recursos")

    # ==========================================================================
    # PASO 1: CARGAR MATRIZ (para métricas)
    # ==========================================================================
    step_start("PASO 1 — Carga de la matriz de expresión")

    with loompy.connect(str(LOOM_PATH)) as ds:
        ex_matrix = pd.DataFrame(
            data=ds[:, :].T, index=ds.ca["CellID"], columns=ds.ra["Gene"]
        )
    rlog(f"- Matriz: {ex_matrix.shape[0]:,} células × {ex_matrix.shape[1]:,} genes")
    rlog(f"- % de ceros: {(ex_matrix.values == 0).mean() * 100:.1f}%")

    tf_names = pd.read_csv(str(TF_PATH), header=None).iloc[:, 0].tolist()
    rlog(f"- Factores de transcripción (hg38): {len(tf_names):,}")

    step_end("PASO 1 — Carga de la matriz de expresión")

    # ==========================================================================
    # PASO 2: GRNBoost2 — pyscenic grn (REAL, arboreto + dask)
    # ==========================================================================
    step_start("PASO 2 — Inferencia de red (pyscenic grn / GRNBoost2)")

    adj_path = OUT_DIR / "adjacencies.tsv"
    run_pyscenic([
        "pyscenic", "grn",
        str(LOOM_PATH),
        str(TF_PATH),
        "--output", str(adj_path),
        "--num_workers", str(N_WORKERS),
        "--seed", str(RANDOM_SEED),
        "--method", "grnboost2",
    ], "GRNBoost2")

    adjacencies = pd.read_csv(str(adj_path), sep="\t")
    rlog(f"- Relaciones TF–gen: {len(adjacencies):,}")
    if len(adjacencies) > 0:
        rlog(f"- Importancia máxima: {adjacencies['importance'].max():.4f}")
        top_tfs = adjacencies.groupby("TF").size().sort_values(ascending=False).head(5)
        rlog("- Top TFs más conectados:")
        for tf, n in top_tfs.items():
            rlog(f"  • {tf}: {n} genes target")

    step_end("PASO 2 — Inferencia de red (pyscenic grn / GRNBoost2)")

    # ==========================================================================
    # PASO 3: cisTarget — pyscenic ctx
    # ==========================================================================
    step_start("PASO 3 — Pruning por motivos (pyscenic ctx / cisTarget)")

    regulons_path = OUT_DIR / "regulons.csv"
    run_pyscenic([
        "pyscenic", "ctx",
        str(adj_path),
        str(RANKING_PATH),
        "--annotations_fname", str(MOTIFS_PATH),
        "--expression_mtx_fname", str(LOOM_PATH),
        "--output", str(regulons_path),
        "--num_workers", str(N_WORKERS),
    ], "cisTarget")

    try:
        df_ctx_check = pd.read_csv(str(regulons_path), index_col=[0, 1], header=[0, 1])
        n_regulons = len(df_ctx_check)
    except Exception:
        n_regulons = 0
    rlog(f"- Regulones identificados: {n_regulons:,}")

    step_end("PASO 3 — Pruning por motivos (pyscenic ctx / cisTarget)")

    # ==========================================================================
    # PASO 4: AUCell — pyscenic aucell
    # ==========================================================================
    step_start("PASO 4 — Scoring de actividad (pyscenic aucell)")

    pyscenic_loom = OUT_DIR / "pyscenic_output.loom"
    run_pyscenic([
        "pyscenic", "aucell",
        str(LOOM_PATH),
        str(regulons_path),
        "--output", str(pyscenic_loom),
        "--num_workers", str(N_WORKERS),
        "--seed", str(RANDOM_SEED),
    ], "AUCell")

    # Extraer la matriz AUC del loom de salida (col attr RegulonsAUC)
    with loompy.connect(str(pyscenic_loom), mode="r") as ds:
        auc_data = ds.ca["RegulonsAUC"]
        cell_ids = ds.ca["CellID"]
    auc_matrix = pd.DataFrame(auc_data, index=cell_ids)
    auc_matrix.index.name = "Cell"

    rlog(f"- Matriz AUCell: {auc_matrix.shape[0]:,} células × {auc_matrix.shape[1]:,} regulones")
    rlog(f"- Score AUC: min={auc_matrix.values.min():.4f}, max={auc_matrix.values.max():.4f}, "
         f"media={auc_matrix.values.mean():.4f}")

    auc_matrix.to_csv(str(OUT_DIR / "auc_matrix.csv"))
    rlog("- Matriz guardada: auc_matrix.csv")

    step_end("PASO 4 — Scoring de actividad (pyscenic aucell)")

    # ==========================================================================
    # PASO 5: VISUALIZACIÓN
    # ==========================================================================
    step_start("PASO 5 — Visualización (scanpy)")

    adata_scenic = sc.AnnData(
        X=auc_matrix.values,
        obs=pd.DataFrame(index=auc_matrix.index.astype(str)),
        var=pd.DataFrame(index=auc_matrix.columns.astype(str)),
    )
    sc.pp.neighbors(adata_scenic, random_state=RANDOM_SEED)
    sc.tl.umap(adata_scenic, random_state=RANDOM_SEED)
    sc.tl.leiden(adata_scenic, random_state=RANDOM_SEED)
    rlog(f"- Clusters Leiden sobre regulones: {adata_scenic.obs['leiden'].nunique()}")

    sc.pl.umap(adata_scenic, color="leiden",
               title="Clusters por actividad de regulones",
               show=False, save="_clusters.png")

    if auc_matrix.shape[1] >= 4:
        top_regs = auc_matrix.mean().sort_values(ascending=False).head(4).index.tolist()
        sc.pl.umap(adata_scenic, color=top_regs, ncols=2,
                   title=[f"Regulón: {r}" for r in top_regs],
                   show=False, save="_top_regulons.png")
        rlog(f"- UMAP de top 4 regulones: {top_regs}")

    adata_scenic.write_h5ad(str(OUT_DIR / "scenic_adata.h5ad"))
    rlog("- AnnData guardado: scenic_adata.h5ad")

    step_end("PASO 5 — Visualización (scanpy)")

except Exception as e:
    error = traceback.format_exc()
    print(f"\n*** ERROR ***\n{error}")

finally:
    write_report(
        ex_matrix=ex_matrix, tf_names=tf_names, adjacencies=adjacencies,
        n_regulons=n_regulons, auc_matrix=auc_matrix, adata_scenic=adata_scenic,
        error=error,
    )
    total = time.time() - _start_time
    print(f"\n{'='*60}")
    print(f"Carpeta de salida: {OUT_DIR.resolve()}")
    print(f"Duración total:    {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*60}")
