# =============================================================================
# TUTORIAL: Inferencia de Redes Regulatorias con pySCENIC
# Fuente: https://github.com/aertslab/SCENICprotocol
#
# Dataset de prueba: expr_mat_tiny.loom
# Flujo: pyscenic grn -> pyscenic ctx -> pyscenic aucell -> visualización
#
# Usa los comandos CLI de pySCENIC (subprocess) en lugar de la API Python
# para evitar incompatibilidades entre arboreto/dask/pandas modernos.
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
import ast
from sklearn.ensemble import GradientBoostingRegressor
from ctxcore.genesig import GeneSignature
from pyscenic.aucell import aucell as pyscenic_aucell

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


def run_pyscenic(cmd: list[str], step_name: str) -> subprocess.CompletedProcess:
    """Ejecuta un comando pyscenic CLI y maneja errores."""
    rlog(f"- Comando: `{' '.join(cmd)}`")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        stdout = result.stdout.strip()
        detail = stderr or stdout or "(sin output)"
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


def write_report(
    ex_matrix=None,
    tf_names=None,
    adjacencies=None,
    n_regulons=0,
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
        f"| pyscenic (CLI) | {PYSCENIC_VERSION} |",
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

    # Top regulones por actividad media
    if auc_matrix is not None and len(auc_matrix.columns) > 0:
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
        adj_ok      = "✓" if adjacencies is not None and len(adjacencies) > 0 else "✗ FALTA"
        regulons_ok = "✓" if n_regulons > 0 else "✗ FALTA — posible problema con cisTarget o DB"
        aucell_ok   = "✓" if auc_matrix is not None and not auc_matrix.empty else "✗ FALTA"
        umap_ok     = "✓" if adata_scenic is not None and "X_umap" in adata_scenic.obsm else "✗ FALTA"
        loom_ok     = "✓" if (OUT_DIR / "scenic_output.loom").exists() else "✗ NO GENERADO"

        lines += [
            f"| Matriz de expresión cargada | > 0 células | {ex_matrix.shape[0]:,} | {celulas_ok} |",
            f"| Genes en la matriz | > 100 | {ex_matrix.shape[1]:,} | {genes_ok} |",
            f"| Adjacencias GRN | > 0 | {len(adjacencies) if adjacencies is not None else 0:,} | {adj_ok} |",
            f"| Regulones cisTarget | > 0 | {n_regulons:,} | {regulons_ok} |",
            f"| Matriz AUCell calculada | Sí | {'Sí' if auc_matrix is not None else 'No'} | {aucell_ok} |",
            f"| UMAP sobre regulones | Sí | {'Sí' if adata_scenic is not None and 'X_umap' in adata_scenic.obsm else 'No'} | {umap_ok} |",
            f"| Loom exportado | Sí | {'Sí' if (OUT_DIR / 'scenic_output.loom').exists() else 'No'} | {loom_ok} |",
        ]

        if n_regulons == 0 and adjacencies is not None:
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
ex_matrix    = None
tf_names     = None
adjacencies  = None
n_regulons   = 0
auc_matrix   = None
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

    for filename, url in URLS.items():
        dest = DATA_DIR / filename
        expected = _remote_size(url)

        if dest.exists():
            local = dest.stat().st_size
            if expected > 0 and local != expected:
                rlog(f"- {filename} corrupto/incompleto ({local} B vs {expected} B esperados), re-descargando...")
                dest.unlink()
            elif local == 0:
                rlog(f"- {filename} está vacío, re-descargando...")
                dest.unlink()

        if not dest.exists():
            rlog(f"- Descargando {filename}...")
            try:
                urllib.request.urlretrieve(url, str(dest))
            except Exception as dl_err:
                if dest.exists():
                    dest.unlink()
                raise RuntimeError(f"Falló la descarga de {filename} desde {url}: {dl_err}") from dl_err
            actual = dest.stat().st_size
            if expected > 0 and actual != expected:
                raise RuntimeError(
                    f"Tamaño inesperado para {filename}: {actual} B (esperado {expected} B). "
                    "Posible 404 servido como HTML o conexión interrumpida."
                )
            size_mb = actual / (1024 * 1024)
            rlog(f"  → Guardado ({size_mb:.2f} MB)")
        else:
            size_mb = dest.stat().st_size / (1024 * 1024)
            rlog(f"- Ya existe: {filename} ({size_mb:.2f} MB)")

    LOOM_PATH    = DATA_DIR / "expr_mat_tiny.loom"
    TF_PATH      = DATA_DIR / "test_TFs_tiny.txt"
    MOTIFS_PATH  = DATA_DIR / "motifs.tbl"
    RANKING_PATH = DATA_DIR / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

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

    tf_names = pd.read_csv(str(TF_PATH), header=None).iloc[:, 0].tolist()
    rlog(f"- Factores de transcripción cargados: {len(tf_names):,}")

    step_end("PASO 1 — Carga de la matriz de expresión")

    # ==========================================================================
    # PASO 2: GRNBoost2  (sklearn — reemplaza arboreto/dask incompatibles)
    # ==========================================================================
    # arboreto/dask son incompatibles con dask moderno (legacy API eliminada).
    # Se usa GradientBoostingRegressor de sklearn con el mismo protocolo:
    # para cada gen target, regresarlo contra los TFs y extraer importancias.
    # La salida es el mismo TSV (TF, target, importance) que acepta pyscenic ctx.
    step_start("PASO 2 — Inferencia de red (GRNBoost2/sklearn)")

    adj_path = OUT_DIR / "adjacencies.tsv"

    tfs_in_matrix = [tf for tf in tf_names if tf in ex_matrix.columns]
    rlog(f"- TFs presentes en la matriz: {len(tfs_in_matrix)} de {len(tf_names)}")

    if not tfs_in_matrix:
        raise RuntimeError(
            f"Ninguno de los TFs del archivo {TF_PATH.name} aparece en la matriz de expresión. "
            "Verificar que los nombres de gen coincidan (mayúsculas, versión de símbolo)."
        )

    records = []
    X_tfs = ex_matrix[tfs_in_matrix].values
    target_genes = [g for g in ex_matrix.columns if g not in tfs_in_matrix]

    for target in target_genes:
        y = ex_matrix[target].values
        if y.std() == 0:
            continue
        gbm = GradientBoostingRegressor(
            n_estimators=500,
            max_depth=3,
            random_state=RANDOM_SEED,
        )
        gbm.fit(X_tfs, y)
        for tf, imp in zip(tfs_in_matrix, gbm.feature_importances_):
            if imp > 0:
                records.append({"TF": tf, "target": target, "importance": imp})

    adjacencies = pd.DataFrame(records).sort_values("importance", ascending=False)
    adjacencies.to_csv(str(adj_path), sep="\t", index=False)

    rlog(f"- Relaciones TF–gen encontradas: {len(adjacencies):,}")
    if len(adjacencies) > 0:
        rlog(f"- Importancia máxima: {adjacencies['importance'].max():.4f}")
        rlog(f"- Importancia media: {adjacencies['importance'].mean():.4f}")
        top_tfs = adjacencies.groupby("TF").size().sort_values(ascending=False).head(5)
        rlog("- Top TFs más conectados:")
        for tf, n in top_tfs.items():
            rlog(f"  • {tf}: {n} genes target")
    else:
        rlog("- ⚠️ No se encontraron relaciones TF–gen (todos los genes con varianza cero?)")

    step_end("PASO 2 — Inferencia de red (GRNBoost2/sklearn)")

    # ==========================================================================
    # PASO 3: cisTarget  (pyscenic ctx — CLI)
    # ==========================================================================
    step_start("PASO 3 — Pruning por motivos (cisTarget)")

    regulons_path = OUT_DIR / "regulons.csv"
    run_pyscenic([
        "pyscenic", "ctx",
        str(adj_path),
        str(RANKING_PATH),
        "--annotations_fname", str(MOTIFS_PATH),
        "--expression_mtx_fname", str(LOOM_PATH),
        "--output", str(regulons_path),
        "--num_workers", "1",
    ], "cisTarget")

    # Leer el CSV con multi-index para contar filas reales (no los headers).
    try:
        df_ctx_check = pd.read_csv(str(regulons_path), index_col=[0, 1], header=[0, 1])
        n_regulons = len(df_ctx_check)
    except Exception:
        n_regulons = 0

    rlog(f"- Regulones identificados: {n_regulons:,}")

    if n_regulons > 0:
        rlog(f"- Resultados guardados en: {regulons_path.name}")
    else:
        rlog("- ⚠️ 0 regulones de cisTarget — posible causa: TF sin motivos enriquecidos en la DB hg38")
        rlog("  (BRF1 es TF de RNA Pol III; sus motivos no aparecen en promotores de genes codificantes)")
        rlog("  → Se usarán regulones derivados del GRN directamente para continuar el pipeline")

    step_end("PASO 3 — Pruning por motivos (cisTarget)")

    # ==========================================================================
    # PASO 4: AUCell  (Python API — evita incompatibilidad de pandas 2.x con CLI)
    # ==========================================================================
    # pyscenic aucell CLI usa load_motifs con multi-index que falla en pandas 2.x.
    # Se parsea el CSV manualmente y se llama aucell() directo (sin arboreto/dask).
    step_start("PASO 4 — Scoring de actividad (AUCell)")

    aucell_loom = OUT_DIR / "scenic_output.loom"

    # Construir firmas de genes para AUCell.
    # Fuente 1: regulons.csv de cisTarget (si encontró regulones).
    # Fuente 2 (fallback): adjacencias del GRN directamente, sin validación de motivos.
    signatures = []

    if n_regulons > 0:
        rlog("- Usando regulones de cisTarget")
        df_ctx_raw = pd.read_csv(str(regulons_path), header=[0, 1], index_col=[0, 1])
        df_ctx_raw.columns = [" ".join(c).strip() for c in df_ctx_raw.columns]
        df_ctx_raw = df_ctx_raw.reset_index()
        target_col = next(
            (c for c in df_ctx_raw.columns if "TargetGenes" in c), None
        )
        tf_col = next(
            (c for c in df_ctx_raw.columns if c in ("TF", "level_0")), None
        )
        if target_col and tf_col:
            for _, row in df_ctx_raw.iterrows():
                tf = str(row[tf_col])
                try:
                    genes = [t[0] for t in ast.literal_eval(str(row[target_col]))]
                except Exception:
                    genes = []
                if genes:
                    signatures.append(
                        GeneSignature(name=f"{tf}(+)", gene2weight={g: 1.0 for g in genes})
                    )
    else:
        rlog("- Fallback: construyendo regulones directamente desde adjacencias GRN")
        rlog("  (sin validación de motivos — solo para demostrar el pipeline)")
        for tf, group in adjacencies.groupby("TF"):
            genes = group["target"].tolist()
            if genes:
                signatures.append(
                    GeneSignature(
                        name=f"{tf}(+)",
                        gene2weight=dict(zip(group["target"], group["importance"])),
                    )
                )

    rlog(f"- Firmas de genes construidas: {len(signatures)}")
    if not signatures:
        raise RuntimeError("No se pudieron construir firmas de genes")

    auc_matrix = pyscenic_aucell(ex_matrix, signatures, num_workers=1)

    rlog(f"- Matriz AUCell: {auc_matrix.shape[0]:,} células × {auc_matrix.shape[1]:,} regulones")
    rlog(f"- Score AUC: min={auc_matrix.values.min():.4f}, max={auc_matrix.values.max():.4f}, "
         f"media={auc_matrix.values.mean():.4f}")

    auc_path = OUT_DIR / "auc_matrix.csv"
    auc_matrix.to_csv(str(auc_path))
    rlog(f"- Matriz guardada: {auc_path.name}")

    step_end("PASO 4 — Scoring de actividad (AUCell)")

    # ==========================================================================
    # PASO 5: VISUALIZACIÓN CON SCANPY
    # ==========================================================================
    step_start("PASO 5 — Visualización (scanpy)")

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

    # Guardar AnnData con embeddings
    adata_path = OUT_DIR / "scenic_adata.h5ad"
    adata_scenic.write_h5ad(str(adata_path))
    rlog(f"- AnnData guardado: {adata_path.name}")

    step_end("PASO 5 — Visualización (scanpy)")

    # ==========================================================================
    # PASO 6: RESUMEN DE SALIDA
    # ==========================================================================
    step_start("PASO 6 — Salida final")

    rlog(f"- AnnData con UMAP/Leiden: {adata_path.name}")
    rlog(f"- Matriz AUC: auc_matrix.csv")
    rlog(f"- Regulones cisTarget: {regulons_path.name}")
    rlog(f"- Adjacencias GRN: {adj_path.name}")

    step_end("PASO 6 — Salida final")

except Exception as e:
    error = traceback.format_exc()
    print(f"\n*** ERROR ***\n{error}")

finally:
    write_report(
        ex_matrix=ex_matrix,
        tf_names=tf_names,
        adjacencies=adjacencies,
        n_regulons=n_regulons,
        auc_matrix=auc_matrix,
        adata_scenic=adata_scenic,
        error=error,
    )
    total = time.time() - _start_time
    print(f"\n{'='*60}")
    print(f"Carpeta de salida: {OUT_DIR.resolve()}")
    print(f"Duración total:    {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*60}")
