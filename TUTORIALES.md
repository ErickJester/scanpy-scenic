# Relación entre los dos tutoriales

Este proyecto une dos tutoriales oficiales que abordan **dos mitades complementarias**
del análisis de single-cell RNA-seq.

| | Tutorial 1 | Tutorial 2 |
|---|---|---|
| **Fuente** | [scanpy — clustering](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html) | [aertslab/SCENICprotocol](https://github.com/aertslab/SCENICprotocol) |
| **Pregunta que responde** | *¿Qué tipos de células hay?* | *¿Qué redes regulatorias las controlan?* |
| **Entrada** | Matriz cruda de conteos (10X `.h5`) | Matriz de expresión (`.loom`) |
| **Método** | QC → normalización → PCA → vecinos → UMAP → Leiden | GRNBoost2 → cisTarget → AUCell |
| **Salida** | Células agrupadas y anotadas (`.h5ad`) | Actividad de regulones por célula (matriz AUC) |
| **Script** | `01_scanpy_clustering.py` | `02_scenic_pipeline.py` |

---

## ¿Cómo se relacionan?

Los dos tutoriales son **secuenciales y complementarios**, no alternativos:

```
   Datos crudos (10X)
          │
          ▼
   ┌──────────────────┐
   │  TUTORIAL 1       │   scanpy: agrupa células en tipos celulares
   │  scanpy clustering│   (Monocytes, Lymphocytes, B Cells, ...)
   └──────────────────┘
          │  adata_clustered.h5ad
          ▼
   ┌──────────────────┐
   │  TUTORIAL 2       │   SCENIC: infiere qué factores de transcripción
   │  SCENIC pipeline  │   están activos en cada célula
   └──────────────────┘
          │  auc_matrix.csv
          ▼
   ┌──────────────────┐
   │  INTEGRACIÓN      │   03: combina ambos → "este tipo celular usa
   │  (este proyecto)  │   estos regulones"
   └──────────────────┘
          │  integrated_output.loom
          ▼
   Biología interpretable
```

scanpy te dice **qué** células tienes. SCENIC te dice **por qué** se comportan así
(qué programas regulatorios las definen). El valor está en cruzar ambos: descubrir
que, por ejemplo, los monocitos están dominados por el regulón `SPI1(+)`.

---

## El problema técnico clave: viven en eras de software distintas

Esta es la parte importante de la relación, y la razón por la que el proyecto
usa **dos entornos conda separados**.

- El tutorial de **scanpy** es la versión *stable actual* (2024+). Usa APIs modernas
  (`sc.pp.scrublet`, `leiden flavor="igraph"`) que **exigen scanpy ≥1.10**, y por
  arrastre `numpy 1.26` + `pandas 2.x`.

- El protocolo **SCENIC** se publicó en *2021*. Su motor de inferencia de redes
  (`GRNBoost2` vía `arboreto` + `dask`) y `pyscenic` se escribieron contra
  `numpy <1.24` (antes de eliminar `np.object`, `np.float`) y `pandas <2.0`
  (antes del cambio en la API de `MultiIndex`).

Estos dos rangos de versiones **son mutuamente incompatibles**: no existe un solo
entorno donde ambos tutoriales corran "tal cual". Si fuerzas a SCENIC al entorno
moderno, fallan `arboreto`/`pyscenic`; si fuerzas a scanpy al entorno viejo, falla
el tutorial de clustering.

### La solución: traspaso por archivos

Los dos pipelines **no comparten una sesión de Python**: se comunican escribiendo y
leyendo archivos en disco.

```
[entorno scanpy-tutorial]  01 ──escribe──▶  adata_clustered.h5ad
                                                   │
[entorno scenic-tutorial]  02 ──escribe──▶  auc_matrix.csv
                                                   │
[cualquier entorno]        03 ──lee ambos──▶ integrated_output.loom
```

Como el handoff es por archivos (`.h5ad`, `.csv`, `.loom`), **cada pipeline puede
correr en el entorno conda que necesita**, en momentos distintos e incluso en
máquinas distintas. Esto es exactamente lo que recomienda el protocolo SCENIC
oficial, que termina exportando un `.loom` integrado para SCope.

---

## Cómo reproducirlo en otra PC

```bash
# --- Tutorial 1: scanpy ---
conda env create -f environment-scanpy.yml
conda activate scanpy-tutorial
python 01_scanpy_clustering.py        # → runs/scanpy_*/adata_clustered.h5ad

# --- Tutorial 2: SCENIC ---
conda env create -f environment-scenic.yml
conda activate scenic-tutorial
python 02_scenic_pipeline.py          # → runs/scenic_*/auc_matrix.csv

# --- Integración (cualquiera de los dos entornos sirve) ---
python 03_integrate_anndata.py        # → runs/integrate_*/integrated_output.loom
```

Cada script deja una carpeta con timestamp en `runs/` que incluye un
`reporte_tecnico.md` con métricas, criterios de validación y los archivos generados.

---

## Nota sobre la implementación actual de GRNBoost2

En el entorno moderno, `02_scenic_pipeline.py` reemplaza GRNBoost2 por una
implementación equivalente con `sklearn.GradientBoostingRegressor` (mismo protocolo:
regresar cada gen contra los TFs y extraer importancias), porque `arboreto` no corre
con dask moderno. **Con el entorno `scenic-tutorial` (versiones fijadas), se puede
usar el GRNBoost2 real del tutorial**, ya que `arboreto` y `dask 2023.3.2` son
compatibles.
