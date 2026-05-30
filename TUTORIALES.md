# Relación entre los dos tutoriales

Este proyecto une dos tutoriales oficiales que cubren **dos mitades complementarias**
del análisis de single-cell RNA-seq, ejecutados **tal cual el ecosistema de 2021**,
sobre **el mismo dataset (PBMC3k)** y en **un único entorno conda**.

| | Tutorial 1 | Tutorial 2 |
|---|---|---|
| **Fuente** | [scanpy — clustering PBMC3k](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html) | [aertslab/SCENICprotocol](https://github.com/aertslab/SCENICprotocol) |
| **Pregunta** | *¿Qué tipos de células hay?* | *¿Qué redes regulatorias las controlan?* |
| **Método** | QC → normalización → PCA → vecinos → UMAP → Leiden | `pyscenic grn` → `pyscenic ctx` → `pyscenic aucell` |
| **Salida** | Células agrupadas y anotadas (`.h5ad`) + loom de counts | Actividad de regulones por célula (matriz AUC) |
| **Script** | `01_scanpy_clustering.py` | `02_scenic_pipeline.py` |

---

## ¿Cómo se relacionan?

Son **secuenciales y complementarios**, no alternativos. scanpy te dice **qué**
células tienes; SCENIC te dice **por qué** se comportan así (qué programas
regulatorios las definen). El protocolo SCENIC está diseñado para correr **sobre
las mismas células** que ya procesó scanpy, y terminar integrando todo en un loom.

```
   PBMC3k (2,700 células)
          │
          ▼
   ┌──────────────────┐
   │  TUTORIAL 1       │   scanpy: QC, clustering, tipos celulares
   │  01_scanpy        │   (CD4 T, Monocytes, B, NK, ...)
   └──────────────────┘
          │  adata_clustered.h5ad
          │  pbmc3k_for_scenic.loom  ◀── counts crudos, mismos barcodes
          ▼
   ┌──────────────────┐
   │  TUTORIAL 2       │   SCENIC: GRNBoost2 → cisTarget → AUCell
   │  02_scenic        │   regulones reales (símbolos HGNC vs DB hg38)
   └──────────────────┘
          │  auc_matrix.csv  ◀── mismos barcodes que scanpy
          ▼
   ┌──────────────────┐
   │  INTEGRACIÓN      │   03: cruza clusters de scanpy con regulones
   │  03_integrate     │   "este tipo celular usa estos regulones"
   └──────────────────┘
          │  integrated_output.loom
          ▼
   Biología interpretable
```

El **puente** es el dataset compartido: `01` exporta `pbmc3k_for_scenic.loom` con
los counts crudos y los **mismos códigos de barras** que su `adata_clustered.h5ad`.
`02` corre SCENIC sobre ese loom, así que la matriz AUC sale indexada por los
mismos barcodes. Por eso `03` puede cruzarlos célula a célula (sin NaN).

---

## El problema técnico que resuelve este proyecto

El tutorial de scanpy del enlace es la versión *actual* (2024) y el protocolo
SCENIC es de *2021*. Sus versiones de software son mutuamente incompatibles:
SCENIC necesita `numpy <1.24` (antes de eliminar `np.object`) y `pandas <2.0`,
mientras que el scanpy moderno arrastra `numpy 1.26` + `pandas 2.x`.

**La solución adoptada aquí: retroceder todo al ecosistema 2021.** Se usa la
versión 2021 del tutorial de scanpy (el workflow clásico de **PBMC3k**) para que
ambos tutoriales convivan en **un solo entorno conda** (`environment.yml`):

- `python 3.10`, `numpy 1.23.5`, `pandas 1.5.3`
- `scanpy 1.9.3`, `anndata 0.8.0`
- `pyscenic 0.12.1`, `arboreto 0.1.6`, `ctxcore 0.2.0`, `dask 2023.3.2`

Con este entorno, `pyscenic grn` corre **GRNBoost2 real** (arboreto + dask), tal
como el protocolo — sin reemplazos ni parches.

---

## Por qué PBMC3k para ambos

El tutorial de scanpy de 2021 usaba PBMC3k; el de SCENIC usaba PBMC10k. Son
muestras distintas (barcodes que no coinciden), así que correr cada tutorial en
su propio dataset haría que la integración no encontrara células en común.

Usando **PBMC3k para toda la cadena** se logra lo mejor de ambos: es el dataset
fiel del tutorial de scanpy, corre rápido en SCENIC, y como es **un único
dataset** la integración tiene biología real. Además, al ser PBMC humano con
símbolos de gen HGNC, `cisTarget` sí encuentra regulones enriquecidos contra la
base de datos hg38 (cosa imposible con datos sintéticos).

---

## Cómo reproducirlo en otra PC

```bash
conda env create -f environment.yml
conda activate scenic-2021

python 01_scanpy_clustering.py    # PBMC3k → clusters + pbmc3k_for_scenic.loom
python 02_scenic_pipeline.py      # GRNBoost2 real → cisTarget → AUCell
python 03_integrate_anndata.py    # une ambos → integrated_output.loom
```

Cada script deja una carpeta con timestamp en `runs/` con un `reporte_tecnico.md`
(métricas, criterios de validación y archivos generados).

### Notas de ejecución

- **`pyscenic grn` es el paso lento.** GRNBoost2 sobre ~2,700 células × ~13,000
  genes con ~1,800 TFs puede tardar de minutos a un par de horas según el número
  de núcleos. Se controla con la variable de entorno `SCENIC_WORKERS`
  (por defecto `min(4, núcleos)`).
- La primera corrida de `02` descarga ~390 MB de recursos cisTarget (lista de TFs,
  `motifs.tbl` y el ranking `feather` de hg38) a `scenic_data/`. Se cachean.
