# Relación entre los dos tutoriales

Este proyecto une dos tutoriales oficiales que cubren **dos mitades complementarias**
del análisis de single-cell RNA-seq, ejecutados sobre **el mismo dataset de médula
ósea humana** (el del enlace actual de scanpy) y en **un único entorno conda**.

| | Tutorial 1 | Tutorial 2 |
|---|---|---|
| **Fuente** | [scanpy — clustering](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html) | [aertslab/SCENICprotocol](https://github.com/aertslab/SCENICprotocol) |
| **Pregunta** | *¿Qué tipos de células hay?* | *¿Qué redes regulatorias las controlan?* |
| **Dataset** | Médula ósea humana (~17,000 células) | El mismo (vía loom exportado por 01) |
| **Método** | QC → normalización → PCA → vecinos → UMAP → Leiden | `pyscenic grn` → `pyscenic ctx` → `pyscenic aucell` |
| **Salida** | Células agrupadas y anotadas (`.h5ad`) + loom de counts | Actividad de regulones por célula (matriz AUC) |
| **Script** | `01_scanpy_clustering.py` | `02_scenic_pipeline.py` |

---

## ¿Cómo se relacionan?

Son **secuenciales y complementarios**, no alternativos. scanpy te dice **qué**
células tienes; SCENIC te dice **por qué** se comportan así (qué programas
regulatorios las definen). El protocolo SCENIC está diseñado para correr **sobre
las mismas células** que ya procesó scanpy y terminar integrando todo en un loom.

```
   Médula ósea (~17,000 células)
          │
          ▼
   ┌──────────────────┐
   │  TUTORIAL 1       │   scanpy: QC, dobletes, clustering, tipos celulares
   │  01_scanpy        │   (Monocytes, Lymphocytes, B Cells, Erythroid, ...)
   └──────────────────┘
          │  adata_clustered.h5ad
          │  bonemarrow_for_scenic.loom  ◀── counts crudos, mismos barcodes
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

El **puente** es el dataset compartido: `01` exporta `bonemarrow_for_scenic.loom`
con los counts crudos y los **mismos códigos de barras** que su
`adata_clustered.h5ad`. `02` corre SCENIC sobre ese loom, así que la matriz AUC
sale indexada por los mismos barcodes. Por eso `03` los cruza célula a célula
(sin NaN).

---

## El problema técnico que resuelve este proyecto

El tutorial de scanpy del enlace es la versión *actual* (2024) y el protocolo
SCENIC es de *2021*. Sus versiones de software chocan: SCENIC necesita
`numpy <1.24` (para que `pyscenic 0.12.1` no rompa con `np.object`) y `pandas <2.0`,
mientras que el scanpy más nuevo arrastra `numpy 1.26` + `pandas 2.x`.

**La solución adoptada aquí: el "punto de compatibilidad".** En lugar de retroceder
todo a 2021 (que obligaría a cambiar el dataset y la API de scanpy), se usa la
versión **más nueva de scanpy que todavía funciona con `numpy 1.23`**, de modo que
el tutorial actual de scanpy (médula ósea, con `sc.pp.scrublet` y
`leiden flavor="igraph"`) corre **tal cual** y al mismo tiempo `pyscenic` ejecuta
**GRNBoost2 real**. Un solo entorno conda (`environment.yml`):

- `python 3.10`, `numpy 1.23.5`, `pandas 1.5.3`
- `scanpy 1.10.4`, `anndata 0.10.8`
- `pyscenic 0.12.1`, `arboreto 0.1.6`, `ctxcore 0.2.0`, `dask 2023.3.2`

`scanpy 1.10` es la pieza clave: es la última serie que conserva las APIs del
tutorial actual **y** tolera `numpy 1.23` (el techo que impone pyscenic).

---

## Por qué un solo dataset (médula ósea para ambos)

Para que la integración tenga sentido, scanpy y SCENIC deben correr sobre **las
mismas células**. `01` procesa la médula ósea y exporta el loom de counts crudos;
`02` corre SCENIC sobre ese loom. Como es médula ósea humana con símbolos de gen
HGNC, `cisTarget` encuentra regulones enriquecidos contra la base de datos hg38
(imposible con datos sintéticos o con genes que no existen en la DB).

---

## Cómo reproducirlo en otra PC

```bash
conda env create -f environment.yml
conda activate scenic-medula

python 01_scanpy_clustering.py    # médula ósea → clusters + bonemarrow_for_scenic.loom
python 02_scenic_pipeline.py      # GRNBoost2 real → cisTarget → AUCell
python 03_integrate_anndata.py    # une ambos → integrated_output.loom
```

Cada script deja una carpeta con timestamp en `runs/` con un `reporte_tecnico.md`.

### Notas de ejecución

- **`pyscenic grn` es el paso lento y pesado.** GRNBoost2 sobre ~17,000 células ×
  ~23,000 genes con ~1,800 TFs puede tardar **horas** y consumir bastante memoria.
  Se controla con `SCENIC_WORKERS` (por defecto `min(4, núcleos)`).
- La primera corrida de `02` descarga ~390 MB de recursos cisTarget hg38
  (lista de TFs, `motifs.tbl` y el ranking `feather`) a `scenic_data/`. Se cachean.
