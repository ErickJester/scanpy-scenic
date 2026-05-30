# scanpy + SCENIC — Análisis de Single-Cell RNA-seq (tutoriales 2021)

Une dos tutoriales oficiales complementarios — **scanpy** (qué tipos de células hay)
y **SCENIC** (qué redes regulatorias las controlan) — ejecutados *tal cual el
ecosistema de 2021*, sobre el **mismo dataset PBMC3k** y en **un único entorno
conda**. Esto permite que `pyscenic grn` corra **GRNBoost2 real** y que la
integración final tenga biología real. Ver [TUTORIALES.md](TUTORIALES.md) para la
explicación de cómo se relacionan.

## Estructura del proyecto

```
scanpy-scenic/
├── 01_scanpy_clustering.py     # Tutorial 1: scanpy PBMC3k → clusters + loom
├── 02_scenic_pipeline.py       # Tutorial 2: pyscenic grn → ctx → aucell (real)
├── 03_integrate_anndata.py     # Integración scanpy + SCENIC → loom final
├── environment.yml             # Entorno conda único (época 2021)
├── requirements.txt            # Alternativa pip (conda recomendado)
├── TUTORIALES.md               # Relación entre los dos tutoriales
└── README.md
```

## Reproducir en cualquier PC

```bash
conda env create -f environment.yml
conda activate scenic-2021

python 01_scanpy_clustering.py    # PBMC3k → adata_clustered.h5ad + pbmc3k_for_scenic.loom
python 02_scenic_pipeline.py      # GRNBoost2 real → cisTarget → AUCell → auc_matrix.csv
python 03_integrate_anndata.py    # une ambos → adata_integrated.h5ad + integrated_output.loom
```

Cada corrida deja una carpeta con timestamp en `runs/` con un `reporte_tecnico.md`
(métricas, criterios de validación, archivos generados). Los pipelines se comunican
por archivos dentro de `runs/`, así que pueden correr por separado.

### Notas

- **`pyscenic grn` es el paso lento** (GRNBoost2 sobre ~2,700 células). Ajustable
  con `SCENIC_WORKERS` (por defecto `min(4, núcleos)`).
- La primera corrida de `02` descarga ~390 MB de recursos cisTarget hg38 a
  `scenic_data/` (se cachean).

## Fuentes oficiales

- scanpy (PBMC3k): https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
- SCENIC: https://github.com/aertslab/SCENICprotocol
