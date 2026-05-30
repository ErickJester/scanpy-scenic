# scanpy + SCENIC — Análisis de Single-Cell RNA-seq (dataset médula ósea)

Une dos tutoriales oficiales complementarios — **scanpy** (qué tipos de células
hay) y **SCENIC** (qué redes regulatorias las controlan) — sobre el **mismo
dataset de médula ósea humana** (~17,000 células, el del enlace actual de scanpy)
y en **un único entorno conda**. SCENIC corre con **GRNBoost2 real** y la
integración final tiene biología real. Ver [TUTORIALES.md](TUTORIALES.md) para la
explicación de cómo se relacionan.

## Estructura del proyecto

```
scanpy-scenic/
├── 01_scanpy_clustering.py     # Tutorial 1: scanpy (médula ósea) → clusters + loom
├── 02_scenic_pipeline.py       # Tutorial 2: pyscenic grn → ctx → aucell (real)
├── 03_integrate_anndata.py     # Integración scanpy + SCENIC → loom final
├── environment.yml             # Entorno conda único (punto de compatibilidad)
├── requirements.txt            # Alternativa pip (conda recomendado)
├── TUTORIALES.md               # Relación entre los dos tutoriales
└── README.md
```

## Reproducir en cualquier PC

```bash
conda env create -f environment.yml
conda activate scenic-medula

python 01_scanpy_clustering.py    # médula ósea → adata_clustered.h5ad + bonemarrow_for_scenic.loom
python 02_scenic_pipeline.py      # GRNBoost2 real → cisTarget → AUCell → auc_matrix.csv
python 03_integrate_anndata.py    # une ambos → adata_integrated.h5ad + integrated_output.loom
```

Cada corrida deja una carpeta con timestamp en `runs/` con un `reporte_tecnico.md`
(métricas, criterios de validación, archivos generados). Los pipelines se
comunican por archivos dentro de `runs/`.

### Notas

- **`pyscenic grn` es el paso lento y pesado.** GRNBoost2 sobre ~17,000 células ×
  ~23,000 genes con ~1,800 TFs puede tardar **horas**. Ajustable con
  `SCENIC_WORKERS` (por defecto `min(4, núcleos)`). Si quieres una corrida rápida
  de prueba, usa la rama con dataset PBMC3k (~2,700 células).
- La primera corrida de `02` descarga ~390 MB de recursos cisTarget hg38 a
  `scenic_data/` (se cachean).

## Fuentes oficiales

- scanpy: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
- SCENIC: https://github.com/aertslab/SCENICprotocol
