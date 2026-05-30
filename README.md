# scanpy + SCENIC — Análisis de Single-Cell RNA-seq

Une dos tutoriales oficiales complementarios: **scanpy** (qué tipos de células hay)
y **SCENIC** (qué redes regulatorias las controlan). Ver
[TUTORIALES.md](TUTORIALES.md) para la explicación de cómo se relacionan.

## Estructura del proyecto

```
scanpy-scenic/
├── 01_scanpy_clustering.py     # Pipeline de scanpy (clustering)
├── 02_scenic_pipeline.py       # Pipeline de pySCENIC (redes regulatorias)
├── 03_integrate_anndata.py     # Integración scanpy + SCENIC → loom final
├── environment-scanpy.yml      # Entorno conda para el tutorial 1 (scanpy)
├── environment-scenic.yml      # Entorno conda para el tutorial 2 (SCENIC)
├── requirements.txt            # Dependencias (entorno único, modo rápido)
├── TUTORIALES.md               # Relación entre los dos tutoriales
└── README.md
```

## Reproducir en cualquier PC (recomendado: dos entornos conda)

Los dos tutoriales fueron escritos para versiones de software incompatibles entre sí
(scanpy actual vs. SCENIC 2021), por lo que cada uno usa su propio entorno conda.
Se comunican por archivos en `runs/`, así que pueden correr por separado o en
máquinas distintas. Detalle completo en [TUTORIALES.md](TUTORIALES.md).

```bash
# --- Tutorial 1: scanpy ---
conda env create -f environment-scanpy.yml
conda activate scanpy-tutorial
python 01_scanpy_clustering.py    # → runs/scanpy_*/adata_clustered.h5ad

# --- Tutorial 2: SCENIC ---
conda env create -f environment-scenic.yml
conda activate scenic-tutorial
python 02_scenic_pipeline.py      # → runs/scenic_*/auc_matrix.csv

# --- Integración (cualquier entorno sirve) ---
python 03_integrate_anndata.py    # → runs/integrate_*/integrated_output.loom
```

Cada corrida deja una carpeta con timestamp en `runs/` con un `reporte_tecnico.md`
(métricas, criterios de validación, archivos generados).

## Alternativa: entorno único con pip

Para correr todo en un solo entorno moderno (sin GRNBoost2 real — usa una
implementación equivalente con sklearn, ver [TUTORIALES.md](TUTORIALES.md)):

```bash
pip install -r requirements.txt
```

## Fuentes oficiales

- scanpy: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
- SCENIC: https://github.com/aertslab/SCENICprotocol
