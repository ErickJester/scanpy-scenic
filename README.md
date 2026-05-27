# scanpy + SCENIC — Análisis de Single-Cell RNA-seq

## Estructura del proyecto

```
scanpy-scenic/
├── 01_scanpy_clustering.py   # Pipeline completo de scanpy
├── 02_scenic_pipeline.py     # Pipeline completo de pySCENIC
├── requirements.txt          # Dependencias
└── README.md
```

## Instalación

```bash
pip install -r requirements.txt
```

## Orden de ejecución

```bash
python 01_scanpy_clustering.py   # genera adata_clustered.h5ad
python 02_scenic_pipeline.py     # genera scenic_output.loom
```

## Fuentes oficiales

- scanpy: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
- SCENIC: https://github.com/aertslab/SCENICprotocol
