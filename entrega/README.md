# Pipeline scRNA-seq: scanpy + SCENIC

Analisis de expresion genica celula a celula sobre el estudio de lupus y
rituximab publicado en Zenodo. Agrupa las celulas por tipo con scanpy y
reconstruye con SCENIC las redes de regulacion que las controlan.

Hay dos formas de correrlo, con el mismo analisis detras:

```
entrega/
├── colab/     el cuaderno, para ejecutar en Google Colab
└── python/    la version de escritorio, para correr en el propio equipo
```

## Cual usar

**Colab** si no se quiere instalar nada. Se sube el cuaderno a
[colab.research.google.com](https://colab.research.google.com) y se ejecuta.
Todo pasa en la nube, pero la sesion se cierra sola tras unas horas de
inactividad y hay un limite de memoria segun el plan.

**Escritorio** si el equipo tiene memoria de sobra, si el analisis va a tardar
varias horas, o si se prefiere no depender de un servicio externo. Necesita
Python y, para el modo `zenodo`, tambien R. Las instrucciones completas estan en
`python/LEEME.txt`.

## Que produce

Las 12 graficas del analisis, las dos figuras del estudio (subtipos de celulas B
y comparacion antes/despues del tratamiento) y tres tablas:

| Archivo | Contenido |
|---|---|
| `DE_completo_<modo>.csv` | Todos los genes evaluados en cada grupo celular |
| `top100_DE_<modo>.csv` | Los mas caracteristicos de cada grupo, hasta 100 |
| `scenic_auc_<modo>.csv` | Actividad de cada regulon en cada celula |

Un regulon es un factor de transcripcion junto con los genes que regula. SCENIC
solo da por bueno un regulon cuando encuentra respaldo para esa relacion en su
base de datos de motivos, y ese filtro es lo que lo separa de una simple
correlacion entre genes. Si un archivo de salida lleva el sufijo `_SIN_VALIDAR`,
significa que ese filtro no dejo nada: las cifras no son interpretables como
resultado de SCENIC.

## Dos modos de datos

- `ejemplo`: medula osea publica. Corre en minutos y sirve para comprobar que
  todo esta bien instalado. No produce resultados reales: son 500 genes,
  demasiado pocos para que SCENIC valide nada.
- `zenodo`: el estudio real, 169.513 celulas. Se puede analizar una muestra
  estratificada por tipo celular en vez del total, para no necesitar tanta
  memoria.

## Que no esta en el repositorio

Las descargas (el archivo del estudio, de 1,6 GB, y las bases de datos de
SCENIC) y los resultados de cada corrida. Son varios gigabytes y se regeneran
ejecutando el pipeline, asi que no tiene sentido versionarlos.
