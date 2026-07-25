# Manual de usuario — Pipeline scRNA-seq: scanpy + SCENIC

> **Para quién es este manual.** Para quien va a **usar** el pipeline, no
> necesariamente a programarlo. No hace falta saber Python. Cada término técnico se
> explica la primera vez que aparece. Si algo falla, la sección 10 tiene la tabla de
> errores y soluciones.

**Producto que se entrega:** `Pipeline_scanpy_SCENIC_Lupus_Colab.ipynb` (versión
1.6.0), un notebook que se ejecuta en Google Colab desde el navegador, sin instalar
nada.

---

## Índice

1. [Qué hace y qué obtienes](#1-qué-hace-y-qué-obtienes)
2. [Antes de empezar](#2-antes-de-empezar)
3. [Puesta en marcha en cinco pasos](#3-puesta-en-marcha-en-cinco-pasos)
4. [La celda de configuración](#4-la-celda-de-configuración)
5. [Los dos modos: cuál usar](#5-los-dos-modos-cuál-usar)
6. [Cómo leer las 12 gráficas](#6-cómo-leer-las-12-gráficas)
7. [Las figuras de lupus](#7-las-figuras-de-lupus-solo-modo-zenodo)
8. [SCENIC: regulones](#8-scenic-qué-es-un-regulón-y-cómo-leerlo)
9. [Cómo saber si el resultado es entregable](#9-cómo-saber-si-el-resultado-es-entregable)
10. [Cuando algo falla](#10-cuando-algo-falla)
11. [Archivos que produce](#11-archivos-que-produce)
12. [Librerías y funciones usadas](#12-librerías-y-funciones-usadas)
13. [Glosario](#13-glosario)

---

## 1. Qué hace y qué obtienes

El pipeline parte de un experimento de **scRNA-seq** (secuenciación de ARN célula por
célula) y responde dos preguntas encadenadas:

1. **¿Qué tipos de células hay en esta muestra?** → lo resuelve **scanpy**
2. **¿Qué "interruptores" internos hacen que cada célula sea lo que es?** → lo resuelve
   **SCENIC**

### La biología mínima, en un párrafo

Imagina que el cuerpo es una ciudad y cada célula es una casa. Todas las casas tienen
la **misma biblioteca** (el ADN), pero cada una **lee libros distintos** según su
oficio. Un músculo lee los libros de músculo; una célula de defensa, los suyos. Cuando
una célula lee un gen, hace una **fotocopia** temporal (el ARN mensajero). El
experimento cuenta esas fotocopias, célula por célula.

Los **factores de transcripción** son los capataces: proteínas que deciden qué libros
se leen. Un capataz más la lista de libros que ordena leer es un **regulón**. scanpy
mira las fotocopias y agrupa células parecidas; SCENIC deduce los capataces detrás de
esas fotocopias.

### Lo que sale al final

- **12 gráficas de análisis** — desde control de calidad hasta trayectorias celulares
- **Una tabla con todos los genes evaluados** de cada grupo de células, marcando cuáles
  son estadísticamente significativos
- **Una tabla de actividad de regulones** por célula (el resultado de SCENIC)
- **Figuras específicas del estudio de lupus** (si eliges ese modo)
- **Un registro de procedencia** que declara si el resultado es apto para entregar

---

## 2. Antes de empezar

| Necesitas | Detalle |
|---|---|
| Una cuenta de Google | El notebook corre en Google Colab, gratis |
| Un navegador | No se instala nada en tu computadora |
| Tiempo | Modo `ejemplo`: ~15 min. Modo `zenodo` (dataset completo): varias horas |
| Espacio en Drive (opcional) | Si quieres conservar los resultados entre sesiones |

**Sobre la memoria (RAM).** El modo `zenodo` descarga un archivo de 1,6 GB con el
estudio completo de Jang *et al.*: **169.513 células**. Por defecto (`N_CELLS_MAX = 0`)
se usan **todas**, y eso **no cabe en el Colab gratuito** — necesitas activar RAM alta
(*Entorno de ejecución → Cambiar tipo de entorno de ejecución*) o Colab Pro/Pro+.

La celda **2.1** comprueba la RAM disponible y **se detiene con una explicación** si
detecta que no alcanza, en vez de arrancar un proceso que probablemente muera sin
aviso claro más adelante. Si solo tienes Colab gratuito, puedes poner un número
concreto en `N_CELLS_MAX` (celda 1, ej. `60000`) para trabajar con una muestra en vez
del dataset completo — ver la nota sobre este trade-off en la sección 4.

**Sobre internet.** El notebook descarga unos 2 GB de bases de datos del genoma humano
la primera vez. Si la conexión se corta a medias, lo detecta y te pide reintentar; no
continúa con un archivo incompleto.

---

## 3. Puesta en marcha en cinco pasos

1. **Abre el notebook en Colab.** Sube `Pipeline_scanpy_SCENIC_Lupus_Colab.ipynb` a
   Google Colab (*Archivo → Subir notebook*).

2. **Elige el modo.** Ve a la **celda 1** (la única que necesitas editar) y pon:

   ```python
   MODO = "ejemplo"     # o "zenodo"
   ```

3. **Ejecuta todo.** *Entorno de ejecución → Ejecutar todas*. Las celdas del modo que
   no elegiste se saltan solas.

4. **Si Colab pide reiniciar**, hazlo y vuelve a *Ejecutar todas*. Es normal: algunas
   librerías necesitan que el entorno arranque de nuevo para cargarse bien.

5. **Lee la celda 6.7 al final.** Te dice si lo que obtuviste es entregable. Esto es
   importante y se explica en la
   [sección 9](#9-cómo-saber-si-el-resultado-es-entregable).

> **Lo más importante que debes saber de entrada.** En modo `"ejemplo"` el notebook
> **se detiene a propósito** en la sección 6.5 con un mensaje explicativo. No es un
> fallo: ese conjunto de datos de demostración es demasiado pequeño para producir
> resultados biológicos reales, y el notebook prefiere decírtelo antes que entregarte
> números sin significado. Para resultados de verdad, usa `MODO = "zenodo"`.

---

## 4. La celda de configuración

Es la **celda 1**, la única que hace falta tocar.

### Elección básica

| Parámetro | Qué hace |
|---|---|
| `MODO` | `"ejemplo"` (demostración) o `"zenodo"` (datos reales de lupus) |
| `N_CELLS_MAX` | Cuántas células conservar del archivo de lupus. **Por defecto `0` = las 169.513 completas** (necesita RAM alta) |

> **Sobre `N_CELLS_MAX` y el trade-off de recortar.** El estudio de Jang *et al.* tiene
> 9 pacientes de lupus seguidos en varios momentos del tratamiento, más 8 controles
> sanos — 169.513 células en total. Si pones un número aquí (por ejemplo `60000`
> porque no tienes RAM alta), el notebook toma una muestra que mantiene las
> proporciones de **tipo celular**, pero **no** las de paciente ni las de momento del
> tratamiento. Para el análisis general (las 12 gráficas) esto no es un problema. Para
> la comparación pre/post-rituximab (sección 7 de este manual), un recorte grande
> puede dejar a algún paciente o momento con pocas células — el notebook tiene un
> mínimo de seguridad (30 células) que para la ejecución si eso ocurre, en vez de
> darte un resultado poco fiable en silencio.

### SCENIC

| Parámetro | Qué hace |
|---|---|
| `SCENIC_DOWNSAMPLE` | `True` = analiza una muestra de células (rápido). `False` = todas (horas) |
| `SCENIC_N_CELLS` | Cuántas células usar si `SCENIC_DOWNSAMPLE = True` |
| `SCENIC_N_GENES` | Cuántos genes variables conservar para el análisis de redes |
| `METODO_GRN` | `"grnboost2"` = el algoritmo real de SCENIC (**déjalo así**). `"sklearn_aprox"` = una aproximación rápida que **no** es comparable con la literatura publicada |
| `PERMITIR_REGULONES_SIN_VALIDAR` | Solo afecta a `MODO="zenodo"`. `False` (recomendado) = si SCENIC no encuentra regulones válidos, el notebook para. `True` = continúa, marcando las salidas como no validadas. **En `MODO="ejemplo"` el notebook siempre continúa**, sin importar este valor — el dataset de demostración nunca puede validar regulones por su tamaño |

### Datos y etiquetas (solo modo `zenodo`)

Todos valen `None` por defecto, lo que significa **"detéctalo tú"**. El notebook busca
la columna adecuada en los metadatos; si no la encuentra o hay varias candidatas,
**para y te muestra la lista real** para que elijas. Solo entonces rellenas estos
valores:

| Parámetro | Qué es |
|---|---|
| `COL_COND` | La columna que dice la condición o diagnóstico del paciente |
| `COL_TIME` | La columna del momento del tratamiento |
| `COL_CTYPE` | La columna con el tipo celular anotado por los autores del estudio |
| `PRE_LABEL` / `POST_LABEL` | Las etiquetas concretas de "antes" y "después" del rituximab |

### Otros

| Parámetro | Qué hace |
|---|---|
| `VERIFICAR_CHECKSUMS` | Comprueba que las bases de datos descargadas son las correctas |
| `ROOT_CLUSTER` | Grupo de partida para la gráfica de trayectoria. `None` = automático |
| `CLUSTER_KEY` | Qué nivel de detalle de agrupación usar en el resto del análisis |

---

## 5. Los dos modos: cuál usar

| | `"ejemplo"` | `"zenodo"` |
|---|---|---|
| **Datos** | Médula ósea humana + conjunto de demostración | Estudio real de lupus (Jang *et al.*, 169.513 células) |
| **Para qué sirve** | Comprobar que el pipeline corre entero | Obtener resultados interpretables |
| **Duración** | ~15 minutos | Varias horas con el dataset completo (por defecto); menos si recortas con `N_CELLS_MAX` |
| **Descarga** | ~400 MB | ~2 GB |
| **RAM** | Cabe en Colab gratuito | **Necesita RAM alta por defecto** (dataset completo) |
| **¿Da biología?** | **No** | Sí |
| **¿Llega hasta el final?** | Sí, pero SCENIC queda marcado `_SIN_VALIDAR` | Sí, resultado apto para entrega |

**Por qué el modo `ejemplo` no da biología.** Su conjunto de datos tiene 500 genes y 20
factores de transcripción. SCENIC necesita muchos más genes para detectar señales
estadísticas fiables. Con tan pocos, nunca valida ningún regulón por motivos de unión —
eso es una limitación del tamaño del dataset, no un fallo del pipeline.

**Aun así, el modo `ejemplo` completa las tres etapas de SCENIC** (GRNBoost2 → cisTarget
→ AUCell) hasta mostrar un UMAP de "regulones", para que puedas ver el mecanismo
completo funcionando de punta a punta. Esa salida queda marcada `_SIN_VALIDAR` en el
nombre del archivo y en pantalla — **no es un resultado biológico real**, es una
ilustración de cómo se ve el flujo. La celda 6.7 nunca la marcará como apta para
entrega, y eso es lo correcto: en modo `ejemplo` es imposible que lo sea.

Úsalo para verificar que todo el flujo funciona en tu entorno, de principio a fin. Para
entregar resultados reales, usa `zenodo`.

---

## 6. Cómo leer las 12 gráficas

Todas se generan en **ambos modos**, en la sección 4 del notebook.

**Gráficas 1 y 2 — Control de calidad.** No todas las "células" del experimento son
células de verdad: hay gotitas vacías y células rotas. Estas gráficas muestran cuántos
genes se detectaron por célula, cuánto ARN total, y qué porcentaje viene de las
mitocondrias. *Mucho ARN mitocondrial = célula moribunda*, y se descarta. Busca nubes
de puntos separadas del grupo principal: son candidatas a basura técnica.

**Gráfica 3 — Selección de genes informativos.** De ~23.000 genes, la mayoría son
iguales en todas las células y no ayudan a distinguir nada. El pipeline elige los 2.000
más variables. Es como fijarse en el rostro y la ropa para distinguir personas, no en
que todas tengan dos pulmones.

**Gráfica 4 — PCA.** Resume miles de mediciones en unas pocas, como la sombra de un
objeto 3D lo resume en 2D sin perder su forma. Si los grupos se separan ya aquí, la
señal biológica es fuerte.

**Gráfica 5 — Grafo de vecinos.** El mapa de parentesco: qué células se parecen a
cuáles.

**Gráfica 6 — Control de calidad sobre el mapa.** Las métricas de calidad proyectadas
sobre el mapa. **Señal de alarma:** si un grupo entero aparece coloreado como "mala
calidad", probablemente es un artefacto técnico y no un tipo celular real.

**Gráfica 7 — Grupos anotados.** El famoso dibujo de nubes de puntos (**UMAP**): cada
punto es una célula, cada nube un tipo celular. Es un **mapa de vecindarios**, no una
foto: sirve para ver quién está cerca de quién, no para medir distancias exactas.

**Gráfica 8 — Marcadores.** Genes "delatores" que identifican tipos celulares
conocidos. Si un grupo enciende `CD14` y `LYZ`, son monocitos. Así los grupos anónimos
reciben nombre biológico.

**Gráficas 9, 10 y 11 — Genes característicos.** Para cada grupo: *¿qué genes están
mucho más encendidos aquí que en el resto?* Es la huella digital de cada grupo. La
gráfica 10 exporta **dos** tablas, derivadas del mismo cálculo:

- `DE_completo_<modo>.csv` — **todos** los genes evaluados por grupo, sin recorte,
  junto con una columna `significativo` que marca cuáles de verdad pasan el corte
  estadístico (`pvals_adj < 0,05`).
- `top100_DE_<modo>.csv` — los genes **realmente significativos** de cada grupo, hasta
  un máximo de 100 (el listado de entrega pedido). El 100 es un **techo, no una cuota**:
  si un grupo tiene 8 genes significativos, el archivo trae 8 filas para ese grupo — no
  se rellena con 92 genes sin significancia real solo para completar el número. Si un
  grupo no tiene ningún gen significativo, simplemente no aparece en este archivo (pero
  sigue en `DE_completo`, por si quieres revisar su ranking igual).

> **Por qué existen los dos.** Un cupo fijo de 100 no tiene respaldo estadístico por sí
> solo — un grupo puede tener 400 genes genuinamente característicos, o solo 8. Por eso
> se conserva la tabla completa como referencia. El listado de 100 es un entregable
> pedido explícitamente, así que también se genera — pero filtrando primero por
> significancia y recortando después, para que "hasta 100" nunca signifique "100 a la
> fuerza, aunque haya que inventarlos".

> Se excluyen a propósito los genes "de mantenimiento" (mitocondriales, ribosomales y
> de hemoglobina). Aparecen como diferenciales en casi todos los grupos por razones
> técnicas, no biológicas, y taparían los marcadores que de verdad importan.

**Gráfica 12 — Trayectoria.** Reconstruye linajes: qué célula parece derivar de cuál,
ordenándolas a lo largo de un proceso de maduración (**pseudotiempo**).

---

## 7. Las figuras de lupus (solo modo `zenodo`)

**Figura 1 — Subtipos de células B.** El lupus es una enfermedad en la que el sistema
inmune ataca al propio cuerpo, y las células B tienen un papel central. Esta figura
muestra sus subtipos en la muestra.

**Figura 2 — Volcano pre vs post-rituximab.** El rituximab es un tratamiento que
elimina células B. Esta gráfica compara la expresión génica **antes y después** del
tratamiento. Cada punto es un gen: a la derecha los que suben tras el tratamiento, a la
izquierda los que bajan, y más arriba cuanto más fiable es la diferencia.

El notebook **detecta solo** cuáles son las etiquetas de "antes" y "después". Si hay
ambigüedad (por ejemplo, varios momentos post-tratamiento), para y te muestra los
valores reales para que elijas en la celda 1.

---

## 8. SCENIC: qué es un regulón y cómo leerlo

SCENIC trabaja en tres pasos. La metáfora: **somos detectives investigando quién manda
en una empresa.**

**Paso 1 — GRNBoost2: ¿quién influye sobre quién?** Si cada vez que el gen X sube
también suben Y y Z, sospechamos que X manda sobre ellos. Produce una lista enorme de
pistas.

> **Ojo:** esto es pura correlación estadística. Como sabe cualquier detective,
> *correlación no es prueba*. Por eso existe el paso 2.

**Paso 2 — cisTarget: descartar coincidencias.** Para que un factor de transcripción
controle de verdad a un gen, tiene que poder **pegarse físicamente** al ADN cerca de
ese gen, en una secuencia con cierta forma llamada **motivo**. cisTarget verifica cada
pista contra una base de datos del genoma humano y **descarta** las que no tienen dónde
pegarse. Lo que sobrevive merece llamarse **regulón**.

**Este paso es el corazón de SCENIC.** Sin él, lo que queda es una simple red de
correlaciones, que es justo lo que SCENIC existe para superar.

**Paso 3 — AUCell: ¿cuán encendido está cada regulón en cada célula?** Da una
calificación de 0 a 1 por célula y regulón. El resultado es una tabla nueva:

| | Regulón SPI1 | Regulón GATA1 | … |
|---|---|---|---|
| Célula 1 | 0,08 | 0,01 | … |
| Célula 2 | 0,00 | 0,12 | … |

Con ella se puede agrupar células **por sus programas de control** en lugar de por sus
genes sueltos.

### Cómo saber si el resultado de SCENIC tiene sentido

El notebook imprime el **rango de las puntuaciones AUC**. Si todas las células tienen
puntuaciones casi idénticas (rango menor a 0,05), te avisa: significa que los regulones
no distinguen unas células de otras, y cualquier agrupación calculada sobre eso sería
ruido, por bonito que se vea el dibujo.

---

## 9. Cómo saber si el resultado es entregable

**Esta es la sección más importante del manual.**

Que el notebook termine sin errores **no basta** para que el resultado sea válido. La
**celda 6.7** imprime un registro de procedencia con un veredicto explícito:

```
  APTO PARA ENTREGA: algoritmo estandar + regulones validados.
```

o bien:

```
  NO APTO PARA ENTREGA como resultado de SCENIC, porque ...
```

Hacen falta **las dos condiciones**:

| Condición | Por qué |
|---|---|
| `METODO_GRN = "grnboost2"` | Es el algoritmo de referencia de SCENIC. La aproximación `sklearn_aprox` no es comparable con la literatura publicada |
| Regulones validados > 0 | El filtro por motivos es lo que separa a SCENIC de una simple red de correlaciones |

**Regla práctica:** si algún archivo lleva el sufijo `_SIN_VALIDAR` en el nombre, **no
es un resultado de SCENIC**. Es una demostración del flujo, y no debe entregarse ni
interpretarse como actividad regulatoria real.

**Guarda la salida de la celda 6.7 junto con cualquier figura o tabla que entregues.**
Un gráfico separado de su contexto pierde exactamente la información que dice si se
puede interpretar.

---

## 10. Cuando algo falla

El notebook está escrito para **fallar con un mensaje que explica la causa**, en vez de
continuar en silencio. Busca aquí el mensaje que te salió:

| Mensaje | Qué pasó | Qué hacer |
|---|---|---|
| `Pediste las 169.513 celulas completas ... pero esta sesion solo tiene` | RAM insuficiente para el dataset completo | Activa RAM alta/Colab Pro, o baja `N_CELLS_MAX` (celda 1) para trabajar con una muestra |
| `Fallo la instalacion de ...` | pip no pudo instalar una librería | Lee el error de pip que aparece encima. Suele resolverse reintentando |
| `Faltan dependencias de SCENIC` | Se instaló, pero el entorno necesita reiniciarse | *Entorno de ejecución → Reiniciar sesión* → *Ejecutar todas* |
| `SHA256 no coincide` | El archivo descargado no es el esperado | Vuelve a ejecutar la celda. Si persiste, el archivo de origen cambió: avisa a quien mantiene el notebook |
| `descarga truncada` | Se cortó la conexión | Vuelve a ejecutar la celda |
| `el servidor devolvio una pagina HTML` | La dirección de descarga está caída o movida | Comprueba la URL a mano; hay que actualizarla |
| `No hay R en este entorno` | Falta R para leer los datos de lupus | Ejecuta `!apt-get install -y r-base` en una celda nueva |
| `convert_rds.R fallo` | La conversión de los datos de lupus falló | Causa habitual: falta memoria. Activa RAM alta o baja `N_CELLS_MAX` |
| `No se encontro ninguna columna para 'ctype'` | Los metadatos no tienen el nombre esperado | El error lista las columnas reales: elige la correcta y ponla en `COL_CTYPE` |
| `El filtro de celulas B encontro N celulas` | La nomenclatura de tipos celulares es distinta | El error lista las etiquetas reales; ajusta la lista de palabras clave en esa celda |
| `No se pudo determinar automaticamente el valor de 'PRE'/'POST'` | Hay varios momentos posibles, o ninguno reconocible | El error lista los valores reales: elige y ponlos en `PRE_LABEL` / `POST_LABEL` |
| `Uno de los grupos tiene solo N celulas` | Muy pocas células para una comparación fiable | Sube `N_CELLS_MAX` o elige otros momentos a comparar |
| `GRNBoost2 fallo` | Incompatibilidad de versiones con dask | El error propone el comando de instalación. **No** cambies a `sklearn_aprox` para entregables |
| `'pyscenic ctx' fallo` | El paso de validación se cayó | Si menciona `np.object`, re-ejecuta la celda 6.0 |
| `cisTarget no valido ningun regulon` | **No es un fallo técnico** | Ver abajo |
| `Ningun factor de transcripcion ... aparece en la matriz` | Los nombres de genes no coinciden con la base de datos | Suele ser símbolos vs identificadores Ensembl |

### El caso especial: "cisTarget no valido ningun regulon"

Este mensaje **no significa que algo se rompió**. Significa que ningún factor de
transcripción tiene motivos de unión enriquecidos entre sus genes candidatos.

- **En modo `ejemplo`: es lo esperado.** El conjunto de demostración es demasiado
  pequeño. Cambia a `zenodo`.
- **En modo `zenodo`:** revisa el diagnóstico que imprime la celda 6.3 — que los
  nombres de genes coincidan con la base de datos, y que la red no sea degenerada.

El notebook **para aquí a propósito**. Continuar sin ese filtro produciría puntuaciones
casi idénticas entre sí — ruido — que luego se grafican como si fueran biología.

---

## 11. Archivos que produce

| Archivo | Qué contiene |
|---|---|
| `DE_completo_<modo>.csv` | Todos los genes evaluados por grupo celular, con columna `significativo` |
| `top100_DE_<modo>.csv` | Genes significativos por grupo, hasta 100 (techo, no cuota — puede haber menos) |
| `scenic_auc_<modo>.csv` | Actividad de cada regulón en cada célula (solo cuando el resultado es válido) |
| `scenic_auc_ejemplo_SIN_VALIDAR.csv` | **Siempre se genera en modo `ejemplo`** — no es un resultado de SCENIC, es una demostración del mecanismo |
| `lupus_rituximab.h5ad` | Los datos de lupus convertidos (solo modo `zenodo`) |

Si conectaste Google Drive al principio, se guardan en `MiUnidad/scanpy_scenic_lupus/`.
Si no, se pierden al cerrar la sesión de Colab — **descárgalos antes de cerrar**.

---

## 12. Librerías y funciones usadas

Solo las esenciales, para quien quiera revisar o adaptar el código.

### Librerías

| Librería | Para qué |
|---|---|
| **scanpy** | Todo el análisis de células: calidad, agrupación, gráficas |
| **pySCENIC** | Descubrimiento de regulones |
| **arboreto** | El algoritmo GRNBoost2 real |
| **anndata** | El formato que guarda la tabla de células con sus anotaciones |
| **pandas** / **numpy** | Manejo de tablas y números |

### Funciones principales

**Análisis con scanpy:**

```python
sc.pp.calculate_qc_metrics()      # métricas de control de calidad
sc.pp.scrublet()                  # detecta "dobletes" (dos células contadas como una)
sc.pp.normalize_total()           # pone todas las células en la misma escala
sc.pp.highly_variable_genes()     # selecciona los genes informativos
sc.tl.pca()                       # reduce dimensiones
sc.pp.neighbors() / sc.tl.umap()  # construye y dibuja el mapa de células
sc.tl.leiden()                    # traza las fronteras de los grupos
sc.tl.rank_genes_groups()         # genes característicos de cada grupo (test de Wilcoxon)
sc.tl.paga() / sc.tl.dpt()        # trayectorias y pseudotiempo
```

**Redes regulatorias con SCENIC:**

```python
arboreto.algo.grnboost2()         # paso 1: infiere la red de influencias
pyscenic ctx                      # paso 2: valida contra motivos de unión (cisTarget)
pyscenic.aucell.aucell()          # paso 3: puntúa la actividad por célula
```

### Nota técnica: compatibilidad con numpy

pySCENIC es de 2022 y usa nombres que numpy eliminó después (`np.object`). El notebook
los **restaura** al arrancar, sin modificar ningún archivo de la librería. Esto evita
que el pipeline se caiga; **no** influye en los resultados.

---

## Aviso sobre los scripts `.py` de esta carpeta

Además del notebook, la carpeta contiene `01_scanpy_clustering.py` y
`02_scenic_pipeline.py`, de una versión anterior del proyecto.

**Estos scripts NO han recibido las correcciones de las versiones 1.2.0 y 1.3.0.** En
concreto:

- `02_scenic_pipeline.py` sustituye GRNBoost2 por una aproximación con scikit-learn, y
  si cisTarget no valida regulones, **los construye igualmente** desde la red sin
  validar, sin detenerse.
- `01_scanpy_clustering.py` exporta la tabla de expresión diferencial con
  `.head(50 * 15)` — un tope fijo de 750 filas aplicado a **la tabla completa de todos
  los clústeres juntos**, no por clúster. Si hay más clústeres o más genes reales de
  los que ese cálculo asumió, los últimos clústeres pueden quedar **incompletos o
  totalmente ausentes** de la tabla, sin ningún aviso. Es una versión más frágil del
  mismo problema que se corrigió en el notebook (ver gráfica 10, sección 6 de este
  manual).

Los resultados guardados en `runs/` fueron generados con esos scripts y **no deben
entregarse**. Para cualquier resultado destinado a un cliente o a publicación, usa el
notebook.

---

## 13. Glosario

**ADN** — El manual de instrucciones de la célula, idéntico en todas las células del
cuerpo.

**Gen** — Un capítulo de ese manual: las instrucciones para fabricar algo concreto.

**Expresión génica** — Que una célula esté *usando* (leyendo) un gen.

**ARN mensajero** — La fotocopia temporal de un gen. Es lo que mide el experimento.

**scRNA-seq** — La técnica que cuenta esas fotocopias célula por célula.

**Factor de transcripción** — Un interruptor que enciende o apaga muchos genes a la
vez.

**Regulón** — Un factor de transcripción **más** la lista de genes que controla.

**Motivo** — La secuencia de ADN a la que un factor de transcripción se pega
físicamente. Es la evidencia que convierte una correlación en un regulón.

**Control de calidad (QC)** — Descartar las "células" que en realidad son basura
técnica.

**Doblete** — Dos células contadas por error como una.

**Normalizar** — Poner todas las células en la misma escala para poder compararlas.

**PCA** — Resumir muchísimas medidas en unas pocas sin perder lo esencial, como una
sombra.

**UMAP** — El dibujo plano de nubes de puntos donde lo parecido queda junto.

**Clúster / Leiden** — Un grupo de células parecidas; Leiden es el algoritmo que traza
esos grupos.

**Genes marcadores** — Genes delatores que identifican un tipo celular.

**Pseudotiempo** — Un orden estimado de las células a lo largo de un proceso de
maduración.

**Volcano** — Gráfica que muestra qué genes cambian entre dos condiciones y con cuánta
fiabilidad.

**GRNBoost2 / cisTarget / AUCell** — Los tres pasos de SCENIC: adivinar influencias,
verificarlas contra el genoma, y medir cuán encendidas están.

**Checksum (SHA256)** — Una huella digital de un archivo, que permite comprobar que se
descargó completo y sin alteraciones.

---

*Manual correspondiente al notebook versión 1.6.0.*
