# Guion del profesor — Explicación completa del proyecto *scanpy + SCENIC*

> **Cómo leer este documento.** Está escrito como si yo, el profesor, estuviera de
> pie frente a ustedes explicando el proyecto de principio a fin, en voz alta, sin
> dar por sentado que saben programar. No hace falta que sepan Python. Cada vez que
> aparezca una palabra técnica, me detendré a explicarla con un ejemplo de la vida
> cotidiana. Pueden leerlo de corrido como una clase, o saltar a la sección que les
> interese usando el índice.

---

## Índice

1. [De qué trata todo esto, en una frase](#1-de-qué-trata-todo-esto-en-una-frase)
2. [La biología mínima que necesitan entender](#2-la-biología-mínima-que-necesitan-entender)
3. [Qué son los datos con los que trabajamos](#3-qué-son-los-datos-con-los-que-trabajamos)
4. [Qué es Python, una librería y un “pipeline”](#4-qué-es-python-una-librería-y-un-pipeline)
5. [Los dos tutoriales y cómo se relacionan](#5-los-dos-tutoriales-y-cómo-se-relacionan)
6. [Tutorial 1 — scanpy, paso a paso](#6-tutorial-1--scanpy-paso-a-paso)
7. [Tutorial 2 — SCENIC, paso a paso](#7-tutorial-2--scenic-paso-a-paso)
8. [La integración — juntar las dos mitades](#8-la-integración--juntar-las-dos-mitades)
9. [La “infraestructura” — qué produce cada programa](#9-la-infraestructura--qué-produce-cada-programa)
10. [El gran problema técnico: las versiones](#10-el-gran-problema-técnico-las-versiones)
11. [Las tres ramas — tres soluciones al mismo problema](#11-las-tres-ramas--tres-soluciones-al-mismo-problema)
12. [Cómo reproducir todo en otra computadora](#12-cómo-reproducir-todo-en-otra-computadora)
13. [Recorrido por cada archivo de la carpeta](#13-recorrido-por-cada-archivo-de-la-carpeta)
14. [Glosario de bolsillo](#14-glosario-de-bolsillo)
15. [Cierre](#15-cierre)

---

## 1. De qué trata todo esto, en una frase

Buenos días. Vamos a empezar por el final, por la foto completa, y después
desmenuzamos cada pieza.

Este proyecto toma **células individuales** de un tejido humano y responde dos
preguntas, una detrás de la otra:

1. **¿Qué tipos de células hay aquí?** (esto lo hace una herramienta llamada
   **scanpy**)
2. **¿Qué “interruptores” internos hacen que cada célula sea lo que es?** (esto lo
   hace una herramienta llamada **SCENIC**)

Y al final, junta las dos respuestas en un solo archivo para poder decir frases
como: *“las células de tipo monocito tienen encendido el interruptor SPI1”*.

Eso es todo. Si entienden esa frase, ya entienden el 80 % del proyecto. El otro
20 % es el **cómo**, y a eso le vamos a dedicar el resto de la clase.

---

## 2. La biología mínima que necesitan entender

No se asusten, esto es breve y lo van a usar todo el tiempo.

### 2.1 La célula y el ADN

Imaginen que el cuerpo es una **ciudad** y cada **célula** es una **casa**. Dentro
de cada casa hay una **biblioteca idéntica**: el **ADN**. En esa biblioteca están
todos los **libros de instrucciones** del organismo. A cada libro lo llamamos
**gen**.

Aquí está lo interesante: **todas las casas tienen la misma biblioteca**, pero no
todas leen los mismos libros. La casa que funciona como “panadería” lee los libros
de hacer pan; la que funciona como “bomberos” lee los libros de apagar incendios.
Tienen el mismo ADN, pero **usan genes distintos**. Por eso una neurona y una
célula de la piel son tan diferentes aunque lleven el mismo manual.

### 2.2 Expresión génica: leer un libro

Cuando una célula “lee” un gen para usarlo, decimos que ese gen se **expresa**. El
proceso físico es: la célula hace una **fotocopia** temporal del libro. Esa
fotocopia se llama **ARN mensajero (ARNm)**.

> **Idea clave:** si yo pudiera contar cuántas fotocopias de cada libro tiene una
> casa en este momento, sabría a qué se dedica esa casa. Muchas fotocopias del
> libro “hemoglobina” → es un glóbulo rojo. Muchas del libro “anticuerpos” → es una
> célula de defensa.

Contar esas fotocopias, casa por casa (célula por célula), es **exactamente** lo
que hace la máquina de laboratorio que produce nuestros datos.

### 2.3 Factores de transcripción: los interruptores

Ahora la segunda mitad. ¿Quién decide qué libros se leen en cada casa? Hay unas
proteínas especiales llamadas **factores de transcripción** (en inglés
*transcription factors*, **TF**). Piensen en ellos como **interruptores de luz** o
**capataces**: un solo factor de transcripción puede encender (o apagar) **decenas
o cientos de genes a la vez**.

Cuando un capataz (un TF) y todos los libros que ese capataz controla trabajan
juntos, a ese conjunto lo llamamos un **regulón**. Un regulón es:

> *“el factor de transcripción X, más toda la lista de genes que X enciende”.*

Si descubrimos qué regulones están activos en una célula, entendemos **el programa
de control** de esa célula, no solo lo que está haciendo, sino **por qué**.

- **scanpy** mira las fotocopias (el ARN) y agrupa células parecidas.
- **SCENIC** deduce los capataces (los regulones) detrás de esas fotocopias.

---

## 3. Qué son los datos con los que trabajamos

### 3.1 La gran tabla

El experimento de laboratorio se llama **secuenciación de ARN de célula
individual** (en inglés *single-cell RNA sequencing*, abreviado **scRNA-seq**).
Su resultado, después de mucho procesamiento de máquina, es **una tabla gigante**:

|                | Gen A | Gen B | Gen C | ... | Gen 36.601 |
|----------------|-------|-------|-------|-----|------------|
| **Célula 1**   | 0     | 5     | 2     | ... | 0          |
| **Célula 2**   | 3     | 0     | 0     | ... | 1          |
| **...**        | ...   | ...   | ...   | ... | ...        |
| **Célula 17.000** | 0  | 1     | 9     | ... | 0          |

Cada **fila** es una célula. Cada **columna** es un gen. Cada **número** es cuántas
fotocopias (cuántas moléculas de ARN) de ese gen se contaron en esa célula. La
mayoría de los números son **cero**, porque cada célula solo usa una fracción de
sus genes. A una tabla así, llena de ceros, los matemáticos la llaman **matriz
dispersa**, y es importante porque nos permite guardarla sin ocupar un disco
entero.

### 3.2 Los dos conjuntos de datos que aparecen en el proyecto

A lo largo del proyecto verán mencionados **dos juegos de células**. Esto importa
porque las distintas versiones del proyecto (las “ramas”, ya llegaremos) usan uno u
otro:

- **Médula ósea humana** (*bone marrow*): unas **17.000 células** y **36.601
  genes**. Es el conjunto que muestra hoy la página oficial de scanpy. Es grande y
  realista.
- **PBMC3k**: unas **2.700 células** de sangre (PBMC = células mononucleares de
  sangre periférica, básicamente glóbulos blancos). Es pequeño y rápido. Era el
  ejemplo “clásico” de scanpy hacia el año 2021.

Los dos vienen en un formato de archivo de la empresa **10X Genomics** (archivos
`.h5`), que es simplemente una forma comprimida de guardar esa tabla gigante.

---

## 4. Qué es Python, una librería y un “pipeline”

Esta sección es para quienes nunca han programado. Si ya saben, sáltenla.

### 4.1 Python

**Python** es un **idioma para darle órdenes a la computadora**. En lugar de hacer
clic en botones, uno escribe **instrucciones en texto**, una debajo de otra, en un
archivo. A ese archivo de instrucciones se le llama **script** o **programa**. Los
archivos de Python terminan en `.py`. En esta carpeta hay tres:
`01_scanpy_clustering.py`, `02_scenic_pipeline.py` y `03_integrate_anndata.py`.

Cuando uno “ejecuta” un script, la computadora lee las instrucciones de arriba
hacia abajo y las va cumpliendo, como una receta de cocina.

### 4.2 Librería (o “paquete”)

Nadie escribe todo desde cero. Si yo quiero hacer un pastel, no cultivo el trigo:
compro harina ya hecha. En programación, esa “harina ya hecha” es una **librería**:
un paquete de instrucciones que alguien más ya escribió y que yo puedo **usar**.

Las librerías estrella de este proyecto son:

- **scanpy** — caja de herramientas para analizar células individuales.
- **pySCENIC** — caja de herramientas para descubrir regulones.
- **numpy** y **pandas** — las dos librerías básicas para manejar tablas y números
  (la “harina y el azúcar” de casi todo en Python científico).
- Otras de apoyo: **anndata** (guarda la tabla de células), **loompy** (otro
  formato de archivo), **matplotlib** y **seaborn** (dibujan las gráficas).

Recuerden estos nombres porque el **gran conflicto** del proyecto (sección 10) será
entre **versiones** de estas librerías.

### 4.3 Pipeline

Un **pipeline** (literalmente “tubería”) es una **cadena de pasos** donde la salida
de uno es la entrada del siguiente, como una **línea de montaje** en una fábrica:
entra la materia prima por un lado, pasa por estaciones, y sale el producto
terminado por el otro. Nuestro proyecto es una línea de montaje de tres estaciones:

```
   datos crudos  →  [01 scanpy]  →  [02 SCENIC]  →  [03 integración]  →  resultado
```

---

## 5. Los dos tutoriales y cómo se relacionan

Todo el proyecto nace de **dos tutoriales oficiales**, dos “recetas” publicadas por
los autores de las herramientas:

| | Tutorial 1 | Tutorial 2 |
|---|---|---|
| **Herramienta** | scanpy | SCENIC |
| **Página oficial** | scanpy.readthedocs.io (clustering) | github.com/aertslab/SCENICprotocol |
| **Pregunta** | ¿Qué tipos de células hay? | ¿Qué regulones las controlan? |
| **Archivo en el proyecto** | `01_scanpy_clustering.py` | `02_scenic_pipeline.py` |

La relación entre ellos es **secuencial**: no son alternativas, son dos mitades de
la misma historia. Primero scanpy **ordena** las células en grupos; después SCENIC
**explica** qué programa interno define a cada grupo. Por eso tiene tanto sentido
correr los dos sobre **las mismas células** y luego unir los resultados.

Guarden esta idea, porque más adelante (sección 11) vamos a ver que hacer que esos
dos tutoriales convivan en la misma computadora resultó ser **el verdadero reto de
ingeniería** del proyecto.

---

## 6. Tutorial 1 — scanpy, paso a paso

Vamos a abrir el primer programa, `01_scanpy_clustering.py`, y a recorrer sus pasos
como si yo los narrara en voz alta. No voy a mostrar código; voy a explicar **qué
hace y por qué** cada paso. (Los números de paso son los mismos que verán como
comentarios dentro del archivo.)

### Paso 1 — Cargar los datos

El programa descarga la tabla gigante de células (la médula ósea, o PBMC3k, según
la versión) y la mete en la memoria de la computadora dentro de un objeto llamado
**AnnData**. Piensen en AnnData como **una hoja de cálculo muy inteligente**: en el
centro tiene la tabla de números (células × genes) y, alrededor, etiquetas y notas
sobre cada fila y cada columna.

### Paso 2 — Control de calidad (QC)

No todas las “células” del experimento son células de verdad. Algunas son **basura
técnica**: gotitas vacías, células rotas, restos. El control de calidad es el
**portero** que decide quién entra y quién no. Se fija en tres señales:

- **Cuántos genes distintos** se detectaron en la célula (muy pocos = probablemente
  una gota vacía).
- **Cuánto ARN total** tiene.
- **Qué porcentaje del ARN viene de las mitocondrias.** Esto es un truco clásico:
  cuando una célula se está muriendo, su contenido se escapa pero el de las
  mitocondrias se queda, así que **mucho ARN mitocondrial = célula moribunda**, y la
  descartamos.

En nuestra corrida real con médula ósea, entraron 17.125 células y, tras el filtro,
quedaron **17.041** (descartamos 84). También se eliminan genes que casi no
aparecen.

### Paso 3 — Detección de “dobletes”

A veces la máquina mete **dos células en la misma gotita** y las cuenta como una.
A ese fantasma lo llamamos **doblete**, y hay que detectarlo porque es una célula
falsa con una mezcla de dos identidades. Una herramienta llamada **Scrublet**
estima la probabilidad de que cada célula sea en realidad un doblete. En nuestra
corrida marcó un 1,2 % de dobletes.

### Paso 4 — Normalización

Unas células son grandes y otras pequeñas, así que naturalmente tienen más o menos
ARN total. Si no corrigiéramos eso, confundiríamos “célula grande” con “célula
distinta”. **Normalizar** es poner a todas las células en la **misma escala**, como
cuando convertimos precios de distintas monedas a dólares para poder compararlos.
Después se aplica una transformación matemática suave (`log`) para que un puñado de
genes hiperactivos no “griten” por encima de todos los demás.

> Antes de normalizar, el programa guarda una **copia de los números crudos** en una
> gaveta llamada `counts`. Esto será **crucial** para el Tutorial 2, que necesita los
> números originales, no los normalizados.

### Paso 5 — Genes altamente variables

De los ~23.000 genes, la mayoría son iguales en todas las células (aburridos, no
ayudan a distinguir tipos). El programa elige los **2.000 genes más variables**, los
que de verdad marcan diferencias. Es como, para distinguir personas, fijarse en el
rostro y la ropa y no en que “todos tienen dos pulmones”.

### Paso 6 — PCA (reducción de dimensiones)

Tenemos 2.000 genes, o sea **2.000 dimensiones**. Nadie puede imaginar un espacio de
2.000 dimensiones. **PCA** (Análisis de Componentes Principales) es un truco para
**resumir** esas 2.000 medidas en unas pocas decenas que capturan lo esencial, igual
que la **sombra** de un objeto 3D sobre una pared lo resume en 2D sin perder su
forma reconocible. Reduce ruido y acelera todo lo que viene después.

### Paso 7 — Grafo de vecinos y UMAP

Ahora el programa calcula, para cada célula, **cuáles son sus células más
parecidas** (sus “vecinas”). Con eso construye un **mapa de parentesco**.

**UMAP** es la técnica que **dibuja ese mapa en una hoja plana** (2D), de modo que
células parecidas queden cerca y células distintas queden lejos. Es el famoso
“dibujo de nubes de puntos” que han visto en los pósters de biología: cada punto es
una célula, cada nube es un tipo celular. Es un **mapa**, no una foto literal:
sirve para ver vecindarios, no para medir distancias exactas.

### Paso 8 — Clustering con Leiden

Hasta aquí tenemos un mapa, pero nadie ha **trazado las fronteras** de los
vecindarios. El algoritmo **Leiden** hace justo eso: agrupa automáticamente las
células en **clústeres** (grupos) según quién está conectado con quién. Es como
mirar un mapa de ciudades de noche desde un avión y dibujar un círculo alrededor de
cada mancha de luces.

Leiden tiene una perilla llamada **resolución**: baja = pocos grupos grandes; alta =
muchos grupos pequeños. El programa prueba tres resoluciones a propósito (en nuestra
corrida dieron **5, 17 y 36** grupos) para que el biólogo elija el nivel de detalle
que necesita.

### Paso 9 — Anotar tipos celulares

Un clúster, por sí solo, es solo “grupo número 3”. Para ponerle **nombre biológico**
(“esto son monocitos”), el programa usa **genes marcadores**: genes que se sabe que
solo encienden ciertos tipos de célula. Si el grupo 3 enciende `CD14` y `LYZ`, es un
monocito. Así, los números anónimos se convierten en nombres: *Lymphocytes,
Monocytes, Erythroid, B Cells…*

### Paso 10 — Expresión diferencial

Por último, el programa pregunta, para cada grupo: *“¿qué genes están mucho más
encendidos aquí que en el resto?”*. Esa lista de genes característicos (una prueba
estadística llamada **Wilcoxon**) es la **huella digital** de cada grupo, y se
guarda en una tabla.

### El resultado del Tutorial 1

Un archivo `adata_clustered.h5ad` que contiene **todo**: las células, sus grupos,
sus nombres, el mapa UMAP y las huellas digitales. Esa es la materia prima del
siguiente tutorial.

---

## 7. Tutorial 2 — SCENIC, paso a paso

Ahora abrimos `02_scenic_pipeline.py`. Cambiamos de pregunta: ya no es *“¿qué
células hay?”* sino *“¿qué interruptores (regulones) las controlan?”*. SCENIC lo
resuelve en tres grandes pasos. Voy a explicarlos con una metáfora única que iremos
ampliando: **somos detectives investigando quién manda en una empresa.**

### Paso A — GRNBoost2: ¿quién influye sobre quién?

Tenemos los niveles de actividad de todos los genes en todas las células. La
pregunta de este paso es: *“cuando el gen-jefe X sube o baja, ¿qué otros genes suben
o bajan con él?”*. Si cada vez que X sube, también suben Y y Z, sospechamos que **X
manda sobre Y y Z**.

La herramienta **GRNBoost2** hace exactamente esa correlación masiva entre los
**factores de transcripción** (los posibles jefes) y todos los demás genes. Produce
una lista enorme de pistas del estilo *“X probablemente influye sobre Y, con fuerza
0,37”*. A esa lista la llamamos **red de regulación génica** (las “adyacencias”).

> **Ojo, detalle importante para después:** este paso es **pura correlación
> estadística**. Como todo detective sabe, *correlación no es prueba*. Por eso existe
> el paso B.

### Paso B — cisTarget: descartar las coincidencias

Que X e Y suban juntos puede ser **casualidad**. El paso B busca **evidencia
física**. Resulta que, para que un factor de transcripción controle de verdad a un
gen, tiene que poder **pegarse** físicamente al ADN cerca de ese gen, en una
secuencia con cierta forma llamada **motivo**.

La herramienta **cisTarget** toma cada pista del paso A y la **verifica** contra una
gran **base de datos** de motivos del genoma humano (los archivos pesados que el
programa descarga, de ~390 MB). Si el factor X **no tiene dónde pegarse** cerca de
Y, se descarta la pista por más que estén correlacionados. Lo que sobrevive a esta
depuración ya merece el nombre de **regulón** de verdad.

> **Aquí ocurrió una lección real del proyecto.** En una de nuestras primeras
> corridas usamos datos **inventados** (genes con nombres falsos como `GENE0001`).
> Como esos genes **no existen** en la base de datos del genoma humano, cisTarget no
> encontró **ningún** motivo y devolvió **cero regulones**. No era un error del
> programa: era que los datos no eran reales. La solución fue usar **datos humanos de
> verdad** (PBMC, médula ósea), donde los genes sí existen en la base de datos.

### Paso C — AUCell: ¿qué tan encendido está cada regulón en cada célula?

Ya sabemos **qué** regulones existen. Falta medir **cuánto** está activo cada uno en
**cada célula**. La herramienta **AUCell** revisa, célula por célula, si los genes de
un regulón están entre los más activos de esa célula, y le pone una **calificación
de 0 a 1** (el “score AUC”). Resultado: una **nueva tabla**, ya no de genes, sino de
**regulones × células**:

|              | Regulón SPI1 | Regulón GATA1 | ... |
|--------------|--------------|---------------|-----|
| Célula 1     | 0,08         | 0,01          | ... |
| Célula 2     | 0,00         | 0,12          | ... |

Esa tabla (`auc_matrix.csv`) es el producto estrella del Tutorial 2. Con ella se
puede dibujar otro UMAP, pero esta vez agrupando células **por sus programas de
control** en lugar de por sus genes sueltos.

### El resultado del Tutorial 2

La tabla de actividad de regulones por célula, más un mapa de esas actividades. En
otras palabras: **la radiografía del “panel de control” de cada célula.**

---

## 8. La integración — juntar las dos mitades

Llega el tercer programa, `03_integrate_anndata.py`. Su trabajo es simple de
enunciar y muy valioso: **pegar** la respuesta del Tutorial 1 (qué tipo es cada
célula) con la del Tutorial 2 (qué regulones tiene encendidos), célula por célula.

¿Cómo sabe que “la célula 5 de scanpy” es “la célula 5 de SCENIC”? Porque cada
célula tiene un **código de barras** único (un identificador de texto). El programa
**cruza** las dos tablas por ese código de barras, como cuando se cruzan dos listas
de invitados por el número de cédula.

El resultado permite, por fin, afirmaciones biológicas completas:

> *“El grupo que scanpy llamó **monocitos** tiene encendido el regulón **SPI1**;
> el grupo de **glóbulos rojos** tiene encendido **GATA1**.”*

Eso es **conocimiento nuevo**: no solo *qué* células hay, sino *qué programa interno*
las define. El programa produce dos archivos finales: `adata_integrated.h5ad` (todo
junto, para seguir analizando) y `integrated_output.loom` (un formato pensado para
visores interactivos).

> **Una advertencia que aprendimos a la fuerza.** El cruce por código de barras solo
> funciona si **ambos tutoriales corrieron sobre las mismas células**. Si scanpy usó
> un conjunto y SCENIC otro distinto, **ningún** código de barras coincide y la
> integración queda vacía (rellena con “N/A”). Esto explica una de las decisiones de
> diseño que verán en las ramas: usar **un solo conjunto de datos para toda la
> cadena**.

---

## 9. La “infraestructura” — qué produce cada programa

Hay un detalle de ingeniería del que estoy orgulloso y que conviene que entiendan,
porque hace al proyecto **serio y reproducible**.

Cada vez que se ejecuta cualquiera de los tres programas, este **no ensucia** la
carpeta tirando archivos por todos lados. En cambio:

1. Crea una **carpeta nueva con la fecha y la hora** en su nombre, dentro de `runs/`.
   Por ejemplo: `runs/scanpy_2026-05-30_13-36-39/`. Así, **cada corrida queda
   archivada por separado** y nunca se pisan los resultados de ayer con los de hoy.
2. Dentro guarda:
   - los **datos de resultado** (las tablas, los `.h5ad`, los `.loom`);
   - una subcarpeta `figures/` con todas las **gráficas** en PNG;
   - y, lo más bonito, un **reporte técnico** en un archivo `reporte_tecnico.md`.

Ese **reporte** es un informe automático, legible por humanos, que dice: en qué
computadora se corrió, qué versión de cada librería se usó, cuánto tardó cada paso,
cuántas células quedaron, cuántos grupos se formaron… y una **tabla de criterios de
validación** con palomitas ✓ o cruces ✗ que permite, de un vistazo, saber si la
corrida salió bien o si algo falló. Es como el **recibo detallado** que entrega un
laboratorio serio junto con los resultados.

Hay otras dos carpetas que conviene nombrar:

- `scenic_data/` — donde se **guardan (caché)** los archivos pesados que SCENIC
  descarga de internet (~390 MB), para no bajarlos cada vez.
- `venv/` — un “Python privado” del proyecto (lo explico en la sección 12).

Estas tres carpetas (`runs/`, `scenic_data/`, `venv/`) están en una **lista de
ignorados** (`.gitignore`) para que **no se suban** al repositorio compartido: son
pesadas y se regeneran solas.

---

## 10. El gran problema técnico: las versiones

Aquí viene la parte que de verdad distingue este proyecto, y la razón de que existan
**tres ramas**. Pongan atención porque es el corazón de la historia.

### 10.1 Qué es una “versión”

Las librerías (scanpy, numpy, pandas…) **cambian con el tiempo**. Sale la versión
1.9, luego la 1.10, la 2.0… A veces, una versión nueva **rompe** cosas que
funcionaban en la vieja: le cambian el nombre a una función, le quitan una opción,
cambian cómo se guarda algo. Es como cuando una aplicación del teléfono se
“actualiza” y de repente el botón que usabas ya no está donde estaba.

### 10.2 El choque concreto

Tenemos dos tutoriales con **edades distintas**:

- El tutorial de **scanpy** que figura en su web es **moderno** (año 2024). Usa
  funciones nuevas que **exigen versiones recientes** de las librerías de base
  (numpy 1.26, pandas 2.x).
- El protocolo de **SCENIC** es de **2021**. Su motor (las piezas `arboreto`,
  `dask`, y el propio `pyscenic`) fue escrito para versiones **viejas** de esas mismas
  librerías. Con las modernas, **se rompe**: por ejemplo, el código de SCENIC pide
  algo llamado `np.object`, que las versiones nuevas de numpy **eliminaron** en 2023.

Aquí está el nudo: **scanpy moderno quiere librerías nuevas; SCENIC quiere librerías
viejas. No se puede tener las dos cosas a la vez en el mismo entorno.** Es como pedir
que un mismo enchufe sea de 110 y de 220 voltios simultáneamente.

> Durante el proyecto vimos estos choques **en vivo**: un error que decía
> `No module named 'pkg_resources'`, otro `module 'numpy' has no attribute 'object'`,
> y una dirección de internet caída (error 404) de donde SCENIC intentaba bajar sus
> datos de ejemplo. Cada uno fue un síntoma del mismo problema de fondo: **piezas de
> épocas distintas que no encajan.**

### 10.3 Las tres maneras de resolverlo

Frente a ese nudo hay tres estrategias razonables, y **cada rama del proyecto
implementa una**:

1. **Parchar el presente.** Quedarse en el mundo moderno y **reemplazar** la pieza de
   SCENIC que no funciona por una equivalente que sí. *(Rama `main`.)*
2. **Volver entero al pasado.** Retroceder **todo** al ecosistema de 2021, de modo
   que SCENIC corra tal cual, original. *(Rama `reproducibilidad-tutoriales`.)*
3. **Buscar el punto medio exacto.** Encontrar la versión **más nueva posible** que
   todavía sea compatible con el SCENIC viejo, para conservar el tutorial moderno de
   scanpy y a la vez correr el SCENIC original. *(Rama `opcion-3-medula-osea`.)*

Vamos a ver cada una con calma.

---

## 11. Las tres ramas — tres soluciones al mismo problema

### 11.0 ¿Qué es una “rama” (branch)?

Antes que nada: el proyecto usa **git**, un sistema que guarda el **historial** de
todos los cambios, como un “control de versiones” con máquina del tiempo. Una
**rama** es una **línea de trabajo paralela**: una copia del proyecto donde puedo
probar una idea **sin afectar** la versión principal. Imaginen un libro de
*“elige tu propia aventura”*: desde un punto común, la historia se bifurca en
caminos distintos. Aquí hay tres caminos:

- `main` — el tronco principal.
- `reproducibilidad-tutoriales` — el camino “volver a 2021”.
- `opcion-3-medula-osea` — el camino “punto de compatibilidad”.

Las tres comparten el mismo origen y la misma idea general; se diferencian en **cómo
resuelven el choque de versiones** y en **qué datos usan**.

### 11.1 Rama `main` — “parchar el presente”

Es la primera versión que funcionó de punta a punta. Su filosofía: **no pelear con
las versiones modernas**. Como el motor original de SCENIC (`arboreto`) no corre con
las librerías nuevas, esta rama lo **sustituye** por una pieza equivalente de otra
librería muy común (`scikit-learn`) que hace un trabajo parecido: deducir qué genes
influyen sobre cuáles.

- **Ventaja:** todo vive en un solo entorno moderno; scanpy corre sin problemas.
- **Desventaja:** el paso de redes **ya no es el GRNBoost2 “oficial”** del tutorial,
  sino un reemplazo. Es decir, **no es SCENIC “tal cual”**.
- **Datos:** médula ósea para scanpy; datos de ejemplo para SCENIC.

Piénsenlo como restaurar un coche antiguo poniéndole un **motor moderno**: arranca y
anda, pero los puristas dirán que “ya no es el original”.

### 11.2 Rama `reproducibilidad-tutoriales` — “volver entero a 2021” (Opción 1)

Aquí la filosofía es la opuesta: **respetar al pie de la letra el protocolo SCENIC
original**, aunque para ello haya que **retroceder todo el entorno** al año 2021.

Como el tutorial moderno de scanpy no existiría en 2021, esta rama usa la **versión
de 2021 del tutorial de scanpy**, que trabajaba con el conjunto pequeño **PBMC3k**.
Y fija (“congela”) todas las librerías en sus versiones de aquella época: Python
3.10, scanpy 1.9, numpy 1.23, pandas 1.5, y el SCENIC original con su **GRNBoost2 de
verdad**.

- **Ventaja:** los **dos** tutoriales corren **originales, sin reemplazos**, en un
  **único entorno** coherente.
- **Ventaja extra:** como scanpy y SCENIC usan **el mismo conjunto PBMC3k**, los
  códigos de barras coinciden y la **integración funciona con biología real**.
- **Desventaja:** el tutorial de scanpy ya no es el que hoy muestra la web (es el de
  2021), y se trabaja con menos células.
- **Datos:** **PBMC3k para toda la cadena**.

Es la restauración **purista**: se le consigue al coche antiguo su **motor original**
y se acepta vivir con la tecnología de su época.

### 11.3 Rama `opcion-3-medula-osea` — “el punto de compatibilidad” (Opción 3)

Esta es la solución **más fina**, un punto medio cuidadosamente calculado. La idea:
*¿cuál es la versión **más nueva** de scanpy que **todavía** se lleva bien con el
SCENIC viejo?* La respuesta resultó ser **scanpy 1.10** junto con **numpy 1.23**.

¿Por qué ese par exacto? Porque:

- numpy debe ser **anterior a la 1.24** (si no, SCENIC se rompe con el famoso
  `np.object`);
- y scanpy 1.10 es la **última** serie que conserva las funciones modernas que usa
  el tutorial **actual** de médula ósea (la detección de dobletes, el clustering
  nuevo) **y** que aún tolera ese numpy 1.23.

Así, esta rama logra lo mejor de los dos mundos: conserva el **tutorial moderno de
scanpy** (médula ósea, ~17.000 células, tal como la web) **y** corre el **SCENIC
original con GRNBoost2 de verdad**, todo en un solo entorno.

- **Ventaja:** tutorial de scanpy **actual** + SCENIC **original**, juntos.
- **Desventaja:** SCENIC sobre 17.000 células es **lento y pesado** (puede tardar
  **horas** y consumir bastante memoria).
- **Datos:** **médula ósea para toda la cadena.** El programa 01 fue ajustado para
  **exportar un archivo puente** (`bonemarrow_for_scenic.loom`) con los números
  crudos y los mismos códigos de barras, de modo que SCENIC corra sobre **las mismas
  células** y la integración vuelva a funcionar.

Es el restaurador **experto** que encuentra la pieza más nueva que aún encaja en el
coche antiguo: máxima modernidad **sin** romper la compatibilidad.

### 11.4 Tabla comparativa de las tres ramas

| | `main` | `reproducibilidad-tutoriales` (Op. 1) | `opcion-3-medula-osea` (Op. 3) |
|---|---|---|---|
| **Filosofía** | Parchar el presente | Volver entero a 2021 | Punto de compatibilidad |
| **Motor de redes (SCENIC)** | Reemplazo con scikit-learn | GRNBoost2 **original** | GRNBoost2 **original** |
| **¿SCENIC “tal cual”?** | No | Sí | Sí |
| **Tutorial de scanpy** | Moderno (médula ósea) | El de 2021 (PBMC3k) | Moderno (médula ósea) |
| **Datos de la cadena** | Médula ósea / ejemplo | **PBMC3k** (uno solo) | **Médula ósea** (uno solo) |
| **scanpy** | 1.11 | 1.9 | 1.10 |
| **numpy** | 1.26 | 1.23 | 1.23 |
| **Entornos** | Uno (moderno) | Uno (2021) | Uno (compatibilidad) |
| **Velocidad de SCENIC** | Rápida | Rápida (pocas células) | Lenta (muchas células) |
| **Integración con biología real** | Limitada | Sí | Sí |

> **¿Cuál es “la buena”?** Ninguna es “la correcta” a secas; cada una sirve para un
> propósito. Si quieren **fidelidad total al protocolo SCENIC** y una corrida rápida,
> usen la **Opción 1**. Si quieren el **tutorial de scanpy tal como está hoy en la
> web** con SCENIC original, usen la **Opción 3**. Si solo quieren algo que **arranque
> en un entorno moderno** sin complicaciones, `main`.

---

## 12. Cómo reproducir todo en otra computadora

“Reproducible” significa que **otra persona, en otra máquina, obtiene los mismos
resultados**. Es uno de los pilares de la ciencia seria, y este proyecto se tomó el
trabajo de lograrlo. Les explico las dos piezas que lo hacen posible.

### 12.1 El entorno aislado

Si yo instalo las librerías “sueltas” en mi computadora, tarde o temprano chocan con
las de otro proyecto. La solución son los **entornos aislados**: una **caja
cerrada** donde viven solo las librerías de **este** proyecto, en sus versiones
exactas, sin molestar a nadie afuera. Hay dos tecnologías para esa caja:

- **venv** — la caja sencilla que ya está en la carpeta (`venv/`). Sirve para correr
  con librerías modernas (estilo `main`).
- **conda** — una caja más poderosa, capaz de instalar también piezas que no son de
  Python. Es la **recomendada** para las ramas de SCENIC, porque varias de sus
  librerías traen componentes delicados que conda maneja mejor.

### 12.2 La “receta” del entorno

Cada rama de SCENIC trae un archivo **`environment.yml`**: la **lista de la compra**
con las versiones exactas de cada librería. Con un solo comando, conda **lee esa
lista y arma la caja idéntica** en cualquier computadora. Eso es lo que garantiza
que en su máquina pase lo mismo que en la mía.

En las ramas de SCENIC, los pasos son:

```bash
conda env create -f environment.yml     # arma la caja a partir de la receta
conda activate <nombre-del-entorno>      # entra en la caja

python 01_scanpy_clustering.py           # tipos de células  → archivos en runs/
python 02_scenic_pipeline.py             # regulones          → archivos en runs/
python 03_integrate_anndata.py           # une ambos          → resultado final
```

(En la rama Opción 1 el entorno se llama `scenic-2021`; en la Opción 3, `scenic-medula`.)

Y recuerden: **cada corrida deja su carpeta con fecha en `runs/` y su
`reporte_tecnico.md`**, así que siempre queda constancia de qué se hizo y cómo salió.

---

## 13. Recorrido por cada archivo de la carpeta

Para cerrar el inventario, aquí está **qué es cada cosa** que verán al abrir la
carpeta. Lo que está marcado “según la rama” cambia entre las tres versiones.

| Archivo o carpeta | Qué es |
|---|---|
| `01_scanpy_clustering.py` | Programa del **Tutorial 1** (scanpy): de células crudas a grupos con nombre. |
| `02_scenic_pipeline.py` | Programa del **Tutorial 2** (SCENIC): de la expresión a los regulones. |
| `03_integrate_anndata.py` | Programa de **integración**: une los resultados de 01 y 02. *(existe en las ramas de trabajo)* |
| `environment.yml` | La **receta** del entorno conda (versiones exactas). *(en las ramas de SCENIC)* |
| `requirements.txt` | Lista de librerías para la caja sencilla (`venv`/pip). |
| `README.md` | La **portada**: resumen e instrucciones rápidas. |
| `TUTORIALES.md` | Documento que explica **la relación entre los dos tutoriales**. *(en las ramas de trabajo)* |
| `GUION_PROFESOR.md` | **Este documento.** |
| `.gitignore` | La lista de cosas que **no** se suben al repositorio (carpetas pesadas). |
| `runs/` | Carpeta donde cada corrida deja sus **resultados, gráficas y reporte**. *(se genera sola)* |
| `scenic_data/` | **Caché** de los archivos pesados que SCENIC descarga. *(se genera sola)* |
| `venv/` | El **Python privado** del proyecto. *(se genera sola)* |

---

## 14. Glosario de bolsillo

- **ADN:** el manual de instrucciones de la célula. Igual en todas las células del
  cuerpo.
- **Gen:** un “capítulo” del manual; las instrucciones para fabricar algo concreto.
- **Expresión génica:** que una célula esté **usando** (leyendo) un gen.
- **ARN mensajero (ARNm):** la **fotocopia** temporal de un gen que la célula usa.
  Es lo que el experimento mide.
- **Factor de transcripción (TF):** un **interruptor** que enciende o apaga muchos
  genes a la vez.
- **Regulón:** un factor de transcripción **más** la lista de genes que controla.
- **scRNA-seq:** la técnica de laboratorio que cuenta las fotocopias de ARN **célula
  por célula**.
- **Matriz / tabla células × genes:** la gran tabla de números que es la materia
  prima de todo.
- **Célula / casa, tipo celular:** un grupo de células que hacen el mismo oficio
  (monocitos, glóbulos rojos, etc.).
- **Control de calidad (QC):** descartar las “células” que en realidad son basura
  técnica.
- **Doblete:** dos células contadas por error como una.
- **Normalizar:** poner a todas las células en la misma escala para poder
  compararlas.
- **PCA:** resumir muchísimas medidas en unas pocas sin perder lo esencial (como una
  sombra).
- **UMAP:** el dibujo plano de nubes de puntos donde lo parecido queda junto.
- **Clúster / Leiden:** un grupo de células parecidas; Leiden es el algoritmo que
  traza esos grupos.
- **Genes marcadores:** genes “delatores” que identifican un tipo celular.
- **GRNBoost2:** el paso de SCENIC que adivina **quién influye sobre quién**.
- **cisTarget:** el paso que **verifica** esas influencias contra el genoma y
  descarta las casualidades.
- **AUCell:** el paso que mide **cuán encendido** está cada regulón en cada célula.
- **Python:** el idioma para darle órdenes a la computadora.
- **Script / programa (`.py`):** un archivo con esas órdenes, ejecutadas en orden.
- **Librería / paquete:** instrucciones ya hechas que uno reutiliza (scanpy, numpy…).
- **Pipeline:** una cadena de pasos tipo línea de montaje.
- **Versión:** la “edición” de una librería; ediciones distintas a veces no encajan
  entre sí.
- **Entorno (venv / conda):** una caja aislada con las librerías de **este** proyecto
  en sus versiones exactas.
- **`environment.yml`:** la receta para armar esa caja idéntica en otra máquina.
- **git / rama (branch):** sistema de historial; una rama es una línea de trabajo
  paralela.
- **AnnData (`.h5ad`) / loom (`.loom`):** formatos de archivo que guardan la tabla de
  células con todas sus anotaciones.
- **Reproducible:** que otra persona, en otra máquina, obtenga los mismos resultados.

---

## 15. Cierre

Recapitulemos la clase en cinco ideas, que es lo que quiero que se lleven:

1. **El propósito.** Partimos de células individuales y respondemos dos preguntas
   encadenadas: *qué* tipos de células hay (scanpy) y *por qué* son así, qué
   regulones las controlan (SCENIC).
2. **La cadena.** Tres programas en línea de montaje: `01` agrupa, `02` descubre
   regulones, `03` los une. Cada uno deja un informe automático con fecha en `runs/`.
3. **El reto real.** Los dos tutoriales nacieron en **épocas distintas** y sus
   librerías **no encajan** entre sí. Resolver ese choque fue el verdadero trabajo de
   ingeniería.
4. **Las tres soluciones (ramas).** `main` parcha el presente (SCENIC con reemplazo);
   `reproducibilidad-tutoriales` vuelve entero a 2021 con PBMC3k (SCENIC original);
   `opcion-3-medula-osea` halla el punto de compatibilidad para conservar el tutorial
   moderno **y** el SCENIC original sobre médula ósea.
5. **La reproducibilidad.** Gracias a los entornos aislados y a las recetas
   `environment.yml`, cualquiera puede rearmar todo en su computadora y obtener los
   mismos resultados.

Si entendieron la **metáfora de la ciudad de casas con bibliotecas idénticas pero
que leen libros distintos según su oficio, y de los capataces que deciden qué libros
se leen**, entonces entendieron la biología. Y si entendieron que **el problema más
difícil no fue la biología sino lograr que herramientas de distintas épocas
trabajaran juntas**, entonces entendieron la ingeniería. Con esas dos ideas, el
resto del proyecto es detalle.

Eso es todo por hoy. Gracias.
