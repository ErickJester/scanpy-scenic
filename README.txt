============================================================
README - Pipeline scanpy + SCENIC (notebook de Colab)
============================================================

Este archivo explica como subir y correr Pipeline_scanpy_SCENIC_Lupus_Colab.ipynb
en Google Colab.

------------------------------------------------------------
1. SUBIR EL NOTEBOOK A COLAB
------------------------------------------------------------

Se entra a https://colab.research.google.com con una cuenta de Google.

En la pantalla de inicio aparece un cuadro para abrir notebooks. Se elige la
pestana "Subir" y se selecciona el archivo Pipeline_scanpy_SCENIC_Lupus_Colab.ipynb
desde la computadora. Tambien se puede arrastrar el archivo directamente a esa
ventana.

El notebook se abre en una pestana nueva, listo para ejecutarse.

------------------------------------------------------------
2. ELEGIR EL MODO
------------------------------------------------------------

En la celda 1, la unica que hace falta tocar, hay una linea:

    MODO = "ejemplo"

Se cambia a "zenodo" para correr sobre el estudio real de lupus, o se deja en
"ejemplo" para la demostracion con medula osea.

------------------------------------------------------------
3. EJECUTAR TODO
------------------------------------------------------------

En el menu de arriba: Entorno de ejecucion > Ejecutar todas.

Si Colab pide reiniciar el entorno a la mitad, es normal, algunas librerias lo
requieren. Se acepta y se vuelve a ejecutar todas las celdas.

------------------------------------------------------------
4. QUE ESPERAR SEGUN EL MODO
------------------------------------------------------------

Modo ejemplo:
  - Corre en unos 15 minutos en el Colab gratuito.
  - No produce biologia real, el dataset es demasiado chico para que SCENIC
    valide regulones, pero muestra el flujo completo de principio a fin,
    incluyendo SCENIC hasta el final, marcado como no validado.

Modo zenodo:
  - Usa el estudio real de lupus: 169.513 celulas por defecto.
  - Esto necesita RAM alta (Colab Pro o el entorno de RAM alta activado). La
    celda 2.1 revisa la RAM disponible y se detiene con una explicacion si no
    alcanza.
  - Si solo se tiene el Colab gratuito, se puede poner un numero menor en
    N_CELLS_MAX (celda 1) para trabajar con una muestra en vez del dataset
    completo.

------------------------------------------------------------
5. COMO SABER SI EL RESULTADO SIRVE PARA ENTREGAR
------------------------------------------------------------

Al final del notebook, la celda de procedencia (seccion 6.7) imprime si la
corrida es apta para entrega o no. Solo lo es en modo zenodo, con GRNBoost2
como algoritmo de red y regulones validados por cisTarget.

Si algun archivo de salida lleva el sufijo _SIN_VALIDAR en el nombre, no es un
resultado de SCENIC. Es una demostracion del mecanismo, no biologia real.

------------------------------------------------------------
6. ARCHIVOS QUE QUEDAN GUARDADOS
------------------------------------------------------------

Si se conecta Google Drive al principio (celda 2.4), los resultados quedan en
Mi unidad/scanpy_scenic_lupus/. Si no, se pierden al cerrar la sesion de Colab,
asi que conviene descargarlos antes de cerrar.

Archivos principales:
  - DE_completo_<modo>.csv y top100_DE_<modo>.csv: genes diferenciales por
    grupo celular.
  - scenic_auc_<modo>.csv: actividad de los regulones por celula.
  - lupus_rituximab.h5ad: los datos convertidos (solo en modo zenodo).

------------------------------------------------------------
7. SI ALGO FALLA
------------------------------------------------------------

El notebook esta escrito para detenerse con un mensaje que explica la causa,
en vez de seguir adelante en silencio. Conviene leer el error completo antes
de reintentar: casi siempre dice que revisar o que parametro ajustar.

============================================================
