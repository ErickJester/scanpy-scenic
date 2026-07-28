# -*- coding: utf-8 -*-
"""
Pipeline scRNA-seq: scanpy + SCENIC. Version de escritorio.

Hace lo mismo que el notebook de Colab, en un solo archivo y sin depender de
Google Colab: control de calidad, agrupamiento de celulas, las 12 graficas,
las dos figuras del estudio y el analisis de regulacion con SCENIC.

Uso basico (todo se configura en el PANEL DE CONTROL de mas abajo):

    python pipeline_scanpy_scenic.py

Tambien se puede cambiar la configuracion desde la linea de comandos sin editar
el archivo:

    python pipeline_scanpy_scenic.py --modo ejemplo
    python pipeline_scanpy_scenic.py --n-cells-max 30000
    python pipeline_scanpy_scenic.py --salida D:\\resultados_lupus

Las graficas se guardan como PNG en la carpeta de resultados, porque aqui no hay
un cuaderno donde mostrarlas. El resto del analisis es identico.

Requisitos: las librerias de Python que instala instalar_dependencias.bat, y
ademas R con el paquete SeuratObject si se usa el modo "zenodo". El script
comprueba las dos cosas antes de empezar y explica que falta.
"""

# ============================================================
#   PANEL DE CONTROL: todo lo configurable del analisis
# ============================================================
# Es lo unico que hay que tocar para adaptar el analisis. Cada valor se puede
# sobreescribir desde la linea de comandos.

# Que datos analizar.
#   "ejemplo" = medula osea publica, para ver el flujo funcionando en minutos.
#   "zenodo"  = el dataset real del estudio, que se descarga de Zenodo.
MODO = "zenodo"

# ---- Modo rapido -------------------------------------------
# En True se saltan las 12 graficas y las dos figuras del estudio, junto con
# los calculos que solo sirven para dibujarlas (PCA, vecinos, UMAP, clustering,
# expresion diferencial y trayectoria). Sirve para comprobar en poco tiempo que
# SCENIC llega hasta el final. Una corrida asi no es entregable: le faltan las
# graficas y las tablas de expresion diferencial.
SALTAR_GRAFICAS = False

# ---- Cuantas celulas analizar (solo MODO="zenodo") ---------
# El dataset trae 169.513 celulas. Con 0 se procesan todas, lo que pide bastante
# memoria: contando con el margen que necesita R para abrir el archivo, conviene
# tener 32 GB. Con un numero menor se toma una muestra estratificada por tipo
# celular, que conserva la proporcion de cada uno. El muestreo solo estratifica
# por esa columna: los grupos que forma la figura 2, que cruza 9 pacientes con
# varios momentos, quedan mas cortos cuanto mas se recorte.
N_CELLS_MAX = 60000

# ---- Modo depuracion (solo MODO="zenodo") ------------------
# Leer el archivo del estudio (1,6 GB) tarda varios minutos cada vez. En True se
# prepara una vez una copia pequena y el resto del analisis trabaja sobre ella,
# lo que permite probar cambios en segundos. La copia sale de recortar el
# archivo original, asi que conserva su estructura real.
# Los resultados de una corrida en este modo son solo para probar, nunca para
# entregar.
RDS_DEBUG = False
N_CELLS_DEBUG = 4000
RDS_DEBUG_FILE = "RTX_zenodo_debug.RDS"

# ---- SCENIC: tamano del analisis ---------------------------
# El primer paso de SCENIC entrena un modelo por gen, y su costo crece rapido con
# el numero de celulas. Con SCENIC_DOWNSAMPLE en True se recorta la entrada a
# SCENIC_N_CELLS antes de empezar. El recorte solo afecta a SCENIC: las graficas
# y las figuras del estudio se calculan sobre el conjunto completo.
SCENIC_DOWNSAMPLE = True
SCENIC_N_CELLS = 2000
SCENIC_N_GENES = 1500

# Algoritmo con el que se reconstruye la red.
#   "grnboost2"     = el metodo oficial de SCENIC
#   "sklearn_aprox" = aproximacion, solo para demostrar el flujo
METODO_GRN = "grnboost2"
# La red no se recorta a mano en ningun momento. SCENIC ya aplica sus propios
# umbrales al formar los grupos de genes, y filtrar antes se apartaria del
# procedimiento publicado.

# ---- Exigir que los regulones esten validados --------------
# SCENIC corre en tres pasos encadenados y este valor decide que pasa si el
# segundo no deja nada: reconstruir la red, filtrarla, y puntuar el resultado
# celula a celula. Es el filtro intermedio el que puede devolver una lista vacia.
# En False el analisis se detiene ahi, que es lo apropiado para un entregable.
# En True sigue al tercer paso y marca todas las salidas como no validadas.
# Con MODO="ejemplo" siempre sigue: ese dataset es demasiado pequeno para pasar
# el filtro y existe solo para recorrer el flujo entero.
PERMITIR_REGULONES_SIN_VALIDAR = False

# ---- Columnas del estudio (None = detectar sola) -----------
# El script busca por su cuenta que columna guarda cada dato. Si no la encuentra,
# muestra las columnas disponibles y se detiene para que se indique aqui cual
# usar.
COL_COND = None     # diagnostico (lupus o control sano)
COL_TIME = None     # momento del tratamiento
# El dataset trae dos columnas de tipo celular con distinto nivel de detalle. La
# primera agrupa todas las celulas B bajo una sola etiqueta, con lo que la figura
# 1 saldria con una unica categoria. La segunda las desglosa en siete, que es lo
# que esa figura necesita.
COL_CTYPE = "Celltype_level2"

# ---- Los dos grupos que compara la figura 2 ----------------
# Se derivan del nombre de cada muestra. Se pueden fijar aqui para comparar otro
# par de los disponibles.
PRE_LABEL = None
POST_LABEL = None

# ---- Integridad de las descargas ---------------------------
# Comprueba que las bases de datos de SCENIC llegaron completas comparando su
# huella digital. Una descarga cortada produciria resultados silenciosamente
# equivocados.
VERIFICAR_CHECKSUMS = True

# ---- Trayectoria celular (grafica 12) ----------------------
# Grupo desde el que arranca el ordenamiento. Se escribe como texto, por ejemplo
# "0". Con None lo elige el propio script.
ROOT_CLUSTER = None

# ---- Agrupamiento de referencia ----------------------------
# Los grupos de celulas que usan la anotacion, la expresion diferencial y la
# trayectoria. El numero es la resolucion: mas alto da mas grupos y mas finos.
CLUSTER_KEY = "leiden_res_0.50"

# ---- Carpetas ----------------------------------------------
# Donde quedan las tablas y las graficas, y donde se guardan las descargas. Se
# separan porque las descargas pesan varios gigabytes y conviene no repetirlas
# entre corridas.
DIR_SALIDA = "resultados"
DIR_DATOS = "datos"
# ============================================================

import argparse
import datetime
import gc
import hashlib
import importlib
import os
import re
import shutil
import site
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
import warnings


# ------------------------------------------------------------------
#   Utilidades generales
# ------------------------------------------------------------------

def _fmt_duracion(segundos):
    m, s = divmod(int(segundos), 60)
    h, m = divmod(m, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


class Heartbeat:
    """Avisa periodicamente que el proceso sigue vivo y cuanto CPU esta usando.

    No hay forma honesta de dar un porcentaje de avance en pasos como la
    conversion en R o el calculo de vecinos de scanpy: esas librerias no
    exponen cuanto les falta. Lo que si se puede dar, y es lo que realmente
    hace falta para distinguir "esta trabajando" de "se colgo", es cuanto
    tiempo lleva cada paso y si el CPU esta activo. Un proceso colgado se
    queda en 0% de CPU indefinidamente; uno que solo tarda sigue consumiendo.
    """

    def __init__(self, intervalo=25):
        self.intervalo = intervalo
        self.paso_actual = "arrancando"
        self.duracion_esperada = None
        self._t_paso = time.time()
        self._t_inicio = time.time()
        self._detener = threading.Event()
        self._hilo = None
        self._proc = None
        try:
            import psutil
            self._proc = psutil.Process()
            self._proc.cpu_percent(interval=None)   # primer llamado, se descarta
        except Exception:
            pass

    def marcar(self, nombre, duracion_esperada=None):
        """Cambia el paso actual y reinicia su cronometro."""
        self.paso_actual = nombre
        self.duracion_esperada = duracion_esperada
        self._t_paso = time.time()

    def _cpu_total(self):
        """Suma el uso de CPU del proceso y de todo lo que haya lanzado (R, dask...)."""
        if self._proc is None:
            return None
        try:
            total = self._proc.cpu_percent(interval=None)
            for hijo in self._proc.children(recursive=True):
                try:
                    total += hijo.cpu_percent(interval=None)
                except Exception:
                    pass
            return total
        except Exception:
            return None

    def _bucle(self):
        while not self._detener.wait(self.intervalo):
            transcurrido_total = time.time() - self._t_inicio
            transcurrido_paso = time.time() - self._t_paso
            cpu = self._cpu_total()
            partes = [f"{_fmt_duracion(transcurrido_total)} transcurridos"]
            partes.append(f"paso actual: {self.paso_actual} ({_fmt_duracion(transcurrido_paso)})")
            if self.duracion_esperada:
                partes.append(f"suele tardar {self.duracion_esperada}")
            if cpu is not None:
                partes.append(f"CPU: {cpu:.0f}%")
                if cpu < 3 and transcurrido_paso > 90:
                    partes.append("<- CPU casi en cero durante mas de un minuto: podria estar colgado")
            print("  [latido] " + " | ".join(partes), flush=True)

    def iniciar(self):
        self._hilo = threading.Thread(target=self._bucle, daemon=True)
        self._hilo.start()

    def detener(self):
        self._detener.set()
        if self._hilo is not None:
            self._hilo.join(timeout=2)


_HEARTBEAT = None    # la activa main(); las funciones de mas abajo la consultan si existe


def titulo(texto, duracion_esperada=None):
    """Separador visible para seguir el avance en una consola larga."""
    print()
    print("=" * 70)
    print("  " + texto)
    print("=" * 70)
    if _HEARTBEAT is not None:
        _HEARTBEAT.marcar(texto, duracion_esperada)


def paso(texto, duracion_esperada=None):
    print()
    print("--- " + texto + " " + "-" * max(0, 66 - len(texto)))
    if duracion_esperada:
        print(f"    (suele tardar {duracion_esperada}; el aviso de mas abajo confirma que sigue vivo)")
    if _HEARTBEAT is not None:
        _HEARTBEAT.marcar(texto, duracion_esperada)


def ejecutar_streaming(cmd, env=None, prefijo="  "):
    """Corre un programa externo mostrando su salida en vivo, linea por linea.

    Los pasos que llaman a R son los que mas tardan y los que menos avisan de
    su progreso. Con subprocess.run() normal, toda su salida queda en un
    buffer y no se ve nada hasta que termina (o hasta que se cuelga): eso es
    justo lo que paso con el filtrado de SCENIC, que se quedo callado 10 horas
    sin que nada distinguiera "trabajando" de "colgado". Mostrando cada linea
    en cuanto R la escribe, un silencio prolongado es una senal real, no una
    posibilidad entre varias.
    Devuelve un objeto con .returncode, .stdout y .stderr, para que el resto
    del codigo que ya revisa esos campos no tenga que cambiar.
    """
    proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, encoding="utf-8", errors="replace", bufsize=1)
    salida_out, salida_err = [], []

    def _leer(flujo, buffer, marca):
        for linea in iter(flujo.readline, ""):
            buffer.append(linea)
            print(f"{prefijo}{marca}{linea.rstrip()}", flush=True)
        flujo.close()

    h_out = threading.Thread(target=_leer, args=(proc.stdout, salida_out, ""))
    h_err = threading.Thread(target=_leer, args=(proc.stderr, salida_err, "! "))
    h_out.start()
    h_err.start()
    proc.wait()
    h_out.join()
    h_err.join()

    class _Resultado:
        pass

    r = _Resultado()
    r.returncode = proc.returncode
    r.stdout = "".join(salida_out)
    r.stderr = "".join(salida_err)
    return r


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blq in iter(lambda: fh.read(chunk), b""):
            h.update(blq)
    return h.hexdigest()


def parece_html(path):
    """Un error 404/403 suele guardarse como pagina HTML con nombre de archivo."""
    with open(path, "rb") as fh:
        inicio = fh.read(512).lstrip().lower()
    return inicio.startswith(b"<!doctype html") or inicio.startswith(b"<html")


def descargar_verificado(nombre, url, destino, esperado, verificar=True):
    """Descarga comprobando que el archivo llego completo y sin alterar."""
    # Un archivo que ya existe puede venir cortado de un intento anterior, asi que
    # se comprueba igual antes de darlo por bueno.
    if os.path.exists(destino):
        if not verificar or esperado is None:
            print(f"  {nombre}: {os.path.getsize(destino)/1e6:.1f} MB (ya estaba)")
            return
        obtenido = sha256(destino)
        if obtenido == esperado:
            print(f"  {nombre}: {os.path.getsize(destino)/1e6:.1f} MB  SHA256 OK")
            return
        print(f"  {nombre}: el archivo local NO coincide con el hash fijado -> se rebaja")
        os.remove(destino)

    # Se descarga con un nombre provisional y se renombra al terminar. Asi una
    # descarga interrumpida nunca queda con el nombre definitivo, haciendose
    # pasar por completa.
    parcial = destino + ".part"
    print(f"  Descargando {nombre} ...")
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            if resp.status != 200:
                raise RuntimeError(f"HTTP {resp.status} al descargar {nombre}")
            declarado = resp.headers.get("Content-Length")
            declarado = int(declarado) if declarado else None
            with open(parcial, "wb") as fh:
                while True:
                    trozo = resp.read(1 << 20)
                    if not trozo:
                        break
                    fh.write(trozo)
    except urllib.error.URLError as e:
        if os.path.exists(parcial):
            os.remove(parcial)
        raise RuntimeError(f"No se pudo descargar {nombre} desde {url}\n  {e}") from e

    real = os.path.getsize(parcial)
    if declarado is not None and real != declarado:
        os.remove(parcial)
        raise RuntimeError(
            f"{nombre}: descarga truncada ({real:,} de {declarado:,} bytes).\n"
            "Vuelve a ejecutar el script; suele ser un corte de red temporal."
        )
    if parece_html(parcial):
        os.remove(parcial)
        raise RuntimeError(
            f"{nombre}: el servidor devolvio una pagina HTML, no el archivo.\n"
            f"URL probablemente caida o movida: {url}"
        )

    obtenido = sha256(parcial)
    if esperado is None:
        print(f"    (sin hash fijado) SHA256 = {obtenido}")
    elif verificar and obtenido != esperado:
        os.remove(parcial)
        raise RuntimeError(
            f"{nombre}: SHA256 no coincide.\n"
            f"  esperado: {esperado}\n  obtenido: {obtenido}\n"
            "O la descarga se corrompio (re-ejecuta), o el archivo de origen cambio.\n"
            "Si el cambio es esperado, actualiza CHECKSUMS con el hash nuevo."
        )

    os.replace(parcial, destino)
    print(f"  {nombre}: {real/1e6:.1f} MB  SHA256 {'OK' if esperado else 'registrado'}")


def buscar_rscript():
    """Devuelve la ruta de Rscript, o None.

    En Windows el instalador de R no siempre lo deja en el PATH, asi que ademas
    de preguntar al sistema se miran las rutas donde suele quedar.
    """
    encontrado = shutil.which("Rscript")
    if encontrado:
        return encontrado
    import glob
    patrones = [
        r"C:\Program Files\R\R-*\bin\x64\Rscript.exe",
        r"C:\Program Files\R\R-*\bin\Rscript.exe",
        r"C:\Program Files (x86)\R\R-*\bin\Rscript.exe",
        os.path.expanduser(r"~\AppData\Local\Programs\R\R-*\bin\x64\Rscript.exe"),
    ]
    candidatos = []
    for p in patrones:
        candidatos.extend(glob.glob(p))
    if candidatos:
        return sorted(candidatos)[-1]      # la version mas reciente
    return None


# ------------------------------------------------------------------
#   Compatibilidad de SCENIC con numpy moderno
# ------------------------------------------------------------------
# SCENIC se publico cuando numpy usaba unos nombres que la version actual ya
# retiro. Aqui se reponen esos nombres para que la libreria funcione sin
# modificar ninguno de sus archivos ni instalar una version antigua de numpy, que
# romperia el resto del entorno.

_MODULO_ALIAS = '''"""Repone en numpy los alias que las versiones 1.24 y 2.0 retiraron.

Lo importa un archivo .pth, asi que corre al arrancar cualquier interprete de
este entorno: por eso no imprime nada ni deja escapar avisos.
No modifica archivos de numpy ni de pySCENIC.
"""
import warnings

try:
    import numpy
except ImportError:
    numpy = None

_TIPOS_124 = (("object", object), ("bool", bool), ("int", int),
              ("float", float), ("str", str), ("complex", complex))

_EXTRA = {"unicode_": "str_", "string_": "bytes_", "round_": "round",
          "product": "prod", "cumproduct": "cumprod", "sometrue": "any",
          "alltrue": "all", "infty": "inf", "Inf": "inf", "Infinity": "inf",
          "NAN": "nan", "NaN": "nan", "float_": "float64",
          "complex_": "complex128", "mat": "asmatrix", "in1d": "isin",
          "row_stack": "vstack", "trapz": "trapezoid"}


def _falta(nombre):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return not hasattr(numpy, nombre)


def _reemplazo(motivo):
    i = motivo.find("`np.")
    if i == -1:
        return None
    j = motivo.find("`", i + 1)
    return motivo[i + 4:j] if j != -1 else None


def restaurar():
    if numpy is None:
        return []
    puestos = []
    for viejo, tipo in _TIPOS_124:
        if _falta(viejo):
            setattr(numpy, viejo, tipo)
            puestos.append(viejo)
    for viejo, nuevo in _EXTRA.items():
        if _falta(viejo) and hasattr(numpy, nuevo):
            setattr(numpy, viejo, getattr(numpy, nuevo))
            puestos.append(viejo)
    for viejo, motivo in getattr(numpy, "__expired_attributes__", {}).items():
        nuevo = _reemplazo(motivo)
        if nuevo and _falta(viejo) and hasattr(numpy, nuevo):
            setattr(numpy, viejo, getattr(numpy, nuevo))
            puestos.append(viejo)
    return puestos


restaurar()
'''


def preparar_alias_numpy():
    """Repone los nombres retirados de numpy, aqui y en los procesos que se lancen.

    Hace falta en dos sitios: en esta sesion y en el programa aparte que SCENIC
    lanza mas adelante, que arranca de cero y no hereda nada de lo que hay aqui.
    Se guarda como archivo para que las dos partes usen exactamente lo mismo.
    """
    import numpy as np

    # El archivo tiene que quedar en la misma carpeta que numpy. Python recorre
    # esas carpetas en orden y ejecuta lo que encuentra en cada una nada mas
    # anadirla, de modo que colocarlo en otra haria que se ejecutase cuando numpy
    # todavia no esta disponible y no serviria de nada.
    carpeta_numpy = os.path.dirname(os.path.dirname(os.path.abspath(np.__file__)))
    candidatas = [carpeta_numpy]
    for extra in (site.getsitepackages() if hasattr(site, "getsitepackages") else []):
        if extra not in candidatas:
            candidatas.append(extra)

    nombre_mod = "_alias_numpy_scenic"
    destino = None
    for carpeta in candidatas:
        try:
            with open(os.path.join(carpeta, nombre_mod + ".py"), "w", encoding="utf-8") as fh:
                fh.write(_MODULO_ALIAS)
            with open(os.path.join(carpeta, "zz_" + nombre_mod + ".pth"), "w", encoding="utf-8") as fh:
                fh.write("import " + nombre_mod + "\n")
            destino = carpeta
            break
        except OSError:
            continue

    if destino is None:
        print("AVISO: no se pudo escribir el archivo de compatibilidad de numpy.")
        print("       Si SCENIC falla mas adelante con 'np.object', reinstala en un")
        print("       entorno virtual propio (ver instalar_dependencias.bat).")
    else:
        print(f"  compatibilidad de numpy instalada en: {destino}")

    # En esta sesion se aplica el mismo archivo que usara el programa aparte, para
    # que los dos se comporten igual.
    sys.path.insert(0, destino or ".")
    try:
        mod = importlib.import_module(nombre_mod)
        puestos = mod.restaurar()
        print(f"  alias repuestos en esta sesion: {len(puestos)}")
    except Exception as e:
        print(f"  AVISO: no se pudieron reponer los alias aqui ({type(e).__name__}: {e})")

    # Se lanza un programa de prueba para confirmar que ahi tambien estan los
    # nombres, en vez de darlo por hecho. Se comprueba uno de cada tanda de
    # retiradas, porque llegaron en versiones distintas de numpy.
    sonda = ("import numpy, warnings; warnings.simplefilter('ignore'); "
             "print(hasattr(numpy,'object') and hasattr(numpy,'unicode_'))")
    try:
        r = subprocess.run([sys.executable, "-c", sonda], capture_output=True, text=True, timeout=120)
        print(f"  subproceso con alias: {r.stdout.strip() or 'sin respuesta'}")
    except Exception as e:
        print(f"  no se pudo comprobar el subproceso ({type(e).__name__})")


# ------------------------------------------------------------------
#   Graficas
# ------------------------------------------------------------------

class Figuras:
    """Guarda cada grafica como PNG conservando la numeracion del cuaderno.

    En el cuaderno cada grafica se dibuja debajo de la celda que la genera, y
    varias de esas celdas producen dos imagenes: la del PCA saca los componentes
    y ademas la varianza explicada, la de anotacion saca los grupos y ademas las
    etiquetas, y asi. Al escribirlas a disco cada imagen es un archivo aparte, de
    modo que numerarlas de corrido las desalinearia de la numeracion del
    cuaderno a partir de la primera que se duplica.
    Por eso el nombre completo, numero incluido, lo decide quien llama, con una
    letra cuando una misma grafica produce varias imagenes: 04a y 04b son las dos
    mitades de la grafica 4.
    """

    def __init__(self, carpeta):
        self.carpeta = carpeta
        self.n = 0
        os.makedirs(carpeta, exist_ok=True)

    def guardar(self, nombre):
        import matplotlib.pyplot as plt
        self.n += 1
        archivo = os.path.join(self.carpeta, f"{nombre}.png")
        plt.savefig(archivo, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"  [figura] {os.path.basename(archivo)}")
        return archivo


# ------------------------------------------------------------------
#   Seleccion de celulas B
# ------------------------------------------------------------------
# Que cuenta como celula B se decide en esta sola funcion, que usan tanto la
# figura de subtipos como el analisis de SCENIC. Teniendo el criterio en un unico
# sitio, las dos partes trabajan siempre sobre el mismo conjunto de celulas.

B_KW = ["naive b", "transitional", "memory b", "switched", "abc",
        "plasmablast", "plasma", "b cell", "b_cell"]

# Etiquetas que designan celulas B sin nombrar ningun subtipo, propias de las
# anotaciones poco detalladas. Se comparan como texto completo y no como
# fragmento: buscar una "b" suelta dentro de cada etiqueta daria por celulas B a
# casi todas las demas.
B_EXACTAS = {"b", "b cells", "bcell", "bcells", "b-cell", "b-cells", "b lymphocyte"}

# Nombres habituales de las columnas que guardan el tipo celular. Sirven para
# proponer una alternativa cuando la columna elegida no encuentra celulas B.
_CLAVES_CTYPE = ("celltype", "cell_type", "cell.type", "annotation", "ident")


def _marcar_B(etiquetas):
    """Mascara booleana de celulas B para una serie de etiquetas."""
    bajas = etiquetas.str.lower().str.strip()
    # Las palabras se tratan como texto literal, de modo que anadir a la lista una
    # etiqueta con puntos o signos no altere la busqueda.
    patron = "|".join(re.escape(k) for k in B_KW)
    return bajas.str.contains(patron, regex=True) | bajas.isin(B_EXACTAS)


def seleccionar_celulas_B(adata, col_ctype, min_celulas=50):
    """Devuelve el subconjunto de celulas B segun la columna de tipo celular.

    Las palabras que busca son las del dataset actual. Con otra nomenclatura el
    filtro no encuentra nada, y en ese caso la funcion para y muestra las
    etiquetas reales en vez de devolver un objeto vacio.
    """
    etiquetas = adata.obs[col_ctype].astype(str)
    es_B = _marcar_B(etiquetas)
    n = int(es_B.sum())

    if n < min_celulas:
        # Cuando no aparecen celulas B, lo habitual es que la columna consultada
        # no sea la adecuada, no que falten palabras en la lista. Los estudios
        # suelen traer varios niveles de anotacion, y la busqueda automatica puede
        # haberse quedado con el menos detallado. Antes de rendirse se prueban las
        # demas columnas y se indica cual funciona.
        alternativas = []
        for c in adata.obs.columns:
            if c == col_ctype or not any(k in c.lower() for k in _CLAVES_CTYPE):
                continue
            try:
                m = int(_marcar_B(adata.obs[c].astype(str)).sum())
            except Exception:
                continue
            if m >= min_celulas:
                alternativas.append((c, m))

        if alternativas:
            sugerencia = (
                "\n\nLa columna elegida no parece ser la correcta. Estas SI encuentran celulas B:\n"
                + "\n".join(f"    COL_CTYPE = {c!r}   ->  {m:,} celulas B" for c, m in alternativas)
                + "\n\nPonla en el PANEL DE CONTROL y vuelve a ejecutar.\n"
                  "No hace falta repetir la conversion de R: el archivo convertido se reutiliza."
            )
        else:
            sugerencia = (
                "\n\nNinguna otra columna de anotacion encuentra celulas B tampoco, asi que\n"
                "el problema si esta en la nomenclatura: ajusta B_KW en este archivo."
            )

        raise ValueError(
            f"El filtro de celulas B encontro {n} celulas (minimo esperado: {min_celulas}).\n"
            f"Se buscaron estas palabras clave en {col_ctype!r}: {B_KW}\n"
            f"Etiquetas reales en {col_ctype!r}:\n"
            + "\n".join(f"    {e!r}: {c:,}" for e, c in etiquetas.value_counts().items())
            + sugerencia
        )

    sub = adata[es_B].copy()
    print(f"Celulas B: {sub.n_obs:,} de {adata.n_obs:,} ({100*sub.n_obs/adata.n_obs:.1f}%)")

    # Encontrar celulas B no garantiza que la anotacion distinga subtipos. Con una
    # etiqueta unica la seleccion es correcta, pero la figura sale con una sola
    # categoria y deja de mostrar lo que su titulo anuncia.
    n_subtipos = sub.obs[col_ctype].astype(str).nunique()
    if n_subtipos < 2:
        print(f"AVISO: {col_ctype!r} solo distingue {n_subtipos} etiqueta(s) dentro de las\n"
              "       celulas B, asi que la figura no mostrara subtipos reales. Si el dataset\n"
              "       tiene una anotacion mas fina (p.ej. un 'level2'), ponla en COL_CTYPE.")
    return sub


# ------------------------------------------------------------------
#   El script de R que traduce el archivo del estudio
# ------------------------------------------------------------------
# El estudio se publica en el formato de Seurat, que es una herramienta de R,
# mientras que este analisis corre en Python con scanpy. Este script traduce de
# uno a otro: exporta la matriz de expresion, los metadatos de cada celula y las
# coordenadas del UMAP a formatos que Python lee.
# Todo lo que el script da por supuesto de la estructura del archivo se comprueba
# antes de usarlo, de modo que un archivo organizado de otra forma se detenga con
# un mensaje que diga que falta, en vez de exportar datos mal.

R_CONVERTIR = r'''
# R solo mira su carpeta personal de paquetes si el directorio ya existe. Se
# anade a mano para no depender de eso: ahi es donde quedaron instalados.
local({
  lib <- Sys.getenv("R_LIBS_USER")
  if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))
})
suppressMessages({library(SeuratObject); library(Matrix)})
n_max <- as.integer(Sys.getenv("N_CELLS_MAX", "60000"))

# El archivo a convertir se recibe desde fuera, para poder apuntar a la copia de
# pruebas sin tocar nada de este script.
fuente <- Sys.getenv("RDS_FILE", "RTX_zenodo.RDS")
salida <- Sys.getenv("DIR_TRABAJO", ".")
cat("Leyendo:", fuente, "\n")
obj <- readRDS(fuente)
if (!inherits(obj, "Seurat"))
  stop("El RDS no contiene un objeto Seurat, sino: ", paste(class(obj), collapse="/"))

cat("=== ESTRUCTURA ===\n"); print(obj)
cat("\n=== COLUMNAS METADATA ===\n"); print(colnames(obj@meta.data))
meta <- obj@meta.data
cat("\n=== GRUPOS ===\n")
for (col in colnames(meta)) {
  v <- meta[[col]]
  if (is.factor(v) || is.character(v) || (is.numeric(v) && length(unique(v)) < 30))
    if (length(unique(v)) <= 30) { cat("\n[", col, "]\n", sep=""); print(table(v)) }
}

# --- La matriz de expresion --------------------------------------------------
# Se lee antes de decidir la muestra, porque el orden de sus columnas es el orden
# real de las celulas y es contra el que se alinea todo lo demas.
ar <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cat("Assay usado:", ar, "de", paste(Assays(obj), collapse=", "), "\n")

# El archivo puede traer varias matrices bajo distintos nombres, y el nombre no
# garantiza el contenido: en este dataset la que se llama de conteos ya viene
# transformada.
# Por eso la eleccion se hace mirando los valores y no el nombre. Los conteos sin
# procesar son enteros; cualquier version transformada tiene decimales. Si
# ninguna resulta de enteros se toma la primera y se avisa, y mas adelante el
# analisis mide en que estado llega la matriz y aplica solo lo que falte.
#
# La lectura va directa a la estructura interna del archivo en lugar de usar las
# funciones habituales, y solo trae un bloque de columnas de cada matriz. Esas
# funciones cargan la matriz entera, y con 169.513 columnas eso agota la memoria
# de la sesion antes de poder decidir nada.
assay_obj <- obj@assays[[ar]]

capas <- tryCatch(names(assay_obj@layers), error=function(e) NULL)
if (is.null(capas) || length(capas) == 0) capas <- c("counts", "data")
cat("Capas del assay:", paste(capas, collapse=", "), "\n")

leer_bloque <- function(nombre, n_celulas=2000L) {
  m <- tryCatch(assay_obj@layers[[nombre]], error=function(e) NULL)
  if (is.null(m))
    m <- tryCatch(slot(assay_obj, nombre), error=function(e) NULL)
  if (is.null(m) || is.null(dim(m))) return(NULL)
  m[, seq_len(min(n_celulas, ncol(m))), drop=FALSE]
}
es_entera <- function(bloque) {
  if (is.null(bloque) || nrow(bloque) == 0 || ncol(bloque) == 0) return(NA)
  v <- if (inherits(bloque, "sparseMatrix")) bloque@x else as.numeric(bloque)
  if (length(v) == 0) return(NA)
  v <- v[seq_len(min(200000L, length(v)))]
  all(abs(v - round(v)) < 1e-8)
}

capa_usada <- NA_character_
for (cp in capas) {
  bloque <- leer_bloque(cp)
  ent <- es_entera(bloque)
  cat("  capa '", cp, "': ",
      if (is.na(ent)) "no legible"
      else if (ent) "valores ENTEROS (conteos crudos)"
      else "valores decimales (ya transformada)", "\n", sep="")
  rm(bloque); gc()
  if (!is.na(ent) && ent) { capa_usada <- cp; break }
}

if (is.na(capa_usada)) {
  capa_usada <- "counts"
  cat("AVISO: ninguna capa del assay trae conteos crudos (enteros).\n",
      "       Se usa 'counts' tal cual viene, ya normalizada en origen.\n",
      "       El analisis lo detecta y NO vuelve a normalizar.\n", sep="")
} else {
  cat("Capa elegida para la matriz: ", capa_usada, "\n", sep="")
}

counts <- tryCatch(assay_obj@layers[[capa_usada]], error=function(e) NULL)
if (is.null(counts))
  counts <- tryCatch(slot(assay_obj, capa_usada), error=function(e) NULL)
if (is.null(counts))
  counts <- tryCatch(GetAssayData(obj, assay=ar, slot=capa_usada),
                     error=function(e) tryCatch(GetAssayData(obj, assay=ar, layer=capa_usada),
                                                error=function(e2) NULL))
rm(assay_obj); gc()
if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0)
  stop("No se pudo leer ninguna matriz de expresion del assay '", ar, "'.")

# --- Emparejar cada celula con su anotacion ----------------------------------
# La matriz de expresion y la tabla de metadatos son dos listas independientes, y
# nada garantiza que traigan las celulas en el mismo orden. Si se emparejan por
# posicion cuando los ordenes difieren, cada fila de metadatos acaba pegada a la
# celula equivocada y todas sus columnas quedan mal asignadas.
# El emparejamiento va por nombre de celula, y se verifica. Es un punto critico
# porque un cruce aqui no produce ningun error visible: el analisis termina bien
# y entrega resultados que parecen validos sin serlo.
if (is.null(colnames(counts)))
  stop("La matriz de conteos no trae nombres de columna (celulas).")
if (ncol(counts) != nrow(meta))
  stop("La matriz tiene ", ncol(counts), " celulas y los metadatos ", nrow(meta), " filas.")

# En este dataset la tabla de metadatos esta numerada del 1 en adelante, sin
# nombres de celula: el identificador real vive en una columna aparte. Se prueba
# primero la numeracion de filas, por si otro dataset si la trae util, y se pasa
# a esa columna cuando no cubre las celulas de la matriz.
clave <- rownames(meta)
origen_clave <- "rownames(meta.data)"
cubre <- !is.null(clave) && !anyDuplicated(clave) && all(colnames(counts) %in% clave)

if (!cubre && "Barcode" %in% colnames(meta)) {
  bc <- as.character(meta$Barcode)
  if (!anyDuplicated(bc) && all(colnames(counts) %in% bc)) {
    clave <- bc; origen_clave <- "la columna 'Barcode'"; cubre <- TRUE
    cat("AVISO: rownames(meta.data) no sirven para alinear (no son identificadores\n",
        "       de celula validos, p.ej.: ", paste(head(rownames(meta), 3), collapse=", "), ").\n",
        "       Se usa la columna 'Barcode' como clave real de celula.\n", sep="")
  }
}

if (!cubre) {
  ejemplo_rn <- if (is.null(rownames(meta))) character(0) else head(rownames(meta), 3)
  stop("No se encontro una clave que alinee los metadatos con la matriz.\n",
       "  rownames(meta.data), ejemplo: ", paste(ejemplo_rn, collapse=", "), "\n",
       "  colnames(counts), ejemplo: ", paste(head(colnames(counts), 3), collapse=", "), "\n",
       "Ni rownames(meta.data) ni la columna 'Barcode' (si existe) cubren todas\n",
       "las celulas de la matriz. Revisa a mano que columna trae el barcode real.")
}

rownames(meta) <- clave
if (identical(colnames(counts), rownames(meta))) {
  cat("Alineacion matriz/metadatos: OK (mismo orden, clave=", origen_clave, ")\n", sep="")
} else {
  cat("Los metadatos se reordenan por nombre de celula (clave=", origen_clave,
      ") para que cada celula lleve su propia anotacion.\n", sep="")
  meta <- meta[colnames(counts), , drop=FALSE]
}

ncells <- ncol(counts)
if (n_max > 0 && n_max < ncells) {
  set.seed(0)
  cc <- grep("celltype|cell_type|cell.type|annotation|ident", colnames(meta),
             ignore.case=TRUE, value=TRUE)
  if (length(cc) > 0) {
    # La muestra se toma dentro de cada tipo celular, de modo que conserve las
    # proporciones del dataset completo.
    grp <- as.character(meta[[cc[1]]]); fr <- n_max / ncells
    idx <- sort(unlist(lapply(split(seq_len(ncells), grp),
                              function(ix) sample(ix, max(1, round(length(ix) * fr))))))
    cat("\n>>> DOWNSAMPLE estratificado por '", cc[1], "': ", ncells, " -> ", length(idx), "\n", sep="")
  } else {
    idx <- sort(sample(ncells, n_max))
    cat("\n>>> DOWNSAMPLE aleatorio:", ncells, "->", length(idx), "\n")
  }
} else {
  idx <- seq_len(ncells); cat("\n>>> TODAS:", ncells, "\n")
}

# Solo se copia la matriz cuando el recorte deja fuera alguna celula. Copiarla
# para quedarse con todas duplicaria varios gigabytes sin cambiar nada, y es
# suficiente para agotar la memoria de la sesion.
if (length(idx) < ncells) counts <- counts[, idx]
meta <- meta[idx, , drop=FALSE]
cat("Matriz de conteos:", nrow(counts), "genes x", ncol(counts), "celulas\n")

# --- Las coordenadas que ya trae el archivo ----------------------------------
# El archivo incluye un UMAP ya calculado. Se exporta para poder contrastarlo con
# el que el analisis calcula por su cuenta. Son dos columnas por celula, asi que
# no pesa. Como el resto de tablas del archivo trae su propio orden de filas, de
# modo que se recorta por nombre y no por posicion.
emb <- NULL; un <- NA_character_
reds <- Reductions(obj)
if (length(reds) == 0) {
  cat("AVISO: el objeto no tiene reducciones; no se exporta el UMAP incluido.\n")
} else {
  con_umap <- reds[grepl("umap", tolower(reds))]
  un <- if (length(con_umap) > 0) con_umap[1] else reds[1]
  e_all <- Embeddings(obj, un)
  if (!all(colnames(counts) %in% rownames(e_all))) {
    cat("AVISO: la reduccion '", un, "' no cubre todas las celulas; no se exporta.\n", sep="")
  } else {
    emb <- e_all[colnames(counts), , drop=FALSE]
  }
  rm(e_all)
}

# Ultima comprobacion de que todo sigue emparejado. Detenerse aqui es preferible
# a escribir en disco unos datos cruzados que nada volveria a revisar.
stopifnot(identical(colnames(counts), rownames(meta)))
if (!is.null(emb)) stopifnot(identical(colnames(counts), rownames(emb)))

# El archivo original se descarta antes de escribir, no despues. Contiene otra
# copia de la matriz que este analisis no usa, y liberarla deja sitio para la
# escritura. La matriz elegida sobrevive porque se guarda aparte.
rm(obj); gc()

Matrix::writeMM(counts, file.path(salida, "counts.mtx"))
write.csv(data.frame(gene=rownames(counts)), file.path(salida, "genes.csv"), row.names=FALSE)
write.csv(data.frame(barcode=colnames(counts)), file.path(salida, "barcodes.csv"), row.names=FALSE)
write.csv(meta, file.path(salida, "metadata.csv"))
if (!is.null(emb)) {
  write.csv(emb, file.path(salida, "umap.csv"))
  cat("Embedding exportado:", un, "(", ncol(emb), "dimensiones )\n")
}

rm(counts, meta); gc()
cat("\nListo.\n")
'''


R_COPIA_PRUEBAS = r'''
# R solo mira su carpeta personal de paquetes si el directorio ya existe. Se
# anade a mano para no depender de eso: ahi es donde quedaron instalados.
local({
  lib <- Sys.getenv("R_LIBS_USER")
  if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))
})
suppressMessages({library(SeuratObject)})

fuente  <- Sys.getenv("RDS_FUENTE", "RTX_zenodo.RDS")
destino <- Sys.getenv("RDS_DESTINO", "RTX_zenodo_debug.RDS")
n_deb   <- as.integer(Sys.getenv("N_CELLS_DEBUG", "4000"))

obj <- readRDS(fuente)
meta <- obj@meta.data

# Los identificadores validos de celula son los de la matriz de expresion. La
# tabla de metadatos de este dataset esta numerada del 1 en adelante, sin nombres
# de celula, asi que consultarla para saber que celulas hay no sirve: devuelve
# numeros de fila que la matriz no reconoce.
ar <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
celdas <- colnames(obj[[ar]])
if (is.null(celdas)) {
  cs0 <- tryCatch(GetAssayData(obj, assay=ar, layer="counts"),
                  error=function(e) tryCatch(GetAssayData(obj, assay=ar, slot="counts"),
                                             error=function(e2) NULL))
  if (is.null(cs0))
    stop("No se pudieron leer los nombres de celula del assay '", ar, "'.")
  celdas <- colnames(cs0); rm(cs0)
}
ncells <- length(celdas)
cat("Objeto completo:", ncells, "celulas\n")
if (ncells != nrow(meta))
  stop("El assay tiene ", ncells, " celulas y meta.data ", nrow(meta), " filas.")

clave <- rownames(meta)
usa_barcode <- FALSE
if (is.null(clave) || anyDuplicated(clave) || !setequal(clave, celdas)) {
  if (!("Barcode" %in% colnames(meta)))
    stop("rownames(meta.data) no identifican las celulas del assay y no hay una ",
         "columna 'Barcode' con la que alinear.")
  bc <- as.character(meta$Barcode)
  if (anyDuplicated(bc) || !setequal(bc, celdas))
    stop("Ni rownames(meta.data) ni la columna 'Barcode' identifican las celulas ",
         "del assay: no se puede muestrear con seguridad.")
  clave <- bc; usa_barcode <- TRUE
  cat("rownames(meta.data) no identifican celulas; se usa la columna 'Barcode'.\n")
}

# Se ponen los metadatos en el mismo orden que la matriz de expresion, y con los
# nombres de celula como identificador. Seurat recorta buscando esos nombres, de
# modo que sin este paso cualquier recorte devolveria un objeto vacio.
obj@meta.data <- meta[match(celdas, clave), , drop=FALSE]
rownames(obj@meta.data) <- celdas

set.seed(0)
cc <- grep("celltype|cell_type|cell.type|annotation|ident", colnames(obj@meta.data),
           ignore.case=TRUE, value=TRUE)
if (length(cc) > 0) {
  grp <- as.character(obj@meta.data[[cc[1]]]); fr <- n_deb / ncells
  # Caso aparte para los tipos celulares con una sola celula: la funcion de
  # muestreo de R interpreta un unico valor como un rango y devolveria una celula
  # distinta de la pedida.
  idx <- sort(unlist(lapply(split(seq_len(ncells), grp), function(ix)
    if (length(ix) == 1) ix else sample(ix, max(2, round(length(ix) * fr))))))
  cat("Muestreo estratificado por '", cc[1], "'\n", sep="")
} else {
  idx <- sort(sample(ncells, min(n_deb, ncells)))
  cat("Muestreo aleatorio\n")
}

celulas <- celdas[idx]
pequeno <- tryCatch(subset(obj, cells = celulas), error = function(e)
  stop("No se pudo subconjuntar el objeto Seurat (", conditionMessage(e), ").\n",
       "Vuelve a RDS_DEBUG = False y trabaja con el archivo completo."))
rm(obj); gc()

# Los metadatos vuelven a quedar como venian en el archivo original: numerados
# por fila y con el identificador de celula solo en su columna. Una copia de
# pruebas mas ordenada que el original no serviria para lo que existe, que es
# comprobar que el analisis maneja bien los datos tal como llegan.
if (usa_barcode) {
  ms <- pequeno@meta.data
  ms <- ms[order(match(rownames(ms), clave)), , drop=FALSE]
  rownames(ms) <- as.character(seq_len(nrow(ms)))
  pequeno@meta.data <- ms
}

cs <- tryCatch(GetAssayData(pequeno, assay=ar, layer="counts"),
               error=function(e) tryCatch(GetAssayData(pequeno, assay=ar, slot="counts"),
                                          error=function(e2) NULL))
ms <- pequeno@meta.data
if (is.null(cs)) stop("El subconjunto no conserva la capa 'counts'.")
if (ncol(cs) != nrow(ms))
  stop("El subconjunto quedo descuadrado: ", ncol(cs), " celulas en counts y ",
       nrow(ms), " filas de metadatos.")

# Se guarda sin comprimir. Ocupa mas espacio, pero se abre mucho mas rapido, que
# es justo lo que se busca en una copia pensada para repetir pruebas.
saveRDS(pequeno, destino, compress=FALSE)
cat("\nCopia de pruebas:", destino, "|", ncol(cs), "celulas x", nrow(cs), "genes\n")
'''


# ------------------------------------------------------------------
#   Carga de datos
# ------------------------------------------------------------------

def cargar_ejemplo():
    """Medula osea publica. Se descarga sola la primera vez."""
    import anndata as ad
    import numpy as np
    import pooch
    import scanpy as sc

    np.random.seed(0)
    EX = pooch.create(path=pooch.os_cache("scverse_tutorials"),
                      base_url="doi:10.6084/m9.figshare.22716739.v1/")
    EX.load_registry_from_doi()
    samples = {"s1d1": "s1d1_filtered_feature_bc_matrix.h5",
               "s1d3": "s1d3_filtered_feature_bc_matrix.h5"}
    adatas = {}
    for sid, fn in samples.items():
        a = sc.read_10x_h5(EX.fetch(fn))
        a.var_names_make_unique()
        adatas[sid] = a
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    print("adata:", adata.shape)
    return adata, None, None, None


def descargar_rds(cfg):
    """Trae el archivo del estudio desde Zenodo si no esta ya en disco."""
    rds = os.path.join(cfg.dir_datos, "RTX_zenodo.RDS")
    url = "https://zenodo.org/records/17868028/files/RTX_zenodo.RDS?download=1"

    def es_rds(path):
        """saveRDS escribe gzip (1f 8b) o RDS sin comprimir ('RDX2'/'RDX3')."""
        with open(path, "rb") as fh:
            m = fh.read(4)
        return m[:2] == b"\x1f\x8b" or m[:3] in (b"RDX", b"RDA")

    def bajar():
        parcial = rds + ".part"
        if _HEARTBEAT is not None:
            _HEARTBEAT.marcar("Descargando el archivo del estudio",
                              "depende de la conexion; suelen ser varios minutos")
        print("Descargando el archivo del estudio (~1.6 GB). Esto tarda un rato...")
        try:
            with urllib.request.urlopen(url, timeout=300) as resp:
                if resp.status != 200:
                    raise RuntimeError(f"HTTP {resp.status} desde Zenodo")
                declarado = resp.headers.get("Content-Length")
                declarado = int(declarado) if declarado else None
                leido = 0
                with open(parcial, "wb") as fh:
                    while True:
                        trozo = resp.read(1 << 22)
                        if not trozo:
                            break
                        fh.write(trozo)
                        leido += len(trozo)
                        if declarado and leido % (1 << 28) < (1 << 22):
                            print(f"  {leido/1e9:.2f} / {declarado/1e9:.2f} GB")
        except urllib.error.URLError as e:
            if os.path.exists(parcial):
                os.remove(parcial)
            raise RuntimeError(f"No se pudo descargar el archivo desde Zenodo:\n  {e}") from e

        real = os.path.getsize(parcial)
        if declarado is not None and real != declarado:
            os.remove(parcial)
            raise RuntimeError(
                f"Descarga truncada: {real:,} de {declarado:,} bytes. Vuelve a ejecutar."
            )
        if not es_rds(parcial):
            os.remove(parcial)
            raise RuntimeError(
                "Lo descargado no es un archivo RDS (puede ser una pagina de error).\n"
                f"Comprueba a mano: {url}"
            )
        os.replace(parcial, rds)

    if not os.path.exists(rds):
        bajar()
    elif not es_rds(rds):
        # Una descarga cortada de un intento anterior queda en disco con el nombre
        # correcto. Se comprueba el formato para no arrastrar ese archivo a medias.
        print("El archivo local esta corrupto -> se descarga de nuevo")
        os.remove(rds)
        bajar()

    print(f"Archivo listo: {os.path.getsize(rds)/1e9:.2f} GB (formato validado)")
    return rds


def preparar_r(cfg):
    """Comprueba que R esta disponible y con los paquetes que hacen falta."""
    rscript = buscar_rscript()
    if rscript is None:
        raise RuntimeError(
            "No se encontro R en este equipo, y el modo 'zenodo' lo necesita para\n"
            "leer el archivo del estudio, que viene en formato de Seurat.\n\n"
            "Que hacer:\n"
            "  1) Instala R desde https://cran.r-project.org/bin/windows/base/\n"
            "  2) Vuelve a ejecutar este script. Los paquetes de R que falten se\n"
            "     instalan solos.\n\n"
            "Si prefieres no instalar R, usa --modo ejemplo: ese modo no lo necesita."
        )
    print(f"R encontrado: {rscript}")

    # Los paquetes se instalan en la carpeta personal del usuario, no en la de R.
    # La de R esta dentro de Archivos de programa y solo se puede escribir como
    # administrador. Cuando R se usa a mano ofrece crear la personal y seguir; al
    # llamarlo desde aqui no hay nadie a quien preguntar, asi que se crea antes.
    codigo = (
        'lib <- Sys.getenv("R_LIBS_USER"); '
        'if (!dir.exists(lib)) dir.create(lib, recursive=TRUE, showWarnings=FALSE); '
        '.libPaths(c(lib, .libPaths())); '
        'for (p in c("SeuratObject","Matrix")) '
        'if (!requireNamespace(p, quietly=TRUE)) '
        'install.packages(p, lib=lib, repos="https://cloud.r-project.org"); '
        'ok <- all(sapply(c("SeuratObject","Matrix"), requireNamespace, quietly=TRUE)); '
        'if (!ok) quit(status=1); cat("paquetes R OK\\n")'
    )
    if _HEARTBEAT is not None:
        _HEARTBEAT.marcar("Preparando paquetes de R",
                          "segundos si ya estan instalados, unos minutos si hay que instalarlos")
    proc = ejecutar_streaming([rscript, "-e", codigo])
    if proc.returncode != 0:
        salida = (proc.stdout or "") + (proc.stderr or "")
        if "not writable" in salida or "unable to install" in salida:
            pista = (
                "\nR no pudo escribir donde guarda sus paquetes. Prueba una de estas:\n"
                "  1) Abre R o RStudio a mano una vez y ejecuta:\n"
                "         install.packages('SeuratObject')\n"
                "     Acepta cuando pregunte si quieres usar una carpeta personal.\n"
                "  2) O ejecuta este script como administrador."
            )
        elif "unable to access index" in salida or "cannot open URL" in salida:
            pista = (
                "\nR no pudo conectarse a su repositorio de paquetes. Comprueba la\n"
                "conexion a internet, y si hay un proxy o un antivirus que este\n"
                "bloqueando el acceso a cloud.r-project.org."
            )
        else:
            pista = "\nRevisa el mensaje de R de aqui arriba: suele decir que falta."
        raise RuntimeError(
            "No se pudieron preparar los paquetes de R.\n"
            "Sin ellos no se puede leer el archivo del estudio." + pista
        )
    return rscript


def convertir_rds(cfg, rscript, rds):
    """Ejecuta el script de R que traduce el archivo a formatos que Python lee."""
    import psutil

    # La memoria disponible se vuelve a medir aqui, justo antes del paso que la
    # consume, y no se reutiliza lo que se midio al arrancar. Un valor heredado
    # puede venir de otro momento y dejar pasar una corrida condenada a quedarse
    # sin memoria.
    ram_gb = psutil.virtual_memory().total / 1e9
    if cfg.n_cells_max == 0 and ram_gb < 20 and not cfg.rds_debug:
        raise RuntimeError(
            f"Pediste las 169.513 celulas completas (N_CELLS_MAX=0) pero este equipo\n"
            f"tiene {ram_gb:.1f} GB de RAM. R se queda sin memoria al leer el archivo\n"
            f"y el sistema corta el proceso.\n\n"
            f"Opciones:\n"
            f"  1) Pon un numero en N_CELLS_MAX, por ejemplo 60000. El muestreo es\n"
            f"     estratificado por tipo celular, asi que conserva las proporciones\n"
            f"     de las poblaciones y sigue siendo valido para el analisis.\n"
            f"  2) Ejecutalo en un equipo con mas memoria."
        )

    fuente = rds
    if cfg.rds_debug:
        destino = os.path.join(cfg.dir_datos, cfg.rds_debug_file)
        if os.path.exists(destino):
            print(f"Ya existe la copia de pruebas ({os.path.getsize(destino)/1e6:.0f} MB): se reutiliza.")
            print("Borrala a mano si quieres regenerarla con otro N_CELLS_DEBUG.")
        else:
            print(f"Generando una copia de pruebas (~{cfg.n_cells_debug:,} celulas).")
            print("Esto lee el archivo completo UNA vez; las pruebas siguientes ya no.")
            script = os.path.join(cfg.dir_trabajo, "hacer_copia_pruebas.R")
            with open(script, "w", encoding="utf-8") as fh:
                fh.write(R_COPIA_PRUEBAS)
            entorno = dict(os.environ, RDS_FUENTE=rds, RDS_DESTINO=destino,
                           N_CELLS_DEBUG=str(cfg.n_cells_debug))
            proc = ejecutar_streaming([rscript, script], env=entorno)
            if proc.returncode != 0:
                raise RuntimeError(
                    f"No se pudo generar la copia de pruebas (codigo {proc.returncode}).\n"
                    "Pon RDS_DEBUG = False para trabajar con el archivo completo."
                )
        fuente = destino

    script = os.path.join(cfg.dir_trabajo, "convert_rds.R")
    with open(script, "w", encoding="utf-8") as fh:
        fh.write(R_CONVERTIR)

    entorno = dict(os.environ,
                   N_CELLS_MAX=str(cfg.n_cells_max),
                   RDS_FILE=fuente,
                   DIR_TRABAJO=cfg.dir_trabajo)
    if _HEARTBEAT is not None:
        _HEARTBEAT.marcar("Conversion en R", "entre 10 y 20 minutos con el archivo completo")
    print("Convirtiendo el archivo del estudio. La salida de R se muestra en vivo, "
          "linea por linea; ahi estan los grupos.\n")
    proc = ejecutar_streaming([rscript, script], env=entorno)

    # Si la conversion falla conviene detenerse aqui. Sin esta comprobacion el
    # error saldria en el paso siguiente, como un archivo que no aparece, y ese
    # mensaje no dice nada de lo que ocurrio en realidad.
    if proc.returncode != 0:
        if proc.returncode < 0 or proc.returncode == 3221225725:
            detalle = (
                "\nEl proceso se corto de golpe, casi siempre por quedarse sin memoria.\n"
                "No es un problema del archivo ni de los paquetes de R, y los mensajes\n"
                "de arriba no van a explicarlo porque R no llego a reaccionar.\n"
                "Reduce N_CELLS_MAX o usa un equipo con mas memoria."
            )
        else:
            detalle = (
                "\nCausas habituales: los paquetes de R no se instalaron bien, o el\n"
                "archivo no tiene la estructura esperada. Los mensajes de arriba\n"
                "deberian distinguir entre las dos."
            )
        raise RuntimeError(f"La conversion fallo (codigo {proc.returncode}).{detalle}")

    obligatorios = ["counts.mtx", "genes.csv", "barcodes.csv", "metadata.csv"]
    faltan = [f for f in obligatorios if not os.path.exists(os.path.join(cfg.dir_trabajo, f))]
    if faltan:
        raise RuntimeError(
            f"R termino sin error pero no genero: {faltan}\n"
            "Probablemente el archivo no tiene la estructura esperada."
        )
    for f in obligatorios + ["umap.csv"]:
        ruta = os.path.join(cfg.dir_trabajo, f)
        if os.path.exists(ruta):
            print(f"  {f}: {os.path.getsize(ruta)/1e6:.1f} MB")
        else:
            print(f"  {f}: no generado (el archivo no traia UMAP; se calculara uno nuevo)")


def reconstruir_adata(cfg):
    """Arma la tabla de trabajo a partir de lo que exporto R."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.io

    d = cfg.dir_trabajo
    X = scipy.io.mmread(os.path.join(d, "counts.mtx")).T.tocsr()
    genes = pd.read_csv(os.path.join(d, "genes.csv"))["gene"].astype(str).values
    bc = pd.read_csv(os.path.join(d, "barcodes.csv"))["barcode"].astype(str).values
    meta = pd.read_csv(os.path.join(d, "metadata.csv"), index_col=0)
    meta.index = meta.index.astype(str)

    # El emparejamiento entre celulas y anotaciones se comprueba, no se da por
    # hecho. Pegar los identificadores sobre los metadatos sin verificar que
    # coinciden dejaria a cada celula con la anotacion de otra, sin ningun aviso,
    # y el analisis entregaria resultados que parecen validos sin serlo.
    if len(meta) != len(bc):
        raise ValueError(f"metadata.csv tiene {len(meta):,} filas y barcodes.csv {len(bc):,}.")
    if not meta.index.is_unique:
        raise ValueError("metadata.csv tiene barcodes duplicados: no se puede alinear con seguridad.")
    if not (meta.index.values == bc).all():
        ausentes = set(bc) - set(meta.index)
        if ausentes:
            raise ValueError(
                f"{len(ausentes):,} barcodes de la matriz no estan en metadata.csv, "
                f"p.ej.: {sorted(ausentes)[:3]}"
            )
        print("AVISO: metadata.csv venia en otro orden que barcodes.csv; se realinea por barcode.")
        meta = meta.loc[bc]

    adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes))
    adata.obs_names = bc

    ruta_umap = os.path.join(d, "umap.csv")
    if os.path.exists(ruta_umap):
        umap = pd.read_csv(ruta_umap, index_col=0)
        umap.index = umap.index.astype(str)
        if not (umap.index.values == bc).all():
            if set(bc) - set(umap.index):
                raise ValueError("umap.csv no cubre todos los barcodes de la matriz.")
            print("AVISO: umap.csv venia en otro orden; se realinea por barcode.")
            umap = umap.loc[bc]
        adata.obsm["X_umap_authors"] = umap.values
        print(f"UMAP incluido en el archivo conservado: {umap.shape[1]} dimensiones")
        del umap
    else:
        print("Sin UMAP en el archivo; se calculara uno nuevo.")
    del X, meta
    gc.collect()
    return adata


def resolver_columnas(adata, cfg):
    """Averigua que columna guarda cada dato. Nunca devuelve None en silencio."""

    def resumen():
        lineas = []
        for c in adata.obs.columns:
            v = adata.obs[c]
            n = v.nunique(dropna=True)
            if n <= 12:
                lineas.append(f"    {c!r:38s} ({n:2d} valores) -> {sorted(map(str, v.dropna().unique()))}")
            else:
                lineas.append(f"    {c!r:38s} ({n:4d} valores, continuo/ID)")
        return "\n".join(lineas)

    def resolver(etiqueta, override, claves, obligatoria):
        if override is not None:
            if override not in adata.obs.columns:
                raise KeyError(
                    f"La columna {override!r} que fijaste para '{etiqueta}' no existe.\n"
                    f"Columnas disponibles:\n{resumen()}"
                )
            print(f"  {etiqueta:6s}: {override!r}  (fijada a mano)")
            return override

        candidatas = [c for c in adata.obs.columns if any(k in c.lower() for k in claves)]
        if len(candidatas) == 1:
            print(f"  {etiqueta:6s}: {candidatas[0]!r}  (detectada sola)")
            return candidatas[0]
        if len(candidatas) > 1:
            print(f"  {etiqueta:6s}: {candidatas[0]!r}  (detectada sola; habia varias "
                  f"{candidatas} -> si no es la correcta fijala en el PANEL DE CONTROL)")
            return candidatas[0]

        msg = (f"No se encontro ninguna columna para '{etiqueta}'.\n"
               f"Este analisis se escribio para el archivo del estudio; otro dataset\n"
               f"puede nombrar sus columnas de otra forma.\n"
               f"Solucion: elige la columna correcta de esta lista y ponla en el\n"
               f"PANEL DE CONTROL como COL_{etiqueta.upper()} = 'nombre_de_columna'.\n"
               f"Columnas disponibles:\n{resumen()}")
        if obligatoria:
            raise KeyError(msg)
        print(f"  {etiqueta:6s}: NO encontrada (opcional)\n{msg}\n")
        return None

    print("Resolviendo columnas de metadatos:")
    col_cond = resolver("cond", cfg.col_cond,
                        ["disease", "condition", "group", "sle", "status", "diagnosis"],
                        obligatoria=False)
    col_time = resolver("time", cfg.col_time,
                        ["time", "visit", "day", "week", "treatment", "point"],
                        obligatoria=False)
    # El tipo celular es imprescindible: sin el no hay forma de saber que celulas
    # son B, y las dos figuras del estudio y todo SCENIC dependen de eso.
    col_ctype = resolver("ctype", cfg.col_ctype,
                         ["celltype", "cell_type", "cell.type", "annotation", "ident"],
                         obligatoria=True)
    return col_cond, col_time, col_ctype


# ------------------------------------------------------------------
#   Control de calidad, normalizacion y las 12 graficas
# ------------------------------------------------------------------

def marcadores_por_modo(adata, modo):
    def present(md):
        out = {k: [g for g in v if g in adata.var_names] for k, v in md.items()}
        return {k: v for k, v in out.items() if v}

    if modo == "ejemplo":
        return present({
            "CD14+ Mono": ["FCN1", "CD14"], "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
            "cDC2": ["CST3", "COTL1", "LYZ", "CLEC10A", "FCER1A"],
            "Erythroblast": ["MKI67", "HBA1", "HBB"],
            "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
            "NK": ["GNLY", "NKG7", "CD247", "TYROBP", "KLRG1"],
            "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
            "Plasma cells": ["MZB1", "HSP90B1", "PRDM1", "IGKC", "JCHAIN"],
            "CD4+ T": ["CD4", "IL7R", "TRBC2"],
            "CD8+ T": ["CD8A", "CD8B", "GZMK", "CCL5", "GZMB"],
            "T naive": ["LEF1", "CCR7", "TCF7"], "pDC": ["IL3RA", "COBLL1", "TCF4"],
        })
    return present({
        "T CD4": ["CD3D", "CD4", "IL7R"], "T CD8": ["CD3D", "CD8A", "GZMK"],
        "B": ["MS4A1", "CD79A", "CD79B"], "NK": ["GNLY", "NKG7", "KLRD1"],
        "Mono": ["CD14", "LYZ", "FCGR3A"], "DC": ["FCER1A", "CST3"],
        "Plasma": ["MZB1", "JCHAIN", "IGHG1"],
    })


def diagnostico_matriz(adata):
    """Comprueba si la matriz trae conteos sin procesar y si quedan genes MT-."""
    import numpy as np

    # Primera comprobacion: si la matriz trae conteos sin procesar, el total que
    # calcula scanpy tiene que coincidir con la columna nCount_RNA que ya venia
    # en los metadatos. Que no coincidan indica que la matriz llega transformada.
    if "nCount_RNA" in adata.obs.columns:
        corr = np.corrcoef(adata.obs["total_counts"], adata.obs["nCount_RNA"])[0, 1]
        print("Correlacion total_counts (scanpy) vs nCount_RNA (del archivo):", corr)
        cols = [c for c in ["total_counts", "nCount_RNA", "n_genes_by_counts", "nFeature_RNA"]
                if c in adata.obs.columns]
        print(adata.obs[cols].describe())
        if corr < 0.99:
            print("\n  total_counts no sigue a nCount_RNA, senal de que la matriz no son\n"
                  "  conteos sin procesar. En este dataset es lo esperado: la matriz viene\n"
                  "  ya normalizada de origen y la version sin procesar no se publico.\n"
                  "  No hace falta tocar nada: la conversion elige la matriz con enteros si\n"
                  "  existe alguna, y el paso siguiente mide el estado real y aplica solo la\n"
                  "  normalizacion que falte.\n"
                  "  Con esta matriz, 'total_counts' no es profundidad de secuenciacion;\n"
                  "  para eso esta nCount_RNA, que si la conserva.")
    else:
        print("'nCount_RNA' no esta en los metadatos; no se puede comparar.")

    # Segunda comprobacion: scanpy calcula el porcentaje mitocondrial buscando los
    # genes cuyo nombre empieza por MT-. Si esa busqueda no encuentra ninguno, la
    # metrica sale cero en todas las celulas sin que nada falle. Se cuenta cuantos
    # hay en la matriz para poder distinguir un cero real de una columna vacia.
    mt = int(adata.var_names.str.startswith("MT-").sum())
    print("\nGenes MT- en la matriz:", mt)
    if "percent.mt" in adata.obs.columns:
        print(adata.obs[["pct_counts_mt", "percent.mt"]].describe())
    if mt == 0:
        print("\n  No hay genes MT- en la matriz, por eso pct_counts_mt sale 0 en todas las\n"
              "  celulas. Ya se quitaron antes de publicar el archivo, asi que el filtro de\n"
              "  calidad por mitocondrial no esta filtrando nada. La columna 'percent.mt' de\n"
              "  los metadatos si conserva el valor que tenian, y es la que hay que usar\n"
              "  para hablar de calidad mitocondrial en este dataset.")


def normalizar_si_hace_falta(adata):
    """Mide en que estado llega la matriz y aplica solo lo que le falte.

    El analisis espera la matriz normalizada y en escala logaritmica, pero no
    todos los datasets llegan en el mismo estado: este viene con los dos pasos ya
    aplicados. Repetirlos transformaria la matriz dos veces y cambiaria la
    seleccion de genes, el PCA y el agrupamiento que salen despues.
    Se distinguen tres casos y hay un cuarto que detiene la ejecucion:
        valores enteros            -> sin procesar: se normaliza y se aplica log1p
        decimales, suma constante  -> normalizada sin log1p: solo log1p
        igual pero tras deshacer   -> ya viene completa: no se toca
        ninguno de los anteriores  -> estado desconocido: se para
    La decision sale de los datos y no de una variable de configuracion, asi que
    sigue valiendo con otro dataset.
    """
    import numpy as np
    import scanpy as sc
    import scipy.sparse as sp

    m = adata.X[:min(200, adata.n_obs)]
    vals = m.data if sp.issparse(m) else np.asarray(m).ravel()
    vals = vals[vals != 0]
    entera = vals.size > 0 and np.allclose(vals, np.round(vals))

    def sumas(mat, deshacer_log):
        mm = mat.copy()
        if sp.issparse(mm):
            if deshacer_log:
                mm.data = np.expm1(mm.data)      # expm1(0)=0: no rompe la dispersion
            return np.asarray(mm.sum(axis=1)).ravel()
        a = np.asarray(mm, dtype=float)
        if deshacer_log:
            a = np.expm1(a)
        return a.sum(axis=1)

    def es_constante(s):
        """True si todas las celulas suman lo mismo: la huella de una normalizacion."""
        s = s[np.isfinite(s)]
        # Se exige que la media sea positiva. Un grupo de celulas todas a cero
        # tiene desviacion cero y pasaria por normalizado, saltandose el paso.
        if s.size == 0 or float(np.mean(s)) <= 0:
            return False
        # El margen no es mas estrecho porque la limpieza de genes del paso
        # anterior descarta genes poco frecuentes, y con ellos se va una parte
        # pequena de cada celula. La suma deja de ser exacta por poco, y el umbral
        # tiene que admitir esa desviacion sin llegar a confundirla con una matriz
        # que nunca se normalizo, que se aparta cien veces mas.
        return float(np.std(s)) / float(np.mean(s)) < 1e-2

    if entera:
        estado = "conteos sin procesar"
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        hecho = "se normaliza (CP10K) y se aplica log1p"
    elif es_constante(sumas(m, True)):
        estado = "ya normalizada y con log1p aplicado en origen"
        adata.layers["normalizada_en_origen"] = adata.X.copy()
        # scanpy deja constancia de haber aplicado el logaritmo, y otras funciones
        # suyas consultan ese registro mas adelante, entre ellas la de expresion
        # diferencial. Como aqui el logaritmo venia de origen, la anotacion se
        # escribe a mano para que esas funciones encuentren lo que esperan.
        adata.uns.setdefault("log1p", {"base": None})
        hecho = "no se toca: ya viene normalizada y logaritmica"
    elif es_constante(sumas(m, False)):
        estado = "normalizada en origen, pero sin log1p"
        adata.layers["normalizada_en_origen"] = adata.X.copy()
        sc.pp.log1p(adata)
        hecho = "solo se aplica log1p"
    else:
        raise RuntimeError(
            "No se pudo determinar como viene la matriz de expresion.\n"
            "No son enteros (luego no son conteos sin procesar) y las celulas tampoco\n"
            "suman lo mismo, ni en crudo ni deshaciendo el logaritmo, asi que no es una\n"
            "normalizacion estandar. Puede estar escalada o venir de otro procesamiento.\n"
            "Se para aqui a proposito: normalizar por encima de una matriz en un estado\n"
            "desconocido produce una seleccion de genes, un PCA y unos grupos celulares\n"
            "que parecen resultados y no lo son.\n"
            "Revisa la salida de la conversion: ahi se lista matriz por matriz cual trae\n"
            "valores enteros."
        )

    print(f"Matriz de entrada: {estado}")
    print(f"  -> {hecho}")


def analisis_scanpy(adata, cfg, fig, col_cond, col_ctype, markers):
    """Control de calidad, agrupamiento y las 12 graficas."""
    import numpy as np
    import scanpy as sc

    paso("Grafica 1: violines de control de calidad")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    if not cfg.saltar_graficas:
        sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                     jitter=0.4, multi_panel=True, show=False)
        fig.guardar("01_violines_qc")
    else:
        print("[SALTAR_GRAFICAS] omitida. Las metricas de calidad si se calcularon.")

    if cfg.modo == "zenodo":
        paso("Diagnostico de la matriz")
        diagnostico_matriz(adata)

    paso("Grafica 2: dispersion de control de calidad")
    if not cfg.saltar_graficas:
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
        fig.guardar("02_scatter_qc")
    else:
        print("[SALTAR_GRAFICAS] omitida.")

    # El dataset del estudio llega ya depurado, asi que solo se aplica una
    # limpieza suave de genes. Los datos de ejemplo vienen sin procesar y
    # necesitan el control de calidad completo.
    if cfg.modo == "ejemplo":
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=3)
        if not cfg.saltar_graficas:
            # Anade la columna 'predicted_doublet' con el resultado de la
            # deteccion. Solo marca: no elimina filas, asi que la matriz que sigue
            # es la misma.
            sc.pp.scrublet(adata, batch_key="sample", random_state=0)
            print("Dobletes:", int(adata.obs["predicted_doublet"].sum()))
    else:
        sc.pp.filter_genes(adata, min_cells=3)

    paso("Grafica 3: seleccion de genes")
    normalizar_si_hace_falta(adata)
    batch = "sample" if "sample" in adata.obs.columns else None
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=batch)
    if not cfg.saltar_graficas:
        sc.pl.highly_variable_genes(adata, show=False)
        fig.guardar("03_genes_variables")
    else:
        print("[SALTAR_GRAFICAS] omitida. La normalizacion y la seleccion si se hicieron.")

    if cfg.saltar_graficas:
        print("\n[SALTAR_GRAFICAS] Se omiten las graficas 4 a 12 y los calculos que solo")
        print("las alimentan. El analisis continua directo con SCENIC.")
        return None

    paso("Grafica 4: componentes principales")
    sc.tl.pca(adata, random_state=0)
    color_group = "sample" if "sample" in adata.obs.columns else (col_cond or "pct_counts_mt")
    # Cada color se empareja por posicion con un par de componentes, asi que la
    # lista se repite para dibujar los mismos datos sobre PC1/PC2 y PC3/PC4.
    sc.pl.pca(adata, color=[color_group, color_group, "pct_counts_mt", "pct_counts_mt"],
              dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)], ncols=2, size=3, show=False)
    fig.guardar("04a_pca")
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    fig.guardar("04b_pca_varianza")

    paso("Grafica 5: grafo de vecinos")
    sc.pp.neighbors(adata, random_state=0)
    sc.tl.umap(adata, random_state=0)
    sc.pl.umap(adata, color=color_group, edges=True, edges_width=0.05,
               title="Grafo de vecinos mas cercanos", show=False)
    fig.guardar("05_grafo_vecinos")

    paso("Agrupamiento de celulas")
    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res,
                     flavor="igraph", n_iterations=2, random_state=0)
        print(f"  resolucion {res}: {adata.obs[f'leiden_res_{res:4.2f}'].nunique()} grupos")

    paso("Grafica 6: metricas de calidad sobre el mapa")
    qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    if "doublet_score" in adata.obs.columns:
        qc_cols.append("doublet_score")
    sc.pl.umap(adata, color=qc_cols, ncols=2, size=3, show=False)
    fig.guardar("06_umap_qc")

    paso("Grafica 7: anotacion de los grupos")
    sc.pl.umap(adata, color=cfg.cluster_key, legend_loc="on data",
               title="Grupos (numeros)", show=False)
    fig.guardar("07a_grupos")
    if cfg.modo == "ejemplo":
        mapa = {"0": "Lymphocytes", "1": "Monocytes", "2": "Erythroid", "3": "B Cells"}
        # El resultado se pasa por texto antes de rellenar los huecos. Al aplicar
        # una correspondencia sobre una columna de categorias, pandas devuelve otra
        # columna de categorias, y ahi no se puede escribir un valor que no este ya
        # en la lista de categorias validas.
        adata.obs["cell_type"] = (adata.obs["leiden_res_0.02"].map(mapa)
                                  .astype(object).fillna("Unknown").astype("category"))
        sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False)
        fig.guardar("07b_anotacion")
    elif col_ctype:
        sc.pl.umap(adata, color=col_ctype, legend_loc="right margin",
                   title="Anotacion incluida en el archivo", show=False)
        fig.guardar("07b_anotacion")

    paso("Grafica 8: marcadores")
    # Cada punto resume dos cosas a la vez: el tamano indica en que fraccion de las
    # celulas del grupo se detecta el gen, y el color, cuanto se expresa.
    sc.pl.dotplot(adata, markers, groupby=cfg.cluster_key, standard_scale="var", show=False)
    fig.guardar("08a_dotplot_marcadores")
    # Los mismos genes sobre el mapa celular, para ver donde se concentra cada uno.
    genes_umap = [g for gs in markers.values() for g in gs][:6]
    sc.pl.umap(adata, color=genes_umap, ncols=3, size=3, show=False)
    fig.guardar("08b_marcadores_umap")

    paso("Expresion diferencial")
    # Se descartan por prefijo de nombre varias familias de genes que aparecen en
    # todas las celulas con valores altos. Coparian las primeras posiciones del
    # ranking sin diferenciar un grupo de otro.
    hk = (adata.var_names.str.startswith(("MT-", "RPS", "RPL", "MRPS", "MRPL")) |
          adata.var_names.str.contains("^HB[^(P)]"))
    adata_de = adata[:, ~hk].copy()
    print(f"Genes de mantenimiento excluidos: {int(hk.sum())} | genes usados: {adata_de.n_vars}")
    sc.tl.rank_genes_groups(adata_de, groupby=cfg.cluster_key, method="wilcoxon")

    paso("Grafica 9: dotplot de genes diferenciales")
    sc.pl.rank_genes_groups_dotplot(adata_de, groupby=cfg.cluster_key,
                                    standard_scale="var", n_genes=5, show=False)
    fig.guardar("09_dotplot_diferenciales")

    paso("Grafica 10: tablas de expresion diferencial")
    de_all = sc.get.rank_genes_groups_df(adata_de, group=None)
    de_all["significativo"] = de_all["pvals_adj"] < 0.05

    # Tabla completa con todos los genes evaluados en cada grupo, resulten
    # significativos o no. Es la referencia para cualquier consulta posterior,
    # porque no aplica ningun recorte.
    out_csv = os.path.join(cfg.dir_salida, f"DE_completo_{cfg.modo}.csv")
    de_all.to_csv(out_csv, index=False)

    # Resumen con los genes mas caracteristicos de cada grupo. Cien es un limite,
    # no una cantidad a completar: primero se descarta lo que no alcanza
    # significancia estadistica y solo despues se recorta la lista. Un grupo con
    # pocos genes significativos aparece con pocos, y no se rellena con ruido.
    top100 = (de_all[de_all["significativo"]]
              .groupby("group", group_keys=False)
              .head(100))
    out_top = os.path.join(cfg.dir_salida, f"top100_DE_{cfg.modo}.csv")
    top100.to_csv(out_top, index=False)

    n_sig = de_all.groupby("group")["significativo"].sum().astype(int)
    print(f"Guardado: {out_csv}  ({len(de_all):,} filas)")
    print(f"Guardado: {out_top}  ({len(top100):,} filas)")
    print("\nGenes significativos por grupo:")
    print(n_sig.to_string())

    incompletos = n_sig[(n_sig < 100) & (n_sig > 0)]
    if len(incompletos) > 0:
        print(f"\n{len(incompletos)} grupo(s) con menos de 100 genes significativos: "
              "el resumen trae solo los que hay, no se rellena.")
    sin_sig = n_sig[n_sig == 0]
    if len(sin_sig) > 0:
        print(f"\nAVISO: {len(sin_sig)} grupo(s) sin ningun gen significativo: "
              f"{list(sin_sig.index)}")
        print("No apareceran en el resumen. Siguen en la tabla completa con su ranking real.")

    sc.pl.rank_genes_groups(adata_de, n_genes=20, sharey=False, show=False)
    fig.guardar("10_rankings_diferenciales")

    paso("Grafica 11: mapa de calor")
    sc.pl.rank_genes_groups_heatmap(adata_de, n_genes=3, groupby=cfg.cluster_key,
                                    standard_scale="var", show_gene_labels=True,
                                    swap_axes=True, figsize=(12, 14), show=False)
    fig.guardar("11_heatmap_diferenciales")

    paso("Grafica 12: trayectoria")
    # Es la ultima grafica y nada de lo que viene despues depende de ella, asi que
    # un fallo aqui se avisa y se sigue. Sin esto, un tropiezo en este punto
    # tiraria tambien el analisis de SCENIC, que es lo que mas tarda.
    try:
        # Mide que grupos celulares estan conectados entre si y con que fuerza.
        sc.tl.paga(adata, groups=cfg.cluster_key)
        sc.pl.paga(adata, color=cfg.cluster_key, title="Conexiones entre grupos", show=False)
        fig.guardar("12a_paga")
        # El mapa se recalcula tomando esas conexiones como posicion de partida, en
        # vez de la inicializacion por defecto, para que el dibujo respete la
        # estructura del grafo anterior.
        sc.tl.umap(adata, init_pos="paga", random_state=0)
        # El pseudotiempo ordena las celulas a lo largo del grafo y necesita un
        # grupo de partida. Se toma el de ROOT_CLUSTER si esta fijado; si no, el
        # menos numeroso.
        root = cfg.root_cluster if cfg.root_cluster is not None else adata.obs[cfg.cluster_key].value_counts().index[-1]
        coincidencias = np.flatnonzero(adata.obs[cfg.cluster_key].astype(str) == str(root))
        if coincidencias.size == 0:
            raise ValueError(
                f"El grupo de partida {root!r} no existe en {cfg.cluster_key!r}. "
                f"Grupos disponibles: {sorted(adata.obs[cfg.cluster_key].astype(str).unique())}"
            )
        adata.uns["iroot"] = int(coincidencias[0])
        sc.tl.dpt(adata)
        print("Grupo de partida del pseudotiempo:", root)
        sc.pl.umap(adata, color=[cfg.cluster_key, "dpt_pseudotime"],
                   legend_loc="on data", size=3, show=False)
        fig.guardar("12b_trayectoria")
    except Exception as e:
        print(f"AVISO: no se pudo calcular la trayectoria ({type(e).__name__}: {e})")
        print("       Suele pasar cuando los grupos quedan sin ninguna conexion entre si.")
        print("       Las demas graficas y el analisis de SCENIC no se ven afectados.")

    return adata_de


# ------------------------------------------------------------------
#   Las dos figuras del estudio
# ------------------------------------------------------------------

def figuras_estudio(adata, cfg, fig, col_time, col_ctype):
    """Subtipos de celulas B, y comparacion antes/despues del tratamiento."""
    import matplotlib.pyplot as plt
    import numpy as np
    import scanpy as sc

    paso("Figura 1: subtipos de celulas B")
    adata_B = seleccionar_celulas_B(adata, col_ctype)
    print(adata_B.obs[col_ctype].value_counts())
    sc.pl.umap(adata_B, color=col_ctype, size=8,
               title="Subtipos de celulas B", show=False)
    fig.guardar("fig1_subtipos_celulas_B")

    paso("Figura 2: comparacion antes y despues del tratamiento")
    celltype_subset = None      # p.ej. 'Memory B' para restringir a un tipo celular

    # La columna de momento tiene mas de dos valores y esta figura compara solo
    # dos, asi que hay que elegir cual va contra el estado previo. Se toma el
    # temprano porque es el unico con los 9 pacientes; al tardio le faltan 3,
    # repartidos entre las dos categorias de la columna Responder. Cambiar este
    # valor por el tardio es valido, contando con esos 3 pacientes menos.
    post_por_defecto = "Early"

    # --- El momento va dentro del identificador de muestra ----------------------
    # No hay columna propia: el identificador combina sujeto y momento separados
    # por un guion bajo. Parte de las muestras no llevan sufijo, y se reconocen por
    # no tener separador. La expresion regular exige ese separador, asi que esas
    # quedan sin valor y el filtro posterior las descarta solo.
    if col_time is None:
        if "sampleID" not in adata.obs.columns:
            raise RuntimeError(
                "No hay columna de momento ni 'sampleID' del que derivarla.\n"
                "Fija COL_TIME en el PANEL DE CONTROL con la columna correcta."
            )
        sid = adata.obs["sampleID"].astype(str)
        adata.obs["timepoint"] = sid.str.extract(r"^[^_]+_(.+)$", expand=False)
        derivados = adata.obs["timepoint"].dropna().unique()
        sin_momento = int(adata.obs["timepoint"].isna().sum())
        if len(derivados) == 0:
            raise RuntimeError(
                "Se intento derivar el momento de 'sampleID' pero ningun valor tiene el\n"
                "formato esperado. Valores reales:\n"
                + "\n".join(f"    {v!r}" for v in sorted(sid.unique())[:20])
                + "\n\nFija COL_TIME en el PANEL DE CONTROL con la columna correcta."
            )
        col_time = "timepoint"
        print(f"Momento derivado de 'sampleID' -> columna {col_time!r}")
        print(f"  momentos encontrados: {sorted(derivados)}")
        print(f"  celulas sin momento: {sin_momento:,}\n")

    valores = adata.obs[col_time].astype(str)
    disponibles = sorted(valores.unique())
    conteos = valores.value_counts()
    lista = "\n".join(f"    {v!r}: {conteos[v]:,} celulas" for v in disponibles)

    def elegir(termino, claves, excluir):
        """Busca UN valor que encaje. Si hay 0 o mas de 1, para y muestra las opciones."""
        cand = [v for v in disponibles
                if any(k in v.lower() for k in claves)
                and not any(x in v.lower() for x in excluir)]
        if len(cand) == 1:
            return cand[0]
        motivo = "no encaja ninguno" if not cand else f"encajan varios: {cand}"
        raise RuntimeError(
            f"No se pudo determinar automaticamente el valor de '{termino}' ({motivo}).\n"
            f"Valores reales de {col_time!r}:\n{lista}\n\n"
            f"Elige el que corresponda y ponlo en el PANEL DE CONTROL:\n"
            f"    PRE_LABEL  = '...'\n"
            f"    POST_LABEL = '...'"
        )

    # Al buscar el momento posterior se descarta cualquier etiqueta que empiece por
    # "pre", para que una palabra como "pretratamiento" no pase por posterior.
    pre = cfg.pre_label if cfg.pre_label is not None else elegir(
        "PRE", ["pre", "baseline", "before", "screening"], ["post"])

    # Buscar la palabra "post" no sirve en este dataset, donde los momentos se
    # llaman temprano y tardio. Se usa el valor elegido mas arriba, y solo si no
    # existe se recurre a la busqueda, que se detendra mostrando los disponibles.
    if cfg.post_label is not None:
        post = cfg.post_label
    elif post_por_defecto in disponibles:
        post = post_por_defecto
    else:
        post = elegir("POST", ["post", "after", "follow"], ["pre-", "pretreat"])

    for etiqueta, valor in (("PRE_LABEL", pre), ("POST_LABEL", post)):
        if valor not in disponibles:
            raise ValueError(f"{etiqueta}={valor!r} no existe en {col_time!r}.\n"
                             f"Valores reales:\n{lista}")
    print(f"Comparacion: {post!r}  vs  {pre!r}  (columna {col_time!r})")

    ad_de = adata
    if celltype_subset and col_ctype:
        ad_de = ad_de[ad_de.obs[col_ctype].astype(str) == celltype_subset]
        print(f"Restringido a {celltype_subset!r}: {ad_de.n_obs:,} celulas")
    ad_de = ad_de[ad_de.obs[col_time].astype(str).isin([pre, post])].copy()

    n_por_grupo = ad_de.obs[col_time].astype(str).value_counts()
    print(n_por_grupo)
    if n_por_grupo.min() < 30:
        raise ValueError(
            f"Uno de los grupos tiene solo {n_por_grupo.min()} celulas.\n"
            "Con tan pocas la comparacion no seria interpretable.\n"
            "Sube N_CELLS_MAX, o elige otros momentos."
        )

    # Se cuentan tambien los pacientes de cada grupo, no solo las celulas. Las
    # celulas de un mismo paciente no son observaciones independientes, asi que el
    # conteo de celulas por si solo no dice si la comparacion esta equilibrada.
    if "patientID" in ad_de.obs.columns:
        pac = ad_de.obs.groupby(ad_de.obs[col_time].astype(str), observed=True)["patientID"].nunique()
        print("\nPacientes por grupo:")
        print(pac.to_string())
        if pac.nunique() > 1:
            solo = {g: set(ad_de.obs.loc[ad_de.obs[col_time].astype(str) == g, "patientID"].unique())
                    for g in pac.index}
            faltan = solo[pre] - solo[post]
            print(f"\nAVISO: los grupos no tienen los mismos pacientes. En {post!r} faltan: "
                  f"{sorted(faltan) if faltan else 'ninguno'}.\n"
                  "       La comparacion deja de estar pareada y parte de la diferencia puede\n"
                  "       venir de que pacientes entran en cada lado.")

    sc.tl.rank_genes_groups(ad_de, groupby=col_time, groups=[post],
                            reference=pre, method="wilcoxon")
    de = sc.get.rank_genes_groups_df(ad_de, group=post)
    print(f"Genes evaluados: {len(de):,}")

    d = de.dropna(subset=["logfoldchanges", "pvals_adj"]).copy()
    d["nlp"] = -np.log10(d["pvals_adj"].clip(lower=1e-300))
    up = (d["logfoldchanges"] > 1) & (d["pvals_adj"] < 0.05)
    dn = (d["logfoldchanges"] < -1) & (d["pvals_adj"] < 0.05)
    plt.figure(figsize=(8, 6))
    plt.scatter(d["logfoldchanges"], d["nlp"], s=6, c="lightgray")
    plt.scatter(d.loc[up, "logfoldchanges"], d.loc[up, "nlp"], s=8, c="#c0392b", label="Up")
    plt.scatter(d.loc[dn, "logfoldchanges"], d.loc[dn, "nlp"], s=8, c="#2e5f9a", label="Down")
    for _, r in d[up | dn].nlargest(15, "nlp").iterrows():
        plt.text(r["logfoldchanges"], r["nlp"], r["names"], fontsize=7)
    plt.axvline(0, color="k", lw=.5)
    plt.axhline(-np.log10(0.05), color="k", ls="--", lw=.5)
    plt.xlabel(f"logFC ({post} - {pre})")
    plt.ylabel("-log10 p")
    plt.legend()
    plt.title("Antes y despues del tratamiento")
    plt.tight_layout()
    fig.guardar("fig2_volcano_pre_post")

    return adata_B


# ------------------------------------------------------------------
#   SCENIC
# ------------------------------------------------------------------

CHECKSUMS = {
    "motifs.tbl": "81eb754118e27e854974301b1400fcf519489f8be5249239671fb288cb501c31",
    # Este nombre no se puede cambiar. SCENIC deduce que tipo de base de datos
    # tiene delante a partir de la terminacion del archivo, y con cualquier otro
    # nombre se detiene nada mas empezar.
    "hg38.genes_vs_motifs.rankings.feather":
        "9c4026a3a8e25fe07cf96749644e2ca028b787410829b30b9932574dc6e78bdb",
    "expr_mat_tiny.loom": "ca57894cc828488d7aeb3ca58ad76a637f265c502688856d45a73d39f9483b4c",
    "allTFs_hg38.txt": None,
}

URLS = {
    "motifs.tbl":
        "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    "hg38.genes_vs_motifs.rankings.feather":
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/"
        "mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    "expr_mat_tiny.loom":
        "https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/expr_mat_tiny.loom",
    # La lista completa de factores de transcripcion humanos, que se usa en los dos
    # modos. La lista reducida que acompana al protocolo de SCENIC hoy contiene un
    # solo factor, y con uno solo no hay red que reconstruir.
    "allTFs_hg38.txt":
        "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
}


def descargar_bases_scenic(cfg):
    carpeta = os.path.join(cfg.dir_datos, "scenic_data")
    os.makedirs(carpeta, exist_ok=True)
    necesarios = ["motifs.tbl", "hg38.genes_vs_motifs.rankings.feather", "allTFs_hg38.txt"]
    if cfg.modo == "ejemplo":
        necesarios.append("expr_mat_tiny.loom")
    print("Las bases de datos de SCENIC pesan varios GB. Solo se descargan una vez.")
    for fn in necesarios:
        descargar_verificado(fn, URLS[fn], os.path.join(carpeta, fn),
                             CHECKSUMS.get(fn), cfg.verificar_checksums)
    print("Bases de datos de SCENIC listas y verificadas.")
    return carpeta


def preparar_matriz_scenic(adata, adata_B, cfg, carpeta, col_ctype):
    """Arma la matriz que entra a SCENIC y la lista de reguladores."""
    import numpy as np
    import pandas as pd

    if cfg.modo == "ejemplo":
        import loompy
        with loompy.connect(os.path.join(carpeta, "expr_mat_tiny.loom")) as ds:
            ex_matrix = pd.DataFrame(ds[:, :].T, index=ds.ca["CellID"], columns=ds.ra["Gene"])
        n_est = 500
    else:
        import scanpy as sc
        if adata_B is None:
            adata_B = seleccionar_celulas_B(adata, col_ctype)
        ad_s = adata_B.copy()
        if cfg.scenic_downsample and ad_s.n_obs > cfg.scenic_n_cells:
            sc.pp.subsample(ad_s, n_obs=cfg.scenic_n_cells, random_state=0)
            print(f"Recorte a {ad_s.n_obs:,} celulas B")
        else:
            print(f"Sin recorte: {ad_s.n_obs:,} celulas B (puede tardar horas)")

        # El pool de genes se recorta contra la base de datos de SCENIC ANTES de
        # elegir los mas variables, no despues. El orden es lo que importa aqui.
        # El filtro intermedio de SCENIC descarta cualquier grupo de genes del que
        # no reconozca al menos el 80% en esa base. Los genes de inmunoglobulina
        # copan las primeras posiciones por variabilidad y no estan en la base, asi
        # que colandose en la seleccion arrastran a casi todos los grupos por
        # debajo del umbral y el paso siguiente devuelve una lista vacia.
        # Se excluyen tambien por prefijo de nombre, ademas de por ausencia en la
        # base, porque no todas sus variantes faltan en ella.
        from ctxcore.rnkdb import FeatherRankingDatabase
        db = FeatherRankingDatabase(
            fname=os.path.join(carpeta, "hg38.genes_vs_motifs.rankings.feather"), name="hg38")
        genes_db = set(db.genes)
        es_ig = ad_s.var_names.str.match(r"^IG[HKL][VDJ]")
        en_db = ad_s.var_names.isin(genes_db)
        print(f"Genes en la matriz de celulas B: {ad_s.n_vars:,} | en la base de SCENIC: "
              f"{en_db.sum():,} ({en_db.mean():.0%}) | inmunoglobulinas excluidas: {es_ig.sum():,}")
        ad_s = ad_s[:, en_db & ~es_ig].copy()
        if ad_s.n_vars < 1000:
            print(f"\nAVISO: solo {ad_s.n_vars} genes quedan tras el filtro. "
                  "Probablemente se validaran pocos regulones o ninguno.")

        sc.pp.highly_variable_genes(ad_s, n_top_genes=min(cfg.scenic_n_genes, ad_s.n_vars - 1))
        ad_s = ad_s[:, ad_s.var.highly_variable].copy()
        Xd = ad_s.X.toarray() if hasattr(ad_s.X, "toarray") else np.asarray(ad_s.X)
        ex_matrix = pd.DataFrame(Xd, index=ad_s.obs_names.astype(str),
                                 columns=ad_s.var_names.astype(str))
        n_est = 200

    all_tfs = pd.read_csv(os.path.join(carpeta, "allTFs_hg38.txt"), header=None).iloc[:, 0].tolist()
    tf_names = [t for t in all_tfs if t in ex_matrix.columns]
    print(f"Matriz: {ex_matrix.shape[0]:,} celulas x {ex_matrix.shape[1]:,} genes | "
          f"reguladores: {len(tf_names)}")

    # SCENIC necesita un numero suficiente de genes para detectar patrones. Con muy
    # pocos, lo esperable es que no valide ningun regulon.
    if ex_matrix.shape[1] < 1000:
        print(f"\nAVISO: solo {ex_matrix.shape[1]} genes en la matriz. Probablemente no se\n"
              "valide ningun regulon. Es lo esperable con los datos de ejemplo.")
    return ex_matrix, tf_names, n_est


def adaptar_from_delayed():
    """Adapta una funcion de dask a como la llaman arboreto y SCENIC.

    Las dos librerias son de 2020 y no se han actualizado. Ambas llaman a la
    misma funcion de la libreria de calculo distribuido, y lo hacen de una forma
    que las versiones actuales ya no admiten: arboreto le pasa una lista vacia y
    SCENIC un generador, mientras que ahora se espera algo cuya longitud se pueda
    consultar de antemano. Cada caso rompe en un sitio distinto y con un mensaje
    que no menciona la causa.
    Se acepta lo que envian y se traduce a lo que la libreria espera. El ajuste
    vive en memoria y solo durante esta ejecucion: no se modifica ningun archivo
    instalado ni se baja la version de nada.
    """
    import dask.dataframe as dd
    import pandas as pd

    def envolver(original):
        if getattr(original, "_adaptada", False):
            return original

        def adaptada(dfs, *args, **kw):
            if not isinstance(dfs, (list, tuple)):
                dfs = list(dfs)
            if len(dfs) == 0:
                meta = kw.get("meta")
                vacio = meta.copy() if meta is not None else pd.DataFrame()
                return dd.from_pandas(vacio, npartitions=1)
            return original(dfs, *args, **kw)

        adaptada._adaptada = True
        return adaptada

    adaptados = []
    for nombre in ("arboreto.core", "pyscenic.prune"):
        try:
            mod = importlib.import_module(nombre)
        except Exception:
            continue
        if hasattr(mod, "from_delayed"):
            mod.from_delayed = envolver(mod.from_delayed)
            adaptados.append(nombre)
    return adaptados


def inferir_red(ex_matrix, tf_names, cfg, carpeta, n_est):
    """Reconstruye la red de regulacion."""
    import pandas as pd

    tfs_in = [t for t in tf_names if t in ex_matrix.columns]
    if not tfs_in:
        raise ValueError(
            "Ningun regulador de la lista aparece en la matriz de expresion.\n"
            "Sin ellos no hay red que inferir. Suele significar que los nombres de gen\n"
            "no coinciden (por ejemplo simbolos frente a identificadores Ensembl)."
        )
    print(f"Reguladores presentes en la matriz: {len(tfs_in)} de {len(tf_names)}")

    # Con menos de dos entradas en la lista de reguladores no hay nada que inferir.
    # Se corta aqui porque el error que llega despues, desde la libreria de calculo
    # distribuido, no menciona la causa.
    if len(tfs_in) < 2:
        raise ValueError(
            f"Solo hay {len(tfs_in)} regulador en la matriz ({tfs_in}).\n"
            "Hacen falta varios para tener algo que comparar."
        )

    t0 = time.time()

    if cfg.metodo_grn == "grnboost2":
        from arboreto.algo import grnboost2
        adaptar_from_delayed()

        def crear_cliente():
            """Devuelve (cliente, cluster), o (None, None) si no arranca aqui.

            El cliente se crea aqui a proposito: si se deja que la libreria monte el
            suyo, le pasa un argumento que la version actual ya no acepta.
            """
            try:
                from distributed import Client, LocalCluster
                n_w = max(1, min(4, os.cpu_count() or 1))
                cluster = LocalCluster(n_workers=n_w, threads_per_worker=1, processes=True)
                print(f"  calculo distribuido: {n_w} procesos")
                return Client(cluster), cluster
            except Exception as e:
                print(f"  no se pudo crear el cliente ({type(e).__name__}: {e})")
                print("  se usa el planificador interno de la libreria")
                return None, None

        cliente, cluster = crear_cliente()
        try:
            adjacencies = grnboost2(expression_data=ex_matrix, tf_names=tfs_in,
                                    client_or_address=cliente if cliente is not None else "local",
                                    seed=42, verbose=True)
        except Exception as e:
            raise RuntimeError(
                f"La reconstruccion de la red fallo: {type(e).__name__}: {e}\n\n"
                "No bajes la version de dask para arreglar esto: encadena mas fallos.\n"
                "Lee el error de arriba: si habla de los datos (pocos reguladores, matriz\n"
                "vacia), el problema esta en la preparacion de la matriz.\n"
                "Como ultimo recurso, --metodo-grn sklearn_aprox es una aproximacion no\n"
                "estandar, valida solo para demostrar el flujo."
            ) from e
        finally:
            if cliente is not None:
                cliente.close()
                cluster.close()
    else:
        # Alternativa que sigue la misma idea, medir que reguladores predicen mejor
        # la expresion de cada gen, pero no es el metodo oficial de SCENIC: le
        # faltan su criterio de parada y su esquema de muestreo. Sirve para
        # recorrer el flujo sin depender de esa libreria.
        from sklearn.ensemble import GradientBoostingRegressor
        print("Aproximacion sklearn (no estandar). Esto puede tardar bastante...")
        X_tfs = ex_matrix[tfs_in].values
        targets = [g for g in ex_matrix.columns if g not in tfs_in]
        registros = []
        for i, tg in enumerate(targets):
            y = ex_matrix[tg].values
            if y.std() == 0:
                continue
            gbm = GradientBoostingRegressor(n_estimators=n_est, max_depth=3, random_state=42)
            gbm.fit(X_tfs, y)
            for tf, imp in zip(tfs_in, gbm.feature_importances_):
                if imp > 0:
                    registros.append({"TF": tf, "target": tg, "importance": float(imp)})
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(targets)}")
        adjacencies = pd.DataFrame(registros)

    if adjacencies.empty:
        raise RuntimeError("La reconstruccion de la red no produjo ninguna conexion.")

    adjacencies = adjacencies.sort_values("importance", ascending=False).reset_index(drop=True)

    # La red se guarda entera. El paso siguiente aplica los umbrales propios de
    # SCENIC al formar los grupos de genes, y filtrar antes se apartaria del
    # procedimiento publicado.
    ruta = os.path.join(carpeta, "adjacencies.tsv")
    adjacencies.to_csv(ruta, sep="\t", index=False)

    por_tf = adjacencies.groupby("TF").size()
    n_posibles = ex_matrix.shape[1] - 1
    print(f"\nMetodo: {cfg.metodo_grn} | {time.time()-t0:.0f}s")
    print(f"Conexiones: {len(adjacencies):,}")
    print(f"Reguladores con genes asociados: {por_tf.size} | genes por regulador: "
          f"mediana {por_tf.median():.0f}, min {por_tf.min()}, max {por_tf.max()}")

    # Densidad de la red resultante. Si cada regulador queda conectado a casi todos
    # los genes, las puntuaciones no estan separando nada y el filtro del paso
    # siguiente recortaria sobre una lista sin orden util.
    if n_posibles > 0:
        densidad = por_tf.median() / n_posibles
        print(f"Densidad: {densidad:.1%}")
        if densidad > 0.9:
            print("\nAVISO: la red conecta casi todo con todo. Los grupos que se formen\n"
                  "despues seran practicamente arbitrarios.")
    return adjacencies


def filtrar_por_motivos(ex_matrix, cfg, carpeta):
    """Deja solo las conexiones con respaldo en la base de datos de SCENIC."""
    import loompy
    import numpy as np
    import pandas as pd

    if cfg.modo == "ejemplo":
        loom = os.path.join(carpeta, "expr_mat_tiny.loom")
    else:
        loom = os.path.join(carpeta, "expr_bcells.loom")
        if os.path.exists(loom):
            os.remove(loom)          # la libreria falla si el archivo ya existe
        loompy.create(loom, ex_matrix.T.values,
                      {"Gene": np.array(ex_matrix.columns)},
                      {"CellID": np.array(ex_matrix.index)})

    salida = os.path.join(carpeta, "regulons.csv")

    # Este paso se hace llamando a SCENIC dentro de este mismo proceso, y no
    # lanzando su programa de linea de comandos aparte.
    # El motivo es que ese programa reparte el trabajo con una libreria de 2020
    # que en Windows con Python 3.12 ya no arranca: los procesos hijos mueren al
    # instante y el principal se queda esperando unos resultados que no van a
    # llegar, sin dar ningun error. Se queda ahi indefinidamente, y ni siquiera
    # bajando a un solo trabajador lo evita. En Linux, que es donde corren los
    # cuadernos en la nube, los procesos se crean de otra forma y no ocurre.
    # Llamandolo desde aqui se puede usar el otro repartidor que trae SCENIC, el
    # mismo que ya funciono al reconstruir la red, y aplicar la adaptacion de
    # from_delayed, que en un proceso aparte no llegaria.
    adaptar_from_delayed()
    from ctxcore.rnkdb import FeatherRankingDatabase
    from pyscenic.cli.utils import save_enriched_motifs
    from pyscenic.prune import prune2df
    from pyscenic.utils import modules_from_adjacencies

    n_trabajadores = max(1, min(4, os.cpu_count() or 1))
    if _HEARTBEAT is not None:
        _HEARTBEAT.marcar("Filtrando por motivos de union", "entre 1 y 5 minutos")
    print(f"Filtrando la red con {n_trabajadores} trabajadores.")
    print("Este paso mostrara una barra de progreso real en cuanto arranque el calculo.\n")

    adjacencies = pd.read_csv(os.path.join(carpeta, "adjacencies.tsv"), sep="\t")
    modulos = list(modules_from_adjacencies(adjacencies, ex_matrix))
    print(f"  grupos de genes formados: {len(modulos):,}")
    if not modulos:
        raise RuntimeError(
            "No se formo ningun grupo de genes a partir de la red.\n"
            "Revisa el diagnostico del paso anterior: la red puede haber quedado\n"
            "demasiado pobre o demasiado densa."
        )

    db = FeatherRankingDatabase(
        fname=os.path.join(carpeta, "hg38.genes_vs_motifs.rankings.feather"),
        name="hg38")
    # Barra de progreso real de dask: cuenta tareas del calculo repartido, no una
    # simulacion. Es el paso que antes se colgaba en silencio, asi que aqui es
    # donde mas hace falta ver movimiento de verdad.
    from dask.diagnostics import ProgressBar
    with ProgressBar():
        df_motifs = prune2df([db], modulos, os.path.join(carpeta, "motifs.tbl"),
                             client_or_address="dask_multiprocessing",
                             num_workers=n_trabajadores)
    save_enriched_motifs(df_motifs, salida)

    # Que el paso termine no significa que haya encontrado algo: queda contar
    # cuantas filas trae el archivo de salida.
    try:
        n_regulons = len(pd.read_csv(salida, index_col=[0, 1], header=[0, 1]))
    except Exception:
        n_regulons = 0

    print(f"\nRegulones validados: {n_regulons}")
    if n_regulons == 0:
        print(
            "\nEl filtrado termino sin error pero no valido ningun regulon. Causas tipicas:\n"
            "  - dataset demasiado pequeno (los datos de ejemplo tienen 500 genes)\n"
            "  - reguladores que no figuran en la base de datos\n"
            "  - nombres de gen que no casan con la base\n"
            "  - la red del paso anterior era demasiado densa o demasiado pobre"
        )
    return n_regulons


def puntuar_celulas(ex_matrix, adjacencies, n_regulons, cfg, carpeta, fig):
    """Calcula cuan activo esta cada regulon en cada celula."""
    import pandas as pd
    from ctxcore.genesig import GeneSignature
    from pyscenic.aucell import aucell as pyscenic_aucell

    signatures = []
    # Los datos de ejemplo tienen 500 genes, demasiado pocos para que ningun
    # regulon supere la validacion. Es una consecuencia de su tamano, no un fallo,
    # asi que ese modo continua siempre: existe para mostrar el procedimiento
    # entero. Con el dataset real la exigencia se mantiene.
    continuar = cfg.permitir_sin_validar or (cfg.modo == "ejemplo")

    if n_regulons > 0:
        validados = True
        # El archivo de salida se lee con las funciones del propio SCENIC en vez de
        # interpretarlo a mano. Su formato tiene dos particularidades que solo ellas
        # resuelven. La primera es como quedan escritos los numeros dentro de las
        # celdas de texto. La segunda es que una misma clave puede ocupar varias
        # filas, y esas funciones las agrupan en una sola entrada uniendo sus genes.
        from pyscenic.transform import df2regulons
        from pyscenic.utils import load_motifs

        df_motifs = load_motifs(os.path.join(carpeta, "regulons.csv"))
        signatures = df2regulons(df_motifs)
        if not signatures:
            raise RuntimeError(
                f"El filtrado reporto {n_regulons} regulones pero no se pudo construir\n"
                "ninguno. Revisa el archivo regulons.csv."
            )
    elif not continuar:
        # Parada deliberada. Sin regulones validados, continuar produciria cifras
        # con aspecto de resultado que no significan nada.
        raise RuntimeError(
            "El filtrado no valido ningun regulon, asi que no hay nada que puntuar.\n"
            "\n"
            "El analisis para aqui a proposito. Construir los regulones sin ese filtro\n"
            "elimina justo el paso que diferencia a SCENIC de una red de coexpresion:\n"
            "el resultado son puntuaciones casi identicas entre si que luego se\n"
            "grafican como si significaran algo.\n"
            "\n"
            "Que hacer:\n"
            "  - Revisa el diagnostico de la red: que los nombres de gen casen con la\n"
            "    base de datos, y que la densidad no sea degenerada.\n"
            "  - Solo para demostrar el flujo: --permitir-sin-validar. Las salidas\n"
            "    quedaran marcadas como no validadas."
        )
    else:
        validados = False
        causa = ("son los datos de ejemplo, que no pueden validar por diseno"
                 if cfg.modo == "ejemplo" and not cfg.permitir_sin_validar
                 else "se pidio expresamente continuar sin validacion")
        print("=" * 70)
        print(f"ATENCION: regulones sin validar. Causa: {causa}")
        print("Esto no es SCENIC completo. Los resultados no son interpretables como")
        print("actividad regulatoria real y no deben entregarse como tal.")
        print("=" * 70)
        for tf, grp in adjacencies.groupby("TF"):
            signatures.append(GeneSignature(name=f"{tf}(+)",
                                            gene2weight=dict(zip(grp["target"], grp["importance"]))))

    sufijo = "" if validados else "_SIN_VALIDAR"
    tam = [len(s.genes) for s in signatures]
    print(f"\nRegulones: {len(signatures)} | genes por regulon: "
          f"mediana {int(pd.Series(tam).median())}, min {min(tam)}, max {max(tam)}")

    auc_matrix = pyscenic_aucell(ex_matrix, signatures, num_workers=1)
    print("Matriz de actividad:", auc_matrix.shape)

    # Estas puntuaciones miden cuan activo esta cada regulon en cada celula. Si
    # todas salen parecidas, los regulones no estan distinguiendo unas celulas de
    # otras y el resultado no dice nada.
    rango = auc_matrix.values.max() - auc_matrix.values.min()
    print(f"Actividad: min {auc_matrix.values.min():.4f} | max {auc_matrix.values.max():.4f} "
          f"| rango {rango:.4f}")
    if rango < 0.05:
        print("\nAVISO: el rango es minusculo; los regulones apenas separan unas celulas\n"
              "de otras. Cualquier agrupamiento calculado sobre esto sera ruido.")

    # --- Grafica de regulones ---
    import scanpy as sc
    asc = sc.AnnData(X=auc_matrix.values,
                     obs=pd.DataFrame(index=auc_matrix.index.astype(str)),
                     var=pd.DataFrame(index=auc_matrix.columns.astype(str)))
    sc.pp.neighbors(asc, random_state=42)
    sc.tl.umap(asc, random_state=42)
    sc.tl.leiden(asc, flavor="igraph", n_iterations=2, random_state=42)

    # El titulo indica si los regulones estan validados. Las figuras acaban en
    # presentaciones e informes separadas del mensaje que las acompanaba en
    # pantalla, y esa advertencia tiene que viajar con la imagen.
    etiqueta = "regulones validados" if validados else "SIN VALIDAR: no interpretable"
    sc.pl.umap(asc, color="leiden",
               title=f"Actividad de regulones ({cfg.modo}): {etiqueta}", show=False)
    fig.guardar("scenic_regulones")

    salida = os.path.join(cfg.dir_salida, f"scenic_auc_{cfg.modo}{sufijo}.csv")
    auc_matrix.to_csv(salida)
    print("Guardado:", salida)
    if not validados:
        print("\nRecuerda: el sufijo _SIN_VALIDAR indica que estas puntuaciones no pasaron\n"
              "el filtro. No las entregues como resultado de SCENIC.")
    return auc_matrix, validados


def resumen_final(cfg, adjacencies, n_regulons, auc_matrix, validados):
    """Deja constancia de con que datos, que metodo y que versiones se ejecuto."""

    def ver(mod):
        try:
            return importlib.import_module(mod).__version__
        except Exception:
            return "no disponible"

    apto = validados and cfg.metodo_grn == "grnboost2"

    print()
    print("=" * 62)
    print("  RESUMEN DEL ANALISIS")
    print("=" * 62)
    print(f"  Fecha              : {datetime.datetime.now():%Y-%m-%d %H:%M}")
    print(f"  Datos              : {cfg.modo}")
    print(f"  Reconstruccion red : {cfg.metodo_grn}"
          f"{'  (metodo estandar de SCENIC)' if cfg.metodo_grn == 'grnboost2' else '  (APROXIMACION NO ESTANDAR)'}")
    print(f"  Conexiones         : {len(adjacencies):,}")
    print(f"  Regulones          : {n_regulons}")
    print(f"  Validados          : {'SI' if validados else 'NO'}")
    print(f"  Matriz de actividad: {auc_matrix.shape[0]:,} celulas x {auc_matrix.shape[1]} regulones")
    print(f"  Rango de actividad : {auc_matrix.values.min():.4f} - {auc_matrix.values.max():.4f}")
    print("-" * 62)
    print(f"  Python {sys.version.split()[0]} | numpy {ver('numpy')} | "
          f"scanpy {ver('scanpy')} | pyscenic {ver('pyscenic')}")
    print("=" * 62)
    if apto:
        print("  Metodo estandar + regulones validados.")
    else:
        motivos = []
        if cfg.metodo_grn != "grnboost2":
            motivos.append("la red no se reconstruyo con el metodo estandar")
        if not validados:
            motivos.append("los regulones no pasaron el filtro")
        print("  NO APTO PARA ENTREGA como resultado de SCENIC, porque "
              + " y ".join(motivos) + ".")
        if cfg.modo == "ejemplo":
            print("  Esto es lo esperado con los datos de ejemplo: 500 genes nunca pueden")
            print("  validar regulones. Sirve para demostrar el mecanismo, no como")
            print("  resultado. Usa --modo zenodo para un analisis real.")
    print("=" * 62)


# ------------------------------------------------------------------
#   Configuracion y arranque
# ------------------------------------------------------------------

class Config:
    """Junta los valores del panel de control con los de la linea de comandos."""

    def __init__(self, args):
        self.modo = args.modo or MODO
        self.saltar_graficas = SALTAR_GRAFICAS if args.saltar_graficas is None else args.saltar_graficas
        self.n_cells_max = N_CELLS_MAX if args.n_cells_max is None else args.n_cells_max
        self.rds_debug = RDS_DEBUG if args.rds_debug is None else args.rds_debug
        self.n_cells_debug = N_CELLS_DEBUG
        self.rds_debug_file = RDS_DEBUG_FILE
        self.scenic_downsample = SCENIC_DOWNSAMPLE
        self.scenic_n_cells = args.scenic_n_cells or SCENIC_N_CELLS
        self.scenic_n_genes = SCENIC_N_GENES
        self.metodo_grn = args.metodo_grn or METODO_GRN
        self.permitir_sin_validar = (PERMITIR_REGULONES_SIN_VALIDAR
                                     if args.permitir_sin_validar is None
                                     else args.permitir_sin_validar)
        self.col_cond = COL_COND
        self.col_time = COL_TIME
        self.col_ctype = args.col_ctype or COL_CTYPE
        self.pre_label = PRE_LABEL
        self.post_label = POST_LABEL
        self.verificar_checksums = VERIFICAR_CHECKSUMS
        self.root_cluster = ROOT_CLUSTER
        self.cluster_key = CLUSTER_KEY
        self.dir_salida = os.path.abspath(args.salida or DIR_SALIDA)
        self.dir_datos = os.path.abspath(args.datos or DIR_DATOS)
        self.dir_trabajo = os.path.join(self.dir_datos, "conversion")

        if self.modo not in ("ejemplo", "zenodo"):
            raise SystemExit(f"MODO invalido: {self.modo!r}. Usa 'ejemplo' o 'zenodo'.")
        if self.metodo_grn not in ("grnboost2", "sklearn_aprox"):
            raise SystemExit(f"METODO_GRN invalido: {self.metodo_grn!r}")

        for d in (self.dir_salida, self.dir_datos, self.dir_trabajo):
            os.makedirs(d, exist_ok=True)

    def mostrar(self):
        print(f"  Datos            : {self.modo}")
        print(f"  Metodo de red    : {self.metodo_grn}")
        print(f"  Celulas maximas  : {'todas' if self.n_cells_max == 0 else f'{self.n_cells_max:,}'}")
        print(f"  Graficas         : {'omitidas' if self.saltar_graficas else 'si'}")
        print(f"  Tipo celular     : {self.col_ctype}")
        print(f"  Resultados en    : {self.dir_salida}")
        print(f"  Descargas en     : {self.dir_datos}")
        if self.saltar_graficas:
            print("\n  AVISO: con las graficas omitidas la corrida no sirve como entregable.")
        if self.rds_debug:
            print("\n  AVISO: modo de pruebas activado. Los resultados no son un entregable.")
        if self.metodo_grn == "sklearn_aprox":
            print("\n  AVISO: el metodo elegido no es el estandar de SCENIC. Los resultados")
            print("         no son comparables con los publicados.")


def leer_argumentos():
    p = argparse.ArgumentParser(
        description="Pipeline scRNA-seq: scanpy + SCENIC (version de escritorio).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Sin argumentos usa los valores del PANEL DE CONTROL, dentro del archivo.")
    p.add_argument("--modo", choices=["ejemplo", "zenodo"],
                   help="'ejemplo' corre en minutos; 'zenodo' es el estudio real.")
    p.add_argument("--n-cells-max", type=int, dest="n_cells_max",
                   help="Celulas a analizar. 0 = todas (pide mucha memoria).")
    p.add_argument("--salida", help="Carpeta donde dejar tablas y graficas.")
    p.add_argument("--datos", help="Carpeta donde guardar las descargas.")
    p.add_argument("--col-ctype", dest="col_ctype", help="Columna con el tipo celular.")
    p.add_argument("--metodo-grn", dest="metodo_grn", choices=["grnboost2", "sklearn_aprox"],
                   help="Metodo para reconstruir la red.")
    p.add_argument("--scenic-n-cells", type=int, dest="scenic_n_cells",
                   help="Celulas que entran a SCENIC.")
    p.add_argument("--saltar-graficas", action="store_true", default=None,
                   dest="saltar_graficas", help="Omite las 12 graficas y va directo a SCENIC.")
    p.add_argument("--rds-debug", action="store_true", default=None, dest="rds_debug",
                   help="Trabaja sobre una copia pequena, solo para probar.")
    p.add_argument("--permitir-sin-validar", action="store_true", default=None,
                   dest="permitir_sin_validar",
                   help="Continua aunque no se valide ningun regulon (marca las salidas).")
    return p.parse_args()


def main():
    global _HEARTBEAT
    cfg = Config(leer_argumentos())
    _HEARTBEAT = Heartbeat()
    _HEARTBEAT.iniciar()

    titulo("PIPELINE scRNA-seq: scanpy + SCENIC")
    cfg.mostrar()
    print("\n  Cada " + str(_HEARTBEAT.intervalo) + " segundos se imprime una linea de "
          "[latido] con el tiempo\n  transcurrido y el uso de CPU. Sirve para distinguir un paso "
          "lento de uno\n  colgado: si el CPU se queda cerca de 0% durante mucho rato, avisa.")

    try:
        # Las graficas se escriben a disco, asi que matplotlib no necesita abrir
        # ninguna ventana. Sin esto, en algunos equipos el script se queda esperando.
        import matplotlib
        matplotlib.use("Agg")
        warnings.filterwarnings("ignore", category=FutureWarning)
        warnings.filterwarnings("ignore", category=UserWarning, module="anndata")

        paso("Comprobando el entorno")
        import psutil
        ram_gb = psutil.virtual_memory().total / 1e9
        print(f"  RAM: {ram_gb:.1f} GB | CPUs: {os.cpu_count()}")
        if cfg.modo == "zenodo" and cfg.n_cells_max == 0 and ram_gb < 20:
            print(f"\n  AVISO: pediste todas las celulas y este equipo tiene {ram_gb:.1f} GB.")
            print("         El paso de conversion se detendra si no alcanza.")
        preparar_alias_numpy()

        import scanpy as sc
        sc.settings.verbosity = 1
        sc.settings.set_figure_params(dpi=100, facecolor="white")
        fig = Figuras(os.path.join(cfg.dir_salida, "graficas"))

        # --- Carga de datos ---
        titulo("1. DATOS")
        if cfg.modo == "ejemplo":
            adata, col_cond, col_time, col_ctype = cargar_ejemplo()
        else:
            rscript = preparar_r(cfg)
            rds = descargar_rds(cfg)
            convertir_rds(cfg, rscript, rds)
            adata = reconstruir_adata(cfg)
            col_cond, col_time, col_ctype = resolver_columnas(adata, cfg)
            print()
            print(adata)
            destino_h5ad = os.path.join(cfg.dir_salida, "datos_convertidos.h5ad")
            adata.write_h5ad(destino_h5ad)
            print("Guardado:", destino_h5ad)
            mtx = os.path.join(cfg.dir_trabajo, "counts.mtx")
            if os.path.exists(mtx):
                os.remove(mtx)
            gc.collect()

        markers = marcadores_por_modo(adata, cfg.modo)
        print("Marcadores:", list(markers.keys()))

        # --- Analisis con scanpy ---
        titulo("2. CONTROL DE CALIDAD, AGRUPAMIENTO Y GRAFICAS",
              "unos minutos con 60.000 celulas; mas si se analizan todas")
        analisis_scanpy(adata, cfg, fig, col_cond, col_ctype, markers)

        # --- Figuras del estudio ---
        adata_B = None
        if cfg.modo == "zenodo" and not cfg.saltar_graficas:
            titulo("3. FIGURAS DEL ESTUDIO", "menos de un minuto")
            adata_B = figuras_estudio(adata, cfg, fig, col_time, col_ctype)

        # --- SCENIC ---
        titulo("4. SCENIC: REDES DE REGULACION")
        carpeta = descargar_bases_scenic(cfg)
        paso("Preparando la matriz", "menos de un minuto")
        ex_matrix, tf_names, n_est = preparar_matriz_scenic(adata, adata_B, cfg, carpeta, col_ctype)
        paso("Reconstruyendo la red", "1 a 3 minutos")
        adjacencies = inferir_red(ex_matrix, tf_names, cfg, carpeta, n_est)
        n_regulons = filtrar_por_motivos(ex_matrix, cfg, carpeta)
        paso("Puntuando celula a celula", "menos de un minuto")
        auc_matrix, validados = puntuar_celulas(ex_matrix, adjacencies, n_regulons, cfg, carpeta, fig)

        # --- Resumen ---
        titulo("5. RESUMEN")
        resumen_final(cfg, adjacencies, n_regulons, auc_matrix, validados)

        print()
        print(f"Se generaron {fig.n} graficas en: {os.path.join(cfg.dir_salida, 'graficas')}")
        print(f"Tablas en: {cfg.dir_salida}")
        print()
    finally:
        # Se detiene aqui tambien si algo fallo arriba, para que no siga
        # imprimiendo latidos de un analisis que ya termino (bien o mal).
        _HEARTBEAT.detener()


if __name__ == "__main__":
    # El calculo distribuido de SCENIC lanza procesos nuevos, y en Windows cada uno
    # vuelve a importar este archivo. Sin esta guarda, cada proceso reiniciaria el
    # analisis entero.
    main()
