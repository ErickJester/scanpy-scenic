@echo off
REM ============================================================
REM   Instalador de dependencias del pipeline scanpy + SCENIC
REM ============================================================
REM Prepara un entorno de Python aislado en la carpeta "entorno" e instala
REM dentro todo lo que el analisis necesita. Se usa un entorno propio para no
REM tocar el Python del sistema: las versiones que pide SCENIC son especificas y
REM podrian romper otros programas.
REM
REM Ejecutar una sola vez, con doble clic o desde la consola.
REM ============================================================

setlocal
cd /d "%~dp0"

echo.
echo ============================================================
echo   Instalacion de dependencias
echo ============================================================
echo.

REM ---- 1. Comprobar que Python esta disponible ---------------
where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: no se encontro Python en este equipo.
    echo.
    echo   Instalalo desde https://www.python.org/downloads/
    echo   IMPORTANTE: marca la casilla "Add Python to PATH" durante la
    echo   instalacion, o este script no podra encontrarlo.
    echo.
    pause
    exit /b 1
)

for /f "tokens=2" %%v in ('python --version 2^>^&1') do set PYVER=%%v
echo Python detectado: %PYVER%
echo.

REM SCENIC y sus dependencias no tienen version compatible con Python 3.13 en
REM el momento de escribir esto. Se avisa en vez de fallar a mitad de la
REM instalacion con un error dificil de interpretar.
echo %PYVER% | findstr /b "3.13 3.14" >nul
if not errorlevel 1 (
    echo AVISO: la version de Python detectada es muy reciente y algunas
    echo        librerias de SCENIC todavia no publican una compatible.
    echo        Si la instalacion falla, instala Python 3.11 o 3.12 y vuelve
    echo        a ejecutar este archivo.
    echo.
)

REM ---- 2. Crear el entorno aislado ---------------------------
if exist entorno\Scripts\python.exe (
    echo El entorno ya existe: se reutiliza.
) else (
    echo Creando el entorno en la carpeta "entorno"...
    python -m venv entorno
    if errorlevel 1 (
        echo.
        echo ERROR: no se pudo crear el entorno.
        pause
        exit /b 1
    )
)
echo.

set PY=entorno\Scripts\python.exe

REM ---- 3. Instalar las librerias -----------------------------
echo Actualizando pip...
"%PY%" -m pip install --upgrade pip --quiet

echo.
echo Instalando la base del entorno...
REM setuptools se fija por debajo de la 81 porque a partir de ahi dejo de incluir
REM pkg_resources, que las librerias de SCENIC todavia usan al arrancar. Sin esta
REM linea la instalacion termina sin errores y el analisis falla despues con un
REM "No module named 'pkg_resources'" que no dice de donde viene. Ademas, desde
REM Python 3.12 los entornos nuevos ya no traen setuptools de serie, asi que hay
REM que pedirlo explicitamente.
"%PY%" -m pip install "setuptools<81" wheel
if errorlevel 1 goto error

echo.
echo Instalando scanpy y el agrupamiento de celulas...
REM Las versiones estan acotadas porque las librerias de analisis unicelular y
REM las de calculo distribuido avanzan a ritmos distintos y no toda combinacion
REM funciona. numpy queda entre 2.0 y 2.1: scanpy necesita la 2, y numba, que
REM scanpy usa por dentro, todavia no admite la 2.1. Para dask se pide una
REM version minima y ninguna maxima.
REM scikit-image lo pide la deteccion de dobletes del modo "ejemplo". En Colab
REM viene preinstalado, asi que alli nunca hace falta nombrarlo; en un equipo
REM limpio no esta, y sin el ese modo se detiene a mitad.
"%PY%" -m pip install "numpy>=2,<2.1" "dask>=2024.1" "distributed>=2024.1" ^
    scanpy leidenalg igraph scikit-misc scikit-image
if errorlevel 1 goto error

echo.
echo Instalando SCENIC...
"%PY%" -m pip install pyscenic==0.12.1
if errorlevel 1 goto error

echo.
echo Instalando el resto...
REM pooch descarga los datos de ejemplo; psutil mide la memoria disponible;
REM matplotlib dibuja las graficas.
"%PY%" -m pip install pooch psutil matplotlib
if errorlevel 1 goto error

REM ---- 4. Comprobar que todo carga ---------------------------
echo.
echo Comprobando que las librerias cargan...
REM Se comprueba que cargan de verdad, no solo que pip dijera que si. Una
REM instalacion puede terminar sin errores y aun asi dejar el entorno en un
REM estado donde alguna libreria no arranca, y entonces el fallo apareceria mas
REM tarde, ya empezado el analisis.
"%PY%" -W ignore -c "import scanpy, pyscenic, ctxcore, loompy, arboreto, dask, distributed, pooch, psutil, matplotlib; print('  todas las librerias cargan correctamente')"
if errorlevel 1 (
    echo.
    echo ------------------------------------------------------------
    echo   AVISO: la instalacion termino pero alguna libreria no carga
    echo ------------------------------------------------------------
    echo   Lee el error de arriba. Si dice "No module named 'pkg_resources'",
    echo   ejecuta esto y vuelve a intentarlo:
    echo.
    echo       entorno\Scripts\python.exe -m pip install "setuptools^<81"
    echo.
    echo   Para cualquier otro error, borra la carpeta "entorno" y ejecuta
    echo   este archivo de nuevo.
    echo ------------------------------------------------------------
    pause
    exit /b 1
)

REM ---- 5. Avisar sobre R -------------------------------------
echo.
where Rscript >nul 2>nul
if errorlevel 1 (
    echo ------------------------------------------------------------
    echo   FALTA R  (solo hace falta para el modo "zenodo"^)
    echo ------------------------------------------------------------
    echo   El archivo del estudio viene en formato de Seurat, que es una
    echo   herramienta de R, asi que hace falta R para leerlo.
    echo.
    echo   Descargalo de: https://cran.r-project.org/bin/windows/base/
    echo   Los paquetes de R que falten se instalan solos la primera vez.
    echo.
    echo   Si solo quieres probar el flujo, el modo "ejemplo" no necesita R:
    echo       ejecutar.bat --modo ejemplo
    echo ------------------------------------------------------------
) else (
    echo R detectado correctamente.
)

echo.
echo ============================================================
echo   Instalacion terminada
echo ============================================================
echo.
echo   Para lanzar el analisis:   ejecutar.bat
echo.
pause
exit /b 0

:error
echo.
echo ============================================================
echo   ERROR durante la instalacion
echo ============================================================
echo.
echo   Revisa el mensaje de pip que aparece mas arriba. Lo mas habitual:
echo     - sin conexion a internet
echo     - version de Python demasiado reciente (usa 3.11 o 3.12^)
echo.
pause
exit /b 1
