@echo off
REM ============================================================
REM   Lanza el pipeline scanpy + SCENIC
REM ============================================================
REM Usa el Python del entorno que preparo instalar_dependencias.bat, no el del
REM sistema: es el unico que tiene las librerias instaladas.
REM
REM Todo lo que se escriba despues del nombre de este archivo se le pasa tal
REM cual al analisis. Por ejemplo:
REM     ejecutar.bat --modo ejemplo
REM     ejecutar.bat --n-cells-max 30000
REM     ejecutar.bat --salida D:\resultados_lupus
REM ============================================================

setlocal
cd /d "%~dp0"

if not exist entorno\Scripts\python.exe (
    echo.
    echo ERROR: falta el entorno con las librerias.
    echo.
    echo   Ejecuta primero:  instalar_dependencias.bat
    echo.
    pause
    exit /b 1
)

entorno\Scripts\python.exe pipeline_scanpy_scenic.py %*
set CODIGO=%errorlevel%

echo.
if %CODIGO% neq 0 (
    echo ============================================================
    echo   El analisis se detuvo (codigo %CODIGO%^)
    echo ============================================================
    echo   Lee el mensaje de arriba: esta escrito para explicar la causa
    echo   y que hacer, no solo para senalar donde fallo.
) else (
    echo ============================================================
    echo   Analisis terminado
    echo ============================================================
)
echo.
pause
exit /b %CODIGO%
