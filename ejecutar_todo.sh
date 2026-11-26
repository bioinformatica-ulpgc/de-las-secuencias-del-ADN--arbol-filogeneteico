#!/usr/bin/env bash
# Ejecuta la práctica completa de principio a fin.
# Uso: bash ejecutar_todo.sh   (con el entorno virtual activado o "python" apuntando a uno con requirements.txt)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python}

echo "== 1. Código proporcionado: GenBank -> COX1 -> alineamiento -> distancias p -> UPGMA =="
$PY proporcionado/procesar_datos.py > /dev/null
$PY proporcionado/generar_animacion.py

echo "== 2. Neighbor-Joining propio =="
$PY nj.py

echo "== 3. Pruebas de validación =="
$PY -m unittest -v

echo "== 4. Comparación NJ propio vs Biopython y UPGMA vs NJ =="
$PY comparar.py > /dev/null && echo "resultados/comparacion.md y resultados/upgma-vs-nj.png generados"

if $PY -c "import typst" 2>/dev/null; then
  echo "== 5. Informe PDF =="
  $PY informe/construir_informe.py
fi
