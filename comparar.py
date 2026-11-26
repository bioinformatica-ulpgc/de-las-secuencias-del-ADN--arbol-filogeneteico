"""
comparar.py - Verificación y comparación de árboles (CÓDIGO PROPIO).

1. NJ propio  vs  NJ de referencia (Bio.Phylo.TreeConstruction.DistanceTreeConstructor.nj)
2. UPGMA (proporcionado) vs NJ propio, con las mismas distancias p
3. Ejemplo sintético aditivo sin reloj molecular (UPGMA falla, NJ acierta)

Genera resultados/comparacion.md y resultados/upgma-vs-nj.png.
Requiere haber ejecutado antes proporcionado/procesar_datos.py y nj.py.
"""
from __future__ import annotations

import io
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from arboles import (biparticiones, describir_split, distancias_patristicas, enraizar_punto_medio, leer_newick,
                     no_triviales, robinson_foulds, upgma)
from nj import ABREV, leer_matriz, neighbor_joining

RAIZ = Path(__file__).resolve().parent
RES = RAIZ / "resultados"
OUTGROUP = "Hylobates_lar"  # el gibón (Hylobatidae) es externo a los grandes simios (Hominidae)


def ab(x: str) -> str:
    return ABREV.get(x.replace("_", " "), x)


def nj_biopython(etiquetas, D):
    """Devuelve (árbol, Newick) del NJ de Biopython. Las comparaciones usan el objeto
    árbol, porque el nº de decimales del Newick escrito depende de la versión."""
    dm = DistanceMatrix(etiquetas, [[float(D[i, j]) for j in range(i + 1)] for i in range(len(etiquetas))])
    arbol = DistanceTreeConstructor().nj(dm)
    h = io.StringIO()
    Phylo.write(arbol, h, "newick")
    return arbol, h.getvalue().strip()


def tabla_splits(arboles: dict[str, str], hojas: set[str]) -> list[str]:
    """Tabla markdown: cada bipartición (interna y terminal) y su longitud en cada árbol."""
    sp = {k: biparticiones(v) for k, v in arboles.items()}
    todas = sorted(set().union(*map(set, sp.values())), key=lambda s: (len(s), sorted(s)))
    filas = ["| Bipartición (arista) | Tipo | " + " | ".join(arboles) + " |",
             "|---|---|" + "---|" * len(arboles)]
    for s in todas:
        tipo = "terminal" if len(s) == 1 or len(s) == len(hojas) - 1 else "interna"
        celdas = [f"{sp[k][s]:.6f}" if s in sp[k] else "—" for k in arboles]
        filas.append(f"| {describir_split(s, hojas, ab)} | {tipo} | " + " | ".join(celdas) + " |")
    return filas


def profundidades(newick: str, punto_medio: bool = False) -> dict[str, float]:
    arbol = enraizar_punto_medio(newick) if punto_medio else leer_newick(newick)
    d = arbol.depths()
    return {t.name: d[t] for t in arbol.get_terminals()}


def grupos_enraizados(newick: str, outgroup: str | None = None, punto_medio: bool = False) -> list[str]:
    arbol = leer_newick(newick)
    if outgroup:
        arbol.root_with_outgroup(outgroup)
    if punto_medio:
        arbol = enraizar_punto_medio(newick)
    out = []
    for c in arbol.get_nonterminals():
        hs = sorted(ab(t.name) for t in c.get_terminals())
        if 1 < len(hs) < len(arbol.get_terminals()):
            out.append("{" + ", ".join(hs) + "}")
    return sorted(out, key=len)


def cuatro_puntos(D: np.ndarray) -> float:
    """Máxima desviación de la condición de los cuatro puntos (0 si la matriz es aditiva)."""
    peor = 0.0
    for i, j, k, l in combinations(range(len(D)), 4):
        s = sorted([D[i, j] + D[k, l], D[i, k] + D[j, l], D[i, l] + D[j, k]])
        peor = max(peor, s[2] - s[1])
    return peor


def ultrametricidad(D: np.ndarray) -> float:
    """Máxima desviación de la condición ultramétrica d(i,k) ≤ max(d(i,j), d(j,k)):
    en cada trío los dos valores mayores deberían ser iguales."""
    peor = 0.0
    for i, j, k in combinations(range(len(D)), 3):
        s = sorted([D[i, j], D[i, k], D[j, k]])
        peor = max(peor, s[2] - s[1])
    return peor


def dibujar_comparacion(nw_upgma: str, nw_nj: str, ruta: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor="white")
    t_up = leer_newick(nw_upgma)
    t_nj = enraizar_punto_medio(nw_nj)  # cae en la rama del gibón: equivale a usarlo de outgroup
    t_nj.ladderize(reverse=True)
    t_up.ladderize(reverse=True)
    for ax, t, tit in [(axes[0], t_up, "UPGMA (proporcionado): enraizado, ultramétrico"),
                       (axes[1], t_nj, "NJ propio enraizado (punto medio = rama del gibón)")]:
        Phylo.draw(t, axes=ax, do_show=False,
                   label_func=lambda c: ab(c.name) if c.is_terminal() else None,
                   branch_labels=lambda c: f"{c.branch_length:.4f}" if c.branch_length else None)
        ax.set_title(tit, fontsize=11)
        ax.set_xlabel("distancia p (diferencias por sitio)")
        ax.set_ylabel("")
        ax.set_xlim(0, 0.095)
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def main() -> None:
    RES.mkdir(exist_ok=True)
    nombres, D = leer_matriz(RAIZ / "proporcionado" / "distancias.csv")
    etiquetas = [n.replace(" ", "_") for n in nombres]
    hojas = set(etiquetas)

    nw_nj = neighbor_joining(nombres, D).newick
    t_bio, nw_bio = nj_biopython(etiquetas, D)
    nw_up = (RAIZ / "proporcionado" / "arbol-upgma.nwk").read_text().strip()
    (RES / "arbol-nj-biopython.nwk").write_text(nw_bio + "\n")

    md = ["# Verificación y comparación de árboles", "",
          "Generado automáticamente por `comparar.py` a partir de `proporcionado/distancias.csv`.", ""]

    # ---------------- 1. NJ propio vs Biopython ----------------
    s_nj, s_bio = biparticiones(nw_nj), biparticiones(t_bio)
    dif = max(abs(s_nj[s] - s_bio[s]) for s in s_nj) if set(s_nj) == set(s_bio) else float("nan")
    md += ["## 1. NJ propio frente a la implementación de referencia (Biopython)", "",
           "```", f"NJ propio   : {nw_nj}", f"NJ Biopython: {nw_bio}", "```", "",
           "Los dos Newick se escriben distinto (orden de los hijos, nombres `Inner` de los nodos "
           "internos), pero describen el mismo árbol no enraizado si tienen las mismas biparticiones:", "",
           *tabla_splits({"NJ propio": nw_nj, "NJ Biopython": t_bio}, hojas), "",
           f"- Distancia de Robinson-Foulds: **{robinson_foulds(nw_nj, t_bio)}** (0 = misma topología)",
           f"- Misma colección de aristas: **{set(s_nj) == set(s_bio)}**",
           f"- Máxima diferencia en longitud de rama: **{dif:.2e}**", ""]

    # ---------------- 2. UPGMA vs NJ ----------------
    n = len(nombres)
    s_up = biparticiones(nw_up)
    P_nj = distancias_patristicas(nw_nj, etiquetas)
    P_up = distancias_patristicas(nw_up, etiquetas)
    iu = np.triu_indices(n, 1)
    prof_up = profundidades(nw_up)
    prof_nj = profundidades(nw_nj, punto_medio=True)
    t_mid = enraizar_punto_medio(nw_nj)
    rama_gibon = next(c for c in t_mid.root.clades if c.is_terminal())
    md += ["## 2. UPGMA frente a NJ (mismas distancias p)", "",
           "```", f"UPGMA: {nw_up}", f"NJ   : {nw_nj}", "```", "",
           "### 2.1 Aristas (biparticiones) y longitudes", "",
           "UPGMA se trata aquí como no enraizado (sus dos ramas de la raíz se suman en una sola arista).", "",
           *tabla_splits({"UPGMA": nw_up, "NJ": nw_nj}, hojas), "",
           f"- Distancia de Robinson-Foulds UPGMA–NJ (no enraizados): **{robinson_foulds(nw_up, nw_nj)}** "
           f"de {2 * (n - 3)} posibles",
           "- Biparticiones internas solo en UPGMA: " + (", ".join(
               describir_split(s, hojas, ab) for s in no_triviales(s_up, n) - no_triviales(s_nj, n)) or "ninguna"),
           "- Biparticiones internas solo en NJ: " + (", ".join(
               describir_split(s, hojas, ab) for s in no_triviales(s_nj, n) - no_triviales(s_up, n)) or "ninguna"),
           "",
           "### 2.2 Agrupaciones como árboles enraizados", "",
           "| Árbol | Clados (grupos con ≥ 2 especies, sin contar la raíz) |", "|---|---|",
           f"| UPGMA (raíz implícita del método) | {'; '.join(grupos_enraizados(nw_up))} |",
           f"| NJ enraizado con el gibón como outgroup | {'; '.join(grupos_enraizados(nw_nj, OUTGROUP))} |",
           f"| NJ enraizado en el punto medio | {'; '.join(grupos_enraizados(nw_nj, None, True))} |",
           "",
           "### 2.3 Distancia de cada hoja a la raíz", "",
           "| Especie | UPGMA | NJ (raíz en el punto medio) |", "|---|---|---|"]
    for e in etiquetas:
        md.append(f"| {ab(e)} | {prof_up[e]:.6f} | {prof_nj[e]:.6f} |")
    md += ["",
           f"UPGMA: todas iguales (rango {max(prof_up.values()) - min(prof_up.values()):.1e}) → ultramétrico. "
           f"NJ: profundidades desiguales (rango {max(prof_nj.values()) - min(prof_nj.values()):.4f}) → "
           "las tasas de cambio no son idénticas en todos los linajes. El punto medio del camino más largo "
           f"cae sobre la rama del gibón ({ab(rama_gibon.name)} a {rama_gibon.branch_length:.4f} de la raíz), "
           "así que el enraizamiento por punto medio coincide con el del outgroup.", "",
           "### 2.4 Ajuste a la matriz de distancias", "",
           "| | UPGMA | NJ |", "|---|---|---|",
           f"| Suma de cuadrados de (d_árbol − d_observada) | {((P_up - D)[iu] ** 2).sum():.3e} | {((P_nj - D)[iu] ** 2).sum():.3e} |",
           f"| Máximo error absoluto | {np.abs(P_up - D)[iu].max():.5f} | {np.abs(P_nj - D)[iu].max():.5f} |",
           "", "Errores por par (d_árbol − d_observada):", "",
           "| Par | d observada | UPGMA | NJ |", "|---|---|---|---|"]
    for i, j in zip(*iu):
        md.append(f"| {ab(etiquetas[i])}–{ab(etiquetas[j])} | {D[i, j]:.5f} | "
                  f"{P_up[i, j] - D[i, j]:+.5f} | {P_nj[i, j] - D[i, j]:+.5f} |")
    md += ["",
           "### 2.5 ¿Cumple la matriz las hipótesis de cada método?", "",
           f"- Desviación máxima de la condición **ultramétrica** (supuesto de UPGMA): {ultrametricidad(D):.5f}",
           f"- Desviación máxima de la condición de **cuatro puntos / aditividad** (caso en que NJ es exacto): {cuatro_puntos(D):.5f}",
           "- Ninguna es 0: la matriz real no es ni ultramétrica ni perfectamente aditiva "
           "(ruido de muestreo, homoplasia, saturación).", ""]

    # ---------------- 3. Ejemplo sintético ----------------
    real = "((A:0.30,B:0.05):0.04,C:0.08,(D:0.07,E:0.10):0.06);"
    nombres_s = list("ABCDE")
    Ds = distancias_patristicas(real, nombres_s)
    nj_s = neighbor_joining(nombres_s, Ds).newick
    up_s = upgma(nombres_s, Ds)
    md += ["## 3. Ejemplo sintético sin reloj molecular", "",
           f"Árbol real (aditivo, no ultramétrico; A acumula 6 veces más cambios que su hermano B): `{real}`", "",
           "```", "Matriz derivada del árbol:",
           *["   " + " ".join(f"{v:5.2f}" for v in fila) for fila in Ds],
           f"NJ   : {nj_s}", f"UPGMA: {up_s}", "```", "",
           f"- NJ recupera el árbol real: RF = **{robinson_foulds(nj_s, real)}**, ramas exactas.",
           f"- UPGMA: RF = **{robinson_foulds(up_s, real)}**. Une B con C porque B, al evolucionar despacio, "
           "parece \"cercano\" a todo; separa A de su verdadero grupo hermano.", ""]

    (RES / "comparacion.md").write_text("\n".join(md), encoding="utf-8")
    dibujar_comparacion(nw_up, nw_nj, RES / "upgma-vs-nj.png")
    print("\n".join(md))
    print(f"\nGuardado: {RES / 'comparacion.md'} y {RES / 'upgma-vs-nj.png'}")


if __name__ == "__main__":
    main()
