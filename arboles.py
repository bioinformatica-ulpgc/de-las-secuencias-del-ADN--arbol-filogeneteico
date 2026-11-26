"""
arboles.py - Utilidades para comparar árboles (CÓDIGO PROPIO).

* biparticiones(): convierte un árbol (Newick) en su conjunto de biparticiones
  ("splits") con la longitud de rama asociada. Dos árboles no enraizados son
  iguales si y solo si tienen las mismas biparticiones, por muy distinto que se
  escriba su Newick (orden de los hijos, punto donde se "cuelga" el árbol...).
* robinson_foulds(): número de biparticiones no triviales que no comparten.
* distancias_patristicas(): distancia a lo largo del árbol entre cada par de hojas.
* upgma(): reescritura, como función reutilizable, del bucle UPGMA de
  proporcionado/procesar_datos.py (misma regla de selección, misma media
  ponderada y mismos desempates). Se usa solo para aplicar UPGMA a matrices de
  prueba; se comprueba en los tests que reproduce exactamente arbol-upgma.nwk.

Bibliotecas: Bio.Phylo (solo para leer Newick y recorrer el árbol), numpy.
"""
from __future__ import annotations

import io
from itertools import combinations

import numpy as np
from Bio import Phylo


def leer_newick(newick):
    """Acepta un Newick (str) o un árbol de Bio.Phylo ya construido."""
    if not isinstance(newick, str):
        return newick
    return Phylo.read(io.StringIO(newick.strip()), "newick")


def biparticiones(newick) -> dict[frozenset, float]:
    """Devuelve {bipartición: longitud}. Cada bipartición se representa por el
    lado que NO contiene a la hoja de referencia (la primera en orden alfabético),
    de modo que la representación no depende de dónde esté la raíz del Newick.
    Si la raíz es bifurcada (árbol enraizado), sus dos ramas forman una única
    arista del árbol no enraizado y sus longitudes se suman."""
    arbol = leer_newick(newick)
    hojas = frozenset(t.name for t in arbol.get_terminals())
    ref = min(hojas)
    res: dict[frozenset, float] = {}
    for clado in arbol.find_clades():
        if clado is arbol.root:
            continue
        lado = frozenset(t.name for t in clado.get_terminals())
        if ref in lado:
            lado = hojas - lado
        if not lado or lado == hojas:
            continue
        res[lado] = res.get(lado, 0.0) + (clado.branch_length or 0.0)
    return res


def no_triviales(splits: dict[frozenset, float], n_hojas: int) -> set[frozenset]:
    """Biparticiones internas (ambos lados con ≥ 2 hojas)."""
    return {s for s in splits if 2 <= len(s) <= n_hojas - 2}


def robinson_foulds(nw1, nw2) -> int:
    s1, s2 = biparticiones(nw1), biparticiones(nw2)
    n = len(leer_newick(nw1).get_terminals())
    return len(no_triviales(s1, n) ^ no_triviales(s2, n))


def distancias_patristicas(newick: str, nombres: list[str]) -> np.ndarray:
    arbol = leer_newick(newick)
    M = np.zeros((len(nombres), len(nombres)))
    for (i, a), (j, b) in combinations(enumerate(nombres), 2):
        M[i, j] = M[j, i] = arbol.distance(a, b)
    return M


def describir_split(s: frozenset, hojas: set[str], abrev=lambda x: x) -> str:
    otro = sorted(hojas - s)
    return "{" + ", ".join(abrev(x) for x in sorted(s)) + "} / {" + ", ".join(abrev(x) for x in otro) + "}"


def upgma(nombres: list[str], D: np.ndarray) -> str:
    """UPGMA con el mismo criterio que procesar_datos.py: se une el par de grupos
    de menor distancia (desempate por nombre), altura h = d/2, rama = h - altura
    del hijo, y d(A∪B, C) = (|A|·d(A,C) + |B|·d(B,C)) / (|A| + |B|)."""
    dists = {frozenset([nombres[i], nombres[j]]): float(D[i, j])
             for i, j in combinations(range(len(nombres)), 2)}
    grupos = {n: {"miembros": [n], "altura": 0.0, "newick": n.replace(" ", "_")} for n in nombres}
    paso = 0
    while len(grupos) > 1:
        clave = min(dists, key=lambda k: (dists[k], tuple(sorted(k))))
        x, y = sorted(clave)
        d = dists[clave]
        h = d / 2
        cx, cy = grupos[x], grupos[y]
        paso += 1
        nuevo = f"G{paso}"
        nx, ny = len(cx["miembros"]), len(cy["miembros"])
        nd = {z: (nx * dists[frozenset([x, z])] + ny * dists[frozenset([y, z])]) / (nx + ny)
              for z in grupos if z not in (x, y)}
        nw = f"({cx['newick']}:{h - cx['altura']:.9f},{cy['newick']}:{h - cy['altura']:.9f})"
        dists = {k: v for k, v in dists.items() if x not in k and y not in k}
        del grupos[x], grupos[y]
        grupos[nuevo] = {"miembros": cx["miembros"] + cy["miembros"], "altura": h, "newick": nw}
        for z, v in nd.items():
            dists[frozenset([nuevo, z])] = v
    return next(iter(grupos.values()))["newick"] + ";"


def enraizar_punto_medio(newick: str):
    """Enraíza en el punto medio del camino más largo entre dos hojas.
    (Se implementa aquí porque Tree.root_at_midpoint() de Biopython 1.8x no deja
    la raíz exactamente a mitad de camino en este árbol.)
    Caso cubierto: el punto medio cae dentro de una rama terminal, que es lo que
    ocurre con los árboles de esta práctica; si no, se recurre a Biopython."""
    arbol = leer_newick(newick)
    hojas = arbol.get_terminals()
    a, b = max(combinations(hojas, 2), key=lambda p: arbol.distance(*p))
    mitad = arbol.distance(a, b) / 2
    for hoja in (a, b):
        if hoja.branch_length is not None and mitad <= hoja.branch_length:
            arbol.root_with_outgroup(hoja.name, outgroup_branch_length=mitad)
            return arbol
    arbol.root_at_midpoint()
    return arbol
