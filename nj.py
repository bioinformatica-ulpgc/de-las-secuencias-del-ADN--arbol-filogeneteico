"""
nj.py - Implementación propia del método Neighbor-Joining (NJ).

CÓDIGO PROPIO (no proporcionado por la profesora).

Lee la matriz de distancias p generada por proporcionado/procesar_datos.py
(distancias.csv), la valida, reconstruye el árbol no enraizado mediante
Neighbor-Joining y exporta:

    arbol-nj.nwk         árbol en formato Newick (no enraizado, trifurcación central)
    iteraciones-nj.csv   resumen de cada iteración (par unido, Q mínimo, ramas, distancias nuevas)
    iteraciones-nj.json  registro completo (matrices D y Q, r(i), ramas) de cada iteración
    traza-nj.txt         la misma traza en texto legible
    arbol-nj.png         representación gráfica (radial no enraizada + filograma)

Notación (la misma que en la teoría de la asignatura):
    r(i)   = suma de las distancias de i al resto de nodos activos
    Q(i,j) = (n - 2) * d(i,j) - r(i) - r(j)          -> se une el par con Q mínimo
    L(i,U) = d(i,j)/2 + [r(i) - r(j)] / [2(n - 2)]
    L(j,U) = d(i,j) - L(i,U)
    d(U,k) = [d(i,k) + d(j,k) - d(i,j)] / 2
    Unión final (3 nodos a, b, c unidos a un nodo central):
    L(a) = [d(a,b) + d(a,c) - d(b,c)] / 2   (y análogamente para b y c)

No se usa ninguna función de biblioteca que implemente Neighbor-Joining.
Bibliotecas usadas: numpy (matrices), csv/json (E/S), matplotlib (dibujo).

Uso:
    python nj.py                                   # usa proporcionado/distancias.csv
    python nj.py --entrada otra.csv --salida dir/  # otra matriz
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parent
TOL = 1e-9  # tolerancia numérica para validar la matriz y para detectar empates en Q


# --------------------------------------------------------------------------- #
# 1. Lectura y validación de la matriz
# --------------------------------------------------------------------------- #
def leer_matriz(ruta: str | Path) -> tuple[list[str], np.ndarray]:
    """Lee un CSV con cabecera [etiqueta, taxón1, taxón2, ...] y una fila por taxón
    (formato de distancias.csv). Devuelve (nombres, matriz)."""
    with open(ruta, newline="", encoding="utf-8") as h:
        filas = [f for f in csv.reader(h) if f]
    cabecera, cuerpo = filas[0][1:], filas[1:]
    nombres_filas = [f[0] for f in cuerpo]
    if nombres_filas != cabecera:
        raise ValueError(
            f"Los nombres de filas {nombres_filas} no coinciden con los de columnas {cabecera}"
        )
    try:
        valores = [[float(x) for x in f[1:]] for f in cuerpo]
    except ValueError as e:
        raise ValueError(f"La matriz contiene un valor no numérico: {e}") from None
    if any(len(f) != len(cabecera) for f in valores):
        raise ValueError("La matriz no es cuadrada: alguna fila tiene un número de columnas distinto")
    return cabecera, np.array(valores, dtype=float)


def validar_matriz(nombres: list[str], D: np.ndarray, tol: float = TOL) -> None:
    """Comprueba que D es cuadrada, simétrica, con diagonal nula y sin valores
    negativos (con tolerancia numérica). Lanza ValueError si algo falla."""
    D = np.asarray(D, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError(f"La matriz no es cuadrada: forma {D.shape}")
    n = D.shape[0]
    if len(nombres) != n:
        raise ValueError(f"Hay {len(nombres)} nombres para una matriz {n}x{n}")
    if len(set(nombres)) != n:
        raise ValueError("Hay nombres de taxones repetidos")
    if n < 2:
        raise ValueError("Se necesitan al menos 2 taxones")
    if not np.all(np.isfinite(D)):
        raise ValueError("La matriz contiene valores no finitos (NaN o infinito)")
    asim = np.abs(D - D.T)
    if asim.max() > tol:
        i, j = np.unravel_index(asim.argmax(), asim.shape)
        raise ValueError(
            f"La matriz no es simétrica: d({nombres[i]},{nombres[j]})={D[i, j]} "
            f"≠ d({nombres[j]},{nombres[i]})={D[j, i]}"
        )
    diag = np.abs(np.diag(D))
    if diag.max() > tol:
        i = int(diag.argmax())
        raise ValueError(f"La diagonal no es nula: d({nombres[i]},{nombres[i]})={D[i, i]}")
    if D.min() < -tol:
        i, j = np.unravel_index(D.argmin(), D.shape)
        raise ValueError(f"Hay distancias negativas: d({nombres[i]},{nombres[j]})={D[i, j]}")


# --------------------------------------------------------------------------- #
# 2. Algoritmo Neighbor-Joining
# --------------------------------------------------------------------------- #
@dataclass
class ResultadoNJ:
    newick: str                      # árbol en formato Newick
    aristas: list[tuple[str, str, float]]  # (nodo_interno, hijo, longitud)
    iteraciones: list[dict] = field(default_factory=list)
    hojas: list[str] = field(default_factory=list)
    centro: str = ""                 # nodo interno al que se unen los 3 últimos nodos


def _etiqueta_newick(nombre: str) -> str:
    """Los espacios no son válidos en Newick sin comillas: se sustituyen por '_'
    (igual que hace procesar_datos.py)."""
    return nombre.replace(" ", "_")


def _elegir_par(Q: np.ndarray, activos: list[str], tol: float) -> tuple[int, int]:
    """Devuelve (i, j) con Q mínimo. Los empates (dentro de la tolerancia) se
    resuelven por orden alfabético de los nombres, como hace el UPGMA de
    procesar_datos.py; así se reproduce también la elección Humano-Chimpancé
    del ejemplo de clase, donde Q(H,C) = Q(G,O) = -32."""
    n = len(activos)
    qmin = min(Q[i, j] for i in range(n) for j in range(i + 1, n))
    candidatos = [
        (tuple(sorted((activos[i], activos[j]))), i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if Q[i, j] <= qmin + tol
    ]
    _, i, j = min(candidatos)
    # que i sea el primero en orden alfabético (solo afecta al orden de impresión)
    if activos[i] > activos[j]:
        i, j = j, i
    return i, j


def neighbor_joining(nombres: list[str], D: np.ndarray, tol: float = TOL) -> ResultadoNJ:
    """Reconstruye un árbol no enraizado por Neighbor-Joining (Saitou y Nei, 1987).

    Devuelve un ResultadoNJ con el Newick, la lista de aristas y el registro
    completo de cada iteración."""
    validar_matriz(nombres, D, tol)
    D = np.array(D, dtype=float)
    activos = list(nombres)                      # nodos que siguen en la matriz
    subarbol = {x: _etiqueta_newick(x) for x in nombres}  # Newick parcial de cada nodo
    aristas: list[tuple[str, str, float]] = []
    iteraciones: list[dict] = []
    contador = 0

    if len(activos) == 2:  # caso trivial: una sola rama entre los dos taxones
        a, b = activos
        newick = f"({subarbol[a]}:{D[0, 1]:.9f},{subarbol[b]}:0.000000000);"
        return ResultadoNJ(newick, [(a, b, float(D[0, 1]))], [], list(nombres), a)

    while len(activos) > 3:
        n = len(activos)
        r = D.sum(axis=1)                                   # r(i)
        Q = (n - 2) * D - r[:, None] - r[None, :]           # Q(i,j)
        np.fill_diagonal(Q, 0.0)                            # la diagonal no se usa
        i, j = _elegir_par(Q, activos, tol)
        a, b = activos[i], activos[j]
        dij = D[i, j]

        # longitudes de las dos ramas nuevas
        Li = dij / 2 + (r[i] - r[j]) / (2 * (n - 2))
        Lj = dij - Li

        contador += 1
        U = f"U{contador}"
        aristas += [(U, a, float(Li)), (U, b, float(Lj))]
        subarbol[U] = f"({subarbol[a]}:{Li:.9f},{subarbol[b]}:{Lj:.9f})"

        # distancias del nuevo nodo U al resto
        resto = [k for k in range(n) if k not in (i, j)]
        dU = {activos[k]: float((D[i, k] + D[j, k] - dij) / 2) for k in resto}

        iteraciones.append({
            "iteracion": contador,
            "n": n,
            "nodos": list(activos),
            "D": D.tolist(),
            "r": {activos[k]: float(r[k]) for k in range(n)},
            "Q": Q.tolist(),
            "par": [a, b],
            "Q_min": float(Q[i, j]),
            "d_par": float(dij),
            "nuevo_nodo": U,
            "ramas": {a: float(Li), b: float(Lj)},
            "distancias_nuevo_nodo": dU,
        })

        # nueva matriz: se eliminan i, j y se añade U al final
        nuevos = [activos[k] for k in resto] + [U]
        Dn = np.zeros((n - 1, n - 1))
        Dn[: n - 2, : n - 2] = D[np.ix_(resto, resto)]
        for idx, k in enumerate(resto):
            Dn[idx, -1] = Dn[-1, idx] = dU[activos[k]]
        activos, D = nuevos, Dn

    # ---- unión final: los 3 nodos restantes se unen a un nodo central ----
    a, b, c = activos
    dab, dac, dbc = D[0, 1], D[0, 2], D[1, 2]
    La = (dab + dac - dbc) / 2
    Lb = (dab + dbc - dac) / 2
    Lc = (dac + dbc - dab) / 2
    centro = f"U{contador + 1}"
    aristas += [(centro, a, float(La)), (centro, b, float(Lb)), (centro, c, float(Lc))]
    iteraciones.append({
        "iteracion": contador + 1,
        "n": 3,
        "nodos": list(activos),
        "D": D.tolist(),
        "union_final": True,
        "nuevo_nodo": centro,
        "ramas": {a: float(La), b: float(Lb), c: float(Lc)},
    })
    newick = f"({subarbol[a]}:{La:.9f},{subarbol[b]}:{Lb:.9f},{subarbol[c]}:{Lc:.9f});"
    return ResultadoNJ(newick, aristas, iteraciones, list(nombres), centro)


# --------------------------------------------------------------------------- #
# 3. Registro legible de las iteraciones
# --------------------------------------------------------------------------- #
ABREV = {
    "Homo sapiens": "Humano",
    "Pan troglodytes": "Chimpancé",
    "Gorilla gorilla": "Gorila",
    "Pongo pygmaeus": "Orangután",
    "Hylobates lar": "Gibón",
}


def _ab(x: str) -> str:
    return ABREV.get(x, x)


def _fmt_matriz(etiquetas: list[str], M, dec: int = 5) -> str:
    ancho = max(10, dec + 5, *(len(_ab(e)) + 1 for e in etiquetas))
    lineas = [" " * ancho + "".join(f"{_ab(e):>{ancho}}" for e in etiquetas)]
    for e, fila in zip(etiquetas, M):
        lineas.append(f"{_ab(e):<{ancho}}" + "".join(f"{v:>{ancho}.{dec}f}" for v in fila))
    return "\n".join(lineas)


def traza_texto(res: ResultadoNJ, dec: int = 5) -> str:
    out = []
    for it in res.iteraciones:
        if it.get("union_final"):
            out.append(f"=== Unión final (quedan 3 nodos) -> nodo central {it['nuevo_nodo']} ===")
            out.append("Matriz D:\n" + _fmt_matriz(it["nodos"], it["D"], dec))
            a, b, c = it["nodos"]
            out.append(f"L(a) = [d(a,b) + d(a,c) - d(b,c)] / 2, etc.")
            for x in (a, b, c):
                out.append(f"  L({_ab(x)}, {it['nuevo_nodo']}) = {it['ramas'][x]:.{dec}f}")
            out.append("")
            continue
        n = it["n"]
        a, b = it["par"]
        out.append(f"=== Iteración {it['iteracion']} (n = {n}) ===")
        out.append("Matriz D:\n" + _fmt_matriz(it["nodos"], it["D"], dec))
        out.append("r(i): " + ", ".join(f"{_ab(k)}={v:.{dec}f}" for k, v in it["r"].items()))
        out.append(f"Matriz Q = (n-2)·d(i,j) - r(i) - r(j):\n" + _fmt_matriz(it["nodos"], it["Q"], dec))
        out.append(f"Par seleccionado (Q mínimo = {it['Q_min']:.{dec}f}): {_ab(a)} + {_ab(b)}"
                   f"  ->  nuevo nodo {it['nuevo_nodo']}")
        out.append(f"  L({_ab(a)}, {it['nuevo_nodo']}) = d/2 + [r({_ab(a)}) - r({_ab(b)})]/[2(n-2)]"
                   f" = {it['ramas'][a]:.{dec}f}")
        out.append(f"  L({_ab(b)}, {it['nuevo_nodo']}) = d - L({_ab(a)}) = {it['ramas'][b]:.{dec}f}")
        out.append(f"Distancias al nuevo nodo d({it['nuevo_nodo']},k) = [d(i,k) + d(j,k) - d(i,j)] / 2:")
        for k, v in it["distancias_nuevo_nodo"].items():
            out.append(f"  d({it['nuevo_nodo']}, {_ab(k)}) = {v:.{dec}f}")
        out.append("")
    negativas = [(p, h, l) for p, h, l in res.aristas if l < 0]
    if negativas:
        out.append("AVISO: ramas con longitud negativa (posible en NJ si la matriz no es aditiva): "
                   + ", ".join(f"{p}-{_ab(h)}={l:.{dec}f}" for p, h, l in negativas))
    out.append("Newick: " + res.newick)
    return "\n".join(out)


def guardar_resultados(res: ResultadoNJ, salida: Path) -> None:
    salida.mkdir(parents=True, exist_ok=True)
    (salida / "arbol-nj.nwk").write_text(res.newick + "\n", encoding="utf-8")
    (salida / "iteraciones-nj.json").write_text(
        json.dumps({"hojas": res.hojas, "newick": res.newick, "aristas": res.aristas,
                    "iteraciones": res.iteraciones}, ensure_ascii=False, indent=2),
        encoding="utf-8")
    (salida / "traza-nj.txt").write_text(traza_texto(res, dec=6) + "\n", encoding="utf-8")
    with open(salida / "iteraciones-nj.csv", "w", newline="", encoding="utf-8") as h:
        w = csv.writer(h)
        w.writerow(["iteracion", "n", "nodo_A", "nodo_B", "Q_min", "d(A,B)", "nuevo_nodo",
                    "rama_A", "rama_B", "rama_C_final", "nuevas_distancias"])
        for it in res.iteraciones:
            if it.get("union_final"):
                a, b, c = it["nodos"]
                w.writerow([it["iteracion"], 3, a, b, "", "", it["nuevo_nodo"],
                            f"{it['ramas'][a]:.9f}", f"{it['ramas'][b]:.9f}",
                            f"{c}={it['ramas'][c]:.9f}", ""])
            else:
                a, b = it["par"]
                w.writerow([it["iteracion"], it["n"], a, b, f"{it['Q_min']:.9f}",
                            f"{it['d_par']:.9f}", it["nuevo_nodo"],
                            f"{it['ramas'][a]:.9f}", f"{it['ramas'][b]:.9f}", "",
                            "; ".join(f"{k}={v:.9f}" for k, v in it["distancias_nuevo_nodo"].items())])


# --------------------------------------------------------------------------- #
# 4. Representación gráfica
# --------------------------------------------------------------------------- #
def coordenadas_radiales(aristas, centro):
    """Disposición 'equal-angle' de un árbol no enraizado: cada subárbol recibe un
    sector angular proporcional a su número de hojas."""
    hijos: dict[str, list[tuple[str, float]]] = {}
    for p, h, l in aristas:
        hijos.setdefault(p, []).append((h, l))

    def n_hojas(x):
        return 1 if x not in hijos else sum(n_hojas(h) for h, _ in hijos[x])

    pos = {centro: (0.0, 0.0)}
    total = n_hojas(centro)

    def colocar(x, ang0, ang1):
        a = ang0
        for h, l in hijos.get(x, []):
            ancho = (ang1 - ang0) * n_hojas(h) / n_hojas(x)
            medio = a + ancho / 2
            x0, y0 = pos[x]
            pos[h] = (x0 + max(l, 0) * math.cos(medio), y0 + max(l, 0) * math.sin(medio))
            colocar(h, a, a + ancho)
            a += ancho

    colocar(centro, 0, 2 * math.pi * (total - 0) / total)
    return pos, hijos


def dibujar(res: ResultadoNJ, ruta_png: Path, titulo: str = "Neighbor-Joining") -> None:
    import io
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo  # solo para dibujar el filograma (no implementa NJ)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), facecolor="white")

    # (a) radial no enraizado
    pos, hijos = coordenadas_radiales(res.aristas, res.centro)
    for p, hs in hijos.items():
        for h, l in hs:
            (x0, y0), (x1, y1) = pos[p], pos[h]
            ax1.plot([x0, x1], [y0, y1], color="black", lw=2)
            # etiqueta desplazada perpendicularmente a la rama para no tapar los nombres
            dx, dy = x1 - x0, y1 - y0
            norma = math.hypot(dx, dy) or 1.0
            ax1.annotate(f"{l:.4f}", ((x0 + x1) / 2, (y0 + y1) / 2),
                         xytext=(-dy / norma * 9, dx / norma * 9), textcoords="offset points",
                         fontsize=8, color="black", ha="center", va="center")
    for h in res.hojas:
        x, y = pos[h]
        ax1.plot(x, y, "o", color="black")
        ax1.annotate(_ab(h), (x, y), xytext=(6, 4), textcoords="offset points", fontsize=11)
    ax1.set_aspect("equal")
    ax1.margins(0.25)
    ax1.axis("off")
    ax1.set_title(f"{titulo}: árbol no enraizado (radial)")

    # (b) filograma: Bio.Phylo lo dibuja "colgado" del nodo central (raíz arbitraria)
    arbol = Phylo.read(io.StringIO(res.newick), "newick")
    Phylo.draw(arbol, axes=ax2, do_show=False,
               label_func=lambda c: _ab(c.name.replace("_", " ")) if c.is_terminal() else None,
               branch_labels=lambda c: f"{c.branch_length:.4f}" if c.branch_length else None)
    ax2.set_title(f"{titulo}: filograma (raíz arbitraria en el nodo central)")
    ax2.set_xlabel("distancia p (diferencias por sitio)")
    ax2.set_ylabel("")
    ax2.set_yticks([])
    fig.tight_layout()
    fig.savefig(ruta_png, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description="Neighbor-Joining propio")
    ap.add_argument("--entrada", default=str(RAIZ / "proporcionado" / "distancias.csv"))
    ap.add_argument("--salida", default=str(RAIZ / "resultados"))
    ap.add_argument("--quiet", action="store_true", help="no imprimir la traza")
    ap.add_argument("--sin-figura", action="store_true")
    args = ap.parse_args()

    nombres, D = leer_matriz(args.entrada)
    validar_matriz(nombres, D)
    print(f"Matriz {D.shape[0]}x{D.shape[0]} leída de {args.entrada}: "
          "cuadrada, simétrica, diagonal nula y sin negativos ✔")
    res = neighbor_joining(nombres, D)
    if not args.quiet:
        print(traza_texto(res))
    salida = Path(args.salida)
    guardar_resultados(res, salida)
    if not args.sin_figura:
        dibujar(res, salida / "arbol-nj.png")
    print(f"Resultados guardados en {salida}/")


if __name__ == "__main__":
    main()
