"""
construir_informe.py - Convierte informe.md en informe.pdf (CÓDIGO PROPIO, auxiliar).

Traduce el subconjunto de Markdown usado en el informe (títulos, párrafos, listas,
tablas, figuras, negrita/cursiva, código y fórmulas LaTeX sencillas) a Typst y lo
compila a PDF con el paquete `typst` de PyPI (pip install typst).

Uso:  python informe/construir_informe.py
"""
import re
from pathlib import Path

DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------- matemáticas
def _llaves(s: str, i: int) -> tuple[str, int]:
    """Devuelve el contenido del grupo {...} que empieza en s[i] y el índice final."""
    assert s[i] == "{"
    nivel, j = 0, i
    while True:
        if s[j] == "{":
            nivel += 1
        elif s[j] == "}":
            nivel -= 1
            if nivel == 0:
                return s[i + 1:j], j + 1
        j += 1


def latex_a_typst(m: str) -> str:
    m = m.replace("{,}", ",")
    # \text{...} -> "..."
    m = re.sub(r"\\text\{([^}]*)\}", r'"\1"', m)
    # \frac{a}{b} -> (a)/(b)  (con llaves anidadas)
    while "\\frac" in m:
        i = m.index("\\frac")
        num, j = _llaves(m, i + 5)
        den, k = _llaves(m, j)
        m = m[:i] + f"({num})/({den})" + m[k:]
    # sub/superíndices con llaves
    out, i = "", 0
    while i < len(m):
        if m[i] in "_^" and i + 1 < len(m) and m[i + 1] == "{":
            cont, j = _llaves(m, i + 1)
            cont = re.sub(r"(?<=[a-z])(?=[a-z])", " ", cont)  # ij -> i j
            out += f"{m[i]}({cont})"
            i = j
        else:
            out += m[i]
            i += 1
    m = out
    for a, b in [(r"\cdot", " dot "), (r"\le", " <= "), (r"\neq", " != "), (r"\cup", " union "),
                 (r"\sum", " sum"), (r"\qquad", " quad "), (r"\,", " ")]:
        m = m.replace(a, b)
    m = re.sub(r"(?<=[A-Za-z])(?=\d)", " ", m)  # U1 -> U 1
    return m.strip()


# ---------------------------------------------------------------- texto en línea
def en_linea(t: str) -> str:
    guard: list[str] = []

    def proteger(s):
        guard.append(s)
        return f"\x00{len(guard) - 1}\x00"

    t = re.sub(r"\$([^$]+)\$", lambda mt: proteger("$" + latex_a_typst(mt.group(1)) + "$"), t)
    t = re.sub(r"`([^`]+)`", lambda mt: proteger("`" + mt.group(1) + "`"), t)
    t = re.sub(r"<(https?://[^>]+)>", lambda mt: proteger(f'#link("{mt.group(1)}")'), t)
    t = t.replace("\\*", proteger("\\*"))
    t = re.sub(r"\*\*(.+?)\*\*", lambda mt: proteger("*" + en_linea_simple(mt.group(1)) + "*"), t)
    t = re.sub(r"\*(.+?)\*", lambda mt: proteger("_" + en_linea_simple(mt.group(1)) + "_"), t)
    t = en_linea_simple(t)
    while "\x00" in t:
        t = re.sub(r"\x00(\d+)\x00", lambda mt: guard[int(mt.group(1))], t)
    return t


def en_linea_simple(t: str) -> str:
    # escapa los caracteres con significado especial en Typst
    return re.sub(r"(?<!\\)([#_@<>~])", r"\\\1", t)


# ---------------------------------------------------------------- bloques
def tabla(lineas: list[str]) -> str:
    filas = [[c.strip() for c in l.strip().strip("|").split("|")] for l in lineas]
    cab, cuerpo = filas[0], filas[2:]
    n = len(cab)
    celdas = [f"[*{en_linea(c)}*]" for c in cab]
    for f in cuerpo:
        celdas += [f"[{en_linea(c)}]" for c in f]
    return (f"#table(columns: {n}, align: left, inset: 4pt, stroke: 0.4pt + luma(170),\n"
            f"  fill: (x, y) => if y == 0 {{ luma(235) }},\n  " + ", ".join(celdas) + ")\n")


def convertir(md: str) -> str:
    meta = {}
    if md.startswith("---"):
        _, fm, md = md.split("---", 2)
        for l in fm.strip().splitlines():
            k, v = l.split(":", 1)
            meta[k.strip()] = v.strip().strip('"')
    out = [
        '#set document(title: "%s")' % meta.get("title", ""),
        '#set page(paper: "a4", margin: (x: 2cm, y: 1.8cm), numbering: "1")',
        '#set text(lang: "es", size: 10pt)',
        "#set par(justify: true, leading: 0.55em)",
        "#set heading(numbering: none)",
        "#show heading.where(level: 1): set text(size: 12.5pt, fill: rgb(\"#18324b\"))",
        "#show raw: set text(size: 8.5pt)",
        "#show figure.caption: set text(size: 8.5pt)",
        "#align(center)[#text(size: 17pt, weight: \"bold\", fill: rgb(\"#18324b\"))[%s] \\ "
        "#text(size: 11pt, style: \"italic\")[%s] \\ #v(2pt) %s \\ #text(size: 9pt)[%s]]"
        % (meta.get("title", ""), meta.get("subtitle", ""), meta.get("author", ""), meta.get("date", "")),
        "",
    ]
    lineas = md.strip().splitlines()
    i = 0
    while i < len(lineas):
        l = lineas[i]
        if l.startswith("|"):
            bloque = []
            while i < len(lineas) and lineas[i].startswith("|"):
                bloque.append(lineas[i])
                i += 1
            out.append(tabla(bloque))
            continue
        if l.startswith("$$"):
            out.append("$ " + latex_a_typst(l.strip().strip("$")) + " $")
        elif m := re.match(r"(#+) (.*)", l):
            out.append("=" * len(m.group(1)) + " " + en_linea(m.group(2)))
        elif m := re.match(r"!\[(.*)\]\((.*?)\)(\{width=(\d+)%\})?", l):
            w = m.group(4) or "100"
            out.append(f'#figure(image("{m.group(2)}", width: {w}%), caption: [{en_linea(m.group(1))}])')
        elif m := re.match(r"\* (.*)", l):
            out.append("- " + en_linea(m.group(1)))
        elif m := re.match(r"\d+\. (.*)", l):
            out.append("+ " + en_linea(m.group(1)))
        else:
            out.append(en_linea(l))
        i += 1
    return "\n".join(out) + "\n"


def main():
    import typst
    fuente = DIR / "informe.typ"
    fuente.write_text(convertir((DIR / "informe.md").read_text(encoding="utf-8")), encoding="utf-8")
    typst.compile(str(fuente), output=str(DIR / "informe.pdf"), root=str(DIR.parent))
    print("Generado", DIR / "informe.pdf")


if __name__ == "__main__":
    main()
