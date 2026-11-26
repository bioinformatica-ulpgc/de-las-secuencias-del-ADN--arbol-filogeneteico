"""
Pruebas de validación de la implementación propia de Neighbor-Joining (CÓDIGO PROPIO).

Ejecutar desde la raíz del repositorio:
    python -m unittest -v
"""
import sys
import unittest
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

from arboles import biparticiones, distancias_patristicas, robinson_foulds, upgma  # noqa: E402
from nj import leer_matriz, neighbor_joining, validar_matriz  # noqa: E402

# Matriz didáctica de la teoría ("Métodos UPGMA y Neighbor-Joining")
NOMBRES_CLASE = ["Humano", "Chimpance", "Gorila", "Orangutan"]
D_CLASE = np.array([
    [0, 2, 6, 10],
    [2, 0, 6, 10],
    [6, 6, 0, 10],
    [10, 10, 10, 0],
], dtype=float)

# Árbol aditivo y NO ultramétrico de 5 taxones (A evoluciona mucho más rápido que B)
NEWICK_ADITIVO = "((A:0.30,B:0.05):0.04,C:0.08,(D:0.07,E:0.10):0.06);"


class TestValidacionMatriz(unittest.TestCase):
    def test_acepta_matriz_correcta(self):
        validar_matriz(NOMBRES_CLASE, D_CLASE)

    def test_tolera_ruido_numerico(self):
        D = D_CLASE.copy()
        D[0, 1] += 1e-12           # asimetría por redondeo
        D[2, 2] = -1e-13           # "cero" negativo por redondeo
        validar_matriz(NOMBRES_CLASE, D)

    def test_rechaza_no_cuadrada(self):
        with self.assertRaisesRegex(ValueError, "cuadrada"):
            validar_matriz(NOMBRES_CLASE, D_CLASE[:3, :])

    def test_rechaza_asimetrica(self):
        D = D_CLASE.copy()
        D[0, 2] = 7
        with self.assertRaisesRegex(ValueError, "simétrica"):
            validar_matriz(NOMBRES_CLASE, D)

    def test_rechaza_diagonal_no_nula(self):
        D = D_CLASE.copy()
        D[1, 1] = 0.5
        with self.assertRaisesRegex(ValueError, "diagonal"):
            validar_matriz(NOMBRES_CLASE, D)

    def test_rechaza_negativos(self):
        D = D_CLASE.copy()
        D[0, 3] = D[3, 0] = -1
        with self.assertRaisesRegex(ValueError, "negativas"):
            validar_matriz(NOMBRES_CLASE, D)

    def test_rechaza_nan(self):
        D = D_CLASE.copy()
        D[0, 3] = D[3, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finitos"):
            validar_matriz(NOMBRES_CLASE, D)

    def test_lee_csv_de_la_practica(self):
        nombres, D = leer_matriz(RAIZ / "proporcionado" / "distancias.csv")
        self.assertEqual(len(nombres), 5)
        validar_matriz(nombres, D)


class TestEjemploDeClase(unittest.TestCase):
    """Reproduce paso a paso el ejemplo resuelto en la teoría."""

    def setUp(self):
        self.res = neighbor_joining(NOMBRES_CLASE, D_CLASE)
        self.it1 = self.res.iteraciones[0]

    def test_r(self):
        self.assertEqual(self.it1["r"], {"Humano": 18, "Chimpance": 18, "Gorila": 22, "Orangutan": 30})

    def test_matriz_Q(self):
        Q = np.array(self.it1["Q"])
        esperado = np.array([
            [0, -32, -28, -28],
            [-32, 0, -28, -28],
            [-28, -28, 0, -32],
            [-28, -28, -32, 0],
        ])
        np.testing.assert_allclose(Q, esperado)

    def test_empate_se_resuelve_como_en_clase(self):
        # Q(H,C) = Q(G,O) = -32; en la teoría se elige Humano-Chimpancé
        self.assertEqual(self.it1["par"], ["Chimpance", "Humano"])

    def test_ramas_y_distancias_nuevas(self):
        self.assertAlmostEqual(self.it1["ramas"]["Humano"], 1)
        self.assertAlmostEqual(self.it1["ramas"]["Chimpance"], 1)
        self.assertAlmostEqual(self.it1["distancias_nuevo_nodo"]["Gorila"], 5)
        self.assertAlmostEqual(self.it1["distancias_nuevo_nodo"]["Orangutan"], 9)

    def test_union_final_2_3_7(self):
        final = self.res.iteraciones[-1]
        self.assertTrue(final.get("union_final"))
        self.assertAlmostEqual(final["ramas"]["U1"], 2)
        self.assertAlmostEqual(final["ramas"]["Gorila"], 3)
        self.assertAlmostEqual(final["ramas"]["Orangutan"], 7)

    def test_particion_HC_GO(self):
        splits = biparticiones(self.res.newick)
        self.assertIn(frozenset({"Gorila", "Orangutan"}), splits)  # = {H,C} | {G,O}

    def test_distancias_patristicas_reproducen_la_matriz(self):
        # la matriz de clase es aditiva: el árbol NJ la reproduce exactamente
        P = distancias_patristicas(self.res.newick, NOMBRES_CLASE)
        np.testing.assert_allclose(P, D_CLASE, atol=1e-9)


class TestArbolAditivo(unittest.TestCase):
    """Con una matriz aditiva, NJ debe recuperar el árbol exacto (topología y
    ramas), aunque las tasas no sean constantes. UPGMA, que supone reloj
    molecular, falla en este caso."""

    def setUp(self):
        self.nombres = list("ABCDE")
        self.D = distancias_patristicas(NEWICK_ADITIVO, self.nombres)

    def test_nj_recupera_topologia_y_ramas(self):
        res = neighbor_joining(self.nombres, self.D)
        self.assertEqual(robinson_foulds(res.newick, NEWICK_ADITIVO), 0)
        s_nj, s_real = biparticiones(res.newick), biparticiones(NEWICK_ADITIVO)
        self.assertEqual(set(s_nj), set(s_real))
        for s in s_real:
            self.assertAlmostEqual(s_nj[s], s_real[s], places=9)

    def test_upgma_falla_sin_reloj_molecular(self):
        self.assertGreater(robinson_foulds(upgma(self.nombres, self.D), NEWICK_ADITIVO), 0)

    def test_tres_taxones(self):
        res = neighbor_joining(["X", "Y", "Z"], np.array([[0, 3, 4], [3, 0, 5], [4, 5, 0]], float))
        self.assertEqual(res.iteraciones[-1]["ramas"], {"X": 1.0, "Y": 2.0, "Z": 3.0})

    def test_dos_taxones(self):
        res = neighbor_joining(["X", "Y"], np.array([[0, 3], [3, 0]], float))
        self.assertEqual(distancias_patristicas(res.newick, ["X", "Y"])[0, 1], 3)


class TestCOX1(unittest.TestCase):
    def setUp(self):
        self.nombres, self.D = leer_matriz(RAIZ / "proporcionado" / "distancias.csv")
        self.res = neighbor_joining(self.nombres, self.D)

    def test_coincide_con_biopython(self):
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        etiquetas = [n.replace(" ", "_") for n in self.nombres]
        dm = DistanceMatrix(etiquetas, [[float(self.D[i, j]) for j in range(i + 1)]
                                        for i in range(len(etiquetas))])
        # se compara el objeto árbol (no su Newick, cuyo nº de decimales depende de la versión)
        ref = DistanceTreeConstructor().nj(dm)
        self.assertEqual(robinson_foulds(self.res.newick, ref), 0)
        s1, s2 = biparticiones(self.res.newick), biparticiones(ref)
        self.assertEqual(set(s1), set(s2))
        for s in s1:
            self.assertAlmostEqual(s1[s], s2[s], places=8)

    def test_punto_medio_equidista_de_las_hojas_mas_alejadas(self):
        from arboles import enraizar_punto_medio
        arbol = enraizar_punto_medio(self.res.newick)
        prof = {t.name: d for t, d in arbol.depths().items() if t.is_terminal()}
        self.assertAlmostEqual(prof["Hylobates_lar"], prof["Pongo_pygmaeus"], places=9)

    def test_upgma_reescrito_reproduce_el_proporcionado(self):
        original = (RAIZ / "proporcionado" / "arbol-upgma.nwk").read_text().strip()
        self.assertEqual(upgma(self.nombres, self.D), original)


if __name__ == "__main__":
    unittest.main()
