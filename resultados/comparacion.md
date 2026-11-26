# Verificación y comparación de árboles

Generado automáticamente por `comparar.py` a partir de `proporcionado/distancias.csv`.

## 1. NJ propio frente a la implementación de referencia (Biopython)

```
NJ propio   : (Homo_sapiens:0.044990272,Pan_troglodytes:0.041909858,(Gorilla_gorilla:0.052285992,(Hylobates_lar:0.078145266,Pongo_pygmaeus:0.074902724):0.019049935):0.008998054);
NJ Biopython: (((Hylobates_lar:0.078145266,Pongo_pygmaeus:0.074902724)Inner1:0.019049935,Gorilla_gorilla:0.052285992)Inner2:0.0089980544,Pan_troglodytes:0.041909858,Homo_sapiens:0.044990272)Inner3:0;
```

Los dos Newick se escriben distinto (orden de los hijos, nombres `Inner` de los nodos internos), pero describen el mismo árbol no enraizado si tienen las mismas biparticiones:

| Bipartición (arista) | Tipo | NJ propio | NJ Biopython |
|---|---|---|---|
| {Humano} / {Gorila, Gibón, Chimpancé, Orangután} | terminal | 0.044990 | 0.044990 |
| {Gibón} / {Gorila, Humano, Chimpancé, Orangután} | terminal | 0.078145 | 0.078145 |
| {Chimpancé} / {Gorila, Humano, Gibón, Orangután} | terminal | 0.041910 | 0.041910 |
| {Orangután} / {Gorila, Humano, Gibón, Chimpancé} | terminal | 0.074903 | 0.074903 |
| {Humano, Chimpancé} / {Gorila, Gibón, Orangután} | interna | 0.008998 | 0.008998 |
| {Gibón, Orangután} / {Gorila, Humano, Chimpancé} | interna | 0.019050 | 0.019050 |
| {Humano, Gibón, Chimpancé, Orangután} / {Gorila} | terminal | 0.052286 | 0.052286 |

- Distancia de Robinson-Foulds: **0** (0 = misma topología)
- Misma colección de aristas: **True**
- Máxima diferencia en longitud de rama: **3.75e-10**

## 2. UPGMA frente a NJ (mismas distancias p)

```
UPGMA: ((((Homo_sapiens:0.043450065,Pan_troglodytes:0.043450065):0.008916991,Gorilla_gorilla:0.052367056):0.020806312,Pongo_pygmaeus:0.073173368):0.002053610,Hylobates_lar:0.075226978);
NJ   : (Homo_sapiens:0.044990272,Pan_troglodytes:0.041909858,(Gorilla_gorilla:0.052285992,(Hylobates_lar:0.078145266,Pongo_pygmaeus:0.074902724):0.019049935):0.008998054);
```

### 2.1 Aristas (biparticiones) y longitudes

UPGMA se trata aquí como no enraizado (sus dos ramas de la raíz se suman en una sola arista).

| Bipartición (arista) | Tipo | UPGMA | NJ |
|---|---|---|---|
| {Humano} / {Gorila, Gibón, Chimpancé, Orangután} | terminal | 0.043450 | 0.044990 |
| {Gibón} / {Gorila, Humano, Chimpancé, Orangután} | terminal | 0.077281 | 0.078145 |
| {Chimpancé} / {Gorila, Humano, Gibón, Orangután} | terminal | 0.043450 | 0.041910 |
| {Orangután} / {Gorila, Humano, Gibón, Chimpancé} | terminal | 0.073173 | 0.074903 |
| {Humano, Chimpancé} / {Gorila, Gibón, Orangután} | interna | 0.008917 | 0.008998 |
| {Gibón, Orangután} / {Gorila, Humano, Chimpancé} | interna | 0.020806 | 0.019050 |
| {Humano, Gibón, Chimpancé, Orangután} / {Gorila} | terminal | 0.052367 | 0.052286 |

- Distancia de Robinson-Foulds UPGMA–NJ (no enraizados): **0** de 4 posibles
- Biparticiones internas solo en UPGMA: ninguna
- Biparticiones internas solo en NJ: ninguna

### 2.2 Agrupaciones como árboles enraizados

| Árbol | Clados (grupos con ≥ 2 especies, sin contar la raíz) |
|---|---|
| UPGMA (raíz implícita del método) | {Chimpancé, Humano}; {Chimpancé, Gorila, Humano}; {Chimpancé, Gorila, Humano, Orangután} |
| NJ enraizado con el gibón como outgroup | {Chimpancé, Humano}; {Chimpancé, Gorila, Humano}; {Chimpancé, Gorila, Humano, Orangután} |
| NJ enraizado en el punto medio | {Chimpancé, Humano}; {Chimpancé, Gorila, Humano}; {Chimpancé, Gorila, Humano, Orangután} |

### 2.3 Distancia de cada hoja a la raíz

| Especie | UPGMA | NJ (raíz en el punto medio) |
|---|---|---|
| Humano | 0.075227 | 0.074660 |
| Chimpancé | 0.075227 | 0.071579 |
| Gorila | 0.075227 | 0.072957 |
| Orangután | 0.075227 | 0.076524 |
| Gibón | 0.075227 | 0.076524 |

UPGMA: todas iguales (rango 0.0e+00) → ultramétrico. NJ: profundidades desiguales (rango 0.0049) → las tasas de cambio no son idénticas en todos los linajes. El punto medio del camino más largo cae sobre la rama del gibón (Gibón a 0.0765 de la raíz), así que el enraizamiento por punto medio coincide con el del outgroup.

### 2.4 Ajuste a la matriz de distancias

| | UPGMA | NJ |
|---|---|---|
| Suma de cuadrados de (d_árbol − d_observada) | 1.334e-04 | 1.266e-04 |
| Máximo error absoluto | 0.00584 | 0.00689 |

Errores por par (d_árbol − d_observada):

| Par | d observada | UPGMA | NJ |
|---|---|---|---|
| Humano–Chimpancé | 0.08690 | +0.00000 | -0.00000 |
| Humano–Gorila | 0.10895 | -0.00422 | -0.00268 |
| Humano–Orangután | 0.14462 | +0.00173 | +0.00332 |
| Humano–Gibón | 0.14916 | +0.00130 | +0.00203 |
| Chimpancé–Gorila | 0.10052 | +0.00422 | +0.00268 |
| Chimpancé–Orangután | 0.14332 | +0.00303 | +0.00154 |
| Chimpancé–Gibón | 0.15499 | -0.00454 | -0.00689 |
| Gorila–Orangután | 0.15110 | -0.00476 | -0.00486 |
| Gorila–Gibón | 0.14462 | +0.00584 | +0.00486 |
| Orangután–Gibón | 0.15305 | -0.00259 | +0.00000 |

### 2.5 ¿Cumple la matriz las hipótesis de cada método?

- Desviación máxima de la condición **ultramétrica** (supuesto de UPGMA): 0.01038
- Desviación máxima de la condición de **cuatro puntos / aditividad** (caso en que NJ es exacto): 0.01816
- Ninguna es 0: la matriz real no es ni ultramétrica ni perfectamente aditiva (ruido de muestreo, homoplasia, saturación).

## 3. Ejemplo sintético sin reloj molecular

Árbol real (aditivo, no ultramétrico; A acumula 6 veces más cambios que su hermano B): `((A:0.30,B:0.05):0.04,C:0.08,(D:0.07,E:0.10):0.06);`

```
Matriz derivada del árbol:
    0.00  0.35  0.42  0.47  0.50
    0.35  0.00  0.17  0.22  0.25
    0.42  0.17  0.00  0.21  0.24
    0.47  0.22  0.21  0.00  0.17
    0.50  0.25  0.24  0.17  0.00
NJ   : (C:0.080000000,(D:0.070000000,E:0.100000000):0.060000000,(A:0.300000000,B:0.050000000):0.040000000);
UPGMA: (A:0.217500000,((B:0.085000000,C:0.085000000):0.030000000,(D:0.085000000,E:0.085000000):0.030000000):0.102500000);
```

- NJ recupera el árbol real: RF = **0**, ramas exactas.
- UPGMA: RF = **2**. Une B con C porque B, al evolucionar despacio, parece "cercano" a todo; separa A de su verdadero grupo hermano.
