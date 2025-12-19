
## Lab - Celda de carga 
Comenzaremos con una simulación teórica de FEM que unifique la ecuación constitutiva con el comportamiento esperado. Se usará para proponer un acondicionamiento adecuado que luego será empleado en físico como un bono en la materia. 

<div style="display: flex; justify-content: center;">
  <img src="https://github.com/mora200217/strain_gauge_sim/blob/main/docs/sim_flow_diagram.png" width="90%"/>
</div>

## Elasticidad lineal y sensores (Gridap.jl + Julia)

Julia es un lenguaje de programación científico. En este curso usaremos el módulo
**Gridap.jl** como herramienta para visualizar y entender cómo la **geometría**
y el **material** influyen en la deformación que mide un sensor
(p. ej. un *strain gauge*).

El objetivo **no es profundizar en FEM**, sino usarlo como apoyo para comprender
el funcionamiento de las **celdas de carga**.

---

## Ley constitutiva (idea física)

En un caso uniaxial, un material elástico obedece la ley de Hooke:

$$
\sigma = E\,\varepsilon
$$

donde:
- $( \sigma )$ es el esfuerzo,
- $( \varepsilon )$ es la deformación,
- $( E )$ es el módulo de Young.

---

## Generalización triaxial

En un sólido real, las deformaciones en una dirección dependen de los esfuerzos
en las otras direcciones (efecto de Poisson):

$$
\varepsilon_x = \frac{1}{E}\left(\sigma_x - \nu(\sigma_y + \sigma_z)\right)
$$

$$
\varepsilon_y = \frac{1}{E}\left(\sigma_y - \nu(\sigma_x + \sigma_z)\right)
$$

$$
\varepsilon_z = \frac{1}{E}\left(\sigma_z - \nu(\sigma_x + \sigma_y)\right)
$$

Aquí \( \nu \) es el **coeficiente de Poisson**.

---

## Forma tensorial (compacta)

Para escribir estas relaciones de forma general y compacta se usan tensores:

$$
\sigma_{ij}
=
\lambda\,\varepsilon_{kk}\,\delta_{ij}
+
2\mu\,\varepsilon_{ij}
$$

donde:
- $( \varepsilon_{ij} )$ es el tensor de deformación,
- $( \sigma_{ij} )$ es el tensor de esfuerzos,
- $( \lambda )$ y $( \mu )$ son constantes del material (parámetros de Lamé).

Esta expresión es solo una forma elegante de escribir las ecuaciones anteriores.

---

## Deformación a partir del desplazamiento

La deformación se obtiene a partir del campo de desplazamientos \( \mathbf{u} \):

$$

\varepsilon_{ij}(\mathbf{u})
=
\frac{1}{2}
\left(
\frac{\partial u_i}{\partial x_j}
+
\frac{\partial u_j}{\partial x_i}
\right)

$$

Esto nos dice **cuánto se estira o comprime el material**, que es exactamente
lo que mide un *strain gauge*.

---

## Idea clave para sensores

- La **fuerza aplicada** genera esfuerzos.
- Los esfuerzos producen **deformaciones**.
- El sensor mide **deformación**, no esfuerzo.
- La **geometría de la celda** controla dónde se concentra la deformación.

Gridap.jl se utiliza aquí como una **herramienta visual y numérica**
para explorar estas relaciones, no como un fin en sí mismo.

## Makefile 
Existe un makefile que facilita el flujo de simulación

Para visualizar el enmallado generado, usaremos
```bash
    make gmsh
```
