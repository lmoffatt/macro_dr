figuras 2 (3 lineas), 3 (1 linea), 4 (0.5 linea)  caption surpassed limit: encogé la figura o reduci el texto.

HECHO 2026-08-26. Las tres entran, medido en el PDF (compilado a `tmp/capfit/`, paginas 12, 14 y 16),
no estimado. Se redujo el texto, NO se encogieron las figuras: estan dibujadas a ancho completo con
tipografia base de 7 pt, que es el piso de eLife, asi que escalarlas hubiera puesto sus rotulos por
debajo del limite.

De donde salieron las lineas, y por que cada corte es seguro:
- Figura 2: las glosas de magnitud y anisotropia (la prosa de esa misma subseccion las da textuales,
  al lado de la ecuacion que las define); "free in every likelihood member and simulated at 1e-4"
  (material de Methods, y ese 1e-4 sin unidades ya venia marcado de la ronda del 2026-08-25);
  "scaled by the reciprocal of the group size" -> "divided by the group size"; y el titulo de source
  data 1, que pasa a una linea.
- Figura 3: la oracion "The channel number is not drawn for either least-squares arm, which cannot
  separate it from the unitary current", que sigue dicha en la caption del suplemento 2 de esa misma
  figura y en la prosa de Results (con la cita de identifiabilidad).
- Figura 4: "and the difference is not cosmetic" (editorializa, no mide) y el titulo del suplemento 4,
  que pasa a una linea.
Ningun numero medido se perdio.

QUEDA ABIERTO, del mismo tipo y NO tocado porque no estaba en la lista: dos paginas de suplementos se
pasan de la caja segun el log de LaTeX, `Figure 5--figure supplement 3` por 121.9 pt (unas diez lineas,
grande) y `Figure 5--figure supplement 2` por 6.7 pt (media linea). Son previas a esta pasada. La
supplement 3 es la de celda mas grande del set (2.7 mm, ocho columnas), asi que ahi la salida puede ser
la figura y no la caption; hay que medirla antes de decidir.
