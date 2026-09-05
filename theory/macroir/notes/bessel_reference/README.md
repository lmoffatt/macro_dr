# bessel_reference — oráculo de la recursión aumentada

2026-08-31. Rol: lo que `papers/1_method/decisions/recompute/mr_vs_ir_boxes.py`
fue para MR/IR — la implementación independiente que fija las fórmulas ANTES
del C++. Los hot paths de producción se validan contra esto.

- `bessel_oracle.py` — simulador exacto del proceso filtrado (canales CTMC +
  filtro por tramos + ruido por Gramiano exacto) + la recursión aumentada
  genérica (nativo K+4 / agrupado K+5 / box K+1 con el mismo código) + 5 gates.
  Correr: `python3 bessel_oracle.py` (~1-2 min). Última línea = veredicto.
- `augmented_recursion_blocks.md` — las ecuaciones de bloques que el C++
  transcribe, con el mapeo a los sitios de qmodel.h.

Compañeros de producción (C++):
- `legacy/acquisition_filter.h` (kernels fase 0) + `tests/math/test_acquisition_filter.cpp`
- `legacy/qdtf_engine.h` (2026-08-31, después): el CUERPO de calc_Qdtf_eig
  (productor espectral de ventana por corrimiento de polos, ensamblaje modal→
  real vía T) y la recursión predict/update, en double, autocontenido.
  Gates: `tests/math/test_qdtf_engine.cpp` (668 aserciones: tablas vs
  cuadratura para nativo/agrupado/box, ancla ν=0 vs ruta Ee/E3 clásica con la
  identidad de colapso en vivo, assembly vs Σ^bnd, invariantes multi-ventana).
  Cruce directo: ventana 1 nativa 0.789411/0.106484 y agrupada
  3.318923/1.214291 == oráculo python (MC-validado) a todos los dígitos.
  El calc_Qdtf_eig de qmodel.h queda como ADAPTADOR fino sobre este motor.

Integración (2026-08-31, tercera pasada):
- `legacy/qdtf_member.h` — el loop del miembro sobre modelo/experimento/registro
  (warmup, sub-intervalos por estado, reads mixtos, box = ancla)
- `include/macrodr/cmd/qdtf_likelihood.h` + comando DSL
  `calc_qdtf_likelihood(model, parameters, experiment, data, n_poles, cutoff)`
- paridad: `tests/macroir/test_qdtf_member.cpp` (box vs miembro av=2)

Compilación standalone (segundos, sin macrodr_core):
`g++ -std=c++20 -O2 tests/math/test_<X>.cpp third_party/catch2/catch_amalgamated.cpp -Ithird_party/catch2 -Ilegacy`
Chequeo semántico de los TU pesados sin build: `g++ -fsyntax-only -std=c++20
-Ilegacy -Iinclude <TU>` (qmodel.h ≈ 3 s; el TU del test ≈ 22 s).

Plan canónico: `../macroir_bessel_plan.md`.
