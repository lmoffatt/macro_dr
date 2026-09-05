---
date: 2026-08-31
status: ecuaciones implementadas en bessel_oracle.py; estado de gates al pie
scope: los bloques de predicción/update de la recursión aumentada (miembro Bessel),
  a escala de conteos, tal como los debe transcribir el C++ (calc_Qdtf_eig + la
  recursión hermana de qmodel.h:4519)
lee: bessel_oracle.py (implementación); macroir_bessel_plan.md (plan);
  ../mr_vs_ir_from_macror.md (la base Σ^bnd del ancla box)
---

# Bloques de la recursión aumentada (miembro Bessel)

## Objetos

Sistema de observación REAL (A, b, read), genérico:

- nativo (K+4): A = bloques modales 2×2 del Bessel ([[−σ,−ω],[ω,−σ]] por par
  conjugado), b = (2Re r, 2Im r), read = c_f (1,0 por par)
- agrupado (K+5): fila cumulador du/dt = c_fᵀφ anexada; read = u/Δ; u, sus
  filas/columnas de covarianza y Cov(N,u) se ponen a CERO al inicio de cada
  ventana grabada (reset legítimo del grabador)
- box (K+1): integrador solo sobre i(t); ES el miembro actual (ancla)

Estado transportado (escala conteos): μ_N (K), S_NN (K×K), m_x (m),
C_Nx (K×m) = Cov(N, x), S_xx (m×m).

Cantidades de ventana (por cuadratura en el oráculo; espectrales con las
divided differences corridas en producción):

- P = e^{QΔ}
- F1[i,j] = E[1{X_Δ=j}·w | X₀=i] ∈ ℝ^m  (primeros PAR-RESUELTOS, por modo),
  con w = ∫₀^Δ e^{A(Δ−s)} b γ_{X_s} ds el aporte de un canal al filtro
- f1[i] = Σ_j F1[i,j]  (primeros condicionados al inicio)
- S2[i] = E[w wᵀ | X₀=i] ∈ ℝ^{m×m}  (segundos SOLO marginalizados al inicio)
- Φ = e^{AΔ}  (descuento del estado viejo; diagonal por bloques 2×2)
- Σ_ξ = S₀·∫₀^Δ e^{As} b bᵀ e^{Aᵀs} ds  (Gramiano del ruido, Van Loan)

## Predicción

    μ_N⁺  = Pᵀ μ_N
    S_NN⁺ = Pᵀ (S_NN − diag μ_N) P + diag(μ_N⁺)          [= σ_pre, qmodel:4655]

    m_x⁺  = Φ m_x + Σ_i μ_N[i]·f1[i]  (+ término de baseline)

    C_Nx⁺[j] = (Pᵀ C_Nx Φᵀ)[j]                            [acarreo del viejo]
             + Σ_i μ_N[i]·(F1[i,j] − P_ij·f1[i])          [dentro de canal]
             + Σ_{i,k} P_ij·S_NN[i,k]·f1[k]               [entre canales]
      — las dos últimas líneas son EXACTAMENTE la forma de gS (qmodel:4637-4645)
        una vez por modo: gS^{(p)} = ḡ^{(p)ᵀ}·SmD·P + μ·gtotal^{(p)}

    S_xx⁺ = Φ S_xx Φᵀ + Φ M1 + (Φ M1)ᵀ + V_w + Σ_ξ
      M1  = Σ_i outer(C_Nx[i], f1[i])                      [Cov(x₀, Σw)]
      V_w = Σ_i μ_N[i]·(S2[i] − f1[i]f1[i]ᵀ)              [E de la Var]
          + Σ_{i,k} S_NN[i,k]·f1[i]f1[k]ᵀ                 [Var de la E]

Todo afín o bilineal en los momentos transportados: propagación exacta.
Los canales son independientes DADOS sus estados de inicio; toda la
correlación entre canales entra por S_NN y C_Nx.

## Update (condicionamiento gaussiano, un escalar z por ventana)

    ẑ = readᵀ m_x⁺        v = readᵀ S_xx⁺ read (+ piso σ_fl²)      δ = z − ẑ
    κ_N = C_Nx⁺·read      κ_x = S_xx⁺·read

    μ_N ← μ_N⁺ + κ_N δ/v        m_x ← m_x⁺ + κ_x δ/v
    S_NN ← S_NN⁺ − κ_Nκ_Nᵀ/v    C_Nx ← C_Nx⁺ − κ_Nκ_xᵀ/v    S_xx ← S_xx⁺ − κ_xκ_xᵀ/v

    logL += −½(log 2πv + δ²/v)

## Ancla box (gate O2)

Con (A,b) = integrador y read = u/Δ, estos bloques deben coincidir con el
ensamblaje de frontera de mr_vs_ir_from_macror.md a escala de conteos:

    Cov(B)_{(ij),(i'j')} = P_ij (S − diag μ)_{ii'} P_{i'j'} + δ_{ii'}δ_{jj'} μ_i P_ij
    ŷ = Σ μ_i P_ij Γ̄_ij     v = Γ̄ᵀCov(B)Γ̄ + Σ μ_i P_ij V̄_ij + S₀/Δ
    κ_j = Σ P_ij (S−diagμ) ḡ⁰ + Σ μ_i P_ij Γ̄_ij

con Γ̄_ij = F1[i,j]/(P_ij Δ) y V̄ la varianza residual del par. La igualdad es
algebraica (verificación simbólica en el diseño: los términos μPΓ̄² se cancelan
entre las dos contabilidades); el gate la chequea numéricamente con la misma
grilla de cuadratura en ambos lados para que solo quede el álgebra.

## Qué consume el C++ y qué no

- primeros: par-resueltos, UNA tabla por modo (F1) — clase E2/Ee corrida
- segundos: SOLO marginalizados al inicio (S2) — colapso e^{(Q−νI)t}𝟙 = e^{−νt}𝟙
- ningún objeto par-resuelto de segundos momentos (IRT excluido por decisión,
  2026-08-31, conversación "E[V̄] residual IRT no se usa")

## Verificación

Gates del oráculo (bessel_oracle.py): O1 ruido del simulador vs formas
cerradas; O2 ancla box vs Σ^bnd (dos ventanas, con update en el medio);
O3 predictivo de ventana 1 vs Monte Carlo (lectura nativa); O4 blancura de
residuales estandarizados + lag-1 positivo del miembro box sobre verdad
filtrada (los números de motivación); O5 ídem para la lectura agrupada (K+5).

ESTADO 2026-08-31 (3): el miembro está INTEGRADO a macro_dr:
`legacy/qdtf_member.h` (loop sobre el registro: warmup de estacionariedad del
filtro, composición de sub-intervalos por estado sin monoide, reads
nativo/agrupado por ventana con un solo layout, baseline como offset por
H(0)=1, NaN = predict-only, S₀ = Current_Noise, Pink como piso,
Proportional NO heredado, N_ch variable rechazado explícitamente) +
`include/macrodr/cmd/qdtf_likelihood.h` + comando DSL
`calc_qdtf_likelihood(model, parameters, experiment, data, n_poles, cutoff)`
(n_poles=0 = box = ancla) + test de paridad
`tests/macroir/test_qdtf_member.cpp` (box vs miembro av=2 a 1e-4 relativo;
tolerancia porque el miembro establecido lleva shrinkage ε_mach·κ_V y α de
confianza que el motor no). Verificado por -fsyntax-only con instanciación
forzada (14-22 s por TU, solo warnings preexistentes); la corrida del test de
paridad requiere el build completo (Luciano).

ESTADO 2026-08-31 (2): estas ecuaciones ya están TRANSCRITAS a C++ en
`legacy/qdtf_engine.h` (productor espectral + predict/update), con sus propios
gates verdes (`tests/math/test_qdtf_engine.cpp`, 668 aserciones) y cruce
numérico exacto contra el oráculo. Detalle numérico hallado ahí: el ensamblaje
modal del Gramiano de ruido cancela ~4 dígitos a las magnitudes de residuos
del Bessel (queda ~1e-8 relativo en double; vigilar PSD si algún día se
suman más polos). Falta SOLO el adaptador en qmodel.h (tipos Qdtf +
MacroR2/dispatch + DSL).

ESTADO 2026-08-31 (1): ALL GATES PASS (corrida en tmp/bessel_oracle_run.txt).
O1: Var 0.18% y Cov₁ 0.008% vs formas cerradas (lectura nativa: R(0)/R(Δ),
NO las formas boxcar — trampa de arnés detectada y corregida en esta fecha).
O2: |Δŷ|, |Δv|, |Δκ| ≤ 2e-14 en dos ventanas (con update en el medio).
O3: media 1.2 SE, varianza 1.4 SE (n=3000). O4: residuales del miembro
mean −0.003, var 1.008, lag1 −0.0000 (SE 0.0019); el box sobre la misma
verdad filtrada: var(r) 0.43 y lag1 +0.596 — los números de motivación del
paper a f_c·Δ = 0.2. O5 (K+5): media 0.3 SE, var 0.1 SE, lag1 +0.002 (SE
0.005). Los kernels C++ pasaron sus 375 aserciones aparte (test_acquisition_
filter.cpp), incluida la fila del apéndice del manuscrito (61.5%/0.641).
