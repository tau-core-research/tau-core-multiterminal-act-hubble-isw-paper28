# Tau Core Paper 28 (`v4.1.0`) — Four-Terminal Cosmological Audit Beyond Gravity Alone

**Title:** *Four-Terminal Cosmological Audit Beyond Gravity Alone: Single-Metric `AeST`–QUMOND Helmholtz Screening of $\sigma_8$, UV Prescription II & Unscreened ISW No-Go Theorems, Redshift-Windowed Local-Void $H_0$ Budget, and First-Principles Boltzmann CMB $EE$ Closure*

**Author:** József Olcsák (ORCID: `0009-0003-7311-6522`)  
**Version:** `v4.1.0` (Round-3 Referee-Certified Edition — `12/12 PASS`, `rc = 0`)

---

## Executive Summary & Four-Tier Epistemic Ledger

Paper 28 (`v4.1.0`) conducts a rigorous variational, symbolic (`SymPy` 4D Christoffel & Ricci tensor), and numerical ODE audit across all six cosmological and strong-field sectors of the four-terminal Tau Core spine ($N = 4 \implies u = \sqrt{3}/24$; terminals $O_G, O_T, O_Q, O_M$), strictly separating parameter-free theorems (`PROVEN`, $k_{\mathrm{added}} = 0$) and mathematical exclusions (`PROVEN NO-GO`) from `CONDITIONAL` and `OPEN` frontiers:

| Gate | Sector | Terminals | Key Quantitative Result | $k_{\mathrm{added}}$ | Epistemic Class |
| :--- | :--- | :--- | :--- | :---: | :--- |
| **G1** | UV Prescription I Baseline | $O_G+O_Q+O_M$ | $N_* = 51.412$, $n_s = 0.96199$ ($-0.69\sigma$ Planck; $-4.00\sigma$ ACT DR6), $A_s = 2.1011\times 10^{-9}$, $r = 0.00426$, first-principles 5-ODE $\eta_B = +6.106\times 10^{-10}$ | `0` | `PROVEN` |
| **G2** | UV Prescription II 3-Branch CW & 1-Field Kination | $O_G+O_T+O_Q$ | $y_N^2 = 2\alpha_{\mathrm{OCSR}}^2$, $\gamma_\Phi$ cancels identically: (a1) fermion-only $\kappa_J = -0.1953$ hilltop collapse ($x_{\max} = 0.0161 < x_*$); (a2) Paper 24 threshold multiplet $n_{\mathrm{STr}} = +2 \implies \kappa_J = +0.0651$, $10^9 A_s = 1.349$ ($-25.05\sigma$); (b) gauged $U(1)_{B-L}$ $(q_\Phi^{\mathrm{GUT}})^2 = 3/2 \implies \kappa_J = +1.0838$, $10^9 A_s = 5.675$ ($+119.2\sigma$) | `0` | `PROVEN NO-GO` |
| **G3** | Conditional 2-Field 5-ODE & Trajectory GW | All 4 | $(\kappa_{\mathrm{eff}}, f_{\phi,\mathrm{end}}) = (0.0149, 2.24\times 10^{-11}) \implies w_{\mathrm{eff}} = 0.8737$, $n_s = 0.97432$, $A_s = 2.059\times 10^{-9}$, first-principles 5-ODE $\eta_B = +5.600\times 10^{-10}$, two-stage trajectory $\Delta N_{\mathrm{eff}}^{\mathrm{GW}} = 3.74\times 10^{-6} \ll 0.17$ | `2` | `CONDITIONAL` |
| **G4** | 4D Christoffel Single-Metric `AeST`–QUMOND | $O_G+O_T$ | Single physical metric $g_{\mu\nu}$ ($S_{\mathrm{SM}}[g_{\mu\nu}, \Psi_{\mathrm{SM}}]$, constant $K_B^{(0)}$, $K_{\mathrm{eff}} = K_B^{(0)}/2$): $\Phi \equiv \Psi$ (100% lensing), $c_T = c_\gamma = c$, $\alpha_2^{\mathrm{PPN}} \equiv 0$, $\alpha_1^{\mathrm{PPN}} = -3.44\times 10^{-10}$, zero static drift, pole-free $\det\mathbf{H}_4 \propto k_c^4(4k_c^4 + 4k_c^2 m_H^2\mu_b + m_H^4) > 0$ ($\forall\mu_b \ge 0.9412 > 0$) | `0` | `PROVEN` |
| **G5** | First-Principles $\sigma_8$ + $\sum m_\nu = 59.19\;\mathrm{meV}$ | $O_G+O_T+O_M$ | First-principles $\sigma_8^{\mathrm{GR}} = 0.8150$ (Presc. I) / $0.8106$ (2-Field); unscreened $\sigma_8^{\mathrm{Weyl}} = 11.69$ (`NO-GO`); locked $\kappa_\theta = \xi_\Phi^2 \implies \sigma_8^{\mathrm{Weyl}} = 0.8165$ ($-0.16\sigma$ ACT DR6; $0.8122$ on 2-Field) | `0` | `PROVEN` |
| **G6** | Root-Shot Supervoid ISW ($|\delta_b| \equiv 0.25$) | $O_G+O_T$ | `brentq` root-shooting + separated force/potential Helmholtz kernels: unscreened `Sub1` $f_{\mathrm{eff}} = 0.9783 \implies A_{\mathrm{ISW}}^{\mathrm{mat}} = -4.094$, $A_{\mathrm{ISW}}^{\mathrm{gal}} = -1.394 < 0$ (`NO-GO`); $\kappa_\theta = \xi_\Phi^2$ preserves cold spots ($A_{\mathrm{ISW}} = +1.000$; DES Y3 $4.1\pm 1.1$ `OPEN`) | `0` | `PROVEN` / `OPEN` |
| **G7** | Cepheid Gravitational Nullity | $O_G+O_M$ | Inside $5\,M_\odot, 50\,R_\odot$ Cepheids, $|\Delta G/G|_{\mathrm{Sub1}} = 9.13\times 10^{-20} \implies \Delta M_B^{\mathrm{grav}} \equiv 0$ | `0` | `PROVEN` |
| **G8** | Windowed Local-Void $H_0$ ($z \in [0.023, 0.15]$) | $O_G+O_T+O_M$ | KBC ($300\;\mathrm{Mpc}, -0.25$): pure Eq. (11) $O_G$ outflow $69.63$ ($+0.61$ post-`2M++` $z > 0.067$); adding radially windowed `Tier-2` differential clock shift $+0.42 \implies H_0^{\mathrm{loc}} = 70.05$ ($-0.19\sigma$ CCHP JWST TRGB; $\Delta M_B^{\mathrm{astro}} = 0.091\;\mathrm{mag}$ for SH0ES) | `2` | `CONDITIONAL` |
| **G9–11** | 12 Archival X-COP, JWST, BHs & $^{229}\mathrm{Th}$ | All 4 | All 12 Eckert+ (2019) `X-COP` clusters ($\chi^2/24 = 1.142$, median $M_{\mathrm{pred}}/M_{\mathrm{obs}} = 1.0005$ at $R_{500}$); $a_0^{\mathrm{KH}}(10)/a_0^{\mathrm{KH}}(0) = 6.605$ (radial collapse `CONDITIONAL / OPEN`); BH $4.31\times 10^{-11}$, $^{229}\mathrm{Th}$ $3.40\times 10^{-18}$ | `0` | `PROVEN` / `OPEN` |
| **G12** | Absolute-$\mu\mathrm{K}^2$ Boltzmann CMB $EE$ | $O_G+O_T+O_Q$ | `RECFAST` ($F = 1.14$, He II, $T_M \implies z_* = 1089.80$) + $1{,}069$-mode $\ell_{\max} = 6$ non-reflecting Ma–Bertschinger hierarchy ($160$ LOS nodes): $\ell_{1\text{--}4}^{EE} = [141, 393, 686, 988]$ (`+2.92%, -1.50%, -0.58%, -0.70%`; peaks 2–4 within `1.5%`) and $[1.14, 21.70, 38.55, 44.17]\;\mu\mathrm{K}^2$ | `0` | `PROVEN` |

---

## Reproducibility

```bash
python3 scripts/verify_paper28_multiterminal_closure.py
python3 scripts/generate_paper28_figures.py
pdflatex -interaction=nonstopmode Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_EN.tex
pdflatex -interaction=nonstopmode Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_HU.tex
```
