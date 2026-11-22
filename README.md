# Tau Core Paper 28 (`v4.2.0`) — Four-Terminal Cosmological Audit Beyond Gravity Alone

**Title:** *Four-Terminal Cosmological Audit Beyond Gravity Alone: Single-Metric `AeST`–QUMOND Helmholtz Screening of $\sigma_8$, UV & ISW No-Go Theorems, `Pantheon+SH0ES` Curvature Audit, LTB Director Expulsion, and Boltzmann CMB $EE$ Closure*

**Author:** József Olcsák (ORCID: `0009-0003-7311-6522`)  
**Version:** `v4.2.0` (Round-4 Referee-Hardened Edition — `12/12 PASS`, `rc = 0`)

---

## Executive Summary & Four-Tier Epistemic Ledger

Paper 28 (`v4.2.0`) conducts a rigorous variational, symbolic (`SymPy` 4D Christoffel & Ricci tensor), and sparse finite-volume PDE/ODE audit across all six cosmological and strong-field sectors of the four-terminal Tau Core spine ($N = 4 \implies u = \sqrt{3}/24$; terminals $O_G, O_T, O_Q, O_M$), strictly separating parameter-free theorems (`PROVEN`, $k_{\mathrm{added}} = 0$) and mathematical exclusions (`PROVEN NO-GO`) from `CONDITIONAL` and `OPEN` frontiers:

| Gate | Sector | Terminals | Key Quantitative Result | $k_{\mathrm{added}}$ | Epistemic Class |
| :--- | :--- | :--- | :--- | :---: | :--- |
| **G1** | UV Prescription I Baseline | $O_G+O_Q+O_M$ | $N_* = 51.412$, $n_s = 0.96199$ ($-0.69\sigma$ Planck; $-4.00\sigma$ ACT DR6), $A_s = 2.1011\times 10^{-9}$, $r = 0.00426$, first-principles 5-ODE $\eta_B = +6.106\times 10^{-10}$ ($-0.35\sigma$) | `0` | `PROVEN` |
| **G2** | UV Presc. II & Pure-Kination No-Go | $O_G+O_T+O_Q$ | $y_N^2 = 2\alpha_{\mathrm{OCSR}}^2$, $\gamma_\Phi$ cancels identically: (a1) fermion $\kappa_J = -0.1953$ hilltop; (a2) Paper 24 $n_{\mathrm{STr}} = +2 \implies \kappa_J = +0.0651$, $10^9 A_s \le 1.645$ across all $f_\phi \in (0,1]$ up to $\Delta N_{\mathrm{eff}}^{\mathrm{GW}} = 0.17$; (b) gauged $U(1)_{B-L}$ $(q_\Phi^{\mathrm{GUT}})^2 = 3/2 \implies 10^9 A_s = 5.675$ ($+119.2\sigma$); pure kination ($\kappa_{\mathrm{eff}}=0$) violates COBE $A_s$ ($+31\sigma$) | `0` | `PROVEN NO-GO` |
| **G3** | 1D COBE-Locked Leptogenesis Ceiling | All 4 | Fixing $10^9 A_s \equiv 2.1000$ ($k_{\mathrm{UV,eff}}=1$), Planck $\eta_B = (6.12\pm 0.04)\times 10^{-10}$ cuts off kination at $\log_{10} f_\phi \ge -9.608$ (5 orders tighter than BBN GWs!), enforcing $n_s \le 0.9729$ ($2\sigma$); minimax $\log_{10} f_\phi = -8.80, \kappa_{\mathrm{eff}} = 0.01121 \implies n_s = 0.97216$ ($-0.61\sigma$ ACT), $\eta_B = +6.095\times 10^{-10}$ ($-0.63\sigma$), $\Delta N_{\mathrm{eff}}^{\mathrm{GW}} = 1.39\times 10^{-8}$ | `1` | `CONDITIONAL` |
| **G4** | 4D Christoffel Single-Metric `AeST`–QUMOND | $O_G+O_T$ | Single physical metric $g_{\mu\nu}$: $\Phi \equiv \Psi$ (100% lensing), $c_T = c_\gamma = c \le c_{S,1}$, $|\theta_{\Phi\varphi}|^2 = 1.85\times 10^{-21} \ll 10^{-15}$, $c_V^2 = 2c^2$, $\alpha_2^{\mathrm{PPN}} \equiv 0$, $\alpha_1^{\mathrm{PPN}} = -3.44\times 10^{-10}$, pole-free $\det\mathbf{H}_4 > 0$, radial EL $\frac{\mathrm{d}}{\mathrm{d}r}[c_2^{\mathrm{eff}}\theta_A] = 0$ | `0` | `PROVEN` |
| **G5** | First-Principles $\sigma_8$ + $\sum m_\nu = 59.19\;\mathrm{meV}$ | $O_G+O_T+O_M$ | First-principles $\sigma_8^{\mathrm{GR}} = 0.8150$ (Presc. I) / $0.8179$ (2-Field); unscreened $\sigma_8^{\mathrm{Weyl}} = 11.69$ (`NO-GO`); locked $\kappa_\theta = \xi_\Phi^2 \implies \sigma_8^{\mathrm{Weyl}} = 0.8165$ ($-0.16\sigma$ ACT DR6) / $0.8195$ ($+0.17\sigma$ on COBE-locked 2-Field) | `0` | `PROVEN` |
| **G6** | Weyl-Shielded Supervoid ISW ($|\delta_b| \equiv 0.25$) | $O_G+O_T$ | `brentq` root-shooting: unscreened `Sub1` $A_{\mathrm{ISW}}^{\mathrm{mat}} = -4.094 < 0$ (`NO-GO`); 4D Maxwell Weyl invariance shields intervening voids from `Tier-2` clock blowup; $\kappa_\theta = \xi_\Phi^2 \implies A_{\mathrm{ISW}} = +1.000$ ($+0.05\sigma$ vs. DESI Legacy linear ISW $0.984\pm 0.349$; DES Y3 $4.1\pm 1.1$ tail `OPEN`) | `0` | `PROVEN` / `OPEN` |
| **G7** | Cepheid Gravitational Nullity | $O_G+O_M$ | Inside $5\,M_\odot, 50\,R_\odot$ Cepheids, $|\Delta G/G|_{\mathrm{Sub1}} = 9.13\times 10^{-20} \implies \Delta M_B^{\mathrm{grav}} \equiv 0$ | `0` | `PROVEN` |
| **G8** | Official `Pantheon+SH0ES` (`1341x1341 STAT+SYS`) | $O_G+O_T+O_M$ | Uncompensated KBC ($\delta_v = -0.25, 300\;\mathrm{Mpc}$) excluded at $\Delta\chi^2 = +56.92$ (`PROVEN NO-GO`); Paper 23 zero-choice `H2a` ($\Delta\chi^2_0 = -1.55$) relaxes $+4.0$ ceiling to $\delta_v = -0.1055 \implies H_0^{\mathrm{loc}} = 69.17$ ($-0.63\sigma$ CCHP HST+JWST TRGB, $-0.29\sigma$ JWST TRGB; $\Delta M_B^{\mathrm{astro}} \ge 0.118\;\mathrm{mag}$ for SH0ES) | `1` | `NO-GO` / `COND.` |
| **G9–11** | 12 X-COP, $z=10$ JWST PDE, BHs & $^{229}\mathrm{Th}$ | All 4 | Exact $F_{\mathrm{OD}}^2-1$ 650-shell mass-compensated LTB PDE proves $\theta_A/(3H) = 6.19\times 10^{-8}$ at $R_{500}$ (99.99% QUMOND across 12 X-COP clusters, $\chi^2/24 = 1.142$, med. $1.0005$) & $2.84\times 10^{-5}$ in $z=10$ JWST disk (99.89% QUMOND, $t_{\mathrm{relax}}/t_{\mathrm{ff}} = 4.5\times 10^{-4}$); BH $4.31\times 10^{-11}$, $^{229}\mathrm{Th}$ $3.40\times 10^{-18}$ | `0` | `PROVEN` |
| **G12** | Absolute-$\mu\mathrm{K}^2$ Boltzmann CMB $EE$ | $O_G+O_T+O_Q$ | `RECFAST` ($F = 1.14$, He II, $T_M \implies z_* = 1089.80$) + $1{,}069$-mode $\ell_{\max} = 6$ non-reflecting Ma–Bertschinger hierarchy ($160$ LOS nodes): $\ell_{1\text{--}4}^{EE} = [141, 393, 686, 988]$ (`+2.92%, -1.50%, -0.58%, -0.70%`; peaks 2–4 within `1.5%`) and $[1.14, 21.70, 38.55, 44.17]\;\mu\mathrm{K}^2$ | `0` | `PROVEN` |

---

## Reproducibility

```bash
python3 scripts/verify_paper28_multiterminal_closure.py
python3 scripts/generate_paper28_figures.py
pdflatex -interaction=nonstopmode Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_EN.tex
pdflatex -interaction=nonstopmode Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_HU.tex
```
