# Tau Core Paper 28 (`v1.0.0`): Four-Terminal Cosmological Closure Beyond Gravity Alone

**Title:** *Four-Terminal Cosmological Closure Beyond Gravity Alone: ACT DR6 Inflation, Multi-Terminal Hubble Tension Budget, ISW Supervoid Sign-Flip Theorem, and Six-Sector Observational Audit*  
**Author:** József Olcsák (ORCID: [0009-0003-7311-6522](https://orcid.org/0009-0003-7311-6522))  
**Version:** `v1.0.0` (2026-11-21)

---

## Core Results (`8/8` Verification Gates Passed, `rc=0`)

1. **Four-Terminal ACT DR6 Inflation Closure (`CONDITIONAL`) & BBN $\Delta N_{\mathrm{eff}}$ Theorem (`PROVEN`):**
   - Coupling Jordan-frame Prescription II running ($\kappa_J = 0.01380$) to a stiff post-inflationary $O_T/O_Q$ kinetic phase ($w_{\mathrm{rh}} = 0.9015$, $\xi_H = 0$) achieves simultaneous $0.00\sigma$ closure on ACT DR6 $n_s = 0.9740$ and $A_s = 2.100 \times 10^{-9}$ ($r = 0.00502$, $N_* = 59.72$) while preserving Paper 24 leptogenesis $\eta_B = 6.106 \times 10^{-10}$ ($-0.35\sigma$).
   - Primordial blue-tilted GWs yield $\Delta N_{\mathrm{eff}}^{\mathrm{GW}} = 8.41 \times 10^{-6}$ ($4.99 \times 10^{-5}$ even at $w_{\mathrm{rh}} = 1$), more than $2 \times 10^4\times$ below the BBN + ACT DR6 limit ($< 0.17$).
2. **Four-Terminal Hubble Tension Budget (`CONDITIONAL`):**
   - Inside a modest observed local underdensity ($\delta_v = -15\%$, $R = 150\;\mathrm{Mpc}$, 50%-screened `OCSR` $\mu_{\mathrm{eff}} = 2.50$): $H_0^{\mathrm{loc}} = 68.22 + 0.42\,(O_T\text{ glob}) + 0.38\,(O_T\text{ loc}) + 3.10\,(O_G) + 1.10\,(O_M/O_Q) = \mathbf{73.22\;\mathrm{km\,s^{-1}\,Mpc^{-1}}}$ ($+0.17\sigma$ vs. SH0ES $73.04 \pm 1.04$; $72.12$ at $\Delta M_B = 0$), with $v_{\mathrm{out}} = 466\;\mathrm{km\,s^{-1}}$ matching CosmicFlows-4.
3. **ISW Supervoid Enhancement & Sign-Flip No-Go Theorem (`PROVEN`):**
   - Exact sign-flip bound: any supervoid model with $\mathrm{d}\ln\mu/\mathrm{d}\ln a > 1 - f(z) = 0.246$ ($z = 0.50$) flips supervoids from cold spots to hot spots.
   - In the EFE-saturated regime ($\mathrm{d}\ln\mu/\mathrm{d}\ln a \approx 0$), the *same* 50%-screened $\mu_{\mathrm{eff}} = 2.50$ factor predicts $\Delta T_{\mathrm{ISW}} = -6.84\;\mu\mathrm{K}$ ($A_{\mathrm{ISW}} = 2.88\times$; full `OCSR`: $4.60\times$), matching DES Y3 ($3.0\text{--}5.2$) with zero new parameters.
4. **Six-Sector Observational Audit (`PROVEN` / `CONDITIONAL`):**
   - **X-COP / Planck SZ Clusters:** $f_b^{\mathrm{Tau}} = 0.15703$ (Planck: $0.1564$), $b_{\mathrm{HSE}} \approx 16\text{--}18\%$.
   - **JWST High-$z$ ($z = 6\text{--}14$):** $O_T$ age shift is only $+0.11\text{--}0.16\%$; nonlinear $O_G$ collapse ($a_0^{\mathrm{KH}}(10) = 6.59\,a_0(0)$) drives early formation.
   - **Black Holes (EHT vs. LIGO/LISA):** $|M_{\mathrm{EHT}}/M_{\mathrm{GW}} - 1| \sim 8.61 \times 10^{-11}$.
   - **Lab Clocks ($\mathrm{Sr}/\mathrm{Yb}$ vs. $^{229}\mathrm{Th}$):** tree-level LPI $= 0$, PPN modulation ceiling $6.80 \times 10^{-18}$.
   - **CMB $EE$ Polarization ($O_Q$):** peaks at $(137.4, 398.6, 688.4, 990.4)$ vs. Planck/ACT $(137, 399, 690, 995)$.

---

## Repository Contents

- `Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_EN.tex` / `.pdf` — English manuscript (`v1.0.0`)
- `Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_HU.tex` / `.pdf` — Hungarian manuscript (`v1.0.0`)
- `scripts/verify_paper28_multiterminal_closure.py` — Standalone verification suite (`python3 scripts/verify_paper28_multiterminal_closure.py`)
- `scripts/generate_paper28_figures.py` — Generates all 3 vector PDF figures in `figures/`
- `paper28_verification_results.json` — Machine-verified numerical ledger
- `observatory/tau_core_cosmic_evolution_observatory.html` — Interactive GPU-accelerated 3D Tau Core Observatory (70,460 GPU quanta, 4 Physical Lenses, live simulation-computed Mollweide CMB Atlas, and `TT`/`EE`/`TE` spectra)
