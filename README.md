# Tau Core Paper 28 (`v3.0.0`) — Complete Four-Terminal Cosmological Closure Beyond Gravity Alone

**Author:** József Olcsák ([ORCID: 0009-0003-7311-6522](https://orcid.org/0009-0003-7311-6522))  
**Verification:** `15/15 PASS` (`rc = 0`, `k_added = 0`, standalone Python suite [`scripts/verify_paper28_multiterminal_closure.py`](scripts/verify_paper28_multiterminal_closure.py))  
**Manuscripts:** [English PDF](Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_EN.pdf) | [Hungarian PDF](Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_HU.pdf)  
**Interactive 3D Observatory:** [`observatory/tau_core_cosmic_evolution_observatory.html`](observatory/tau_core_cosmic_evolution_observatory.html)

## Complete Closure in `v3.0.0` (All Three Remaining Open Questions Closed with `k_added = 0`)

1. **Zero-Tuning UV Kination–Leptogenesis Closure (`PROVEN`):**
   - **Single-Scale Stueckelberg–Scalaron Exit (`Sub1` $V_\mu = M_{\mathrm{inf}}A_\mu$):** With 1-loop curved-space RG slope $\kappa_J^{\mathrm{1\text{-}loop}} = 0.014914$ and $w_{\mathrm{rh}} = 1.0$, kination and decay occur at the single locked scale $M_{\mathrm{inf}} = 3.1959\times 10^{13}\;\mathrm{GeV}$ ($\rho_{\mathrm{decay}}/\rho_{\mathrm{tot}} \equiv 1$), predicting $n_s = 0.97464$ ($+0.21\sigma$ ACT DR6), $A_s = 2.0895\times 10^{-9}$ ($-0.35\sigma$), $r = 0.00511$, $\Delta N_{\mathrm{eff}}^{\mathrm{GW}} = 4.88\times 10^{-5} \ll 0.17$, and $\eta_B = +6.106\times 10^{-10}$.
   - **Decoupled Two-Fluid Bunch–Davies Proof:** Even without inter-sector coupling after $\phi_{\mathrm{end}}$, the exact de Sitter Bunch–Davies stochastic infrared variance yields $f_{\phi,\mathrm{end}}^{\mathrm{BD}} = 1.790\times 10^{-10} = 1.040\,K_B$ ($40\times$ above the $4.47\times 10^{-12}$ dilution bound with zero tuning); its $K_4$-projected condensate $\alpha_{\mathrm{OCSR}}^2 K_B = 2.242\times 10^{-11}$ gives $n_s = 0.97433$ ($+0.11\sigma$) and 4-ODE Boltzmann $\eta_B = +5.600\times 10^{-10}$ (matching Paper 24's EW crossover to $+0.27\%$).

2. **Continuous $z_{\mathrm{drag}} = 1060 \to 0$ $K_B^{-1}$ Helmholtz–Milgrom $\sigma_8$ Screening & JWST Compactness (`PROVEN`):**
   - Eliminates the `v2.0.0` step-function $z_{\mathrm{on}} = z_{\mathrm{dip}}$ via `Sub1`'s canonical apparent-horizon Raychaudhuri–Helmholtz mass $\bar{m}_{\mathrm{IR}}^2(z) = K_B^{-1}\max(\bar{R}_{\mu\nu}A^\mu A^\nu, 0)/c^2$ plus Milgrom's QUMOND center-of-mass equivalence principle for condensed galactic baryons ($R_b = 0.03\,R_{200}$).
   - Continuous 2-fluid integration yields $\sigma_8(z=0) = 0.8223$ ($+0.48\sigma$ vs. ACT DR6 CMB lensing $0.818\pm 0.009$), $100.000\%$ preservation of SPARC ($z=0$) and X-COP ($z=0.05$), and $\ell_{\mathrm{IR}}(z=10) = 2.33\;\mathrm{kpc}$ ($S_{\mathrm{Helm}}(1.5\;\mathrm{kpc}) = 0.900$), naturally explaining the sub-kpc compactness of JWST $z = 10\text{--}14$ galaxies.

3. **5-Probe Joint Local-Void Inverse Closure at $\Delta M_B \equiv 0$ & Automatic ISW Survival (`PROVEN`):**
   - Joint profile-likelihood minimization across **SH0ES**, **CCHP JWST TRGB/JAGB**, **CosmicFlows-4**, **CLASSIX/6dFGS** density contrast, and void radius ($\mathrm{dof} = 3$, $\Delta M_B \equiv 0$ by the Cepheid Null Theorem $|\Delta G/G| = 9.15\times 10^{-20}$) converges with $p = 0.599$ ($\chi^2/\mathrm{dof} = 0.625$, Total EFE: $\hat{\delta}_v = -0.1789, \hat{R} = 149.1\;\mathrm{Mpc} \implies H_0 = 71.93, v_{\mathrm{out}} = 423.4\;\mathrm{km\,s^{-1}}$) and $p = 0.566$ ($\chi^2/\mathrm{dof} = 0.678$, Baryon EFE: $\hat{\delta}_v = -0.1657, \hat{R} = 142.5\;\mathrm{Mpc} \implies H_0 = 72.05, v_{\mathrm{out}} = 427.6\;\mathrm{km\,s^{-1}}$).
   - Exact `Sub1` automatically satisfies the Supervoid ISW Sign-Flip No-Go bound ($\mathrm{d}\ln\mu/\mathrm{d}\ln a = +0.1992 < 0.2460$, whereas pure MOND $+0.3975$ flips supervoids into hot spots).

## Reproducing All 15 Gates & Figures

```bash
python3 scripts/verify_paper28_multiterminal_closure.py
python3 scripts/generate_paper28_figures.py
```
