# Tau Core Paper 28 (`v2.0.0`) — Four-Terminal Cosmological Closure Beyond Gravity Alone

**Author:** József Olcsák ([ORCID: 0009-0003-7311-6522](https://orcid.org/0009-0003-7311-6522))  
**Verification:** `12/12 PASS` (`rc = 0`, standalone Python suite [`scripts/verify_paper28_multiterminal_closure.py`](scripts/verify_paper28_multiterminal_closure.py))  
**Manuscripts:** [English PDF](Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_EN.pdf) | [Hungarian PDF](Tau_Core_Paper_28_Multi_Terminal_ACT_Hubble_ISW_HU.pdf)  
**Interactive 3D Observatory:** [`observatory/tau_core_cosmic_evolution_observatory.html`](observatory/tau_core_cosmic_evolution_observatory.html)

## What Tightened in `v2.0.0` (All Three Conditional Sectors Narrowed from First Principles)

1. **Zero-Parameter (`k_added = 0`) ACT DR6 Inflation & Two-Fluid Dilution Bound:**
   - Replaced the `v1.0.0` 2-parameter $(\kappa_J, w_{\mathrm{rh}})$ root-fit with the exact 1-loop curved-space RG slope for Paper 24's complex $\mathrm{Spin}(10)$ $B-L$ seesaw radial modulus ($y_{\mathrm{eff}}^2 = \alpha_{\mathrm{OCSR}}^2/2$, $n_{\mathrm{eff}}^{\mathrm{STr}} = 2$, $d_S = 3$): $\kappa_J^{\mathrm{1\text{-}loop}} = 0.014914$.
   - Paired with pure $O_T$ director kination ($w_{\mathrm{rh}} = 1.0$, $k_{\mathrm{added}} = 0$): $N_* = 60.175$, $n_s = 0.97464$ ($+0.21\sigma$ vs. ACT DR6), $A_s = 2.0895\times 10^{-9}$ ($-0.35\sigma$), $r = 0.00511$, and $\Delta N_{\mathrm{eff}}^{\mathrm{GW}} = 4.88\times 10^{-5} \ll 0.17$ (`PROVEN`).
   - Proved the **Two-Fluid Kination–Leptogenesis Dilution Bound**: independent scalaron + kination fluids require $1 - f_{\mathrm{kin,end}} \ge \sqrt{\rho_{\mathrm{rh}}/\rho_{\mathrm{end}}} = 4.47\times 10^{-12}$ to avoid diluting $\eta_B = 6.106\times 10^{-10}$, pinning the non-tuned branch to a coupled quintessential-kination exit (`CONDITIONAL`).

2. **Exact `Sub1` Local-Void Hubble Budget (`Delta M_B = 0`, Zero Ad-Hoc 50% Factor) & Cepheid Null Theorem:**
   - Eliminated the `v1.0.0` ad-hoc 50% screening factor completely via exact `Sub1` Eq. (82–83) baryon-only QUMOND + GR aether-dust ($\mu_{\mathrm{Sub1}} = (1-f_b) + f_b\,\mu_{\mathrm{QUMOND}}^{\mathrm{Gauss}}$, $f_b = 0.15703$) with cosmic external field $g_{\mathrm{ext}}(R, z)$ integrated from the CMB-normalized linear power spectrum $P(k, z)$.
   - Proved the **Cepheid / SN Ia Gravitational Null Theorem**: $|\Delta G/G|_{\mathrm{Cepheid}} = 9.15\times 10^{-20} \implies \Delta M_B^{\mathrm{grav}} \equiv 0$.
   - Audited five published local-void catalogs strictly at $\Delta M_B = 0$ and $z_{\mathrm{on}} = z_{\mathrm{dip}} = 0.6316$: **CLASSIX X-ray** ($\delta_v = -0.20, R = 140\;\mathrm{Mpc}$, Total EFE) hits $H_0^{\mathrm{loc}} = 72.30\;\mathrm{km\,s^{-1}\,Mpc^{-1}}$ ($-0.71\sigma$ vs. SH0ES) and $v_{\mathrm{out}} = 442.0\;\mathrm{km\,s^{-1}}$ ($+0.64\sigma$ vs. CosmicFlows-4) simultaneously; **Baseline** ($\delta_v = -0.15, R = 150\;\mathrm{Mpc}$, Baryonic EFE) hits $H_0^{\mathrm{loc}} = 71.72$ and $v_{\mathrm{out}} = 405.4\;\mathrm{km\,s^{-1}}$ ($-0.38\sigma$ vs. CF4).

3. **Linear $\sigma_8$ Screening No-Go, Automatic `Sub1` ISW Sign-Flip Survival & Exact JWST `Sub1` Boost:**
   - Proved the **Linear $\sigma_8$ Screening No-Go Theorem**: unscreened linear QUMOND from $z_{\mathrm{on}}=1000$ drives $\sigma_8 = 13.075$ (and $1.458$ at $z_{\mathrm{on}}=10$), whereas locking activation to `Sub1`'s Cartan inflection $z_{\mathrm{dip}} = 0.6316$ bounds $\sigma_8 = 0.8493$ ($+4.7\%$).
   - Proved **Automatic `Sub1` ISW Sign-Flip Survival**: $\mathrm{d}\ln\mu_{\mathrm{Sub1}}/\mathrm{d}\ln a = +0.1992 < 0.2460$ (linear) and $< 0.2963$ (nonlinear void), whereas pure unscreened MOND ($+0.3974$) flips supervoids into hot spots ($-3.93\times$), establishing a sharp dynamic ($0.50\text{--}0.86\times$) vs. EFE-saturated ($2.20\text{--}2.61\times$) ISW discriminator.
   - Replaced the ad-hoc JWST $\mu=1.5$ ansatz with exact `Sub1` $\mu_{\mathrm{Sub1}}(z=10) = 1.612\text{--}2.153$ alongside full six-sector closure (X-COP clusters, EHT/LIGO BHs, $\mathrm{Sr}/^{229}\mathrm{Th}$ clocks, CMB $EE$ peaks).

## Reproducing All Results & Figures

```bash
python3 scripts/verify_paper28_multiterminal_closure.py
python3 scripts/generate_paper28_figures.py
```
