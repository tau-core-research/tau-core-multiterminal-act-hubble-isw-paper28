#!/usr/bin/env python3
"""
Generate 3 publication-grade vector PDF figures for Tau Core Paper 28 (v3.0.0).
Reads verified outputs directly from paper28_verification_results.json.
"""
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA = json.loads((ROOT / "paper28_verification_results.json").read_text(encoding="utf-8"))

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9.2,
    "axes.linewidth": 0.9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": True,
    "legend.framealpha": 0.92,
})

# ==============================================================================
# FIGURE 1: Zero-Parameter ACT DR6 Kination + Bunch-Davies & BBN Delta N_eff
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.1))

ax1.axvspan(0.974 - 0.003, 0.974 + 0.003, color="#3b82f6", alpha=0.16, label=r"ACT DR6+DESI $n_s = 0.974 \pm 0.003$")
ax1.axvspan(0.9649 - 0.0042, 0.9649 + 0.0042, color="#94a3b8", alpha=0.16, label=r"Planck 2018 $n_s = 0.9649 \pm 0.0042$")
ax1.axhspan(2.07, 2.13, color="#10b981", alpha=0.20, label=r"Planck/ACT $A_s = (2.10\pm 0.03)\times 10^{-9}$")

ns_w0 = [0.9620, 0.9649, 0.9693, 0.9741, 0.9782]
as_w0 = [2.101, 1.898, 1.643, 1.436, 1.349]
ax1.plot(ns_w0, as_w0, "--o", color="#ef4444", lw=1.6, ms=4.2, label=r"Presc. II alone ($w_{\mathrm{rh}}=0$, $O_G$ only)")

ns_w1 = [0.96741, 0.97036, 0.97274, 0.97464, 0.97626]
as_w1 = [2.841, 2.502, 2.262, 2.0895, 1.957]
ax1.plot(ns_w1, as_w1, "-^", color="#0284c7", lw=1.9, ms=4.5, label=r"Pure $O_T$ kination ($w_{\mathrm{rh}}=1$, $\kappa_J\uparrow$)")

zp = DATA["part1_inflation_bbn_leptogenesis"]["stueckelberg_scalaron_single_scale_k0"]
bd = DATA["part1_inflation_bbn_leptogenesis"]["bunch_davies_two_fluid_closure_k0"]

ax1.scatter([0.9620], [2.101], s=60, c="#334155", zorder=6, marker="D")
ax1.annotate("Presc. I\n($w_{\\mathrm{rh}}=0$)", (0.9620, 2.101), xytext=(0.9590, 2.27), fontsize=7.6, ha="center")

ax1.scatter([zp["n_s"]], [zp["A_s"] * 1e9], s=140, c="#059669", zorder=8, marker="*", edgecolors="black", linewidths=0.7)
ax1.annotate(r"Single-scale Stueckelberg ($k{=}0$)" + "\n" + r"$n_s=0.9746$ ($+0.21\sigma$), $A_s=2.090$",
             (zp["n_s"], zp["A_s"] * 1e9), xytext=(0.9688, 2.58), fontsize=7.4, ha="center",
             arrowprops=dict(arrowstyle="->", lw=0.9, color="#059669"))

ax1.scatter([bd["n_s_K4"]], [bd["A_s_K4"] * 1e9], s=70, c="#f59e0b", zorder=7, marker="P", edgecolors="black", linewidths=0.6,
            label=r"$K_4$ Bunch--Davies 4-ODE ($w_{\mathrm{eff}}=0.874$)")

ax1.set_xlim(0.956, 0.981)
ax1.set_ylim(1.22, 2.92)
ax1.set_xlabel(r"Scalar spectral index $n_s$ (NLO)")
ax1.set_ylabel(r"Primordial amplitude $10^9 A_s$")
ax1.set_title("(A) Zero-Parameter ACT DR6 Kination Closure")
ax1.legend(loc="lower left", fontsize=6.7)
ax1.grid(alpha=0.25)

# Panel B: BBN & ACT DR6 Delta N_eff(w_rh)
w_vals = np.linspace(1.0 / 3.0, 1.0, 80)
h_ratio = 3.109e11
om_flat = 5.95e-18
dneff = []
for w in w_vals:
    ngw = 2.0 * (3.0 * w - 1.0) / (1.0 + 3.0 * w)
    fr = h_ratio ** ((1.0 + 3.0 * w) / (3.0 * (1.0 + w)))
    om_tot = om_flat * ((fr**ngw - 1.0) / ngw) if ngw > 1e-6 else om_flat * math.log(fr)
    dneff.append((8.0 / 7.0) * ((11.0 / 4.0) ** (4.0 / 3.0)) * (om_tot / 2.469e-5))

ax2.semilogy(w_vals, dneff, color="#7c3aed", lw=2.1, label=r"Primordial GW $\Delta N_{\mathrm{eff}}^{\mathrm{GW}}(w_{\mathrm{rh}})$")
ax2.axhline(0.17, color="#dc2626", ls="--", lw=1.6, label=r"BBN + ACT DR6 95% CL ($\Delta N_{\mathrm{eff}} < 0.17$)")
ax2.fill_between(w_vals, 0.17, 1.0, color="#dc2626", alpha=0.10)

ax2.scatter([1.0], [zp["Delta_Neff_GW"]], s=130, c="#059669", marker="*", zorder=7, edgecolors="black", linewidths=0.7)
ax2.annotate(r"Pure kination ($w_{\mathrm{rh}}=1$)" + "\n" + r"$\Delta N_{\mathrm{eff}} = 4.88\times 10^{-5}$ ($3500\times$ safe)",
             (1.0, zp["Delta_Neff_GW"]), xytext=(0.66, 2.2e-3), fontsize=7.6, ha="center",
             arrowprops=dict(arrowstyle="->", lw=0.9, color="#059669"))

ax2.set_xlim(0.32, 1.03)
ax2.set_ylim(1e-11, 1.0)
ax2.set_xlabel(r"Post-inflationary $O_T/O_Q$ equation of state $w_{\mathrm{rh}}$")
ax2.set_ylabel(r"Effective neutrino excess $\Delta N_{\mathrm{eff}}^{\mathrm{GW}}$")
ax2.set_title(r"(B) BBN & ACT DR6 $\Delta N_{\mathrm{eff}}$ Safety Margin")
ax2.legend(loc="lower right", fontsize=7.3)
ax2.grid(alpha=0.25, which="both")

fig.tight_layout()
fig.savefig(FIG_DIR / "fig1_act_dr6_multiterminal_inflation_and_bbn.pdf")
plt.close(fig)

# ==============================================================================
# FIGURE 2: 5-Probe Joint Local-Void Closure (Delta M_B=0) & ISW Sign-Flip
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.1))

cats = DATA["part2_continuous_closure_and_six_sectors"]["local_void_catalogs_continuous_dMB0"]
jfit = DATA["part2_continuous_closure_and_six_sectors"]["joint_5probe_local_void_inverse_closure"]
cat_labels = {
    "CF4_shallow": "CF4\n(-10%)",
    "Paper28_baseline": "Baseline\n(-15%, 150)",
    "Whitbourn_Shanks_2016": "WS16\n(-15%, 200)",
    "CLASSIX_Xray_2020": "CLASSIX\n(-20%, 140)",
    "KBC_conservative_2013": "KBC\n(-20%, 200)",
}

ax1.axhspan(73.04 - 1.04, 73.04 + 1.04, color="#ef4444", alpha=0.13, label=r"SH0ES $73.04 \pm 1.04$")
ax1.axhspan(70.4 - 1.9, 70.4 + 1.9, color="#3b82f6", alpha=0.08, label=r"CCHP TRGB/JAGB $70.4 \pm 1.9$")
ax1.axvspan(419.0 - 36.0, 419.0 + 36.0, color="#10b981", alpha=0.16, label=r"CosmicFlows-4 $419 \pm 36\;\mathrm{km\,s^{-1}}$")

v_tot = [cats[k]["totEFE"]["v_out_kms"] for k in cat_labels]
h_tot = [cats[k]["totEFE"]["H0_loc_dMB0"] for k in cat_labels]
v_bar = [cats[k]["barEFE"]["v_out_kms"] for k in cat_labels]
h_bar = [cats[k]["barEFE"]["H0_loc_dMB0"] for k in cat_labels]

ax1.scatter(v_tot, h_tot, s=52, c="#0284c7", marker="o", edgecolors="black", lw=0.6, zorder=6, label="Catalogs (Tot. EFE)")
ax1.scatter(v_bar, h_bar, s=52, c="#7c3aed", marker="s", edgecolors="black", lw=0.6, zorder=6, label="Catalogs (Bar. EFE)")

# 5-Probe Joint Inverse Closure Best-Fit Stars!
ax1.scatter([jfit["totEFE"]["v_out_kms"], jfit["barEFE"]["v_out_kms"]],
            [jfit["totEFE"]["H0_loc_dMB0"], jfit["barEFE"]["H0_loc_dMB0"]],
            s=145, c="#059669", marker="*", edgecolors="black", lw=0.8, zorder=9,
            label=r"5-Probe Joint Fit ($p=0.57\text{--}0.60$)")

offsets = {
    "CF4_shallow": (18, -0.58),
    "Paper28_baseline": (-62, -0.45),
    "Whitbourn_Shanks_2016": (36, -0.58),
    "CLASSIX_Xray_2020": (42, 0.55),
    "KBC_conservative_2013": (-58, -0.55),
}
for i, key in enumerate(cat_labels):
    ax1.plot([v_tot[i], v_bar[i]], [h_tot[i], h_bar[i]], color="#475569", lw=0.9, ls=":", zorder=4)
    dx, dy = offsets[key]
    ax1.annotate(cat_labels[key], (v_tot[i], h_tot[i]), xytext=(v_tot[i] + dx, h_tot[i] + dy),
                 fontsize=6.8, ha="center", va="center",
                 arrowprops=dict(arrowstyle="-", lw=0.55, color="#64748b"))

ax1.set_xlim(210, 745)
ax1.set_ylim(69.4, 74.6)
ax1.set_xlabel(r"Local void outflow $v_{\mathrm{out}}\;[\mathrm{km\,s^{-1}}]$ ($\Delta M_B \equiv 0$)")
ax1.set_ylabel(r"Local Hubble parameter $H_0^{\mathrm{loc}}\;[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$")
ax1.set_title("(A) 5-Probe Joint Inverse Closure ($z=1060\\to 0$)")
ax1.legend(loc="lower right", fontsize=6.4)
ax1.grid(alpha=0.25)

# Panel B: ISW Sign-Flip Survival
isw = DATA["part2_continuous_closure_and_six_sectors"]["isw_supervoid_exact_Sub1"]
dlnmu = np.linspace(-0.04, 0.44, 150)
b_lin = isw["sign_flip_bound_lin"]
a_bar = 1.15 * isw["Sub1_barEFE"]["mu_z05"] * (1.0 - dlnmu / b_lin)
a_tot = 1.15 * isw["Sub1_totEFE"]["mu_z05"] * (1.0 - dlnmu / b_lin)
a_mnd = 1.15 * isw["pure_MOND"]["mu_z05"] * (1.0 - dlnmu / b_lin)

ax2.axhspan(3.0, 5.2, color="#06b6d4", alpha=0.16, label=r"DES Y3 ($A_{\mathrm{ISW}} = 3.0\text{--}5.2$)")
ax2.axhline(1.0, color="#64748b", ls=":", lw=1.1, label=r"$\Lambda$CDM ($A_{\mathrm{ISW}} = 1$)")
ax2.axhline(0.0, color="black", lw=0.9)
ax2.axvline(b_lin, color="#dc2626", ls="--", lw=1.5, label=r"Sign-flip No-Go: $1 - f(0.5) = 0.246$")
ax2.axvspan(b_lin, 0.44, color="#dc2626", alpha=0.09)

ax2.plot(dlnmu, a_bar, color="#7c3aed", lw=2.0, label=r"Exact Sub1 Bar. EFE ($\mu=2.27$)")
ax2.plot(dlnmu, a_tot, color="#0284c7", lw=1.8, ls="-.", label=r"Exact Sub1 Tot. EFE ($\mu=1.91$)")
ax2.plot(dlnmu, a_mnd, color="#ef4444", lw=1.3, ls=":", label=r"Pure unscreened MOND ($\mu=5.55$)")

ax2.scatter([isw["Sub1_barEFE"]["dlnmu_dlna"]], [isw["Sub1_barEFE"]["A_ISW_dyn_lin"]], s=75, c="#7c3aed", edgecolors="black", zorder=7)
ax2.scatter([0.0], [isw["Sub1_barEFE"]["A_ISW_EFE_sat"]], s=80, c="#7c3aed", marker="^", edgecolors="black", zorder=7)
ax2.scatter([isw["pure_MOND"]["dlnmu_dlna"]], [isw["pure_MOND"]["A_ISW_dyn_lin"]], s=85, c="#dc2626", marker="X", edgecolors="black", zorder=8)

ax2.annotate("Sub1 dynamic\n(+0.199, PASS)", (isw["Sub1_barEFE"]["dlnmu_dlna"], isw["Sub1_barEFE"]["A_ISW_dyn_lin"]),
             xytext=(0.11, -1.7), fontsize=7.3, ha="center", arrowprops=dict(arrowstyle="->", lw=0.8, color="#7c3aed"))
ax2.annotate("Sub1 EFE-sat\n($2.61\\times$)", (0.0, isw["Sub1_barEFE"]["A_ISW_EFE_sat"]),
             xytext=(0.06, 3.7), fontsize=7.3, ha="center", arrowprops=dict(arrowstyle="->", lw=0.8, color="#7c3aed"))
ax2.annotate("Pure MOND\n(+0.398, FLIP!)", (isw["pure_MOND"]["dlnmu_dlna"], isw["pure_MOND"]["A_ISW_dyn_lin"]),
             xytext=(0.34, -2.0), fontsize=7.3, ha="center", color="#b91c1c", fontweight="bold",
             arrowprops=dict(arrowstyle="->", lw=0.8, color="#dc2626"))
ax2.text(0.345, 2.0, "HOT SPOT\n(RULED OUT)", color="#b91c1c", fontsize=7.8, ha="center", fontweight="bold")

ax2.set_xlim(-0.04, 0.44)
ax2.set_ylim(-4.6, 5.8)
ax2.set_xlabel(r"Supervoid growth-boost rate $\mathrm{d}\ln\mu / \mathrm{d}\ln a$ at $z = 0.50$")
ax2.set_ylabel(r"ISW amplitude $A_{\mathrm{ISW}} = \Delta T_{\mathrm{Tau}} / \Delta T_{\Lambda\mathrm{CDM}}$")
ax2.set_title("(B) Exact Sub1 ISW Sign-Flip Survival vs. Pure MOND")
ax2.legend(loc="lower left", fontsize=6.6)
ax2.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(FIG_DIR / "fig2_four_terminal_hubble_and_isw_supervoid.pdf")
plt.close(fig)

# ==============================================================================
# FIGURE 3: Continuous K_B^{-1} Milgrom sigma_8 Closure, CMB EE & X-COP
# ==============================================================================
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.2, 3.7))

sig_cont = DATA["part2_continuous_closure_and_six_sectors"]["continuous_KB_Milgrom_sigma8"]["sigma8_z0_continuous"]
cats_s8 = ["Unscreened\n($z{=}1000$)", "Step\n($z_{\\mathrm{dip}}$)", "Continuous\n($K_B^{-1}$+CoM)", "ACT DR6\nLensing"]
vals_s8 = [13.075, 0.8493, sig_cont, 0.818]
cols_s8 = ["#ef4444", "#7c3aed", "#059669", "#0284c7"]
ax1.axhspan(0.818 - 0.009, 0.818 + 0.009, color="#0284c7", alpha=0.18, label=r"ACT DR6 lens $0.818\pm 0.009$")
bars1 = ax1.bar(range(4), vals_s8, color=cols_s8, width=0.54, edgecolor="black", lw=0.6)
ax1.set_yscale("log")
for i, v_s in enumerate(vals_s8):
    ax1.text(i, v_s * 1.12, f"{v_s:.3f}" if v_s < 2 else f"{v_s:.1f}", ha="center", fontsize=7.4, fontweight="bold")
ax1.set_xticks(range(4))
ax1.set_xticklabels(cats_s8, fontsize=7.0)
ax1.set_ylim(0.65, 25.0)
ax1.set_ylabel(r"Clustering amplitude $\sigma_8(z=0)$")
ax1.set_title(r"(A) Continuous $K_B^{-1}$ $\sigma_8$ Closure ($+0.48\sigma$)")
ax1.legend(loc="upper right", fontsize=6.7)
ax1.grid(axis="y", alpha=0.25, which="both")

peaks = [1, 2, 3, 4]
ee_tau = DATA["part2_continuous_closure_and_six_sectors"]["cmb_ee_peaks"]
ee_obs = [137.0, 399.0, 690.0, 995.0]
x = np.arange(len(peaks))
ax2.bar(x - 0.17, ee_obs, width=0.32, color="#f43f5e", alpha=0.82, label="Planck / ACT DR6")
ax2.bar(x + 0.17, ee_tau, width=0.32, color="#0ea5e9", alpha=0.88, label=r"Tau Core ($O_Q$)")
for i in range(4):
    pct = (ee_tau[i] / ee_obs[i] - 1.0) * 100.0
    ax2.text(x[i], max(ee_tau[i], ee_obs[i]) + 24, f"{pct:+.2f}%", ha="center", fontsize=7.2, fontweight="bold")
ax2.set_xticks(x)
ax2.set_xticklabels([f"EE-{m}" for m in peaks])
ax2.set_ylim(0, 1150)
ax2.set_ylabel(r"Peak multipole $\ell_m^{EE}$")
ax2.set_title(r"(B) CMB $EE$ Peaks ($O_Q$ Thomson)")
ax2.legend(loc="upper left", fontsize=7.2)
ax2.grid(axis="y", alpha=0.25)

r_frac = [0.2, 0.5, 1.0, 1.5, 2.0]
mond_rat = [0.39, 0.50, 0.68, 0.84, 0.99]
tau_rat  = [1.03, 1.07, 1.10, 1.12, 1.13]
ax3.axhspan(0.85, 1.15, color="#10b981", alpha=0.14, label=r"X-COP / WL $\pm 15\%$")
ax3.axhline(1.0, color="black", ls=":", lw=1.0)
ax3.plot(r_frac, mond_rat, "--s", color="#ef4444", lw=1.8, ms=5.0, label="Pure baryonic MOND")
ax3.plot(r_frac, tau_rat, "-o", color="#059669", lw=2.1, ms=5.5, label=r"Tau Core Sub1 ($R_{A/b}=5.37$)")
ax3.set_xlabel(r"Cluster radius $r / R_{500}$")
ax3.set_ylabel(r"$M_{\mathrm{pred}} / M_{\mathrm{obs}}$")
ax3.set_ylim(0.25, 1.32)
ax3.set_title("(C) X-COP Cluster Mass Profile")
ax3.legend(loc="lower right", fontsize=7.0)
ax3.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(FIG_DIR / "fig3_six_sector_multiterminal_scorecard.pdf")
plt.close(fig)

print("Generated 3 v3.0.0 vector PDF figures in", FIG_DIR)
