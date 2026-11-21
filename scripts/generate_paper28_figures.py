#!/usr/bin/env python3
"""
Generate 3 publication-grade vector PDF figures for Tau Core Paper 28 (v1.0.0).
"""
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9.5,
    "axes.linewidth": 0.9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": True,
    "legend.framealpha": 0.92,
})

# ==============================================================================
# FIGURE 1: ACT DR6 4-Terminal Inflation (n_s, A_s) & BBN Delta N_eff(w_rh)
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.0))

# Panel A: (n_s, A_s) trajectory from Presc I -> Presc II -> 4-terminal ACT fixed point
# Observational bands
ax1.axvspan(0.974 - 0.003, 0.974 + 0.003, color="#3b82f6", alpha=0.16, label=r"ACT DR6+DESI $n_s = 0.974 \pm 0.003$")
ax1.axvspan(0.9649 - 0.0042, 0.9649 + 0.0042, color="#94a3b8", alpha=0.18, label=r"Planck 2018 $n_s = 0.9649 \pm 0.0042$")
ax1.axhspan(2.07, 2.13, color="#10b981", alpha=0.20, label=r"Planck/ACT $A_s = (2.10\pm 0.03)\times 10^{-9}$")

# Presc II curve at w_rh = 0 (only gravity/scalar terminal -> A_s drops!)
ns_w0 = [0.9620, 0.9649, 0.9693, 0.9741, 0.9782]
as_w0 = [2.101, 1.898, 1.643, 1.436, 1.349]
ax1.plot(ns_w0, as_w0, "--o", color="#ef4444", lw=1.7, ms=4.5, label=r"Presc. II alone ($w_{\mathrm{rh}}=0$, $O_G$ only)")

# 4-Terminal trajectory at fixed kappa_J = 0.0138 as w_rh goes 0 -> 1/3 -> 0.9015 -> 1.0
ns_4t = [0.9688, 0.9718, 0.9740, 0.9746]
as_4t = [1.66, 1.88, 2.100, 2.18]
ax1.plot(ns_4t, as_4t, "-s", color="#0284c7", lw=2.0, ms=5.0, label=r"4-Terminal ($O_G+O_T+O_Q+O_M$, $w_{\mathrm{rh}}\uparrow$)")

ax1.scatter([0.9620], [2.101], s=75, c="#334155", zorder=6, marker="D")
ax1.annotate("Presc. I\n($w_{\\mathrm{rh}}=0$)", (0.9620, 2.101), xytext=(0.9585, 2.24), fontsize=8, ha="center")
ax1.scatter([0.9740], [2.100], s=110, c="#059669", zorder=7, marker="*", edgecolors="black", linewidths=0.7)
ax1.annotate("4-Term. Fixed Point\n($w_{\\mathrm{rh}}=0.902, \\kappa_J=0.0138$)", (0.9740, 2.100),
             xytext=(0.9708, 2.34), fontsize=8, ha="center", arrowprops=dict(arrowstyle="->", lw=0.9))

ax1.set_xlim(0.956, 0.981)
ax1.set_ylim(1.22, 2.52)
ax1.set_xlabel(r"Scalar spectral index $n_s$ (NLO)")
ax1.set_ylabel(r"Primordial amplitude $10^9 A_s$")
ax1.set_title("(A) Multi-Terminal ACT DR6 Inflation Closure")
ax1.legend(loc="lower left", fontsize=7.2)
ax1.grid(alpha=0.25)

# Panel B: BBN & ACT DR6 Delta N_eff(w_rh)
w_vals = np.linspace(1.0 / 3.0, 1.0, 80)
h_ratio = 3.109e11
om_flat = 6.097e-18
dneff = []
for w in w_vals:
    ngw = 2.0 * (3.0 * w - 1.0) / (1.0 + 3.0 * w)
    fr = h_ratio ** ((1.0 + 3.0 * w) / (3.0 * (1.0 + w)))
    om_tot = om_flat * ((fr**ngw - 1.0) / ngw) if ngw > 1e-6 else om_flat * math.log(fr)
    dneff.append((8.0 / 7.0) * ((11.0 / 4.0) ** (4.0 / 3.0)) * (om_tot / 2.469e-5))

ax2.semilogy(w_vals, dneff, color="#7c3aed", lw=2.2, label=r"Tau Core primordial GW $\Delta N_{\mathrm{eff}}^{\mathrm{GW}}(w_{\mathrm{rh}})$")
ax2.axhline(0.17, color="#dc2626", ls="--", lw=1.6, label=r"BBN + ACT DR6 95% CL ($\Delta N_{\mathrm{eff}} < 0.17$)")
ax2.scatter([0.9015], [8.719e-6], s=100, c="#059669", marker="*", zorder=6, edgecolors="black", linewidths=0.7)
ax2.annotate(r"ACT Fixed Point" + "\n" + r"$\Delta N_{\mathrm{eff}} = 8.7\times 10^{-6}$",
             (0.9015, 8.719e-6), xytext=(0.64, 1.8e-4), fontsize=8,
             arrowprops=dict(arrowstyle="->", lw=0.9))
ax2.fill_between(w_vals, 0.17, 1.0, color="#dc2626", alpha=0.10)
ax2.set_xlim(0.32, 1.02)
ax2.set_ylim(1e-11, 1.0)
ax2.set_xlabel(r"Post-inflationary $O_T/O_Q$ equation of state $w_{\mathrm{rh}}$")
ax2.set_ylabel(r"Effective neutrino excess $\Delta N_{\mathrm{eff}}^{\mathrm{GW}}$")
ax2.set_title(r"(B) BBN & ACT DR6 $\Delta N_{\mathrm{eff}}$ Safety Margin")
ax2.legend(loc="lower right", fontsize=7.5)
ax2.grid(alpha=0.25, which="both")

fig.tight_layout()
fig.savefig(FIG_DIR / "fig1_act_dr6_multiterminal_inflation_and_bbn.pdf")
plt.close(fig)

# ==============================================================================
# FIGURE 2: 4-Terminal Hubble Budget & ISW Supervoid Sign-Flip Theorem
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.0))

# Panel A: 4-Terminal H0 waterfall
labels = ["ACT+Planck\n(CMB)", r"$+O_T$ glob." + "\n(H2-vol)", r"$+O_T$ loc." + "\n(Void-Clk)", r"$+O_G$ (50%)" + "\n(OCSR)", r"$+O_M/O_Q$" + "\n($\\Delta M_B$)", "Tau Core\n4-Term."]
bottoms = [0, 68.22, 68.64, 69.02, 72.17, 0]
heights = [68.22, 0.42, 0.38, 3.15, 1.10, 73.27]
colors = ["#475569", "#8b5cf6", "#06b6d4", "#f59e0b", "#ec4899", "#10b981"]

ax1.axhspan(73.04 - 1.04, 73.04 + 1.04, color="#ef4444", alpha=0.14, label=r"SH0ES $73.04 \pm 1.04$")
ax1.axhline(73.04, color="#ef4444", ls="--", lw=1.2)
ax1.axhspan(69.85 - 1.0, 69.85 + 1.0, color="#3b82f6", alpha=0.10, label=r"JWST TRGB/JAGB $\sim 70$")
bars = ax1.bar(range(len(labels)), heights, bottom=bottoms, color=colors, width=0.58, edgecolor="black", linewidth=0.6)
for i, (b, h) in enumerate(zip(bottoms, heights)):
    val_txt = f"{h:.2f}" if i in (0, 5) else f"+{h:.2f}"
    ax1.text(i, b + h + 0.16, val_txt, ha="center", va="bottom", fontsize=7.8, fontweight="bold")
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels, fontsize=7.6)
ax1.set_ylim(66.5, 75.5)
ax1.set_ylabel(r"Hubble constant $H_0\;[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$")
ax1.set_title("(A) Four-Terminal Hubble Budget ($\\delta_v=-15\\%, 150\\,\\mathrm{Mpc}$)")
ax1.legend(loc="upper left", fontsize=7.5)
ax1.grid(axis="y", alpha=0.25)

# Panel B: ISW Supervoid A_ISW vs d(ln mu)/d(ln a) Sign-Flip Theorem
dlnmu = np.linspace(-0.08, 0.36, 120)
f_gr_05 = 0.754
a_isw_50 = 2.50 * 1.15 * (1.0 - dlnmu / (1.0 - f_gr_05))
a_isw_100 = 4.00 * 1.15 * (1.0 - dlnmu / (1.0 - f_gr_05))

ax2.axhspan(3.0, 5.2, color="#06b6d4", alpha=0.16, label=r"DES Y3 / Kovács+ ($A_{\mathrm{ISW}} = 3.0\text{--}5.2$)")
ax2.axhline(1.0, color="#64748b", ls=":", lw=1.2, label=r"$\Lambda$CDM ($A_{\mathrm{ISW}} = 1$)")
ax2.axhline(0.0, color="black", lw=0.9)
ax2.axvline(1.0 - f_gr_05, color="#dc2626", ls="--", lw=1.5, label=r"Sign-flip No-Go: $1 - f(0.5) = 0.246$")
ax2.axvspan(1.0 - f_gr_05, 0.36, color="#dc2626", alpha=0.09)

ax2.plot(dlnmu, a_isw_50, color="#0284c7", lw=2.1, label=r"Tau Core 50%-screened ($\mu=2.50$)")
ax2.plot(dlnmu, a_isw_100, color="#7c3aed", lw=1.8, ls="-.", label=r"Tau Core full OCSR ($\mu=4.00$)")
ax2.scatter([0.0], [2.87], s=85, c="#0284c7", edgecolors="black", zorder=6)
ax2.scatter([0.0], [4.60], s=85, c="#7c3aed", edgecolors="black", zorder=6)
ax2.text(0.295, 1.8, "HOT SPOT\n(RULED OUT)", color="#b91c1c", fontsize=8, ha="center", fontweight="bold")

ax2.set_xlim(-0.08, 0.36)
ax2.set_ylim(-2.2, 6.2)
ax2.set_xlabel(r"Supervoid growth-boost rate $\mathrm{d}\ln\mu / \mathrm{d}\ln a$ at $z = 0.50$")
ax2.set_ylabel(r"ISW amplitude $A_{\mathrm{ISW}} = \Delta T_{\mathrm{Tau}} / \Delta T_{\Lambda\mathrm{CDM}}$")
ax2.set_title("(B) Supervoid ISW Anomaly & Sign-Flip Theorem")
ax2.legend(loc="lower left", fontsize=7.0)
ax2.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(FIG_DIR / "fig2_four_terminal_hubble_and_isw_supervoid.pdf")
plt.close(fig)

# ==============================================================================
# FIGURE 3: 6-Sector Multi-Terminal Scorecard (EE Peaks & Cluster Profiles)
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 3.8))

# Panel A: CMB EE polarization peaks (O_Q)
peaks = [1, 2, 3, 4]
ee_tau = [137.4, 398.6, 688.4, 990.4]
ee_obs = [137.0, 399.0, 690.0, 995.0]
x = np.arange(len(peaks))
ax1.bar(x - 0.17, ee_obs, width=0.32, color="#f43f5e", alpha=0.82, label="Planck PR3 / ACT DR6")
ax1.bar(x + 0.17, ee_tau, width=0.32, color="#0ea5e9", alpha=0.88, label=r"Tau Core ($u=\sqrt{3}/24, O_Q$)")
for i in range(4):
    pct = (ee_tau[i] / ee_obs[i] - 1.0) * 100.0
    ax1.text(x[i], max(ee_tau[i], ee_obs[i]) + 22, f"{pct:+.2f}%", ha="center", fontsize=7.8, fontweight="bold")
ax1.set_xticks(x)
ax1.set_xticklabels([f"EE-{m}" for m in peaks])
ax1.set_ylim(0, 1140)
ax1.set_ylabel(r"Peak multipole $\ell_m^{EE}$")
ax1.set_title(r"(A) CMB $EE$ Polarization Peaks ($O_Q$ Thomson Readout)")
ax1.legend(loc="upper left", fontsize=8)
ax1.grid(axis="y", alpha=0.25)

# Panel B: X-COP Cluster radial mass ratio (Pure MOND vs Tau Core Sub1)
r_frac = [0.2, 0.5, 1.0, 1.5, 2.0]
mond_rat = [0.39, 0.50, 0.68, 0.84, 0.99]
tau_rat  = [1.03, 1.07, 1.10, 1.12, 1.13]
ax2.axhspan(0.85, 1.15, color="#10b981", alpha=0.14, label=r"X-COP / WL $\pm 15\%$ envelope")
ax2.axhline(1.0, color="black", ls=":", lw=1.0)
ax2.plot(r_frac, mond_rat, "--s", color="#ef4444", lw=1.9, ms=5.5, label="Pure baryonic MOND/RAR (no dust)")
ax2.plot(r_frac, tau_rat, "-o", color="#059669", lw=2.2, ms=6.0, label=r"Tau Core Sub1 ($R_{A/b}=5.3681$ + QUMOND)")
ax2.set_xlabel(r"Cluster radius $r / R_{500}$ ($M_{500}=6\times 10^{14}\,M_\odot$)")
ax2.set_ylabel(r"Predicted / Observed Total Mass $M_{\mathrm{pred}} / M_{\mathrm{obs}}$")
ax2.set_ylim(0.25, 1.32)
ax2.set_title("(B) X-COP Cluster Multi-Terminal Mass Profile")
ax2.legend(loc="lower right", fontsize=7.6)
ax2.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(FIG_DIR / "fig3_six_sector_multiterminal_scorecard.pdf")
plt.close(fig)

print("Generated 3 vector PDF figures in", FIG_DIR)
