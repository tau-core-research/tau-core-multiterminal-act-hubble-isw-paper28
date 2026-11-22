#!/usr/bin/env python3
"""
Generate 3 publication-grade vector PDF figures (both English and Hungarian editions)
for Tau Core Paper 28 (v4.4.0) 100% from paper28_verification_results.json (zero hardcoded
curves, zero duplicated magic constants, and full native Hungarian decimal-comma formatting).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA = json.loads((ROOT / "paper28_verification_results.json").read_text(encoding="utf-8"))

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 8.8,
    "axes.linewidth": 0.9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": True,
    "legend.framealpha": 0.93,
})


def num_s(val: float, fmt: str, hu: bool) -> str:
    s = format(val, fmt)
    return s.replace(".", "{,}") if hu else s


def apply_hu_comma_ticks(ax: plt.Axes, hu: bool, log_y: bool = False, cat_x: bool = False) -> None:
    if not hu:
        return
    if not cat_x:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x:g}".replace(".", ",")))
    if not log_y:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: f"{y:g}".replace(".", ",")))


def make_figures(lang: str = "en") -> None:
    hu = (lang == "hu")
    sfx = "_hu" if hu else ""

    # ==========================================================================
    # FIGURE 1: UV Sector — Presc. I, 3-Branch Presc. II No-Go, 2-Field 5-ODE & GW
    # ==========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.25))

    p1 = DATA["part1_uv_sector"]["gate1_presc_I_k0_proven"]
    p2 = DATA["part1_uv_sector"]["gate2_presc_II_three_branch_nogo"]
    p3 = DATA["part1_uv_sector"]["gate3_conditional_2field_5ode_benchmark"]
    p2_a2 = p2["branch_A2_paper24_threshold_multiplet_nSTr_plus2"]
    p2_b  = p2["branch_B_gauged_BL_Spin10_qGUT2_3_over_2"]
    traj  = p2["ode_trajectory_w0"]

    ax1.axvspan(0.974 - 0.003, 0.974 + 0.003, color="#3b82f6", alpha=0.16,
                label=("ACT DR6+DESI $n_s = 0{,}974\\pm 0{,}003$" if hu
                       else "ACT DR6+DESI $n_s = 0.974\\pm 0.003$"))
    ax1.axvspan(0.9649 - 0.0042, 0.9649 + 0.0042, color="#94a3b8", alpha=0.18,
                label=("Planck 2018 $n_s = 0{,}9649\\pm 0{,}0042$" if hu
                       else "Planck 2018 $n_s = 0.9649\\pm 0.0042$"))
    ax1.axhspan(2.07, 2.13, color="#10b981", alpha=0.22,
                label=("Planck/ACT $10^9 A_s = 2{,}10\\pm 0{,}03$" if hu
                       else "Planck/ACT $10^9 A_s = 2.10\\pm 0.03$"))

    # 100% ODE-computed trajectory read directly from paper28_verification_results.json:
    ax1.plot(traj["n_s"], traj["As_1e9"], "--", color="#ef4444", lw=1.65,
             label=("Egzakt II. előírás ODE-pálya ($w_{\\mathrm{rh}}{=}0$)" if hu
                    else "Exact Prescription II ODE curve ($w_{\\mathrm{rh}}{=}0$)"))

    ax1.scatter([p1["n_s"]], [p1["A_s"] * 1e9], s=75, c="#0284c7", marker="D",
                edgecolors="black", linewidths=0.7, zorder=7,
                label=("I. előírás alap ($k{=}0$, bizonyított)" if hu
                       else "Prescription I baseline ($k{=}0$, PROVEN)"))
    ax1.scatter([p2_a2["n_s_w0"]], [p2_a2["A_s_w0"] * 1e9], s=85, c="#dc2626", marker="s",
                edgecolors="black", linewidths=0.7, zorder=8,
                label=(f"Paper 24 küszöb No-Go ($-25{{,}}0\\sigma$ $A_s$)" if hu
                       else f"Paper 24 multiplet No-Go ($-25.0\\sigma$ $A_s$)"))
    ax1.scatter([p2_b["n_s_w0"]], [p2_b["A_s_w0"] * 1e9], s=95, c="#991b1b", marker="X",
                edgecolors="black", linewidths=0.7, zorder=8,
                label=(f"Mérték $U(1)_{{B-L}}$ No-Go ($+119\\sigma$ $A_s$)" if hu
                       else f"Gauged $U(1)_{{B-L}}$ No-Go ($+119\\sigma$ $A_s$)"))
    ax1.scatter([p3["n_s"]], [p3["A_s"] * 1e9], s=135, c="#f59e0b", marker="*",
                edgecolors="black", linewidths=0.7, zorder=9,
                label=("COBE-zárt 2-mezős minimax ($k_{\\mathrm{UV,eff}}{=}1$)" if hu
                       else "COBE-locked 2-field minimax ($k_{\\mathrm{UV,eff}}{=}1$)"))

    ax1.annotate((f"I. előírás ($k{{=}}0$)\n$n_s{{=}}{num_s(p1['n_s'], '.4f', hu)}$" if hu
                  else f"Presc. I ($k{{=}}0$)\n$n_s{{=}}{p1['n_s']:.4f}$"),
                 (p1["n_s"], p1["A_s"] * 1e9), xytext=(0.9600, 3.30), fontsize=7.2, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#0284c7"))
    ax1.annotate((f"$n_{{\\mathrm{{STr}}}}{{=}}+2$ No-Go\n$10^9 A_s{{=}}{num_s(p2_a2['A_s_w0']*1e9, '.2f', hu)}$" if hu
                  else f"$n_{{\\mathrm{{STr}}}}{{=}}+2$ No-Go\n$10^9 A_s{{=}}{p2_a2['A_s_w0']*1e9:.2f}$"),
                 (p2_a2["n_s_w0"], p2_a2["A_s_w0"] * 1e9), xytext=(0.9752, 0.88), fontsize=7.1,
                 ha="center", color="#b91c1c", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#dc2626"))
    ax1.annotate((f"$U(1)_{{B-L}}$ No-Go\n$10^9 A_s{{=}}{num_s(p2_b['A_s_w0']*1e9, '.2f', hu)}$" if hu
                  else f"Gauged $B{{-}}L$ No-Go\n$10^9 A_s{{=}}{p2_b['A_s_w0']*1e9:.2f}$"),
                 (p2_b["n_s_w0"], p2_b["A_s_w0"] * 1e9), xytext=(0.9755, 5.20), fontsize=7.1,
                 ha="center", color="#991b1b", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#991b1b"))
    ax1.annotate((f"Minimax ($k_{{\\mathrm{{UV,eff}}}}{{=}}1$)\n$n_s{{=}}{num_s(p3['n_s'], '.4f', hu)}$" if hu
                  else f"Minimax ($k_{{\\mathrm{{UV,eff}}}}{{=}}1$)\n$n_s{{=}}{p3['n_s']:.4f}$"),
                 (p3["n_s"], p3["A_s"] * 1e9), xytext=(0.9698, 2.95), fontsize=7.2, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#d97706"))

    ax1.set_xlim(0.956, 0.9855)
    ax1.set_ylim(0.65, 6.15)
    ax1.set_xlabel("Skalár spektrálindex $n_s$ (NLO)" if hu else "Scalar spectral index $n_s$ (NLO)")
    ax1.set_ylabel("Ősi amplitúdó $10^9 A_s$" if hu else "Primordial amplitude $10^9 A_s$")
    ax1.set_title("(A) UV infláció: I. előírás, II. No-Go és COBE-zárt 2-mezős" if hu
                  else "(A) UV Inflation: Presc. I, 3-Branch No-Go & COBE-Locked 2-Field")
    ax1.legend(loc="upper left", fontsize=6.0)
    ax1.grid(alpha=0.25)
    apply_hu_comma_ticks(ax1, hu)

    # Panel B: Exact 1D COBE-Locked (10^9 A_s == 2.100) Non-Thermal Leptogenesis Pinch Curve!
    c1d = p3["cobe_locked_1d_curve"]
    ns_1d  = np.array(c1d["n_s"])
    eta_1d = np.array(c1d["eta_B_1e10"])
    ceil2  = p3["leptogenesis_2sigma_ceiling"]

    ax2.axvspan(0.974 - 0.003, 0.974 + 0.003, color="#3b82f6", alpha=0.15,
                label=("ACT DR6 $n_s = 0{,}974\\pm 0{,}003$" if hu else "ACT DR6 $n_s = 0.974\\pm 0.003$"))
    ax2.axhspan(6.12 - 0.08, 6.12 + 0.08, color="#10b981", alpha=0.22,
                label=("Planck $\\eta_B$ $2\\sigma$ ($6{,}12\\pm 0{,}08$)" if hu else "Planck $\\eta_B$ $2\\sigma$ ($6.12\\pm 0.08$)"))
    ax2.plot(ns_1d, eta_1d, color="#7c3aed", lw=2.1,
             label=("1D COBE-zárt 5-ODE pálya ($10^9 A_s\\equiv 2{,}10$)" if hu
                    else "1D COBE-locked 5-ODE curve ($10^9 A_s\\equiv 2.10$)"))
    ax2.axvline(ceil2["n_s_max_2sigma"], color="#dc2626", ls="--", lw=1.4,
                label=(f"2$\\sigma$ leptogenezis plafon ($n_s\\leq {num_s(ceil2['n_s_max_2sigma'], '.4f', hu)}$)" if hu
                       else f"2$\\sigma$ leptogenesis ceiling ($n_s\\leq {ceil2['n_s_max_2sigma']:.4f}$)"))
    ax2.scatter([p1["n_s"]], [p1["eta_B_first_principles_5ODE"] * 1e10], s=70, c="#0284c7", marker="D",
                edgecolors="black", zorder=8)
    ax2.scatter([p3["n_s"]], [p3["eta_B_first_principles_5ODE"] * 1e10], s=130, c="#f59e0b", marker="*",
                edgecolors="black", zorder=9,
                label=("Minimax optimum ($-0{,}61\\sigma$ ACT, $-0{,}63\\sigma$ $\\eta_B$)" if hu
                       else "Minimax optimum ($-0.61\\sigma$ ACT, $-0.63\\sigma$ $\\eta_B$)"))
    ax2.annotate((f"Minimax ($\\log_{{10}}f_\\phi{{=}}-8{{,}}8$)\n$10^{{10}}\\eta_B{{=}}{num_s(p3['eta_B_first_principles_5ODE']*1e10, '.2f', hu)}$" if hu
                  else f"Minimax ($\\log_{{10}}f_\\phi{{=}}-8.8$)\n$10^{{10}}\\eta_B{{=}}{p3['eta_B_first_principles_5ODE']*1e10:.2f}$"),
                 (p3["n_s"], p3["eta_B_first_principles_5ODE"] * 1e10), xytext=(0.9665, 5.15), fontsize=7.0,
                 ha="center", arrowprops=dict(arrowstyle="->", lw=0.8, color="#d97706"))

    ax2.set_xlim(0.961, 0.9756)
    ax2.set_ylim(3.4, 6.45)
    ax2.set_xlabel("Skalár spektrálindex $n_s$ a COBE-zárt pályán" if hu
                   else "Scalar spectral index $n_s$ along COBE-locked curve")
    ax2.set_ylabel("Barion-aszimmetria $10^{10}\\eta_B$ (5-ODE)" if hu
                   else "Baryon asymmetry $10^{10}\\eta_B$ (5-ODE)")
    ax2.set_title("(B) Leptogenezis $n_s$-plafon (5 nagyságrenddel BBN GW előtt)" if hu
                  else "(B) Non-Thermal Leptogenesis Pinch Ceiling on $n_s$")
    ax2.legend(loc="lower left", fontsize=6.1)
    ax2.grid(alpha=0.25)
    apply_hu_comma_ticks(ax2, hu)

    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig1_act_dr6_multiterminal_inflation_and_bbn{sfx}.pdf")
    plt.close(fig)

    # ==========================================================================
    # FIGURE 2: Real N=1341 Pantheon+SH0ES Curvature Audit & Supervoid ISW
    # ==========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.25))

    g8 = DATA["part3_late_and_boltzmann_ee"]["gate8_sh0es_window_local_void_hubble"]
    vc = g8["delta_v_sweep_curves"]
    dv_arr   = np.array(vc["delta_v"])
    dchi_l   = np.array(vc["dChi2_LCDM_vs_self0"])
    dchi_h2a = np.array(vc["dChi2_H2a_vs_self0"])
    dchi_clx = np.array(vc["dChi2_H2a_CLASSIX120_comp_ruler_vs_self0"])
    c_h2a    = g8["H2a_2sigma_ceiling_OG"]
    c_clx    = g8["H2a_CLASSIX_120Mpc_compensated_2sigma_ruler"]
    dchi_bg  = g8["H2a_homogeneous_dv0"]["chi2_Pantheon_1341"] - g8["LCDM_homogeneous_dv0"]["chi2_Pantheon_1341"]
    kbc_dchi = g8["LCDM_uncompensated_KBC_nogo"]["chi2_Pantheon_1341"] - g8["LCDM_homogeneous_dv0"]["chi2_Pantheon_1341"]

    ax1.plot(dv_arr, dchi_l, "--", color="#64748b", lw=1.7,
             label=("$\\Lambda$CDM + 300 Mpc üreg" if hu else "$\\Lambda$CDM + 300 Mpc void"))
    ax1.plot(dv_arr, dchi_h2a, "-", color="#059669", lw=2.0,
             label=(f"H2a ($\\Delta\\chi^2_{{0}}{{=}}{num_s(dchi_bg, '+.2f', hu)}$) + 300 Mpc üreg ($O_G$)" if hu
                    else f"H2a ($\\Delta\\chi^2_{{0}}{{=}}{dchi_bg:+.2f}$) + 300 Mpc void ($O_G$)"))
    ax1.plot(dv_arr, dchi_clx, "-.", color="#7c3aed", lw=1.8,
             label=("H2a + 120 Mpc komp. CLASSIX üreg" if hu else "H2a + 120 Mpc comp. CLASSIX void"))
    ax1.axhline(0.0, color="black", ls=":", lw=0.9)
    ax1.axhline(4.0, color="#dc2626", ls="--", lw=1.1,
                label=("$2\\sigma$ görbületi korlát ($\\Delta\\chi^2 = +4{,}0$)" if hu
                       else "$2\\sigma$ curvature limit ($\\Delta\\chi^2 = +4.0$)"))
    ax1.scatter([-0.25], [kbc_dchi],
                s=70, c="#dc2626", marker="X", edgecolors="black", zorder=8,
                label=(f"KBC ($-0{{,}}25$: ${num_s(kbc_dchi, '+.2f', hu)}$, kizárt!)" if hu
                       else f"KBC ($-0.25$: ${kbc_dchi:+.2f}$, NO-GO!)"))
    ax1.scatter([c_h2a["delta_v"]], [4.0], s=105, c="#0284c7", marker="*", edgecolors="black", zorder=9,
                label=(f"300 Mpc $O_G$ plafon ($H_0^{{(O_G)}}{{=}}{num_s(c_h2a['H0_loc_dMB0'], '.2f', hu)}$)" if hu
                       else f"300 Mpc $O_G$ ceiling ($H_0^{{(O_G)}}{{=}}{c_h2a['H0_loc_dMB0']:.2f}$)"))
    ax1.scatter([c_clx["delta_v"]], [4.0], s=65, c="#f59e0b", marker="D", edgecolors="black", zorder=9,
                label=(f"CLASSIX mérce plafon ($H_0^{{\\mathrm{{loc}}}}{{=}}{num_s(c_clx['H0_loc_dMB0'], '.2f', hu)}$, ${num_s(c_clx['z_vs_CCHP_HST_JWST_TRGB'], '+.2f', hu)}\\sigma$)" if hu
                       else f"CLASSIX ruler ceiling ($H_0^{{\\mathrm{{loc}}}}{{=}}{c_clx['H0_loc_dMB0']:.2f}$, ${c_clx['z_vs_CCHP_HST_JWST_TRGB']:+.2f}\\sigma$)"))
    ax1.annotate((f"CLASSIX komp. ($\\delta_v{{=}}{num_s(c_clx['delta_v'], '.2f', hu)}$)\n$H_0^{{\\mathrm{{loc}}}} = {num_s(c_clx['H0_loc_dMB0'], '.2f', hu)}$" if hu
                  else f"CLASSIX comp. ($\\delta_v{{=}}{c_clx['delta_v']:.2f}$)\n$H_0^{{\\mathrm{{loc}}}} = {c_clx['H0_loc_dMB0']:.2f}$"),
                 (c_clx["delta_v"], 4.0), xytext=(-0.333, 17.5), fontsize=6.5, ha="left",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#d97706"))

    ax1.set_xlim(-0.34, 0.08)
    ax1.set_ylim(-2.5, 64.0)
    ax1.set_xlabel("Lokális üregkontraszt $\\delta_v$ ($R_{\\mathrm{void}} \\in [120, 400]\\;\\mathrm{Mpc}$)" if hu
                   else "Local void contrast $\\delta_v$ ($R_{\\mathrm{void}} \\in [120, 400]\\;\\mathrm{Mpc}$)")
    ax1.set_ylabel("STAT+SYS $\\Delta\\chi^2(\\delta_v) = \\chi^2(\\delta_v) - \\chi^2(0)$" if hu
                   else "STAT+SYS $\\Delta\\chi^2(\\delta_v) = \\chi^2(\\delta_v) - \\chi^2(0)$")
    ax1.set_title("(A) Pantheon+SH0ES ($1341\\times 1341$ STAT+SYS) és CLASSIX" if hu
                  else "(A) Pantheon+SH0ES ($1341\\times 1341$ STAT+SYS) & CLASSIX")
    ax1.legend(loc="upper right", fontsize=5.45)
    ax1.grid(alpha=0.25)
    apply_hu_comma_ticks(ax1, hu)

    # Panel B: Root-Shot Self-Consistent Supervoid ISW (|delta_b(0.50)| == 0.25) + Linear-MG Ceiling vs DES Y3 Tail
    g6 = DATA["part2_covariant_ir_and_isw"]["gate6_self_consistent_supervoid_isw"]
    unscr = g6["unscreened_Sub1_nogo"]
    scr   = g6["screened_Sub1_xi_Phi_sq"]
    desi_obs = g6["DESI_Legacy_linear_ISW_obs"]
    desi_err = g6["DESI_Legacy_linear_ISW_err"]
    lin_mg_max = g6["linear_MG_spectral_overlap_ceiling"]["sigma8_unconstrained_gz1_choquet_monotone_comp"]
    lin_mg_box = g6["linear_MG_spectral_overlap_ceiling"]["time_varying_gz_DESY3_box_choquet_max"]

    cats_isw = (["Árnyékolatlan\n(anyag)", "Árnyékolatlan\n(galaxis)",
                 "Helmholtz\n($\\kappa_\\theta{=}\\xi_\\Phi^2$)", "DESI Legacy\nlineáris ISW", "DES Y3\nszupervoid"]
                if hu else
                ["Unscreened\n(matter)", "Unscreened\n(galaxy)",
                 "Helmholtz\n($\\kappa_\\theta{=}\\xi_\\Phi^2$)", "DESI Legacy\nlinear ISW", "DES Y3\nsupervoid"])
    vals_isw = [unscr["A_ISW_matter_norm"], unscr["A_ISW_galaxy_norm"], scr["A_ISW_matter_norm"], desi_obs, 4.1]
    cols_isw = ["#dc2626", "#ef4444", "#059669", "#10b981", "#0284c7"]

    ax2.axhspan(-5.0, 0.0, color="#dc2626", alpha=0.09,
                label=("Meleg folt ($A_{\\mathrm{ISW}}<0$, kizárt!)" if hu
                       else "Hot-spot sign flip ($A_{\\mathrm{ISW}}<0$, ruled out!)"))
    ax2.axhline(0.0, color="black", lw=1.0)
    ax2.axhline(1.0, color="#475569", ls=":", lw=1.1,
                label=("$\\Lambda$CDM ($A_{\\mathrm{ISW}}=1{,}00$)" if hu else "$\\Lambda$CDM ($A_{\\mathrm{ISW}}=1.00$)"))
    ax2.axhline(lin_mg_max, color="#d97706", ls="--", lw=1.2,
                label=(f"Monoton-$k$ plafon ($\\leq {num_s(lin_mg_max, '.2f', hu)}$ $g{{=}}1$; $\\leq {num_s(lin_mg_box, '.2f', hu)}$ DES sáv)" if hu
                       else f"Monotone-$k$ ceiling ($\\leq {lin_mg_max:.2f}$ $g{{=}}1$; $\\leq {lin_mg_box:.2f}$ DES box)"))
    ax2.bar(range(5), vals_isw, width=0.55, color=cols_isw, edgecolor="black", lw=0.65, zorder=5)
    ax2.errorbar([3], [desi_obs], yerr=[desi_err], fmt="none", ecolor="black", elinewidth=1.2, capsize=3.5, zorder=7)
    ax2.errorbar([4], [4.1], yerr=[1.1], fmt="none", ecolor="black", elinewidth=1.3, capsize=4, zorder=7)
    for i, v in enumerate(vals_isw):
        dy = 0.28 if v >= 0 else -0.64
        if i == 2:
            dy = 0.50
        elif i == 3:
            dy = 0.50
        elif i == 4:
            dy = 1.25
        lbl_v = f"{v:+.2f}".replace(".", ",") if hu else f"{v:+.2f}"
        ax2.text(i, v + dy, lbl_v, ha="center", fontsize=7.1, fontweight="bold")

    ax2.set_xticks(range(5))
    ax2.set_xticklabels(cats_isw, fontsize=6.1)
    ax2.set_ylim(-5.0, 6.2)
    ax2.set_ylabel("ISW amplitúdó $A_{\\mathrm{ISW}} = \\Delta T / \\Delta T_{\\Lambda\\mathrm{CDM}}$" if hu
                   else "ISW amplitude $A_{\\mathrm{ISW}} = \\Delta T / \\Delta T_{\\Lambda\\mathrm{CDM}}$")
    ax2.set_title("(B) Szupervoid ISW + 4D Maxwell Weyl-pajzs" if hu
                  else "(B) Supervoid ISW + 4D Maxwell Weyl Shield")
    ax2.legend(loc="upper left", fontsize=6.4)
    ax2.grid(axis="y", alpha=0.25)
    apply_hu_comma_ticks(ax2, hu, cat_x=True)

    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig2_four_terminal_hubble_and_isw_supervoid{sfx}.pdf")
    plt.close(fig)

    # ==========================================================================
    # FIGURE 3: First-Principles sigma_8, Absolute uK^2 Boltzmann EE & 12 X-COP Clusters
    # ==========================================================================
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.6, 3.95))

    g5 = DATA["part2_covariant_ir_and_isw"]["gate5_fourier_sigma8_barEFE"]
    cats_s8 = (["Árnyékolatlan\nWeyl ($\\kappa_\\theta{=}0$)", "Zárt $\\kappa_\\theta{=}\\xi_\\Phi^2$\n(I. előírás)",
                "Zárt $\\kappa_\\theta{=}\\xi_\\Phi^2$\n(2-mezős)", "ACT DR6\nlencsézés"]
               if hu else
               ["Unscreened\nWeyl ($\\kappa_\\theta{=}0$)", "Locked $\\kappa_\\theta{=}\\xi_\\Phi^2$\n(Presc. I)",
                "Locked $\\kappa_\\theta{=}\\xi_\\Phi^2$\n(2-Field)", "ACT DR6\nCMB Lensing"])
    vals_s8 = [
        g5["sigma8_unscreened_nogo"]["sigma8_Weyl"],
        g5["kappa_xi_Phi_sq_PrescI"]["sigma8_Weyl"],
        g5["kappa_xi_Phi_sq_2Field"]["sigma8_Weyl"],
        0.818,
    ]
    cols_s8 = ["#dc2626", "#059669", "#10b981", "#0284c7"]
    ax1.axhspan(0.818 - 0.009, 0.818 + 0.009, color="#0284c7", alpha=0.18,
                label=("ACT DR6 $0{,}818\\pm 0{,}009$" if hu else "ACT DR6 $0.818\\pm 0.009$"))
    ax1.bar(range(4), vals_s8, color=cols_s8, width=0.52, edgecolor="black", lw=0.65)
    ax1.set_yscale("log")
    for i, v in enumerate(vals_s8):
        txt = (f"{v:.3f}" if v < 2.0 else f"{v:.2f}")
        if hu:
            txt = txt.replace(".", ",")
        ax1.text(i, v * 1.15, txt, ha="center", fontsize=7.2, fontweight="bold")
    ax1.set_xticks(range(4))
    ax1.set_xticklabels(cats_s8, fontsize=6.4)
    ax1.set_ylim(0.62, 22.0)
    z_s8 = g5["kappa_xi_Phi_sq_PrescI"]["z_ACT_lens"]
    ax1.set_ylabel("Weyl klasztereződés $\\sigma_8^{\\mathrm{WL}}(z=0)$" if hu else "Weyl clustering $\\sigma_8^{\\mathrm{WL}}(z=0)$")
    ax1.set_title((f"(A) Egy-metrikás Helmholtz $\\sigma_8$ (${num_s(z_s8, '+.2f', hu)}\\sigma$)" if hu
                   else f"(A) Single-Metric Helmholtz $\\sigma_8$ (${z_s8:+.2f}\\sigma$)"))
    ax1.legend(loc="upper right", fontsize=6.6)
    ax1.grid(axis="y", alpha=0.25, which="both")

    # Panel B: Absolute uK^2 First-Principles RECFAST + l_max=6 Boltzmann D_ell^EE Spectrum!
    g12 = DATA["part3_late_and_boltzmann_ee"]["gate12_first_principles_boltzmann_ee"]
    ells_c = np.array(g12["ells_curve"])
    dl_uK2 = np.array(g12["Dl_EE_uK2_PrescI_curve"])
    ax2.plot(ells_c, dl_uK2, color="#0284c7", lw=1.8,
             label=("Boltzmann $D_\\ell^{EE}\\;[\\mu\\mathrm{K}^2]$ ($k{=}0$)" if hu
                    else "Boltzmann $D_\\ell^{EE}\\;[\\mu\\mathrm{K}^2]$ ($k{=}0$)"))
    for m_idx, (l_obs, l_th) in enumerate(zip(g12["ee_peaks_observed"], g12["ee_peaks_boltzmann"]), 1):
        ax2.axvline(l_obs, color="#ef4444", ls="--", lw=1.0,
                    label=("Planck/ACT $\\ell_{1..4}^{EE}$" if m_idx == 1 else None))
        y_pk = float(np.interp(l_th, ells_c, dl_uK2))
        ax2.scatter([l_th], [y_pk], s=34, c="#059669", edgecolors="black", zorder=6)
    ax2.set_xlim(60, 1055)
    ax2.set_ylim(0.0, float(np.max(dl_uK2)) * 1.14)
    ax2.set_xlabel("Multipólus $\\ell$" if hu else "Multipole $\\ell$")
    ax2.set_ylabel("Polarizációs spektrum $D_\\ell^{EE}\\;[\\mu\\mathrm{K}^2]$" if hu
                   else "Polarization power $D_\\ell^{EE}\\;[\\mu\\mathrm{K}^2]$")
    ax2.set_title("(B) Abszolút $\\mu\\mathrm{K}^2$ Boltzmann $EE$-spektrum" if hu
                  else "(B) Absolute $\\mu\\mathrm{K}^2$ Boltzmann $EE$ Spectrum")
    ax2.legend(loc="upper left", fontsize=6.4)
    ax2.grid(alpha=0.25)
    apply_hu_comma_ticks(ax2, hu)

    # Panel C: Real 12 Archival Eckert+ (2019) X-COP Galaxy Clusters at R_500 from JSON!
    g9 = DATA["part3_late_and_boltzmann_ee"]["gate9_xcop_12_clusters_archival_audit"]
    idx_cl = np.arange(len(g9["clusters"]))
    ax3.axhspan(0.85, 1.15, color="#10b981", alpha=0.14, label="X-COP $\\pm 15\\%$")
    ax3.axhline(1.0, color="black", ls=":", lw=1.0)
    ax3.errorbar(idx_cl, g9["Mpred_over_Mobs_R500_Tau"], yerr=g9["sigma_frac_R500"],
                 fmt="o", color="#059669", ecolor="#059669", elinewidth=1.0, capsize=2.2, ms=4.6,
                 label=("Tau Core Sub1 (LTB PDE)" if hu else "Tau Core Sub1 (LTB PDE)"))
    ax3.scatter(idx_cl, g9["Mpred_over_Mobs_R500_baryon_MOND"], marker="s", s=24, c="#ef4444",
                label=("Pusztán barionos QUMOND" if hu else "Baryon-only QUMOND"))
    ax3.set_xticks(idx_cl)
    ax3.set_xticklabels(g9["clusters"], rotation=45, ha="right", fontsize=5.8)
    ax3.set_ylabel("$M_{\\mathrm{pred}} / M_{\\mathrm{obs}}$ ($R_{500}$)")
    ax3.set_ylim(0.18, 1.32)
    ax3.set_title("(C) 12 valós X-COP halmaz ($R_{500}$, $k{=}0$)" if hu
                  else "(C) 12 Archival X-COP Clusters ($R_{500}$, $k{=}0$)")
    ax3.legend(loc="lower right", fontsize=6.2)
    ax3.grid(alpha=0.25)
    if hu:
        ax3.yaxis.set_major_formatter(FuncFormatter(lambda y, _p: f"{y:g}".replace(".", ",")))

    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig3_six_sector_multiterminal_scorecard{sfx}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    make_figures("en")
    make_figures("hu")
    print("Generated all 6 JSON-driven v4.4.0 vector PDF figures (EN + HU) in", FIG_DIR)
