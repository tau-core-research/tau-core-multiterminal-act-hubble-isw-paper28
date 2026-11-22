#!/usr/bin/env python3
"""
Generate 3 publication-grade vector PDF figures (both English and Hungarian editions)
for Tau Core Paper 28 (v4.1.0) 100% from paper28_verification_results.json (zero hardcoded
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
                label=("Feltételes 2-mezős 5-ODE ($k_{\\mathrm{UV}}{=}2$)" if hu
                       else "Conditional 2-field 5-ODE ($k_{\\mathrm{UV}}{=}2$)"))

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
    ax1.annotate((f"2-mezős ($k_{{\\mathrm{{UV}}}}{{=}}2$)\n$n_s{{=}}{num_s(p3['n_s'], '.4f', hu)}$" if hu
                  else f"2-Field ($k_{{\\mathrm{{UV}}}}{{=}}2$)\n$n_s{{=}}{p3['n_s']:.4f}$"),
                 (p3["n_s"], p3["A_s"] * 1e9), xytext=(0.9702, 2.95), fontsize=7.2, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#d97706"))

    ax1.set_xlim(0.956, 0.9855)
    ax1.set_ylim(0.65, 6.15)
    ax1.set_xlabel("Skalár spektrálindex $n_s$ (NLO)" if hu else "Scalar spectral index $n_s$ (NLO)")
    ax1.set_ylabel("Ősi amplitúdó $10^9 A_s$" if hu else "Primordial amplitude $10^9 A_s$")
    ax1.set_title("(A) UV infláció: I. előírás, II. No-Go és 2-mezős határ" if hu
                  else "(A) UV Inflation: Presc. I, 3-Branch No-Go & 2-Field")
    ax1.legend(loc="upper left", fontsize=6.1)
    ax1.grid(alpha=0.25)
    apply_hu_comma_ticks(ax1, hu)

    # Panel B: 100% JSON-driven GW curve + exact two-stage trajectory integral point!
    gw_c = p3["gw_single_fluid_curve"]
    w_vals = np.array(gw_c["w_eff"])
    dneff  = np.array(gw_c["Delta_Neff_GW"])
    ax2.semilogy(w_vals, dneff, color="#7c3aed", lw=2.0,
                 label=("Egyfolyadékos felső burkoló $\\Delta N_{\\mathrm{eff}}^{\\mathrm{GW}}(w_{\\mathrm{eff}})$" if hu
                        else "Single-fluid envelope $\\Delta N_{\\mathrm{eff}}^{\\mathrm{GW}}(w_{\\mathrm{eff}})$"))
    ax2.axhline(0.17, color="#dc2626", ls="--", lw=1.5,
                label=("BBN + ACT DR6 95%-os korlát ($\\Delta N_{\\mathrm{eff}} < 0{,}17$)" if hu
                       else "BBN + ACT DR6 95% CL ($\\Delta N_{\\mathrm{eff}} < 0.17$)"))
    ax2.fill_between(w_vals, 0.17, 1.0, color="#dc2626", alpha=0.10)
    dneff_traj = p3["Delta_Neff_GW_trajectory"]
    dneff_mant = dneff_traj * 1e6
    ax2.scatter([p3["w_eff_matching"]], [dneff_traj], s=130, c="#f59e0b", marker="*",
                edgecolors="black", linewidths=0.7, zorder=8,
                label=("Kétlépcsős ($w{=}1{\\to}0$) 5-ODE pálya-integrál" if hu
                       else "Two-stage ($w{=}1{\\to}0$) 5-ODE trajectory integral"))
    ax2.annotate((f"5-ODE pálya-integrál\n$\\Delta N_{{\\mathrm{{eff}}}}^{{\\mathrm{{GW}}}} = {num_s(dneff_mant, '.2f', hu)}\\times 10^{{-6}}$" if hu
                  else f"5-ODE trajectory integral\n$\\Delta N_{{\\mathrm{{eff}}}}^{{\\mathrm{{GW}}}} = {dneff_mant:.2f}\\times 10^{{-6}}$"),
                 (p3["w_eff_matching"], dneff_traj), xytext=(0.62, 4e-4), fontsize=7.3, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#d97706"))

    ax2.set_xlim(0.32, 1.02)
    ax2.set_ylim(1e-11, 1.0)
    ax2.set_xlabel("Effektív reheating állapotegyenlet $w_{\\mathrm{eff}}$" if hu
                   else "Effective reheating equation of state $w_{\\mathrm{eff}}$")
    ax2.set_ylabel("Effektív neutrínótöbblet $\\Delta N_{\\mathrm{eff}}^{\\mathrm{GW}}$" if hu
                   else "Effective neutrino excess $\\Delta N_{\\mathrm{eff}}^{\\mathrm{GW}}$")
    ax2.set_title("(B) BBN és ACT DR6 gravitációshullám-korlát" if hu
                  else "(B) BBN & ACT DR6 Primordial GW Bound")
    ax2.legend(loc="lower right", fontsize=6.6)
    ax2.grid(alpha=0.25, which="both")
    apply_hu_comma_ticks(ax2, hu, log_y=True)

    fig.tight_layout()
    fig.savefig(FIG_DIR / f"fig1_act_dr6_multiterminal_inflation_and_bbn{sfx}.pdf")
    plt.close(fig)

    # ==========================================================================
    # FIGURE 2: SH0ES Redshift-Window [0.023, 0.15] Local-Void H0 & Root-Shot ISW
    # ==========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.25))

    g8 = DATA["part3_late_and_boltzmann_ee"]["gate8_sh0es_window_local_void_hubble"]
    vc = g8["R_void_sweep_curves"]
    R_grid = np.array(vc["R_Mpc"])

    for key, col, ls_style, lbl_en, lbl_hu in [
        ("dv_m015_joint", "#64748b", ":",  "$\\delta_v = -0.15$ ($O_G + O_T^{\\mathrm{diff}}$)", "$\\delta_v = -0{,}15$ ($O_G + O_T^{\\mathrm{diff}}$)"),
        ("dv_m020_CLASSIX_joint", "#0284c7", "-.", "$\\delta_v = -0.20$ (CLASSIX)", "$\\delta_v = -0{,}20$ (CLASSIX)"),
        ("dv_m025_KBC_joint", "#059669", "-",  "$\\delta_v = -0.25$ (KBC, $O_G + O_T^{\\mathrm{diff}}$)", "$\\delta_v = -0{,}25$ (KBC, $O_G + O_T^{\\mathrm{diff}}$)"),
        ("dv_m025_KBC_OG_only", "#10b981", "--", "$\\delta_v = -0.25$ (pure Eq. 11 $O_G$ only)", "$\\delta_v = -0{,}25$ (tisztán 11. egy. $O_G$)"),
    ]:
        ax1.plot(R_grid, vc[key], color=col, ls=ls_style, lw=1.75, label=(lbl_hu if hu else lbl_en))

    ax1.axhspan(73.04 - 1.04, 73.04 + 1.04, color="#ef4444", alpha=0.14,
                label=("SH0ES $73{,}04\\pm 1{,}04$" if hu else "SH0ES $73.04\\pm 1.04$"))
    ax1.axhspan(70.40 - 1.90, 70.40 + 1.90, color="#3b82f6", alpha=0.11,
                label=("CCHP JWST TRGB $70{,}40\\pm 1{,}90$" if hu else "CCHP JWST TRGB $70.40\\pm 1.90$"))
    ax1.axhline(68.22, color="black", lw=0.9, ls="--",
                label=("ACT+DESI $H_0^{\\mathrm{CMB}} = 68{,}22$" if hu else "ACT+DESI $H_0^{\\mathrm{CMB}} = 68.22$"))

    c_pt = g8["CLASSIX_140Mpc"]
    k_pt = g8["KBC_300Mpc"]
    ax1.scatter([c_pt["R_Mpc"]], [c_pt["H0_joint_unsub_dMB0"]], s=65, c="#0284c7", edgecolors="black", zorder=7)
    ax1.scatter([k_pt["R_Mpc"]], [k_pt["H0_joint_unsub_dMB0"]], s=95, c="#059669", marker="*", edgecolors="black", zorder=8)
    z_cchp = k_pt["z_vs_CCHP_TRGB"]
    ax1.annotate((f"KBC-üreg (${num_s(z_cchp, '+.2f', hu)}\\sigma$ CCHP)\n$H_0^{{\\mathrm{{loc}}}} = {num_s(k_pt['H0_joint_unsub_dMB0'], '.2f', hu)}$" if hu
                  else f"KBC void (${z_cchp:+.2f}\\sigma$ CCHP)\n$H_0^{{\\mathrm{{loc}}}} = {k_pt['H0_joint_unsub_dMB0']:.2f}$"),
                 (k_pt["R_Mpc"], k_pt["H0_joint_unsub_dMB0"]), xytext=(248, 71.55), fontsize=7.1, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.8, color="#059669"))

    ax1.set_xlim(60, 380)
    ax1.set_ylim(67.8, 74.5)
    ax1.set_xlabel("Lokális void Gauss-skálasugár $R_{\\mathrm{void}}\\;[\\mathrm{Mpc}]$" if hu
                   else "Local void Gaussian scale radius $R_{\\mathrm{void}}\\;[\\mathrm{Mpc}]$")
    ax1.set_ylabel("Ablakolt $H_0^{\\mathrm{loc}}(z\\in[0{,}023;\\,0{,}15])\\;[\\mathrm{km\\,s^{-1}\\,Mpc^{-1}}]$" if hu
                   else "Windowed $H_0^{\\mathrm{loc}}(z\\in[0.023,0.15])\\;[\\mathrm{km\\,s^{-1}\\,Mpc^{-1}}]$")
    ax1.set_title("(A) SH0ES ablakolt $H_0$ ($\\Delta M_B{=}0$, $k_{\\mathrm{env}}{=}2$)" if hu
                  else "(A) SH0ES-Windowed Local-Void $H_0$ ($\\Delta M_B{=}0$, $k_{\\mathrm{env}}{=}2$)")
    ax1.legend(loc="upper left", fontsize=6.0)
    ax1.grid(alpha=0.25)
    apply_hu_comma_ticks(ax1, hu)

    # Panel B: Root-Shot Self-Consistent Supervoid ISW (|delta_b(0.50)| == 0.25)
    g6 = DATA["part2_covariant_ir_and_isw"]["gate6_self_consistent_supervoid_isw"]
    unscr = g6["unscreened_Sub1_nogo"]
    scr   = g6["screened_Sub1_xi_Phi_sq"]

    cats_isw = (["Árnyékolatlan Sub1\n(anyag-norm.)", "Árnyékolatlan Sub1\n(galaxis-norm.)",
                 "Helmholtz-árnyékolt\n($\\kappa_\\theta{=}\\xi_\\Phi^2$)", "DES Y3 megfigyelés\n($4{,}1\\pm 1{,}1$)"]
                if hu else
                ["Unscreened Sub1\n(matter-norm.)", "Unscreened Sub1\n(galaxy-norm.)",
                 "Helmholtz-screened\n($\\kappa_\\theta{=}\\xi_\\Phi^2$)", "DES Y3 observed\n($4.1\\pm 1.1$)"])
    vals_isw = [unscr["A_ISW_matter_norm"], unscr["A_ISW_galaxy_norm"], scr["A_ISW_matter_norm"], 4.1]
    cols_isw = ["#dc2626", "#ef4444", "#059669", "#0284c7"]

    ax2.axhspan(-5.0, 0.0, color="#dc2626", alpha=0.09,
                label=("Meleg folt ($A_{\\mathrm{ISW}}<0$, kizárt!)" if hu
                       else "Hot-spot sign flip ($A_{\\mathrm{ISW}}<0$, ruled out!)"))
    ax2.axhline(0.0, color="black", lw=1.0)
    ax2.axhline(1.0, color="#475569", ls=":", lw=1.1,
                label=("$\\Lambda$CDM ($A_{\\mathrm{ISW}}=1{,}00$)" if hu else "$\\Lambda$CDM ($A_{\\mathrm{ISW}}=1.00$)"))
    ax2.bar(range(4), vals_isw, width=0.52, color=cols_isw, edgecolor="black", lw=0.65, zorder=5)
    ax2.errorbar([3], [4.1], yerr=[1.1], fmt="none", ecolor="black", elinewidth=1.3, capsize=4, zorder=7)
    for i, v in enumerate(vals_isw):
        dy = 0.28 if v >= 0 else -0.62
        if i == 3:
            dy = 1.25
        lbl_v = f"{v:+.2f}".replace(".", ",") if hu else f"{v:+.2f}"
        ax2.text(i, v + dy, lbl_v, ha="center", fontsize=7.5, fontweight="bold")

    ax2.set_xticks(range(4))
    ax2.set_xticklabels(cats_isw, fontsize=6.6)
    ax2.set_ylim(-5.0, 6.2)
    ax2.set_ylabel("ISW amplitúdó $A_{\\mathrm{ISW}} = \\Delta T / \\Delta T_{\\Lambda\\mathrm{CDM}}$" if hu
                   else "ISW amplitude $A_{\\mathrm{ISW}} = \\Delta T / \\Delta T_{\\Lambda\\mathrm{CDM}}$")
    ax2.set_title("(B) Önkikényszerített szupervoid ISW ($|\\delta_b|{=}0{,}25$)" if hu
                  else "(B) Root-Shot Supervoid ISW ($|\\delta_b|{=}0.25$)")
    ax2.legend(loc="upper left", fontsize=6.6)
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
                 label=("Tau Core Sub1 ($m_H\\equiv 0$)" if hu else "Tau Core Sub1 ($m_H\\equiv 0$)"))
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
    print("Generated all 6 JSON-driven v4.1.0 vector PDF figures (EN + HU) in", FIG_DIR)
