# -*- coding: utf-8 -*-
"""Figures for the revised manuscript (PNG 600 dpi + editable SVG). Author: RNV."""
import os, sys
ROOT = (sys.argv[1].rstrip("/\\") + "/") if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/") + "/"
RAW = ROOT + "data/raw/"; OUTD = ROOT + "outputs/"; os.makedirs(OUTD + "figures/", exist_ok=True)
import pandas as pd, numpy as np, json, sys, warnings, itertools
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Ellipse
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
import statsmodels.formula.api as smf
warnings.filterwarnings("ignore")
DATA = OUTD; FIG = OUTD + "figures/"
import os; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "axes.labelsize": 10, "axes.titlesize": 10})
TEXT = "#1a1a1a"; TEXT2 = "#555555"; GRID = "#dcdcdc"
CAT = {"L": "#2a6fb0", "M": "#e08a2e", "H": "#4c9f70"}  # colour-blind safe (blue / orange / green)
CATLAB = {"L": "Low disturbance", "M": "Medium disturbance", "H": "High disturbance"}
LABS = list("ABCDEFGHI")
df = pd.read_csv(DATA + "individuals_clean.csv"); frag = pd.read_csv(DATA + "fragments.csv"); R = json.load(open(DATA + "results.json"))
fl = pd.read_csv(DATA + "fragment_level.csv")
def save(fig, name):
    fig.savefig(FIG + name + ".png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG + name + ".svg", bbox_inches="tight", facecolor="white"); plt.close(fig); print("saved", name)

# ---------------------------------------------------------------- Fig 1: causal structure diagram
def fig1():
    fig, ax = plt.subplots(figsize=(8.6, 5.0)); ax.set_xlim(0, 118); ax.set_ylim(-4, 66); ax.axis("off")
    def box(x, y, w, h, title, sub="", fc="#f4f3ef", ec="#666666", tfs=8.4):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=ec, lw=1.0))
        if sub:
            ax.text(x + w/2, y + h*0.70, title, ha="center", va="center", fontsize=tfs, color=TEXT, weight="bold")
            ax.text(x + w/2, y + h*0.27, sub, ha="center", va="center", fontsize=6.6, color=TEXT2, linespacing=1.1)
        else: ax.text(x + w/2, y + h/2, title, ha="center", va="center", fontsize=tfs, color=TEXT, weight="bold")
    def arrow(p1, p2, ls="-", color="#555555", lw=1.1, rad=0.0):
        ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=11, lw=lw, color=color, ls=ls, connectionstyle=f"arc3,rad={rad}"))
    for y, t in [(57, "Process"), (42.5, "Landscape" + "\n" + "pattern" + "\n" + "(measured)"), (27.5, "Local" + "\n" + "mechanisms" + "\n" + "(not measured)"), (10.5, "Phenotype" + "\n" + "(measured)")]:
        ax.text(0.5, y, t, ha="left", va="center", fontsize=7.2, color=TEXT2)
    ax.plot([17, 17], [2, 64], color=GRID, lw=0.8)
    box(31, 53, 48, 8, "Deforestation and land-use change", "historical and ongoing process in the Atlantic Forest", fc="#e9e8e3")
    box(92, 53, 24, 8, "Elevation", "40–1,000 m (covariate)", fc="#e9e8e3")
    box(20, 38, 24, 9, "Fragment area", "4.1–3,600 ha (log10)", fc="#e3eefb", ec=CAT["L"])
    box(47, 38, 26, 9, "Isolation", "distance to nearest forest," + "\n" + "0.05–15 km (log10)", fc="#e3eefb", ec=CAT["L"])
    box(76, 38, 40, 9, "Disturbance level", "ordinal L / M / H; confounded with area" + "\n" + "in this set of remnants (Spearman ρ = −0.94)", fc="#e3eefb", ec=CAT["L"])
    box(22, 23, 34, 9, "Floral and nesting resources", "quantity and quality around nests", fc="#fff3e6", ec=CAT["M"])
    box(62, 23, 40, 9, "Larval provisioning by females", "philopatric nesting; season and climate", fc="#fff3e6", ec=CAT["M"])
    box(30, 6, 62, 9, "Adult male body size (wing centroid size, CS)", "fixed at the end of larval development; measured in 211 males", fc="#e8f5ee", ec=CAT["H"])
    arrow((46, 53), (33, 47.2)); arrow((55, 53), (59, 47.2)); arrow((66, 53), (90, 47.2))
    arrow((31, 38), (34, 32.2)); ax.text(28.5, 35.3, "H1 (+)", fontsize=7, style="italic", color=TEXT2, ha="center")
    arrow((57, 38), (46, 32.2)); ax.text(44.5, 35.0, "H2 (−)", fontsize=7, style="italic", color=TEXT2, ha="center")
    arrow((94, 38), (86, 32.2)); ax.text(95.5, 35.3, "H3 (−)", fontsize=7, style="italic", color=TEXT2, ha="center")
    arrow((78, 38), (54, 32.2), ls=(0, (3, 2)), color="#9a9890"); ax.text(78.6, 36.2, "alt. (+)", fontsize=6.2, style="italic", color="#9a9890", ha="left", va="center")
    arrow((56, 27.5), (62, 27.5), color=CAT["M"]); arrow((40, 23), (52, 15.2)); arrow((80, 23), (70, 15.2))
    arrow((106, 53), (90, 10.5), ls=(0, (1, 2)), color="#9a9890", rad=-0.35); ax.text(111, 24, "temperature," + "\n" + "phenology", fontsize=6.3, style="italic", color="#9a9890", ha="center")
    arrow((18.5, 10.5), (30, 10.5), ls=(0, (3, 2)), color="#9a9890"); ax.text(18.5, 16.8, "dispersal between" + "\n" + "nearby fragments", fontsize=6.3, style="italic", color="#9a9890", ha="left", va="center")
    arrow((97, 21), (97, 15.2), ls=(0, (1, 2)), color="#9a9890"); ax.text(103, 18, "season of" + "\n" + "capture" + "\n" + "(covariate)", fontsize=6.3, style="italic", color="#9a9890", ha="left", va="center")
    ax.text(59, -3.5, "Solid arrows: hypothesised paths examined as associations between fragment attributes and CS (H1–H3), with fragment as a random effect." + "\n" +
            "Dashed and dotted arrows: alternative or confounding paths not tested with the present data (except season, included as a covariate).", ha="center", va="bottom", fontsize=6.6, color=TEXT2)
    save(fig, "Fig1_conceptual_diagram")

# ---------------------------------------------------------------- Fig 4: box-plots of CS per fragment
def fig4():
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    cld = R["cld_MW"]; rng = np.random.default_rng(1)
    for i, f in enumerate(LABS):
        y = df.cs[df.frag == f].values; c = CAT[frag.disturbance[frag.frag == f].iloc[0]]
        ax.boxplot(y, positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=c, alpha=0.25, edgecolor=c, lw=1.1), medianprops=dict(color=c, lw=1.8), whiskerprops=dict(color=c, lw=1), capprops=dict(color=c, lw=1))
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(y)), y, s=9, color=c, alpha=0.75, lw=0, zorder=3)
        ax.text(i, 318.5, cld[f], ha="center", va="bottom", fontsize=8.5, color=TEXT)
        ax.text(i, 238.5, f"n = {len(y)}", ha="center", va="top", fontsize=7, color=TEXT2)
    ax.set_xticks(range(9)); ax.set_xticklabels([f"{f}\n{a:g} ha" for f, a in zip(LABS, frag.area_ha)], fontsize=8)
    ax.set_ylim(232, 324); ax.set_xlim(-0.6, 8.6); ax.set_ylabel("Wing centroid size (CS)"); ax.set_xlabel("Forest fragment (ordered by decreasing area)")
    ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    hs = [plt.Line2D([], [], marker="s", ls="", ms=7, color=CAT[k], alpha=0.6, label=CATLAB[k]) for k in "LMH"]
    ax.legend(handles=hs, loc="lower left", frameon=False, fontsize=7.5, ncol=3, bbox_to_anchor=(0, -0.3))
    save(fig, "Fig4_boxplots_CS")

# ---------------------------------------------------------------- Fig 5: heatmap D2 (upper) + geographic distance (lower)
def fig5():
    D2 = pd.read_csv(DATA + "mahalanobis_D2_matrix.csv", index_col=0); geo = pd.read_csv(DATA + "geo_distances_table2.csv", index_col=0)
    hot = pd.read_csv(DATA + "shape_pairwise_hotelling_mahalanobis.csv"); sig = set(hot.pair[hot.p_holm < 0.05])
    cm_b = LinearSegmentedColormap.from_list("b", ["#deebf7", "#9ecae1", "#4292c6", "#08519c"]); cm_o = LinearSegmentedColormap.from_list("o", ["#fee6ce", "#fdae6b", "#e6550d", "#8c2d04"])
    fig, ax = plt.subplots(figsize=(6.6, 5.6)); k = 9
    mmin, mmax = np.nanmin(D2.values), np.nanmax(D2.values); gmax = geo.values.max()
    for i, r in enumerate(LABS):
        for j, c in enumerate(LABS):
            x, y = j, k - 1 - i
            if i == j:
                ax.add_patch(Rectangle((x, y), 1, 1, fc="#f2f1ed", ec="white", lw=2)); ax.text(x+.5, y+.5, r, ha="center", va="center", fontsize=10, weight="bold", color=TEXT)
            elif i < j:
                v = D2.loc[r, c]; t = (v - mmin) / (mmax - mmin); ax.add_patch(Rectangle((x, y), 1, 1, fc=cm_b(t), ec="white", lw=2))
                ax.text(x+.5, y+.5, f"{v:.1f}" + ("*" if f"{r}-{c}" in sig else ""), ha="center", va="center", fontsize=8, color="white" if t > 0.55 else TEXT)
            else:
                v = geo.loc[r, c]; t = v / gmax; ax.add_patch(Rectangle((x, y), 1, 1, fc=cm_o(t), ec="white", lw=2))
                ax.text(x+.5, y+.5, f"{v:g}" if v < 10 else f"{v:.0f}", ha="center", va="center", fontsize=8, color="white" if t > 0.55 else TEXT)
    ax.set_xlim(0, k); ax.set_ylim(0, k); ax.set_aspect("equal"); ax.set_xticks(np.arange(k)+.5); ax.set_xticklabels(LABS); ax.set_yticks(np.arange(k)+.5); ax.set_yticklabels(LABS[::-1])
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]; ax.xaxis.tick_top()
    import matplotlib.cm as mcm
    fig.subplots_adjust(left=0.06, right=0.80, top=0.95, bottom=0.05)
    for sm, pos, lab in [(mcm.ScalarMappable(cmap=cm_b, norm=plt.Normalize(mmin, mmax)), [0.84, 0.55, 0.025, 0.37], "Squared Mahalanobis distance (D²)\nwing shape — upper triangle"),
                         (mcm.ScalarMappable(cmap=cm_o, norm=plt.Normalize(0, gmax)), [0.84, 0.08, 0.025, 0.37], "Geographic distance (km)\nlower triangle")]:
        sm.set_array([]); cb = fig.colorbar(sm, cax=fig.add_axes(pos)); cb.set_label(lab, fontsize=7.5, color=TEXT2); cb.ax.tick_params(labelsize=7, colors=TEXT2); cb.outline.set_visible(False)
    save(fig, "Fig5_heatmap_D2_geo")

# ---------------------------------------------------------------- Fig 6: fragment area vs CS (LMM) + season
def fig6():
    m = smf.mixedlm("cs ~ log_area", df, groups=df["frag"]).fit(reml=True, method=["bfgs", "powell"], maxiter=2000)
    b = m.params[["Intercept", "log_area"]].values; V = m.cov_params().loc[["Intercept", "log_area"], ["Intercept", "log_area"]].values
    xs = np.linspace(0.5, 3.65, 100); X = np.column_stack([np.ones_like(xs), xs]); yhat = X @ b; se = np.sqrt(np.einsum("ij,jk,ik->i", X, V, X))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.6), gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax.fill_between(xs, yhat - 1.96*se, yhat + 1.96*se, color="#bbbbbb", alpha=0.35, lw=0); ax.plot(xs, yhat, color=TEXT, lw=1.4)
    for _, r in fl.iterrows():
        c = CAT[r.disturbance]; ax.errorbar(r.log_area, r.mean_cs, yerr=r.se, fmt="none", ecolor=c, elinewidth=1, capsize=2)
        ax.scatter(r.log_area, r.mean_cs, s=18 + 4.5*r.n, color=c, edgecolor="white", lw=0.8, zorder=4)
        ax.annotate(r.frag, (r.log_area, r.mean_cs), xytext=(6, 5), textcoords="offset points", fontsize=8, color=TEXT)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["10", "100", "1,000"]); ax.set_xlabel("Fragment area (ha, log scale)"); ax.set_ylabel("Mean wing centroid size (± SE)")
    ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.text(-0.13, 1.02, "a", transform=ax.transAxes, fontsize=11, weight="bold", va="bottom")
    ax.text(0.98, 0.04, f"LMM slope = {b[1]:.1f} CS units per 10-fold area\n95% CI {b[1]-1.96*np.sqrt(V[1,1]):.1f}–{b[1]+1.96*np.sqrt(V[1,1]):.1f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.3, color=TEXT2)
    hs = [plt.Line2D([], [], marker="o", ls="", ms=6, color=CAT[k], label=CATLAB[k]) for k in "LMH"]; ax.legend(handles=hs, loc="upper left", frameon=False, fontsize=7)
    # panel b: season within fragments sampled in both seasons
    sw = R["season_within"]; fr = [f for f in LABS if f in sw]; rng = np.random.default_rng(2)
    for i, f in enumerate(fr):
        for j, (s, mk) in enumerate([("cold-dry", "o"), ("hot-wet", "^")]):
            y = df.cs[(df.frag == f) & (df.season == s)].values; x = i + (-0.17 if j == 0 else 0.17)
            ax2.scatter(x + rng.uniform(-0.06, 0.06, len(y)), y, s=8, color="#888888" if j == 0 else "#c0392b", alpha=0.55, lw=0, marker=mk)
            ax2.errorbar(x, y.mean(), yerr=stats.sem(y) if len(y) > 1 else 0, fmt=mk, color="#444444" if j == 0 else "#c0392b", ms=5, capsize=2, zorder=4)
        ax2.text(i, 319.5, f"n = {sw[f]['n_dry']}/{sw[f]['n_wet']}", ha="center", fontsize=6.6, color=TEXT2)
    ax2.set_xticks(range(len(fr))); ax2.set_xticklabels(fr); ax2.set_xlabel("Fragment"); ax2.set_ylabel("Wing centroid size (CS)"); ax2.set_ylim(238, 324)
    hs = [plt.Line2D([], [], marker="o", ls="", ms=5, color="#444444", label="cold-dry season (Apr–Sep)"), plt.Line2D([], [], marker="^", ls="", ms=5, color="#c0392b", label="hot-wet season (Oct–Mar)")]
    ax2.legend(handles=hs, loc="lower left", bbox_to_anchor=(0, 0.04), frameon=False, fontsize=6.8); ax2.yaxis.grid(True, color=GRID, lw=0.6); ax2.set_axisbelow(True); ax2.text(-0.18, 1.02, "b", transform=ax2.transAxes, fontsize=11, weight="bold", va="bottom")
    fig.tight_layout(); save(fig, "Fig6_area_LMM_season")

# ---------------------------------------------------------------- Fig S1: CVA scatter with 90% ellipses
def figS1():
    cv = pd.read_csv(DATA + "cva_scores.csv"); fig, ax = plt.subplots(figsize=(6.4, 5.6))
    cols = dict(zip(LABS, ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#17becf"])); mk = dict(zip(LABS, "osD^vP*Xh"))
    for f in LABS:
        s = cv[cv.frag == f]; ax.scatter(s.CV1, s.CV2, s=16, color=cols[f], marker=mk[f], alpha=0.75, lw=0, label=f"{f} (n = {len(s)})")
        cov = np.cov(s.CV1, s.CV2); vals, vecs = np.linalg.eigh(cov); ang = np.degrees(np.arctan2(*vecs[:, 1][::-1])); w, h = 2*np.sqrt(vals*stats.chi2.ppf(0.90, 2))
        ax.add_patch(Ellipse((s.CV1.mean(), s.CV2.mean()), w, h, angle=ang, fc="none", ec=cols[f], lw=1.0, alpha=0.9))
    p = R["CVA"]["prop"]; ax.set_xlabel(f"Canonical variate 1 ({p[0]*100:.0f}% of among-group variance)"); ax.set_ylabel(f"Canonical variate 2 ({p[1]*100:.0f}%)")
    ax.axhline(0, color=GRID, lw=0.6); ax.axvline(0, color=GRID, lw=0.6); ax.legend(frameon=False, fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.14), title="Fragment (n)", title_fontsize=7.5)
    save(fig, "FigS1_CVA_scatter")

# ---------------------------------------------------------------- Fig S2: regression diagnostics (OLS individual + LMM)
def figS2():
    d = pd.read_csv(DATA + "ols_individual_diagnostics.csv"); mmd = pd.read_csv(DATA + "mm_area_diagnostics.csv")
    fig, axs = plt.subplots(2, 3, figsize=(8.4, 5.4)); axs = axs.ravel()
    axs[0].scatter(d.fitted, d.resid, s=8, color="#2a6fb0", alpha=0.6, lw=0); axs[0].axhline(0, color=TEXT2, lw=0.8); axs[0].set_xlabel("Fitted values"); axs[0].set_ylabel("Residuals"); axs[0].set_title("a  OLS (individual level): residuals vs fitted")
    stats.probplot(d.resid, dist="norm", plot=axs[1]); axs[1].get_lines()[0].set(markersize=3, color="#2a6fb0", alpha=0.6); axs[1].get_lines()[1].set(color=TEXT2); axs[1].set_title("b  OLS: normal Q–Q")
    axs[2].scatter(d.fitted, np.sqrt(np.abs(d.std_resid)), s=8, color="#2a6fb0", alpha=0.6, lw=0); axs[2].set_xlabel("Fitted values"); axs[2].set_ylabel("√|standardized residuals|"); axs[2].set_title("c  OLS: scale–location")
    axs[3].scatter(d.leverage, d.std_resid, s=8, color="#2a6fb0", alpha=0.6, lw=0); axs[3].set_xlabel("Leverage"); axs[3].set_ylabel("Standardized residuals"); axs[3].set_title("d  OLS: residuals vs leverage")
    axs[4].scatter(mmd.fitted, mmd.resid, s=8, color="#4c9f70", alpha=0.6, lw=0); axs[4].axhline(0, color=TEXT2, lw=0.8); axs[4].set_xlabel("Fitted values (fixed + random)"); axs[4].set_ylabel("Residuals"); axs[4].set_title("e  LMM (area): residuals vs fitted")
    stats.probplot(mmd.resid, dist="norm", plot=axs[5]); axs[5].get_lines()[0].set(markersize=3, color="#4c9f70", alpha=0.6); axs[5].get_lines()[1].set(color=TEXT2); axs[5].set_title("f  LMM: normal Q–Q")
    for a in axs: a.title.set_fontsize(8); a.title.set_ha("left"); a.title.set_position((0, 1))
    fig.tight_layout(); save(fig, "FigS2_regression_diagnostics")

# ---------------------------------------------------------------- Fig S3: CS by collection date within fragment
def figS3():
    d = df.copy(); d["date"] = pd.to_datetime(d.date); fig, axs = plt.subplots(3, 3, figsize=(8.4, 6.6), sharey=True); axs = axs.ravel(); rng = np.random.default_rng(3)
    for ax, f in zip(axs, LABS):
        s = d[d.frag == f]; dates = sorted(s.date.unique())
        for i, dt_ in enumerate(dates):
            y = s.cs[s.date == dt_].values; c = "#c0392b" if pd.Timestamp(dt_).month in [10, 11, 12, 1, 2, 3] else "#444444"
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(y)), y, s=9, color=c, alpha=0.6, lw=0); ax.errorbar(i, np.median(y), fmt="_", color=c, ms=14, mew=1.6)
        ax.set_xticks(range(len(dates))); ax.set_xticklabels([pd.Timestamp(x).strftime("%b\n%Y") for x in dates], fontsize=6.5)
        t = R["temporal_within_frag"].get(f); ax.set_title(f"{f}  (n = {len(s)})" + (f"  KW p = {t['p']:.2f}" if t else "  single survey"), fontsize=8, loc="left")
        ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    axs[3].set_ylabel("Wing centroid size (CS)"); hs = [plt.Line2D([], [], marker="o", ls="", ms=5, color="#444444", label="cold-dry season"), plt.Line2D([], [], marker="o", ls="", ms=5, color="#c0392b", label="hot-wet season")]
    fig.legend(handles=hs, loc="lower center", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.01)); fig.tight_layout(rect=(0, 0.03, 1, 1)); save(fig, "FigS3_CS_by_survey_date")

# ---------------------------------------------------------------- Fig 3: wing with numbered landmarks (detect white dots and match to mean shape)
def fig3():
    from PIL import Image; from scipy import ndimage
    im = Image.open(RAW + "wing_photo.png").convert("L"); A = np.array(im)
    mask = A > 235; lab, n = ndimage.label(mask); sizes = ndimage.sum(mask, lab, range(1, n+1)); cents = ndimage.center_of_mass(mask, lab, range(1, n+1))
    cand = [(c, s) for c, s in zip(cents, sizes) if 15 < s < 400]
    # keep blobs that are roughly round (bbox aspect)
    objs = ndimage.find_objects(lab); keep = []
    for idx, (c, s) in enumerate(zip(cents, sizes)):
        if not (15 < s < 400): continue
        sl = objs[idx]; h = sl[0].stop - sl[0].start; w = sl[1].stop - sl[1].start
        if 0.6 < h/w < 1.6 and s > 0.5*h*w: keep.append((c[1], c[0], s))
    pts = np.array([(x, y) for x, y, s in keep]); print("white blobs detected:", len(pts))
    ms = np.load(DATA + "meanshape.npy")  # 12x2, unit CS
    if len(pts) < 12: print("landmark detection failed"); return None
    # choose 12 points best matching mean shape by Procrustes (try all subsets is too many; use iterative: align using Hungarian after initial similarity fit)
    from scipy.optimize import linear_sum_assignment
    def fit(P, Q):  # similarity transform mapping Q (12x2) onto P (12x2), returns rmsd
        Pm = P - P.mean(0); Qm = Q - Q.mean(0); nP = np.linalg.norm(Pm); nQ = np.linalg.norm(Qm); u, s, vt = np.linalg.svd(Qm.T @ Pm); Rm = u @ vt
        if np.linalg.det(Rm) < 0: u[:, -1] *= -1; Rm = u @ vt
        sc = s.sum() / nQ**2; Qa = sc * Qm @ Rm + P.mean(0); return np.sqrt(((Qa - P)**2).sum() / len(P)), Qa
    best = None
    for refl in [np.array([[1, 0], [0, 1]]), np.array([[-1, 0], [0, 1]]), np.array([[1, 0], [0, -1]])]:
        M = ms @ refl
        # initial: scale mean shape to point cloud extent and align centroids, then assignment, then refine
        Pc = pts.copy(); Mc = M - M.mean(0); Mc = Mc / np.ptp(Mc[:, 0]) * np.ptp(Pc[:, 0]) + Pc.mean(0)
        for it in range(6):
            cost = np.linalg.norm(Pc[:, None, :] - Mc[None, :, :], axis=2); ri, ci = linear_sum_assignment(cost)
            sel = np.empty((12, 2)); sel[ci] = Pc[ri]
            rmsd, Qa = fit(sel, M); Mc = Qa
        if best is None or rmsd < best[0]: best = (rmsd, sel.copy())
    rmsd, sel = best; print("landmark match rmsd (px):", round(rmsd, 2))
    pass
    fig, ax = plt.subplots(figsize=(6.8, 2.2)); ax.imshow(A, cmap="gray", vmin=0, vmax=255); ax.axis("off")
    for i, (x, y) in enumerate(sel):
        ax.add_patch(plt.Circle((x, y), 9, fc="#ffd400", ec="black", lw=0.6, zorder=3)); ax.text(x, y, str(i+1), ha="center", va="center", fontsize=5.6, weight="bold", zorder=4)
    fig.savefig(FIG + "Fig3_wing_landmarks.png", dpi=600, bbox_inches="tight", pad_inches=0.02, facecolor="white"); fig.savefig(FIG + "Fig3_wing_landmarks.svg", bbox_inches="tight", pad_inches=0.02); plt.close(fig); print("saved Fig3_wing_landmarks")
    return sel

if __name__ == "__main__":
    fig1(); fig4(); fig5(); fig6(); figS1(); figS2(); figS3()
