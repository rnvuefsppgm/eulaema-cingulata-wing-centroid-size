# -*- coding: utf-8 -*-
"""Reanalysis of Eulaema cingulata wing centroid size data (JIC revision). Author: RNV.
Inputs: individuals_full.csv (built by 01_build_individual_dataset.py from data/raw/*.xlsx). Outputs: tables (CSV) + results.json"""
import os, sys
ROOT = (sys.argv[1].rstrip("/\\") + "/") if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/") + "/"
RAW = ROOT + "data/raw/"; OUTD = ROOT + "outputs/"; os.makedirs(OUTD + "figures/", exist_ok=True)
import pandas as pd, numpy as np, json, itertools, warnings, datetime as dt, openpyxl, sys
from scipy import stats
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.stattools import durbin_watson
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.multivariate.manova import MANOVA
warnings.filterwarnings("ignore")
rng = np.random.default_rng(2026)
OUT = OUTD; D = RAW; R = {}
LABS = list("ABCDEFGHI")
# ---------------- 0. dataset (rebuild dates without duplicates)
df = pd.read_csv(OUT + "individuals_full.csv").drop(columns=["date", "year", "month", "season"]).drop_duplicates("id_AI").reset_index(drop=True)
ws = openpyxl.load_workbook(D + "Data_coleta_El_cingulata_ morfometria.xlsx", data_only=True)["Plan1"]
dates = [(r[0], r[1]) for r in ws.iter_rows(values_only=True) if r[0] is not None and isinstance(r[1], dt.datetime)]
dd = pd.DataFrame(dates, columns=["id_orig", "date"])
dupl = dd[dd.id_orig.duplicated(keep=False)].sort_values("id_orig")
R["date_duplicates"] = dupl.assign(date=dupl.date.dt.strftime("%Y-%m-%d")).to_dict("records")
dd = dd.drop_duplicates("id_orig", keep="first")
df = df.merge(dd, on="id_orig", how="left")
df["date"] = pd.to_datetime(df.date); df["year"] = df.date.dt.year; df["month"] = df.date.dt.month
df["season"] = np.where(df.month.isin([10, 11, 12, 1, 2, 3]), "hot-wet", "cold-dry")
assert len(df) == 211 and df.date.notna().all()
frag = pd.DataFrame({"frag": LABS,
    "area_ha": [3600, 1200, 876, 630, 141, 200, 55, 6.6, 4.1], "isolation_km": [0.15, 1.70, 0.05, 15.0, 0.05, 0.57, 0.47, 0.09, 0.09],
    "elevation_m": [150, 40, 790, 400, 1000, 80, 60, 200, 200], "disturbance": list("LLMMMMHHH"),
    "lat": [-(22+25/60+44/3600), -(21+24/60+45/3600), -(22+5/60+51/3600), -(21+47/60+34/3600), -(22+4/60+43/3600), -(21+20/60+10/3600), -(21+19/60+10/3600), -(21+23/60+10/3600), -(21+23/60+30/3600)],
    "lon": [-(42+2/60+4/3600), -(41+5/60+10/3600), -(42+5/60+38/3600), -(41+26/60+30/3600), -(42+5/60+45/3600), -(41+13/60+30/3600), -(41+7/60+15/3600), -(41+55/60+15/3600), -(41+55/60+30/3600)]})
frag["log_area"] = np.log10(frag.area_ha); frag["log_iso"] = np.log10(frag.isolation_km * 1000); frag["elev100"] = frag.elevation_m / 100
frag["dist_ord"] = frag.disturbance.map({"L": 1, "M": 2, "H": 3})
df = df.merge(frag, on="frag")
df.to_csv(OUT + "individuals_clean.csv", index=False); frag.to_csv(OUT + "fragments.csv", index=False)
# ---------------- 1. descriptive (Table S2) and sampling records (Table S5)
def desc(g):
    q1, q3 = np.percentile(g, [25, 75])
    return pd.Series(dict(N=len(g), Mean=g.mean(), SD=g.std(ddof=1), SE=g.std(ddof=1)/np.sqrt(len(g)), Variance=g.var(ddof=1), Median=np.median(g), Q1=q1, Q3=q3, IQR=q3-q1, Min=g.min(), Max=g.max()))
t2 = df.groupby("frag").cs.apply(desc).unstack(); t2.loc["Total"] = desc(df.cs); t2 = t2.round(2); t2["N"] = t2.N.astype(int)
t2.to_csv(OUT + "cs_descriptive_by_fragment.csv"); R["tableS2"] = t2.reset_index().to_dict("records")
sub = pd.DataFrame({"frag": LABS + ["Total"], "Mean_sub": [294.76, 286.64, 284.89, 295.43, 292.73, 287.61, 273.42, 273.69, 267.09, 286.3], "SD_sub": [8.36, 10.38, 9.91, 9.36, 8.14, 6.56, 16.46, 18.11, 10.66, 13.75], "Min_sub": [275.06, 258.85, 288.09, 270, 268.95, 274.75, 249.12, 243.76, 243.76, 243.76], "Max_sub": [311.05, 304.41, 310.93, 315.08, 306.37, 297.14, 297.46, 312.58, 281.8, 315.08], "SE_sub": [1.53, 1.86, 1.81, 1.71, 1.49, 1.37, 4.10, 3.86, 3.07, 0.92]}).set_index("frag")
t2[["N", "Mean", "SD", "SE", "Min", "Max"]].join(sub).to_csv(OUT + "cs_descriptive_vs_submitted_table.csv")
t5 = df.groupby(["frag", "date", "season", "year"]).size().reset_index(name="males_analysed").sort_values(["frag", "date"])
t5["date"] = t5.date.dt.strftime("%d/%m/%Y"); t5.to_csv(OUT + "capture_records.csv", index=False)
R["tableS5"] = t5.to_dict("records")
R["season_by_frag"] = df.groupby(["frag", "season"]).size().unstack(fill_value=0).to_dict("index")
R["year_by_frag"] = df.groupby(["frag", "year"]).size().unstack(fill_value=0).to_dict("index")
# ---------------- 2. Kruskal-Wallis + post hoc (Mann-Whitney/Bonferroni and Dunn/Holm)
groups = [g.cs.values for _, g in df.groupby("frag")]
H, p = stats.kruskal(*groups); R["KW"] = dict(H=H, df=8, p=p, eps2=(H - 8) / (len(df) - 9))
lev = stats.levene(*groups, center="median"); R["levene"] = dict(W=lev.statistic, p=lev.pvalue)
R["shapiro_by_frag"] = df.groupby("frag").cs.apply(lambda g: stats.shapiro(g).pvalue).round(3).to_dict()
aov = stats.f_oneway(*groups); R["ANOVA_cs"] = dict(F=aov.statistic, p=aov.pvalue)
pairs = list(itertools.combinations(LABS, 2)); m = len(pairs)
mw = {}
for a, b in pairs:
    u, pu = stats.mannwhitneyu(df.cs[df.frag == a], df.cs[df.frag == b], alternative="two-sided")
    mw[(a, b)] = dict(U=u, p_raw=pu, p_bonf=min(1, pu * m))
df["rank"] = df.cs.rank(); N = len(df)
ties = df.cs.value_counts(); Cties = 1 - ((ties**3 - ties).sum()) / (N**3 - N)
dunn = {}
for a, b in pairs:
    ra, rb = df["rank"][df.frag == a], df["rank"][df.frag == b]
    z = (ra.mean() - rb.mean()) / np.sqrt((N * (N + 1) / 12) * Cties * (1 / len(ra) + 1 / len(rb)))
    dunn[(a, b)] = dict(z=z, p_raw=2 * stats.norm.sf(abs(z)))
def holm(p):
    p = np.array(p, dtype=float); mm_ = int(np.isfinite(p).sum()); adj = np.full(len(p), np.nan)
    order = np.argsort(np.where(np.isfinite(p), p, np.inf)); run = 0
    for i, ix in enumerate(order[:mm_]):
        run = max(run, p[ix] * (mm_ - i)); adj[ix] = min(1, run)
    return adj
adj = holm([dunn[k]["p_raw"] for k in pairs])
for k, a in zip(pairs, adj): dunn[k]["p_holm"] = a
ph = pd.DataFrame([dict(pair=f"{a}-{b}", MW_U=mw[(a, b)]["U"], MW_p_raw=mw[(a, b)]["p_raw"], MW_p_bonf=mw[(a, b)]["p_bonf"], Dunn_z=dunn[(a, b)]["z"], Dunn_p_raw=dunn[(a, b)]["p_raw"], Dunn_p_holm=dunn[(a, b)]["p_holm"]) for a, b in pairs])
ph.to_csv(OUT + "cs_pairwise_posthoc_tests.csv", index=False)
R["posthoc_sig_MW_bonf"] = ph.pair[ph.MW_p_bonf < 0.05].tolist(); R["posthoc_sig_Dunn_holm"] = ph.pair[ph.Dunn_p_holm < 0.05].tolist()
def cld(labs, sigpairs, centre):
    order = sorted(labs, key=lambda l: -centre[l]); groupsL = []
    for l in order:
        placed = False
        for gset in groupsL:
            if all((min(l, o), max(l, o)) not in sigpairs for o in gset): gset.add(l); placed = True
        if not placed: groupsL.append({l})
    groupsL = [g for i, g in enumerate(groupsL) if not any(g < h for j, h in enumerate(groupsL) if i != j)]
    letters = {l: "" for l in labs}
    for i, g in enumerate(groupsL):
        for l in sorted(g): letters[l] += chr(97 + i)
    return letters
med = df.groupby("frag").cs.median().to_dict()
R["cld_MW"] = cld(LABS, {tuple(p.split("-")) for p in R["posthoc_sig_MW_bonf"]}, med)
R["cld_Dunn"] = cld(LABS, {tuple(p.split("-")) for p in R["posthoc_sig_Dunn_holm"]}, med)
# ---------------- 3. individual-level OLS (original approach) with diagnostics
ols = smf.ols("cs ~ log_area + log_iso + elev100", data=df).fit()
Xd = sm.add_constant(df[["log_area", "log_iso", "elev100"]])
bp = het_breuschpagan(ols.resid, Xd); vif = {c: variance_inflation_factor(Xd.values, i) for i, c in enumerate(Xd.columns) if c != "const"}
R["OLS_ind"] = dict(F=ols.fvalue, df_num=int(ols.df_model), df_den=int(ols.df_resid), p=ols.f_pvalue, R2=ols.rsquared, adjR2=ols.rsquared_adj,
    coef=ols.params.to_dict(), se=ols.bse.to_dict(), pvals=ols.pvalues.to_dict(), ci=ols.conf_int().T.to_dict(),
    shapiro_W=stats.shapiro(ols.resid).statistic, shapiro_p=stats.shapiro(ols.resid).pvalue, BP_LM=bp[0], BP_p=bp[1], DW=durbin_watson(ols.resid), VIF=vif,
    r_cs_logarea=stats.pearsonr(df.cs, df.log_area)[0], r_p=stats.pearsonr(df.cs, df.log_area)[1], rho_cs_logarea=stats.spearmanr(df.cs, df.log_area)[0])
R["OLS_ind_single"] = {v: dict(r=stats.pearsonr(df.cs, df[v])[0], p=stats.pearsonr(df.cs, df[v])[1], R2=stats.pearsonr(df.cs, df[v])[0]**2) for v in ["log_area", "log_iso", "elev100"]}
alt_area = {"A": 3600, "B": 1200, "C": 876, "D": 600, "E": 180, "F": 150, "G": 45, "H": 5, "I": 2}
df["log_area_sheet"] = np.log10(df.frag.map(alt_area))
ols2 = smf.ols("cs ~ log_area_sheet + log_iso + elev100", data=df).fit(); R["OLS_ind_sheetareas"] = dict(F=ols2.fvalue, R2=ols2.rsquared, p=ols2.f_pvalue)
infl = ols.get_influence()
pd.DataFrame(dict(fitted=ols.fittedvalues, resid=ols.resid, std_resid=infl.resid_studentized_internal, leverage=infl.hat_matrix_diag, cooks=infl.cooks_distance[0], frag=df.frag)).to_csv(OUT + "ols_individual_diagnostics.csv", index=False)
# ---------------- 4. mixed models (fragment random intercept)
def mm(formula, reml=True):
    return smf.mixedlm(formula, df, groups=df["frag"]).fit(reml=reml, method=["bfgs", "powell", "nm"], maxiter=2000)
m0 = mm("cs ~ 1"); var_f = float(m0.cov_re.iloc[0, 0]); var_e = float(m0.scale); icc = var_f / (var_f + var_e)
R["MM_null"] = dict(intercept=m0.params["Intercept"], var_frag=var_f, var_resid=var_e, ICC=icc)
def mmres(m_):
    ci = m_.conf_int(); keep = [k for k in m_.params.index if k != "Group Var"]
    return dict(coef={k: m_.params[k] for k in keep}, se={k: m_.bse[k] for k in keep}, p={k: m_.pvalues[k] for k in keep}, ci={k: [ci.loc[k, 0], ci.loc[k, 1]] for k in keep},
                var_frag=float(m_.cov_re.iloc[0, 0]), var_resid=float(m_.scale), loglik=m_.llf, aic=m_.aic)
ma = mm("cs ~ log_area"); R["MM_area"] = mmres(ma); R["MM_area"]["prop_between_var_explained"] = 1 - float(ma.cov_re.iloc[0, 0]) / var_f
R["MM_iso"] = mmres(mm("cs ~ log_iso")); R["MM_elev"] = mmres(mm("cs ~ elev100")); R["MM_full"] = mmres(mm("cs ~ log_area + log_iso + elev100"))
R["MM_dist"] = mmres(mm("cs ~ dist_ord")); R["MM_area_season"] = mmres(mm("cs ~ log_area + season")); R["MM_area_year"] = mmres(mm("cs ~ log_area + C(year)"))
R["MM_area_dist"] = mmres(mm("cs ~ log_area + dist_ord"))
l0 = mm("cs ~ 1", reml=False); l1 = mm("cs ~ log_area", reml=False)
def lrt(big, small, dfd): lr = 2 * (big.llf - small.llf); return dict(LR=lr, df=dfd, p=stats.chi2.sf(lr, dfd))
R["LRT_area"] = lrt(l1, l0, 1); R["LRT_season"] = lrt(mm("cs ~ log_area + season", reml=False), l1, 1); R["LRT_year"] = lrt(mm("cs ~ log_area + C(year)", reml=False), l1, 1)
R["LRT_iso_elev_given_area"] = lrt(mm("cs ~ log_area + log_iso + elev100", reml=False), l1, 2); R["LRT_dist_given_area"] = lrt(mm("cs ~ log_area + dist_ord", reml=False), l1, 1)
fixed = ma.params["Intercept"] + ma.params["log_area"] * df.log_area; vf = np.var(fixed)
R["MM_area"]["R2_marginal"] = vf / (vf + float(ma.cov_re.iloc[0, 0]) + float(ma.scale)); R["MM_area"]["R2_conditional"] = (vf + float(ma.cov_re.iloc[0, 0])) / (vf + float(ma.cov_re.iloc[0, 0]) + float(ma.scale))
R["MM_area"]["resid_shapiro_p"] = stats.shapiro(ma.resid).pvalue
re_ = pd.Series({k: float(v.iloc[0]) for k, v in ma.random_effects.items()}); R["MM_area"]["ranef_shapiro_p"] = stats.shapiro(re_).pvalue; R["MM_area"]["ranef"] = re_.round(2).to_dict()
pd.DataFrame(dict(fitted=ma.fittedvalues, resid=ma.resid, frag=df.frag, log_area=df.log_area)).to_csv(OUT + "mm_area_diagnostics.csv", index=False)
# ---------------- 5. fragment-level regression (n = 9), conservative inference (t with 7 df)
fl = df.groupby("frag").agg(mean_cs=("cs", "mean"), sd=("cs", "std"), n=("cs", "size"), median=("cs", "median")).reset_index().merge(frag, on="frag")
fl["se"] = fl.sd / np.sqrt(fl.n)
def flfit(formula, w=None):
    mdl = smf.wls(formula, fl, weights=w).fit() if w is not None else smf.ols(formula, fl).fit()
    ci = mdl.conf_int(); ex = mdl.model.exog
    return dict(coef=mdl.params.to_dict(), se=mdl.bse.to_dict(), t=mdl.tvalues.to_dict(), p=mdl.pvalues.to_dict(), ci={k: [ci.loc[k, 0], ci.loc[k, 1]] for k in mdl.params.index}, R2=mdl.rsquared, adjR2=mdl.rsquared_adj, df_resid=int(mdl.df_resid),
                shapiro_p=stats.shapiro(mdl.resid).pvalue, BP_p=het_breuschpagan(mdl.resid, ex)[1], DW=durbin_watson(mdl.resid))
R["FL_area_ols"] = flfit("mean_cs ~ log_area"); R["FL_area_wls"] = flfit("mean_cs ~ log_area", w=fl.n)
R["FL_iso"] = flfit("mean_cs ~ log_iso"); R["FL_elev"] = flfit("mean_cs ~ elev100"); R["FL_dist"] = flfit("mean_cs ~ dist_ord"); R["FL_full"] = flfit("mean_cs ~ log_area + log_iso + elev100")
R["FL_area_dist"] = flfit("mean_cs ~ log_area + dist_ord")
Xf = sm.add_constant(fl[["log_area", "log_iso", "elev100"]]); R["FL_full"]["VIF"] = {c: variance_inflation_factor(Xf.values, i) for i, c in enumerate(Xf.columns) if c != "const"}
R["FL_pred_corr"] = fl[["log_area", "log_iso", "elev100", "dist_ord"]].corr(method="spearman").round(2).to_dict()
R["FL_spearman"] = {v: dict(rho=stats.spearmanr(fl.mean_cs, fl[v])[0], p=stats.spearmanr(fl.mean_cs, fl[v])[1]) for v in ["log_area", "log_iso", "elev100", "dist_ord"]}
R["FL_median_spearman_area"] = dict(rho=stats.spearmanr(fl["median"], fl.log_area)[0], p=stats.spearmanr(fl["median"], fl.log_area)[1])
fl.to_csv(OUT + "fragment_level.csv", index=False)
loo = {}
for f in LABS:
    mdl = smf.ols("mean_cs ~ log_area", fl[fl.frag != f]).fit(); loo[f] = dict(slope=mdl.params["log_area"], p=mdl.pvalues["log_area"], R2=mdl.rsquared)
R["FL_area_LOO"] = loo
# ---------------- 6. temporal variation within fragments
temp = {}
for f, g in df.groupby("frag"):
    if g.date.nunique() > 1:
        gg = [x.cs.values for _, x in g.groupby("date")]; h_, pv = stats.kruskal(*gg)
        temp[f] = dict(n_dates=int(g.date.nunique()), H=h_, df=int(g.date.nunique() - 1), p=pv, by_date={d.strftime("%d/%m/%Y"): dict(n=int(len(x)), mean=float(x.cs.mean()), median=float(x.cs.median())) for d, x in g.groupby("date")})
R["temporal_within_frag"] = temp
both = df.groupby("frag").season.nunique(); both = both[both > 1].index.tolist()
R["season_within"] = {f: dict(n_dry=int(((df.frag == f) & (df.season == "cold-dry")).sum()), n_wet=int(((df.frag == f) & (df.season == "hot-wet")).sum()),
    mean_dry=float(df.cs[(df.frag == f) & (df.season == "cold-dry")].mean()), mean_wet=float(df.cs[(df.frag == f) & (df.season == "hot-wet")].mean()),
    MW_p=float(stats.mannwhitneyu(df.cs[(df.frag == f) & (df.season == "cold-dry")], df.cs[(df.frag == f) & (df.season == "hot-wet")]).pvalue)) for f in both}
# ---------------- 7. shape: GPA, PCA, MANOVA/CVA, Hotelling, Mahalanobis, Mantel, allometry, Box's M, Mardia
X = df[[f"{a}{i}" for i in range(1, 13) for a in ("x", "y")]].values.reshape(-1, 12, 2).astype(float)
def gpa(X, iters=200):
    X = X - X.mean(axis=1, keepdims=True); cs = np.sqrt((X**2).sum(axis=(1, 2))); X = X / cs[:, None, None]
    ref = X[0].copy()
    for _ in range(iters):
        al = []
        for s in X:
            u, _, vt = np.linalg.svd(s.T @ ref); rmat = u @ vt
            if np.linalg.det(rmat) < 0: u[:, -1] *= -1; rmat = u @ vt
            al.append(s @ rmat)
        al = np.array(al); new = al.mean(axis=0); new /= np.sqrt((new**2).sum())
        if np.abs(new - ref).max() < 1e-12: ref = new; break
        ref = new
    return al, ref, cs
P, meanshape, cs_raw = gpa(X)
np.save(OUT + "meanshape.npy", meanshape)
R["cs_check_corr_with_sheet"] = float(np.corrcoef(cs_raw * 1e-5, df.cs)[0, 1]); R["cs_check_max_abs_diff"] = float(np.abs(cs_raw * 1e-5 - df.cs).max())
Pf = P.reshape(len(P), 24); grp = df.frag.values
pd.DataFrame(Pf, columns=[f"p{a}{i}" for i in range(1, 13) for a in ("x", "y")]).assign(frag=grp, id=df.id_AI.values).to_csv(OUT + "procrustes_coords.csv", index=False)
Pc = Pf - Pf.mean(axis=0); U, Sv, Vt = np.linalg.svd(Pc, full_matrices=False); ev = Sv**2 / (len(Pc) - 1); prop = ev / ev.sum()
k20 = 20; PCs = Pc @ Vt[:k20].T
R["PCA"] = dict(var_prop_first10=prop[:10].round(4).tolist(), cum95_npc=int(np.searchsorted(np.cumsum(prop), 0.95) + 1), cum_first20=float(np.cumsum(prop)[19]))
pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(k20)]).assign(frag=grp).to_csv(OUT + "shape_PCs.csv", index=False)
dfm = pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(k20)]); dfm["frag"] = grp
mv = MANOVA.from_formula(" + ".join(dfm.columns[:k20]) + " ~ frag", data=dfm).mv_test()
tab = mv.results["frag"]["stat"]; R["MANOVA_20PC"] = {idx: dict(value=float(tab.loc[idx, "Value"]), F=float(tab.loc[idx, "F Value"]), df_num=float(tab.loc[idx, "Num DF"]), df_den=float(tab.loc[idx, "Den DF"]), p=float(tab.loc[idx, "Pr > F"])) for idx in tab.index}
G = LABS; ng = {g: int((grp == g).sum()) for g in G}; cent = {g: PCs[grp == g].mean(axis=0) for g in G}
W = sum(((PCs[grp == g] - cent[g]).T @ (PCs[grp == g] - cent[g])) for g in G) / (len(PCs) - len(G))
B = sum(ng[g] * np.outer(cent[g] - PCs.mean(axis=0), cent[g] - PCs.mean(axis=0)) for g in G) / (len(G) - 1)
Winv = np.linalg.inv(W); evals, evecs = np.linalg.eig(Winv @ B); idx = np.argsort(-evals.real); evals = evals.real[idx]; evecs = evecs.real[:, idx]
# scale CV axes to unit within-group variance
cvv = evecs[:, :2]; cvv = cvv / np.sqrt(np.diag(cvv.T @ W @ cvv)); CV = PCs @ cvv
R["CVA"] = dict(eig=evals[:8].round(4).tolist(), prop=(evals[:8] / evals[:8].sum()).round(3).tolist())
pd.DataFrame(dict(CV1=CV[:, 0], CV2=CV[:, 1], frag=grp)).to_csv(OUT + "cva_scores.csv", index=False)
D2 = pd.DataFrame(np.nan, index=G, columns=G); hot = []; p_dim = k20; dfe = len(PCs) - len(G)
for a, b in itertools.combinations(G, 2):
    d = cent[a] - cent[b]; d2 = float(d @ Winv @ d); D2.loc[a, b] = D2.loc[b, a] = d2
    na, nb = ng[a], ng[b]; T2 = (na * nb / (na + nb)) * d2; F = T2 * (dfe - p_dim + 1) / (dfe * p_dim); pv = stats.f.sf(F, p_dim, dfe - p_dim + 1)
    Sp = ((PCs[grp == a] - cent[a]).T @ (PCs[grp == a] - cent[a]) + (PCs[grp == b] - cent[b]).T @ (PCs[grp == b] - cent[b])) / (na + nb - 2)
    dfe2 = na + nb - 2
    if dfe2 - p_dim + 1 > 0 and np.linalg.cond(Sp) < 1e12:
        d2p = float(d @ np.linalg.solve(Sp, d)); T2b = (na * nb / (na + nb)) * d2p; F2 = T2b * (dfe2 - p_dim + 1) / (dfe2 * p_dim); pv2 = stats.f.sf(F2, p_dim, dfe2 - p_dim + 1)
    else: d2p = pv2 = np.nan
    hot.append(dict(pair=f"{a}-{b}", D2_pooled_all=d2, T2=T2, F=F, p_raw=pv, D2_pairpooled=d2p, p_pairpooled_raw=pv2))
hot = pd.DataFrame(hot); hot["p_holm"] = holm(hot.p_raw); hot["p_pairpooled_holm"] = holm(hot.p_pairpooled_raw)
hot.to_csv(OUT + "shape_pairwise_hotelling_mahalanobis.csv", index=False); D2.round(3).to_csv(OUT + "mahalanobis_D2_matrix.csv")
R["Hotelling_sig_pooled_holm"] = hot.pair[hot.p_holm < 0.05].tolist(); R["Hotelling_sig_pairpooled_holm"] = hot.pair[hot.p_pairpooled_holm < 0.05].tolist(); R["Hotelling_sig_pooled_raw"] = hot.pair[hot.p_raw < 0.05].tolist()
R["D2_range"] = dict(min=float(np.nanmin(D2.values)), max=float(np.nanmax(D2.values)), HI=float(D2.loc["H", "I"]), argmin=hot.loc[hot.D2_pooled_all.idxmin(), "pair"], argmax=hot.loc[hot.D2_pooled_all.idxmax(), "pair"], pairpooled_nan=hot.pair[hot.D2_pairpooled.isna()].tolist())
geo_sub = {("B","A"):150,("C","A"):37,("C","B"):127.7,("D","A"):95,("D","B"):56,("D","C"):80,("E","A"):38,("E","B"):127.9,("E","C"):0.17,("E","D"):74,("F","A"):148,("F","B"):16,("F","C"):128,("F","D"):56.5,("F","E"):125,("G","A"):156,("G","B"):11,("G","C"):133,("G","D"):62.5,("G","E"):135,("G","F"):10.5,("H","A"):116,("H","B"):86,("H","C"):79,("H","D"):67,("H","E"):78,("H","F"):72.5,("H","G"):83,("I","A"):115,("I","B"):87,("I","C"):78,("I","D"):66,("I","E"):77,("I","F"):73,("I","G"):84,("I","H"):0.7}
geo = pd.DataFrame(0.0, index=G, columns=G)
for (a, b), v in geo_sub.items(): geo.loc[a, b] = geo.loc[b, a] = v
def hav(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, [la1, lo1, la2, lo2]); a = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2; return 6371 * 2 * np.arcsin(np.sqrt(a))
hv = pd.DataFrame([[hav(frag.lat[i], frag.lon[i], frag.lat[j], frag.lon[j]) for j in range(9)] for i in range(9)], index=G, columns=G)
R["geo_check_max_abs_diff_km"] = float(np.abs(hv.values - geo.values).max()); R["geo_check_pairs_diff_gt10"] = [(a, b, round(geo.loc[a, b], 1), round(hv.loc[a, b], 1)) for a, b in itertools.combinations(G, 2) if abs(geo.loc[a, b] - hv.loc[a, b]) > 10]
hv.round(1).to_csv(OUT + "geo_distances_haversine.csv"); geo.to_csv(OUT + "geo_distances_table2.csv")
def mantel(A, Bm, nperm=9999):
    iu = np.triu_indices(len(A), 1); x = A[iu]; y = Bm[iu]; r = np.corrcoef(x, y)[0, 1]; cnt = 0
    for _ in range(nperm):
        pi = rng.permutation(len(A)); cnt += abs(np.corrcoef(x, Bm[np.ix_(pi, pi)][iu])[0, 1]) >= abs(r)
    return dict(r=float(r), p=(cnt + 1) / (nperm + 1))
R["Mantel_D2_geo"] = mantel(D2.values, geo.values); R["Mantel_D2_loggeo"] = mantel(D2.values, np.log10(geo.values + 1))
dcs = np.abs(np.subtract.outer(fl.mean_cs.values, fl.mean_cs.values)); R["Mantel_dCS_geo"] = mantel(dcs, geo.values); R["Mantel_D2_dCS"] = mantel(D2.values, dcs)
def permanova(Y, g, nperm=4999):
    from scipy.spatial.distance import pdist, squareform
    Dm = squareform(pdist(Y)); n = len(Y)
    def ss(gv):
        tot = (Dm[np.triu_indices(n, 1)]**2).sum() / n; w = 0
        for lv in set(gv):
            ix = np.where(gv == lv)[0]; w += (Dm[np.ix_(ix, ix)][np.triu_indices(len(ix), 1)]**2).sum() / len(ix)
        return tot, w
    tot, w = ss(g); a = len(set(g)); F0 = ((tot - w) / (a - 1)) / (w / (n - a)); cnt = 0
    for _ in range(nperm):
        t, w2 = ss(rng.permutation(g)); cnt += ((t - w2) / (a - 1)) / (w2 / (n - a)) >= F0
    return dict(pseudoF=F0, p=(cnt + 1) / (nperm + 1), R2=(tot - w) / tot, df_num=a - 1, df_den=n - a)
R["PERMANOVA_shape"] = permanova(Pf, grp)
def boxM(Y, g):
    G_ = sorted(set(g)); p_ = Y.shape[1]; n = len(Y); k = len(G_); Sp = 0; lnS = 0; c = 0
    for lv in G_:
        y = Y[g == lv]; ni = len(y); S = np.cov(y.T); Sp = Sp + (ni - 1) * S; lnS += (ni - 1) * np.log(np.linalg.det(S)); c += 1 / (ni - 1)
    Sp = Sp / (n - k); M = (n - k) * np.log(np.linalg.det(Sp)) - lnS
    Cc = (2 * p_**2 + 3 * p_ - 1) / (6 * (p_ + 1) * (k - 1)) * (c - 1 / (n - k)); chi = M * (1 - Cc); dfb = p_ * (p_ + 1) * (k - 1) / 2
    return dict(M=M, chi2=chi, df=dfb, p=stats.chi2.sf(chi, dfb))
R["BoxM"] = {f"first{kk}PC": boxM(PCs[:, :kk], grp) for kk in [2, 4, 6, 8]}
def mardia(Y):
    n, p_ = Y.shape; Yc = Y - Y.mean(axis=0); S = np.cov(Yc.T, bias=True); Si = np.linalg.inv(S); Dm = Yc @ Si @ Yc.T
    b1 = (Dm**3).sum() / n**2; b2 = (np.diag(Dm)**2).mean()
    skew_stat = n * b1 / 6; df1 = p_ * (p_ + 1) * (p_ + 2) / 6; kurt_z = (b2 - p_ * (p_ + 2)) / np.sqrt(8 * p_ * (p_ + 2) / n)
    return dict(b1=b1, skew_chi2=skew_stat, skew_df=df1, skew_p=stats.chi2.sf(skew_stat, df1), b2=b2, kurt_z=kurt_z, kurt_p=2 * stats.norm.sf(abs(kurt_z)))
resid = np.vstack([PCs[grp == g] - cent[g] for g in G]); R["Mardia_pooled_resid_20PC"] = mardia(resid); R["Mardia_pooled_resid_8PC"] = mardia(resid[:, :8])
R["Mardia_by_group_4PC"] = {g: mardia(PCs[grp == g][:, :4]) for g in G}
R["shape_univariate_shapiro_resid_frac_sig"] = float(np.mean([stats.shapiro(resid[:, j]).pvalue < 0.05 for j in range(20)]))
R["sample_size"] = dict(n_total=211, p_shape_dims=20, n_per_group=ng, pooled_df=dfe)
def proc_anova(Y, x, nperm=4999):
    Xm = np.column_stack([np.ones(len(x)), x]); beta = np.linalg.lstsq(Xm, Y, rcond=None)[0]; fit = Xm @ beta; sst = ((Y - Y.mean(axis=0))**2).sum(); ssr = ((fit - Y.mean(axis=0))**2).sum()
    F0 = (ssr / 1) / ((sst - ssr) / (len(x) - 2)); cnt = 0
    for _ in range(nperm):
        Xp = np.column_stack([np.ones(len(x)), rng.permutation(x)]); ssrp = ((Xp @ np.linalg.lstsq(Xp, Y, rcond=None)[0] - Y.mean(axis=0))**2).sum(); cnt += (ssrp / 1) / ((sst - ssrp) / (len(x) - 2)) >= F0
    return dict(R2=ssr / sst, F=F0, p=(cnt + 1) / (nperm + 1))
lcs = np.log10(df.cs.values); R["allometry_overall"] = proc_anova(Pf, lcs)
Yw = np.vstack([Pf[grp == g] - Pf[grp == g].mean(axis=0) for g in G]); xw = np.concatenate([lcs[grp == g] - lcs[grp == g].mean() for g in G]); R["allometry_within"] = proc_anova(Yw, xw)
dfm["lcs"] = lcs; mv2 = MANOVA.from_formula(" + ".join(dfm.columns[:k20]) + " ~ lcs + frag", data=dfm).mv_test()
R["MANCOVA_frag_given_lcs"] = dict(wilks=float(mv2.results["frag"]["stat"].loc["Wilks' lambda", "Value"]), p=float(mv2.results["frag"]["stat"].loc["Wilks' lambda", "Pr > F"]), wilks_lcs=float(mv2.results["lcs"]["stat"].loc["Wilks' lambda", "Value"]), p_lcs=float(mv2.results["lcs"]["stat"].loc["Wilks' lambda", "Pr > F"]))
# shape-group effect size: proportion of Procrustes variance among groups
R["shape_eta2_groups"] = float(((np.vstack([np.tile(Pf[grp == g].mean(axis=0), (ng[g], 1)) for g in G]) - Pf.mean(axis=0))**2).sum() / ((Pf - Pf.mean(axis=0))**2).sum())
# ---------------- save
def conv(o):
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (pd.Timestamp, dt.datetime)): return o.strftime("%Y-%m-%d")
    if isinstance(o, dict): return {str(k): conv(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [conv(v) for v in o]
    return o
json.dump(conv(R), open(OUT + "results.json", "w"), indent=1)
r2 = lambda d, n=3: {k: (round(v, n) if isinstance(v, (float, np.floating)) else v) for k, v in d.items()}
print("Table S2 recomputed:\n", t2[["N", "Mean", "SD", "SE", "Median", "IQR", "Min", "Max"]].to_string())
print("date duplicates:", R["date_duplicates"])
print("KW:", r2(R["KW"], 4), "Levene p", round(lev.pvalue, 4), "ANOVA", r2(R["ANOVA_cs"], 4), "Shapiro by frag", R["shapiro_by_frag"])
print("MW-Bonf sig:", R["posthoc_sig_MW_bonf"]); print("Dunn-Holm sig:", R["posthoc_sig_Dunn_holm"]); print("CLD MW:", R["cld_MW"]); print("CLD Dunn:", R["cld_Dunn"])
o = R["OLS_ind"]; print("OLS ind: F(%d,%d)=%.2f p=%.2e R2=%.3f coef=%s p=%s SW p=%.3g BP p=%.3g DW=%.2f VIF=%s r=%.3f rho=%.3f" % (o["df_num"], o["df_den"], o["F"], o["p"], o["R2"], r2(o["coef"], 2), r2(o["pvals"], 4), o["shapiro_p"], o["BP_p"], o["DW"], r2(o["VIF"], 2), o["r_cs_logarea"], o["rho_cs_logarea"]))
print("OLS single:", {k: r2(v) for k, v in R["OLS_ind_single"].items()}); print("OLS sheet areas:", r2(R["OLS_ind_sheetareas"]))
print("MM null ICC=%.3f varF=%.1f varE=%.1f" % (icc, var_f, var_e))
for k in ["MM_area", "MM_iso", "MM_elev", "MM_full", "MM_dist", "MM_area_dist", "MM_area_season", "MM_area_year"]:
    v = R[k]; print(k, "coef", r2(v["coef"], 2), "se", r2(v["se"], 2), "p", r2(v["p"], 4), "varF", round(v["var_frag"], 1), "AIC", round(v["aic"], 1))
print("MM_area CI", {k: [round(x, 2) for x in v] for k, v in R["MM_area"]["ci"].items()}, "R2m", round(R["MM_area"]["R2_marginal"], 3), "R2c", round(R["MM_area"]["R2_conditional"], 3), "between var expl", round(R["MM_area"]["prop_between_var_explained"], 3), "resid SW p", round(R["MM_area"]["resid_shapiro_p"], 3), "ranef SW p", round(R["MM_area"]["ranef_shapiro_p"], 3))
print("LRT:", {k: r2(R[k], 4) for k in ["LRT_area", "LRT_season", "LRT_year", "LRT_iso_elev_given_area", "LRT_dist_given_area"]})
for k in ["FL_area_ols", "FL_area_wls", "FL_iso", "FL_elev", "FL_dist", "FL_full", "FL_area_dist"]:
    v = R[k]; print(k, "coef", r2(v["coef"], 2), "se", r2(v["se"], 2), "p", r2(v["p"], 4), "R2", round(v["R2"], 3), "SW p", round(v["shapiro_p"], 3), "BP p", round(v["BP_p"], 3), "DW", round(v["DW"], 2))
print("FL area CI", [round(x, 2) for x in R["FL_area_ols"]["ci"]["log_area"]], "VIF", r2(R["FL_full"]["VIF"], 2)); print("FL spearman", {k: (round(v["rho"], 2), round(v["p"], 3)) for k, v in R["FL_spearman"].items()}); print("pred corr", R["FL_pred_corr"])
print("LOO slopes:", {k: (round(v["slope"], 2), round(v["p"], 3), round(v["R2"], 2)) for k, v in loo.items()})
print("temporal within:", {k: (v["n_dates"], round(v["p"], 3)) for k, v in temp.items()}); print("season within:", {k: (v["n_dry"], v["n_wet"], round(v["mean_dry"], 1), round(v["mean_wet"], 1), round(v["MW_p"], 3)) for k, v in R["season_within"].items()})
print("CS check corr", R["cs_check_corr_with_sheet"], "maxdiff", R["cs_check_max_abs_diff"]); print("PCA", R["PCA"])
print("MANOVA Wilks:", r2(R["MANOVA_20PC"]["Wilks' lambda"], 4), "Pillai p", round(R["MANOVA_20PC"]["Pillai's trace"]["p"], 5)); print("CVA eig prop", R["CVA"]["prop"][:4], "eta2 shape", round(R["shape_eta2_groups"], 3))
print("Hotelling pooled Holm sig:", R["Hotelling_sig_pooled_holm"]); print("Hotelling pooled raw sig:", R["Hotelling_sig_pooled_raw"]); print("Hotelling pair-pooled Holm sig:", R["Hotelling_sig_pairpooled_holm"]); print("D2:", r2(R["D2_range"])); print(D2.round(2).to_string())
print("Mantel D2-geo", r2(R["Mantel_D2_geo"], 4), "log", r2(R["Mantel_D2_loggeo"], 4), "dCS-geo", r2(R["Mantel_dCS_geo"], 4), "D2-dCS", r2(R["Mantel_D2_dCS"], 4)); print("geo check", round(R["geo_check_max_abs_diff_km"], 1), R["geo_check_pairs_diff_gt10"])
print("PERMANOVA", r2(R["PERMANOVA_shape"], 4)); print("BoxM", {k: (round(v["chi2"], 1), int(v["df"]), round(v["p"], 4)) for k, v in R["BoxM"].items()})
print("Mardia pooled 20PC", r2(R["Mardia_pooled_resid_20PC"], 4)); print("Mardia pooled 8PC", r2(R["Mardia_pooled_resid_8PC"], 4))
print("Mardia by group 4PC p(skew,kurt):", {g: (round(v["skew_p"], 3), round(v["kurt_p"], 3)) for g, v in R["Mardia_by_group_4PC"].items()}); print("univ shapiro frac sig", R["shape_univariate_shapiro_resid_frac_sig"])
print("allometry overall", r2(R["allometry_overall"], 4), "within", r2(R["allometry_within"], 4)); print("MANCOVA", r2(R["MANCOVA_frag_given_lcs"], 5))
