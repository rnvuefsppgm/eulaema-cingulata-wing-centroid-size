# -*- coding: utf-8 -*-
"""Complementary analyses: pairwise permutation tests for shape, D2 small-sample note, abundance not available. Author: RNV."""
import os, sys
ROOT = (sys.argv[1].rstrip("/\\") + "/") if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/") + "/"
RAW = ROOT + "data/raw/"; OUTD = ROOT + "outputs/"; os.makedirs(OUTD + "figures/", exist_ok=True)
import pandas as pd, numpy as np, json, itertools, sys
DATA = OUTD; rng = np.random.default_rng(7)
P = pd.read_csv(DATA + "procrustes_coords.csv"); Y = P[[c for c in P.columns if c.startswith("p")]].values; g = P.frag.values; LABS = list("ABCDEFGHI")
R = json.load(open(DATA + "results.json"))
# pairwise permutation test on Procrustes distance between group means (distribution-free), Holm
def pd_mean(a, b): return np.linalg.norm(Y[g == a].mean(0) - Y[g == b].mean(0))
res = []
for a, b in itertools.combinations(LABS, 2):
    ix = np.where((g == a) | (g == b))[0]; lab = g[ix].copy(); d0 = pd_mean(a, b); cnt = 0
    Ysub = Y[ix]
    for _ in range(4999):
        lp = rng.permutation(lab); d = np.linalg.norm(Ysub[lp == a].mean(0) - Ysub[lp == b].mean(0)); cnt += d >= d0
    res.append(dict(pair=f"{a}-{b}", procrustes_dist=d0, p_perm=(cnt + 1) / 5000))
res = pd.DataFrame(res); p = res.p_perm.values; order = np.argsort(p); adj = np.empty(len(p)); run = 0; m = len(p)
for i, ix in enumerate(order): run = max(run, p[ix] * (m - i)); adj[ix] = min(1, run)
res["p_perm_holm"] = adj
hot = pd.read_csv(DATA + "shape_pairwise_hotelling_mahalanobis.csv").merge(res, on="pair"); hot.to_csv(DATA + "shape_pairwise_hotelling_mahalanobis.csv", index=False)
R["perm_pairwise_sig_holm"] = res.pair[res.p_perm_holm < 0.05].tolist(); R["perm_pairwise_sig_raw"] = res.pair[res.p_perm < 0.05].tolist()
# agreement between methods
sets = {"hotelling_pooled_holm": set(R["Hotelling_sig_pooled_holm"]), "perm_holm": set(R["perm_pairwise_sig_holm"]), "pairpooled_holm": set(R["Hotelling_sig_pairpooled_holm"])}
R["pairwise_shape_agreement"] = {"both_hotelling_and_perm": sorted(sets["hotelling_pooled_holm"] & sets["perm_holm"]), "only_hotelling": sorted(sets["hotelling_pooled_holm"] - sets["perm_holm"]), "only_perm": sorted(sets["perm_holm"] - sets["hotelling_pooled_holm"])}
# D2 bias note: expected inflation ~ p*(1/na+1/nb) under H0
ng = pd.Series(g).value_counts().to_dict(); R["D2_H0_expectation"] = {f"{a}-{b}": 20 * (1/ng[a] + 1/ng[b]) for a, b in itertools.combinations(LABS, 2)}
# fragment mean CS vs D2 to A (largest) etc. not needed. Save
json.dump(R, open(DATA + "results.json", "w"), indent=1)
print("perm sig Holm:", R["perm_pairwise_sig_holm"]); print("agreement:", R["pairwise_shape_agreement"])
print("D2 H0 expectation H-I:", round(R["D2_H0_expectation"]["H-I"], 2), "A-B:", round(R["D2_H0_expectation"]["A-B"], 2))
