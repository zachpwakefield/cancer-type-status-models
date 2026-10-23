"""D-77 Phase 2 (scripts/62 compute/report): the statistics, the permutation schemes, the guards and the
claim-level verdict, checked on synthetic data against brute-force implementations of the
pre-registered definitions (config/followup.yaml: msi_rising_proof).

Run: cd tests && PYTHONPATH=../scripts:../src python -m unittest test_proof_stats
"""
from __future__ import annotations

import importlib.util
import os
import re
import sys
import unittest
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import torch  # noqa: E402,F401  (load torch's OpenMP runtime before numpy's; see scripts/_common.py)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
PRE = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["msi_rising_proof"]


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), REPO / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = load_script("62_msi_rising_proof.py")


# ---------------------------------------------------------------- synthetic data and brute force

def synth(seed=1, cancers=(4, 1, 3, 2), n_tum=(30, 25, 40, 12), n_exon=9, miss=0.25):
    """Per cancer an adj matrix (NaN = not measured) and observed MMRd labels."""
    rng = np.random.default_rng(seed)
    As, ys = [], []
    for nd, n in zip(cancers, n_tum):
        A = rng.normal(0, 0.05, (n, n_exon)) + rng.normal(0, 0.02, n_exon)
        A[rng.random((n, n_exon)) < miss] = np.nan
        y = np.zeros(n, bool)
        y[rng.choice(n, nd, replace=False)] = True
        A[y] += 0.01
        As.append(A)
        ys.append(y)
    return As, ys


def brute_exon(As, ys, e, min_mmr=5):
    """delta and its studentized value for one exon, straight from the definitions."""
    meas = [np.isfinite(A[:, e]) for A in As]
    use = [(m & y).any() and (m & ~y).any() for m, y in zip(meas, ys)]
    N = sum(int((m & y).sum()) for m, y, u in zip(meas, ys, use) if u)
    if N < min_mmr:
        return np.nan, np.nan, N
    num, u_all, Kc, vp = 0.0, [], 0, 0.0
    for A, m, y, u in zip(As, meas, ys, use):
        if not u:
            continue
        a = A[:, e]
        d, p = a[m & y], a[m & ~y]
        num += (d - p.mean()).sum()
        uu = d / N
        u_all.append(uu - uu.mean())
        Kc += 1
        v = p * len(d) / (N * len(p))
        if len(p) > 1:
            vp += len(p) * v.var(ddof=1)
    delta = num / N
    if N - Kc < 1:
        return delta, np.nan, N
    s2d = np.square(np.concatenate(u_all)).sum() / (N - Kc)
    return delta, delta / np.sqrt(N * s2d + vp), N


def brute_class(As, ys, min_mmr=5, weights=None):
    """RP2a's studentized class statistic from explicit per-tumour contributions (N = every MMRd
    tumour), with the Welch-Satterthwaite df."""
    E = As[0].shape[1]
    ex = [brute_exon(As, ys, e, min_mmr) for e in range(E)]
    valid = np.array([np.isfinite(d) for d, _, _ in ex])
    Ns = np.array([n for _, _, n in ex], float)
    w = np.where(valid, 1 / valid.sum(), 0) if weights is None else weights(valid[None, :], None)[0]
    T, U, Vc, comps = 0.0, [], 0.0, []
    for A, y in zip(As, ys):
        M_ = np.isfinite(A)
        nD, nP = (M_ & y[:, None]).sum(0), (M_ & ~y[:, None]).sum(0)
        ok = (nD > 0) & (nP > 0)
        u = np.array([sum(w[e] * A[t, e] / Ns[e] for e in range(E) if valid[e] and ok[e] and M_[t, e]) for t in np.flatnonzero(y)])
        v = np.array([sum(w[e] * A[t, e] * nD[e] / (Ns[e] * nP[e]) for e in range(E) if valid[e] and ok[e] and M_[t, e])
                      for t in np.flatnonzero(~y)])
        if y.any():
            T += u.sum() - v.sum()
            U.append(u - u.mean())
            if len(v) > 1:
                comps.append((len(v) * v.var(ddof=1), len(v) - 1))
    N, K = sum(int(y.sum()) for y in ys), sum(bool(y.any()) for y in ys)
    s2d = np.square(np.concatenate(U)).sum() / (N - K)
    se2 = N * s2d + sum(c for c, _ in comps)
    df = se2 ** 2 / ((N * s2d) ** 2 / (N - K) + sum(c ** 2 / d for c, d in comps))
    mean_delta = float(np.sum([w[e] * ex[e][0] for e in range(E) if valid[e]]))
    return T, se2, N, K, df, mean_delta


class TestStudentized(unittest.TestCase):
    """inference.studentized and inference.interval against brute force."""

    def test_exon_engine_matches_definition(self):
        As, ys = synth(seed=3)
        D, T, N = M.exon_engine(As, [y[:, None] for y in ys])
        for e in range(As[0].shape[1]):
            d, t, n = brute_exon(As, ys, e)
            self.assertEqual(N[0, e], n)
            if np.isfinite(d):
                self.assertAlmostEqual(float(D[0, e]), d, places=6)
                self.assertTrue(np.isfinite(T[0, e]) == np.isfinite(t))
                if np.isfinite(t):
                    self.assertAlmostEqual(float(T[0, e]), t, places=4)
            else:
                self.assertTrue(np.isnan(D[0, e]))

    def test_exon_engine_drops_cancers_without_a_measured_mmrp_reference(self):
        As, ys = synth(seed=4)
        As[0][~ys[0], 0] = np.nan                       # cancer 0: no MMRp measured at exon 0
        D, T, N = M.exon_engine(As, [y[:, None] for y in ys])
        d, t, n = brute_exon(As, ys, 0)
        self.assertEqual(N[0, 0], n)
        self.assertEqual(n, sum(int((np.isfinite(A[:, 0]) & y).sum()) for A, y in zip(As[1:], ys[1:])))

    def test_class_stat_se_against_brute_force_with_n_minus_k(self):
        As, ys = synth(seed=5, n_exon=12)
        out = M.class_stat(As, [y[:, None] for y in ys])
        T, se2, N, K, df, md = brute_class(As, ys)
        self.assertAlmostEqual(out["T"][0], T, places=10)
        self.assertAlmostEqual(out["T"][0], md, places=10)       # the statistic is the mean delta
        self.assertAlmostEqual(out["SE2"][0], se2, places=10)
        self.assertEqual((out["N"][0], out["K"][0]), (N, K))
        self.assertEqual(out["comps"]["df_mmrd"], N - K)          # 10 - 4 here; 50 - 16 = 34 held out
        lo, hi, df_ = M.interval(out["T"][0], out["SE2"][0], out["comps"])
        self.assertAlmostEqual(df_, df, places=8)
        from scipy.stats import t as tdist
        self.assertAlmostEqual(hi, T + tdist.ppf(0.975, df) * np.sqrt(se2), places=10)
        self.assertAlmostEqual(lo, T - tdist.ppf(0.975, df) * np.sqrt(se2), places=10)

    def test_class_stat_gene_weights(self):
        As, ys = synth(seed=6, n_exon=8)
        genes = np.array(["a", "a", "a", "b", "c", "c", "d", "e"])
        fn = M.gene_weights(genes)
        out = M.class_stat(As, [y[:, None] for y in ys], fn)
        T, se2, *_ = brute_class(As, ys, weights=fn)
        self.assertAlmostEqual(out["T"][0], T, places=10)
        self.assertAlmostEqual(out["SE2"][0], se2, places=10)
        ex = [brute_exon(As, ys, e) for e in range(8)]
        g = pd.DataFrame({"g": genes, "d": [d for d, _, _ in ex]}).dropna().groupby("g").d.mean()
        self.assertAlmostEqual(out["T"][0], g.mean(), places=10)

    def test_permuted_columns_equal_a_run_on_the_permuted_labels(self):
        As, ys = synth(seed=7)
        g = np.concatenate([np.full(len(y), k) for k, y in enumerate(ys)])
        Y = M.label_perms(np.concatenate(ys), g, 5, seed=0)
        Ys = [Y[g == k] for k in range(len(ys))]
        out = M.class_stat(As, Ys)
        for k in (1, 4):
            ref = M.class_stat(As, [y[:, k:k + 1] for y in Ys])
            self.assertAlmostEqual(out["z"][k], ref["z"][0], places=10)
        D, T, N = M.exon_engine(As, Ys)
        D3, T3, _ = M.exon_engine(As, [y[:, 3:4] for y in Ys])
        np.testing.assert_allclose(D[3], D3[0], rtol=1e-6)

    def test_arm_stat_against_brute_force(self):
        rng = np.random.default_rng(8)
        g = np.repeat([0, 1, 2, 3], [40, 30, 25, 20])
        v = rng.normal(0, 1, len(g))
        v[rng.random(len(g)) < 0.15] = np.nan
        y = np.zeros(len(g), bool)
        for c, k in zip(range(4), [5, 3, 0, 2]):
            y[rng.choice(np.flatnonzero(g == c), k, replace=False)] = True
        out = M.arm_stat(v, y[:, None], g)
        ok = np.isfinite(v)
        valid = [c for c in range(4) if (ok & y & (g == c)).any() and (ok & ~y & (g == c)).any()]
        N = sum(int((ok & y & (g == c)).sum()) for c in valid)
        T, ssd, vp = 0.0, 0.0, 0.0
        for c in valid:
            d, p = v[ok & y & (g == c)], v[ok & ~y & (g == c)]
            T += (d - p.mean()).sum() / N
            u = d / N
            ssd += np.square(u - u.mean()).sum()
            vv = p * len(d) / (N * len(p))
            vp += len(p) * vv.var(ddof=1)
        se2 = N * ssd / (N - len(valid)) + vp
        self.assertAlmostEqual(out["T"][0], T, places=12)
        self.assertAlmostEqual(out["SE2"][0], se2, places=12)
        self.assertEqual(out["N"][0], N)


class TestPermutationSchemes(unittest.TestCase):
    def test_label_perms_keep_counts_within_group_and_are_seeded(self):
        rng = np.random.default_rng(0)
        g = rng.integers(0, 5, 200)
        y = rng.random(200) < 0.2
        A = M.label_perms(y, g, 50, seed=0)
        B = M.label_perms(y, g, 50, seed=0)
        C = M.label_perms(y, g, 20, seed=0)
        np.testing.assert_array_equal(A, B)
        np.testing.assert_array_equal(A[:, :21], C)            # the count model's 200 = the first 200 of 2,000
        np.testing.assert_array_equal(A[:, 0], y)
        for k in range(5):
            np.testing.assert_array_equal(A[g == k].sum(0), y[g == k].sum())
        self.assertFalse(np.array_equal(A[:, 1], y))

    def test_constants_match_the_block(self):
        inf = {k: " ".join(str(v).split()) for k, v in PRE["inference"].items()}
        self.assertIn("2,000 permutations, seed 0", inf["label_perm"])
        self.assertEqual(M.N_PERM, 2000)
        self.assertIn("10,000 permutations, seed 0", inf["class_perm"])
        self.assertEqual(M.N_CLASS_PERM, 10_000)
        self.assertIn("200 permutations, seed 0", " ".join(" ".join(PRE["sensitivity"]).split()))
        self.assertEqual(M.N_COUNT_PERM, 200)
        self.assertEqual((M.SEED, M.DRY_SEED), (0, 999))
        self.assertTrue(PRE["sesoi"].startswith(f"{M.SESOI} PSI"))
        self.assertEqual((M.MIN_CLASS, M.MIN_TUMOUR, M.MIN_RP3D, M.MIN_RP3E, M.MIN_SHAP), (20, 10, 15, 40, 10))
        self.assertIn("1,000 draws (seed 0)", " ".join(PRE["instrument_check"]["design"].split()))
        self.assertEqual(M.N_DRAWS, 1000)

    def test_perm_p_and_iut(self):
        self.assertAlmostEqual(M.perm_p(2.0, np.array([1, -3, 2, 0.5])), 3 / 5)
        self.assertAlmostEqual(M.perm_p(-2.0, np.array([1, -3, np.nan])), 2 / 3)
        self.assertEqual(M.iut(0.01, 0.2), 0.2)
        self.assertEqual(M.iut(0.3, np.nan), 0.3)
        self.assertTrue(np.isnan(M.iut(np.nan, None)))

    def test_bh_matches_statsmodels(self):
        from statsmodels.stats.multitest import multipletests
        p = np.array([0.001, 0.04, 0.03, 0.2, np.nan, 0.012])
        q = M.bh(p)
        ok = np.isfinite(p)
        np.testing.assert_allclose(q[ok], multipletests(p[ok], method="fdr_bh")[1])
        self.assertTrue(np.isnan(q[~ok]).all())

    def test_class_perm_within_coarsened_strata(self):
        rng = np.random.default_rng(2)
        lab = np.array(["0-0:0"] * 12 + ["1-2"] * 9 + ["3-3:1"] * 7)
        first = np.zeros(len(lab), bool)
        first[[0, 1, 2, 12, 13, 21]] = True
        v = rng.normal(0, 1, len(lab))
        v[first] += 0.5
        obs, null = M.class_perm(v, first, lab, 4000, seed=0)
        self.assertAlmostEqual(obs, M.mean_by_strata(v[None, :], first, ~first, lab)[0], places=12)
        w = {s: first[lab == s].sum() for s in set(lab)}
        ref = sum(w[s] * (v[(lab == s) & first].mean() - v[(lab == s) & ~first].mean()) for s in w) / sum(w.values())
        self.assertAlmostEqual(obs, ref, places=12)
        # permutations stay within strata: shifting every exon of one stratum leaves the observed
        # statistic and its null unchanged (a permutation across strata would carry the shift)
        self.assertLess(abs(null.mean()), 0.05)
        v2 = v.copy()
        v2[lab == "3-3:1"] += 100
        obs2, null2 = M.class_perm(v2, first, lab, 4000, seed=0)
        self.assertAlmostEqual(obs2, obs, places=8)
        np.testing.assert_allclose(null2, null, atol=1e-8)

    def test_class_perm_exhaustive_single_stratum(self):
        from itertools import combinations
        v = np.array([0.3, 1.2, -0.4, 0.8, 0.1])
        first = np.array([1, 1, 0, 0, 0], bool)
        lab = np.array(["s"] * 5)
        obs, null = M.class_perm(v, first, lab, 20000, seed=1)
        exact = [v[list(c)].mean() - np.delete(v, list(c)).mean() for c in combinations(range(5), 2)]
        p_exact = np.mean(np.abs(exact) >= abs(obs) - 1e-12)
        self.assertAlmostEqual(M.perm_p(obs, null), p_exact, delta=0.01)

    def test_partial_spearman_and_in_cell_permutation(self):
        from scipy.stats import rankdata
        rng = np.random.default_rng(4)
        n = 60
        L, o = rng.integers(20, 40, n).astype(float), rng.integers(3, 15, n).astype(float)
        g = rng.normal(0, 1, n) + 0.1 * L
        x = 0.5 * rankdata(g) + rng.normal(0, 5, n)
        Z = np.column_stack([rankdata(L), rankdata(o)])
        X = np.column_stack([np.ones(n), Z])
        rx = rankdata(x) - X @ np.linalg.lstsq(X, rankdata(x), rcond=None)[0]
        rg = rankdata(g) - X @ np.linalg.lstsq(X, rankdata(g), rcond=None)[0]
        self.assertAlmostEqual(M.partial_spearman(x, g, Z), rx @ rg / np.sqrt(rx @ rx * rg @ rg), places=10)
        cells = np.array([f"{a}:{b}" for a, b in zip(np.digitize(L, [20, 25]), np.digitize(o, [3, 5]))])
        idx = M.within_perm_index(cells, 200, seed=0)
        for k in range(200):
            self.assertTrue((cells[idx[k]] == cells).all())
            self.assertEqual(sorted(idx[k]), list(range(n)))
        obs, null = M.partial_spearman_perm(x, g, Z, cells, 200, seed=0)
        for k in (0, 7, 199):
            self.assertAlmostEqual(null[k], M.partial_spearman(x, g[idx[k]], Z), places=10)

    def test_freedman_lane_against_statsmodels_and_within_cells(self):
        import statsmodels.api as sm
        rng = np.random.default_rng(5)
        n = 120
        cells = np.repeat(["a", "b", "c"], 40)
        X = np.column_stack([np.ones(n), rng.normal(0, 1, n), (cells == "b").astype(float), rng.random(n) < 0.3])
        y = X[:, 1] * 0.5 + rng.standard_t(4, n)
        b, t = M.ols_hc3(X, y, 3)
        f = sm.OLS(y, X).fit(cov_type="HC3")
        self.assertAlmostEqual(b, f.params[3], places=10)
        self.assertAlmostEqual(t, f.tvalues[3], places=8)
        b2, t2, null = M.freedman_lane(X, y, 3, cells, 300, seed=0)
        self.assertEqual((b2, t2), (b, t))
        idx = M.within_perm_index(cells, 300, seed=0)
        Xr = np.delete(X, 3, axis=1)
        fit = Xr @ np.linalg.lstsq(Xr, y, rcond=None)[0]
        for k in (0, 123):
            ys = fit + (y - fit)[idx[k]]
            self.assertAlmostEqual(null[k], M.ols_hc3(X, ys, 3)[1], places=8)
        self.assertTrue(all((cells[idx[k]] == cells).all() for k in range(300)))

    def test_freedman_lane_is_calibrated_under_the_null(self):
        rng = np.random.default_rng(6)
        ps = []
        for r in range(60):
            n = 80
            cells = np.repeat(np.arange(4), 20)
            X = np.column_stack([np.ones(n), rng.normal(0, 1, n), rng.random(n) < 0.25])
            y = X[:, 1] + rng.normal(0, 1, n)
            _, t, null = M.freedman_lane(X, y, 2, cells, 199, seed=r)
            ps.append(M.perm_p(t, null))
        self.assertLess(np.mean(np.array(ps) < 0.05), 0.15)
        self.assertGreater(np.mean(ps), 0.35)


class TestDesigns(unittest.TestCase):
    def test_band_adjust_neighbour_rule(self):
        rng = np.random.default_rng(1)
        E = 60
        band = rng.integers(0, 13, E)
        band[:3] = -1                                   # outside the universe
        bmask = rng.random(E) < 0.7
        dev = rng.normal(0, 0.1, (5, E)).astype(np.float32)
        dev[rng.random((5, E)) < 0.2] = np.nan
        dev[:, band < 0] = np.nan
        adj, ga = M.band_adjust(dev, band, bmask & (band >= 0), min_n=4)
        for t in range(5):
            for e in range(3, E):
                if not np.isfinite(dev[t, e]):
                    continue
                pool = bmask & (band >= 0) & np.isfinite(dev[t])
                same = pool & (band == band[e])
                nb = pool & (np.abs(band - band[e]) <= 1)
                if same.sum() >= 4:
                    ref = dev[t, same].mean()
                elif nb.sum() >= 4:
                    ref = dev[t, nb].mean()
                else:
                    self.assertTrue(np.isnan(adj[t, e]))
                    continue
                self.assertAlmostEqual(float(adj[t, e]), float(dev[t, e] - ref), places=5)
            self.assertAlmostEqual(ga[t], float(np.nanmean(np.where(bmask & (band >= 0), dev[t], np.nan))), places=5)

    def test_mirrored_instrument_draws(self):
        rng = np.random.default_rng(3)
        n = 400
        up = np.zeros(n, bool)
        up[:30] = True
        dn = ~up & (rng.random(n) < 0.8)
        base = np.where(up, rng.uniform(0, 0.15, n), rng.uniform(0.8, 1.0, n))
        terc = rng.integers(0, 3, n)
        band = M.band_of(base)
        band_m = band.copy()
        band_m[up] = M.band_of(1 - base[up])            # the mirrored baseline of each R_up exon
        self.assertTrue((band_m[up] >= 10).all())
        lab, kept = M.coarsen(band_m, terc, up, dn, per_first=1)
        draws = M.coarsened_draws(lab, up & kept, dn & kept, 1, 200, np.random.default_rng(0))[:, :, 0]
        fi = np.flatnonzero(up & kept)
        for d in range(200):
            got = draws[d][draws[d] >= 0]
            self.assertEqual(len(got), len(set(got)))                # without replacement within a draw
            self.assertTrue(dn[got].all())
            for k, e in enumerate(fi):
                if draws[d, k] >= 0:
                    self.assertEqual(lab[draws[d, k]], lab[e])       # same mirrored band x tertile stratum
        self.assertFalse(np.array_equal(draws[0], draws[1]))
        again = M.coarsened_draws(lab, up & kept, dn & kept, 1, 200, np.random.default_rng(0))[:, :, 0]
        np.testing.assert_array_equal(draws, again)

    def test_instrument_draw_statistic_is_rp2a_on_the_drawn_exons(self):
        As, ys = synth(seed=21, n_exon=14)
        draws = np.array([[0, 3, 5, 7, 9], [1, 2, 4, 6, 8], [10, 11, 12, 13, -1]])
        cols = np.arange(14)
        Ys = [np.repeat(y[:, None], 3, axis=1) for y in ys]
        out = M.class_stat(As, Ys, M.draw_weights(draws, cols), block=2)     # blocks smaller than the draws
        for d in range(3):
            sel = draws[d][draws[d] >= 0]
            ref = M.class_stat([A[:, sel] for A in As], [y[:, None] for y in ys])
            self.assertAlmostEqual(out["T"][d], ref["T"][0], places=12)
            self.assertAlmostEqual(out["z"][d], ref["z"][0], places=10)

    def test_rp2c_cut_from_b_sample_permutations(self):
        rng = np.random.default_rng(22)
        E, P = 120, 301
        tab = pd.DataFrame({"cls": ["R_up"] * 10 + ["B_far"] * 100 + ["R_down"] * 10,
                            "b_sample": [False] * 10 + [True] * 100 + [False] * 10})
        Tm = rng.normal(0, 1, (P, E)).astype(np.float32)
        Tm[:, 115:] = 50.0                                # non-B exons never set the cut
        Tm[0, :10] = 6.0
        elig = np.ones(E, bool)
        elig[105:110] = False                             # ineligible B exons never set the cut
        Tm[:, 105:110] = 40.0
        band, terc = np.zeros(E, int), np.zeros(E, int)
        ctx = M.XCtx(M.SPECS["primary"], "held_out", np.arange(E), tab, Tm, Tm, None, [], elig, band, terc,
                     np.zeros(E), np.ones(E))
        r = M.rp2c(ctx)
        bs = tab.b_sample.to_numpy() & elig
        self.assertAlmostEqual(r["cut"], float(np.percentile(Tm[1:, bs].astype(float).ravel(), 97.5)), places=5)
        self.assertEqual(r["n_sets"], 10)
        self.assertEqual(r["share_up"], 1.0)
        self.assertLess(r["p_class"], 0.01)
        self.assertLess(r["p_label"], 0.02)

    def test_contrast_label_and_class_permutation(self):
        rng = np.random.default_rng(23)
        E, P = 200, 401
        first = np.zeros(E, bool)
        first[:40] = True
        ctrl = ~first
        band = rng.integers(0, 4, E)
        terc = rng.integers(0, 3, E)
        D = rng.normal(0, 0.01, (P, E))
        D[0, :40] += 0.01
        D[5, 3] = np.nan                                   # an exon losing its delta in one column
        ctx = M.XCtx(M.SPECS["primary"], "held_out", np.arange(E), pd.DataFrame(index=range(E)), D, D, None, [],
                     np.ones(E, bool), band, terc, np.zeros(E), np.ones(E))
        r = M.contrast(ctx, first, ctrl, min_n=20)
        lab, kept = M.coarsen(band, terc, first, ctrl)
        use = kept
        stat = M.mean_by_strata(D[:, use], first[use], ctrl[use], lab[use])
        self.assertAlmostEqual(r["estimate"], stat[0], places=12)
        self.assertAlmostEqual(r["p_label"], M.perm_p(stat[0], stat[1:]), places=12)
        self.assertLess(r["p_class"], 0.01)
        self.assertTrue(r["testable"])
        self.assertFalse(M.contrast(ctx, first & (np.arange(E) < 15), ctrl, min_n=20)["testable"])

    def test_rp6c_held_out_draws(self):
        rng = np.random.default_rng(0)
        project = np.repeat(["A", "B", "C"], [40, 30, 30])
        mmr = np.zeros(100, bool)
        mmr[[0, 1, 2, 3, 4, 5, 40, 41, 42, 43, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82]] = True
        prof = {p: np.flatnonzero((project == p) & ~mmr) for p in "ABC"}
        ia, ib = M.rp6c_draws(np.flatnonzero(mmr), project, prof, rng)
        self.assertEqual(len(ia), 50)
        for a, b in zip(ia, ib):
            self.assertEqual(len(set(a)), 20)
            self.assertTrue(mmr[a].all() and not mmr[b].any())
            self.assertTrue((project[a] == project[b]).all())         # each with an MMRp tumour of its cancer
            self.assertEqual(len(set(b)), 20)                         # no MMRp tumour twice within a draw
        ia2, ib2 = M.rp6c_draws(np.flatnonzero(mmr), project, prof, np.random.default_rng(0))
        self.assertTrue(all((x == y).all() for x, y in zip(ib, ib2)))

    def test_instrument_verdict_rule(self):
        null = np.random.default_rng(0).normal(0, 1, 2000)
        crit = np.percentile(null, 2.5)
        v = M.instrument_verdict(np.full(100, crit - 0.5), null)
        self.assertTrue(v["passes"])
        self.assertAlmostEqual(v["critical_value"], crit)
        z = np.concatenate([np.full(79, crit - 1), np.full(21, crit + 1)])
        self.assertFalse(M.instrument_verdict(z, null)["passes"])            # 79% below: fails
        z = np.concatenate([np.full(80, crit - 1), np.full(20, 5.0)])
        self.assertTrue(M.instrument_verdict(z, null)["passes"])

    def test_rp2c_cut_and_rise_shares(self):
        rng = np.random.default_rng(9)
        P, E = 201, 40
        T = rng.normal(0, 1, (P, E))
        T[0, :5] = 4.0                                  # first exons rise in the observed column
        bs = np.zeros(E, bool)
        bs[5:] = True
        cut = float(np.nanpercentile(T[1:, bs], 97.5))
        self.assertAlmostEqual(cut, np.percentile(T[1:, 5:].ravel(), 97.5))
        sets = np.column_stack([np.arange(5), np.arange(5, 30).reshape(5, 5)])
        rise = T > cut
        st = M.rise_share(rise, sets)
        self.assertAlmostEqual(st[0], 1.0 - rise[0, 5:30].mean())
        null = M.rise_class_perm(rise[0], sets, 5000, seed=0)
        self.assertLess(M.perm_p(st[0], null), 0.01)
        self.assertLess(abs(null.mean()), 0.05)

    def test_host_gene_change_recomputed_per_label_column(self):
        rng = np.random.default_rng(11)
        Gs = [rng.normal(5, 1, (n, 4)) for n in (20, 15, 30)]
        ys = [np.isin(np.arange(n), rng.choice(n, k, replace=False)) for n, k in ((20, 3), (15, 0), (30, 4))]
        g = np.repeat(np.arange(3), [20, 15, 30])
        Y = M.label_perms(np.concatenate(ys), g, 6, seed=0)
        Ys = [Y[g == k] for k in range(3)]
        Us = [np.array([True, True, False]), np.array([True, False, True]), np.array([False, True, True])]
        gi = np.array([0, 2, 3])
        hc = M.gene_change(Gs, Ys, Us, gi)
        for k in (0, 5):
            for e in range(3):
                num, den = 0.0, 0
                for G, Yc, U in zip(Gs, Ys, Us):
                    y = Yc[:, k]
                    if U[e] and y.any():
                        num += (G[y, gi[e]] - G[~y, gi[e]].mean()).sum()
                        den += y.sum()
                self.assertAlmostEqual(hc[k, e], num / den, places=10)
        self.assertFalse(np.allclose(hc[0], hc[5]))

    def test_matched_draws_same_band_and_tertile(self):
        rng = np.random.default_rng(12)
        fb, ft = rng.integers(0, 4, 15), rng.integers(0, 3, 15)
        pb, pt = rng.integers(0, 4, 300), rng.integers(0, 3, 300)
        d = M.matched_draws(fb, ft, pb, pt, 50, np.random.default_rng(0))
        self.assertTrue((pb[d] == fb[None, :]).all())
        self.assertTrue((pt[d] == ft[None, :]).all())
        for k in range(50):
            for key in set(zip(fb, ft)):
                j = [i for i in range(15) if (fb[i], ft[i]) == key]
                self.assertEqual(len(set(d[k, j])), len(j))            # without replacement within a set

    def test_auroc_matrix(self):
        from sklearn.metrics import roc_auc_score
        rng = np.random.default_rng(13)
        S = rng.integers(0, 6, (80, 3)).astype(float)
        S[rng.random((80, 3)) < 0.1] = np.nan
        Y = rng.random((80, 2)) < 0.3
        A = M.auroc_matrix(S, Y)
        for i in range(3):
            for j in range(2):
                f = np.isfinite(S[:, i])
                self.assertAlmostEqual(A[i, j], roc_auc_score(Y[f, j], S[f, i]), places=10)

    def test_quasi_binomial_against_statsmodels(self):
        import statsmodels.api as sm
        rng = np.random.default_rng(14)
        n = 300
        c = rng.integers(0, 3, n)
        X = np.column_stack([np.eye(3)[c], rng.random(n), rng.random(n) < 0.2])
        tot = rng.integers(5, 60, n)
        off = np.log(rng.uniform(0.5, 2, n))
        eta = X @ np.array([-1, 0, 0.5, 0.3, 0.4]) + off
        y = rng.binomial(tot, 1 / (1 + np.exp(-eta)))
        b, cov, scale, conv = M.quasi_binomial(y, tot, off, X)
        f = sm.GLM(y / tot, X, family=sm.families.Binomial(), var_weights=tot, offset=off).fit(scale="X2")
        self.assertTrue(conv)
        np.testing.assert_allclose(b, f.params, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(np.sqrt(np.diag(cov)), f.bse, rtol=1e-5)
        self.assertAlmostEqual(scale, f.scale, places=6)

    def test_count_model_drops_cancers_without_inclusion_reads(self):
        rng = np.random.default_rng(15)
        n = 240
        cidx = np.repeat([0, 1, 2], 80)
        I = rng.binomial(30, 0.4, (n, 1)).astype(float)
        Sj = (30 - I)
        I[cidx == 2] = 0                                 # cancer 2: no inclusion read at this exon
        R = np.ones((n, 1))
        ga = rng.random(n)
        y = np.zeros(n, bool)
        y[[1, 5, 90, 99, 170]] = True
        Y = M.label_perms(y, cidx, 3, seed=0)
        U = np.ones((3, 1), bool)
        B, SE, _ = M._glm_chunk(I, Sj, R, U, cidx, ga, Y, 3)
        keep = cidx < 2
        X = np.column_stack([np.eye(2)[cidx[keep]], ga[keep], Y[keep, 0]])
        b, cov, _, _ = M.quasi_binomial(I[keep, 0], (I + Sj)[keep, 0], np.zeros(keep.sum()), X)
        self.assertAlmostEqual(B[0, 0], b[-1], places=6)
        self.assertAlmostEqual(SE[0, 0], np.sqrt(cov[-1, -1]), places=6)

    def test_inverse_variance_class_means(self):
        beta = np.array([[0.2, 0.1], [0.4, -0.1], [0.0, 0.3], [1.0, 0.0]])
        se = np.array([[0.1, 0.1], [0.2, 0.2], [0.1, 0.1], [0.5, 0.5]])
        first = np.array([True, True, False, False])
        m = M.ivw_contrast(beta, se, first, None, None)
        w = 1 / se[:2, 0] ** 2
        self.assertAlmostEqual(m[0], (w * beta[:2, 0]).sum() / w.sum())
        lab = np.array(["s", "s", "s", "s"])
        d = M.ivw_contrast(beta, se, first, ~first, lab)
        wc = 1 / se[2:, 0] ** 2
        self.assertAlmostEqual(d[0], m[0] - (wc * beta[2:, 0]).sum() / wc.sum())
        null = M.ivw_class_perm(beta[:, 0], se[:, 0], first, ~first, lab, 500, seed=0)
        self.assertTrue(np.isfinite(null).all())

    def test_rp6c_null_contrasts_and_draws(self):
        rng = np.random.default_rng(16)
        tab = pd.DataFrame({"file_id": [f"f{i}" for i in range(130)],
                            "label": ["A"] * 40 + ["B"] * 30 + ["C"] * 25 + ["MSI"] * 20 + ["D"] * 15})
        tab["sample15"] = [f"s{i}" for i in range(130)]
        tasks = {"X": tab, "Y": tab.assign(label=["P"] * 60 + ["Q"] * 70)}
        mmr = {"s0", "s1", "s50"}
        d71 = set(tab.file_id) - {"f60", "f61", "f62", "f63", "f64", "f65"}   # not D-71 tumours
        out = M.null_contrasts(tasks, mmr, d71)
        pairs = [(c, a, b) for c, a, b, *_ in out]
        self.assertEqual(pairs, [("X", "A", "B"), ("X", "A", "C"), ("X", "B", "C"), ("Y", "P", "Q")])
        fa = dict(((c, a, b), x) for c, a, b, x, _ in out)
        self.assertNotIn("f0", fa[("X", "A", "B")])      # MMRd removed
        self.assertEqual(M.RP6C_NULL, 57)
        self.assertIn("the 57 pairs of classes", " ".join(next(t["claim"] for t in PRE["tests"] if t["id"] == "RP6c").split()))
        A = rng.normal(0, 1, (100, 6))
        A[rng.random((100, 6)) < 0.2] = np.nan
        ia = [rng.choice(50, 20, replace=False) for _ in range(5)]
        ib = [50 + rng.choice(50, 20, replace=False) for _ in range(5)]
        ref = []
        for a, b in zip(ia, ib):
            ok = (np.isfinite(A[a]).sum(0) >= 10) & (np.isfinite(A[b]).sum(0) >= 10)
            ref.append(np.mean(np.nanmean(A[a][:, ok], 0) - np.nanmean(A[b][:, ok], 0)))
        self.assertAlmostEqual(M.contrast_draw_stat(A, ia, ib), np.mean(ref), places=12)

    def test_rp7_pair_rule_and_statistic(self):
        s = pd.DataFrame({"file_id": ["t1", "t2", "t3", "n2b", "n2a", "n3", "n9"],
                          "case_id": ["c1", "c2", "c3", "c2", "c2", "c3", "c9"],
                          "role": ["tumour"] * 3 + ["normal"] * 4})
        p = M.pair_rows(s)
        self.assertEqual(dict(zip(p.case_id, p.normal_file_id)), {"c2": "n2a", "c3": "n3"})   # first normal in sort order
        self.assertEqual(dict(zip(p.case_id, p.normal_row)), {"c2": 4, "c3": 5})
        arm = np.array([0.1, 0.3, -0.2, 0.0, 0.05, np.nan, 0.2, 0.4])
        g = np.array(["A"] * 5 + ["B"] * 3)
        y = np.array([1, 0, 0, 0, 0, 1, 0, 0], bool)
        st = M.pair_stat(arm, y[:, None], g)
        self.assertAlmostEqual(st[0], 0.1 - np.mean([0.3, -0.2, 0.0, 0.05]))   # B's MMRd pair has no arm
        Y = M.label_perms(y, g, 10, seed=0)
        self.assertTrue((Y[g == "A"].sum(0) == 1).all())

    def test_shap_direction_from_a_toy_array(self):
        from scipy.stats import spearmanr
        rng = np.random.default_rng(17)
        X = rng.normal(0, 1, (30, 4))
        S = np.column_stack([X[:, 0], -X[:, 1] ** 3, rng.normal(0, 1, 30), np.zeros(30)])
        rows = np.arange(0, 30, 2)
        d = M.shap_directions(X, S, rows)
        self.assertAlmostEqual(d[0], 1.0)
        self.assertAlmostEqual(d[1], -1.0)
        self.assertAlmostEqual(d[2], spearmanr(X[rows, 2], S[rows, 2]).statistic, places=10)
        self.assertTrue(np.isnan(d[3]))                  # no spread in SHAP
        jd = {"classes": ["CIN", "MSI", "GS"]}
        task = pd.DataFrame({"label": ["CIN", "MSI", "GS", "MSI"], "y": [0, 1, 2, 1]})
        self.assertEqual(M.msi_class_index(jd, task), 1)
        with self.assertRaises(AssertionError):
            M.msi_class_index({"classes": ["MSI", "CIN"]}, task)

    def test_rpd1_permutation_null_is_calibrated_where_binomial_is_not(self):
        """A-22: with correlated inputs (one latent factor drives every input's value, another every
        input's SHAP, independently), the binomial test on the share of positive directions rejects far
        too often; the within-fold row permutation does not."""
        from scipy.stats import binomtest
        rej_b, rej_p = [], []
        for r in range(30):
            rng = np.random.default_rng(100 + r)
            folds = []
            for c in range(3):
                for k in range(3):
                    n, p_ = 10, 40
                    z, w = rng.normal(0, 1, (n, 1)), rng.normal(0, 1, (n, 1))
                    X = z + 0.5 * rng.normal(0, 1, (n, p_))
                    Sv = w + 0.5 * rng.normal(0, 1, (n, p_))
                    folds.append({"cohort": f"c{c}", "X": X, "S": Sv, "fids": np.array([f"e{j}" for j in range(p_)]),
                                  "arm": np.array(["R_down"] * p_)})
            out = M.rpd1_test(folds, n_perm=99, seed=r)
            kk, nn = M.direction_shares(folds)[0]["R_down"]
            rej_b.append(binomtest(kk, nn, 0.5).pvalue < 0.05)
            rej_p.append(out["R_down"]["p_class"] < 0.05)
        self.assertGreater(np.mean(rej_b), 0.3)
        self.assertLess(np.mean(rej_p), 0.15)

    def test_rpd1_shares_match_fold_averaged_directions(self):
        rng = np.random.default_rng(3)
        folds = [{"cohort": "A", "X": rng.normal(0, 1, (8, 3)), "S": rng.normal(0, 1, (8, 3)),
                  "fids": np.array(["e1", "e2", "e3"]), "arm": np.array(["R_up", "R_up", "R_down"])},
                 {"cohort": "A", "X": rng.normal(0, 1, (9, 2)), "S": rng.normal(0, 1, (9, 2)),
                  "fids": np.array(["e1", "e3"]), "arm": np.array(["R_up", "R_down"])}]
        sh, t = M.direction_shares(folds)
        d = [M.shap_directions(f["X"], f["S"], np.arange(len(f["X"]))) for f in folds]
        e1 = (d[0][0] + d[1][0]) / 2
        self.assertAlmostEqual(float(t.set_index("feature_id").direction["e1"]), e1, places=12)
        self.assertEqual(sh["R_up"], (int(e1 > 0) + int(d[0][1] > 0), 2))
        self.assertEqual(sh["R_down"][1], 1)

    def test_separated_fits_dropped(self):
        B = np.array([[0.2, 12.0], [-11.0, 0.1], [0.3, 0.4]])
        SE = np.array([[0.1, 0.2], [0.5, 20.0], [15.0, 0.3]])
        np.testing.assert_array_equal(M.separated(B, SE), [[False, True], [True, True], [True, False]])
        self.assertEqual(M.SEP_LIMIT, 10.0)

    def test_committed_descriptives_blank_cells_under_five_mmrd(self):
        d = pd.DataFrame({"item": ["per-cancer RP2b"] * 3 + ["other"], "group": ["ACC", "PRAD", "BRCA", ""],
                          "value": [0.1, 0.2, 0.3, 0.4], "n": ["3 MMRd", "4 MMRd", "11 MMRd", "5"],
                          "note": [np.nan, "x", "", "keep"], "p": [0.5, 0.6, 0.7, np.nan],
                          "ci_low": [np.nan] * 4, "ci_high": [np.nan] * 4, "n_mmrd": [3, 4, 11, np.nan]})
        b = M.blank_small(d)
        self.assertTrue(b.value[:2].isna().all() and b.p[:2].isna().all())
        self.assertEqual(list(b.value[2:]), [0.3, 0.4])
        self.assertEqual(list(b.note), ["blanked (< 5 MMRd tumours)", "x; blanked (< 5 MMRd tumours)", "", "keep"])
        self.assertEqual(list(b.n), list(d.n))                 # rows and their n stay
        self.assertEqual(d.value[0], 0.1)                      # the raw table is not changed

    def test_matched_draws_log_the_rule_used(self):
        fb, ft = np.array([0, 0, 1, 2, 4, 4, 5]), np.array([0, 0, 1, 0, 0, 0, 1])
        pb, pt = np.array([0, 0, 0, 1, 1, 3, 4, 5]), np.array([0, 1, 1, 0, 2, 0, 1, 1])
        log = {}
        d = M.matched_draws(fb, ft, pb, pt, 20, np.random.default_rng(0), log)
        self.assertEqual(log, {"band x tertile": 1, "tertiles merged": 3, "neighbouring bands": 1,
                               "tertiles merged, with replacement": 2})
        self.assertTrue((pb[d[:, 3]] != 2).all())
        self.assertTrue((d[:, 4] == d[:, 5]).all())            # one candidate, drawn with replacement

    def test_m4_guard(self):
        M.m4_guard(0.9004, real=True)
        M.m4_guard(0.8996, real=True)
        with self.assertRaises(SystemExit):
            M.m4_guard(0.8985, real=True)
        with self.assertRaises(SystemExit):
            M.m4_guard(0.95, real=True)
        M.m4_guard(0.5, real=False)                       # the dry run reports, never enforces


class TestRunGuards(unittest.TestCase):
    def test_real_labels_need_approval_and_paths_are_separated(self):
        with self.assertRaises(SystemExit):
            M.Run("real")
        with self.assertRaises(SystemExit):
            M.Run("observed")
        dry = M.Run("shuffled")
        self.assertIn("dryrun", dry.out.parts)
        self.assertIn("dryrun", dry.local.parts)
        for bad in [M.OUT / "tests.csv", M.PROOF / "tests" / "x.npz", M.PROOF / "x.npy"]:
            self.assertFalse(dry.allowed(bad), bad)
        self.assertTrue(dry.allowed(M.OUT / "dryrun" / "tests.csv"))
        with self.assertRaises(RuntimeError):
            dry.path("../tests.csv")
        real = M.Run("real", approved=True)
        self.assertFalse(real.allowed(M.OUT / "dryrun" / "tests.csv"))
        self.assertTrue(real.allowed(M.OUT / "tests.csv"))
        self.assertTrue(str(dry.local).startswith(str(M.m51.DATA)))     # per-sample outputs: data.nosync
        self.assertTrue(str(real.local).startswith(str(M.m51.DATA)))

    def test_per_sample_outputs_are_written_locally(self):
        src = (REPO / "scripts" / "62_msi_rising_proof.py").read_text()
        calls = re.findall(r'run\.path\(\s*(f?"[^"]+")([^)]*)\)', src)
        self.assertTrue(calls)
        for name, rest in calls:
            if re.search(r"\.(parquet|npz|npy|pkl|json)\"$", name) and name.strip('f"') not in (
                    "instrument.json",) and "timing" not in name and "verdict" not in name:
                self.assertIn("local=True", rest, name)
        for name in ["tumours.parquet", "pairs.parquet", "rpd1_directions.parquet", "count_model.npz"]:
            self.assertRegex(src, rf'run\.path\("{re.escape(name)}", local=True\)')

    def test_shuffle_keeps_counts_within_cancer_and_pairs(self):
        rng = np.random.default_rng(0)
        t = pd.DataFrame({"project": np.repeat(["A", "B", "C"], [50, 40, 30]),
                          "has_pair": rng.random(120) < 0.2})
        t["mmr"] = rng.random(120) < 0.1
        t["hyper"] = ~t.mmr & (rng.random(120) < 0.1)
        t["n_hp"] = np.where(t.mmr, 20, 0)
        t["n_snv"] = rng.integers(0, 900, 120)
        for c in ["mantis", "msisensor"]:
            t[c] = rng.random(120)
        for c in ["call_mantis", "call_msisensor"]:
            t[c] = np.where(rng.random(120) < 0.1, "MSI-H", "MSS")
        s = M.shuffle_labels(t, 999)
        pd.testing.assert_frame_equal(s, M.shuffle_labels(t, 999))
        self.assertFalse(s.mmr.equals(t.mmr))
        for keys, g in t.groupby(["project", "has_pair"]):
            h = s.loc[g.index]
            self.assertEqual(int(h.mmr.sum()), int(g.mmr.sum()))
            self.assertEqual(sorted(h.n_snv), sorted(g.n_snv))
        self.assertTrue(((s.n_hp == 20) == s.mmr).all())               # label columns move together
        self.assertTrue((~(s.hyper & s.mmr)).all())


class TestVerdict(unittest.TestCase):
    """msi_rising_proof.verdict as code, on hand-made tables."""

    RP = ["RP2a", "RP3a", "RP3e", "RP4a", "RP4b", "RP4c", "RP4d", "RP4e", "RP5a", "RP5b", "RP6a", "RP6c"]

    def table(self, **over):
        rows = {i: ("supported", 0.01, 0.02) for i in self.RP}
        rows.update({"RPC2a": ("supported", -0.01, np.nan), "RP3c": ("supported", 0.2, np.nan),
                     "RP7": ("not resolved", 0.004, np.nan)})
        rows["RP2a"] = ("supported", 0.006, 0.009)
        rows.update(over)
        return pd.DataFrame([{"id": k, "verdict": v, "estimate": e, "ci_high": u} for k, (v, e, u) in rows.items()])

    PASS, FAIL = {"passes": True}, {"passes": False}

    def test_established(self):
        v = M.claim_verdict(self.table(), self.PASS)
        self.assertEqual(v["verdict"], "established")
        self.assertEqual(v["branch_point_test"], "RP3e")

    def test_uninterpretable_takes_precedence(self):
        for rpc in ["not resolved", "against", "not testable"]:
            v = M.claim_verdict(self.table(RPC2a=(rpc, 0.0, np.nan), RP2a=("against", -0.01, 0.0)), self.PASS)
            self.assertEqual(v["verdict"], "uninterpretable")

    def test_refuted_clause_0_needs_no_instrument(self):
        v = M.claim_verdict(self.table(RP2a=("against", -0.004, -0.001)), self.FAIL)
        self.assertEqual(v["verdict"], "refuted")
        self.assertTrue(v["refuted_clauses"]["(0) RP2a against"])

    def test_refuted_clause_i_needs_the_instrument(self):
        t = self.table(RP2a=("not resolved", 0.0005, 0.002))
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "refuted")
        v = M.claim_verdict(t, self.FAIL)
        self.assertEqual(v["verdict"], "partial")
        self.assertIn("instrument cannot see a near-silent falling class", v["reasons"])
        t = self.table(RP2a=("not resolved", -0.001, 0.004))     # estimate <= 0, bound above SESOI
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "refuted")
        t = self.table(RP2a=("not resolved", 0.002, 0.004))      # a rise of interest not excluded
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "partial")

    def test_refuted_clause_ii_and_iii(self):
        self.assertEqual(M.claim_verdict(self.table(RP3a=("not resolved", -0.001, np.nan)), self.PASS)["verdict"], "refuted")
        t = self.table(RP3e=("not resolved", -0.1, np.nan), RP3c=("not resolved", 0.0, np.nan))
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "refuted")
        t = self.table(RP3e=("not resolved", -0.1, np.nan), RP3c=("not resolved", 0.05, np.nan))
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "partial")
        t = self.table(RP3e=("not testable", np.nan, np.nan), RP3c=("not resolved", -0.05, np.nan))
        v = M.claim_verdict(t, self.PASS)
        self.assertEqual(v["verdict"], "refuted")
        self.assertEqual(v["branch_point_test"], "RP3c")
        t = self.table(RP3e=("not testable", np.nan, np.nan), RP3c=("supported", 0.2, np.nan))
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "established")

    def test_refuted_beats_established(self):
        t = self.table(RP3a=("supported", 0.0, np.nan))           # supported needs a sign; estimate 0 refutes
        self.assertEqual(M.claim_verdict(t, self.PASS)["verdict"], "refuted")

    def test_partial_names_unmet_conditions(self):
        v = M.claim_verdict(self.table(RP2a=("supported", 0.002, 0.0025)), self.PASS)
        self.assertEqual(v["verdict"], "partial")
        self.assertIn("rise below the smallest effect of interest", v["reasons"])
        v = M.claim_verdict(self.table(RP7=("against", -0.01, np.nan)), self.PASS)
        self.assertEqual(v["verdict"], "partial")
        self.assertTrue(any(r.startswith("RP7") for r in v["reasons"]))
        v = M.claim_verdict(self.table(RP5b=("not resolved", 0.01, np.nan), RP6c=("not testable", np.nan, np.nan)), self.PASS)
        self.assertEqual(v["verdict"], "partial")
        self.assertTrue({"RP5b not resolved (needed: supported)", "RP6c not testable (needed: supported)"} <= set(v["reasons"]))
        v = M.claim_verdict(self.table(RP7=("not resolved", 0.0, np.nan)), self.PASS)
        self.assertEqual(v["verdict"], "partial")                 # RP7's estimate must be > 0

    def test_test_verdicts(self):
        self.assertEqual(M.test_verdict(0.1, 0.01, 1), "supported")
        self.assertEqual(M.test_verdict(-0.1, 0.01, 1), "against")
        self.assertEqual(M.test_verdict(-0.1, 0.01, -1), "supported")
        self.assertEqual(M.test_verdict(0.1, 0.2, 1), "not resolved")
        self.assertEqual(M.test_verdict(0.1, np.nan, 1), "not resolved")
        self.assertEqual(M.test_verdict(0.1, 0.01, 1, testable=False), "not testable")


if __name__ == "__main__":
    unittest.main()
