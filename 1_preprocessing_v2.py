"""
1_preprocessing_v2.py  —  statistically corrected replacement for the inferential
sections (1.3-1.5, 1.7-1.8) of 1_preprocessing.py.

WHY
---
The original scans test individual MD frames pooled across replicates (n~900-1900) when the
experimental unit is the trajectory (3-4 replicates/system). That pseudoreplication inflates
significance by ~8 orders of magnitude (see prototype_replicate_level_test.py). This module
fixes it.

METHOD (the corrected statistics)
---------------------------------
1. Collapse each trajectory to ONE replicate-level summary over its EQUILIBRATED frames
   (SUMMARY = 'median' [default, robust] or 'mean' [thermodynamic] for energies; MEDIAN always
   for BC). This makes the replicate the unit and kills the pseudoreplication.
2. Compare conditions with WELCH's t-test on those replicate-level summaries (unequal-variance,
   n=3-4 per group), then Benjamini-Hochberg FDR across the tested family. Report the effect
   size (Delta) and its 95% CI.
3. Keep the original |Delta mean| > 2 kJ/mol effect-size gate for the energy scans.
4. BC (bounded [0,1], ~14% zeros, skew ~3.9): replicate summary = MEDIAN; residue-level omnibus
   = Kruskal-Wallis on replicate medians; per-contrast dBC = Welch on replicate medians. A
   non-significant dBC cell is set to NaN so the existing notebook renders it as "n.s." gray.
5. The 2.83x row-duplication from the many-to-many merge (1_preprocessing.py:246) is fixed by
   de-duplicating the energy substrate up front.

SUMMARY-STATISTIC SENSITIVITY
   Both energy hit lists are emitted every run: the primary (SUMMARY) under the canonical name,
   and the alternative under a '_<other>sens' name. The primary *_FULL table carries the
   alternative's p/q and a `robust_to_summary` flag. This matters because the B*38:01 APO state
   is genuinely floppy (median within-trajectory SD ~19 kJ/mol vs ~0.6 elsewhere), so mean and
   median diverge there (B*38:01 loaded-vs-free: ~643 hits by mean, ~305-315 by median) while the
   stable systems are unaffected. The gap is itself a readout of that apo instability.

EQUILIBRATION
   Energies now use the same equilibrated window as BC (frames with index >= EQUIL_MIN_FRAME;
   BC used Frame > 20 -> 219 of 239 frames). The saved df_cons_pairs_saved keeps ALL frame
   columns so the notebook's polymorphic KDEs are unchanged; only the inference uses the window.

WHY NOT EMPIRICAL-BAYES MODERATION (limma): MD pair-energy variances span ~4e6 across pairs and
   BC ~2e9, so the EB prior df collapses to <1 and moderation degenerates to anti-conservative
   equal-variance Student's t. The fitted prior df is reported at run time; we use Welch.

INPUTS (produced unchanged by the current pipeline; no /mnt/MORIA or PDBs needed):
    df_cons_pairs_saved_2025_07_10.csv   annotated per-(pair,system,repeat) VdW energies
    df_bc_equil_2025_07_10.csv           per-frame betweenness centrality

OUTPUTS (written to OUTPUT_DIR, date-stamped DATE):
    df_cons_pairs_saved_<DATE>.csv            de-duplicated substrate (all frame columns kept)
    df_sig_aff_pairs_loaded_<DATE>.csv        corrected 1.4 hit list (primary summary)
    df_sig_aff_pairs_loaded_<DATE>_<alt>sens.csv   sensitivity hit list (other summary)
    df_sig_aff_pairs_loaded_<DATE>_FULL.csv   every gated contrast + effect/CI/q + sensitivity cols
    df_sig_aff_pairs_b39_<DATE>.csv           corrected 1.5 hit list (+ _<alt>sens, _FULL)
    df_bc_equil_cons_resids_<DATE>.csv        BC + reproducibility flag
    df_bc_sigaff_resids_<DATE>.csv            corrected 1.8 (Kruskal omnibus + Welch dBC)
    df_bc_change_{omnibus,contrasts}_<DATE>_FULL.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.special import digamma, polygamma

# ----------------------------- configuration --------------------------------
# The v1 annotated-energy substrate lives in the backup folder (it is a v1 output consumed here
# as input; v2 de-duplicates it). df_bc_equil is a valid intermediate kept in the repo root.
IN_CONS_PAIRS = "backup_v1_faulty_stats/df_cons_pairs_saved_2025_07_10.csv"
IN_BC_EQUIL   = "df_bc_equil_2025_07_10.csv"
OUT_MAIN      = "."                   # corrected primary + sensitivity outputs -> repo root
OUT_SUPP      = "v2_supplementary"    # per-contrast FULL / omnibus / contrast detail tables
DATE          = "2026_07_27"         # date stamp for the corrected outputs
SUMMARY       = "median"             # per-replicate energy summary: 'median' (robust) or 'mean'
EFFECT_GATE   = 2.0                  # |Delta| kJ/mol, matches 1_preprocessing.py:368
ALPHA         = 0.05
MIN_REPS      = 3                    # require >=3 replicate trajectories per group
EQUIL_MIN_FRAME = 20                 # keep energy frames with index >= this (aligns to BC Frame>20)
WRITE_BC_EQUIL_CONS = True

ALLELES  = ["B3801", "B3901", "B3906"]
PEPTIDES = ["A", "AL", "ALL"]
SUMMARY_COL = {"median": "rep_med", "mean": "rep_mean"}
ALT = "mean" if SUMMARY == "median" else "median"


# ======================= statistics core ====================================
def benjamini_hochberg(pvals):
    p = np.asarray(pvals, float)
    q = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    pm = p[m]; n = pm.size
    if n == 0:
        return q
    o = np.argsort(pm); r = pm[o]
    raw = np.minimum.accumulate((r * n / np.arange(1, n + 1))[::-1])[::-1]
    qm = np.empty(n); qm[o] = np.clip(raw, 0, 1); q[m] = qm
    return q


def welch(a, b):
    """Welch t-test. Returns (effect=mean(a)-mean(b), p, se, df, ci_low, ci_high)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = a.size, b.size
    va, vb = a.var(ddof=1), b.var(ddof=1)
    effect = a.mean() - b.mean()
    se = np.sqrt(va / na + vb / nb)
    if not np.isfinite(se) or se == 0:
        return effect, 1.0, 0.0, np.nan, effect, effect
    with np.errstate(divide="ignore", invalid="ignore"):
        df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    t = effect / se
    p = 2 * st.t.sf(abs(t), df)
    tcrit = st.t.ppf(0.975, df)
    return effect, p, se, df, effect - tcrit * se, effect + tcrit * se


def _trigamma_inverse(x):
    if not np.isfinite(x) or x <= 0:
        return np.inf
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    y = 0.5 + 1.0 / x
    for _ in range(60):
        tri = polygamma(1, y)
        dif = tri * (1 - tri / x) / polygamma(2, y)
        y += dif
        if -dif / y < 1e-8:
            break
    return y


def eb_prior_df0(s2, dfres):
    """Fitted prior d.f. of the limma variance model; df0<1 => moderation degenerate (reported only)."""
    s2 = np.asarray(s2, float); dfres = np.asarray(dfres, float)
    ok = np.isfinite(s2) & (s2 > 0) & np.isfinite(dfres) & (dfres > 0)
    s2, dfres = s2[ok], dfres[ok]
    if s2.size < 2:
        return np.inf
    e = np.log(s2) - digamma(dfres / 2.0) + np.log(dfres / 2.0)
    evar = e.var(ddof=1) - np.mean(polygamma(1, dfres / 2.0))
    return 2.0 * _trigamma_inverse(evar) if evar > 0 else np.inf


# ======================= energy scans (1.4, 1.5) ============================
def load_energy_substrate(path):
    df = pd.read_csv(path, low_memory=False)
    fc = [c for c in df.columns if c.isdigit()]
    n_before = len(df)
    df = df.drop_duplicates(subset=["resnum_chain12", "system", "repeat"], keep="first").reset_index(drop=True)
    fc_eq = [c for c in fc if int(c) >= EQUIL_MIN_FRAME]          # equilibrated window (align to BC)
    Xeq = df[fc_eq].to_numpy(float)
    df["rep_mean"] = np.nanmean(Xeq, axis=1)
    df["rep_med"] = np.nanmedian(Xeq, axis=1)
    print(f"[energy] {path}: de-dup {n_before} -> {len(df)} ({n_before/len(df):.2f}x); "
          f"frames {len(fc)} -> {len(fc_eq)} equilibrated (index >= {EQUIL_MIN_FRAME}); "
          f"reps/system {df.groupby('system')['repeat'].nunique().min()}-"
          f"{df.groupby('system')['repeat'].nunique().max()}")
    return df, fc


def build_contrasts(df, kind):
    out = []
    for pair, dfp in df.groupby("resnum_chain12", sort=False):
        systems = set(dfp["system"].unique())
        if kind == "loaded_vs_free":
            for allele in dfp["allele"].unique():
                if allele in systems:
                    for pep in PEPTIDES:
                        if f"{allele}_{pep}" in systems:
                            out.append((pair, allele, pep, f"{allele}_{pep}", allele))
        elif kind == "b39_vs_b38":
            for pep in PEPTIDES:
                ref = f"B3801_{pep}"
                if ref in systems:
                    for allele in ("B3901", "B3906"):
                        if f"{allele}_{pep}" in systems:
                            out.append((pair, allele, pep, f"{allele}_{pep}", ref))
    return out


def analyse_energy(df, kind, summary):
    col = SUMMARY_COL[summary]
    gm = {k: g[col].to_numpy(float) for k, g in df.groupby(["resnum_chain12", "system"])}
    rows, s2s, dfs = [], [], []
    for pair, allele, pep, test_sys, ref_sys in build_contrasts(df, kind):
        a = gm.get((pair, test_sys)); b = gm.get((pair, ref_sys))
        if a is None or b is None or len(a) < MIN_REPS or len(b) < MIN_REPS:
            continue
        na, nb = len(a), len(b)
        effect, p, se, dfw, lo, hi = welch(a, b)
        rows.append(dict(allele=allele, peptide=pep, system=test_sys, ref_system=ref_sys,
                         pair=pair, change_type="attractive" if effect < 0 else "repulsive",
                         delta_mean=effect, p_value=p, ci_low=lo, ci_high=hi, welch_df=dfw,
                         n_ref_rep=nb, n_test_rep=na))
        s2s.append(((na-1)*a.var(ddof=1) + (nb-1)*b.var(ddof=1)) / (na+nb-2))
        dfs.append(na + nb - 2)
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["gated"] = res.delta_mean.abs() >= EFFECT_GATE
    res["q_value"] = np.nan
    res.loc[res.gated, "q_value"] = benjamini_hochberg(res.loc[res.gated, "p_value"])
    res.attrs["eb_df0"] = eb_prior_df0(np.array(s2s)[res.gated.values], np.array(dfs)[res.gated.values])
    return res


HIT_COLS = ["allele", "peptide", "system", "pair", "change_type", "p_value", "q_value",
            "delta_mean", "ci_low", "ci_high", "n_ref_rep", "n_test_rep"]


def run_energy_scan(df, kind, mainp, suppp, base):
    """Run the scan under both summaries; write primary hit list, sensitivity hit list, FULL."""
    prim = analyse_energy(df, kind, SUMMARY)
    alt = analyse_energy(df, kind, ALT)
    df0 = prim.attrs.get("eb_df0", float("nan"))
    key = ["allele", "peptide", "pair"]
    a = alt.set_index(key)
    prim = prim.join(a[["delta_mean", "p_value", "q_value", "gated"]]
                     .rename(columns={"delta_mean": f"delta_{ALT}sens", "p_value": f"p_{ALT}sens",
                                      "q_value": f"q_{ALT}sens", "gated": f"gated_{ALT}sens"}), on=key)
    prim["robust_to_summary"] = (prim.gated & (prim.q_value < ALPHA)
                                 & prim[f"gated_{ALT}sens"].fillna(False) & (prim[f"q_{ALT}sens"] < ALPHA))
    prim_hits = prim[prim.gated & (prim.q_value < ALPHA)][HIT_COLS].reset_index(drop=True)
    alt_hits = alt[alt.gated & (alt.q_value < ALPHA)][HIT_COLS].reset_index(drop=True)
    prim_hits.to_csv(mainp(f"{base}_{DATE}.csv"), index=False)
    alt_hits.to_csv(mainp(f"{base}_{DATE}_{ALT}sens.csv"), index=False)
    prim.to_csv(suppp(f"{base}_{DATE}_FULL.csv"), index=False)
    print(f"[{base}] primary({SUMMARY})={len(prim_hits)} hits, sensitivity({ALT})={len(alt_hits)}; "
          f"robust-to-summary={int(prim.robust_to_summary.sum())}; "
          f"gated={int(prim.gated.sum())}; EB df0={df0:.2f} (moderation not used)")
    return prim_hits, alt_hits


# ======================= BC scan (1.7, 1.8) =================================
def analyse_bc(path):
    print(f"[bc] loading {path} (3.9M rows; ~30s) ...")
    bc = pd.read_csv(path, low_memory=False)
    repmed = (bc.groupby(["Residue_id", "system", "repeat"])["BC"].median()
                .rename("rep_med").reset_index())
    present = (repmed.groupby(["Residue_id", "system"])["repeat"].size()
               .rename("n_rep").reset_index())
    present["reproducible"] = present["n_rep"] >= MIN_REPS
    gm = {(r, s): g["rep_med"].to_numpy(float) for (r, s), g in repmed.groupby(["Residue_id", "system"])}
    usable = present[present.reproducible].groupby("Residue_id")["system"].size()
    cons_res = usable[usable >= 2].index.tolist()
    print(f"[bc] {len(cons_res)} reproducible residues (>= {MIN_REPS} reps in >= 2 systems)")

    all_sys = ALLELES + [f"{a}_{p}" for a in ALLELES for p in PEPTIDES]
    om = []
    for r in cons_res:
        groups = [gm[(r, s)] for s in all_sys if (r, s) in gm and len(gm[(r, s)]) >= MIN_REPS]
        if len(groups) < 2:
            continue
        try:
            p = st.kruskal(*groups).pvalue
        except ValueError:
            p = np.nan
        om.append((r, p))
    om = pd.DataFrame(om, columns=["Residue_id", "p_omni"])
    om["q_omni"] = benjamini_hochberg(om["p_omni"])
    sig_res = om[om.q_omni < ALPHA]["Residue_id"].tolist()
    print(f"[bc] omnibus (Kruskal + FDR): {len(sig_res)}/{len(om)} residues change somewhere")

    crows = []
    for r in sig_res:
        for a in ALLELES:
            free = gm.get((r, a))
            if free is None or len(free) < MIN_REPS:
                continue
            for p in PEPTIDES:
                load = gm.get((r, f"{a}_{p}"))
                if load is None or len(load) < MIN_REPS:
                    continue
                effect, pp, se, dfw, lo, hi = welch(load, free)
                crows.append(dict(Residue_id=r, allele=a, peptide=p, diff=effect,
                                  ref_bc=float(np.mean(free)), p_value=pp, ci_low=lo, ci_high=hi))
    cdf = pd.DataFrame(crows)
    cdf["q_value"] = benjamini_hochberg(cdf["p_value"])
    print(f"[bc] per-contrast dBC (Welch + FDR): {int((cdf.q_value < ALPHA).sum())}/{len(cdf)} significant")

    om_idx = om.set_index("Residue_id")
    wide = []
    for r in sig_res:
        row = {"Residue_id": r, "p_value": om_idx.at[r, "q_omni"]}
        for a in ALLELES:
            free = gm.get((r, a))
            row[f"median_bc_{a}_nopep"] = float(np.mean(free)) if free is not None else np.nan
            for p in PEPTIDES:
                sub = cdf[(cdf.Residue_id == r) & (cdf.allele == a) & (cdf.peptide == p)]
                if len(sub):
                    s = sub.iloc[0]; sig = s.q_value < ALPHA
                    row[f"diff_bc_{a}_{p}"] = s["diff"] if sig else np.nan
                    row[f"diff_bc_raw_{a}_{p}"] = s["diff"]
                    row[f"q_bc_{a}_{p}"] = s.q_value
                    row[f"pt_diff_bc_{a}_{p}"] = (s["diff"] / abs(s.ref_bc)
                                                  if sig and s.ref_bc != 0 else np.nan)
                else:
                    for c in (f"diff_bc_{a}_{p}", f"diff_bc_raw_{a}_{p}", f"q_bc_{a}_{p}", f"pt_diff_bc_{a}_{p}"):
                        row[c] = np.nan
        wide.append(row)
    return pd.DataFrame(wide), om, cdf, present, bc


# ============================== main ========================================
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    os.makedirs(OUT_SUPP, exist_ok=True)
    mainp = lambda f: os.path.join(OUT_MAIN, f)
    suppp = lambda f: os.path.join(OUT_SUPP, f)

    print("=" * 78)
    print(f"1_preprocessing_v2  —  corrected statistics layer  (SUMMARY={SUMMARY}, DATE={DATE})")
    print("=" * 78)

    edf, fc = load_energy_substrate(IN_CONS_PAIRS)
    # de-duplicated substrate, ALL frame columns kept (notebook polymorphic KDEs unchanged)
    edf.drop(columns=["rep_mean", "rep_med"]).to_csv(mainp(f"df_cons_pairs_saved_{DATE}.csv"), index=False)

    run_energy_scan(edf, "loaded_vs_free", mainp, suppp, "df_sig_aff_pairs_loaded")
    run_energy_scan(edf, "b39_vs_b38", mainp, suppp, "df_sig_aff_pairs_b39")

    wide, om, cdf, present, bc = analyse_bc(IN_BC_EQUIL)
    wide.to_csv(mainp(f"df_bc_sigaff_resids_{DATE}.csv"), index=False)
    om.to_csv(suppp(f"df_bc_change_omnibus_{DATE}_FULL.csv"), index=False)
    cdf.to_csv(suppp(f"df_bc_change_contrasts_{DATE}_FULL.csv"), index=False)
    if WRITE_BC_EQUIL_CONS:
        flag = present[present.reproducible][["Residue_id", "system"]].copy()
        flag["p_value"] = 1.0
        bc.merge(flag, on=["Residue_id", "system"], how="left") \
          .to_csv(mainp(f"df_bc_equil_cons_resids_{DATE}.csv"), index=False)
        print(f"[bc] wrote df_bc_equil_cons_resids_{DATE}.csv")

    print(f"\nCorrected outputs (date {DATE}): primary in {OUT_MAIN}/, detail in {OUT_SUPP}/. "
          f"Primary summary = {SUMMARY}; '{ALT}' emitted as sensitivity.")
