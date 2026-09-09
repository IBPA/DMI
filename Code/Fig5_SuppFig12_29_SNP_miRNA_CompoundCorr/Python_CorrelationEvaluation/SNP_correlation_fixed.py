import os
import time
import gzip
import csv
import argparse
import warnings
from multiprocessing import cpu_count, get_context
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

warnings.filterwarnings("ignore")

SNP_FILE = "SNP_ad_freq.csv"
COMPOUND_FILE = "compound_conc.csv"
PCA_TSNE_FILE = "compound_conc.csv"
COMPOUND_SUBSET_FILE = None
EXPECTED_N_COMPOUNDS: Optional[int] = 5310
REQUIRE_SUBSET = True

OUTPUT_DIR = "correleation_results_snp"

MIN_SAMPLES = 10

SIG_ABS_R = 0.5
Q_THRESHOLD = 0.05
SIG_MIN_OVERLAP_DEFAULT = 30

THRESHOLDS = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

N_CORES = max(1, int(cpu_count() * 0.8))
DTYPE_FLOAT = np.float32
CSV_FLOAT_FORMAT = "%.6g"
ROW_CHUNK = 250

FLAGGED_SNP = "GK000026.2_43287817_G_A"
FLAGGED_COMPOUND = "DMD303043"

DEFAULT_EXACT_QVALUE_MAX_PAIRS = 50_000_000
DEFAULT_HIST_BINS = 200_000

_G_SNP: Optional[np.ndarray] = None
_G_CMP: Optional[np.ndarray] = None
_G_MIN: int = 10


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--quick-test", action="store_true")
    p.add_argument("--test-snps", type=int, default=200)
    p.add_argument("--test-compounds", type=int, default=200)
    p.add_argument("--test-seed", type=int, default=123)
    p.add_argument("--no-parallel", action="store_true")
    p.add_argument("--qmode", choices=["auto", "exact", "hist", "skip"], default="auto")
    p.add_argument("--hist-bins", type=int, default=DEFAULT_HIST_BINS)
    p.add_argument("--sig-min-overlap", type=int, default=SIG_MIN_OVERLAP_DEFAULT)
    p.add_argument("--no-tests", action="store_true")
    p.add_argument("--test-k", type=int, default=8)
    p.add_argument("--test-seed2", type=int, default=999)
    p.add_argument("--strict-tests", action="store_true")
    return p.parse_args()


def _init_worker(snp_arr: np.ndarray, cmp_arr: np.ndarray, min_samples: int) -> None:
    global _G_SNP, _G_CMP, _G_MIN
    _G_SNP = snp_arr
    _G_CMP = cmp_arr
    _G_MIN = int(min_samples)


def _compute_for_snp(snp_idx: int) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    assert _G_SNP is not None and _G_CMP is not None
    snp = _G_SNP[:, snp_idx]
    cmp = _G_CMP
    min_n = _G_MIN
    n_compounds = cmp.shape[1]

    pcc = np.full(n_compounds, np.nan, dtype=DTYPE_FLOAT)
    spr = np.full(n_compounds, np.nan, dtype=DTYPE_FLOAT)
    p_pcc = np.full(n_compounds, np.nan, dtype=DTYPE_FLOAT)
    p_spr = np.full(n_compounds, np.nan, dtype=DTYPE_FLOAT)
    n_ov = np.zeros(n_compounds, dtype=np.uint16)

    snp_mask = ~np.isnan(snp)
    if int(snp_mask.sum()) < min_n:
        return snp_idx, pcc, spr, p_pcc, p_spr, n_ov

    snp0 = snp[snp_mask]
    if np.std(snp0) == 0:
        return snp_idx, pcc, spr, p_pcc, p_spr, n_ov

    for j in range(n_compounds):
        y = cmp[:, j]
        mask = snp_mask & ~np.isnan(y)
        n = int(mask.sum())
        n_ov[j] = n
        if n < min_n:
            continue
        x = snp[mask]
        yy = y[mask]
        if np.std(yy) == 0:
            continue
        try:
            r1, pv1 = pearsonr(x, yy)
            r2, pv2 = spearmanr(x, yy)
            pcc[j] = DTYPE_FLOAT(r1)
            p_pcc[j] = DTYPE_FLOAT(pv1)
            spr[j] = DTYPE_FLOAT(r2)
            p_spr[j] = DTYPE_FLOAT(pv2)
        except Exception:
            pass

    return snp_idx, pcc, spr, p_pcc, p_spr, n_ov


def _read_numeric_matrix(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").astype(DTYPE_FLOAT)
    return df


def _load_subset_ids(path: str) -> List[str]:
    if path.lower().endswith(".txt"):
        with open(path, "r", encoding="utf-8") as f:
            return [ln.strip() for ln in f if ln.strip()]
    sub = pd.read_csv(path, header=None)
    return sub.iloc[:, 0].astype(str).tolist()


def _write_subset_ids(path: str, ids: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for x in ids:
            x = str(x).strip()
            if x:
                f.write(x + "\n")


def _pca_header_compound_ids(pca_tsne_path: str) -> List[str]:
    header_df = pd.read_csv(pca_tsne_path, nrows=0)
    cols = list(header_df.columns)
    if len(cols) < 2:
        raise ValueError(f"Bad PCA/TSNE header: {pca_tsne_path}")
    out = []
    for c in cols[1:]:
        c = str(c).strip()
        if not c:
            continue
        if c.lower().startswith("unnamed"):
            continue
        out.append(c)
    return out


def ensure_subset_file(subset_path: str, pca_tsne_path: str, expected_n: Optional[int]) -> None:
    if os.path.exists(subset_path):
        return
    if not os.path.exists(pca_tsne_path):
        raise FileNotFoundError(f"Missing subset file and PCA/TSNE file: {subset_path} / {pca_tsne_path}")
    ids = _pca_header_compound_ids(pca_tsne_path)
    if expected_n is not None and len(ids) != expected_n:
        raise ValueError(f"Auto subset has {len(ids)} IDs, expected {expected_n}")
    _write_subset_ids(subset_path, ids)
    print(f"✓ Auto-generated {subset_path} from {pca_tsne_path} (n={len(ids)})")


def _gzip_csv_writer(path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return gzip.open(path, "wt", newline="")


def _save_matrix_csv_gz(memmap_arr: np.memmap, row_ids: pd.Index, col_ids: pd.Index, out_path: str) -> None:
    with _gzip_csv_writer(out_path) as f:
        f.write("," + ",".join(map(str, col_ids)) + "\n")
        n_rows = memmap_arr.shape[0]
        for i0 in tqdm(range(0, n_rows, ROW_CHUNK), desc=f"Saving {os.path.basename(out_path)}"):
            i1 = min(n_rows, i0 + ROW_CHUNK)
            chunk = memmap_arr[i0:i1, :]
            for r, rid in enumerate(row_ids[i0:i1]):
                row = chunk[r, :]
                row_str = ",".join("" if not np.isfinite(v) else (CSV_FLOAT_FORMAT % float(v)) for v in row)
                f.write(str(rid) + "," + row_str + "\n")


def _save_matrix_uint_csv_gz(memmap_arr: np.memmap, row_ids: pd.Index, col_ids: pd.Index, out_path: str) -> None:
    with _gzip_csv_writer(out_path) as f:
        f.write("," + ",".join(map(str, col_ids)) + "\n")
        n_rows = memmap_arr.shape[0]
        for i0 in tqdm(range(0, n_rows, ROW_CHUNK), desc=f"Saving {os.path.basename(out_path)}"):
            i1 = min(n_rows, i0 + ROW_CHUNK)
            chunk = memmap_arr[i0:i1, :]
            for r, rid in enumerate(row_ids[i0:i1]):
                row = chunk[r, :]
                row_str = ",".join(str(int(v)) for v in row)
                f.write(str(rid) + "," + row_str + "\n")


def plot_pair_scatter(
    snp_id: str,
    compound_id: str,
    snp_df: pd.DataFrame,
    compound_df: pd.DataFrame,
    sample_order: List[str],
    out_png: str,
) -> Dict[str, float]:
    x = snp_df.loc[snp_id, sample_order]
    y = compound_df.loc[compound_id, sample_order]
    mask = x.notna() & y.notna()
    xv = x[mask].astype(float).values
    yv = y[mask].astype(float).values
    n = len(xv)

    if n >= 3 and np.std(xv) > 0 and np.std(yv) > 0:
        r_p, p_p = pearsonr(xv, yv)
        r_s, p_s = spearmanr(xv, yv)
    else:
        r_p, p_p, r_s, p_s = np.nan, np.nan, np.nan, np.nan

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.figure()
    plt.scatter(xv, yv)
    plt.xlabel(snp_id)
    plt.ylabel(compound_id)
    plt.title(f"N={n} | PCC={r_p:.4g} (p={p_p:.3g}) | Spearman={r_s:.4g} (p={p_s:.3g})")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return {"N": float(n), "PCC": float(r_p), "PCC_p": float(p_p), "SPR": float(r_s), "SPR_p": float(p_s)}


def _estimate_pi0_storey(p: np.ndarray, lambdas: np.ndarray = None) -> float:
    if lambdas is None:
        lambdas = np.arange(0.05, 0.95, 0.05, dtype=np.float64)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return 1.0
    vals = []
    for lam in lambdas:
        vals.append(np.mean(p > lam) / max(1e-12, (1.0 - lam)))
    pi0 = float(min(vals))
    return float(min(max(pi0, 0.0), 1.0))


def qvalues_storey_exact(p: np.ndarray) -> np.ndarray:
    p = p.astype(np.float64, copy=False)
    m = p.size
    pi0 = _estimate_pi0_storey(p)
    order = np.argsort(p)
    p_sorted = p[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)
    q_sorted = (pi0 * m * p_sorted) / ranks
    for i in range(m - 2, -1, -1):
        if q_sorted[i] > q_sorted[i + 1]:
            q_sorted[i] = q_sorted[i + 1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return q.astype(np.float32, copy=False)


def compute_qvalues_exact_from_memmap(pvals_mm: np.memmap, qvals_mm: np.memmap, label: str) -> None:
    t0 = time.time()
    flat = np.array(pvals_mm.reshape(-1), dtype=np.float64, copy=True)
    nan_mask = ~np.isfinite(flat)
    flat[nan_mask] = 1.0
    q = qvalues_storey_exact(flat)
    q[nan_mask] = np.nan
    qvals_mm[:] = q.reshape(pvals_mm.shape)
    qvals_mm.flush()
    print(f"✓ Exact {label} q-values computed in {(time.time() - t0)/60:.2f} minutes")


def compute_qvalues_hist_from_memmap(pvals_mm: np.memmap, qvals_mm: np.memmap, label: str, bins: int) -> None:
    t0 = time.time()
    bins = int(bins)
    if bins < 2000:
        raise ValueError("hist-bins must be >= 2000")

    hist = np.zeros(bins, dtype=np.int64)
    m = 0
    n_rows = pvals_mm.shape[0]

    for i0 in tqdm(range(0, n_rows, ROW_CHUNK), desc=f"[{label}] Histogram pass"):
        i1 = min(n_rows, i0 + ROW_CHUNK)
        chunk = np.array(pvals_mm[i0:i1, :], copy=False)
        flat = chunk.reshape(-1)
        finite = flat[np.isfinite(flat)]
        if finite.size == 0:
            continue
        finite = np.clip(finite.astype(np.float64, copy=False), 0.0, 1.0)
        idx = np.minimum((finite * bins).astype(np.int64), bins - 1)
        hist += np.bincount(idx, minlength=bins)
        m += int(finite.size)

    if m == 0:
        qvals_mm[:] = np.nan
        qvals_mm.flush()
        print(f"⚠ [{label}] No finite p-values; q-values set to NaN.")
        return

    lambdas = np.arange(0.05, 0.95, 0.05, dtype=np.float64)
    tail = np.cumsum(hist[::-1])[::-1]
    pi0_vals = []
    for lam in lambdas:
        k = min(int(lam * bins), bins - 1)
        k2 = min(k + 1, bins - 1)
        tail_count = int(tail[k2])
        pi0_vals.append((tail_count / m) / max(1e-12, (1.0 - lam)))
    pi0 = float(min(pi0_vals))
    pi0 = float(min(max(pi0, 0.0), 1.0))

    cdf = np.cumsum(hist).astype(np.float64)
    ranks = np.maximum(cdf, 1.0)
    p_rep = (np.arange(bins, dtype=np.float64) + 0.5) / bins
    q_bin = (pi0 * m * p_rep) / ranks
    for k in range(bins - 2, -1, -1):
        if q_bin[k] > q_bin[k + 1]:
            q_bin[k] = q_bin[k + 1]
    q_bin = np.clip(q_bin, 0.0, 1.0).astype(np.float32)

    for i0 in tqdm(range(0, n_rows, ROW_CHUNK), desc=f"[{label}] Assign q-values"):
        i1 = min(n_rows, i0 + ROW_CHUNK)
        chunk = np.array(pvals_mm[i0:i1, :], copy=False)
        out = np.empty_like(chunk, dtype=np.float32)

        flat = chunk.reshape(-1)
        out_flat = out.reshape(-1)
        finite_mask = np.isfinite(flat)
        out_flat[~finite_mask] = np.nan
        if np.any(finite_mask):
            pvals = np.clip(flat[finite_mask].astype(np.float64, copy=False), 0.0, 1.0)
            idx = np.minimum((pvals * bins).astype(np.int64), bins - 1)
            out_flat[finite_mask] = q_bin[idx]

        qvals_mm[i0:i1, :] = out

    qvals_mm.flush()
    print(f"✓ [{label}] Hist q-values computed in {(time.time() - t0)/60:.2f} minutes (pi0={pi0:.4g}, bins={bins})")


def choose_subset(index: pd.Index, k: int, seed: int) -> pd.Index:
    k = min(int(k), len(index))
    rng = np.random.default_rng(int(seed))
    pos = rng.choice(np.arange(len(index)), size=k, replace=False)
    return pd.Index(index[pos])


def ensure_included(base: pd.Index, must_include: Optional[str], k: int) -> pd.Index:
    if must_include is None:
        return base[:k]
    if must_include in base:
        return base[:k]
    merged = pd.Index([must_include]).append(base)
    seen = set()
    out = []
    for x in merged:
        if x not in seen:
            out.append(x)
            seen.add(x)
        if len(out) >= k:
            break
    return pd.Index(out)


def post_run_tests(
    out_dir: str,
    snp_df: pd.DataFrame,
    cmp_df: pd.DataFrame,
    common_samples: List[str],
    subset_path: str,
    pca_tsne_path: str,
    pcc_mm: np.memmap,
    spr_mm: np.memmap,
    p_pcc_mm: np.memmap,
    p_spr_mm: np.memmap,
    n_ov_mm: np.memmap,
    sig_pairs_path: str,
    args,
) -> Tuple[bool, str]:
    lines = []
    ok = True

    def check(cond: bool, ok_msg: str, fail_msg: str):
        nonlocal ok
        if cond:
            lines.append(f"[OK]   {ok_msg}")
        else:
            lines.append(f"[FAIL] {fail_msg}")
            ok = False

    lines.append("=" * 80)
    lines.append("POST-RUN TEST REPORT (paste this back to me)")
    lines.append("=" * 80)

    check(os.path.exists(subset_path), "Subset file exists", f"Subset file missing: {subset_path}")
    if os.path.exists(subset_path):
        ids = [str(x).strip() for x in _load_subset_ids(subset_path)]
        ids = [x for x in ids if x]
        check(len(ids) == EXPECTED_N_COMPOUNDS,
              f"Subset has {len(ids)} IDs (expected {EXPECTED_N_COMPOUNDS})",
              f"Subset has {len(ids)} IDs (expected {EXPECTED_N_COMPOUNDS})")

        if os.path.exists(pca_tsne_path):
            pca_ids = _pca_header_compound_ids(pca_tsne_path)
            check(len(pca_ids) == EXPECTED_N_COMPOUNDS,
                  f"PCA/TSNE header has {len(pca_ids)} compound IDs",
                  f"PCA/TSNE header has {len(pca_ids)} compound IDs (expected {EXPECTED_N_COMPOUNDS})")
            check(pca_ids == ids,
                  "Subset list matches PCA/TSNE header order exactly",
                  "Subset list does NOT match PCA/TSNE header order exactly")

        cmp_ids = set(cmp_df.index.astype(str))
        inter = len(set(ids) & cmp_ids)
        expected_here = len(cmp_df.index)
        check(inter == expected_here,
              f"Subset IDs present in current run compound matrix: {inter}/{expected_here}",
              f"Subset intersection mismatch: {inter}/{expected_here}")

    check(len(common_samples) > 0, f"Common samples = {len(common_samples)}", "No common samples found")
    check(list(snp_df.columns) == common_samples,
          "SNP sample order preserved after alignment",
          "SNP sample order mismatch after alignment")
    check(list(cmp_df.columns) == common_samples,
          "Compound sample columns match SNP sample order",
          "Compound sample order mismatch")

    expected_files = [
        "pcc_correlation_matrix.csv.gz",
        "spearman_correlation_matrix.csv.gz",
        "pcc_pvalues.csv.gz",
        "spearman_pvalues.csv.gz",
        "n_overlap_matrix.csv.gz",
        "threshold_analysis.csv.gz",
        "significant_pairs.csv.gz",
        "summary_statistics.csv.gz",
    ]
    for fn in expected_files:
        path = os.path.join(out_dir, fn)
        check(os.path.exists(path), f"Found output file {fn}", f"Missing output file {fn}")

    thr_path = os.path.join(out_dir, "threshold_analysis.csv.gz")
    if os.path.exists(thr_path):
        try:
            thr = pd.read_csv(thr_path)
            check(len(thr) == 6, "threshold_analysis has 6 rows", f"threshold_analysis rows={len(thr)} (expected 6)")
            req = {"Threshold", "Pairs_Passing_Both", "Percentage_of_All_Pairs"}
            check(req.issubset(set(thr.columns)), "threshold_analysis has required columns",
                  f"threshold_analysis missing columns: {req - set(thr.columns)}")
        except Exception as e:
            check(False, "", f"threshold_analysis not readable: {e}")

    sum_path = os.path.join(out_dir, "summary_statistics.csv.gz")
    if os.path.exists(sum_path):
        try:
            s = pd.read_csv(sum_path)
            check(len(s) == 1, "summary_statistics has 1 row", f"summary_statistics rows={len(s)} (expected 1)")
        except Exception as e:
            check(False, "", f"summary_statistics not readable: {e}")

    sig_path = sig_pairs_path
    if os.path.exists(sig_path):
        try:
            sig = pd.read_csv(sig_path)
            required_cols = {"SNP", "Compound", "N_overlap", "PCC_r", "PCC_p", "PCC_q", "Spearman_r", "Spearman_p", "Spearman_q"}
            check(required_cols.issubset(set(sig.columns)),
                  "significant_pairs has required columns",
                  f"significant_pairs missing columns: {required_cols - set(sig.columns)}")
            if len(sig) > 0 and "N_overlap" in sig.columns:
                bad = int((sig["N_overlap"] < int(args.sig_min_overlap)).sum())
                check(bad == 0,
                      "All significant pairs satisfy sig-min-overlap",
                      f"{bad} significant pairs violate sig-min-overlap")
            lines.append(f"[INFO] significant_pairs rows: {len(sig)}")
        except Exception as e:
            check(False, "", f"significant_pairs not readable: {e}")

    nov = np.array(n_ov_mm[:], copy=False).reshape(-1)
    nov_pos = nov[nov > 0]
    if nov_pos.size > 0:
        mn, med, mx = int(nov_pos.min()), int(np.median(nov_pos)), int(nov_pos.max())
        lines.append(f"[INFO] N_overlap (positive) min/median/max: {mn}/{med}/{mx}")
        check(mx <= len(common_samples), "N_overlap max <= #samples", f"N_overlap max {mx} > #samples {len(common_samples)}")
        nov_comp = nov[nov >= MIN_SAMPLES]
        if nov_comp.size > 0:
            mn2, med2, mx2 = int(nov_comp.min()), int(np.median(nov_comp)), int(nov_comp.max())
            lines.append(f"[INFO] N_overlap (computed >= {MIN_SAMPLES}) min/median/max: {mn2}/{med2}/{mx2}")
        else:
            lines.append(f"[INFO] No pairs with N_overlap >= {MIN_SAMPLES}")
    else:
        check(False, "", "No overlap values recorded")

    rng = np.random.default_rng(int(args.test_seed2))
    K = int(max(1, args.test_k))
    n_snps, n_cmp = pcc_mm.shape

    candidates = []
    for _ in range(4000):
        i = int(rng.integers(0, n_snps))
        j = int(rng.integers(0, n_cmp))
        rp = float(pcc_mm[i, j])
        rs = float(spr_mm[i, j])
        if np.isfinite(rp) and np.isfinite(rs):
            candidates.append((i, j))
            if len(candidates) >= K:
                break

    check(len(candidates) >= 1, f"Found {len(candidates)} finite entries to spot-check", "Could not find finite entries to spot-check")

    tol_r = 2e-3
    tol_p = 5e-6

    for idx, (i, j) in enumerate(candidates[:K], start=1):
        snp_id = str(snp_df.index[i])
        cmp_id = str(cmp_df.index[j])

        x = snp_df.loc[snp_id, common_samples]
        y = cmp_df.loc[cmp_id, common_samples]
        mask = x.notna() & y.notna()
        xv = x[mask].astype(float).values
        yv = y[mask].astype(float).values

        if len(xv) < MIN_SAMPLES or np.std(xv) == 0 or np.std(yv) == 0:
            check(False, "", f"Spot-check {idx}: recompute invalid but matrix finite for ({snp_id}, {cmp_id})")
            continue

        rp2, pp2 = pearsonr(xv, yv)
        rs2, ps2 = spearmanr(xv, yv)

        rp1 = float(pcc_mm[i, j])
        rs1 = float(spr_mm[i, j])
        pp1 = float(p_pcc_mm[i, j])
        ps1 = float(p_spr_mm[i, j])
        n1 = int(n_ov_mm[i, j])

        cond = (
            abs(rp1 - rp2) <= tol_r and
            abs(rs1 - rs2) <= tol_r and
            abs(pp1 - pp2) <= tol_p and
            abs(ps1 - ps2) <= tol_p and
            n1 == int(len(xv))
        )
        check(cond,
              f"Spot-check {idx}: OK ({snp_id} vs {cmp_id})",
              f"Spot-check {idx}: MISMATCH ({snp_id} vs {cmp_id}) | "
              f"rP {rp1:.4g} vs {rp2:.4g}, rS {rs1:.4g} vs {rs2:.4g}, "
              f"pP {pp1:.3g} vs {pp2:.3g}, pS {ps1:.3g} vs {ps2:.3g}, "
              f"N {n1} vs {len(xv)}")

    if FLAGGED_SNP in snp_df.index and FLAGGED_COMPOUND in cmp_df.index:
        x = snp_df.loc[FLAGGED_SNP, common_samples]
        y = cmp_df.loc[FLAGGED_COMPOUND, common_samples]
        mask = x.notna() & y.notna()
        used = [s for s, m in zip(common_samples, mask.values) if bool(m)]
        lines.append(f"[INFO] Flagged pair overlap N={len(used)}. Samples used:")
        lines.append("       " + ", ".join(used))
    else:
        lines.append("[INFO] Flagged pair not present in this run subset")

    lines.append("=" * 80)
    lines.append(f"ALL PASSED: {ok}")
    lines.append("=" * 80)

    return ok, "\n".join(lines)


if __name__ == "__main__":
    args = parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    out_dir = OUTPUT_DIR
    if args.quick_test:
        out_dir = os.path.join(OUTPUT_DIR, "quick_test")
    os.makedirs(out_dir, exist_ok=True)

    plot_dir = os.path.join(out_dir, "validation_plots")
    os.makedirs(plot_dir, exist_ok=True)

    print("=" * 80)
    print(f"Working directory: {os.getcwd()}")
    print(f"Cores: {N_CORES} / {cpu_count()}")
    print(f"MIN_SAMPLES (compute): {MIN_SAMPLES}")
    print(f"SIG_MIN_OVERLAP (significant list): {args.sig_min_overlap}")
    print(f"QMODE: {args.qmode}")
    if args.quick_test:
        print(f"QUICK TEST: {args.test_snps} SNPs × {args.test_compounds} compounds (seed={args.test_seed})")
    print("=" * 80)

    if REQUIRE_SUBSET:
        ensure_subset_file(COMPOUND_SUBSET_FILE, PCA_TSNE_FILE, EXPECTED_N_COMPOUNDS)

    if not os.path.exists(SNP_FILE):
        raise FileNotFoundError(f"Missing SNP file: {SNP_FILE}")
    if not os.path.exists(COMPOUND_FILE):
        raise FileNotFoundError(f"Missing compound file: {COMPOUND_FILE}")

    print("\nLoading SNP matrix...")
    snp_df = _read_numeric_matrix(SNP_FILE)
    print(f"SNP shape: {snp_df.shape} (rows=SNPs, cols=samples)")

    print("\nLoading compound matrix...")
    cmp_df = _read_numeric_matrix(COMPOUND_FILE)
    print(f"Compound shape (raw): {cmp_df.shape} (rows=samples, cols=compounds)")
    cmp_df = cmp_df.T
    print(f"Compound shape (transposed): {cmp_df.shape} (rows=compounds, cols=samples)")

    if COMPOUND_SUBSET_FILE is not None:
        if not os.path.exists(COMPOUND_SUBSET_FILE):
            if REQUIRE_SUBSET:
                raise FileNotFoundError(f"Missing COMPOUND_SUBSET_FILE: {COMPOUND_SUBSET_FILE}")
        keep = [str(x).strip() for x in _load_subset_ids(COMPOUND_SUBSET_FILE)]
        keep = [x for x in keep if x]
        before = cmp_df.shape[0]
        cmp_df = cmp_df.loc[cmp_df.index.intersection(pd.Index(keep))]
        after = cmp_df.shape[0]
        print(f"Compounds: {before:,} → {after:,} after subset filter")
        if EXPECTED_N_COMPOUNDS is not None and after != EXPECTED_N_COMPOUNDS:
            raise ValueError(f"Expected {EXPECTED_N_COMPOUNDS} compounds after filtering, got {after}")

    common_samples = [s for s in snp_df.columns if s in cmp_df.columns]
    if len(common_samples) == 0:
        raise ValueError("No common sample IDs between SNP and compound matrices")
    snp_df = snp_df[common_samples]
    cmp_df = cmp_df[common_samples]
    print(f"\nCommon samples: {len(common_samples)}")

    if args.quick_test:
        print("\n" + "=" * 80)
        print("QUICK TEST MODE ENABLED")
        print("=" * 80)
        snp_keep = choose_subset(snp_df.index, args.test_snps, args.test_seed)
        cmp_keep = choose_subset(cmp_df.index, args.test_compounds, args.test_seed + 1)
        if FLAGGED_SNP in snp_df.index:
            snp_keep = ensure_included(snp_keep, FLAGGED_SNP, args.test_snps)
        if FLAGGED_COMPOUND in cmp_df.index:
            cmp_keep = ensure_included(cmp_keep, FLAGGED_COMPOUND, args.test_compounds)
        snp_df = snp_df.loc[snp_keep]
        cmp_df = cmp_df.loc[cmp_keep]
        print(f"Quick test shapes: SNPs={snp_df.shape[0]:,}, Compounds={cmp_df.shape[0]:,}, Samples={snp_df.shape[1]:,}")

    snp_arr = snp_df.T.values.astype(DTYPE_FLOAT, copy=False)
    cmp_arr = cmp_df.T.values.astype(DTYPE_FLOAT, copy=False)

    n_samples, n_snps = snp_arr.shape
    _, n_compounds = cmp_arr.shape
    total_pairs = int(n_snps) * int(n_compounds)

    print("\nAligned shapes:")
    print(f"  samples × SNPs      : {snp_arr.shape}")
    print(f"  samples × compounds : {cmp_arr.shape}")
    print(f"Total pairs: {total_pairs:,}")

    mm_dir = os.path.join(out_dir, "memmap")
    os.makedirs(mm_dir, exist_ok=True)

    def mm(path, dtype, shape):
        return np.memmap(path, dtype=dtype, mode="w+", shape=shape)

    pcc_mm = mm(os.path.join(mm_dir, "pcc.dat"), np.float32, (n_snps, n_compounds))
    spr_mm = mm(os.path.join(mm_dir, "spr.dat"), np.float32, (n_snps, n_compounds))
    p_pcc_mm = mm(os.path.join(mm_dir, "p_pcc.dat"), np.float32, (n_snps, n_compounds))
    p_spr_mm = mm(os.path.join(mm_dir, "p_spr.dat"), np.float32, (n_snps, n_compounds))
    n_ov_mm = mm(os.path.join(mm_dir, "n_overlap.dat"), np.uint16, (n_snps, n_compounds))

    pcc_mm[:] = np.nan
    spr_mm[:] = np.nan
    p_pcc_mm[:] = np.nan
    p_spr_mm[:] = np.nan
    n_ov_mm[:] = 0

    print("\nComputing correlations...")
    t0 = time.time()
    indices = list(range(n_snps))

    if args.no_parallel:
        _init_worker(snp_arr, cmp_arr, MIN_SAMPLES)
        for i in tqdm(indices, desc="SNPs (serial)"):
            snp_idx, pcc_row, spr_row, pcc_p_row, spr_p_row, n_ov_row = _compute_for_snp(i)
            pcc_mm[snp_idx, :] = pcc_row
            spr_mm[snp_idx, :] = spr_row
            p_pcc_mm[snp_idx, :] = pcc_p_row
            p_spr_mm[snp_idx, :] = spr_p_row
            n_ov_mm[snp_idx, :] = n_ov_row
    else:
        ctx = get_context("spawn")
        with ctx.Pool(
            processes=N_CORES,
            initializer=_init_worker,
            initargs=(snp_arr, cmp_arr, MIN_SAMPLES),
            maxtasksperchild=500,
        ) as pool:
            for snp_idx, pcc_row, spr_row, pcc_p_row, spr_p_row, n_ov_row in tqdm(
                pool.imap(_compute_for_snp, indices, chunksize=25),
                total=len(indices),
                desc="SNPs",
            ):
                pcc_mm[snp_idx, :] = pcc_row
                spr_mm[snp_idx, :] = spr_row
                p_pcc_mm[snp_idx, :] = pcc_p_row
                p_spr_mm[snp_idx, :] = spr_p_row
                n_ov_mm[snp_idx, :] = n_ov_row

    pcc_mm.flush()
    spr_mm.flush()
    p_pcc_mm.flush()
    p_spr_mm.flush()
    n_ov_mm.flush()
    elapsed = time.time() - t0
    print(f"\n✓ Correlations done in {elapsed/60:.2f} minutes ({elapsed/3600:.2f} hours)")

    finite_pcc = int(np.isfinite(pcc_mm[:]).sum())
    finite_spr = int(np.isfinite(spr_mm[:]).sum())
    nov_all = n_ov_mm[:]
    nov_pos = nov_all[nov_all > 0]
    med_all = int(np.median(nov_all))
    if nov_pos.size > 0:
        med_pos = int(np.median(nov_pos))
    else:
        med_pos = 0
    print("\nSanity checks:")
    print(f"  Finite PCC entries      : {finite_pcc:,}")
    print(f"  Finite Spearman entries : {finite_spr:,}")
    print(f"  N_overlap median (all)  : {med_all}")
    print(f"  N_overlap median (>0)   : {med_pos}")
    print(f"  N_overlap min/max       : {int(nov_all.min())} / {int(nov_all.max())}")
    print(f"  N_overlap >= MIN_SAMPLES: {int((nov_all >= MIN_SAMPLES).sum()):,}")

    print("\nThreshold analysis (|r|>=t in BOTH PCC & Spearman)...")
    thresh_counts = {t: 0 for t in THRESHOLDS}
    valid_count = 0

    for i0 in tqdm(range(0, n_snps, ROW_CHUNK), desc="Scanning thresholds"):
        i1 = min(n_snps, i0 + ROW_CHUNK)
        pcc_chunk = pcc_mm[i0:i1, :]
        spr_chunk = spr_mm[i0:i1, :]
        valid = np.isfinite(pcc_chunk) & np.isfinite(spr_chunk)
        valid_count += int(valid.sum())
        abs_p = np.abs(pcc_chunk)
        abs_s = np.abs(spr_chunk)
        for t in THRESHOLDS:
            both = (abs_p >= t) & (abs_s >= t) & valid
            thresh_counts[t] += int(both.sum())

    thresh_rows = []
    for t in THRESHOLDS:
        c = thresh_counts[t]
        pct = (c / total_pairs) * 100.0 if total_pairs else 0.0
        thresh_rows.append({"Threshold": t, "Pairs_Passing_Both": c, "Percentage_of_All_Pairs": pct})

    thresh_df = pd.DataFrame(thresh_rows)
    thresh_out = os.path.join(out_dir, "threshold_analysis.csv.gz")
    thresh_df.to_csv(thresh_out, index=False, compression="gzip")
    print(f"✓ Saved: {thresh_out}")

    q_pcc_mm = mm(os.path.join(mm_dir, "q_pcc.dat"), np.float32, (n_snps, n_compounds))
    q_spr_mm = mm(os.path.join(mm_dir, "q_spr.dat"), np.float32, (n_snps, n_compounds))
    q_pcc_mm[:] = np.nan
    q_spr_mm[:] = np.nan

    qmode = args.qmode
    if qmode == "auto":
        if args.quick_test or total_pairs <= DEFAULT_EXACT_QVALUE_MAX_PAIRS:
            qmode = "exact"
        else:
            qmode = "hist"

    q_done = False
    if qmode == "skip":
        print("\nSkipping q-values (--qmode skip).")
    elif qmode == "exact":
        print("\nComputing q-values (EXACT)...")
        compute_qvalues_exact_from_memmap(p_pcc_mm, q_pcc_mm, "PCC")
        compute_qvalues_exact_from_memmap(p_spr_mm, q_spr_mm, "Spearman")
        q_done = True
    elif qmode == "hist":
        print("\nComputing q-values (HIST)...")
        compute_qvalues_hist_from_memmap(p_pcc_mm, q_pcc_mm, "PCC", bins=args.hist_bins)
        compute_qvalues_hist_from_memmap(p_spr_mm, q_spr_mm, "Spearman", bins=args.hist_bins)
        q_done = True
    else:
        raise ValueError(f"Unknown qmode: {qmode}")

    sig_out = os.path.join(out_dir, "significant_pairs.csv.gz")
    sig_min_overlap = int(args.sig_min_overlap)

    print(f"\nExtracting significant pairs to: {sig_out}")
    print(f"Significance rule: |PCC|>= {SIG_ABS_R}, |Spearman|>= {SIG_ABS_R}, "
          f"{'q<=0.05' if q_done else 'p<=0.05'}, N_overlap >= {sig_min_overlap}")

    sig_n = 0
    with gzip.open(sig_out, "wt", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "SNP", "Compound", "N_overlap",
                "PCC_r", "PCC_p", "PCC_q",
                "Spearman_r", "Spearman_p", "Spearman_q",
                "Direction", "Average_abs_r",
            ],
        )
        writer.writeheader()

        for i0 in tqdm(range(0, n_snps, ROW_CHUNK), desc="Scanning significant pairs"):
            i1 = min(n_snps, i0 + ROW_CHUNK)
            pcc_chunk = pcc_mm[i0:i1, :]
            spr_chunk = spr_mm[i0:i1, :]
            ppc_chunk = p_pcc_mm[i0:i1, :]
            psr_chunk = p_spr_mm[i0:i1, :]
            nov_chunk = n_ov_mm[i0:i1, :].astype(np.int32, copy=False)

            valid = np.isfinite(pcc_chunk) & np.isfinite(spr_chunk) & (nov_chunk >= sig_min_overlap)

            if q_done:
                qpc_chunk = q_pcc_mm[i0:i1, :]
                qsr_chunk = q_spr_mm[i0:i1, :]
                valid = valid & np.isfinite(qpc_chunk) & np.isfinite(qsr_chunk)
                mask = (
                    valid &
                    (np.abs(pcc_chunk) >= SIG_ABS_R) &
                    (np.abs(spr_chunk) >= SIG_ABS_R) &
                    (qpc_chunk <= Q_THRESHOLD) &
                    (qsr_chunk <= Q_THRESHOLD)
                )
            else:
                mask = (
                    valid &
                    (np.abs(pcc_chunk) >= SIG_ABS_R) &
                    (np.abs(spr_chunk) >= SIG_ABS_R) &
                    (ppc_chunk <= 0.05) &
                    (psr_chunk <= 0.05)
                )
                qpc_chunk = np.full_like(ppc_chunk, np.nan, dtype=np.float32)
                qsr_chunk = np.full_like(psr_chunk, np.nan, dtype=np.float32)

            if not mask.any():
                continue

            coords = np.argwhere(mask)
            for (ri, cj) in coords:
                snp_id = str(snp_df.index[i0 + ri])
                cmp_id = str(cmp_df.index[cj])
                pr = float(pcc_chunk[ri, cj])
                sr = float(spr_chunk[ri, cj])

                writer.writerow({
                    "SNP": snp_id,
                    "Compound": cmp_id,
                    "N_overlap": int(nov_chunk[ri, cj]),
                    "PCC_r": pr,
                    "PCC_p": float(ppc_chunk[ri, cj]),
                    "PCC_q": float(qpc_chunk[ri, cj]) if q_done else np.nan,
                    "Spearman_r": sr,
                    "Spearman_p": float(psr_chunk[ri, cj]),
                    "Spearman_q": float(qsr_chunk[ri, cj]) if q_done else np.nan,
                    "Direction": "Positive" if pr > 0 else "Negative",
                    "Average_abs_r": (abs(pr) + abs(sr)) / 2.0,
                })
                sig_n += 1

    print(f"✓ Significant pairs written: {sig_n:,}")

    print("\nSaving matrices (CSV.GZ, streamed)...")
    _save_matrix_csv_gz(pcc_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "pcc_correlation_matrix.csv.gz"))
    _save_matrix_csv_gz(spr_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "spearman_correlation_matrix.csv.gz"))
    _save_matrix_csv_gz(p_pcc_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "pcc_pvalues.csv.gz"))
    _save_matrix_csv_gz(p_spr_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "spearman_pvalues.csv.gz"))
    if q_done:
        _save_matrix_csv_gz(q_pcc_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "pcc_qvalues.csv.gz"))
        _save_matrix_csv_gz(q_spr_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "spearman_qvalues.csv.gz"))
    _save_matrix_uint_csv_gz(n_ov_mm, snp_df.index, cmp_df.index, os.path.join(out_dir, "n_overlap_matrix.csv.gz"))

    summary = {
        "Total_SNPs": int(n_snps),
        "Total_Compounds": int(n_compounds),
        "Total_Pairs": int(total_pairs),
        "Common_Samples": int(len(common_samples)),
        "Valid_Correlations": int(valid_count),
        "CPU_Cores_Used": int(N_CORES),
        "Min_Samples_Per_Correlation": int(MIN_SAMPLES),
        "Sig_abs_r": float(SIG_ABS_R),
        "Sig_q_threshold": float(Q_THRESHOLD) if q_done else None,
        "Sig_min_overlap": int(sig_min_overlap),
        "Significant_Pairs": int(sig_n),
        "Computation_Time_Minutes": float(elapsed / 60.0),
        "Quick_Test_Mode": bool(args.quick_test),
        "Q_Mode_Used": qmode,
        "Histogram_Bins": int(args.hist_bins) if qmode == "hist" else None,
    }
    pd.DataFrame([summary]).to_csv(os.path.join(out_dir, "summary_statistics.csv.gz"), index=False, compression="gzip")

    print("\nGenerating validation scatter plots...")
    flagged_plot = os.path.join(plot_dir, f"scatter_FLAGGED__{FLAGGED_SNP}__{FLAGGED_COMPOUND}.png")
    if FLAGGED_SNP in snp_df.index and FLAGGED_COMPOUND in cmp_df.index:
        stats = plot_pair_scatter(FLAGGED_SNP, FLAGGED_COMPOUND, snp_df, cmp_df, common_samples, flagged_plot)
        print(f"✓ Flagged scatter saved: {os.path.basename(flagged_plot)} | N={int(stats['N'])}")
    else:
        print("⚠ Flagged pair not present in this run subset")

    flagged_debug_path = os.path.join(out_dir, "flagged_pair_debug.csv")
    if FLAGGED_SNP in snp_df.index and FLAGGED_COMPOUND in cmp_df.index:
        x = snp_df.loc[FLAGGED_SNP, common_samples]
        y = cmp_df.loc[FLAGGED_COMPOUND, common_samples]
        mask = x.notna() & y.notna()
        dbg = pd.DataFrame({
            "Sample": [s for s, m in zip(common_samples, mask.values) if bool(m)],
            "SNP_value": x[mask].astype(float).values,
            "Compound_value": y[mask].astype(float).values,
        })
        dbg.to_csv(flagged_debug_path, index=False)
        print(f"✓ Flagged pair debug saved: {os.path.basename(flagged_debug_path)} (N={len(dbg)})")

    if not args.no_tests:
        print("\nRunning post-run tests...")
        all_ok, report = post_run_tests(
            out_dir=out_dir,
            snp_df=snp_df,
            cmp_df=cmp_df,
            common_samples=common_samples,
            subset_path=COMPOUND_SUBSET_FILE,
            pca_tsne_path=PCA_TSNE_FILE,
            pcc_mm=pcc_mm,
            spr_mm=spr_mm,
            p_pcc_mm=p_pcc_mm,
            p_spr_mm=p_spr_mm,
            n_ov_mm=n_ov_mm,
            sig_pairs_path=sig_out,
            args=args,
        )
        print(report)
        report_path = os.path.join(out_dir, "post_run_test_report.txt")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"✓ Test report saved to: {report_path}")
        if args.strict_tests and not all_ok:
            raise SystemExit(2)

    print("\n" + "=" * 80)
    print("✓ DONE — outputs generated.")
    print(f"Output folder: {out_dir}/")
    if args.quick_test:
        print("Quick test finished. Rerun without --quick-test for the full run.")
    print("=" * 80)
