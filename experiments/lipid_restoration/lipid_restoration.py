#!/usr/bin/env python3
"""Lipid representability-restoration experiment (AGILE / LANTERN panel).

Compares a lossy encoder against a lossless one on held-out interaction contrasts:
  radius-2 Morgan count fingerprints  -> merges 1100 states into 920 classes,
                                         rank(CH)=740, encoding floor 0.3150033
  radius-8 Morgan count fingerprints  -> separates all 1100 states, floor exactly 0

Training rows are identical in both arms, so training support is held constant by
construction and only representability changes. Evaluation is 5-fold grouped by
head group A, so each fold predicts tail-by-tail interactions of unseen head groups.

Self-check runs first and ABORTS unless it reproduces the published diagnostic:
920 classes, rank(CH)=740, floor 0.3150033148, recorded contrast energy 12.425324.

Usage:
    python3 lipid_restoration.py /path/to/lantern.csv [--seeds 20] [--max-iter 3000]

Needs: rdkit, numpy, scikit-learn. The published diagnostic used rdkit 2026.03.1.
"""
import argparse, csv, json, re, sys, warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as fpg
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

EXPECT = dict(classes=920, rank=740, floor=0.3150033148, energy=12.425324)

def say(*a):
    print(*a); sys.stdout.flush()

def load(path):
    rows = list(csv.DictReader(open(path)))
    pat = re.compile(r'^A(\d+)C(\d+)B(\d+)$')
    rec = {}
    for r in rows:
        m = pat.match(r['Lipid'])
        if not m: sys.exit('unparsable lipid id: %r' % r['Lipid'])
        rec[tuple(int(x) for x in m.groups())] = (r['SMILES'], float(r['Target']))
    return rec

def contrasts(rec):
    """C x B interaction contrasts within each A level, anchored at the first C and B."""
    A = sorted({k[0] for k in rec}); C = sorted({k[1] for k in rec}); B = sorted({k[2] for k in rec})
    if len(rec) != len(A)*len(C)*len(B): sys.exit('panel is not a complete factorial')
    keys = sorted(rec); idx = {k: i for i, k in enumerate(keys)}
    c0, b0 = C[0], B[0]
    M = np.zeros((len(A)*(len(C)-1)*(len(B)-1), len(keys)))
    r = 0
    for a in A:
        for c in C[1:]:
            for b in B[1:]:
                for k, s in (((a,c,b),1), ((a,c,b0),-1), ((a,c0,b),-1), ((a,c0,b0),1)):
                    M[r, idx[k]] += s
                r += 1
    return keys, M, (len(A), len(C), len(B))

def fingerprints(rec, keys, radius, nbits=2048):
    gen = fpg.GetMorganGenerator(radius=radius, fpSize=nbits)
    out = []
    for k in keys:
        mol = Chem.MolFromSmiles(rec[k][0])
        if mol is None: sys.exit('rdkit could not parse %s' % (k,))
        out.append(np.asarray(gen.GetCountFingerprintAsNumPy(mol), dtype=np.float64))
    return np.vstack(out)

def classes_of(X):
    u = {}
    for i, v in enumerate(map(tuple, X)): u.setdefault(v, []).append(i)
    return u

def floor_of(X, Cm, I, q):
    u = classes_of(X)
    H = np.zeros((X.shape[0], len(u)))
    for j, ix in enumerate(u.values()):
        for i in ix: H[i, j] = 1.0
    M = Cm @ H
    res = I - M @ np.linalg.pinv(M) @ I
    return float(res @ res) / q, int(np.linalg.matrix_rank(M)), len(u)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('panel'); ap.add_argument('--seeds', type=int, default=20)
    ap.add_argument('--max-iter', type=int, default=3000)
    ap.add_argument('--hidden', type=int, default=64)
    ap.add_argument('--folds', type=int, default=5)
    ap.add_argument('--out', default='lipid_restoration_results.json')
    a = ap.parse_args()

    rec = load(a.panel)
    keys, Cm, (nA, nC, nB) = contrasts(rec)
    q = Cm.shape[0]
    y = np.array([rec[k][1] for k in keys])
    I = Cm @ y
    energy = float(I @ I) / q
    say('panel %s : states=%d  A=%d C=%d B=%d  q=%d' % (a.panel, len(rec), nA, nC, nB, q))

    X2 = fingerprints(rec, keys, 2); X8 = fingerprints(rec, keys, 8)
    f2, r2, k2 = floor_of(X2, Cm, I, q)
    f8, r8, k8 = floor_of(X8, Cm, I, q)
    say('radius 2 : classes=%d  rank(CH)=%d  floor=%.10f  (%.4f%% of energy)'
        % (k2, r2, f2, 100*f2/energy))
    say('radius 8 : classes=%d  rank(CH)=%d  floor=%.10f' % (k8, r8, f8))
    say('recorded contrast energy = %.6f' % energy)

    bad = []
    if k2 != EXPECT['classes']: bad.append('classes %d != %d' % (k2, EXPECT['classes']))
    if r2 != EXPECT['rank']:    bad.append('rank %d != %d' % (r2, EXPECT['rank']))
    if abs(f2 - EXPECT['floor']) > 5e-9:  bad.append('floor %.10f != %.10f' % (f2, EXPECT['floor']))
    if abs(energy - EXPECT['energy']) > 5e-6: bad.append('energy %.6f != %.6f' % (energy, EXPECT['energy']))
    if bad:
        say('\nSELF-CHECK FAILED: ' + '; '.join(bad))
        say('Setup does not match the published diagnostic. Stopping; do not use these fits.')
        sys.exit(1)
    say('self-check OK: reproduces the published diagnostic exactly\n')
    if a.seeds == 0:
        say('--seeds 0, stopping after the self-check'); return

    def cmse(pred):
        d = I - Cm @ pred
        return float(d @ d) / q

    Aof = np.array([k[0] for k in keys]); Alev = sorted(set(Aof.tolist()))
    Z2 = StandardScaler().fit_transform(X2); Z8 = StandardScaler().fit_transform(X8)

    say('=== ridge reference (checks whether the fitted contrast vector is ~0) ===')
    for nm, Z in (('radius 2', Z2), ('radius 8', Z8)):
        p = Ridge(alpha=10.0).fit(Z, y).predict(Z)
        say('  %s in-sample contrast MSE %.6f   (= energy %.6f means zero fitted contrasts)'
            % (nm, cmse(p), energy))

    capped = [0]
    def heldout(Z, seed):
        rng = np.random.default_rng(seed)
        pred = np.zeros_like(y)
        for ch in np.array_split(rng.permutation(Alev), a.folds):
            te = np.isin(Aof, ch)
            m = MLPRegressor(hidden_layer_sizes=(a.hidden,), max_iter=a.max_iter,
                             random_state=seed, early_stopping=False, tol=1e-7,
                             n_iter_no_change=25, learning_rate_init=3e-3)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', ConvergenceWarning)
                m.fit(Z[~te], y[~te])
            if m.n_iter_ >= a.max_iter: capped[0] += 1
            pred[te] = m.predict(Z[te])
        return cmse(pred)

    say('\n=== held-out MLP, %d-fold grouped by head group, %d paired seeds ==='
        % (a.folds, a.seeds))
    h2, h8 = [], []
    for s in range(a.seeds):
        x, z = heldout(Z2, s), heldout(Z8, s)
        h2.append(x); h8.append(z)
        say('  seed %3d : radius2 %.6f   radius8 %.6f   diff %+.6f' % (s, x, z, x - z))
    h2 = np.array(h2); h8 = np.array(h8); d = h2 - h8
    say('\n  radius2 mean %.6f   radius8 mean %.6f' % (h2.mean(), h8.mean()))
    say('  paired diff (radius2 - radius8): mean %+.6f  min %+.6f  max %+.6f'
        % (d.mean(), d.min(), d.max()))
    say('  radius8 better in %d of %d seeds' % (int((d > 0).sum()), a.seeds))
    rng = np.random.default_rng(0)
    bs = np.array([rng.choice(d, size=a.seeds, replace=True).mean() for _ in range(20000)])
    lo, hi = np.percentile(bs, 2.5), np.percentile(bs, 97.5)
    say('  bootstrap 95%% CI on mean diff: [%+.6f, %+.6f]  %s'
        % (lo, hi, 'excludes zero' if lo * hi > 0 else 'includes zero'))
    say('  fits that hit the iteration cap: %d of %d  %s'
        % (capped[0], 2 * a.seeds * a.folds,
           '<-- RAISE --max-iter, results not safe to quote' if capped[0] else '(all converged)'))
    json.dump(dict(panel=a.panel, q=q, energy=energy, floor_r2=f2, floor_r8=f8,
                   rank_r2=r2, classes_r2=k2, seeds=a.seeds, max_iter=a.max_iter,
                   hidden=a.hidden, folds=a.folds, heldout_r2=h2.tolist(),
                   heldout_r8=h8.tolist(), paired_mean=float(d.mean()),
                   ci=[float(lo), float(hi)], capped_fits=capped[0]),
              open(a.out, 'w'), indent=1)
    say('\nwrote %s' % a.out)

if __name__ == '__main__':
    main()
