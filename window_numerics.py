"""
Are the "undecided windows" a real phenomenon or an artefact of the proof?

Everything is exact: the state space M(s_0) = {gamma in {0,1}^d : |gamma| <= s_0}
is enumerated, the tempered targets pi_lambda are evaluated at every state, the
flip-and-exchange kernel P_lambda is assembled, and its spectral gap is obtained
by diagonalisation.  No Monte Carlo.

Four questions, one per panel of the figure:

  (A) Is Gap(P_lambda) monotone in lambda, and does it dip where k_lambda > 0?
      -> if the true gap is flat across a window, the window is an artefact.

  (B) k_lambda, Theta_lambda and Delta^C_lambda along the path.

  (C) Are the block partition functions c(J) = pibar_lambda(J) / pi_lambda(v_J)
      nearly constant in J?  This is the load-bearing guess behind the claim
      that the {2(2s*+1)}^6 in Proposition 7.1 is removable.  We report
      spread(lambda) = max_J c(J) / min_J c(J).

  (D) The exact hypercube congestion varrho_lambda of eq:rholam, against the
      bound {2(2s*+1)}^6 exp{6 k_lambda Delta^C_lambda} of Proposition 7.1,
      and against the congestion of the same ensemble when segment (ii) is
      allowed to use the EXCHANGE move of eq:neigh instead of flipping the
      coordinates of J xor J' one at a time.

Outputs
-------
window_numerics.pdf   four-panel figure
window_numerics.txt   the numbers

Run:  python3 window_numerics.py
"""

import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG_SEED = 20240824

N_OBS = 100          # n
D = 14               # p
S0 = 4               # s_0
S_STAR = 3           # s* = |S|, influential set S = {0,1,2}
SIGMA0 = 1.0
KAPPA = 2.0
ALPHA = 1.5
G_PRIOR = D ** (2 * ALPHA)

S_SET = tuple(range(S_STAR))
SC_SET = tuple(range(S_STAR, D))
THR = KAPPA * np.log(D)              # kappa log p
MARGIN = np.log(2 * S_STAR)          # the log(2s*) of def:cpm

# three separated coefficient magnitudes -> three separated windows
BETA_S = np.array([1.20, 0.70, 0.45])


# --------------------------------------------------------------------------
# design, data, marginal likelihood
# --------------------------------------------------------------------------

def build_design(rng, rho):
    """Equicorrelated influential columns; null columns residualised against
    span(X_S) so that Assumption 6.1 is comfortable.  ||X_j||^2 = n."""
    Z = rng.standard_normal((N_OBS, D))
    Sigma = (1 - rho) * np.eye(S_STAR) + rho * np.ones((S_STAR, S_STAR))
    L = np.linalg.cholesky(Sigma)
    X = np.empty((N_OBS, D))
    X[:, :S_STAR] = Z[:, :S_STAR] @ L.T
    Q, _ = np.linalg.qr(X[:, :S_STAR])
    Xnull = Z[:, S_STAR:]
    X[:, S_STAR:] = Xnull - Q @ (Q.T @ Xnull)
    X /= np.linalg.norm(X, axis=0, keepdims=True) / np.sqrt(N_OBS)
    return X


def enumerate_states(d, s0):
    out = []
    for k in range(s0 + 1):
        out.extend(itertools.combinations(range(d), k))
    return out


def log_marginal(states, X, y):
    yy = float(y @ y)
    out = np.empty(len(states))
    for i, x in enumerate(states):
        if len(x) == 0:
            r2 = 0.0
        else:
            Q, _ = np.linalg.qr(X[:, list(x)])
            r2 = float(np.sum((Q.T @ y) ** 2)) / yy
        out[i] = (-0.5 * len(x) * np.log1p(G_PRIOR)
                  - 0.5 * N_OBS * np.log1p(G_PRIOR * (1.0 - r2)))
    return out


# --------------------------------------------------------------------------
# kernel and exact spectral gap  (as in binary_numerics.py)
# --------------------------------------------------------------------------

def neighbour_structure(states):
    index = {s: i for i, s in enumerate(states)}
    rows, cols, prop = [], [], []
    for i, x in enumerate(states):
        xs, k = set(x), len(x)
        for j in range(D):
            if j in xs:
                y = tuple(sorted(xs - {j}))
            else:
                if k + 1 > S0:
                    continue
                y = tuple(sorted(xs | {j}))
            rows.append(i); cols.append(index[y]); prop.append(1.0 / (4 * D))
        if 0 < k < D:
            p = 1.0 / (4.0 * k * (D - k))
            for a in x:
                for b in range(D):
                    if b in xs:
                        continue
                    y = tuple(sorted((xs - {a}) | {b}))
                    rows.append(i); cols.append(index[y]); prop.append(p)
    return np.asarray(rows), np.asarray(cols), np.asarray(prop)


def spectral_gap(logpi, nb):
    rows, cols, prop = nb
    m = logpi.shape[0]
    delta = logpi[cols] - logpi[rows]
    off_sym = prop * np.exp(-0.5 * np.abs(delta))
    off_k = prop * np.exp(np.minimum(delta, 0.0))
    diag = 1.0 - np.bincount(rows, weights=off_k, minlength=m)
    Sym = sp.coo_matrix((off_sym, (rows, cols)), shape=(m, m)).tocsr()
    Sym = Sym + sp.diags(diag)
    v = np.exp(0.5 * (logpi - logpi.max()))
    v /= np.linalg.norm(v)
    lin = spla.LinearOperator(
        (m, m), matvec=lambda z: Sym @ z - v * float(v @ z), dtype=float)
    lam2 = float(spla.eigsh(lin, k=1, which="LA", tol=1e-11,
                            maxiter=50000)[0][0])
    return 1.0 - lam2


# --------------------------------------------------------------------------
# quantities of the theory
# --------------------------------------------------------------------------

def psi_extremes(logm, index):
    """psi_j^- , psi_j^+ : min / max of psi_j(v) over v subset of S \\ {j}."""
    lo, hi = {}, {}
    for j in S_SET:
        rest = [i for i in S_SET if i != j]
        vals = []
        for r in range(len(rest) + 1):
            for ctx in itertools.combinations(rest, r):
                a = index[tuple(sorted(ctx))]
                b = index[tuple(sorted(ctx + (j,)))]
                vals.append(logm[b] - logm[a])
        lo[j], hi[j] = min(vals), max(vals)
    return lo, hi


def decided_sets(lam, lo, hi):
    """c^-_lambda, c^+_lambda, C_lambda of def:cpm."""
    cminus = tuple(j for j in S_SET if lam * lo[j] >= THR + MARGIN)
    cplus = tuple(j for j in S_SET if lam * hi[j] > THR - MARGIN)
    C = tuple(j for j in cplus if j not in cminus)
    return cminus, cplus, C


def delta_C(lam, logm, index, cminus, C):
    """Delta^C_lambda of eq:DeltaC: the spread of psi_j over contexts
    c^-_lambda u U with U subset of C_lambda \\ {j}."""
    best = 0.0
    for j in C:
        others = [i for i in C if i != j]
        vals = []
        for r in range(len(others) + 1):
            for U in itertools.combinations(others, r):
                ctx = tuple(sorted(set(cminus) | set(U)))
                vals.append(logm[index[tuple(sorted(set(ctx) | {j}))]]
                            - logm[index[ctx]])
        best = max(best, max(vals) - min(vals))
    return best


def block_measure(states, logpi, C):
    """pibar_lambda(J) = pi_lambda(B_J), B_J = {gamma : gamma n C = J}."""
    Cs = set(C)
    keys = [tuple(sorted(Cs & set(x))) for x in states]
    subsets = [tuple(sorted(J))
               for r in range(len(C) + 1)
               for J in itertools.combinations(sorted(C), r)]
    pos = {J: i for i, J in enumerate(subsets)}
    w = np.exp(logpi - logpi.max())
    bar = np.zeros(len(subsets))
    for k, wi in zip(keys, w):
        bar[pos[k]] += wi
    bar /= bar.sum()
    return subsets, bar


def swap_congestion(subsets, bar, C):
    """As hypercube_congestion, but segment (ii) pairs off J1 \\ J2 with
    J2 \\ J1 using exchange moves before flipping the remainder.  The chain
    already has this move; eq:edge covers it with 1/(4 p s_0) instead of
    1/(4 p)."""
    if not C:
        return 1.0
    pos = {J: i for i, J in enumerate(subsets)}
    edge = {}
    for J1 in subsets:
        for J2 in subsets:
            if J1 == J2:
                continue
            A = sorted(set(J1) - set(J2)); B = sorted(set(J2) - set(J1))
            m = min(len(A), len(B))
            steps = [("swap", a, b) for a, b in zip(A[:m], B[:m])]
            steps += [("flip", a, None) for a in A[m:]]
            steps += [("flip", b, None) for b in B[m:]]
            cur = set(J1)
            for kind, a, b in steps:
                nxt = (cur - {a}) | {b} if kind == "swap" else cur ^ {a}
                e = (tuple(sorted(cur)), tuple(sorted(nxt)))
                edge[e] = edge.get(e, 0.0) + bar[pos[J1]] * bar[pos[J2]]
                cur = nxt
    return max(f / min(bar[pos[a]], bar[pos[b]])
               for (a, b), f in edge.items())


def hypercube_congestion(subsets, bar, C):
    """varrho_lambda of eq:rholam, exactly: paths flip J1 x-or J2 in
    increasing coordinate order."""
    if not C:
        return 1.0
    pos = {J: i for i, J in enumerate(subsets)}
    edge = {}
    for J1 in subsets:
        for J2 in subsets:
            if J1 == J2:
                continue
            cur = set(J1)
            for j in sorted(set(J1) ^ set(J2)):
                nxt = cur ^ {j}
                e = (tuple(sorted(cur)), tuple(sorted(nxt)))
                edge[e] = edge.get(e, 0.0) + bar[pos[J1]] * bar[pos[J2]]
                cur = nxt
    worst = 0.0
    for (a, b), flow in edge.items():
        worst = max(worst, flow / min(bar[pos[a]], bar[pos[b]]))
    return worst


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def run(rho, states, index, nb, X, y, log):
    logm = log_marginal(states, X, y)
    sizes = np.array([len(x) for x in states], dtype=float)
    logprior = -KAPPA * np.log(D) * sizes
    lo, hi = psi_extremes(logm, index)

    grid = np.linspace(1e-3, 1.0, 240)
    rec = {k: [] for k in ("gap", "k", "theta", "dC", "spread", "rho_ex",
                           "rho_bd", "rho_sw", "barmin")}
    for lam in grid:
        logpi = lam * logm + logprior
        logpi -= logpi.max()
        rec["gap"].append(spectral_gap(logpi, nb))

        cminus, cplus, C = decided_sets(lam, lo, hi)
        rec["k"].append(len(C))
        rec["theta"].append(max([np.exp(lam * hi[j] - THR) for j in C],
                                default=1.0))
        dC = delta_C(lam, logm, index, cminus, C) if C else 0.0
        rec["dC"].append(dC)

        subsets, bar = block_measure(states, logpi, C)
        # c(J) = pibar(J) / pi(v_J),  v_J = c^- u J
        pin = np.exp(logpi - logpi.max()); pin /= pin.sum()
        cJ = []
        for J in subsets:
            vJ = tuple(sorted(set(cminus) | set(J)))
            if len(vJ) <= S0:
                cJ.append(bar[subsets.index(J)] / pin[index[vJ]])
        rec["spread"].append(max(cJ) / min(cJ) if cJ else 1.0)

        rec["rho_ex"].append(hypercube_congestion(subsets, bar, C))
        rec["rho_sw"].append(swap_congestion(subsets, bar, C))
        rec["barmin"].append(bar.min())
        rec["rho_bd"].append((2 * (2 * S_STAR + 1)) ** 6
                             * np.exp(6 * len(C) * dC))

    for k in rec:
        rec[k] = np.asarray(rec[k])
    rec["grid"] = grid

    win = rec["k"] > 0
    log.append(f"\n=== rho = {rho} ===")
    log.append(f"psi_j^-  : {[round(lo[j],2) for j in S_SET]}")
    log.append(f"psi_j^+  : {[round(hi[j],2) for j in S_SET]}")
    log.append(f"kappa log p = {THR:.3f},  margin log(2s*) = {MARGIN:.3f}")
    log.append(f"measure of {{k_lambda > 0}} = {win.mean():.3f}")
    log.append(f"max k_lambda = {rec['k'].max()},  "
               f"max Delta^C = {rec['dC'].max():.4f}")
    log.append(f"gap: min {rec['gap'].min():.3e} at lambda="
               f"{grid[rec['gap'].argmin()]:.3f}, "
               f"max {rec['gap'].max():.3e}, "
               f"ratio {rec['gap'].max()/rec['gap'].min():.2f}")
    if win.any():
        log.append(f"gap inside windows : min {rec['gap'][win].min():.3e}, "
                   f"mean {rec['gap'][win].mean():.3e}")
        log.append(f"gap outside windows: min {rec['gap'][~win].min():.3e}, "
                   f"mean {rec['gap'][~win].mean():.3e}")
    d = np.diff(rec["gap"])
    log.append(f"monotone in lambda? increasing frac = "
               f"{(d > 0).mean():.3f}  (1.0 = increasing, 0.0 = decreasing)")
    log.append(f"c(J) spread max_J/min_J : max over lambda = "
               f"{rec['spread'].max():.4f}")
    log.append(f"varrho exact, flip-only paths : max = {rec['rho_ex'].max():.4g}")
    log.append(f"varrho exact, exchange allowed : max = {rec['rho_sw'].max():.4g}")
    log.append(f"  ratio flip/exchange at worst lambda : "
               f"{(rec['rho_ex']/rec['rho_sw']).max():.3g}")
    log.append(f"lightest block min_J pibar(J) : {rec['barmin'].min():.3e}")
    log.append(f"varrho bound (Prop 7.1) : max = {rec['rho_bd'].max():.3e}"
               f"   [ (2(2s*+1))^6 = {(2*(2*S_STAR+1))**6:.3e} ]")
    log.append(f"bound / exact at the worst lambda : "
               f"{(rec['rho_bd']/np.maximum(rec['rho_ex'],1)).max():.3e}")
    return rec


def main():
    noise = np.random.default_rng(RNG_SEED + 7).standard_normal(N_OBS)
    states = enumerate_states(D, S0)
    index = {s: i for i, s in enumerate(states)}
    nb = neighbour_structure(states)
    log = [f"states |M(s_0)| = {len(states)}, n = {N_OBS}, p = {D}, "
           f"s_0 = {S0}, s* = {S_STAR}, kappa = {KAPPA}, alpha = {ALPHA}",
           f"beta*_S = {BETA_S}"]

    results = {}
    for rho in (0.0, 0.5, 0.85, 0.95):
        X = build_design(np.random.default_rng(RNG_SEED), rho)
        y = X[:, :S_STAR] @ BETA_S + SIGMA0 * noise
        results[rho] = run(rho, states, index, nb, X, y, log)

    txt = "\n".join(log)
    open("window_numerics.txt", "w").write(txt + "\n")
    print(txt)

    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    cols = {0.0: "tab:blue", 0.5: "tab:orange", 0.85: "tab:red",
            0.95: "tab:purple"}
    for rho, r in results.items():
        g = r["grid"]
        ax[0, 0].semilogy(g, r["gap"], color=cols[rho], label=f"$\\rho={rho}$")
        ax[0, 1].plot(g, r["dC"], color=cols[rho])
        ax[0, 1].plot(g, r["k"] / 10.0, color=cols[rho], ls=":", lw=1)
        ax[1, 0].plot(g, r["spread"], color=cols[rho])
        ax[1, 1].semilogy(g, np.maximum(r["rho_ex"], 1e-3), color=cols[rho])
        ax[1, 1].semilogy(g, np.maximum(r["rho_sw"], 1e-3), color=cols[rho],
                          ls="--", lw=1.4)
        for a in ax.ravel():
            a.fill_between(g, *a.get_ylim(), where=r["k"] > 0,
                           color=cols[rho], alpha=0.07, lw=0)

    ax[0, 0].set_title("(A) exact $\\mathrm{Gap}(\\mathbf{P}_\\lambda)$; "
                       "shading = $k_\\lambda>0$")
    ax[0, 0].set_xlabel("$\\lambda$"); ax[0, 0].legend(fontsize=8)
    ax[0, 1].set_title("(B) $\\Delta^C_\\lambda$ (solid), "
                       "$k_\\lambda/10$ (dotted)")
    ax[0, 1].set_xlabel("$\\lambda$")
    ax[1, 0].set_title("(C) block ratio spread "
                       "$\\max_J c(J)/\\min_J c(J)$")
    ax[1, 0].set_xlabel("$\\lambda$")
    ax[1, 1].set_title("(D) exact $\\varrho_\\lambda$: flip-only paths "
                       "(solid) vs exchange allowed (dashed)")
    ax[1, 1].set_xlabel("$\\lambda$")
    fig.tight_layout()
    fig.savefig("window_numerics.pdf")
    print("\nwrote window_numerics.pdf, window_numerics.txt")


if __name__ == "__main__":
    main()
