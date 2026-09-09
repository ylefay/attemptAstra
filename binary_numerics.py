"""
Exact computations for the numerical illustration of the Bayesian
variable-selection section (Section `sec:binary`).

Everything is exact: the state space X = {x in {0,1}^d : |x| <= s_0} is
enumerated, the tempered targets pi_lambda are evaluated at every state, the
flip-and-exchange kernel K_lambda is assembled, and its spectral gap is
obtained by diagonalisation.  No Monte Carlo is involved.

Outputs
-------
binary_numerics.pdf   three-panel figure
binary_numerics.txt   the numbers quoted in the text

Run:  python3 binary_numerics.py
"""

import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG_SEED = 20240824

# --------------------------------------------------------------------------
# 1. Design and data
# --------------------------------------------------------------------------

N_OBS = 100          # n
D = 14               # d
S0 = 4               # s_0, truncation level
S_STAR = 3           # s* = |S|; influential set is S = {0, 1, 2}
RHO = 0.3            # correlation among the influential columns
SIGMA0 = 1.0
KAPPA = 2.0          # sparsity prior exponent
ALPHA = 1.5          # g = d^{2 alpha}
G_PRIOR = D ** (2 * ALPHA)

S_SET = tuple(range(S_STAR))
SC_SET = tuple(range(S_STAR, D))


def build_design(rng):
    """Columns of S are equicorrelated; columns of S^c are orthogonal to
    span(X_S), so that the proxy energy Xi of eq:binary_proxy vanishes.
    All columns are normalised to ||X_j||^2 = n."""
    Z = rng.standard_normal((N_OBS, D))
    Sigma = (1 - RHO) * np.eye(S_STAR) + RHO * np.ones((S_STAR, S_STAR))
    L = np.linalg.cholesky(Sigma)
    X = np.empty((N_OBS, D))
    X[:, :S_STAR] = Z[:, :S_STAR] @ L.T
    # residualise the null columns against span(X_S)
    Q, _ = np.linalg.qr(X[:, :S_STAR])
    Xnull = Z[:, S_STAR:]
    X[:, S_STAR:] = Xnull - Q @ (Q.T @ Xnull)
    X /= np.linalg.norm(X, axis=0, keepdims=True) / np.sqrt(N_OBS)
    return X


# --------------------------------------------------------------------------
# 2. State space and marginal likelihood
# --------------------------------------------------------------------------

def enumerate_states(d, s0):
    states = []
    for k in range(s0 + 1):
        states.extend(itertools.combinations(range(d), k))
    return states


def log_marginal(states, X, y):
    """log m(x) of eq:binary_m, up to an additive constant."""
    yy = float(y @ y)
    out = np.empty(len(states))
    for i, x in enumerate(states):
        if len(x) == 0:
            r2 = 0.0
        else:
            Xx = X[:, list(x)]
            Q, _ = np.linalg.qr(Xx)
            r2 = float(np.sum((Q.T @ y) ** 2)) / yy
        out[i] = (-0.5 * len(x) * np.log1p(G_PRIOR)
                  - 0.5 * N_OBS * np.log1p(G_PRIOR * (1.0 - r2)))
    return out


# --------------------------------------------------------------------------
# 3. Kernel and spectral gap
# --------------------------------------------------------------------------

def neighbour_structure(states):
    """Return, for each ordered pair (x, x') with x' a *proposal* from x
    inside X, the row/col index and the proposal probability of the lazy
    flip-and-exchange kernel K1-K4.  The proposal law is symmetric, so a
    single array of probabilities suffices."""
    index = {s: i for i, s in enumerate(states)}
    rows, cols, prop = [], [], []
    for i, x in enumerate(states):
        xs = set(x)
        k = len(x)
        # (K2) single flip: probability 1/2 * 1/2 * 1/d
        for j in range(D):
            if j in xs:
                y = tuple(sorted(xs - {j}))
            else:
                if k + 1 > S0:
                    continue          # proposal outside X, rejected
                y = tuple(sorted(xs | {j}))
            rows.append(i); cols.append(index[y]); prop.append(1.0 / (4 * D))
        # (K3) exchange: probability 1/2 * 1/2 * 1/(|x| (d - |x|))
        if 0 < k < D:
            p = 1.0 / (4.0 * k * (D - k))
            for a in x:
                for b in range(D):
                    if b in xs:
                        continue
                    y = tuple(sorted((xs - {a}) | {b}))
                    rows.append(i); cols.append(index[y]); prop.append(p)
    return (np.asarray(rows), np.asarray(cols), np.asarray(prop))


def spectral_gap(logpi, nb):
    """Gap(K_lambda) = 1 - lambda_2, computed from the symmetrised kernel

        Sym(x, y) = prop(x, y) exp(-|log pi(x) - log pi(y)| / 2),

    which is similar to K and numerically stable however small pi gets."""
    rows, cols, prop = nb
    m = logpi.shape[0]
    delta = logpi[cols] - logpi[rows]
    off_sym = prop * np.exp(-0.5 * np.abs(delta))
    off_k = prop * np.exp(np.minimum(delta, 0.0))       # true K(x, y)
    diag = 1.0 - np.bincount(rows, weights=off_k, minlength=m)
    Sym = sp.coo_matrix((off_sym, (rows, cols)), shape=(m, m)).tocsr()
    Sym = Sym + sp.diags(diag)
    # deflate the top eigenpair, whose eigenvector is sqrt(pi)
    v = np.exp(0.5 * (logpi - logpi.max()))
    v /= np.linalg.norm(v)
    lin = spla.LinearOperator(
        (m, m), matvec=lambda z: Sym @ z - v * float(v @ z), dtype=float)
    lam2 = float(spla.eigsh(lin, k=1, which="LA",
                            tol=1e-10, maxiter=20000)[0][0])
    return 1.0 - lam2


# --------------------------------------------------------------------------
# 4. Quantities of the theory
# --------------------------------------------------------------------------

def psi_bounds(logm, index):
    """psi_j^- and psi_j^+ of eq:binary_spread: min and max over x subset of
    S \\ {j} of log{m(x u {j}) / m(x)}."""
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


def null_exponent(logm, states, index):
    """The smallest a_0 for which Assumption ass:binary_null holds, namely
    a_0 = alpha + max psi_j(x) / log d over x in X with |x| < s_0 and
    j in S^c \\ x.  Also returns the sparsity-budget slack
    min(kappa, kappa + alpha - a_0) - 1 of eq:binary_budget_endpoints."""
    best = -np.inf
    for i, x in enumerate(states):
        if len(x) >= S0:
            continue
        xs = set(x)
        for j in SC_SET:
            if j in xs:
                continue
            best = max(best, logm[index[tuple(sorted(xs | {j}))]] - logm[i])
    a0 = ALPHA + best / np.log(D)
    slack = min(KAPPA, KAPPA + ALPHA - a0) - 1.0
    lam_c = 1.0 if a0 <= ALPHA else min(1.0, (KAPPA - 1.0) / (a0 - ALPHA))
    return a0, slack, lam_c


def margin_failures(lam, lo, hi):
    """C_lambda of def:binary_margin and its cardinality k_lambda."""
    thr = KAPPA * np.log(D)
    marg = np.log(2 * S_STAR)
    bad = [j for j in S_SET
           if not (lam * lo[j] >= thr + marg or lam * hi[j] <= thr - marg)]
    return bad


def variance_function(logm, logprior, lam):
    """V(lambda) = Var_{pi_lambda}[log m] of eq:binary_V."""
    w = lam * logm + logprior
    w -= w.max()
    p = np.exp(w); p /= p.sum()
    mean = float(p @ logm)
    return float(p @ (logm - mean) ** 2)


def schedule(logm, logprior, grid):
    """The schedule of prop:binary_horizon(ii).  Delta_t is the largest
    admissible step, found by bisection; the constraint
    Delta^2 max{V(l) : l in [cur, cur + 2 Delta]} <= log 2 is evaluated by
    sampling V on the look-ahead interval."""
    V = np.array([variance_function(logm, logprior, l) for l in grid])

    def admissible(cur, step, npts=33):
        pts = np.linspace(cur, min(cur + 2 * step, 2.0), npts)
        vmax = max(variance_function(logm, logprior, t) for t in pts)
        return step ** 2 * vmax <= np.log(2)

    lams, cur = [0.0], 0.0
    while cur < 1.0 - 1e-9 and len(lams) < 1000:
        hi_s = 1.0 - cur
        if admissible(cur, hi_s):
            step = hi_s
        else:
            lo_s = 0.0
            for _ in range(40):
                mid = 0.5 * (lo_s + hi_s)
                if admissible(cur, mid):
                    lo_s = mid
                else:
                    hi_s = mid
            step = max(lo_s, 1e-6)
        cur = min(1.0, cur + step)
        lams.append(cur)
    return np.array(lams), grid, V


# --------------------------------------------------------------------------
# 6. Truncation experiment (the slow-mixing example of YWJ, Appendix B)
# --------------------------------------------------------------------------

def neighbour_structure_gen(states, d, s0):
    """Same as neighbour_structure, with (d, s0) passed explicitly."""
    index = {s: i for i, s in enumerate(states)}
    rows, cols, prop = [], [], []
    for i, x in enumerate(states):
        xs = set(x)
        k = len(x)
        for j in range(d):
            if j in xs:
                y = tuple(sorted(xs - {j}))
            else:
                if k + 1 > s0:
                    continue
                y = tuple(sorted(xs | {j}))
            rows.append(i); cols.append(index[y]); prop.append(1.0 / (4 * d))
        if 0 < k < d:
            pr = 1.0 / (4.0 * k * (d - k))
            for a in x:
                for b in range(d):
                    if b in xs:
                        continue
                    y = tuple(sorted((xs - {a}) | {b}))
                    rows.append(i); cols.append(index[y]); prop.append(pr)
    return (np.asarray(rows), np.asarray(cols), np.asarray(prop))


def log_marginal_gen(states, X, y, n, g):
    yy = float(y @ y)
    out = np.empty(len(states))
    for i, x in enumerate(states):
        if len(x) == 0:
            r2 = 0.0
        else:
            Q, _ = np.linalg.qr(X[:, list(x)])
            r2 = min(float(np.sum((Q.T @ y) ** 2)) / yy, 1.0)
        out[i] = (-0.5 * len(x) * np.log1p(g)
                  - 0.5 * n * np.log1p(g * (1.0 - r2)))
    return out


def truncation_experiment(rng, d=12, alpha=1.5, kappa=2.0):
    """Section B of Yang, Wainwright and Jordan (2016): n = d, y = w pure
    noise, and no truncation.  The saturated model is then a trap, because
    R^2 = 1 there.  We sweep the truncation level s_0 and record Gap(K_1)."""
    n = d
    g = d ** (2 * alpha)
    X = rng.standard_normal((n, d))
    X /= np.linalg.norm(X, axis=0, keepdims=True) / np.sqrt(n)
    y = rng.standard_normal(n)          # y = w: the null model is true
    out = []
    for s0 in range(1, d + 1):
        states = enumerate_states(d, s0)
        nb = neighbour_structure_gen(states, d, s0)
        logm = log_marginal_gen(states, X, y, n, g)
        sizes = np.array([len(t) for t in states])
        logpi = logm - kappa * np.log(d) * sizes
        out.append((s0, len(states), spectral_gap(logpi, nb)))
    return out


# --------------------------------------------------------------------------
# 5. Experiment
# --------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(RNG_SEED)
    X = build_design(rng)
    w = SIGMA0 * rng.standard_normal(N_OBS)

    states = enumerate_states(D, S0)
    index = {s: i for i, s in enumerate(states)}
    nb = neighbour_structure(states)
    sizes = np.array([len(s) for s in states])
    logprior = -KAPPA * np.log(D) * sizes
    x_star_idx = index[S_SET]
    print(f"|X| = {len(states)}   (d={D}, s_0={S0})")

    # ---- signal-strength sweep -------------------------------------------
    # snr = n C_beta^2 / (sigma_0^2 log d), the scale of the beta-min condition
    snr_grid = np.concatenate([np.linspace(0.5, 12, 24), np.linspace(14, 40, 9)])
    lam_grid = np.linspace(0.0, 1.0, 41)

    gap1, gapmin, kmax, postmass, spread = [], [], [], [], []
    a0s, slacks, lamcs = [], [], []
    curves = {}
    highlight = {}
    for snr in snr_grid:
        cbeta = np.sqrt(snr * SIGMA0 ** 2 * np.log(D) / N_OBS)
        beta = np.zeros(D); beta[list(S_SET)] = cbeta
        y = X @ beta + w
        logm = log_marginal(states, X, y)
        lo, hi = psi_bounds(logm, index)
        spread.append(max(hi[j] - lo[j] for j in S_SET))

        gaps = []
        ks = []
        for lam in lam_grid:
            logpi = lam * logm + logprior
            gaps.append(spectral_gap(logpi, nb))
            ks.append(len(margin_failures(lam, lo, hi)))
        gaps = np.array(gaps)
        curves[round(float(snr), 3)] = (gaps, np.array(ks), lo, hi)

        a0_e, slack_e, lamc_e = null_exponent(logm, states, index)
        a0s.append(a0_e); slacks.append(slack_e); lamcs.append(lamc_e)
        p = np.exp(logm + logprior - (logm + logprior).max())
        p /= p.sum()
        postmass.append(float(p[x_star_idx]))
        gap1.append(gaps[-1]); gapmin.append(gaps.min()); kmax.append(max(ks))

    a0s = np.array(a0s); slacks = np.array(slacks); lamcs = np.array(lamcs)
    gap1 = np.array(gap1); gapmin = np.array(gapmin)
    kmax = np.array(kmax); postmass = np.array(postmass)
    spread = np.array(spread)

    # three representative signal strengths for panel (a): below, inside and
    # above the window in which the posterior mass moves from 0 to 1
    picks = [round(float(snr_grid[np.argmin(np.abs(snr_grid - t))]), 3)
             for t in (1.5, 12.0, 40.0)]

    # ---- schedule for the strongest signal --------------------------------
    cbeta = np.sqrt(picks[-1] * SIGMA0 ** 2 * np.log(D) / N_OBS)
    beta = np.zeros(D); beta[list(S_SET)] = cbeta
    y = X @ beta + w
    logm = log_marginal(states, X, y)
    fine = np.linspace(0.0, 2.0, 801)
    lams, vgrid, V = schedule(logm, logprior, fine)

    # ---- truncation sweep -------------------------------------------------
    trunc = truncation_experiment(np.random.default_rng(RNG_SEED + 1))

    # ---- dense cross-check of the sparse eigensolver -----------------------
    logpi_chk = 1.0 * logm + logprior
    rows, cols, prop = nb
    delta = logpi_chk[cols] - logpi_chk[rows]
    Kd = np.zeros((len(states), len(states)))
    Kd[rows, cols] = prop * np.exp(np.minimum(delta, 0.0))
    Kd[np.arange(len(states)), np.arange(len(states))] = 1.0 - Kd.sum(axis=1)
    piv = np.exp(logpi_chk - logpi_chk.max()); piv /= piv.sum()
    Sd = np.sqrt(piv)[:, None] * Kd / np.sqrt(piv)[None, :]
    ev = np.linalg.eigvalsh(0.5 * (Sd + Sd.T))
    gap_dense = 1.0 - ev[-2]
    gap_sparse = spectral_gap(logpi_chk, nb)

    # ---- figure -----------------------------------------------------------
    plt.rcParams.update({"font.size": 9.5, "font.family": "serif",
                         "axes.linewidth": 0.7, "lines.linewidth": 1.3})
    fig, ax2 = plt.subplots(2, 2, figsize=(6.7, 4.8))
    ax = [ax2[0, 0], ax2[0, 1], ax2[1, 0], ax2[1, 1]]

    # (a) gap along the path
    colours = ["#8a8a8a", "#c1121f", "#1b4965"]
    styles = ["-.", "--", "-"]
    for c, st, snr in zip(colours, styles, picks):
        gaps, ks, lo, hi = curves[snr]
        ax[0].semilogy(lam_grid, gaps, color=c, ls=st,
                       label=rf"$nC_\beta^2/(\sigma_0^2\log d)={snr:.0f}$")
    # shade where the margin condition (M_lambda) fails, for the strongest signal
    gaps, ks, lo, hi = curves[picks[-1]]
    thr, marg = KAPPA * np.log(D), np.log(2 * S_STAR)
    bad_iv = []
    for j in S_SET:
        if lo[j] <= 0:
            continue
        left, right = (thr - marg) / hi[j], (thr + marg) / lo[j]
        if right > 0 and left < 1:
            bad_iv.append((max(left, 0.0), min(right, 1.0)))
            ax[0].axvspan(bad_iv[-1][0], bad_iv[-1][1], color="#1b4965",
                          alpha=0.10, lw=0)
    ax[0].axhline(1 / (64 * D * S0), color="0.35", ls=":", lw=0.9)
    ax[0].text(0.02, 1 / (64 * D * S0) * 1.35, r"$1/(64ds_0)$",
               color="0.35", fontsize=8)
    ax[0].set_ylim(1e-4, 8e-2)
    ax[0].set_xlabel(r"temperature $\lambda$")
    ax[0].set_ylabel(r"$\mathrm{Gap}(K_\lambda)$")
    ax[0].set_title("(a) gap along the tempering path", fontsize=9.5)
    ax[0].legend(frameon=False, fontsize=7.5, loc="lower right")

    # (b) sweep through the signal-strength window
    ax[1].semilogy(snr_grid, gap1, color="#c1121f", marker="o", ms=3.0,
                   lw=1.0, label=r"$\mathrm{Gap}(K_1)$")
    ax[1].semilogy(snr_grid, gapmin, color="#1b4965", marker="s", ms=2.2,
                   lw=1.0, ls="--",
                   label=r"$\min_{\lambda}\mathrm{Gap}(K_\lambda)$")
    ax[1].axhline(1 / (64 * D * S0), color="0.35", ls=":", lw=0.9)
    ax[1].set_xlabel(r"$nC_\beta^2/(\sigma_0^2\log d)$")
    ax[1].set_ylabel("spectral gap")
    ax[1].set_title("(b) sweep of the signal strength", fontsize=9.5)
    ax[1].legend(frameon=False, fontsize=7.5, loc="lower right")
    ax[1].set_ylim(1e-4, 8e-2)
    axb = ax[1].twinx()
    axb.plot(snr_grid, postmass, color="0.55", lw=1.0, ls=":")
    axb.set_ylabel(r"$\pi(x^\star\mid y)$", color="0.45", fontsize=9,
                   labelpad=1)
    axb.tick_params(axis="y", colors="0.45", labelsize=8)
    axb.set_ylim(-0.03, 1.03)

    # (c) truncation sweep: the slow-mixing example of YWJ, Appendix B
    s0s = np.array([t[0] for t in trunc])
    gtr = np.array([t[2] for t in trunc])
    ax[2].semilogy(s0s, gtr, color="#1b4965", marker="o", ms=3.0)
    ax[2].set_xlabel(r"truncation level $s_0$  ($d=n=12$, $y=w$)")
    ax[2].set_ylabel(r"$\mathrm{Gap}(K_1)$")
    ax[2].set_title("(c) removing the truncation", fontsize=9.5)
    ax[2].set_xticks(s0s[::2])

    # (d) variance function and schedule
    sel = vgrid <= 1.0
    ax[3].plot(vgrid[sel], V[sel], color="#1b4965")
    ax[3].plot(lams, np.interp(lams, vgrid, V), "o", ms=3.0, color="#c1121f")
    ax[3].set_xlabel(r"temperature $\lambda$")
    ax[3].set_ylabel(r"$V(\lambda)=\mathrm{Var}_{\pi_\lambda}[\log m]$")
    ax[3].set_title(rf"(d) variance function and schedule, $T+1={len(lams)-1}$",
                    fontsize=9.5)

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    axb.spines[["top"]].set_visible(False)
    fig.tight_layout()
    fig.savefig("binary_numerics.pdf")

    # ---- numbers ----------------------------------------------------------
    with open("binary_numerics.txt", "w") as fh:
        def emit(line):
            print(line); fh.write(line + "\n")
        emit(f"seed={RNG_SEED}  n={N_OBS} d={D} s0={S0} s*={S_STAR} "
             f"rho={RHO} kappa={KAPPA} alpha={ALPHA} |X|={len(states)}")
        emit(f"1/(64 d s0) = {1/(64*D*S0):.3e}")
        emit("snr   Gap(K_1)     min_l Gap    k_max  Delta_int  pi(x*|y)"
             "    a_0    slack  lambda_c")
        for i, snr in enumerate(snr_grid):
            emit(f"{snr:5.1f} {gap1[i]:.4e}  {gapmin[i]:.4e}   "
                 f"{kmax[i]:d}     {spread[i]:7.3f}   {postmass[i]:.4f}"
                 f"  {a0s[i]:6.3f} {slacks[i]:7.3f}  {lamcs[i]:.3f}")
        emit(f"a_0 over the sweep: min {a0s.min():.3f} max {a0s.max():.3f}; "
             f"budget slack min {slacks.min():.3f}; "
             f"lambda_c min {lamcs.min():.3f}")
        emit(f"Gap(K_0) = {curves[picks[0]][0][0]:.6e},  "
             f"1/(4d) = {1/(4*D):.6e}")
        emit(f"gap range over all lambda and all snr: "
             f"[{gapmin.min():.4e}, {max(c[0].max() for c in curves.values()):.4e}]")
        emit(f"ratio min_l Gap / Gap(K_1), max over sweep: "
             f"{np.max(gapmin/gap1):.3f}")
        emit(f"min over sweep of min_l Gap: {gapmin.min():.4e}")
        emit(f"min over sweep of Gap(K_1):  {gap1.min():.4e}")
        for snr in picks:
            gaps, ks, lo, hi = curves[snr]
            emit(f"snr={snr}: psi- ={[round(lo[j],2) for j in S_SET]} "
                 f"psi+ ={[round(hi[j],2) for j in S_SET]} "
                 f"Delta_int={max(hi[j]-lo[j] for j in S_SET):.2f} "
                 f"k_1={ks[-1]} k_max={max(ks)}")
        emit(f"margin-failure intervals at snr={picks[-1]}: "
             f"{[(round(a,3), round(b,3)) for a, b in bad_iv]}")
        emit(f"schedule for snr={picks[-1]}: T+1={len(lams)-1}, "
             f"steps={np.round(np.diff(lams), 4).tolist()}")
        emit(f"sqrt(max V / log 2) = {np.sqrt(V.max()/np.log(2)):.2f}")
        emit(f"cross-check at lambda=1: sparse {gap_sparse:.6e} "
             f"dense {gap_dense:.6e}")
        emit("truncation sweep (d=n=12, y=w, kappa=2, alpha=1.5):")
        emit("  s_0   |X|      Gap(K_1)")
        for s0, sz, gp in trunc:
            emit(f"  {s0:3d} {sz:6d}   {gp:.4e}")
        emit(f"Gap(K_1) at s_0=4 / at s_0=d : "
             f"{trunc[3][2]/trunc[-1][2]:.3e}")


if __name__ == "__main__":
    main()
