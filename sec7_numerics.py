"""
Exact numerics for Section `sec:numerics`: the two-level ensemble of
`def:eroute`, its dip `Gamma_lambda`, its tie index `T_lambda`, and the bound
`eq:gapT2` against the exact spectral gap.

Everything is exact.  The state space X = {gamma in {0,1}^p : |gamma| <= s_0}
is enumerated, pi^(lambda) is evaluated at every state, the flip-and-exchange
kernel is assembled, its gap is obtained by diagonalisation, and the
congestion, the dip and the tie index of the ensemble are *computed*, not
bounded.  No Monte Carlo, no fitted constant.

The design and likelihood machinery is imported unchanged from
`window_numerics.py`; the cancelling design, the greedy chain, the E-route and
the three quantities of `def:dipA` are built here.

Outputs
-------
sec7_numerics.pdf   the four-panel figure of Section `sec:numerics`
sec7_numerics.txt   the numbers quoted in the text
sec7_sstar.pdf      the s*-sweep figure           (--sstar)
sec7_sstar.txt      the numbers of the s*-sweep   (--sstar)

Run:  python3 sec7_numerics.py           (the six-design experiment)
      python3 sec7_numerics.py --sstar   (the s*-sweep, at s_0 = 6)
"""

import itertools
import sys
from math import comb

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import window_numerics as W          # read-only reuse
import binary_numerics as B          # read-only reuse: the (d, s_0)-generic
                                     # kernel and marginal likelihood

N, P, S0 = W.N_OBS, W.D, W.S0
KAPPA, GP, SIG = W.KAPPA, W.G_PRIOR, W.SIGMA0
LOGP = np.log(P)
# alpha is *defined* by eq:alphadef, so that alpha log p = (1/2) log(1+g)
# exactly; W.ALPHA = 1.5 is the exponent used to build g = p^{2 alpha}, and
# the two differ by (1/2) log(1 + 1/g), which eq:wlam is sensitive to.
ALPHA = np.log1p(GP) / (2 * LOGP)
LAMS = np.linspace(0.0, 1.0, 41)


# --------------------------------------------------------------------------
# 1.  Designs
# --------------------------------------------------------------------------

def design_equicorrelated(rng, rho):
    """The four designs of `window_numerics.py`: s* = 3, separated
    coefficients, influential columns equicorrelated at rho."""
    X = W.build_design(rng, rho)
    beta = np.zeros(P)
    beta[:W.S_STAR] = W.BETA_S
    return X, beta, tuple(range(W.S_STAR))


def design_cancelling(rng, r, b):
    """The cancelling pair: X_0, X_1 correlated at r and beta* = (b, -b).
    Neither column explains anything alone, the pair explains everything, so
    the greedy gains increase and (B2) of `ass:BB` fails."""
    Z = rng.standard_normal((N, P))
    Sigma = np.array([[1.0, r], [r, 1.0]])
    L = np.linalg.cholesky(Sigma)
    X = np.empty((N, P))
    X[:, :2] = Z[:, :2] @ L.T
    Q, _ = np.linalg.qr(X[:, :2])
    X[:, 2:] = Z[:, 2:] - Q @ (Q.T @ Z[:, 2:])
    X /= np.linalg.norm(X, axis=0, keepdims=True) / np.sqrt(N)
    beta = np.zeros(P)
    beta[0], beta[1] = b, -b
    return X, beta, (0, 1)


# --------------------------------------------------------------------------
# 2.  The residual D_gamma of eq:Ddef, and the exponential family eq:wlam
# --------------------------------------------------------------------------

def log_D(X, y, subset):
    """log D_gamma = log{ g^{-1}||Y||^2 + ||(I - Phi_gamma) Y||^2 }."""
    yy = float(y @ y)
    if len(subset) == 0:
        rss = yy
    else:
        Q, _ = np.linalg.qr(X[:, list(subset)])
        rss = yy - float(np.sum((Q.T @ y) ** 2))
    return np.log(yy / GP + max(rss, 0.0))


def check_expfam(states, logm, X, y):
    """Verify eq:wlam: lambda logm + logprior and log w_lambda differ by a
    constant.  Returns the largest deviation over a lambda grid."""
    sizes = np.array([len(s) for s in states])
    logD = np.array([log_D(X, y, s) for s in states])
    worst = 0.0
    for lam in (0.0, 0.37, 1.0):
        a = lam * logm - KAPPA * LOGP * sizes
        bb = -0.5 * lam * N * logD - (KAPPA + lam * ALPHA) * sizes * LOGP
        d = a - bb
        worst = max(worst, float(d.max() - d.min()))
    return worst


# --------------------------------------------------------------------------
# 3.  The greedy chain of def:greedy and the conditions of ass:BB
# --------------------------------------------------------------------------

def greedy_chain(X, y, S_SET):
    """P_0 = () subset P_1 subset ... subset P_{s*} = gamma*, chosen by
    eq:greedychain, with the gains g_m of eq:gains."""
    chain, cur = [tuple()], []
    remaining = list(S_SET)
    while remaining:
        j = min(remaining, key=lambda k: log_D(X, y, tuple(sorted(cur + [k]))))
        cur = sorted(cur + [j])
        remaining.remove(j)
        chain.append(tuple(cur))
    gains = [0.5 * N * (log_D(X, y, chain[m]) - log_D(X, y, chain[m + 1]))
             - ALPHA * LOGP for m in range(len(chain) - 1)]
    return chain, np.array(gains)


def check_B1(X, y, S_SET, chain):
    """(B1) exchange descent: the matched exchanges of (E1) never increase D.
    Returns the largest increase of log D along any matched exchange."""
    worst = 0.0
    for k in range(len(S_SET) + 1):
        for A in itertools.combinations(S_SET, k):
            cur = set(A)
            a = sorted(set(A) - set(chain[k]))
            bb = sorted(set(chain[k]) - set(A))
            for ai, bi in zip(a, bb):
                before = log_D(X, y, tuple(sorted(cur)))
                cur = (cur - {ai}) | {bi}
                after = log_D(X, y, tuple(sorted(cur)))
                worst = max(worst, after - before)
    return worst


# --------------------------------------------------------------------------
# 4.  The E-route of def:eroute and the three quantities of def:dipA
# --------------------------------------------------------------------------

def neighbours(A, S_SET):
    """The neighbours of A inside 2^{gamma*}: flips and exchanges."""
    A = set(A)
    out = [tuple(sorted(A ^ {j})) for j in S_SET]
    out += [tuple(sorted((A - {a}) | {b}))
            for a in A for b in S_SET if b not in A]
    return out


def V_pot(A, chain, m_lam):
    """The potential V of eq:Vpot."""
    k = len(A)
    return len(set(A) - set(chain[k])) + abs(k - m_lam)


def ascent_map(subs, wl, S_SET, chain, m_lam):
    """F_lambda of eq:ascent: the w_lambda-largest admissible neighbour,
    admissible meaning of strictly smaller V.  wl is w_lambda on 2^{gamma*}
    (proportional to pi^(lambda) there, by eq:wlam)."""
    M = {}
    for A in subs:
        v = V_pot(A, chain, m_lam)
        cand = [z for z in neighbours(A, S_SET)
                if V_pot(z, chain, m_lam) < v]
        M[A] = max(cand, key=lambda z: wl[z]) if cand else A
    return M


def orbit(A, M):
    out = [A]
    while M[out[-1]] != out[-1]:
        out.append(M[out[-1]])
    return out


def route(A, Ap, orbs):
    """R_{A,A'}: the orbit of A, then the reverse of the orbit of A',
    spliced at their first common state."""
    oa, ob = orbs[A], orbs[Ap]
    pos = {z: j for j, z in enumerate(ob)}
    for i, z in enumerate(oa):
        if z in pos:
            return oa[:i + 1] + ob[:pos[z]][::-1]
    raise RuntimeError("orbits do not meet")


def route_quantities(subsets, pibar, wl, chain, m_lam, s_star, S_SET):
    """Gamma_lambda of eq:dip, the exact congestion varrho_lambda of
    eq:varrhoA, the longest route L_lambda, and the cut ratio C_lambda of
    eq:cut.  Also returns the crude precedent ratio max_z x(z)/pibar(z),
    which is the estimate rem:tie discards."""
    pb = pibar
    M = ascent_map(subsets, wl, S_SET, chain, m_lam)
    fixed = [A for A in subsets if M[A] == A]
    orbs = {A: orbit(A, M) for A in subsets}
    gamma_lam, demand, longest = 1.0, {}, 0
    for A in subsets:
        for Ap in subsets:
            if A == Ap:
                continue
            R = route(A, Ap, orbs)
            longest = max(longest, len(R) - 1)
            src = min(pb[A], pb[Ap])
            for z in R:
                gamma_lam = max(gamma_lam, src / pb[z])
            for z, zp in zip(R[:-1], R[1:]):
                demand[(z, zp)] = demand.get((z, zp), 0.0) + pb[A] * pb[Ap]
    varrho = (max(d / min(pb[z], pb[zp]) for (z, zp), d in demand.items())
              if demand else 1.0)
    prec = {A: 0.0 for A in subsets}
    for A in subsets:
        for z in orbs[A]:
            prec[z] += pb[A]
    moving = [z for z in subsets if M[z] != z]
    cut = max((prec[z] * (1 - prec[z]) / pb[z] for z in moving), default=0.0)
    crude = max(prec[z] / pb[z] for z in subsets)
    uphill = all(wl[M[A]] >= wl[A] for A in subsets)
    return (gamma_lam, varrho, longest, cut, crude, uphill, len(fixed))


def a0_value(states, X, y, S_SET):
    """The smallest a_0 for which ass:twosided holds:
    a_0 = alpha + max{ psi_j(gamma) : |gamma| < s_0, j in S^c \\ gamma } / log p.
    ass:twosided also asks kappa + alpha >= 2 + a_0."""
    best = -np.inf
    SC = [j for j in range(P) if j not in S_SET]
    for g in states:
        if len(g) >= S0:
            continue
        dg = log_D(X, y, g)
        for j in SC:
            if j in g:
                continue
            gj = tuple(sorted(set(g) | {j}))
            best = max(best, 0.5 * N * (dg - log_D(X, y, gj)) - ALPHA * LOGP)
    return ALPHA + best / LOGP, best


def sparse_eigenvalue(X, s):
    """nu of ass:B: the smallest eigenvalue of n^{-1} X_gamma^T X_gamma over
    all |gamma| <= s.  Exact, by enumeration."""
    G = X.T @ X / N
    best = np.inf
    for k in range(1, s + 1):
        for idx in itertools.combinations(range(P), k):
            ii = list(idx)
            best = min(best, float(np.linalg.eigvalsh(G[np.ix_(ii, ii)])[0]))
    return best


def cross_congestion(subsets, pibar, order):
    """Exact congestion of the crossover route of def:cross: set the
    coordinates to their A'-values in the greedy order."""
    dem = {}
    for A in subsets:
        for Ap in subsets:
            if A == Ap:
                continue
            cur = set(A)
            for j in order:
                if (j in Ap) != (j in cur):
                    nxt = set(cur)
                    nxt.discard(j) if j in cur else nxt.add(j)
                    z, zp = tuple(sorted(cur)), tuple(sorted(nxt))
                    dem[(z, zp)] = dem.get((z, zp), 0.0) + pibar[A]*pibar[Ap]
                    cur = nxt
    return (max(d / min(pibar[z], pibar[zp]) for (z, zp), d in dem.items())
            if dem else 1.0)


def cut_dependence(subsets, pibar, order):
    """M_lambda of eq:Mdef: the largest ratio pibar(a)pibar(b)/pibar(a u b)
    over the prefix/suffix cuts of the greedy order."""
    k = len(order); out = 1.0
    for t in range(k + 1):
        pre, suf = set(order[:t]), set(order[t:])
        ma, mb, mj = {}, {}, {}
        for A in subsets:
            a = tuple(sorted(set(A) & pre)); b = tuple(sorted(set(A) & suf))
            ma[a] = ma.get(a, 0.0) + pibar[A]
            mb[b] = mb.get(b, 0.0) + pibar[A]
            mj[(a, b)] = mj.get((a, b), 0.0) + pibar[A]
        for (a, b), v in mj.items():
            if v > 0:
                out = max(out, ma[a] * mb[b] / v)
    return out


def product_defect(subsets, pibar, S_SET):
    """r_lambda of eq:pdefect: the Chebyshev distance from log pibar to the
    additive functions, i.e. min over product nu of osc{log(pibar/nu)}.
    An exact linear programme: minimise t subject to
    |log pibar(A) - c - sum_{j in A} a_j| <= t/2."""
    from scipy.optimize import linprog
    k = len(S_SET)
    y = np.array([np.log(max(pibar[A], 1e-300)) for A in subsets])
    M = np.zeros((len(subsets), k + 1))
    M[:, 0] = 1.0
    for i, A in enumerate(subsets):
        for jj, j in enumerate(S_SET):
            M[i, jj + 1] = 1.0 if j in A else 0.0
    # variables (c, a_1..a_k, t);  M x - y <= t/2,  y - M x <= t/2
    nv = k + 2
    Aub = np.vstack([np.hstack([M, -0.5*np.ones((len(subsets), 1))]),
                     np.hstack([-M, -0.5*np.ones((len(subsets), 1))])])
    bub = np.concatenate([y, -y])
    c = np.zeros(nv); c[-1] = 1.0
    r = linprog(c, A_ub=Aub, b_ub=bub,
                bounds=[(None, None)]*(nv-1) + [(0, None)], method="highs")
    return float(r.x[-1]) if r.success else np.inf


# --------------------------------------------------------------------------
# 5.  One design, all temperatures
# --------------------------------------------------------------------------

def analyse(name, X, beta, S_SET, rng):
    s_star = len(S_SET)
    y = X @ beta + SIG * rng.standard_normal(N)
    states = W.enumerate_states(P, S0)
    index = {s: i for i, s in enumerate(states)}
    nb = W.neighbour_structure(states)
    logm = W.log_marginal(states, X, y)
    sizes = np.array([len(s) for s in states])
    logprior = -KAPPA * LOGP * sizes

    chain, gains = greedy_chain(X, y, S_SET)
    b1 = check_B1(X, y, S_SET, chain)
    # (B2) holds iff every increment of the gains is <= 0
    b2 = float(np.max(np.diff(gains))) if len(gains) > 1 else -np.inf
    subsets = [tuple(sorted(A)) for k in range(s_star + 1)
               for A in itertools.combinations(S_SET, k)]
    # block index: gamma -> gamma n gamma*
    blk = [tuple(sorted(set(s) & set(S_SET))) for s in states]

    nu = sparse_eigenvalue(X, S0 + s_star)
    a0, psimax_null = a0_value(states, X, y, S_SET)
    cbeta = float(np.min(np.abs(beta[list(S_SET)])))
    rec = {"name": name, "s_star": s_star, "chain": chain, "gains": gains,
           "nu": nu, "a0": a0, "psimax_null": psimax_null,
           "betamin_lhs": N * nu ** 2 * cbeta ** 2 / (SIG ** 2 * LOGP),
           "B1_slack": b1, "B2_slack": b2,
           "expfam": check_expfam(states, logm, X, y),
           "lam": [], "gap": [], "Gamma": [], "cut": [], "crude": [],
           "varrho": [], "m_lam": [], "bT1": [], "bT2": [], "funnel": [],
           "uphill": [], "nfix": [], "cross": [], "rdef": [],
           "Mdep": []}

    for lam in LAMS:
        logpi = lam * logm + logprior
        gap = W.spectral_gap(logpi, nb)
        w = np.exp(logpi - logpi.max())
        pi = w / w.sum()
        pibar = {A: 0.0 for A in subsets}
        for i, A in enumerate(blk):
            pibar[A] += pi[i]
        phi = np.array([lam * logm[index[c]] + logprior[index[c]]
                        for c in chain])
        m_lam = int(np.argmax(phi))
        wl = {A: pi[index[A]] for A in subsets}
        G, vr, L, cut, crude, up, nfix = route_quantities(
            subsets, pibar, wl, chain, m_lam, s_star, S_SET)
        funnel = max(pibar[A] / pi[index[A]] for A in subsets)
        rec["lam"].append(lam);      rec["gap"].append(gap)
        rec["Gamma"].append(G);      rec["cut"].append(cut)
        rec["crude"].append(crude);  rec["varrho"].append(vr)
        rec["m_lam"].append(m_lam);  rec["funnel"].append(funnel)
        rec["uphill"].append(up);    rec["nfix"].append(nfix)
        order = [c for m in range(s_star) for c in chain[m+1]
                 if c not in chain[m]]
        rec["cross"].append(cross_congestion(subsets, pibar, order))
        rec["Mdep"].append(cut_dependence(subsets, pibar, order))
        rec["rdef"].append(product_defect(subsets, pibar, S_SET))
        rec["bT1"].append(1.0 / (4 * P * np.exp(1 / P)
                                 * max(1.0, s_star * vr) * (2 * S0 + L)))
        rec["bT2"].append(1.0 / (8 * P * (2 * S0 + 4 * s_star)
                                 * max(1.0, 2 * s_star * cut)))
        # the theorem pays the smaller of the two flows
    for k in ("lam", "gap", "Gamma", "cut", "crude", "varrho", "m_lam",
              "bT1", "bT2", "funnel", "uphill", "nfix", "cross", "rdef",
              "Mdep"):
        rec[k] = np.asarray(rec[k])
    return rec


# --------------------------------------------------------------------------
# 6.  Figure and tables
# --------------------------------------------------------------------------

def figure(recs, fname):
    plt.rcParams.update({"font.size": 9.5, "font.family": "serif",
                         "axes.linewidth": 0.7, "lines.linewidth": 1.3})
    fig, ax2 = plt.subplots(2, 2, figsize=(6.7, 4.9))
    ax = [ax2[0, 0], ax2[0, 1], ax2[1, 0], ax2[1, 1]]
    eq = [r for r in recs if r["name"].startswith("rho")]
    ca = [r for r in recs if r["name"].startswith("cancel")]
    cols = ["#8a8a8a", "#4a7ba7", "#c1121f", "#1b4965"]
    dark = ["#7a1f1f", "#000000"]

    # (a) the centre moves with the temperature
    for c, st, r in zip(cols, ["-", "--", "-.", ":"], eq):
        ax[0].step(r["lam"], r["m_lam"], where="post", color=c, ls=st,
                   label=r["name"])
    ax[0].set_xlabel(r"temperature $\lambda$")
    ax[0].set_ylabel(r"$m_\lambda$")
    ax[0].set_yticks(range(0, 4))
    ax[0].set_title(r"(a) the centre $P_{m_\lambda}$ moves", fontsize=9.5)
    ax[0].set_ylim(-0.25, 3.6)
    ax[0].legend(frameon=False, fontsize=7, loc="upper left")

    # (b) the cut ratio, against the true deficit
    for c, r in zip(cols, eq):
        ax[1].semilogy(r["lam"], r["cut"], color=c)
    for c, r in zip(dark, ca):
        ax[1].semilogy(r["lam"], r["cut"], color=c, ls="-",
                       label=r["name"])
        ax[1].semilogy(r["lam"], (1 / (4 * P)) / r["gap"], color=c, ls=":",
                       lw=1.0)
    ax[1].set_xlabel(r"temperature $\lambda$")
    ax[1].set_ylabel(r"$C_\lambda$")
    ax[1].set_title(r"(b) the cut ratio (dotted: "
                    r"$(4p)^{-1}/\mathrm{Gap}$)", fontsize=9.5)
    ax[1].legend(frameon=False, fontsize=7, loc="upper left")

    # (c) the bound against the truth
    show = [eq[2], ca[-1]]
    for c, r in zip(["#1b4965", "#c1121f"], show):
        ax[2].semilogy(r["lam"], r["gap"], color=c, label=r["name"])
        ax[2].semilogy(r["lam"], r["bT2"], color=c, ls="--")
        ax[2].semilogy(r["lam"], r["bT1"], color=c, ls=":")
    ax[2].set_xlabel(r"temperature $\lambda$")
    ax[2].set_ylabel("spectral gap")
    ax[2].set_title(r"(c) truth (---), (T2) (- -), (T1) ($\cdots$)",
                    fontsize=9.5)
    ax[2].legend(frameon=False, fontsize=7, loc="lower left")

    # (d) the tie index does not degrade
    for c, st, r in zip(cols, ["-", "--", "-.", ":"], eq):
        ax[3].semilogy(r["lam"], r["cut"], color=c, ls=st)
        ax[3].semilogy(r["lam"], r["crude"], color=c, ls=":", lw=0.8)
    for c, r in zip(dark, ca):
        ax[3].semilogy(r["lam"], r["cut"], color=c)
    ax[3].axhline(2 ** 3, color="0.35", ls="--", lw=0.9)
    ax[3].text(0.02, 2 ** 3 * 1.15, r"$2^{s^\ast}$", color="0.35",
               fontsize=8)
    ax[3].set_xlabel(r"temperature $\lambda$")
    ax[3].set_ylabel(r"$C_\lambda$")
    ax[3].set_title(r"(d) the cut ratio (dotted: the count it replaces)",
                    fontsize=9.5)

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fname)


def main():
    recs = []
    for rho in (0.0, 0.5, 0.85, 0.95):
        rng = np.random.default_rng(W.RNG_SEED)
        X, beta, S_SET = design_equicorrelated(rng, rho)
        recs.append(analyse(f"rho={rho}", X, beta, S_SET,
                            np.random.default_rng(W.RNG_SEED + 1)))
    for r, b in ((0.90, 3.0), (0.95, 3.0)):
        rng = np.random.default_rng(W.RNG_SEED + 7)
        X, beta, S_SET = design_cancelling(rng, r, b)
        recs.append(analyse(f"cancel r={r}", X, beta, S_SET,
                            np.random.default_rng(W.RNG_SEED + 1)))

    figure(recs, "sec7_numerics.pdf")

    with open("sec7_numerics.txt", "w") as fh:
        def emit(s):
            print(s); fh.write(s + "\n")
        emit(f"n={N} p={P} s_0={S0} kappa={KAPPA} alpha={ALPHA} "
             f"g=p^(2alpha)={GP:.0f}  |X|={len(W.enumerate_states(P,S0))}")
        emit(f"lambda grid: {len(LAMS)} points on [0,1];  1/(4p) = "
             f"{1/(4*P):.4e}")
        emit("")
        emit("design        s*  eq:wlam dev   (B1) slack   (B2) slack   "
             "greedy gains g_m")
        for r in recs:
            emit(f"{r['name']:<13s} {r['s_star']:2d}  {r['expfam']:.2e}     "
                 f"{r['B1_slack']:+.3e}  {r['B2_slack']:+.3e}   "
                 f"{np.round(r['gains'],2).tolist()}")
        emit("")
        emit("design         max C_lam  max varrho  max funnel  "
             "min Gap    min bT2   min bT1   worst bT2/Gap  worst bT1/Gap")
        for r in recs:
            emit(f"{r['name']:<13s} {r['cut'].max():9.3f} "
                 f"{r['varrho'].max():10.3e} "
                 f"{r['funnel'].max():10.4f} {r['gap'].min():.3e} "
                 f"{r['bT2'].min():.3e} {r['bT1'].min():.3e} "
                 f"{(r['bT2']/r['gap']).max():13.3e} "
                 f"{(r['bT1']/r['gap']).max():13.3e}")
        emit("")
        emit("")
        emit("the two flows of th:gapA: the ascent (cut ratio C_lambda) and "
             "the crossover (congestion, and product defect r_lambda):")
        emit("design         max C_lam   max crossover     max M_lam  "
             "M/crossover  max r_lam   max min{C,M}")
        for r in recs:
            mn = np.minimum(r["cut"], r["Mdep"]).max()
            emit(f"  {r['name']:<13s} {r['cut'].max():9.3f} "
                 f"{r['cross'].max():15.4e} {r['Mdep'].max():13.4e} "
                 f"{(r['Mdep']/np.maximum(r['cross'],1e-300)).max():11.4f} "
                 f"{r['rdef'].max():10.3f} {mn:14.3f}")
        emit("")
        emit("the ascent F_lambda: is every step uphill, and is the fixed "
             "point unique?  (prop:dip predicts yes iff (B1) and (B2) hold)")
        for r in recs:
            emit(f"  {r['name']:<13s} uphill at "
                 f"{int(r['uphill'].sum())}/{len(r['lam'])} temperatures;  "
                 f"unique fixed point at "
                 f"{int((r['nfix']==1).sum())}/{len(r['lam'])};  "
                 f"max crude/cut ratio {(r['crude']/np.maximum(r['cut'],1e-300)).max():.1f}")
        emit("")
        emit("ass:twosided (Assumption 6.1): a_0 = alpha + max_{j in S^c} "
             "psi_j / log p, and it asks kappa + alpha >= 2 + a_0, "
             f"here {KAPPA + ALPHA:.4f} >= 2 + a_0:")
        for r in recs:
            ok = (KAPPA + ALPHA) >= 2 + r["a0"]
            emit(f"  {r['name']:<13s} max_{{S^c}} psi_j = "
                 f"{r['psimax_null']:+8.4f}   a_0 = {r['a0']:7.4f}   "
                 f"kappa_lambda >= {KAPPA - max(0.0, r['a0'] - ALPHA):.4f}   "
                 f"{'HOLDS' if ok else 'FAILS'}")
        emit("")
        emit("the regime: ass:betamin asks "
             "n nu^2 C_beta^2/(sigma_0^2 log p) >= 128(L+Ltilde+alpha+kappa),"
             " which is at least 832 under ass:B and ass:C.")
        for r in recs:
            emit(f"  {r['name']:<13s} nu = {r['nu']:.4f}   "
                 f"n nu^2 C_beta^2/(sigma_0^2 log p) = "
                 f"{r['betamin_lhs']:.3e}   "
                 f"shortfall factor {832/r['betamin_lhs']:.3e}")
        emit("")
        emit(f"e^(1/p) = {np.exp(1/P):.4f}   "
             f"largest funnel constant over all designs and temperatures: "
             f"{max(r['funnel'].max() for r in recs):.4f}")
        emit(f"validity: min over all designs and temperatures of "
             f"Gap - bT2 = "
             f"{min((r['gap']-r['bT2']).min() for r in recs):.3e}  "
             f"(must be >= 0)")
        emit(f"validity: min over all designs and temperatures of "
             f"Gap - bT1 = "
             f"{min((r['gap']-r['bT1']).min() for r in recs):.3e}  "
             f"(must be >= 0)")
        emit("")
        emit("does the cut ratio track the true deficit?  "
             "C_lambda against (4p)^{-1}/Gap(lambda), at lambda = 1:")
        for r in recs:
            defic = (1 / (4 * P)) / r["gap"][-1]
            emit(f"  {r['name']:<13s} C_1 = {r['cut'][-1]:9.3e}   "
                 f"deficit = {defic:9.3e}   ratio = "
                 f"{r['cut'][-1]/defic:6.3f}")
        emit("")
        emit("per-temperature detail, cancelling r=0.95 "
             "(lambda, m_lam, C_lam, varrho, Gap, bT1, bT2):")
        r = recs[-1]
        for i in range(0, len(LAMS), 4):
            emit(f"  {r['lam'][i]:.3f}  {r['m_lam'][i]:d}  "
                 f"{r['cut'][i]:.3e}  "
                 f"{r['varrho'][i]:.3e}  {r['gap'][i]:.3e}  "
                 f"{r['bT1'][i]:.3e}  {r['bT2'][i]:.3e}")
    return recs



# --------------------------------------------------------------------------
# 7.  The s*-sweep: how the gap, and the bound, depend on the number of
#     influential covariates.  Here s_0 = 6, so that s* can run to 5.
# --------------------------------------------------------------------------

S0_BIG = 6
SSTARS = (1, 2, 3, 4, 5)
PROFILE = np.array([1.20, 1.00, 0.85, 0.70, 0.45])
LAMS_BIG = np.linspace(0.0, 1.0, 21)


def design_equicorrelated_gen(rng, rho, s_star):
    """As design_equicorrelated, for an arbitrary number s* of influential
    columns: they are equicorrelated at rho, the rest are residualised
    against their span, and every column is normalised to ||X_j||^2 = n."""
    Z = rng.standard_normal((N, P))
    Sigma = (1 - rho) * np.eye(s_star) + rho * np.ones((s_star, s_star))
    L = np.linalg.cholesky(Sigma)
    X = np.empty((N, P))
    X[:, :s_star] = Z[:, :s_star] @ L.T
    Q, _ = np.linalg.qr(X[:, :s_star])
    X[:, s_star:] = Z[:, s_star:] - Q @ (Q.T @ Z[:, s_star:])
    X /= np.linalg.norm(X, axis=0, keepdims=True) / np.sqrt(N)
    beta = np.zeros(P)
    beta[:s_star] = PROFILE[:s_star]
    return X, beta, tuple(range(s_star))


def log_D_big(X, y, subset):
    """log D_gamma of eq:Ddef; identical to log_D, kept separate only because
    the sweep uses its own truncation level."""
    return log_D(X, y, subset)


def analyse_sstar(rho, s_star, states, nb, rng):
    """One (design, s*) pair: the exact gap, the dip, the tie index and the
    two bounds of th:gapA, at every temperature of LAMS_BIG."""
    X, beta, S_SET = design_equicorrelated_gen(
        np.random.default_rng(W.RNG_SEED + 100 * s_star), rho, s_star)
    y = X @ beta + SIG * rng.standard_normal(N)
    index = {st: i for i, st in enumerate(states)}
    logm = B.log_marginal_gen(states, X, y, N, GP)
    sizes = np.array([len(st) for st in states])
    logprior = -KAPPA * LOGP * sizes

    chain, gains = greedy_chain(X, y, S_SET)
    b1 = check_B1(X, y, S_SET, chain)
    b2 = float(np.max(np.diff(gains))) if len(gains) > 1 else -np.inf
    subsets = [tuple(sorted(A)) for k in range(s_star + 1)
               for A in itertools.combinations(S_SET, k)]
    blk = [tuple(sorted(set(st) & set(S_SET))) for st in states]

    out = {"rho": rho, "s_star": s_star, "gains": gains,
           "B1_slack": b1, "B2_slack": b2,
           "gap": [], "Gamma": [], "cut": [], "crude": [], "varrho": [],
           "bT1": [], "bT2": [], "uphill": []}
    for lam in LAMS_BIG:
        logpi = lam * logm + logprior
        out["gap"].append(W.spectral_gap(logpi, nb))
        w = np.exp(logpi - logpi.max())
        pi = w / w.sum()
        pibar = {A: 0.0 for A in subsets}
        for i, A in enumerate(blk):
            pibar[A] += pi[i]
        phi = np.array([lam * logm[index[c]] + logprior[index[c]]
                        for c in chain])
        m_lam = int(np.argmax(phi))
        wl = {A: pi[index[A]] for A in subsets}
        G, vr, L, cut, crude, up, nfix = route_quantities(
            subsets, pibar, wl, chain, m_lam, s_star, S_SET)
        out["Gamma"].append(G); out["cut"].append(cut)
        out["crude"].append(crude); out["varrho"].append(vr)
        out["uphill"].append(up)
        out["bT1"].append(1.0 / (4 * P * np.exp(1 / P)
                                 * max(1.0, s_star * vr) * (2 * S0_BIG + L)))
        out["bT2"].append(1.0 / (8 * P * (2 * S0_BIG + 4 * s_star)
                                 * max(1.0, 2 * s_star * cut)))
    for k in ("gap", "Gamma", "cut", "crude", "varrho", "bT1", "bT2",
              "uphill"):
        out[k] = np.asarray(out[k])
    return out


def figure_sstar(rows, fname):
    plt.rcParams.update({"font.size": 9.5, "font.family": "serif",
                         "axes.linewidth": 0.7, "lines.linewidth": 1.3})
    fig, ax2 = plt.subplots(2, 2, figsize=(6.7, 4.9))
    ax = [ax2[0, 0], ax2[0, 1], ax2[1, 0], ax2[1, 1]]
    cols = ["#8a8a8a", "#4a7ba7", "#c1121f", "#1b4965"]
    styles = ["-", "--", "-.", ":"]
    rhos = sorted({r["rho"] for r in rows})
    ss = np.array(SSTARS)

    def get(rho, k):
        return next(r for r in rows if r["rho"] == rho and r["s_star"] == k)

    def series(rho, key, how):
        return np.array([how(get(rho, k)[key]) for k in SSTARS])

    # (a) the exact gap, in units of 1/(4p)
    for c, st, rho in zip(cols, styles, rhos):
        v = 4 * P * series(rho, "gap", np.min)
        ax[0].semilogy(ss, v, color=c, ls=st, lw=1.1, label=rf"$\rho={rho}$")
        b2 = np.array([get(rho, k)["B2_slack"] for k in SSTARS])
        ok = b2 <= 0
        ax[0].plot(ss[ok], v[ok], "o", color=c, ms=4.0)
        ax[0].plot(ss[~ok], v[~ok], "o", color=c, ms=4.5, mfc="none")
    ax[0].axhline(1.0, color="0.35", ls=":", lw=0.9)
    ax[0].set_xlabel(r"$s^\ast$")
    ax[0].set_ylabel(r"$4p\,\min_\lambda\,\mathrm{Gap}(\mathbf{P}_\lambda)$")
    ax[0].set_title(r"(a) the exact gap; open marker: (B2) fails",
                    fontsize=9.5)
    ax[0].set_xticks(ss)
    ax[0].legend(frameon=False, fontsize=7, loc="lower left", ncol=2)

    # (b) the dip against the true deficit
    for c, st, rho in zip(cols, styles, rhos):
        ax[1].semilogy(ss, series(rho, "cut", np.max), color=c, ls=st,
                       marker="o", ms=3.0, label=rf"$\rho={rho}$")
        ax[1].semilogy(ss, 1.0 / (4 * P * series(rho, "gap", np.min)),
                       color=c, ls=":", lw=0.9, marker="s", ms=2.2)
    ax[1].set_xlabel(r"$s^\ast$")
    ax[1].set_ylabel(r"$\max_\lambda\,C_\lambda$")
    ax[1].set_title(r"(b) the cut ratio (dotted: "
                    r"$(4p\,\mathrm{Gap})^{-1}$)", fontsize=9.5)
    ax[1].set_xticks(ss)

    # (c) the certificate
    for c, st, rho in zip(cols, styles, rhos):
        ax[2].semilogy(ss, series(rho, "bT2", np.min), color=c, ls=st,
                       marker="o", ms=3.0)
    ref = (series(rhos[-1], "bT2", np.min)[0]
           * (1 * (2 * S0_BIG + 4)) / (ss * (2 * S0_BIG + 4 * ss)))
    ax[2].semilogy(ss, ref, color="0.35", ls="--", lw=0.9)
    ax[2].text(2.6, ref[1] * 1.5,
               r"$\propto\{s^\ast(2s_0+4s^\ast)\}^{-1}$", color="0.35",
               fontsize=8)
    ax[2].set_xlabel(r"$s^\ast$")
    ax[2].set_ylabel(r"$\min_\lambda$ bound")
    ax[2].set_title(r"(c) the certificate decays polynomially",
                    fontsize=9.5)
    ax[2].set_xticks(ss)

    # (d) the tie index against the binomial count it is
    for c, st, rho in zip(cols, styles, rhos):
        ax[3].semilogy(ss, series(rho, "cut", np.max), color=c, ls=st,
                       marker="o", ms=3.0)
    ax[3].semilogy(ss, [2.0 ** k for k in SSTARS], color="0.35", ls="--",
                   lw=0.9)
    ax[3].text(3.3, 2.0 ** 3.6, r"$2^{s^\ast}$", color="0.35", fontsize=8)
    ax[3].set_xlabel(r"$s^\ast$")
    ax[3].set_ylabel(r"$\max_\lambda\,C_\lambda$")
    ax[3].set_title(r"(d) the cut ratio does not grow with $s^\ast$",
                    fontsize=9.5)
    ax[3].set_xticks(ss)

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fname)


def main_sstar():
    states = W.enumerate_states(P, S0_BIG)
    nb = B.neighbour_structure_gen(states, P, S0_BIG)
    rows = []
    for rho in (0.0, 0.5, 0.85, 0.95):
        for s_star in SSTARS:
            rows.append(analyse_sstar(
                rho, s_star, states, nb,
                np.random.default_rng(W.RNG_SEED + 1)))
    figure_sstar(rows, "sec7_sstar.pdf")

    with open("sec7_sstar.txt", "w") as fh:
        def emit(t):
            print(t); fh.write(t + "\n")
        emit(f"s*-sweep: n={N} p={P} s_0={S0_BIG} kappa={KAPPA} "
             f"alpha={ALPHA:.5f}  |X|={len(states)}  "
             f"{len(LAMS_BIG)} temperatures on [0,1]")
        emit(f"beta*_S = the first s* entries of "
             f"{np.round(PROFILE,2).tolist()};  1/(4p) = {1/(4*P):.4e}")
        emit("")
        emit("rho    s*  (B1)     (B2)      min Gap    Gap(1)   "
             "4p*minGap  max C_lam  max varrho  min bound  "
             "bound/Gap  uphill")
        for r in rows:
            emit(f"{r['rho']:<5.2f} {r['s_star']:2d}  "
                 f"{r['B1_slack']:+.1e} {r['B2_slack']:+.1e}  "
                 f"{r['gap'].min():.3e}  {r['gap'][-1]:.3e}  "
                 f"{4*P*r['gap'].min():8.4f}  {r['cut'].max():9.3f}  "
                 f"{r['varrho'].max():10.3f}  "
                 f"{r['bT2'].min():.3e}  "
                 f"{(r['bT2']/r['gap']).max():.3e}  "
                 f"{int(r['uphill'].sum())}/{len(S.LAMS_BIG) if False else len(r['gap'])}")
        emit("")
        emit("greedy gains g_m, by design:")
        for r in rows:
            emit(f"  rho={r['rho']:.2f} s*={r['s_star']}: "
                 f"{np.round(r['gains'],2).tolist()}")
        emit("")
        g0 = min(4 * P * r["gap"].min() for r in rows)
        g1 = max(4 * P * r["gap"].min() for r in rows)
        emit(f"4p * min_lambda Gap over the whole sweep: [{g0:.4f}, "
             f"{g1:.4f}] -- the truth stays within a factor "
             f"{g1/g0:.2f} of 1/(4p) as s* runs from 1 to 5")
        r1 = [r for r in rows if r["s_star"] == 1]
        r5 = [r for r in rows if r["s_star"] == 5]
        emit(f"min bound at s*=1: {min(r['bT2'].min() for r in r1):.3e};  "
             f"at s*=5: {min(r['bT2'].min() for r in r5):.3e};  "
             f"ratio {min(r['bT2'].min() for r in r1)/min(r['bT2'].min() for r in r5):.1f}")
        emit(f"validity over the sweep: min(Gap - bound) = "
             f"{min((r['gap']-r['bT2']).min() for r in rows):.3e} "
             f"(must be >= 0)")
    return rows


if __name__ == "__main__":
    if "--sstar" in sys.argv:
        main_sstar()
    else:
        main()
