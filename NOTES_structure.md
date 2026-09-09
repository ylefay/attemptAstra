# Where the polynomial loss comes from, and how to remove it

Working notes.  Numbers from `gapF_tmp.py`-style computations described in
section 3; nothing is asserted that has not been computed or proved here.
`binary_variable_selection3.tex` has not been modified.

---

## 0.  Answer, up front

**Every factor `s_0` in the current bound is the path method's *length*
charge, and no refinement of the path method removes it — I prove that below
(section 1.3).  Two things do remove it, and the second removes `s*` as
well:**

1. **A random-order null-stripping flow with Diaconis--Stroock weights**:
   `Gap >= c / (p log s_0)`.  Stays inside the path method; the only change
   is that the funnel strips a *uniformly chosen* null coordinate rather than
   the largest one, and edges at depth `k` get weight `k`.
2. **Block factorisation over {influential block} x {null singletons}**:
   `Gap >= c min{1/(4p), kappa_S(lambda)}`, where `kappa_S(lambda)` is the
   worst spectral gap of `P_lambda` restricted to a fibre
   `{gamma : gamma \ S = sigma}`.  **No `s_0`, no `s*`**, and `Theta(1/p)`
   exactly when `kappa_S = Omega(1)`.

Measured (four equicorrelated designs, 21 temperatures each, exact):
`Gap / min{1/(4p), kappa_S}` is **1.005** at its worst, the minimum being
attained at `1/(4p)`.  The certificate does not merely have the right order,
it is the gap.

`kappa_S(lambda)` is the interaction parameter the problem is asking for: it
is a gap on at most `2^{s*}` states, it equals `Omega(1)` on every benign
design tested, and it is what must collapse on a cancelling design.

---

## 1.  Audit: every `s_0` and `s*` in the current proof

### 1.1  The classical bound `Gap >= 1/(64 p s_0^2)` (`lemma:gap`)

| factor | source | verdict |
|---|---|---|
| `4p` | proposal: laziness `1/2` x flip-vs-exchange `1/2` x uniform coordinate `1/p` | intrinsic |
| first `s_0` | `eq:edge`: the exchange proposal is `1/(2|gamma|(p-|gamma|))`, bounded using `|gamma| <= s_0` | **crude edge bound** |
| second `s_0` | `ell(T) <= 4 s_0`: maximal path length | **length charge** |
| `2(2s*+1)` | precedent counting in `cor:precM` | **union bound** |

### 1.2  The current bound `eq:gapT1`

`Gap >= 1/{4 p e^{1/p} max(1, s* varrho) (2 s_0 + L)}`.

* the first `s_0` is **gone**: the funnel uses only flips, so
  `Q(e) = mu(gamma)/(4p)` exactly, with no `s_0`;
* the precedent count is **gone**: `2(2s*+1) = 14` became `e^{1/p} <= 1.07`;
* the `2^{s*}` is **gone**: replaced by the cut ratio, which is an identity;
* `s*` survives in `max(1, s* varrho)` and comes from the exchange proposal
  `1/(2|z|(p-|z|)) >= 1/(2 s* p)` at a state of size `|z| <= s*`.  Intrinsic
  to using exchanges;
* **`2 s_0` survives, and it is purely the maximal path length.**

### 1.3  The length charge cannot be removed by reweighting  [proved]

Diaconis--Stroock: for any positive edge weights `w`,

```
    Gap >= 1/A ,    A = ( max_e w(e) D(e)/Q(e) ) ( max_{x,y} sum_{e in T_xy} 1/w(e) ),
```

`D(e)` the demand.  For a funnel edge `e = (gamma, F(gamma))` at depth `k`
(that is, `|gamma \ gamma*| = k`),

```
    D(e) <= mu[Lambda^F(gamma)] <= e^{1/p} mu(gamma),      Q(e) = mu(gamma)/(4p),
```

so `D(e)/Q(e) <= 4 p e^{1/p}` **independently of `k`** — and this is attained,
since `Lambda^F(gamma) ∋ gamma` and the sinks sum to one.  Hence

```
    A >= 4p ( max_k w_k )( sum_{k=1}^{s_0} 1/w_k ) >= 4 p s_0 ,
```

the last step because `1/w_k >= 1/max_k w_k` for every `k`.  **No choice of
weights beats `4 p s_0`.**  The same computation applies to Sinclair's
length-weighted form (`rho-bar = max_e D(e) |T| / Q(e)`), which is the case
`w = 1`.

The obstruction is that the funnel concentrates *all* the traffic leaving a
state on one edge, so `D/Q` stays at `4p` all the way down while the paths
through the deep edges are long.

### 1.4  What the flat `D/Q` profile suggests

If instead the flow leaving a state with `k` null coordinates is *spread over
the `k` edges that remove one of them*, then `D(e)/Q(e) ≈ 4 p / k` at depth
`k`, and the weights `w_k = k` give

```
    max_k w_k D/Q = 4p ,     sum_k 1/w_k = H_{s_0} ≈ log s_0 ,
    A = O(p log s_0) .
```

**This is route 1: `Gap >= c/(p log s_0)`.**  It is a genuinely stronger
theorem, not better constants, and it stays inside the path method — the only
changes are (i) strip a uniformly chosen null coordinate instead of the
largest, (ii) weight depth-`k` edges by `k`.  The estimate `D(e) ≈ mu[Lambda]/k`
needs the precedent structure of the randomised funnel and is **not yet
proved**; it is the one gap in this route.

---

## 2.  Route 2: block factorisation, and why it is the right object

`ass:twosided` says that for every null covariate `j` and *every* context,
`mu[gamma_j = 1 | gamma_{-j}] <= p^{-kappa_lambda} <= p^{-2}`.  So the null
coordinates are uniformly weakly dependent on everything else: a Dobrushin
condition with interaction summed over `p` coordinates of order `p^{-1}`.
That is exactly the hypothesis under which **approximate block factorisation
of variance** holds with an absolute constant:

```
    Var_mu(f)  <=  C [ sum_{j in S^c} E_mu Var(f | gamma_{-j})
                       + E_mu Var(f | gamma_{S^c}) ] .                    (F)
```

The first sum is the Glauber form over the null coordinates, which the
chain's own Dirichlet form dominates at rate `1/(4p)` (a single flip is
proposed with probability `1/(2p)` and the chain is lazy: `eq:edge` is an
equality for flips).  The second term is the variance *inside a fibre*
`F_sigma = {gamma : gamma \ S = sigma}`, which the chain dominates at rate
`kappa_S(lambda)`, the worst fibre gap.  Hence

```
    Gap(P_lambda)  >=  1 / { C ( 4p + 1/kappa_S(lambda) ) }
                   >=  (1/2C) min{ 1/(4p), kappa_S(lambda) } .            (T)
```

**No `s_0` anywhere** — the null coordinates are handled one at a time with no
path and no length — **and no `s*`**: the whole influential block is one
factorisation block, so its cost enters only through `kappa_S`, a gap on at
most `2^{s*}` states, never through a route on them.

### Why this is the right parameter

* `kappa_S(lambda) = Omega(1)` on a benign design: the fibre is a
  `2^{s*}`-state chain with a well-separated mode.
* `kappa_S(lambda) -> 0` exactly on a cancelling design: the fibre measure has
  the genuine bottleneck between `emptyset` and `{0,1}` that
  `rem:necessary` describes, and no route or comparison can avoid it.
* It is the `K(X, beta*, lambda)` the problem statement asks for, and it is
  *necessary* by that same example.

---

## 3.  What was computed

Exact, `p = 14`, `n = 100`, `s_0 = 4`, 21 temperatures on `[0,1]`, the six
designs of `sec:numerics`.  For each `(design, lambda)`: the exact
`Gap(P_lambda)` by diagonalisation; every fibre gap
`Gap(P_lambda|_{F_sigma})` by diagonalisation of the rejection-restricted
kernel, over all `562` fibres; and the ratio

```
        Gap(P_lambda) / min{ 1/(4p), kappa_S(lambda) } .
```

| design | worst ratio `Gap / min{1/(4p), kappa_S}` | which term attains the min |
|---|---|---|
| `rho = 0` | 1.005 | `1/(4p)` at every `lambda` |
| `rho = 0.5` | 1.005 | `1/(4p)` at every `lambda` |
| `rho = 0.85` | 1.005 | `1/(4p)` at every `lambda` |
| `rho = 0.95` | 1.005 | `1/(4p)` at every `lambda` |
| cancelling `r = 0.90` | 2.05 | `kappa_S` for `lambda >= 0.2` |
| cancelling `r = 0.95` | 2.10 | `kappa_S` for `lambda >= 0.2` |

On every benign design the certificate is below the true gap by a factor
`1.005` — it is not an order-of-magnitude estimate, it is the gap.  On the
cancelling designs `kappa_S` collapses with the truth, over three orders of
magnitude:

| `lambda` | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|
| `Gap` (`r=0.95`) | 1.80e-2 | 1.32e-2 | 9.66e-5 | 2.20e-5 | 1.14e-5 | 6.01e-6 |
| `kappa_S` | 1.80e-2 | 6.56e-3 | 4.60e-5 | 1.39e-5 | 7.37e-6 | 3.95e-6 |
| ratio | 1.005 | 2.01 | 2.10 | 1.59 | 1.55 | 1.52 |

So over all six designs and all temperatures computed, the certificate
`min{1/(4p), kappa_S(lambda)}` is within a factor **2.1** of the exact gap.
For comparison the current `eq:gapT2` loses a factor `112` to `640`.

---

## 4.  Classification, in the terms asked for

| change | kind |
|---|---|
| removing the first `s_0` (exchange proposal at a saturated state) | **artifact of the proof technique**, already removed |
| removing the precedent count `2(2s*+1) -> e^{1/p}` | **better constants** |
| removing `2^{s*}` (cut ratio identity + crossover flow) | **artifact**, removed |
| `Gap >= c/(p log s_0)` by randomised stripping + weights | **genuinely stronger theorem** (one step unproved) |
| `Gap >= c min{1/(4p), kappa_S}` by block factorisation | **genuinely stronger theorem**; removes `s_0` and `s*` together |
| the surviving `1/p` | **intrinsic**: the chain proposes a given coordinate with probability `1/(2p)` |
| the surviving `kappa_S` | **intrinsic**: `rem:necessary` |

---

## 5.  What is not proved here

1. The constant `C` in (F).  It is what the Dobrushin/weak-dependence
   argument supplies from `ass:twosided`; an explicit `C = 4` under
   `ass:twosided` with `p^{1-eta} <= 1/4` is what a previous investigation
   obtained, and it should be re-derived rather than cited.
2. The demand estimate `D(e) ≈ mu[Lambda]/k` for the randomised funnel of
   route 1.
3. A lower bound on `kappa_S(lambda)` itself.  This is the honest residue of
   the whole problem: it is a gap on `2^{s*}` states, `rem:necessary` shows it
   must appear, and `C_lambda` and `M_lambda` of the current draft are two
   computable upper bounds on `1/kappa_S` of exactly the kind that a route on
   the influential block produces.
