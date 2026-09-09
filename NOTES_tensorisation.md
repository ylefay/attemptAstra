# Does the variance of `pi^(lambda)` factorise approximately?

Working notes.  Every number below is produced by `tensor_numerics.py`:
`python3 tensor_numerics.py` writes `tensor_numerics.txt`, `--quick` uses fewer
temperatures, and `--screen` runs the random-design study of section 6.  All of
it is exact: `M(s_0)` is enumerated, `pi^(lambda)` is evaluated at every state,
and every constant is the exact extremum of a Rayleigh quotient obtained by
dense diagonalisation.  Nothing in
`binary_variable_selection3.tex` or `bvs3_app.tex` has been touched.

---

## 0. Answer, up front

**`Theta(1/p)` is provable by approximate tensorisation, but not over single
coordinates: single-site tensorisation is false here, and what is true is a
*block* factorisation in which the `s*` influential coordinates are lumped into
one block.  That factorisation holds with constant `4` under Assumption 6.1
alone (with `p^{1-eta} <= 1/4` in place of `<= 1/2`), and it gives**

```
                                       1                    theta_S
    Gap(P_lambda)   >=   ---------------------------   =   -------------- ,
                          4 { 4p + 1/kappa_S(lambda) }      16 p (1+theta_S)

    kappa_S(lambda) = min over sigma subset of S^c of the spectral gap of
                      P_lambda restricted to the fibre {gamma : gamma n S^c = sigma},
    theta_S         = 4 p kappa_S .
```

**This is `Theta(1/p)` if and only if `theta_S = Omega(1)`, which is a
hypothesis about a chain on at most `2^{s*}` states.  Assumptions 3.1--3.5 do
not deliver it and it is false for the cancelling pair.  The missing hypothesis
is therefore exactly a conductance bound for the conditional law of `gamma_S`
given `gamma_{S^c}` — a hypothesis of the same kind as Assumption 6.3, but
quantitative and about the fibre measure on `2^S` rather than the block measure
on `2^{W_lambda}`.  It is not a new kind of hypothesis, and it is not one that
Assumptions 3.1--3.5 can be sharpened into: section 6 exhibits a collapsed
design with a *larger* `nu` and a *smaller* `omega(X)` than a benign one.**

The value of the route is that the certificate *tracks* the truth.  Over the
126 (design, temperature) pairs computed, the proved bound is between `1/23`
and `1/4.6` of the exact gap — **including** the cancelling pair, where both
the truth and the certificate fall by three orders of magnitude together.  For
comparison, `eq:gapM0` is a factor `2.2e2` below the truth at `s* = 3`, the
two-level route of `NOTES_flow.md` a factor `52` at `s* = 5` where its
routing hypothesis holds, and `eq:gapM` of the previous draft a factor `10^{30}`
and more where it does not.

Acceptance test, as required: at `p=14, s_0=4, n=100`, `beta* = (3,-3)`,
`r = 0.95`, the exact single-site tensorisation constant is `C_atv = 2238`
(and `702` at `r = 0.9`), against `1.0` to `2.0` for the benign designs.  The
constant does degrade, by three orders of magnitude, exactly where it must.

---

## 1. Setting

Throughout, `mu = pi_n^{(lambda)}(. | Y)` on `M(s_0)`, and

```
    D_gamma  = g^{-1}||Y||^2 + ||(I - Phi_gamma)Y||^2  =  ||Y||^2 (g^{-1}+1-R^2_gamma),
    g_j(gamma) = ||(Phi_{gamma u j} - Phi_gamma)Y||^2 ,   D_{gamma u j} = D_gamma - g_j(gamma),
    psi_j(gamma) = -(1/2)log(1+g) - (n/2) log(1 - g_j(gamma)/D_gamma) .
```

The single-site conditional of the target is a logistic function of the
one-step log-ratio of `eq:onestepG`:

```
                                { sigmoid( lambda psi_j(gamma_{-j}) - kappa log p ),  |gamma_{-j}| < s_0,
    mu[ gamma_j = 1 | gamma_{-j} ] = 
                                { 0,                                                  |gamma_{-j}| = s_0.       (1)
```

Two elementary facts are used repeatedly.

**(F1)  The chain's Dirichlet form dominates the Glauber form.**  For a single
flip, `eq:edge` is an equality: `Q_lambda(gamma, gamma^j) = (4p)^{-1}
min{mu(gamma), mu(gamma^j)}`.  Since `min(a,b) >= ab/(a+b)`,

```
    sum_{j in A} E[ Var(f | gamma_{-j}) ]  <=  4 p  E_{P_lambda}(f,f)          (2)
```

for every set `A` of coordinates, `E_{P_lambda}` being the Dirichlet form of
the lazy chain of `eq:MH`.  In particular a single-site tensorisation with
constant `C` would give `Gap(P_lambda) >= 1/(4 p C)`.

**(F2)  The influence of `i` on `j` is the mixed second difference of
`log L_n`.**  For `gamma` with `|gamma| <= s_0 - 2` and `i != j` outside
`gamma`,

```
    delta_ij(gamma) := psi_j(gamma u {i}) - psi_j(gamma)
                     = (n/2) log{ D_{gamma u i} D_{gamma u j} / ( D_gamma D_{gamma u i u j} ) } .  (3)
```

The prior term `kappa log p` and the Occam factor `(1/2)log(1+g)` cancel
identically, so `delta_ij` is a functional of the design and the data alone,
free of `kappa`, `g` and `lambda`.  It is symmetric in `i, j`.
(`tensor_numerics.py` checks (3) against a direct evaluation of
`psi_j(gamma u i) - psi_j(gamma)`: agreement to `2.1e-13` on every design.)

By (1), the Dobrushin influence of `j` on `i`,
`R_ij = max { |mu(gamma_i=1|.) - mu(gamma_i=1|.')| : the two contexts differ
only at j }`, obeys the two bounds

```
    R_ij  <=  (lambda/4) max_gamma |delta_ij(gamma)|          (sigmoid' <= 1/4),   (4a)
    R_ij  <=  max_gamma  mu[ gamma_i = 1 | gamma_{-i} ] .                          (4b)
```

`(4b)` needs no smoothness and is unaffected by the truncation, which can only
*lower* a conditional probability to `0`.  `(4a)` is stated for the pairs of
contexts at which both conditionals are non-degenerate; where the truncation
makes one of them degenerate only `(4b)` is available.  Nothing below uses
`(4a)` in a proof — it appears only in section 2.2, to say what a Dobrushin
argument would have to control.

---

## 2. Single-site tensorisation is false, and the reason is explicit

### 2.1 The size of `delta_ij` in the paper's own quantities

Write `G_ij(gamma) = ||(Phi_{gamma u {i,j}} - Phi_gamma)Y||^2`, so that
`D_{gamma u {i,j}} = D_gamma - G_ij(gamma)`.  Dropping the two factors
`1 - g_i/D_gamma <= 1` and `1 - g_j/D_gamma <= 1` from (3),

```
    delta_ij(gamma)  <=  (n/2) log{ D_gamma / D_{gamma u {i,j}} }
                     <=  (n/2) . G_ij(gamma) / D_{gamma u {i,j}} .                 (5)
```

**Proposition A (pair amplification).**  Let `u = (I-Phi_gamma)X_i`,
`v = (I-Phi_gamma)X_j`, and `r = <u,v>/(||u|| ||v||)`.  Then

```
    G_ij(gamma) = { t_i^2 - 2 r t_i t_j + t_j^2 } / (1 - r^2),
    t_i = <u, (I-Phi_gamma)Y>/||u||,   g_i(gamma) = t_i^2,     (same for j)
```

and under Assumption 3.2 at level `s = s_0 + s*`,

```
    1 - |r|  >=  nu,     hence     G_ij(gamma)  <=  2 { g_i(gamma) + g_j(gamma) } / nu .   (6)
```

*Proof.*  By YWJ's Lemma 6 (`eq:underfit` in `bvs3_app.tex`) applied to the
nested pair `gamma subset gamma u {i,j}` of `M(s)`, the matrix
`M = n^{-1} X_{ {i,j} }^T (I - Phi_gamma) X_{ {i,j} }` satisfies
`M >= nu I_2`; its diagonal entries are `n^{-1}||(I-Phi_gamma)X_i||^2 <= 1`
because `||X_i||^2 = n`.  Hence the correlation matrix
`[[1,r],[r,1]] = diag(M)^{-1/2} M diag(M)^{-1/2}` has smallest eigenvalue
`>= nu / lambda_max(diag M) >= nu`, i.e. `1 - |r| >= nu` and
`1 - r^2 >= nu`.  The displayed formula for `G_ij` is the rank-two projection
written in the basis `(u/||u||, v/||v||)`, and
`t_i^2 - 2 r t_i t_j + t_j^2 <= (|t_i|+|t_j|)^2 <= 2(t_i^2+t_j^2)`.  ∎

Combining (5), (6) and the residual floor `D >= n sigma_0^2 / 8`
(`lemma:residual`, at level `s_0 + 2`):

```
    delta_ij(gamma)  <=  8 { g_i(gamma) + g_j(gamma) } / ( nu sigma_0^2 ) ,        (7)
```

and, crudely, using `g_i <= D_gamma` and `lemma:Dupper`,

```
    max_{i,j,gamma} delta_ij  <=  32 ||X_S beta*_S||^2 / (nu sigma_0^2)  +  128 n / nu .   (8)
```

**This is of order `n/nu`, and it is not slack.**  In the cancelling-pair
design (`beta* = (b,-b)`, `b = 3`, `r = 0.9`, `n = 100`) the exact maximum of
`delta_ij` over `i, j in S` is `+44.2`, and its *minimum* is `+41.5`: the
log-gain of either covariate is larger by more than `41` in every context that
already contains its partner.  Proposition A's bound `Xi <= 2/nu` is checked
in `tensor_numerics.txt`: `Xi = 8.48` against `2/nu = 21.7` there, and
`Xi <= 1.32` against `2/nu` between `3.4` and `47` on the equicorrelated
designs.

### 2.2 What a Dobrushin argument would need, and why it cannot have it

Dobrushin's condition `alpha = max_i sum_{j != i} R_ij < 1` implies (by path
coupling in the Hamming metric, and the Lipschitz-contraction argument of
section 3.2 below) `Gap(Glauber) >= (1-alpha)/p`, hence single-site
tensorisation with `C = 1/(1-alpha)`.  Splitting the row sum by (4a)-(4b),

```
    alpha  <=  p . max_{i in S^c, gamma} mu[gamma_i=1|gamma_{-i}]
               +  (lambda/4) max_{i in S} sum_{j != i} max_gamma |delta_ij(gamma)| .
```

The first term is `<= p^{1-eta} <= 1/2` by Assumption 6.1 — free.  The second
requires `lambda max_{i in S, j} |delta_ij| <~ 1/p`.  By (8) that asks for
`lambda n / nu <~ 1/p`, i.e. for the target to be a product measure to within
`O(1/p)` in the exponent.  **No hypothesis of the paper is anywhere near
this, and none should be: the log-gain of an influential covariate genuinely
depends on which of its partners are already in the model, by an amount of
order `n`, and that is the phenomenon `sub:twosided` and `rem:third` are
about.**

### 2.3 What is measured

`C_atv = sup_f Var(f) / sum_j E[Var(f | gamma_{-j})]` is exactly
`1/(p Gap(Glauber))` and is computed exactly.  Maxima over `lambda in [0,1]`:

| design | `nu` | `max_{i,j in S} delta_ij` | `Xi_S` | `max_lambda C_atv` | `min_lambda 4p Gap` |
|---|---|---|---|---|---|
| equicorr `rho=0`    | 0.594 | +12.70 | 0.974 | 1.2 | 1.005 |
| equicorr `rho=0.5`  | 0.383 | -2.60  | 0.788 | 2.0 | 0.931 |
| equicorr `rho=0.85` | 0.120 | -6.32  | 0.715 | 27.0 | 0.907 |
| equicorr `rho=0.95` | 0.042 | -3.50  | 0.689 | 95.2 | 1.005 |
| **cancelling `r=0.9`**  | 0.092 | **+44.24** | **8.48** | **701.9** | **0.0016** |
| **cancelling `r=0.95`** | 0.046 | **+26.52** | **10.04** | **2237.9** | **0.0006** |

`C_atv` is not `O(1)`, it is not bounded in `nu`, and it is the right size in
the cancelling pair.  **Single-site approximate tensorisation of variance is
false for this target.**  (Entropy is worse: `C_ent >= C_atv` always.)

Two remarks on the table.  First, `Xi_S < 1` in every equicorrelated design —
the *explained energy* `Y^T Phi_gamma Y` is submodular on pairs inside `S`
there — while `Xi_S ~ 10` in the cancelling pair.  Second, submodularity of
`log L_n` itself is **not** the discriminating property: at `rho = 0` the
mixed second difference reaches `+12.7` and the chain is nonetheless as fast
as it ever gets.  The reason is structural: `-(n/2) log D` is a *convex*
increasing function of the explained energy, so `log L_n` acquires increasing
returns even when the explained energy has diminishing ones.  What separates
the two regimes in every experiment is `Xi_S`, not the sign of `delta_ij`.

---

## 3. What does hold: block factorisation with `S` as one block

### 3.1 The statement

Let `V_0 = S` and `V_j = {j}` for `j in S^c`, so `N = p - s* + 1` blocks.

> **Theorem T.**  Let `lambda in [0,1]`, assume Assumption 6.1 with an exponent
> `eta` satisfying `p^{1-eta} <= 1/4`, and work on the good event.  Then for
> every `f : M(s_0) -> R`,
>
> ```
>     Var_mu(f)  <=  4 { sum_{j in S^c} E[ Var(f | gamma_{-j}) ]
>                        +  E[ Var(f | gamma_{S^c}) ] } ,                   (T1)
> ```
>
> and consequently, with `kappa_S(lambda)` the worst fibre gap defined below,
>
> ```
>     Gap(P_lambda)  >=  1 / { 4 ( 4p + 1/kappa_S(lambda) ) } .             (T2)
> ```

Here

```
    kappa_S(lambda) = min over sigma subset of S^c, |sigma| <= s_0, of
                      inf_f  E^sigma_{P_lambda}(f,f) / Var_{mu(. | gamma_{S^c} = sigma)}(f),
```

`E^sigma_{P_lambda}` being the part of the chain's own Dirichlet form carried
by the edges that stay inside the fibre `{gamma : gamma n S^c = sigma}`,
normalised by the fibre mass.  A fibre has at most `2^{s*}` states, and the
edges available inside it are single flips of `S`-coordinates (rate `1/(4p)`)
and exchanges between two `S`-coordinates (rate `1/(4|gamma|(p-|gamma|))
>= 1/(4 p s_0)`).

The only hypothesis is Assumption 6.1, and only through

```
    eps_1 := max_{j in S^c}  max_{gamma_{-j}}  mu[ gamma_j = 1 | gamma_{-j} ]
          <= p^{-kappa_lambda} <= p^{-eta} ,     required:  p eps_1 <= 1/4.     (9)
```

The paper already assumes `p^{1-eta} <= 1/2` in `cor:precM` and `th:gapM`;
`<= 1/4` is the same condition with `log 4` in place of `log 2` in
`eq:etamin`.  Measured, `p eps_1` lies between `0.045` and `0.071` on every
design and temperature — a factor `4` inside the requirement.

### 3.2 Proof sketch

*Step 1 (block dynamics).*  Let `P_blk = N^{-1} sum_k P_{V_k}`, where `P_{V_k}`
resamples `gamma_{V_k}` from `mu( . | gamma_{V_k^c})`.  Its Dirichlet form is
`E_blk(f) = N^{-1} sum_k E[ Var(f | gamma_{V_k^c}) ]`, so a bound
`Gap(P_blk) >= (1-alpha)/N` is exactly (T1) with constant `1/(1-alpha)`.

*Step 2 (weighted path coupling).*  Give block `V_0 = S` the weight `c_0 = 1/2`
and each singleton the weight `1`, and let
`d(gamma,gamma') = sum_k c_k 1[gamma_{V_k} != gamma'_{V_k}]`.
`d` **is** the path metric of the "differ in one block" adjacency on `M(s_0)`:
`M(s_0)` is downward closed, so performing the deletions before the additions
takes any `gamma` to any `gamma'` one block at a time without the size ever
exceeding `max(|gamma|,|gamma'|) <= s_0`.  For `gamma, gamma'` differing only
in block `k`, couple as follows.  With probability `1/N` the chain picks block
`k`; the two conditionals then coincide, because the states agree off `V_k`, so
the update can be made identical and the discrepancy is removed (`-c_k`).  For
each `l != k`, with probability `1/N` the chain picks block `l`; the two
conditionals are then within `R_lk` in total variation, so a maximal coupling
leaves them disagreeing with probability at most `R_lk`, adding `c_l`.  Hence

```
    E[ Delta d ]  <=  N^{-1} ( - c_k + sum_{l != k} c_l R_lk ),
    alpha := max_k  c_k^{-1} sum_{l != k} c_l R_lk .
```

*Step 3 (the entries).*  By (4b), `R_{j,k} <= eps_1` for every `j in S^c` and
every `k` — including `k = S`, and including the case where the truncation
makes a conditional degenerate.  For `k = {j}` and `l = S` use the trivial
`R_{S,j} <= 1`.  The two column sums are then

```
    k = S:      sum_{j in S^c} 1 . R_{j,S}          <=  p eps_1        (need < c_0 = 1/2),
    k = {j}:    c_0 R_{S,j} + sum_{i in S^c} R_{i,j} <=  1/2 + p eps_1  (need < 1),
```

so `p eps_1 <= 1/4` gives `alpha <= max{ (1/4)/(1/2), 3/4 } = 3/4` and
`Gap(P_blk) >= 1/(4N)`.  (Contraction of `W_1` in a path metric at rate `rho`
bounds every non-trivial eigenvalue by `1-rho`, because
`||P f||_Lip <= (1-rho) ||f||_Lip` and an eigenfunction of a non-constant
eigenvalue has positive Lipschitz seminorm.)  This is (T1).

*Step 4 (back to the chain).*  The first term of (T1) is at most
`4p E_{P_lambda}(f)` by (F1)/(2).  For the second, the definition of
`kappa_S` gives `Var_{mu(.|sigma)}(f) <= kappa_S^{-1} E^sigma(f)` on every
fibre; multiplying by `mu(fibre_sigma)` and summing over `sigma` turns the
right-hand side into `kappa_S^{-1} sum_sigma sum_{e inside fibre_sigma}
Q_lambda(e) (grad f)^2`, and the fibres being disjoint this is at most
`kappa_S^{-1} E_{P_lambda}(f)`.  (T2) follows.  The two families of edges used
are themselves disjoint — a flip of `j in S^c` leaves its fibre, an edge inside
a fibre does not touch `S^c` — so no edge is charged twice, though the argument
does not need this.  ∎

### 3.3 Reading the constant

`kappa_S` is proportional to `1/p` by construction — the fibre edges carry the
chain's own rates — so write `theta_S(lambda) = 4 p kappa_S(lambda)`, a pure
shape constant of the conditional law of `gamma_S`.  Then (T2) reads

```
    Gap(P_lambda)  >=  theta_S / { 16 p (1 + theta_S) } .                   (T3)
```

Three regimes.

* `theta_S = Omega(1)`: `Gap = Omega(1/p)`, the observed order, with an
  absolute constant.  Measured, `theta_S >= 0.35` on all four equicorrelated
  designs at all temperatures.
* `theta_S >= c/s_0`, which is the generic benign case (the binding fibre
  constraint is an exchange at a saturated state, of rate
  `1/(4 s_0 (p-s_0))`): `Gap = Omega(1/(p s_0))`.  This is `eq:gapM0` without
  the `(2s*+1)` and without needing `k_lambda = 0`.
* Under a quantitative version of Assumption 6.3 — unimodality of the *fibre*
  measure on `2^S` rather than of the block measure on `2^{W_lambda}`, with the
  ascent routing of `lemma:ascent` and its congestion `<= 2^{s*}` and length
  `<= 2^{s*}` — Sinclair's bound on the fibre gives
  `kappa_S >= 1/{4 s_0 (p-s_0) 4^{s*}}` and hence
  `Gap = Omega( 1 / (p s_0 4^{s*}) )`: the same order as `eq:gapM`, better by a
  factor `~ (s*)^2`, but *not* an asymptotic improvement.  (This third bullet
  is an estimate I have not checked numerically; the first two are.)

So the route does **not** beat `th:gapM` when `kappa_S` is estimated by
canonical paths.  What it does is *localise* the whole difficulty: the
`p - s*` null coordinates are handled optimally and once and for all, and
everything that can go wrong is confined to a chain on `2^{s*}` states.

---

## 4. Numerical check

`p = 14, n = 100, s_0 = 4, kappa = 2, alpha = 1.5`, `21` temperatures, six
designs (`tensor_numerics.txt`).  Over all `126` (design, temperature) pairs:

| quantity | range | claim |
|---|---|---|
| `C_blk` (exact) | `[1.000, 1.023]` | `<= 4` proved in Theorem T |
| `p eps_1` | `[0.045, 0.071]` | `<= 1/4` required by (9) |
| Dobrushin `alpha` (exact, weighted) | `[0.532, 0.554]` | `<= 3/4` proved in Step 3 |
| `cert4 / Gap`, `cert4 = 1/{4(4p+1/kappa_S)}` | `[0.044, 0.219]` | must be `<= 1` |
| `C_atv` (exact) | `[1.00, 2237.94]` | unbounded — section 2 |
| marginal boundedness `b` | `[2.1e-15, 5.1e-3]` | see section 5 |
| `lambda_max(Psi)` (spectral independence) | `[0.000, 1.000]` | see section 5 |

Per design, worst over `lambda`:

| design | `min Gap` | `4p min Gap` | `min theta_S` | `max C_atv` | `max C_blk` | `min cert4/Gap` |
|---|---|---|---|---|---|---|
| equicorr `rho=0`    | 1.795e-2 | 1.005 | 0.458 | 1.2 | 1.014 | 0.076 |
| equicorr `rho=0.5`  | 1.662e-2 | 0.931 | 0.350 | 2.0 | 1.010 | 0.047 |
| equicorr `rho=0.85` | 1.619e-2 | 0.907 | 0.350 | 27.0 | 1.012 | 0.044 |
| equicorr `rho=0.95` | 1.795e-2 | 1.005 | 0.583 | 95.2 | 1.005 | 0.053 |
| **cancelling `r=0.9`**  | **2.91e-5** | **0.0016** | **0.0014** | 701.9 | 1.018 | 0.085 |
| **cancelling `r=0.95`** | **1.12e-5** | **0.0006** | **0.0004** | 2237.9 | 1.023 | 0.076 |

Read the last three columns together.  **`C_blk` does not move** — the block
factorisation is just as good in the collapsed case as in the benign one,
because the collapse is *inside* the `S` block, which the factorisation does
not attempt to break.  `theta_S` moves by a factor `10^3`, in lock step with
the true gap, and the certificate follows it.  That is the sense in which the
constant "degrades in the cancelling case": not `C_blk`, but the second term
`1/kappa_S` of (T2), which is the term the theorem leaves open.

**Acceptance test.**  The cancelling pair was reproduced independently before
anything else was computed (`beta* = (3,-3)`, null columns residualised against
`span(X_S)`, `40` temperatures): `min Gap = 3.13e-5` at `r = 0.9` and
`1.14e-5` at `r = 0.95`, i.e. factors `571` and `1565` below `1/(4p)`, the
minimum at `lambda = 1`.  These are the same phenomenon as `rem:gapM`'s
`1.6e-5 / 1120` and `2.7e3`; the small discrepancy is the design realisation
(which columns are residualised, and the noise draw), not the effect.

**Scaling in `p`.**  Separately, with `s_0 = 3`, `s* = 2`, equicorrelated at
`rho in {0, 0.85}` and `p = 8, 10, 12, 14, 16` (scratch script, same code
paths):

```
    4p Gap      in [0.90, 1.29]     -- the gap really is Theta(1/p)
    theta_S     in [0.42, 1.04]     -- and theta_S does not drift with p
    cert4/Gap   in [0.058, 0.129]   -- the certificate is a stable fraction of it
```

So (T3) is `Theta(1/p)` with a constant that does not degrade as `p` grows, in
the regime where `theta_S = Omega(1)`.

---

## 5. The hard constraint `|gamma| <= s_0`

**Not fatal, and in this argument free.**  Three separate points.

1. **It does not obstruct the block factorisation.**  `M(s_0)` is downward
   closed (it is the family of independent sets of the uniform matroid
   `U_{s_0,p}`).  Two consequences are used, both in Step 2 of section 3.2:
   the weighted Hamming distance is still the *path* metric of the
   block adjacency, because deletions can be performed before additions; and
   the constraint can only push a conditional probability `mu[gamma_i=1|.]`
   *down*, to `0`, so the bound `(4b)` — the only bound the argument uses on
   the `S^c` rows — is untouched.  The influence the constraint creates is
   negative and is absorbed for free.
2. **Empirically it costs nothing.**  At `lambda = 0.001` the target is
   essentially the truncated product prior, and the exact single-site constant
   is `C_atv = 1.000` on every design.  The reason is that under Assumption 6.1
   the number of null coordinates in the model is stochastically dominated by a
   sum of Bernoulli's with `p eps_1 <= 1/4`, so
   `mu[ |gamma n S^c| >= m ] <= (p eps_1)^m <= 4^{-m}` and the truncation binds
   with probability at most `4^{-(s_0-s*)}`; Assumption 3.4 already forces
   `s_0 > s*`.
3. **Down-up / matroid results do apply formally, and fail for three separate
   reasons.**  Adding `s_0` dummy elements homogenises `M(s_0)` into the bases
   of a rank-`s_0` matroid, and the move set of the down-up walk on that
   complex coincides with the flips and exchanges of `eq:neigh` (remove a real
   element and add a dummy = deletion, remove a dummy and add a real =
   addition, real for real = exchange); the two chains differ only in that the
   down-up walk resamples the added element from the full conditional whereas
   `eq:MH` proposes one uniform candidate and Metropolises.  But:
   * *(a)*  The log-concavity theorems (Anari--Liu--Oveis Gharan--Vinzant, and
     the `0`-spectral-independence of matroids that follows) require the
     weights to be a **product** measure on the matroid, or more generally the
     generating polynomial to be log-concave.  `log L_n(Y|.)` is not additive,
     and nothing in Assumptions 3.1--3.5 makes the polynomial log-concave; the
     cancelling pair shows it is not.
   * *(b)*  Every quantitative local-to-global theorem that turns spectral
     independence into an `Omega(1/#sites)` gap (Chen--Liu--Vigoda and
     successors) also requires **marginal boundedness**: every non-degenerate
     marginal bounded below by a constant `b`.  Here `b` is measured down to
     `2.1e-15`, and this is not an artefact — the prior charges `p^{-kappa}`
     per covariate and Theorem 3.1 says the posterior concentrates, so `b` is
     at best `p^{-kappa}` and generally exponentially small in `n`.  The
     hypothesis fails by a margin no constant can absorb.
   * *(c)*  Even granting (a) and (b), the down-up walk's gap is `O(1/s_0)`
     intrinsically, and comparing it to the one-candidate Metropolis chain of
     `eq:MH` costs a further factor `p`.  The best this route could yield is
     `Omega(1/(p s_0))`, which (T2) already gives without any of it.

   **A fourth, decisive, empirical point:** spectral independence *holds* here
   and is blind to the obstruction.  `lambda_max(Psi)`, `Psi` the pairwise
   influence matrix, is at most `1.000` on every design and every temperature —
   including the cancelling pair, where it is `0.998--1.000` while the gap is
   three orders of magnitude below `1/(4p)`.  A framework whose only input is
   `lambda_max(Psi) <= eta` cannot distinguish the fast case from the collapsed
   one, so it cannot be the right framework for this model.

---

## 6. What is missing, precisely

Theorem T reduces the problem, with no loss beyond a factor `4`, to bounding

```
    kappa_S(lambda) = min_sigma Gap( P_lambda restricted to {gamma n S^c = sigma} ),
```

a spectral gap of a chain on at most `2^{s*}` states.  **Assumptions 3.1--3.5
do not bound it.**  The exact place where they stop is
`lemma:forward`\ref{it:fs_a}: it produces a gain floor for the *greedy*
addition only,

```
    ||(Phi_{gamma u j_gamma} - Phi_gamma) X_S beta*_S||^2  >=  n nu^2 betabar(gamma* \ gamma),
```

which certifies that *one* direction out of every underfitted state is good.
A conductance bound needs that *every* pair of `S`-configurations be joined by
a route that never dips.  The cancelling pair is precisely a design where
`lemma:forward` is **tight and still useless**.  At `gamma = {}` it certifies a
gain of `A^2 = nu ||X_S beta*_S||^2 / s*`, i.e. a `nu/s* = 0.0462` share of the
signal energy at `r = 0.9`; the best single covariate delivers a share of
`0.0461`, so the lemma is attained to three digits.  And that share is not
enough to beat the Occam factor: measured,

```
   design                psi_j(empty)        psi_j(S \ j)     kappa log p = 5.28
   equicorr rho = 0.95   84.3, 83.4, 83.1    -1.2, -1.9, -2.4
   cancelling r = 0.9    -0.63, -3.68        40.9, 37.8
```

At `r = 0.9` *both* influential covariates have a **negative** log-gain in the
empty context and a log-gain above `37` once the partner is in.  The
equicorrelated design has exactly the opposite ordering, which is why it is
fast.  `lemma:forward` cannot see the difference: it bounds the gain from
below, and the bound is met.  So no strengthening of `nu` alone can close the
gap; the missing hypothesis is genuinely about the *joint* behaviour of pairs.

The sharpest candidate I found, and the one that separates the regimes in every
experiment, is the **pair amplification** of Proposition A:

```
    Xi  :=  max over gamma, i != j in S  of
            ||(Phi_{gamma u {i,j}} - Phi_gamma)Y||^2 / { g_i(gamma) + g_j(gamma) } .
```

`Xi <= 2/nu` always (Proposition A, from Assumption 3.2 alone); `Xi <= 1` is
"the explained energy is submodular on pairs inside `S`".  Measured
`Xi_S = 0.97, 0.79, 0.71, 0.69` on the four equicorrelated designs and
`8.48, 10.04` on the two cancelling pairs.

**But `Xi_S` is not a reliable predictor, and I am not proposing it as the
hypothesis.**  `python3 tensor_numerics.py --screen` draws 60 random designs at
`p=12, s_0=3, s*=3` (random `3x3` correlation matrix for `X_S`, random signs
and magnitudes for `beta*`, null columns residualised) and computes `Xi_S`,
`theta_S` and `4p Gap` exactly:

| `Xi_S` bucket | # | min `theta_S` | median `theta_S` | min `4p Gap` | max `C_atv` |
|---|---|---|---|---|---|
| `[0, 1.05)`   | 8  | 0.311 | 0.444 | 0.373 | 2690 |
| `[1.05, 1.3)` | 19 | 0.146 | 0.444 | 0.274 | 725 |
| `[1.3, 2.0)`  | 15 | 0.010 | 0.444 | 0.016 | 4316 |
| `[2.0, inf)`  | 18 | 0.014 | 0.445 | 0.015 | 173 |

(`C_atv` reaching `4316` in a bucket whose median `theta_S` is `0.444` is a
second, independent demonstration that the single-site constant is the wrong
object: it is large in designs whose gap is perfectly healthy.)

Every one of the 8 designs whose gap collapsed (`4p Gap < 0.2`) had
`Xi_S >= 1.40`, and no design with `Xi_S <= 1.05` collapsed; but the rank
correlation between `Xi_S` and `theta_S` over all 60 is only `0.07`, because
most designs are fast whatever `Xi_S` does.  So the evidence supports
"`Xi_S` small is *sufficient*" on a sample of 8, and refutes "`Xi_S` large is
informative".  That is not enough to build a hypothesis on.

**I could not prove that `Xi_S <= 1` implies `kappa_S >= c/(p s_0)`,
and I do not believe the implication in that clean form**: `Xi_S <= 1` bounds
the second difference of the *explained energy*, whereas the barrier is a
second difference of `log L_n`, and the two differ by the convexity of
`-(n/2) log D` — at `rho = 0` the explained energy is submodular (`Xi_S <= 1`)
and `delta_ij` is nevertheless `+12.7`.  What `Xi_S <= 1` does give, from (3)
and `log(1+x) <= x`, is the sharper bound

```
    delta_ij(gamma)  <=  (n/2) . g_i(gamma) g_j(gamma)
                             / { D_gamma ( D_gamma - g_i(gamma) - g_j(gamma) ) } ,
```

the denominator being positive because `D_gamma - g_i - g_j >= D_gamma - G_ij
= D_{gamma u {i,j}} > 0` under `Xi_S <= 1`.  This is still `O(n)` and therefore
still useless for a Dobrushin argument, but it is the exact statement of how
much increasing return survives when the energies do not amplify: it is
second order in the individual gains, whereas (5)-(7) is first order.

So the honest list of what would have to be assumed, in decreasing strength:

1. `theta_S(lambda) = 4 p kappa_S(lambda) >= c` for an absolute `c`.  Gives
   `Gap = Omega(1/p)`, the observed order.  This is the assumption; it is
   directly checkable by enumeration on `2^{s*}` states, it holds at `>= 0.35`
   in every benign experiment, and it fails at `4e-4` in the cancelling pair.
2. Assumption 6.3 (unimodal fibres) plus the ascent routing of
   `lemma:ascent` applied inside a fibre.  Gives `Gap = Omega(1/(p s_0 4^{s*}))`
   — same order as `th:gapM`, not better.
3. Nothing beyond Assumptions 3.1--3.5.  Gives nothing, and here is the sharp
   form of that statement.  `nu` and `omega(X)` **cannot** separate the two
   regimes, because the collapsed design is the *better* one on both:

   | design | `nu` | `omega(X)` | `4 nu^{-2} omega + 1` | `min Gap` |
   |---|---|---|---|---|
   | equicorr `rho=0.85` | 0.120 | 1.53 | 424 | 1.62e-2 |
   | equicorr `rho=0.95` | 0.042 | 1.83 | 4102 | 1.80e-2 |
   | **cancelling `r=0.9`** | **0.092** | **0.82** | **388** | **2.91e-5** |
   | **cancelling `r=0.95`** | 0.046 | 0.91 | 1695 | 1.12e-5 |

   The cancelling pair at `r = 0.9` has a *larger* restricted eigenvalue and a
   *smaller* `omega(X)` than the equicorrelated design at `rho = 0.85`, and
   both quantities enter `ass:B` and `ass:D` monotonically — yet its gap is
   `550` times smaller.  Its signal is also well inside `ass:A`'s cap
   (`||X_S beta*_S||^2 = 166` against `g sigma_0^2 log p = 7242`) and its null
   columns are orthogonal to `span(X_S)`, so `Ltilde = 0`.  (At `p = 14`,
   `n = 100` neither family meets the numerical constants of `ass:C` and
   `ass:betamin` — those are asymptotic and no `p = 14` example meets them,
   the benign designs included; the comparison above is between designs that
   stand or fall together on that score.)  **A hypothesis that excludes the
   cancelling pair must therefore be of a different kind from `ass:B` and
   `ass:D`: it has to constrain pairs of covariates jointly, not the spectrum
   of `X_gamma`.**

---

## 7. What I would take from this

* **Do not claim, and do not conjecture, `Gap = Theta(1/p)` unconditionally.**
  Section 6's `rem:gapM` is already careful about this via `ass:unimodal`; the
  present notes add that the *tensorisation* route runs into exactly the same
  wall, at exactly the same place, and cannot be the escape `rem:gapM` gestures
  at ("A bound of that order would have to come from approximate tensorisation
  rather than from canonical paths").  It would be more accurate to say that
  approximate tensorisation removes the `p - s*` null coordinates optimally and
  leaves the `2^{s*}`-state problem untouched.
* **Theorem T is worth stating**, not because its constant beats `eq:gapM` —
  it does not, once `kappa_S` is bounded by canonical paths — but because it is
  the only bound in the file that is within a small constant factor of the
  truth on every design tested, and because it isolates the obstruction
  exactly.  It also has no `n`, no `lambda`, no `s*` and no `4^{s*}` in it: all
  of those sit inside `kappa_S`.
* **The exchange move is load bearing here too.**  Inside a fibre, two
  configurations of `S` of the same size are joined only by an exchange; delete
  the exchange move and `kappa_S` collapses with the flip-only gaps of
  `exchange_test.txt`.
* If anything is to be added to the draft, the cheapest useful item is the
  identity (3) together with Proposition A: they say, in the paper's own
  quantities and in one line each, exactly what "context dependence of
  `psi_j`" is and how large Assumption 3.2 permits it to be.  `rem:third` and
  `sub:twosided` currently describe that phenomenon qualitatively.


---

## 8. Reproducing

```
    python3 tensor_numerics.py            ->  tensor_numerics.txt   (~10 min)
    python3 tensor_numerics.py --quick    ->  9 temperatures instead of 21
    python3 tensor_numerics.py --screen   ->  the 60 random designs of section 6
```

`tensor_numerics.py` imports `window_numerics` for `enumerate_states`,
`log_marginal` and the model constants (`p = 14, n = 100, s_0 = 4, kappa = 2,
alpha = 1.5`) and builds everything else itself, so `s*` and `beta*` are free
parameters; the two cancelling designs are `s* = 2`, `beta* = (3, -3)` at
`rho = 0.9` and `0.95`.  The `p`-scaling table of section 4 and the
`nu`/`omega(X)` table of section 6 were produced by short scratch scripts using
the same functions; they are quoted, not re-derived, and can be regenerated by
calling `restricted_eig`, `fibre_gap` and `form_gap` at other `(p, s_0, s*)`.
