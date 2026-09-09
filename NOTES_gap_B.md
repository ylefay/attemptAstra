# A lambda-dependent gap bound: the undecided window as a single block

Working notes.  Every number below is produced by one of
`gapB_numerics.py` (-> `gapB_numerics.txt`, `--big` -> `gapB_numerics_big.txt`),
`gapB_checks.py` (-> `gapB_checks.txt`), `gapB_adversarial.py`
(-> `gapB_adversarial.txt`) and `gapB_witness.py` (-> `gapB_witness.txt`), all
of which import `window_numerics` and `flow_numerics` read-only and build the
rest in `gapB_core.py`.  Nothing in `binary_variable_selection3.tex`,
`bvs3_app.tex` or any pre-existing `.py` has been touched.  Nothing is
asserted that has not been computed.  Statements are tagged
**[proved]**, **[computed]** or **[conjecture]**.

---

## 0.  Answer, up front

**Section 6 should not be a canonical-path argument at all.**  The right
object is not a transition function with a centre that moves with `lambda`,
but a *partition of the coordinates* that moves with `lambda`: at temperature
`lambda` the chain's `p` coordinates split into those the temperature has
already decided — every null covariate, and every influential covariate whose
log-gain has not yet reached the sparsity penalty — and the handful
`m^+_lambda ⊆ gamma*` it has not.  The decided coordinates are pinned, in
*every* context, to within `p^{-2}` or `1/(8s*)` of `0`; a weighted
Dobrushin/path-coupling argument over the partition
`{ m^+_lambda } u { singletons off m^+_lambda }` then factorises the variance
with an absolute constant, and the *whole* residual difficulty is the spectral
gap of the conditional law of `gamma_{m^+_lambda}` given everything else — a
chain on at most `2^{s*}` states.  This gives **[proved]**

```
    Gap(P_lambda)  >=  (1 - alpha_lambda) . min{ 1/(4p) ,  kappa^+_lambda }
                   >=  (1/3) min{ 1/(4p), kappa^+_lambda }       (p >= 8),
```

with `alpha_lambda` an explicit absolute constant `<= 0.6404` and
`kappa^+_lambda` the worst fibre gap.  Besides `p >= 8` and the good event on
which it is stated, it needs **no hypothesis beyond Assumption 6.1**
(`ass:twosided`) — not `ass:unimodal`, not Assumptions 3.1--3.6, no hypothesis
on the modal structure of the target.  In (B) there is no `2^{k_lambda}`, no
`Theta_lambda`, no `(1+k Theta)^k` and no `4^{s*}`; `m^-_lambda`, the greedy
chain of `def:chain` and the transition function `G` of `def:G` do not appear
at all.

**What it costs.**  `kappa^+_lambda` is bounded below by canonical paths on the
fibre, giving a fully explicit `B(lambda)`.  Over 231 exactly-computed
(design, temperature) pairs — the four equicorrelated designs at `s* = 3` and
at `s* = 5`, and three cancelling pairs — the ratio bound/truth is

| form of the bound | `B/Gap` range | worst factor lost |
|---|---|---|
| `B(lambda)` with the exact `kappa^+_lambda` | `[0.0795, 0.4628]` | 12.6 |
| `B(lambda)` with the canonical-path bound on `kappa^+_lambda` | `[0.0061, 0.4617]` | 164 |
| `eq:gapM` (`th:gapM`, = the two-level bound of `NOTES_flow`) | `[1.9e-6, 8.3e-1]` | 5.3e5 |
| `eq:gapM0` (only where the window is empty) | `[0.0019, 0.0062]` | 531 |
| the block-tensorisation certificate of `NOTES_tensorisation` (`s*=3` only) | `[0.0437, 0.2186]` | 22.9 |
| the v2 bound `1/{32 p s_0(2s*+1) max(1,Theta)(1+kTheta)^k}` | `[0, 6.3e-3]` | unbounded |

**The quantity that degrades, and only there.**  It is the *shape congestion*
`varrho_lambda` of the fibre ensemble — equivalently `1/(4p kappa^+_lambda)`, a
bottleneck of the conditional law of the undecided covariates.  **[computed]**

| design | `varrho_lambda` over `{m^+_lambda != 0}` | `4p kappa^+_lambda` |
|---|---|---|
| equicorr `rho=0`, `s*=3` | `[0.443, 1.429]` | `[0.458, 1.128]` |
| equicorr `rho=0.5` | `[0.490, 1.429]` | `[0.350, 1.364]` |
| equicorr `rho=0.85` | `[0.502, 1.428]` | `[0.350, 1.512]` |
| equicorr `rho=0.95` | `[0.502, 0.986]` | `[0.583, 1.930]` |
| equicorr `rho=0`, `s*=5` | `[0.337, 1.601]` | `[0.208, 0.995]` |
| equicorr `rho=0.85`, `s*=5` | `[0.906, 1.143]` | `[0.292, 1.701]` |
| **cancelling `r=0.9`, `b=3`** | **`[0.649, 183.9]`** | **`[0.0014, 0.782]`** |
| **cancelling `r=0.95`, `b=3`** | **`[0.574, 628.4]`** | **`[0.0004, 0.870]`** |
| cancelling `r=0.9`, `b=1` (benign) | `[0.499, 0.500]` | `[1.0008, 1.0014]` |

Three orders of magnitude, at the cancelling pair, at the temperatures where
the true gap falls by three orders of magnitude, and nowhere else.

**A sharp negative result on the side.**  `eq:gapM` of `th:gapM` is *false*
without `ass:unimodal`.  `gapB_witness.py` exhibits a design (`p=14`, `n=100`,
`s_0=4`, `s*=2`, a cancelling pair at `r=0.98`, `b=4`) on which
`ass:twosided` holds with `a_0 = 0.668`, the good event
`A_n n B_n n C_n n D_n` holds, `ass:unimodal` fails (two local maxima), and at
`lambda = 1`

```
        Gap(P_1) = 3.4688e-06   <   eq:gapM = 9.3006e-06        (factor 2.681)
```

Theorem B holds there with a margin of 3.6.  So `ass:unimodal` is load-bearing
in `th:gapM` and cannot be quietly dropped; over the 120-design search of
`gapB_adversarial.py` the violation reaches a factor 2.68 and occurs in 24 of
the 60 cancelling-*pair* designs, and in none of the 60 cancelling-*triple*
designs, where `eq:gapM` carries an extra `2^{k}` with `k=3`. **[computed]**

**What I could not do.**  I could not bound `kappa^+_lambda` from
`(X, beta*, sigma_0, g, kappa)` by anything better than a fixed power of `p`
(Corollary B4 below, `p^{-(k+1)(kappa + lambda alpha)}`, which is 8 to 16
orders of magnitude looser than the truth), and I do not believe a better
assumption-free bound exists — section 6.  And Theorem B's residual is a
minimum over *all* fibres, including those in which the truncation
`|gamma| <= s_0` binds; those are the binding ones in 16--20 of 21
temperatures on the equicorrelated designs, they cost a factor 1.0--3.2
**[computed]**, and I have no a priori bound for them at all — section 6.2.
Consequently Theorem B does **not** subsume `th:mixing`: the two are
complementary and one should take the maximum.

---

## 1.  Notation and the two definitions the theorem needs

Write `L := log p`, `mu := pi_n^(lambda)(. | Y)`, and use `alpha log p` for
`(1/2)log(1+g)` throughout (`eq:alphadef`; guidance item 3 adopted, it makes
every threshold below transparent).  As in `eq:psidef`,

```
    psi_j(gamma) = log L_n(Y|gamma u {j}) - log L_n(Y|gamma),   |gamma| < s_0,
    log { mu(gamma u j)/mu(gamma) } = lambda psi_j(gamma) - kappa log p .    (eq:onestepG)
```

**Definition 1 (extremes over all contexts).**
```
    psi_j^-  :=  min { psi_j(gamma) : gamma in M(s_0), |gamma| < s_0, j not in gamma },
    psi_j^+  :=  max { same } .
```
This differs from `eq:psipm`, which restricts the context to `gamma ⊆ gamma*`.
The enlargement is needed because the Dobrushin argument is a worst-case
statement over contexts, and it is nearly free: measured, `psi_j^-` is
*identical* to the `eq:psipm` value on every design and temperature computed,
and `psi_j^+` exceeds it by `0.48` to `4.13` (e.g. `rho=0.85`:
`6.335 / 78.089` against `6.335 / 74.713`).
**[computed]**

**Definition 2 (decided and undecided, margin `log(8s*)`).**
```
    m^-_lambda := { j in gamma* : lambda psi_j^- >= kappa log p + log(8 s*) },
    m^+_lambda := { j in gamma* : lambda psi_j^+ >  kappa log p - log(8 s*) },
    W_lambda   := m^+_lambda \ m^-_lambda,  k_lambda := |W_lambda| .
```
This is `def:cpm` with `log(2s*)` replaced by `log(8s*)` and without the
"longest prefix of the greedy chain" restriction (dropping the prefix rule can
only enlarge `m^-` and shrink `m^+`, hence only shrink the window).  The
larger margin is what buys the factor `1/8` in `delta-bar` below; it costs
nothing else.

**Definition 3 (fibres and the fibre gap).**  For `sigma ⊆ [p] \ m^+_lambda`
with `|sigma| <= s_0`, let `F_sigma := { gamma in M(s_0) : gamma \ m^+_lambda
= sigma }`; the `F_sigma` partition `M(s_0)`.  Write `mu_sigma := mu(. |
F_sigma)`, let `E^sigma(f) := mu(F_sigma)^{-1} sum_{e ⊆ F_sigma} Q_lambda(e)
(grad f)^2` be the part of the chain's own Dirichlet form carried by the edges
with *both* endpoints in `F_sigma`, and

```
    kappa_sigma := inf_{f non-constant} E^sigma(f) / Var_{mu_sigma}(f),
    kappa^+_lambda := min over sigma with |F_sigma| >= 2 of kappa_sigma ,
```
with `kappa^+_lambda := +infinity` when `m^+_lambda = 0`.  A fibre has at most
`2^{|m^+_lambda|} <= 2^{s*}` states and the edges available inside it are flips
of `m^+_lambda`-coordinates (rate `1/(4p)`) and exchanges between two
`m^+_lambda`-coordinates (rate `>= 1/(4 p s_0)`).

---

## 2.  The theorem

> **Theorem B.**  Let `lambda in [0,1]` and `p >= 8`, assume Assumption 6.1
> (`ass:twosided`) with `kappa_lambda = kappa - lambda(a_0 - alpha) >= 2`, and
> work on the good event.  Put
> ```
>     delta-bar := p^{1-kappa_lambda} + s*/(1 + 8 s*)      ( <= 1/p + 1/8 <= 1/4 ),
>     c_0       := ( sqrt(delta-bar^2 + 4 delta-bar) - delta-bar ) / 2,
>     alpha_lambda := c_0 + delta-bar .
> ```
> Then `alpha_lambda <= 0.6404 < 1` and
> ```
>     Gap(P_lambda)  >=  (1 - alpha_lambda) . min{ 1/(4p),  kappa^+_lambda }.   (B)
> ```
> In particular `Gap(P_lambda) >= (1/3) min{ 1/(4p), kappa^+_lambda }`, and if
> `m^+_lambda = 0` then `Gap(P_lambda) >= (1 - delta-bar)/(4p) >= 3/(16 p)`.

At `p = 14`, `s* = 3`: `delta-bar = 0.1914`, `c_0 = 0.3522`,
`1 - alpha_lambda = 0.4564`.  As `p -> infinity`, `1 - alpha_lambda -> 0.5785`.

### 2.1  Proof

Throughout, `Gap` is `1 - lambda_2`, which by `lemma:kernel` is the absolute
spectral gap because `P_lambda` is non-negative definite; and it suffices to
prove the Poincaré inequality `Var_mu(f) <= C E_{P_lambda}(f,f)` with
`C = (1-alpha_lambda)^{-1} max{4p, 1/kappa^+_lambda}`.

**Step 0 (the pinned conditionals).**  Fix `j` not in `m^+_lambda` and a
context `gamma_{-j}`; write `gamma` for the state with that context and
`gamma_j = 0`.  Then
```
    mu[ gamma_j = 1 | gamma_{-j} ] = { 1 + exp( kappa log p - lambda psi_j(gamma) ) }^{-1}
```
if `|gamma| < s_0`, and `= 0` if `|gamma| = s_0`.

*(i) `j in S^c`.*  **[Assumption 6.1 used here, and only here]**
`psi_j(gamma) <= (a_0 - alpha) log p` for every `gamma in M(s_0)` with
`|gamma| < s_0`, so `kappa log p - lambda psi_j(gamma) >= kappa_lambda log p`
and `mu[gamma_j = 1 | .] <= p^{-kappa_lambda} =: eps_j`.

*(ii) `j in gamma* \ m^+_lambda`.*  By Definition 2,
`lambda psi_j^+ <= kappa log p - log(8 s*)`, so
`mu[gamma_j = 1 | .] <= 1/(1 + 8 s*) =: eps_j`.

In both cases the truncation only replaces the value by `0`, so the bound
survives it.  Hence
```
    delta := sum_{j not in m^+_lambda} sup_{contexts} mu[gamma_j = 1|gamma_{-j}]
          <= p . p^{-kappa_lambda} + s* / (1 + 8 s*)  =  delta-bar .        (1)
```

**Step 1 (block dynamics).**  Blocks `V_0 := m^+_lambda` and `V_j := {j}` for
`j not in m^+_lambda`; `N` = the number of blocks (`<= p`).  Let `P_{V}` be the
kernel that redraws `gamma_V` from `mu(. | gamma_{V^c})` and
`P_blk := N^{-1} sum_k P_{V_k}`.  Each `P_{V_k}` is a conditional expectation,
hence a self-adjoint idempotent on `L^2(mu)`; so `P_blk` is `mu`-reversible and
non-negative definite, and its Dirichlet form is
```
    E_blk(f) = N^{-1} sum_k E_mu[ Var(f | gamma_{V_k^c}) ] .                (2)
```

**Step 2 (weighted path coupling).**  Give `V_0` weight `c_0` and every
singleton weight `1`, and set
`d(gamma,gamma') := sum_k c_k 1[gamma_{V_k} != gamma'_{V_k}]`.

*(a) `d` is the path metric of the "differ in exactly one block" graph on
`M(s_0)`.*  Any block path must touch each differing block at least once, so
its length is at least `d`.  Conversely `M(s_0)` is downward closed, so
performing first every block move that does not increase `|gamma|` and then
every move that does takes `gamma` to `gamma'` one block at a time: in the
first phase the size never exceeds `|gamma| <= s_0`, in the second it never
exceeds `|gamma'| <= s_0`.  (This is the point at which the hard constraint
`|gamma| <= s_0` has to be handled, and it is handled for free; the same
observation is `NOTES_tensorisation` section 5 item 1.)

*(b) Contraction.*  Let `gamma, gamma'` differ only in block `k`.  Pick a block
uniformly.  With probability `1/N` it is `k`: the two states agree off `V_k`,
so the two conditionals of `gamma_{V_k}` are *identical* (the truncation
constrains both in the same way), the update can be made identical, and the
discrepancy is removed, gaining `-c_k`.  With probability `1/N` it is
`l != k`: a maximal coupling of the two conditionals of `gamma_{V_l}` leaves
them disagreeing with probability at most
`R_{lk} := TV( mu(gamma_{V_l} in . | gamma_{V_l^c}), mu(gamma_{V_l} in . |
gamma'_{V_l^c}) )`, costing `+c_l`.  So
```
    E[ Delta d ] <= N^{-1} ( -c_k + sum_{l != k} c_l R_{lk} ) ,
    alpha := max_k  c_k^{-1} sum_{l != k} c_l R_{lk} .
```

*(c) The entries.*  For a singleton `l`, Step 0 says the conditional law of
`gamma_l` is Bernoulli with parameter in `[0, eps_l]` in *every* context, so
`R_{lk} <= eps_l` for every `k`, `V_0` included.  For `l = V_0` use the trivial
`R_{V_0,k} <= 1`.  Hence, by (1),
```
    column k = V_0 :  c_0^{-1} sum_{l singleton} eps_l  <=  delta-bar / c_0 ,
    column k = {j} :  c_0 . 1 + sum_{l != j} eps_l      <=  c_0 + delta-bar ,
```
so `alpha <= max{ delta-bar/c_0, c_0 + delta-bar }`.  The stated `c_0` is the
root of `c^2 + delta-bar c - delta-bar = 0`, which equalises the two, giving
`alpha <= c_0 + delta-bar = alpha_lambda`.  With `p >= 8` and
`kappa_lambda >= 2`, `delta-bar <= 1/4`, `c_0 <= 0.39039` and
`alpha_lambda <= 0.64039`, so `1 - alpha_lambda >= 0.3596 > 1/3`.
If `m^+_lambda = 0` there is no `V_0` block and `alpha <= delta-bar <= 1/4`.

*(d) From contraction to a gap.*  `W_1`-contraction at rate `rho` in a path
metric gives `||P_blk f||_Lip <= (1-rho)||f||_Lip`; an eigenfunction of a
non-constant eigenvalue has positive Lipschitz seminorm, so every non-trivial
eigenvalue of `P_blk` is at most `1-rho`, and `P_blk` being non-negative
definite, `Gap(P_blk) >= rho = (1 - alpha_lambda)/N`.  With (2),
```
    Var_mu(f)  <=  (1 - alpha_lambda)^{-1} sum_k E_mu[ Var(f|gamma_{V_k^c}) ] .   (3)
```

**Step 3 (each term against the chain's own form).**

*(i) Singletons.*  For a single flip, `eq:MH` and `eq:lazy` give
`Q_lambda(gamma, gamma^j) = (4p)^{-1} min{mu(gamma), mu(gamma^j)}` *exactly*.
The conditional given `gamma_{-j}` is carried by the two states
`{gamma, gamma^j}`, so with `a := mu(gamma)`, `b := mu(gamma^j)`,
`Var(f|gamma_{-j}) = ab(grad f)^2/(a+b)^2` and, since `ab/(a+b) <= min(a,b)`,
```
    sum_{j not in m^+} E_mu[Var(f|gamma_{-j})]
        = sum_{j not in m^+} sum_{flip pairs} { ab/(a+b) } (grad f)^2
        <= 4p . E^out(f) ,                                                  (4)
```
`E^out` being the part of `E_{P_lambda}` carried by flips of coordinates
outside `m^+_lambda`.

*(ii) The block.*  By Definition 3,
`Var_{mu_sigma}(f) <= kappa_sigma^{-1} E^sigma(f)`, and
`mu(F_sigma) E^sigma(f) = sum_{e ⊆ F_sigma} Q_lambda(e)(grad f)^2`, so
```
    E_mu[ Var(f|gamma_{V_0^c}) ] = sum_sigma mu(F_sigma) Var_{mu_sigma}(f)
        <= (kappa^+_lambda)^{-1} E^in(f) ,                                  (5)
```
`E^in` being the part of `E_{P_lambda}` carried by edges inside fibres.

*(iii) Disjointness.*  A flip of `j not in m^+_lambda` changes
`gamma \ m^+_lambda`, hence leaves its fibre; an edge inside a fibre changes
only `m^+_lambda`-coordinates.  So the two edge families are disjoint and
`E^out(f) + E^in(f) <= E_{P_lambda}(f)`.  Adding (4) and (5) and using
`x A + y B <= max(x,y)(A+B)` for non-negative `A,B`,
```
    sum_k E_mu[Var(f|.)]  <=  max{ 4p, 1/kappa^+_lambda } . E_{P_lambda}(f) .
```
Substituting into (3) gives the Poincaré inequality, hence (B).  `∎`

### 2.2  Corollaries

> **B1 (empty window).** `m^+_lambda = 0`  =>  `Gap(P_lambda) >= 3/(16 p)`.
> **[proved]**  (Step 2(c), no `V_0` block.)

> **B2 (one possibly-influential covariate).**  `|m^+_lambda| <= 1`  =>  every
> fibre has at most two states, joined by a single flip, so
> `kappa_sigma = 1/(4p max(nu_0,nu_1)) >= 1/(4p)` and
> `Gap(P_lambda) >= 1/(12 p)`.  **[proved]**
> Verified exactly on the `b=1` cancelling design: `4p kappa^+_lambda` between
> `1.0008` and `1.0014` at `lambda in {0.7,0.8,0.9,1.0}` (`gapB_checks.txt`,
> item E). **[computed]**

> **B3 (canonical paths on the fibre).**  For each fibre take a spanning tree
> `R_sigma` of `(F_sigma, flips and exchanges)` rooted at the
> `mu_sigma`-largest state and route each ordered pair along the tree path.
> By Sinclair's Corollary 6' / Diaconis--Stroock,
> `kappa_sigma >= 1/(rho_sigma ell_sigma)` with
> `rho_sigma = max_e Q^sigma(e)^{-1} sum_{T ∋ e} mu_sigma(A) mu_sigma(B)`
> and `Q^sigma(e) = c_e min(mu_sigma(u), mu_sigma(v))`, `c_e = 1/(4p)` for a
> flip and `>= 1/(4 p s_0)` for an exchange.  Writing
> `R_lambda := max_sigma rho_sigma ell_sigma = 4 p s_0 varrho_lambda ell_lambda`,
> ```
>     Gap(P_lambda) >= (1 - alpha_lambda) / ( 4p . max{1, s_0 varrho_lambda ell_lambda} ).
> ```
> `varrho_lambda` is the *shape congestion*; for a tree ensemble
> `varrho_lambda = max_{e in R} 2 mu_sigma(C_e) mu_sigma(C_e^c) /
> min(mu_sigma(u), mu_sigma(v))`, a Cheeger ratio along the tree.  **[proved]**
> Two trees are computed and the better taken: the maximum-bottleneck spanning
> tree (Prim on edge weight `min(mu_u, mu_v)`; by the max-spanning-tree
> property its paths maximise the minimum mass encountered) and the monotone
> tree `parent(A) = A u {min(A^dagger \ A)}` corrected by exchanges.

> **B4 (a priori, uncapped fibres only).**  If `|sigma| <= s_0 - |m^+_lambda|`
> then every subset of `m^+_lambda` is a legal state of `F_sigma`, so the route
> `A -> A u B -> B` (add the elements of `B \ A` in increasing order, then
> delete those of `A \ B`) stays in the fibre and has at most
> `|A xor B| <= |m^+_lambda|` edges.  Since `psi_j(gamma) >= -alpha log p`
> always (adding a column never decreases the fit), every addition costs at
> most `p^{kappa+lambda alpha}`, so `mu_sigma(A) <= mu_sigma(u)
> p^{|u \ A|(kappa+lambda alpha)}` for `A ⊆ u`.  Both legs of the route through
> an edge `e = (u, u u {j})` have `u ⊇ A` (resp. `u ⊇ B`), so the demand across
> `e` is at most `2 sum_{A ⊆ u} mu_sigma(A) <= 2 mu_sigma(u)
> (1 + p^{kappa+lambda alpha})^{|u|}`, while
> `min(mu_sigma(u), mu_sigma(u u {j})) >= mu_sigma(u) p^{-(kappa+lambda alpha)}`.
> Hence
> ```
>     varrho_lambda <= 2 (1 + p^{kappa + lambda alpha})^{|m^+_lambda| + 1},
>     ell_lambda    <= |m^+_lambda| .
> ```
> **[proved]**  This is affine in `lambda` in the exponent, contains no
> `2^{k_lambda}` and no `Theta_lambda`, and is a *fixed power of `p`* — but it
> is `p^{-(s*+1)(kappa+alpha)}` in the worst case, which is weaker than
> `p^{-s*}`, and numerically it is 8 to 16 orders of magnitude above the exact
> congestion (`gapB_checks.txt`, item D: bound `1.6e8`--`2.2e16` against exact
> `0.73`--`628`).  **[computed]**  It is the honest assumption-free content and
> it is not the operative bound.

> **B5 (the cold end).**  Let
> `lambda_cold := (kappa log p - log(8 s*)) / max_{j in gamma*} psi_j^+`.  For
> `lambda <= lambda_cold`, `m^+_lambda = 0` and B1 applies.  On the good event
> `psi_j(gamma) <= -alpha log p + 4||Y||^2/sigma_0^2` for every `j` and every
> `gamma` (this is `eq:over1b`--`eq:over2`, which use only that
> `Phi_{gamma u j} - Phi_gamma` is a rank-one projection and
> `lemma:residual`), and on `C_n`, `ass:A`,
> `||Y||^2 <= 2||X_S beta*_S||^2 + 4 Ltilde sigma_0^2 log p + 6 n sigma_0^2`, so
> ```
>     lambda_cold >= (kappa log p - log 8s*) /
>                    ( 8||X_S beta*_S||^2/sigma_0^2 + 16 Ltilde log p + 24 n ) .
> ```
> **[proved]**  Measured `lambda_cold` is `0.024`--`0.101` on the equicorrelated
> and cancelling designs (`0.660` on the benign `b=1` cancelling pair), and the
> proved lower bound is `7.7e-4`--`3.5e-3` — right in form, 30 to 100 times too
> small (`gapB_checks.txt`, item F).  **[computed]**

> **B6 (decided window, uncapped fibres).**  If `W_lambda = 0` and
> `|sigma| <= s_0 - |m^+_lambda|`, the map `A -> A u {min(m^-_lambda \ A)}` is a
> valid transition function on `F_sigma` with a single fixed point, every step
> increasing `mu_sigma` by a factor at least `8 s*` (Definition 2), at most
> `s*` immediate precedents per state and at most `s*` steps.  Telescoping,
> `mu_sigma[Lambda(A)] <= (8/7) mu_sigma(A)`, every edge is an increasing single
> flip, so `rho_sigma <= (64/7) p` and `ell_sigma <= 2 s*`, whence
> `kappa_sigma >= 7/(128 p s*)` and
> ```
>     Gap(P_lambda)  >=  7 / (384 p s*) .
> ```
> **[proved, uncapped fibres only]**  (Demand across an edge is counted in both
> orientations, the convention the scripts use; counting directed edges as
> `eq:rhoreduce` does would give `7/(192 p s*)`.)  This is the analogue of
> `eq:gapM0 = 1/(32 p s_0 (2s*+1))` and is better than it by a factor
> `32 s_0 (2s*+1) . 7 / (384 s*) ~ 1.17 s_0`.

---

## 3.  Where the lambda-dependence enters, and what it costs

Exactly two places, and both are explicit.

1. **The partition moves with `lambda`.**  `V_0 = m^+_lambda` grows from `0`
   (below `lambda_cold`, Corollary B5) to a subset of `gamma*`.  This is
   guidance item 2 in its strongest form: the construction has *no centre at
   all* on `M(s_0)`; the only ensemble left is inside a fibre of at most
   `2^{s*}` states, and there its root is the fibre mode, which moves with
   `lambda`.  Nothing is directed at `gamma*`, and `prop:sharp` — which says
   that no ensemble directed at `gamma*` can work below `lambda_star` — has
   nothing to bite on.

2. **The fibre measure is a Gibbs measure with `lambda` in front of the
   energy.**  On `F_sigma`,
   ```
       log mu_sigma(A) = lambda { log L_n(Y|sigma u A) - log L_n(Y|sigma) }
                          - kappa |A| log p + const ,      A ⊆ m^+_lambda ,
   ```
   *affine in `lambda`* up to the normalisation.  **[proved]**  Consequently
   every log-ratio inside a fibre, `log mu_sigma(A) - log mu_sigma(w)`, is an
   affine function of `lambda`; so the barrier that drives Corollary B3 is a
   maximum of affine functions of `lambda`, hence convex and piecewise linear.
   `log(1/kappa^+_lambda)` itself is not exactly of that form — it also carries
   the log-partition function of the fibre, which is a log-sum-exp of affine
   functions **[proved]** — but it is dominated by the barrier whenever the
   fibre measure is concentrated, which is why `varrho_lambda` at the
   cancelling pair grows monotonically along the temperature grid from `0.65`
   to `184` (`r=0.9`) and `0.57` to `628` (`r=0.95`) **[computed]**.  This is
   what makes the bound genuinely `lambda`-dependent rather than a
   `lambda`-uniform bound with a threshold: it is the same formula at every
   temperature and it interpolates.

**The cost.**  `B(lambda)` loses a factor at most `12.6` to the truth in its
exact form and at most `164` in its path form, uniformly over the 231 pairs
computed.  Decomposing the loss at `rho = 0.85`, `lambda = 1`, `s* = 3`:
`4p Gap = 0.907`; `1 - alpha_lambda = 0.4564` (the exact Dobrushin `alpha` is
`0.357`, and the exact factorisation constant is `1.0072` against the proved
`2.191`, so this factor is `2.2` and almost all of it is slack in the
Dobrushin bound, not in the factorisation); `4p kappa^+_lambda = 0.350` (a real
factor 2.6, of which 1.9 is the capped fibres, section 6.2); and
`varrho_lambda ell_lambda s_0 = 1.428 x 3 x 4 = 17.1` in the path form.

---

## 4.  Numeric validation

`p = 14`, `n = 100`, `kappa = 2`, `alpha = 1.5`, `g = p^{2 alpha}`,
`sigma_0 = 1`, 21 temperatures on `[0.001, 1]`, `|M(s_0)| = 1471` (`s_0=4`) and
`6476` (`s_0=6`).  Designs as in `flow_numerics.py`.  Full tables in
`gapB_numerics.txt` and `gapB_numerics_big.txt`.

### 4.1  Bound against truth, per design (min over `lambda` of each column)

| design | `min Gap` | `B` exact form | `B/Gap` range | `B` path form | `B/Gap` range | `eq:gapM` | ratio range | cert of `NOTES_tensorisation` | ratio range |
|---|---|---|---|---|---|---|---|---|---|
| equicorr `rho=0`, `s*=3` | 1.795e-2 | 3.73e-3 | `[0.203,0.454]` | 4.75e-4 | `[0.025,0.454]` | 4.43e-6 | `[2.3e-4,1.3e-3]` | 1.40e-3 | `[0.076,0.125]` |
| equicorr `rho=0.5` | 1.662e-2 | 2.85e-3 | `[0.115,0.454]` | 4.75e-4 | `[0.019,0.454]` | 1.66e-6 | `[6.9e-5,1.3e-3]` | 1.16e-3 | `[0.047,0.125]` |
| equicorr `rho=0.85` | 1.619e-2 | 2.85e-3 | `[0.108,0.454]` | 4.76e-4 | `[0.018,0.454]` | 1.66e-6 | `[6.2e-5,1.3e-3]` | 1.16e-3 | `[0.044,0.125]` |
| equicorr `rho=0.95` | 1.795e-2 | 4.75e-3 | `[0.154,0.454]` | 6.89e-4 | `[0.022,0.454]` | 1.66e-6 | `[5.4e-5,1.3e-3]` | 1.64e-3 | `[0.053,0.125]` |
| **cancelling `r=0.9`, `b=3`** | **2.909e-5** | 1.17e-5 | `[0.171,0.461]` | 3.73e-6 | `[0.035,0.461]` | 9.30e-6 | `[7.8e-4,**0.320**]` | 6.36e-6 | `[0.085,0.219]` |
| **cancelling `r=0.95`, `b=3`** | **1.122e-5** | 3.68e-6 | `[0.142,0.462]` | 1.09e-6 | `[0.032,0.462]` | 9.30e-6 | `[5.5e-4,**0.829**]` | 1.99e-6 | `[0.076,0.178]` |
| cancelling `r=0.9`, `b=1` | 1.777e-2 | 8.22e-3 | `[0.458,0.463]` | 4.11e-3 | `[0.231,0.461]` | 4.96e-5 | `[2.8e-3,2.8e-3]` | 2.18e-3 | `[0.123,0.125]` |
| equicorr `rho=0`, `s*=5` | 5.956e-3 | 1.68e-3 | `[0.133,0.451]` | 1.60e-4 | `[0.009,0.451]` | 1.81e-7 | `[1.0e-5,3.5e-4]` | — | — |
| equicorr `rho=0.5`, `s*=5` | 1.795e-2 | 2.36e-3 | `[0.102,0.451]` | 1.69e-4 | `[0.006,0.451]` | 5.76e-8 | `[2.1e-6,3.5e-4]` | — | — |
| equicorr `rho=0.85`, `s*=5` | 1.669e-2 | 2.36e-3 | `[0.093,0.450]` | 1.95e-4 | `[0.008,0.450]` | 5.76e-8 | `[1.9e-6,3.5e-4]` | — | — |
| equicorr `rho=0.95`, `s*=5` | 1.797e-2 | 2.36e-3 | `[0.080,0.450]` | 1.80e-4 | `[0.006,0.450]` | 5.76e-8 | `[1.9e-6,3.5e-4]` | — | — |

The v2 bound `1/{32 p s_0(2s*+1) max(1,Theta_lambda)(1+k Theta_lambda)^{k}}`
attains, over the same grids, `8.99e-23`, `4.79e-29`, `7.07e-127`,
`1.84e-143`, `1.15e-51`, `4.67e-28`, `1.12e-4`, `1.98e-44`, `2.67e-51`,
`8.05e-238`, `0` — it is not on the same scale as anything else here.

`eq:gapM0` applies only where `k_lambda = 0`, where it reads `7.97e-5`
(`s*=3`), `1.12e-4` (`s*=2`) and `3.38e-5` (`s*=5`) against a truth of
`1.795e-2` — factors `225`, `161` and `531`.  Corollary B1 gives `8.15e-3`
(`s* = 3, 2`) and `8.09e-3` (`s* = 5`) at the same temperatures, a factor `2.2`.

### 4.2  The exactness of the factorisation (`gapB_checks.txt`, item A)

The proof's one inequality that could be badly lossy is (3).  Its exact
constant
`C_blk = sup_f Var(f)/{sum_j E[Var(f|gamma_{-j})] + E[Var(f|gamma_{V_0^c})]}`
was computed by dense generalised diagonalisation on all 1471 states, at four
temperatures and six designs:
```
    max over the 24 cases   C_blk = 1.0227      proved bound 1/(1-alpha) = 2.191
    exact weighted Dobrushin alpha in [0.049, 0.404]   proved  alpha <= 0.5436
```
**[computed]**  The factorisation is essentially lossless; the slack is in the
Dobrushin estimate, and the largest single contribution to it is the trivial
`R_{V_0,j} <= 1`.

### 4.3  The degrading quantity

Section 0's table.  Restated as a claim: over 231 (design, temperature) pairs,
`varrho_lambda <= 1.61` and `4 p kappa^+_lambda >= 0.21` on every
equicorrelated design at `s* = 3` and `s* = 5` and on the benign `b=1`
cancelling pair, while on the two `b=3` cancelling pairs `varrho_lambda` rises
to `184` and `628` and `4 p kappa^+_lambda` falls to `1.4e-3` and `4.0e-4`, at
the same temperatures at which `4p Gap` falls to `1.63e-3` and `6.28e-4`.
**[computed]**

### 4.4  The counterexample to `eq:gapM` (`gapB_witness.txt`)

Design: `p=14`, `n=100`, `s_0=4`, `s*=2`; `X_0, X_1` correlated at `r=0.98`,
the other twelve columns residualised against their span;
`beta* = (4,-4,0,...,0)`; `sigma_0 = 1`.  Verified exactly:

```
    ass:twosided:  max_{j in S^c, gamma} psi_j(gamma) = -2.1946 = (a_0-alpha) log p,
                   a_0 = 0.6684,  kappa+alpha = 3.5 >= 2+a_0 = 2.668,  kappa_1 = 2.832 >= 2
    good event:    A_n holds with L = 1.821  (ass:B asks L >= 4/nu, so A_n is not binding)
                   B_n: 1.1816 <= 8 ;  C_n: 0.0854 <= 1/2 ;  D_n: 0.0677 <= 8.278
    def:cpm:       m^-_1 = 0, m^+_1 = W_1 = {0,1}, k_1 = 2
    pibar_1 on 2^W = ( 4.013e-5, 4.638e-9, 4.640e-9, 1.000 )  -> TWO local maxima
                     for flips and exchanges, so ass:unimodal FAILS
    Gap(P_1)  = 3.4688e-06        eq:gapM = 9.3006e-06        ratio 2.6812
    Theorem B = 9.5038e-07 (exact form, ratio 0.274) / 1.5843e-07 (path form, 0.046)
```
**[computed]**  `th:gapM` is stated under `ass:twosided`, the good event and
`ass:unimodal`; the first two hold here and the third does not, so this is a
counterexample to `eq:gapM` with `ass:unimodal` deleted — not to `th:gapM` as
stated.  Over the 120-design search of `gapB_adversarial.py` (60 cancelling
pairs and 60 cancelling triples: five seeds, `r in {0.9,0.95,0.98}`,
`b in {2,3,4,6}`, 60 temperatures each) 24 designs violate it, all of them
cancelling pairs, worst factor `2.681`; Theorem B holds at every one of them,
with ratio `0.042`--`0.33` in its exact form.

### 4.5  Reproducing

```
    python3 gapB_numerics.py           ->  gapB_numerics.txt        (~80 s)
    python3 gapB_numerics.py --big     ->  gapB_numerics_big.txt    (~120 s)
    python3 gapB_checks.py             ->  gapB_checks.txt          (~40 s)
    python3 gapB_adversarial.py        ->  gapB_adversarial.txt     (~290 s)
    python3 gapB_witness.py            ->  gapB_witness.txt         (~10 s)
```

---

## 5.  The guidance, item by item

1. **"Start at the cold end."**  Survives, but *not* in the form "perturb the
   prior".  A Dobrushin/comparison continuation from `lambda = 0` needs
   `lambda max_{i in S} sum_j |delta_ij| <~ 1/p`, i.e.
   `lambda <~ nu/(n p)` by `NOTES_tensorisation` (8) — very much shorter than
   the interval that Corollary B5 already covers for free, and it gives nothing
   beyond it.  The correct cold-end statement is combinatorial, not
   perturbative: *below `lambda_cold` every influential covariate is decided
   out*, the target is a product measure up to Dobrushin-negligible
   interactions, and `Gap >= 3/(16p)`.  **[proved]**  But `lambda_cold` is
   `0.024`--`0.101` **[computed]** against a `lambda_star` that can be as large
   as `1/2`, so the cold segment is a *small* fraction of what has to be
   covered, and "continuation upward in temperature" is not how the rest gets
   done — it gets done by the same theorem at every `lambda`.

2. **"The centre must move with `lambda`."**  Survives, and is superseded:
   there is no centre.  See section 3.

3. **`alpha log p`.**  Adopted throughout.

4. **"Is the exponent 3 of `cor:critical` required?"**  No — and it buys
   nothing structural.  **[proved, `gapB_checks.txt` item C]**  Replacing `3`
   by `q` changes `lambda_star` to
   `max{(kappa+q)log p/(Psi - alpha log p), 2q log p/Psi}`; the first term is
   `<= (kappa+q)/(2(L+Ltil)+alpha+2kappa)`, which tends to `1/2` as
   `kappa -> infinity` *for every fixed `q`*, so no `q` lowers the a priori
   bound on `lambda_star` below `1/2`.  What the exponent `3` does buy is the
   concentration rate of `prop:conc` (`1-O(p^{-1})` rather than
   `1-O(s* p^{-1})`) and the freedom, in `cor:prec`, from any condition
   relating `s*` to `p`: with a common exponent `2` one needs
   `(s*+1) p^{-1} <= 1/2`.  Theorem B uses neither `cor:critical` nor
   `cor:prec`.

5. **"`cor:below1` proves more than it states."**  Confirmed.  **[proved]**
   Under `ass:B` (`L nu >= 4`, `nu <= 1`, so `L >= 4`), `ass:C`
   (`alpha >= 1/2`, `kappa >= 2`, `kappa + alpha >= 8(L+Ltil)+2`, hence
   `kappa >= 34 - alpha`), `ass:betamin` and `ass:nsize` (which give
   `Psi >= 2 Q log p`, `Q := L+Ltil+alpha+kappa`),
   ```
       lambda_star <= max{ (kappa+3)/(2(L+Ltil)+alpha+2kappa),  3/Q }
                   <= max{ (kappa+3)/(8+alpha+2kappa),  3/38 }
                   <  1/2 ,
   ```
   the first term being `1/2 - (alpha/2 + 5/2)/(8+alpha+2kappa) < 1/2` and
   increasing in `kappa`, so `1/2` is the exact supremum and is not attained
   (grid search reaches `0.499999`).  There is no dependence on `n`, `p`, `s*`,
   `nu` or `C_beta`.
   *What it means, taken seriously, and the caveat.*  Under Assumptions
   3.1--3.6 the region Section 6 has to cover is the fixed initial segment
   `[0, 1/2]`, and on it the likelihood is damped by at least a factor two.
   This does buy something: on `[0,1/2]` every log-barrier
   `lambda x (barrier in psi units)` is halved, so the worst-case degradation
   is square-rooted.  **But it is not where the difficulty is.**  In every
   design in which the exact gap varies appreciably — the cancelling pairs,
   where it falls by three orders of magnitude — its minimum over the grid is
   at `lambda = 1` **[computed: `gapB_checks.txt` item B, argmin `lambda` = 1
   in all nine cases]**; on the equicorrelated designs `4p Gap` stays in
   `[0.91, 1.75]` and has no meaningful minimum.  So the hard temperatures lie
   in `[lambda_star, 1]`, exactly where `th:mixing` already applies — *under
   Assumptions 3.1--3.6*.  Those assumptions are what
   exclude the cancelling pair: `ass:betamin` asks
   `n nu^2 C_beta^2 >= 128(L+Ltil+alpha+kappa) sigma_0^2 log p`, whose right
   side is at least `2534` at `p=14`, and the left side is `7.66` for the
   cancelling pair at `r=0.9` and `0.29` for the equicorrelated design at
   `rho=0.85` **[computed]** — no `p=14` design meets it, the benign ones
   included.  So the `p = 14` experiments can test the *architecture* of a
   bound but cannot test the assumption set, and the honest reading of
   `lambda_star < 1/2` is: it partitions the temperature axis correctly, but it
   does not localise the hard part of the problem, because the hard part is
   excluded by hypothesis rather than by temperature.

---

## 6.  Obstructions and open points

**6.1  `kappa^+_lambda` is a conductance of the target, and that is
unavoidable.**  The brief asks for a bound whose only extra hypothesis is
checkable from `(X, beta*, sigma_0, g, kappa)` and is not about the
conductance of the target.  Theorem B needs *no* extra hypothesis — but the
price is that the conductance appears inside `B(lambda)` as
`kappa^+_lambda`.  This is forced.  **[computed]** the cancelling pair at
`r=0.95` has `4p Gap = 6.28e-4`, so no bound of the form `c/(p . poly(s_0,s*))`
with an absolute `c` can be correct; **[proved]** the only quantity Theorem B
leaves free is `kappa^+_lambda`, so it is the one that has to fall.  It does,
and by the right amount **[computed]**: over `lambda in [0,1]`,
`4p kappa^+_lambda` falls by a factor `559` and `4p Gap` by `617` at `r=0.9`;
`2175` and `1600` at `r=0.95`.
The mitigating fact is that `kappa^+_lambda` is a gap on at most `2^{s*}`
states and is computed by enumeration; it is a *small, explicitly localised*
residual, not the original problem.

**6.2  The capped fibres: the one place the truncation is not free.**  If
`|sigma| > s_0 - |m^+_lambda|` the fibre `F_sigma` is
`{A ⊆ m^+_lambda : |A| <= s_0 - |sigma|}`, additions inside it can be illegal,
and both Corollary B4 (which routes upward) and Corollary B6 (which adds the
decided covariates) break.  Two facts, both computed:
* the capped fibres *are* the binding ones: the argmin of `kappa_sigma` is a
  capped fibre at 16, 19, 19, 19 of 21 temperatures on the four `s*=3`
  equicorrelated designs and at 13--20 of 21 at `s*=5`;
* but they cost only a factor `1.00`--`3.22`: `4p kappa^+_lambda` computed over
  all `sigma` against `4p kappa^+_lambda` computed over uncapped `sigma` only
  is `0.458/1.068`, `0.350/1.118`, `0.350/0.653`, `0.583/1.052` at
  `lambda = 1`, and `1.00` exactly on the cancelling pairs (there
  `|m^+_lambda| = 2 <= s_0 - |sigma|` for the binding `sigma`).
I have **no a priori bound** for a capped fibre.  The obstruction is real: in a
capped fibre two configurations of the same size are joined only by exchanges,
and an exchange that removes an undecided covariate has no a priori cost bound,
whereas the funnel of `lemma:funnel` avoids the problem altogether by deleting
every non-`gamma*` covariate first — at the price of the factor `s_0 s*` that
Theorem B removes.  **This is the price of the block architecture, and it is
the one structural respect in which `th:gapM`'s architecture is better.**
A fix would need a variance factorisation that is allowed to ignore fibres of
small mass; the mass of the capped fibres is at most `delta-bar^{s_0-s*+1}`
**[conjecture: the stochastic-domination step of `NOTES_tensorisation` section
5 item 2 is stated there but I did not verify it]**, but approximate
tensorisation as used here takes a *minimum* over fibres and is blind to their
mass.

**6.3  Theorem B does not subsume `th:mixing`.**  Under Assumptions 3.1--3.6
and `lambda >= lambda_star` one has `m^+_lambda = gamma*` **[proved: by
`th:drop`, every `j in gamma*` occurs as the greedy addition `j_{gamma_r}` at
some point of the chain of `def:chain`, so `psi_j^+ >= Psi_r - alpha log p >=
(2Q - alpha) log p`, and `lambda(2Q-alpha) log p >= (kappa+3) log p > kappa log
p - log 8s*`]**, but `m^-_lambda` need not be `gamma*` (`lemma:forward` bounds
the *greedy* gain, not every gain — `NOTES_tensorisation` section 6), so the
window need not shrink and `kappa^+_lambda` is not bounded by Assumptions
3.1--3.6.  On *uncapped* fibres one does get `kappa_sigma >= 1/(32 p s*)` from
`cor:critical` (the greedy addition stays inside the fibre and drops by
`p^{-3}`; with at most `s*` immediate precedents per state the precedent mass
ratio is `(1 - s* p^{-3})^{-1} <= 2`, so an increasing flip edge has congestion
at most `16 p` and the routes have length at most `2 s*`), hence
`Gap >= 1/(96 p s*)` — better than `1/(64 p s_0^2)` by a factor
`(2/3) s_0^2 / s*` — but 6.2 blocks it.  Practically: report
`max{ Theorem B, eq:tau's 1/(64 p s_0^2) for lambda >= lambda_star }`.

**6.4  The a priori congestion bound is the weak link.**  Corollary B4 is
`p^{-(k_lambda+1)(kappa + lambda alpha)}` up to constants, which fails the
brief's non-triviality test (`B(lambda) >= p^{-s*}`) whenever
`(k_lambda+1)(kappa + lambda alpha) > s*` — always, since `kappa + alpha >= 5/2`
by `ass:C`.  I do not think this can be repaired without a hypothesis: the base
`p^{kappa + lambda alpha}` per undecided covariate is the exact cost of the
first addition into an empty context, `exp(kappa log p - lambda psi_j(sigma))`,
and `psi_j(sigma)` can be as low as `-alpha log p` **[proved: `psi_j >= -alpha
log p` with equality iff the column is orthogonal to `Y`]**; the cancelling
pair realises `psi_j(0) = -0.63` against `-alpha log p = -3.96` **[computed]**,
so the base is approached but not attained in the designs tested.

**6.5  Things I did not try.**  (a) A fractional (multicommodity) second level
inside a fibre — `NOTES_flow` section 3 shows the LP optimum is `1.00`--`2.27`
where the best single route is `2.9e5`, so it would tighten the path form of
`B(lambda)` but, as there, yields no theorem.  (b) Optimising the block weight
`c_0` per design; the exact `alpha` is `0.049`--`0.404` against the proved
`0.5436`, so at most a factor `1.6` is available.  (c) Sharpening
`R_{V_0,j} <= 1`; this single entry contributes `c_0 = 0.352` of the `0.544`,
and any bound on the influence of the window block on a null coordinate would
improve `1 - alpha_lambda` towards `1 - delta-bar = 0.81`, i.e. a further
factor `1.8`.  (d) The case `p < 8`, where `delta-bar > 1/4` and the constants
must be redone.
