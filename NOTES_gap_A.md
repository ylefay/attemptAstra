# A lambda-dependent spectral gap bound: the best-subset residual profile

Working notes, attempt A.  Every number below is produced by one of

    gapA_probe.py       -> gapA_probe.txt        (exploration; four route families)
    gapA_validate.py    -> gapA_validate.txt     (the theorem, validated per temperature)
    gapA_lambdastar.py  -> gapA_lambdastar.txt   (brief items 4 and 5)
    gapA_adversarial.py                          (how far the cancelling pair goes)
    gapA_table.py                                (compact tables, read back from the .txt)

All of it is exact: `M(s_0)` is enumerated, `pi^(lambda)` evaluated at every
state, the flip-and-exchange kernel assembled, the gap obtained by
diagonalisation, and every congestion *computed* rather than bounded.  The
design/likelihood machinery is `window_numerics.py` and `flow_numerics.py`,
reused unchanged.  Nothing in `binary_variable_selection3.tex`, `bvs3_app.tex`
or any existing `.py` has been touched.

Throughout, `proved`, `verified numerically` and `conjectured` are kept apart
and labelled.

---

## 0. Answer, up front

**There is a lambda-dependent bound of the required shape, and the quantity
that carries the lambda-dependence is the *best-subset residual profile* of the
influential block — not the modal structure of the target.  Two things go away:
`ass:unimodal` is replaced by a lambda-free hypothesis about
`(X, beta*, sigma_0, g)` alone, and one of the two exponentials of `eq:gapM`
(the one in the path *length*) is removed, so `4^{s*}` becomes `2^{s*}`.  One
thing does not go away: a factor bounded only by `2^{s*}`, which is a
*tie-counting* cost of the path method and not a bottleneck cost, remains.  I
could not remove it, and section 9 says exactly where it sits.**

The construction is in two levels and has no `W_lambda`, no `m^{-}_lambda`,
no `Theta_lambda` and no margin `log(2 s*)`.

* **Level 1** funnels the null coordinates away: `F(gamma)` deletes
  `max(gamma \ gamma*)`.  This uses `ass:twosided` and nothing else, has orbits
  of at most `s_0` single flips, and costs a factor `e^{1/p} <= 1.07` in mass
  (measured: `1.047` to `1.063`), against the `2(2s*+1) = 14` of
  `lemma:funnel`\ref{it:fn_iii}.
* **Level 2** is a chain on `2^{gamma*}`, at most `2^{s*}` states, on which the
  target is *exactly* the two-parameter exponential family

        pi^(lambda)(A)  proportional to  D_A^{-lambda n/2} p^{-(kappa + lambda alpha)|A|},
        D_A := g^{-1}||Y||^2 + ||(I - Phi_A)Y||^2 ,

  so the whole second level is a function of the `2^{s*}` residuals `D_A` and of
  the single tilt direction `(-lambda n/2, -(kappa + lambda alpha) log p)`,
  which rotates with the temperature.  The route is: *align by exchanges to the
  greedy chain, then walk the greedy chain to the level `m_lambda` that the
  temperature selects.*  **The centre `P_{m_lambda}` moves with `lambda`, from
  the empty model at `lambda = 0` to `gamma*` at `lambda = 1`** — this is item 2
  of the brief, realised.

The bound is

        Gap(P_lambda)  >=  1 / { 8 p (2 s_0 + 4 s*) max(1, 2 s* T_lambda) },

`T_lambda = T_lambda(Gamma_lambda)` a *tie index* and `Gamma_lambda` the *dip*
of the route, both explicit and both computed exactly below.  `Gamma_lambda` is
the quantity that must degrade, and it does: it equals `1` at every temperature
on all four equicorrelated designs at `s* = 3`, and rises to `367`, `1256`,
`3645` on the cancelling pairs, at exactly the temperatures at which the true
gap falls, tracking the true deficit `(1/4p)/Gap(lambda)` to within a factor
`0.6` to `1.7` over three orders of magnitude.  `T_lambda` does **not** degrade
there (it stays near `2`), which is the numerical proof that the two factors
measure different things.

`Gamma_lambda = 1` for every `lambda in [0,1]` is guaranteed by one lambda-free
design hypothesis, **(B)**: *the matched exchanges toward the greedy chain do
not increase the residual, and the greedy marginal gains are non-increasing*.
Clause (B2) — diminishing returns of the best-subset regression path — is
exactly what the cancelling pair violates: there `g_0 = -0.63` and
`g_1 = +37.8`, the second variable being worth 60 times the first.

Against the truth: the bound is a **fixed factor** `755`–`1013` below the exact
gap at `s* = 3`, `258`–`266` on the cancelling pairs, and `6350`–`7920` at
`s* = 5`, *uniformly in lambda, and in particular on the cancelling pairs,
where the truth itself falls to `3.0e-04` of `1/(4p)`*.
For comparison, on the same runs `eq:gapM` bottoms out at `10^{-127}`
(`s* = 3`) and `10^{-267}` (`s* = 5`), and the two-level bound (T-U) of
`NOTES_flow`, whose hypothesis fails on the cancelling family, is **false**
there: at `(r,b) = (0.97,5)`, `(0.97,8)`, `(0.99,8)` it exceeds the true gap by
`1.19`, `1.36`, `1.71`.

Two side answers to the brief.  **`lambda_star < 1/2` always** (proved, and
`1/2` is the exact supremum over the admissible parameter region), and this is
unchanged if the exponent `3` of `cor:critical`\ref{it:crit_under} is lowered
to `2` — item 4 buys a constant and nothing structural, item 5 is true but
does not help, because the cold end reaches only `lambda ~ 0.03`–`0.06` while
the window is non-empty on 92%–96% of `[0,1]` on the four equicorrelated
designs.

---

## 1. The reduction, and the identity that makes it work

### 1.1 The target in the residual variable

**Proved.**  Write `D_gamma := g^{-1}||Y||^2 + ||(I - Phi_gamma)Y||^2`, the
regularised residual, so that `||Y||^2 (g^{-1} + 1 - R^2_gamma) = D_gamma`.
Then `1 + g(1 - R^2_gamma) = g D_gamma / ||Y||^2`, and `eq:Lik` together with
`eq:prior` and `eq:alphadef` gives, for every `gamma in M(s_0)` and every
`lambda >= 0`,

        pi_n^{(lambda)}(gamma | Y)  =  w_lambda(gamma) / sum_{gamma'} w_lambda(gamma'),
        w_lambda(gamma) := exp{ -(lambda n/2) log D_gamma - (kappa + lambda alpha)|gamma| log p }.   (1)

This is an identity, not an estimate: the whole tempering path is a
two-parameter exponential family with sufficient statistic
`(log D_gamma, |gamma|)` and natural parameter
`(-lambda n/2, -(kappa + lambda alpha) log p)`.  Raising the temperature
rotates the tilt from "prefer few variables" towards "prefer small residual",
and nothing else about the target changes.  In this notation

        psi_j(gamma) = -alpha log p + (n/2) log( D_gamma / D_{gamma u j} )              (2)

is `eq:psidef`, and `eq:onestepG` is the increment of `log w_lambda`.

Everything below is written in `w_lambda`; `alpha log p` replaces
`(1/2)log(1+g)` everywhere, as item 3 of the brief asks.

### 1.2 The two levels

Put `S := gamma*` and, for `A subset S`,

        B_A := { gamma in M(s_0) : gamma n gamma* = A },        pibar_lambda(A) := pi^(lambda)(B_A).

The blocks partition `M(s_0)` into at most `2^{s*}` classes.  Note the
difference from `def:Mlambda`: the blocks are indexed by *all* of `gamma*`, not
by an undecided window, and their representative `v_A = A` is a pure state.
There is no `m^{-}_lambda`, no `m^{+}_lambda`, no `Theta_lambda`, and no
temperature enters the definition of the partition.

---

## 2. The construction

### 2.1 Level 1: the null funnel

**Definition.**  `F(gamma) := gamma \ {max(gamma \ gamma*)}` if
`gamma notsubset gamma*`, and `F(gamma) := gamma` otherwise.

`F` deletes only coordinates of `S^c`, so it preserves `gamma n gamma*`, maps
`B_A` into itself, and has `A` as its unique fixed point there; its orbits are
single flips and have at most `|gamma \ gamma*| <= s_0` edges.

### 2.2 Level 2: the greedy chain and the E-route

**Definition (greedy chain and its gains).**  Let `j_1, ..., j_{s*}` be defined
by `j_{m+1} := argmin_{j in S \ P_m} D_{P_m u {j}}`, `P_m := {j_1,...,j_m}`
(ties lexicographic), and put

        g_m := (n/2) log( D_{P_m} / D_{P_{m+1}} ) - alpha log p  =  psi_{j_{m+1}}(P_m),
        phi_lambda(m) := log w_lambda(P_m),      phi_lambda(m+1) - phi_lambda(m) = lambda g_m - kappa log p,
        m_lambda := min argmax_{0 <= m <= s*} phi_lambda(m).                          (3)

`g_m` is the marginal log-gain of the greedy (best-subset) regression path; it
is a function of `(X, Y, g)` and of nothing else.  `m_lambda` is the size the
temperature selects, and it is the only place where `lambda` enters the
construction.  `m_0 = 0` and, when `max_m lambda g_m > kappa log p` at
`lambda = 1`, `m_1 = s*`.

**Definition (the E-route).**  For `A subset S` the *orbit* of `A` is

  * **(E1) align.**  Write `A \ P_{|A|} = {a_1 < ... < a_q}` and
    `P_{|A|} \ A = {b_1 < ... < b_q}` (the two have the same size).  Perform
    the exchanges `a_i -> b_i`, `i = 1..q`, in this order.  The orbit is now at
    `P_{|A|}`.  These are `N_2` moves between states of size `|A| <= s*`.
  * **(E2) climb or descend the chain.**  From `P_{|A|}` add `j_{l+1}` while
    `l < m_lambda`, or delete `j_l` while `l > m_lambda`.  These are `N_1`
    moves.

The route `R_{A,A'}` is the orbit of `A` followed by the reversed orbit of
`A'`, truncated at their first common state.  The whole canonical path from
`x in B_A` to `y in B_{A'}` is: the `F`-orbit of `x`, then `R_{A,A'}`, then the
reversed `F`-orbit of `y`.

The orbit has at most `q <= s*` exchanges and at most `s*` flips, so

        L_lambda := max route length  <=  4 s* .                                       (4)

### 2.3 The two quantities the bound carries

        Gamma_lambda := max over ordered pairs (A,A') and states z on R_{A,A'} of
                        min{ pibar_lambda(A), pibar_lambda(A') } / pibar_lambda(z)   >= 1,     (5)

the **dip** of the ensemble, and

        T_lambda(c)  := max_{A subset S}  pibar_lambda{ B : pibar_lambda(B) <= c pibar_lambda(A) }
                                          / pibar_lambda(A) ,
        T_lambda     := T_lambda(Gamma_lambda) ,                                       (6)

the **tie index**.  Both are computed by one pass over `M(s_0)` and `4^{s*}`
route evaluations; both are `>= 1`; and `T_lambda(c) <= c 2^{s*}` always.
Replacing `pibar_lambda` by `w_lambda` in (5)-(6) — which makes them functions
of the `2^{s*}` numbers `D_A` alone — costs a factor `e^{1/p}` each way, by
Lemma 1 below.

---

## 3. The theorem

> **Theorem A (proved).**  Let `p >= 2`, `s* <= s_0`, `lambda in [0,1]`.
> Assume `ass:twosided` (Assumption 6.1), so that `kappa_lambda >= 2` by
> `eq:kappalam`, and work on the good event.  Let `varrho_lambda` and
> `L_lambda` be the congestion and the longest route of *any* family
> `{R_{A,A'}}` of paths on `(2^{gamma*}, flips and exchanges)`, the congestion
> normalised as in `def:flowens`.  Then
>
>     Gap(P_lambda)  >=  1 / { 4 p e^{1/p} max(1, s* varrho_lambda) (2 s_0 + L_lambda) }.   (T1)
>
> For the E-route of section 2.2,
>
>     varrho_lambda <= 2 T_lambda(Gamma_lambda),      L_lambda <= 4 s*,
>
> so, using `e^{1/p} <= 2`,
>
>     Gap(P_lambda)  >=  1 / { 8 p (2 s_0 + 4 s*) max(1, 2 s* T_lambda(Gamma_lambda)) }.     (T2)
>
> Unconditionally `T_lambda(Gamma) <= Gamma 2^{s*}`, hence
>
>     Gap(P_lambda)  >=  1 / { 16 p s* 2^{s*} (2 s_0 + 4 s*) Gamma_lambda }.                 (T3)

> **Proposition B (proved).**  Suppose, on the good event,
>
>   * **(B1) exchange descent.**  For every `A subset gamma*`, the matched
>     exchanges of step (E1) do not increase `D`:
>     `D_{A_i} <= D_{A_{i-1}}` at every step, `A_0 = A`, `A_q = P_{|A|}`.
>   * **(B2) diminishing greedy gains.**  `g_0 >= g_1 >= ... >= g_{s*-1}`.
>
> Then `Gamma_lambda <= e^{1/p} <= 2` for **every** `lambda in [0,1]`, and (T2)
> reads `Gap(P_lambda) >= 1/{8 p (2 s_0 + 4 s*) max(1, 2 s* T_lambda(2))}`,
> with `T_lambda(2) <= 2^{s*+1}`.

Note (B1) implies that `P_m` is a global minimiser of `D` over the `m`-subsets
of `gamma*`, and (B2) is equivalent to convexity of `m -> log D_{P_m}`.  Both
clauses are statements about the `2^{s*}` numbers `{D_A : A subset gamma*}`,
that is about `(X, beta*, sigma_0, g)` on the good event; **neither mentions
`lambda`, `kappa`, the target, its modes or its conductance.**

> **Proposition C (the cold end; proved, no hypothesis beyond the good event).**
> Let `Psi^{max} := max{ psi_j(gamma) : gamma in M(s_0), |gamma| < s_0, j notin gamma }`
> and
>
>     lambda_cold := ( kappa - 1 - log 2 / log p ) log p / Psi^{max} .
>
> For every `lambda in [0, lambda_cold]` the map `gamma -> gamma \ {max gamma}`
> is a valid transition function with centre the null model and every step
> `pi^(lambda)`-increasing, and
>
>     Gap(P_lambda)  >=  1 / ( 4 p e^{1/2} . 2 s_0 )  >=  1 / (14 p s_0).
>
> On the good event, `Psi^{max} <= -alpha log p + 8 ||X_S beta*_S||_2^2/sigma_0^2
> + 16(L + Ltilde) log p`, by the computation of `prop:a0` with the restriction
> `j in S^c` dropped.

> **Proposition D (`lambda_star`; proved, and verified numerically).**  Under
> `ass:A`--`ass:nsize`, and with the exponent `3` of `eq:underconc` replaced by
> an arbitrary `q > 0` (so that `eq:lambdastar` becomes
> `max{(kappa+q)log p / (Psi - alpha log p), 2q log p / Psi}`),
>
>     lambda_star  <  1/2 ,
>
> and `1/2` is the exact supremum over the region cut out by
> `ass:B`, `ass:C`, `ass:betamin`, `ass:nsize`, for **every** `q`.  In
> particular `lambda_star` is bounded by an absolute numerical constant with no
> dependence on `n, p, s*, nu, C_beta`, and lowering the exponent buys a
> constant only.

---

## 4. Proofs

Throughout, `mu := pi_n^{(lambda)}(. | Y)`, and the assumption used is flagged
at the point of use.

### 4.1 Lemma 1 (the funnel)

> **Lemma 1.**  On the good event and under `ass:twosided`: (i) every step of
> `F` is `mu`-increasing, with `mu(gamma)/mu(F(gamma)) <= p^{-kappa_lambda}
> <= p^{-2}`; (ii) `mu[Lambda^F(gamma)] <= e^{1/p} mu(gamma)` for every
> `gamma`; (iii) in particular `pi^(lambda)(A) <= pibar_lambda(A)
> <= e^{1/p} pi^(lambda)(A)` for every `A subset gamma*`.

*Proof.*  (i)  Let `j = max(gamma \ gamma*) in S^c` and `gamma' = gamma \ {j}`,
so `|gamma'| < s_0` and `j in S^c \ gamma'`.  By `eq:onestepG`,
`mu(gamma)/mu(gamma') = exp{lambda psi_j(gamma') - kappa log p}`.
**[`ass:twosided` used here, and only here]**  `psi_j(gamma') <=
(a_0 - alpha) log p`, so the ratio is at most
`p^{lambda(a_0 - alpha) - kappa} = p^{-kappa_lambda} <= p^{-2}` by
`eq:kappalam`.

(ii)  `bar gamma in Lambda^F(gamma)` at distance `k` means
`bar gamma = gamma u sigma` with `sigma subset S^c`, `|sigma| = k` and
`min sigma > max(gamma \ gamma*)`; there are at most `C(p,k) <= p^k/k!` of
them, and telescoping (i) along the orbit gives
`mu(bar gamma) <= p^{-2k} mu(gamma)`.  Hence
`mu[Lambda^F(gamma)]/mu(gamma) <= sum_k (p^k/k!) p^{-2k} = e^{1/p}`.

(iii)  `B_A = Lambda^F(A)` by the first paragraph of section 2.1, and
`A in B_A`. `[]`

Measured (`gapA_validate.txt`): `max_A pibar_lambda(A)/pi^(lambda)(A)` is
`1.057` (`s* = 3`), `1.063` (`s* = 2`), `1.047` (`s* = 5`), against
`e^{1/14} = 1.0742`.

### 4.2 Proof of (T1)

The ensemble uses two disjoint kinds of edge: *funnel* edges
`(bar gamma, F(bar gamma))`, which join two states of one block, and *block*
edges `(A, A')` with `A != A'`, which join two pure states of different blocks.
So `rho(T)` is the larger of the two maxima, not their sum.

*Funnel edges.*  Such an edge is a single flip and is increasing by Lemma 1(i),
so by `eq:MH` (proposal `1/(2p)`, laziness `1/2`)
`Q_lambda(e) = mu(bar gamma)/(4p)` exactly — the factor `s_0` of `eq:edge` is
not needed, `F` using only flips.  The forward orientation can be used only in
the first leg of a path, whose source must then lie in `Lambda^F(bar gamma)`,
and the reverse orientation only in the third leg, whose sink must; either way
the demand is at most `mu[Lambda^F(bar gamma)] <= e^{1/p} mu(bar gamma)` by
Lemma 1(ii).  Congestion `<= 4 p e^{1/p}`.

*Block edges.*  `e = (A, A')` is a flip or an exchange between two states of
size at most `s*`.  For an exchange the proposal is
`{2 |A| (p - |A|)}^{-1} >= {2 s* p}^{-1}`, so
`Q_lambda(e) >= min{mu(A), mu(A')}/(4 p s*)`, and
`mu(A) >= e^{-1/p} pibar_lambda(A)` by Lemma 1(iii).  The demand across `e` is
`sum over pairs (x,y) whose route uses e of mu(x)mu(y)
 = sum over pairs (A_1,A_2) whose second-level route uses e of
   pibar_lambda(A_1) pibar_lambda(A_2)
 <= varrho_lambda min{pibar_lambda(A), pibar_lambda(A')}` by the definition of
`varrho_lambda`.  Congestion `<= 4 p s* e^{1/p} varrho_lambda`.

Hence `rho(T) <= 4 p e^{1/p} max{1, s* varrho_lambda}`.  A path is one
`F`-orbit (`<= s_0` edges), one second-level route (`<= L_lambda`), and one
reversed `F`-orbit (`<= s_0`), so `ell(T) <= 2 s_0 + L_lambda`.  The paths are
simple: the two `F`-orbits lie in the distinct blocks `B_A`, `B_{A'}`, the
second-level route visits pure states of pairwise distinct blocks, and when
`A = A'` the two orbits are truncated at their first common state.
`eq:sinclair` with `lemma:kernel` gives (T1). `[]`

### 4.3 Proof of `varrho_lambda <= 2 T_lambda(Gamma_lambda)`

Fix a directed block edge `e = (z, z')` and put `m := min{pibar(z), pibar(z')}`.
Both `z` and `z'` lie on every route that uses `e`, so by the definition (5) of
`Gamma_lambda`, every pair `(A,A')` using `e` has
`min{pibar(A), pibar(A')} <= Gamma_lambda m`.  Splitting according to which of
the two attains the minimum, and using `sum_{A'} pibar(A') <= 1`,

    sum over pairs using e of pibar(A) pibar(A')
        <=  2 sum_{A : pibar(A) <= Gamma_lambda m} pibar(A)
        =   2 pibar{ B : pibar(B) <= Gamma_lambda m }.

Since `m` is the mass of an actual state, dividing by `m` and maximising gives
`varrho_lambda <= 2 T_lambda(Gamma_lambda)`, which is (6).  The trivial bound
`T_lambda(c) <= c 2^{s*}` holds because each of the at most `2^{s*}` states in
the set contributes at most `c pibar(A)`. `[]`

The length bound `L_lambda <= 4 s*` is (4).  This proves (T2) and (T3).

### 4.4 Proof of Proposition B

Under (B1), each exchange of step (E1) does not increase `D` and does not
change `|A|`, so by (1) it does not decrease `w_lambda`, for every
`lambda >= 0`.  Under (B2) the increments
`phi_lambda(m+1) - phi_lambda(m) = lambda g_m - kappa log p` are non-increasing
in `m`, so `phi_lambda` is concave, hence unimodal with maximum at `m_lambda`;
therefore each step of (E2), which moves `m` towards `m_lambda`, does not
decrease `phi_lambda`, that is does not decrease `w_lambda`.  So `w_lambda` is
non-decreasing along every orbit.  A route is an initial segment of the orbit
of `A` followed by a reversed initial segment of the orbit of `A'`, so every
state `z` on it has `w_lambda(z) >= min{w_lambda(A), w_lambda(A')}`.  By
Lemma 1(iii), `pibar_lambda` and `w_lambda/sum w_lambda` agree to within a
factor `e^{1/p}`, so
`min{pibar(A),pibar(A')}/pibar(z) <= e^{1/p}`, which is
`Gamma_lambda <= e^{1/p}`. `[]`

*Remark.*  (B2) alone controls (E2); (B1) alone controls (E1).  Numerically
(B1) is the mild one: it holds on 11 of the 12 designs tested, failing only at
`rho = 0.85, s* = 5`.  (B2) is the discriminating one.

### 4.5 Proof of Proposition C

At `lambda <= lambda_cold`, for every `gamma` with `|gamma| < s_0` and every
`j notin gamma`, `lambda psi_j(gamma) - kappa log p <= -(1 + log2/log p) log p`,
so `mu(gamma u {j}) <= (p^{-1}/2) mu(gamma)`: every addition is unfavourable by
a definite factor.  The map `G_0(gamma) = gamma \ {max gamma}` is then
increasing, its orbits are `<= s_0` flips, its unique fixed point is the null
model, and its precedents at distance `k` number at most `C(p,k) <= p^k/k!`, so
`mu[Lambda^{G_0}(gamma)] <= sum_k (p^k/k!)(p^{-1}/2)^k mu(gamma)
= e^{1/2} mu(gamma)`.  As in section 4.2, `rho <= 4 p e^{1/2}` and
`ell <= 2 s_0`.  The bound on `Psi^{max}` is the display in the proof of
`prop:a0` with `j in S^c` relaxed to `j notin gamma`: only `lemma:residual`,
`lemma:noise`\ref{it:noise_ii} and `ass:A` are used there, none of which needs
`j` to be uninfluential, and `||Pi X_S beta*_S||^2 <= ||X_S beta*_S||^2`
replaces the appeal to `lemma:proxy`. `[]`

### 4.6 Proof of Proposition D

By `eq:floorlower`, `Psi >= 2 Q log p` with `Q = L + Ltilde + alpha + kappa`,
so with a general exponent `q`

    lambda_star  <=  max{ (kappa + q)/(2Q - alpha),  q/Q } .

`ass:B` and `ass:C` give `L >= 4`, `Ltilde >= 0`, `alpha >= 1/2`,
`kappa >= 2`, so `Q >= 6.5` and `q/Q <= q/6.5`.  For the first term, write
`B := L + Ltilde >= 4`; then `2Q - alpha = 2B + alpha + 2 kappa` and

    (kappa+q)/(2B + alpha + 2kappa)  =  1/2  -  (2B + alpha - 2q) / (2(2B+alpha+2kappa))  <  1/2

whenever `2B + alpha > 2q`, which holds for `q <= 4` since `2B + alpha >= 8.5`.
The expression increases in `kappa` and tends to `1/2`, so `1/2` is the
supremum and is not attained.  Hence `lambda_star < 1/2`.  `[]`

**Verified numerically** (`gapA_lambdastar.txt`): exhaustive search over the
admissible region gives `sup = 0.499933, 0.499963, 0.499993` for `q = 2, 3, 4`
respectively, attained in the limit `kappa -> infinity` with `alpha = 1/2`,
`L + Ltilde = 4`.

---

## 5. Where the lambda-dependence enters, and what it costs

Exactly three places, and they are all in the second level.

1. **The centre.**  `m_lambda = argmax_m phi_lambda(m)` of (3).  It is `0` at
   `lambda = 0` (the prior mode) and `s*` at `lambda = 1` under the hypotheses
   of `cor:below1`, and it moves monotonically in between when (B2) holds,
   because `phi_lambda(m+1) - phi_lambda(m) = lambda g_m - kappa log p` is
   increasing in `lambda`.  This is the "staircase" of `def:chain`, but read
   off the *best-subset* residual profile rather than off the worst-case floor
   `Psi_r`, and it is the whole of item 2 of the brief: no ensemble directed at
   `gamma*` can work below `lambda_star` (`prop:sharp`), and this one is
   directed at `P_{m_lambda}`.
2. **The dip `Gamma_lambda`.**  It is `1` whenever `phi_lambda` is quasiconcave
   and (E1) descends; it is the exponential of the deepest valley the route has
   to cross.  It costs nothing when (B) holds and it costs exactly the true
   bottleneck when (B) fails — see section 6.3.
3. **The tie index `T_lambda`.**  It is the effective number of states of
   comparable mass; it is largest at `lambda -> 0`, where all `C(s*,m)` subsets
   of a given size tie under the prior, and it decreases as the temperature
   concentrates the target.  Measured maxima: `2.96`–`3.51` at `s* = 3` and
   `9.87`–`16.4` at `s* = 5`, against `max_m C(s*,m) = 3` and `10`.

The price paid, relative to `eq:gapM`, is: nothing in the hypotheses (one
lambda-free design hypothesis replaces `ass:unimodal`), a gain of one
exponential in the length (`4 s_0 + 2^{k+1}` becomes `2 s_0 + 4 s*`), and no
gain in the surviving exponential, since `T_lambda <= 2^{s*}` and `s*` replaces
`k_lambda <= s*`.  Where `k_lambda < s*`, `eq:gapM` has the smaller exponent
and the two bounds should be taken together (their maximum is a valid bound).

---

## 6. The numbers

`p = 14`, `n = 100`, `kappa = 2`, `alpha = 1.5`, 13 temperatures on
`lambda in [0.001, 1]`; `s* = 3` uses `s_0 = 4`, `beta*_S = (1.20, 0.70, 0.45)`
and the equicorrelated designs of `window_numerics.build_design`; `s* = 5` uses
`s_0 = 6`, `beta*_S = (1.20, 1.00, 0.85, 0.70, 0.45)`; the cancelling pairs are
`X_1 = r X_0 + sqrt(1-r^2) Z`, `beta* = (b, -b, 0, ...)`, null columns
residualised, as in `flow_numerics.adversarial`.  `1/(4p) = 1.7857e-02`.

### 6.1 The bound against the truth, and against the alternatives

`B(lambda)` is (T2).  `1/(rho ell)` is the exact congestion and length of the
ensemble that (T2) analyses, so the gap between the two columns is the cost of
the estimates and not of the construction.  `min` is over the temperature grid.

| design | min Gap | min B | worst B/Gap | max Gamma | max T | min 1/(rho l) | min eq:gapM | min (T-U) | min cert |
|---|---|---|---|---|---|---|---|---|---|
| equicorr rho=0.00, s*=3 | 1.795e-02 | 2.512e-05 | 1.33e-03 | 1 | 2.96 | 5.191e-04 | 8.99e-23 | 4.43e-06 | 1.402e-03 |
| equicorr rho=0.50, s*=3 | 1.669e-02 | 2.118e-05 | 9.87e-04 | 1 | 3.51 | 5.162e-04 | 4.79e-29 | 1.66e-06 | 1.157e-03 |
| equicorr rho=0.85, s*=3 | 1.619e-02 | 2.480e-05 | 1.02e-03 | 1 | 3.00 | 5.125e-04 | 7.07e-127 | 1.66e-06 | 1.158e-03 |
| equicorr rho=0.95, s*=3 | 1.795e-02 | 2.467e-05 | 1.16e-03 | 1 | 3.02 | 5.112e-04 | 1.84e-143 | 1.66e-06 | 1.644e-03 |
| cancelling r=0.90 b=1 | 1.777e-02 | 6.965e-05 | 3.88e-03 | 1 | 2.00 | 8.474e-04 | 1.12e-04 | 4.96e-05 | 2.177e-03 |
| **cancelling r=0.90 b=3** | **2.909e-05** | **3.792e-07** | **3.88e-03** | **367** | **368** | 4.846e-06 | 1.15e-51 | 9.30e-06 | 6.359e-06 |
| **cancelling r=0.95 b=3** | **1.122e-05** | **1.110e-07** | **3.76e-03** | **1256** | **1257** | 1.418e-06 | 4.67e-28 | 9.30e-06 | 1.994e-06 |
| cancelling r=0.99 b=3 | 1.785e-02 | 6.960e-05 | 3.88e-03 | 1 | 2.01 | 8.468e-04 | 1.12e-04 | 4.96e-05 | 2.206e-03 |
| equicorr rho=0.00, s*=5 | 3.696e-03 | 1.699e-06 | 1.57e-04 | 3.58 | 16.4 | 1.015e-04 | 3.20e-59 | 1.81e-07 | -- |
| equicorr rho=0.50, s*=5 | 1.483e-02 | 2.826e-06 | 1.57e-04 | 1 | 9.87 | 1.017e-04 | 3.37e-76 | 5.77e-08 | -- |
| equicorr rho=0.85, s*=5 | 1.405e-02 | 1.788e-06 | 1.26e-04 | 4.26 | 15.6 | 1.012e-04 | 5.83e-215 | 5.77e-08 | -- |
| equicorr rho=0.95, s*=5 | 1.795e-02 | 2.812e-06 | 1.45e-04 | 1 | 9.92 | 1.011e-04 | 3.73e-267 | 5.77e-08 | -- |

Read four things off it.

* **`B(lambda)/Gap(lambda)` is a fixed factor.**  It lies in
  `[9.9e-4, 3.9e-3]` at `s* = 3` and in `[1.26e-4, 8.0e-4]` at `s* = 5`, over
  every design and every temperature, *including the two cancelling designs on
  which the true gap falls by factors 614 and 1592*.  The bound therefore
  tracks the truth, which is what the brief asks for.
* **The construction is much better than the estimates.**  `1/(rho ell)` for
  the actual ensemble is a factor `34.6`–`38.9` below the truth at `s* = 3`,
  `21.2` on the cancelling pairs and `176`–`178` at `s* = 5`; (T2) is a further
  factor `20.7`–`24.4` (`s* = 3`), `12.2`–`12.8` (cancelling) and `36`–`60`
  (`s* = 5`) below that.  Both losses are in `varrho <= 2 T` and in the
  constants, none in the routes.
* **Against `eq:gapM`:** better by `2.8e17`, `4.4e23`, `3.5e121`, `1.3e138` on
  the four `s* = 3` designs and by `5.3e52`, `8.4e69`, `3.1e208`, `7.5e260` at
  `s* = 5`.  The two exceptions are the two cancelling designs whose window is
  empty at every temperature (`r = 0.9, b = 1` and `r = 0.99, b = 3`), where
  `eq:gapM` reduces to `eq:gapM0` and is `1.6` times larger than (T2).  Against
  `eq:gapM0` itself (`= 7.97e-05` at `s* = 3`, `3.38e-05` at `s* = 5`), which
  is available only at the `4.5%`–`8%` of temperatures where the window is
  empty, (T2) is a factor `3.2` (`s* = 3`) and `12`–`20` (`s* = 5`) worse — but
  holds at *every* temperature.  Against the two-level bound (T-U) of
  `NOTES_flow`: better by `5.7`, `12.8`, `14.9`, `14.9` at `s* = 3` and by
  `9.4`, `49`, `31`, `49` at `s* = 5`, and, unlike (T-U), valid without
  `ass:unimodal`.  On the cancelling pairs (T-U) is *larger* than (T2) by
  `25`–`85`, which is precisely the regime in which it is not a bound at all
  (section 6.4).
* **Against the block-tensorisation certificate of `NOTES_tensorisation`:**
  worse by a factor `17` to `67`.  That certificate is `1/{4(4p + 1/kappa_S)}`
  with `kappa_S` the worst fibre gap, computed here independently
  (`fibre_gap_min` in `gapA_validate.py`); it reproduces their numbers
  (`min cert/Gap = 0.076, 0.048, 0.044, 0.053, 0.085, 0.076` against their
  `0.076, 0.047, 0.044, 0.053, 0.085, 0.076`).  It is sharper — by a factor `17` to `67` — and it leaves `kappa_S`
  unbounded; (T2) bounds everything it names.

### 6.2 Hypothesis (B), design by design

| design | (B1) exchange descent | (B2) diminishing gains | greedy gains `g_m` |
|---|---|---|---|
| equicorr rho=0.00, s*=3 | yes | yes | 32.17, 15.79, 7.24 |
| equicorr rho=0.50, s*=3 | yes | yes | 52.69, 17.32, 3.21 |
| equicorr rho=0.85, s*=3 | yes | yes | 74.71, 7.65, -0.81 |
| equicorr rho=0.95, s*=3 | yes | yes | 84.26, 1.64, -2.35 |
| **cancelling r=0.90 b=1** | yes | **no** | **-1.83, +1.17** |
| **cancelling r=0.90 b=3** | yes | **no** | **-0.63, +37.83** |
| **cancelling r=0.95 b=3** | yes | **no** | **-1.86, +20.66** |
| **cancelling r=0.99 b=3** | yes | **no** | **-2.85, +0.89** |
| equicorr rho=0.00, s*=5 | yes | no | 18.32, 13.97, 13.82, 22.09, 6.11 |
| equicorr rho=0.50, s*=5 | yes | yes | 40.74, 17.62, 16.70, 13.00, 6.10 |
| equicorr rho=0.85, s*=5 | **no** | no | 65.58, 14.03, 12.70, 4.64, 6.09 |
| equicorr rho=0.95, s*=5 | yes | no | 77.94, 19.21, 4.14, 4.86, -0.98 |

`kappa log p = 5.278`.  **(B2) is exactly the cancelling-pair test**: there the
first variable is worthless alone (`g_0 < 0`) and the second is worth `+37.8`
once the first is in.  (B2) also fails for three `s* = 5` designs, so it is
*sufficient and not necessary* — but it fails there only mildly
(`Gamma_lambda <= 4.26`), and at `rho = 0, s* = 5` its failure is real: the
true gap dips to `3.70e-03` at `lambda = 0.334`, a factor `4.83`, and
`Gamma_lambda = 3.58` there.

### 6.3 `Gamma_lambda` degrades exactly where it must, and only there

`deficit(lambda) := (1/4p)/Gap(lambda)` is the true loss.  Excerpt
(`gapA_table.py`; the full table is in the script output):

| design | lambda | Gap | deficit | Gamma | Gamma/deficit |
|---|---|---|---|---|---|
| cancelling r=0.9 b=3 | 0.001 | 1.795e-02 | 1.00 | 1 | 1.00 |
| cancelling r=0.9 b=3 | 0.168 | 7.346e-03 | 2.43 | 2.89 | 1.19 |
| cancelling r=0.9 b=3 | 0.251 | 3.927e-04 | 45.5 | 67.4 | 1.48 |
| cancelling r=0.9 b=3 | 0.334 | 8.761e-05 | 204 | 242 | 1.19 |
| cancelling r=0.9 b=3 | 0.667 | 4.397e-05 | 406 | 298 | 0.73 |
| cancelling r=0.9 b=3 | 1.000 | 2.909e-05 | 614 | 367 | 0.60 |
| cancelling r=0.95 b=3 | 0.251 | 1.485e-02 | 1.20 | 1 | 0.83 |
| cancelling r=0.95 b=3 | 0.334 | 4.109e-03 | 4.35 | 5.08 | 1.17 |
| cancelling r=0.95 b=3 | 0.584 | 5.046e-05 | 354 | 580 | 1.64 |
| cancelling r=0.95 b=3 | 1.000 | 1.122e-05 | 1592 | 1256 | 0.79 |
| equicorr rho=0, s*=5 | 0.334 | 3.696e-03 | 4.83 | 3.58 | 0.74 |
| equicorr rho=0, s*=5 | 1.000 | 1.796e-02 | 0.99 | 1 | 1.00 |

`Gamma_lambda / deficit(lambda)` stays in `[0.60, 1.71]` on the cancelling
family across three orders of magnitude, and in `[0.43, 4.14]` on the `s* = 5`
equicorrelated designs.  On the four `s* = 3` equicorrelated designs
`Gamma_lambda = 1` at every temperature and `deficit <= 1.11`.  **So
`Gamma_lambda` is a faithful, computable proxy for the true deficit, it is
temperature-selective (it is `1` below `lambda = 0.17` even on the worst
cancelling design), and it is what the brief asked to identify.**  The tie
index `T_lambda`, by contrast, is `2.00`–`2.01` on the benign cancelling
designs and rises with `Gamma` only because it is evaluated at
`c = Gamma_lambda`; evaluated at `c = 1` its maximum over the grid is `2.28`,
`2.65`, `2.50`, `2.14`, `2.43`, `2.45`, `2.00`, `2.66` over the eight `(r,b)`
pairs of section 6.4 — that is, `T_lambda(1) <= 2.66` on the whole cancelling
family, while `Gamma_lambda` reaches `3645`.

### 6.4 The cancelling family, and the failure of (T-U)

`gapA_adversarial.py`, 25 temperatures:

| r | b | min Gap | min Gap x 4p | min B | worst B/Gap | max Gamma | (T-U) at argmin | **(T-U)/Gap** | modes |
|---|---|---|---|---|---|---|---|---|---|
| 0.90 | 3 | 2.909e-05 | 0.0016 | 3.79e-07 | 3.8e-03 | 367 | 9.30e-06 | 0.320 | 2 |
| 0.90 | 8 | 2.627e-05 | 0.0015 | 3.08e-07 | 3.2e-03 | 451 | 9.30e-06 | 0.354 | 2 |
| 0.95 | 3 | 1.122e-05 | 0.0006 | 1.11e-07 | 3.8e-03 | 1256 | 9.30e-06 | 0.829 | 2 |
| 0.95 | 8 | 9.648e-06 | 0.0005 | 9.35e-08 | 2.7e-03 | 1491 | 9.30e-06 | 0.964 | 2 |
| **0.97** | **5** | 7.817e-06 | 0.0004 | 6.95e-08 | 3.2e-03 | 2007 | 9.30e-06 | **1.190** | 2 |
| **0.97** | **8** | 6.818e-06 | 0.0004 | 5.88e-08 | 2.7e-03 | 2370 | 9.30e-06 | **1.364** | 2 |
| **0.99** | **8** | 5.442e-06 | 0.0003 | 3.83e-08 | 3.3e-03 | 3645 | 9.30e-06 | **1.709** | 2 |

Two things.  First, `B(lambda)/Gap(lambda)` never leaves
`[2.7e-03, 3.9e-03]` while the true gap falls to `3.0e-04` of `1/(4p)`, a
factor `3281`.  Second, the
formula (T-U) of `NOTES_flow`, evaluated where its hypothesis `ass:unimodal`
fails (two local maxima in every row), **exceeds the true gap by up to `1.71`**:
it is not merely unproved there, it is false.  This is the sharp version of
"`ass:unimodal` cannot be dropped".

---

## 7. What (B) is, and whether it is the right hypothesis

**(B) is checkable from `(X, beta*, sigma_0, g)`**: it is two statements about
the residual sums of squares `D_A` of the `2^{s*}` submodels of `gamma*`, i.e.
about the best-subset regression path on the influential set.  It does not
mention `lambda`, `kappa`, `pi^(lambda)`, its modes or its conductance, and it
is the same hypothesis at every temperature.  In the paper's own language,
(B2) says: reading `g_m = psi_{j_{m+1}}(P_m)` from (2), *the greedy log-gains
along the greedy chain are non-increasing*, which is `lemma:forward`\ref{it:fs_a}
strengthened from "there is a floor" to "the floors are ordered".

**It is exactly the hypothesis the cancelling pair kills.**  `NOTES_tensorisation`
section 6 records the same phenomenon from the other side:
`psi_j(empty) = -0.63, -3.68` and `psi_j(S \ j) = +40.9, +37.8` at `r = 0.9`,
i.e. the log-gains *increase* along the chain, which is `not (B2)`.  And it is
of a different kind from `ass:B` and `ass:D`: the cancelling pair at `r = 0.9`
has a *larger* `nu` and a *smaller* `omega(X)` than the equicorrelated design at
`rho = 0.85` (their table), so no sharpening of `nu` or `omega(X)` can separate
them, while `g_0 = -0.63 < g_1 = 37.8` separates them at a glance.

**Two honest qualifications.**

1. *(B) is sufficient, not necessary.*  It fails on three of the four `s* = 5`
   equicorrelated designs whose gap is healthy.  What the theorem uses is
   `Gamma_lambda`, which is computed; (B) is one checkable condition forcing
   `Gamma_lambda = 1`.  A weaker sufficient condition would be
   "`phi_lambda` is quasiconcave for the `lambda` of interest", but that is
   temperature-dependent and closer to a statement about the target.
2. *(B) is still a landscape hypothesis, one level down.*  It is a hypothesis
   about the shape of the map `A -> D_A` on `2^{gamma*}`, not about the modal
   structure of `pi^(lambda)` on `2^{W_lambda}`.  I regard the move as real —
   it is lambda-free, it is a standard statistical object, and it is *not*
   needed for the theorem, only for the corollary `Gamma_lambda = 1` — but I do
   not claim to have escaped landscape hypotheses altogether, and section 9
   explains why nobody can.

---

## 8. The cold end, and what `lambda_star < 1/2` is worth

**Verified numerically** (`gapA_lambdastar.txt`), at `p = 14`, `n = 100`,
`kappa = 2`:

| design | max psi | min greedy gain | lambda_cold | lambda_* (empirical) | ratio | meas{k_lambda>0} |
|---|---|---|---|---|---|---|
| equicorr rho=0.00 | 45.76 | 7.24 | 0.0577 | 1.823 | 31.6 | 0.920 |
| equicorr rho=0.50 | 54.86 | 3.21 | 0.0481 | 4.106 | 85.4 | 0.935 |
| equicorr rho=0.85 | 78.09 | -0.81 | 0.0338 | inf | inf | 0.950 |
| equicorr rho=0.95 | 88.39 | -2.35 | 0.0299 | inf | inf | 0.955 |
| cancelling r=0.90 b=3 | 43.74 | -0.63 | 0.0603 | inf | inf | 0.905 |
| cancelling r=0.95 b=3 | 24.78 | -1.86 | 0.1065 | inf | inf | 0.830 |

`lambda_cold` here is `(kappa-1)log p / max psi`, the exact analogue of
Proposition C, and the last column is the fraction of `[0,1]` on which the
window of `def:cpm` is non-empty.

**So the answer to items 4 and 5 of the brief is: yes and no.**  Yes,
`lambda_star` is bounded by an absolute constant, and the constant is `1/2`,
sharp, with no dependence on `n, p, s*, nu, C_beta`, and unchanged when the
exponent `3` is lowered to `2` (Proposition D).  No, this does not buy
anything structural.  The two ends do not meet: `lambda_cold` is of order
`kappa log p / (n ||beta*||^2/sigma_0^2)` and `lambda_star` of order
`(kappa+3) log p / (n nu^2 C_beta^2 / sigma_0^2)`, so
`lambda_star / lambda_cold >= const . ||beta*||^2/(nu^2 C_beta^2)`, which is at
least `s*/nu^2` — measured `31.6` and `85.4` on the two designs where
`lambda_star` is finite at `p = 14`.  The interval
`[lambda_cold, lambda_star]`, which carries 83%–96% of the path in these
experiments, is the whole problem, and damping the likelihood by a factor `2`
does not touch it, because `psi ~ n` while `kappa log p ~ log p`: at
`lambda = 1/2`, `lambda max_j psi_j` is still `4.3` to `8.4` times
`kappa log p` on the six designs of the table.

The exponent question separately: the change from `3` to `2` in
`cor:critical`\ref{it:crit_under} is *legitimate* — `cor:prec` needs
`a = (s*+1)p^{-eta+1} <= 1/2`, and with `eta = 2` this is
`(s*+1) <= p/2`, which `ass:D` already gives with room
(`s* + 1 <= s <= n/(256 log p) <= p/(256 log p)`), and `prop:conc` then reads
`1 - O(s*/p)` instead of `1 - O(1/p)`, still `>= 3/4`.  It lowers
`lambda_star` by the factor `(kappa+2)/(kappa+3)` and no more; the supremum
`1/2` is unmoved (Proposition D, verified for `q = 2, 3, 4`).

---

## 9. Obstructions and open points

1. **The `2^{s*}` in `T_lambda` is not removed, and I know why.**
   `T_lambda(1) = max_A pibar{B : pibar(B) <= pibar(A)}/pibar(A)` is the
   effective number of states of comparable mass, and it is genuinely
   `max_m C(s*,m)` in the degenerate limit `lambda -> 0`, where the prior ties
   every subset of a given size: measured `2.96` at `s* = 3`
   (`max_m C(3,m) = 3`) and `9.89` at `s* = 5` (`max_m C(5,m) = 10`).  This is a
   cost of *counting*, not of a bottleneck: the tied states are joined by
   exchange edges of full weight and the chain crosses between them instantly
   (the true gap at `lambda -> 0` is `1.795e-02 = 1/(4p)`).  The classical
   device that makes ties free is Sinclair's *encoding*, which for the
   fixed-order flip route on a product measure gives `varrho <= 1` with no
   counting at all.  **I could not make an encoding work here**, for a reason I
   can state exactly: the encoding needs the map `(A,A') -> (z,eta)`,
   `eta = A xor A' xor z`, to be injective at a fixed edge, and it is, for
   *flip* routes; but the flip route is catastrophic on collinear designs
   (measured `varrho = 1.07e27` at `rho = 0.85, s* = 3` and `4.2e42` at
   `rho = 0.95, s* = 5`), because it must pass through the empty combination of
   a collinear group.  The exchange route repairs that but destroys the
   injectivity: at an exchange edge, `(z, eta)` determines `A u A'`, `A n A'`
   and the *sizes* of `A \ A'` and `A' \ A`, but not which elements of
   `z xor eta` belong to which side, and the ambiguity is `C(|z xor eta|, q)`.
   **Open:** an encoding for matched-exchange routes on `2^{gamma*}`, or a
   proof that the multiplicity is compensated by the masses.  This is the one
   place where a poly(`s*`) bound is within reach.
2. **No bound of this shape can be unconditional.**  The cancelling pair makes
   the exact gap fall by `3281` at `p = 14` (section 6.4) while every hypothesis
   of Section 3 that is monotone in the design gets *better*, so some quantity
   must degrade.  In (T2) it is `Gamma_lambda`, and section 6.3 shows it
   degrades there, tracks the size of the collapse, and equals `1` everywhere
   else.  I regard this as settled rather than open.
3. **(B) is sufficient, not necessary, and I have no sharp characterisation.**
   The exact condition for `Gamma_lambda = 1` at a given `lambda` is
   quasiconcavity of `phi_lambda` plus (B1); the exact condition for
   `Gamma_lambda = 1` at *every* `lambda` is (B1) plus quasiconcavity of
   `m -> -theta log D_{P_m} - c m` for every `theta, c > 0`, which is (B2).
   So (B) is sufficient for the E-route to be dip-free at all temperatures
   simultaneously, and (B2) is also *necessary* for `phi_lambda` to be
   quasiconcave at every `lambda >= 0`: if `g_a < g_b` for some `a < b`, any
   `lambda` with `lambda g_a < kappa log p < lambda g_b` makes increment `a`
   negative and increment `b` positive, which is a valley.  Restricted to
   `lambda in [0,1]` that argument needs `g_b > kappa log p`, so on `[0,1]`
   (B2) is sufficient and only partially necessary — which is exactly the slack
   that the three `s* = 5` designs violating (B2) with a healthy gap exhibit.  **Conjectured, not proved:** a quantitative version,
   `Gamma_lambda <= exp{max_{a<m<b}(min(phi(a),phi(b)) - phi(m))^+}` with the
   maximum over the *layer optima* rather than over the greedy chain, would be
   both necessary and sufficient up to (B1); I did not prove that the layer
   optima are reachable without (B1).
4. **The within-layer clause (B1) is the weakest part.**  It is a hypothesis
   about a *local search* succeeding, which is uncomfortable.  It holds on 11
   of the 12 designs tested and its failure at `rho = 0.85, s* = 5` costs a dip
   of only `4.26`.  A cleaner replacement would be a statement about pairs of
   covariates in `gamma*` — the pair-amplification `Xi` of
   `NOTES_tensorisation` Proposition A is the obvious candidate — but
   `tensor_numerics.py --screen` already shows `Xi` is not a reliable
   predictor, and I did not find an implication either way.
5. **The bound is `Theta(1/p)` only up to `2^{s*}`, and the truth is
   `Theta(1/p)` flat.**  Measured, `4p Gap` lies in `[0.787, 1.728]` on every
   benign design and every temperature, the single exception being the one
   genuine dip (`rho = 0`, `s* = 5`, `lambda = 0.334`, where it is `0.207`);
   (T2) delivers `1/(8p(2s_0+4s*)2s*T)`.
   Closing the remaining `poly(s_0, s*) . 2^{s*}` is the same problem as item 1,
   plus the `s_0` that `NOTES_flow` item 7 already identifies as an artefact of
   the path method.
6. **What (T2) does *not* need, and what that buys.**  It does not need
   `cor:critical`, `cor:prec`, `prop:conc`, `lambda_star`, `Psi_r`, `W_lambda`,
   `Theta_lambda`, or `ass:unimodal`; Section 3 enters only through
   `prop:a0`, which supplies `a_0` for `ass:twosided`, and through
   `s* <= s_0`.  In particular (T2) is a *single* statement covering
   `[0,1]`, whereas `sec:gap` covers `[lambda_star, 1]` and `sec:alltemp`
   covers `[0,1]` at the price of `ass:unimodal`.  The concentration statement
   `prop:conc` is lost below `lambda_star` — as it must be, `pi^(0)(gamma*)`
   being `p^{-kappa s*}` — so the mixing-time corollary `eq:tau` must replace
   it.  `lemma:minmass` is available at every `lambda`, and combined with
   `1 <= |M(s_0)| max_gamma pi^(lambda)(gamma)` and `|M(s_0)| <= (p+1)^{s_0}`
   it gives
   `log[1/min_gamma pi^(lambda)] <= 2{kappa s_0 log p + lambda(n+s_0)alpha log p}
   + s_0 log(p+1)`, which fed into `eq:sandwich` together with (T2) yields a
   mixing-time bound of the same order as `eq:tau` at every temperature.  I did
   not optimise the constant.
7. **A negative I checked and discarded.**  Continuation in `lambda` (bounding
   `Gap(P_lambda)` by `Gap(P_{lambda'})` for nearby temperatures, via a
   Diaconis--Saloff-Coste comparison) gives nothing: the comparison constant is
   `exp{|lambda - lambda'| osc(log L_n)}` and `osc(log L_n) = Theta(n)`, so
   reaching `lambda_star` from `lambda = 0` needs `Theta(n)` steps and costs
   `e^{Theta(n)}`.  The temperature must be handled globally, which is what
   `m_lambda` does.

---

## 10. Reproducing

    python3 gapA_probe.py [--quick]        ->  gapA_probe.txt      (~3 min)
    python3 gapA_validate.py --quick --big ->  gapA_validate.txt   (~25 min)
    python3 gapA_lambdastar.py             ->  gapA_lambdastar.txt (~2 min)
    python3 gapA_adversarial.py                                     (~5 min)
    python3 gapA_table.py                                           (reads gapA_validate.txt)

`gapA_probe.py` and `gapA_validate.py` import `window_numerics` for the design,
the marginal likelihood and the model constants, and `flow_numerics` for the
state space, the exact gap, `decided_sets` (for the `eq:gapM` column) and
`two_level` (for the exact congestion of the finished ensemble).  Neither is
modified.
