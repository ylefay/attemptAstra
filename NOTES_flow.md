# Can a multicommodity flow remove the commitment penalty?

Working notes.  All numbers below are produced by `flow_numerics.py`:
`python3 flow_numerics.py` writes `flow_numerics.txt` (p = 14, s_0 = 4,
s* = 3, 40 temperatures, plus the adversarial designs) and
`python3 flow_numerics.py --big` writes `flow_numerics_big.txt`
(p = 14, s_0 = 6, s* = 5, 20 temperatures).  Nothing is asserted that has not
been computed.  Nothing in `binary_variable_selection3.tex` or `bvs3_app.tex`
has been touched.

---

## 0. Answer, up front

**No — not in the form the brief proposes, and not in the regime the paper
cares about.  Splitting the demand over the orders in which the undecided
coordinates are flipped is worth a factor between 1.7 and 8.5.  It is not what
removes `(1 + k_lambda Theta_lambda)^{k_lambda}`.**  (There is one regime where
the *optimal* flow does buy something no single path can; it is described at
the end of this section and in section 3, and it is exactly the regime in
which I can prove nothing.)

Splitting the demand between two states over the *orders* in which the
undecided coordinates are flipped — the textbook cure for a commitment
penalty — was computed exactly against the deterministic fixed-order path on
the same graph.  Over every design and every temperature tested, the ratio
(fixed order) / (uniform random order) is

| design | fixed order H1 | uniform random order H4 | ratio |
|---|---|---|---|
| rho = 0.00 | 5.55 | 3.23 | 1.7 |
| rho = 0.50 | 36.6 | 19.7 | 1.9 |
| rho = 0.85 | 1.07e27 | 4.92e26 | 2.2 |
| rho = 0.95 | 8.35e33 | 4.19e33 | 2.0 |
| s* = 5, rho = 0 | 519.9 | 142.8 | 3.6 |
| s* = 5, rho = 0.85 | 1.62e35 | 1.90e34 | 8.5 |

A congestion of `4.9e26` is not better than one of `1.07e27` in any sense that
matters.  Fractionality is a bounded-factor device; the penalty is a power of
`p`.

**What does remove the penalty is a different pair of changes, both of which
are compatible with a single canonical path:**

1. **Do not contract the residual classes.**  Keep the block/funnel
   construction that was in `binary_variable_selection2.tex`
   (`Section "A gap at every temperature"`, `lemma:funnel`,
   `def:flowens`): funnel each block
   `B_J = {gamma : gamma cap W_lambda = J}` onto its own representative
   `v_J = m^-_lambda u J`, so that there are `2^{k_lambda}` fixed points
   instead of one.  This is where "no commitment" actually happens.  It is
   necessary but **not sufficient**: with the increasing-order flip route of
   `def:flowens` the congestion is still `1.07e27` at rho = 0.85.
2. **Route the second level uphill, and let it use the exchange move.**  Two
   routes achieve this; both are single paths, no flow needed:
   - **H2**, pair `J1 \ J2` with `J2 \ J1` and use exchanges, then flip the
     remainder (this is `swap_congestion` in `window_numerics.py`);
   - **H3**, steepest ascent of `pibar_lambda` on `2^{W_lambda}` to its mode,
     over flips *and* exchanges — i.e. **route through the mode**, item 2 of
     the brief.

   Exact hypercube congestion `varrho_lambda` (max over 40 temperatures):

   | routing | rho=0 | rho=0.5 | rho=0.85 | rho=0.95 |
   |---|---|---|---|---|
   | H1 flips, increasing order (single path)  | 5.55 | 36.6 | 1.07e27 | 8.35e33 |
   | H2 exchange pairing (single path)         | 5.55 | 2.04 | **1.00** | **1.00** |
   | H3 ascent to the mode (single path)       | **1.00** | **1.00** | **1.01** | **1.00** |
   | H4 uniform flip order (FLOW)              | 3.23 | 19.7 | 4.92e26 | 4.19e33 |
   | H5 random matching + order (FLOW)         | 2.87 | 1.00 | 1.00 | 1.00 |
   | H6 optimal multicommodity flow, LP (FLOW) | **1.00** | **1.00** | **1.00** | **1.00** |

   H3 is the only *single* route that is uniformly `~1`; H2 is uniformly `~1`
   except in the orthogonal design, where it costs 5.55.  The optimal flow H6
   equals 1 everywhere, so nothing better than `varrho = 1` exists and H3
   attains it — **there is no room left for fractionality here.**

So the honest summary is: the commitment penalty is removed — the gap bound
improves by 17 to 137 orders of magnitude on the four `s* = 3` designs and by
52 to 259 orders on the four `s* = 5` designs — but multicommodity flow is not
the instrument.  The instruments are
the block decomposition of the state space and the exchange move, which the
chain already has and which the current `M_lambda` deliberately never uses
(`lemma:Mvalid`\ref{it:mv_i}).

**One caveat in the other direction.**  There *is* a regime where fractionality
does something no single path can, and it is exactly the regime in which the
ascent route H3 becomes unavailable: when the block measure is not unimodal.
At `p = 14, s_0 = 6, s* = 5`, rho = 0.85, `pibar_lambda` has two local maxima,
H3 does not exist, the best single route (H2) has congestion `2.9e5`, and the
optimal flow has congestion `1.00`.  See section 3.  But I have no analytic
bound on the optimal flow in that regime, so it buys numbers and not a
theorem.

---

## 1. Why the current bound loses: a sharper diagnosis than "near-neutrality"

The draft motivates `W_lambda` as the set of *nearly neutral* covariates, for
which `lambda psi_j ~ kappa log p`, and says a nearly neutral coordinate is
easy for the chain.  That is true but it is not the case that hurts.
`def:cpm` puts `j` in `W_lambda` when

    lambda psibar_j > kappa log p - log(2 s*)      (not certainly out)
    lambda psiunderbar_j < kappa log p + log(2 s*) (not certainly in),

i.e. the test uses `psiunderbar_j` (min over contexts) on one side and
`psibar_j` (max over contexts) on the other, and `Theta_lambda` is then
defined from `psibar_j`.  So a covariate with `psiunderbar_j` just below the
threshold and `psibar_j` far above it is classified undecided and charged
`Theta_lambda = exp(lambda psibar_j - kappa log p)`, which is enormous.  It is
not nearly neutral at all: it is *always worth having*, but by an amount that
depends violently on which of its collinear partners are already in.

The numbers say this is the dominant effect, and it is present even at rho = 0:

| design | `psiunderbar_j` | `psibar_j` | max `Theta_lambda` | max `max(1,Theta)(1+k Theta)^k` |
|---|---|---|---|---|
| rho=0.00 | 31.9, 6.04, 4.44 | 43.6, 18.6, 7.24 | 6.05e5 | 8.87e17 |
| rho=0.50 | 26.5, 10.6, 3.21 | 52.7, 33.0, 33.2 | 1.29e12 | 1.67e24 |
| rho=0.85 | 6.34, 1.33, -0.81 | 74.7, 67.6, 67.0 | 1.43e30 | 1.13e122 |
| rho=0.95 | -1.17, -1.89, -2.35 | 84.3, 83.4, 83.1 | 2.00e34 | 4.32e138 |

(`kappa log p = 5.28`.  The last column is the exact maximum over the 40
temperatures of the factor `eq:gapM` divides by; the maxima of `Theta_lambda`
and of `k_lambda` are not attained at the same temperature, so it is not the
product of the neighbouring entries.)

A second, independent defect follows: rule
`M:undec` **deletes** `min(gamma cap W_lambda)`.  When the covariate is in
fact always worth having, that step goes downhill by `Theta_lambda`, and the
fixed point `m^-_lambda` of `M_lambda` is a state of tiny mass.  The exact
Sinclair congestion of the ensemble `M_lambda` actually generates is

| design | exact `rho(T_lambda)` | gap bound `1/(rho ell)` | true gap |
|---|---|---|---|
| rho=0.00 | 4.20e6 | 2.38e-8 | 1.79e-2 |
| rho=0.50 | 7.69e4 | 1.62e-6 | 1.66e-2 |
| rho=0.85 | 1.52e32 | 8.23e-34 | 1.62e-2 |
| rho=0.95 | 4.95e35 | 2.52e-37 | 1.80e-2 |

so the loss is in the *ensemble*, not only in the estimates used to bound it:
`Theorem "A gap at every temperature"` gives `7.1e-127` at rho = 0.85 and the
ensemble it analyses really does have congestion `1.5e32`.  Both are hopeless;
the bound is about 93 orders of magnitude looser than the object it bounds,
which is the `(1+k Theta)^k` counting on top of the `Theta` per step.

---

## 2. The construction I propose

It is the two-level ensemble of `binary_variable_selection2.tex`
(`def:flowens`) with segment (ii) replaced.  Notation as in
`binary_variable_selection3.tex`.

**Level 1 (funnel), unchanged.**  Let `F_lambda` be `M_lambda` with rule
`M:undec` deleted, i.e. rules `M:null`, `M:del`, `M:add` only, and
`F_lambda(gamma) = gamma` when none applies.  No rule touches a coordinate of
`W_lambda`, so `F_lambda` maps `B_J` into itself and its unique fixed point
there is `v_J = m^-_lambda u J`.  All of `lemma:Mvalid` applies verbatim with
the fourth term of `eq:VM` deleted; the orbit has at most `s_0 + 2 s* <= 2 s_0`
edges, every one a single flip, and `th:dropM` gives ratios `p^{-eta}` and
`1/(2s*)` — **the `Theta_lambda` case never arises**.  The precedent
computation in the proof of `cor:precM`, with the `m`-sum deleted, gives

    mu[Lambda^F(gamma)] <= 2 (2 s* + 1) mu(gamma)                      (F1)

and in particular, since `B_J = Lambda^F(v_J)`,

    mu(v_J) >= pibar_lambda(J) / {2 (2 s* + 1)},                       (F2)

which is `lemma:funnel`\ref{it:fn_iii} of v2.  Measured: the constant actually
needed in (F2), `max_J pibar_lambda(J)/mu(v_J)`, is 1.16, 1.37, 1.16, 1.20 at
rho = 0, 0.5, 0.85, 0.95, against the bound `2(2s*+1) = 14`; and the funnel
edges of the finished ensemble have congestion 58 to 67, against the bound
`8p(2s*+1) = 784` proved below.

**Level 2 (block routing), new.**  On `2^{W_lambda}` with vertices `v_J` the
chain offers two kinds of edge, and both are in `N(gamma)` of `eq:neigh`:

* a **flip** between `v_J` and `v_{J u {j}}`;
* an **exchange** between `v_{J u {j}}` and `v_{J u {j'}}` — these two states
  have the same size and differ by one swap, so this is exactly `N_2`.

Let `A_lambda` be steepest ascent of `pibar_lambda` over *both* kinds of edge:
`A_lambda(J)` is the neighbour of largest `pibar_lambda` if it exceeds
`pibar_lambda(J)`, and `J` otherwise.  Route `J1 -> J2` by the orbit of `J1`
under `A_lambda`, spliced with the reversed orbit of `J2` at their first common
state.  This is `route through the mode`.

**Assumption (U).**  `pibar_lambda` has exactly one local maximum on
`2^{W_lambda}` for the flip-and-exchange neighbourhood.

Under (U) the ascent has a unique fixed point and is a valid transition
function; the congestion estimate is (F3) in the proof below.  Measured:
`max_J pibar[Lambda^A(J)]/pibar(J)` is 1.97, 2.06, 2.73, 3.14 at rho = 0, 0.5,
0.85, 0.95, against the bound `2^k = 4, 8, 8, 8`; the resulting `varrho` is
1.00, 1.00, 1.01, 1.00; the ascent orbit is at most 3 steps long, against the
bound `2^k - 1 = 7`.

### Theorem (what this gives)

State it for an arbitrary flow at the second level, and then specialise; the
flow version of Sinclair's bound (`Sinclair1992`, Corollary 6') is what
`eq:sinclair` becomes when a pair's demand is split over several routes, with
`ell` the longest route carrying positive flow.

Let `f_lambda` be any flow on the graph `(2^{W_lambda}, flips + exchanges)`
that routes a demand `pibar_lambda(J1) pibar_lambda(J2)` between every ordered
pair, let

    varrho_lambda := max over directed edges (J,J') of
                     [ flow of f_lambda across (J,J') ]
                     / min{ pibar_lambda(J), pibar_lambda(J') }

be its congestion in the normalisation of `def:rholam`, and let `L_lambda` be
its longest route.

> **Theorem (T).**  Let `lambda in [0,1]`, assume `ass:twosided` with an
> exponent `eta` satisfying `p^{1-eta} <= 1/2`, and work on the good event.
> Then, for every such flow,
>
>     Gap(P_lambda) >= 1 / { 8 p (2s*+1) max(1, s* varrho_lambda)
>                                (4 s_0 + 2 L_lambda) }.

> **Corollary (U).**  If in addition `pibar_lambda` has exactly one local
> maximum on `2^{W_lambda}` for the flip-and-exchange neighbourhood, then
> steepest ascent to that maximum is an admissible single route with
> `varrho_lambda <= 2^{k_lambda}` and `L_lambda <= 2^{k_lambda} - 1`, so
>
>     Gap(P_lambda) >= 1 / { 8 p (2s*+1) max(1, s* 2^{k_lambda})
>                                (4 s_0 + 2^{k_lambda+1}) }.            (T-U)
>
> With `k_lambda <= s*` this is
> `1/{ 8 p (2s*+1) s* 2^{s*} (4 s_0 + 2^{s*+1}) }`, a fixed rational function of
> `p, s_0, s*` times `4^{s*}` — **with no `p` in the exponent**.

*Proof of (T), with the constants.*  The ensemble uses two disjoint kinds of
edge (a funnel edge joins two states of the same block, a block edge joins two
different blocks, so no edge is of both kinds and the congestion is the maximum
of the two, not the sum).

*Funnel edges* `e = (gamma, F_lambda(gamma))`, a single flip.  The ensemble is
memoryless inside the block, so `e` can be traversed forwards only by a pair
whose source is a precedent of `gamma`, and backwards only by a pair whose
sink is; either way the demand it carries is at most
`mu[Lambda^F(gamma)] . 1 <= 2(2s*+1) mu(gamma)` by (F1).  The step is uphill,
so `Q_lambda(e) >= mu(gamma)/(4p)` by `eq:edge` *without* the factor `s_0`,
`F_lambda` using only flips.  Congestion `<= 8 p (2s*+1)`.

*Block edges* `e = (v_J, v_{J'})`, a flip or an exchange between two states of
size at most `s*`.  The demand is `<= varrho_lambda min{pibar(J), pibar(J')}`
by the definition of `varrho_lambda`, and
`Q_lambda(e) >= min{mu(v_J), mu(v_{J'})}/(4 p s*)
            >= min{pibar(J), pibar(J')} / {4 p s* . 2(2s*+1)}`
by `eq:edge` — with `s*` rather than `s_0`, legitimate because `|v_J| <= s*` —
and (F2).  Congestion `<= 8 p s* (2s*+1) varrho_lambda`.

Hence `rho <= 8 p (2s*+1) max{1, s* varrho_lambda}`.  A route is two funnel
orbits (`<= 2 s_0` edges each) and one second-level route (`<= L_lambda`
edges), and the spliced version of the latter is at most `2 L_lambda`, so
`ell <= 4 s_0 + 2 L_lambda`.  Sinclair's Corollary 6', with `lemma:kernel` for
the identification with the absolute spectral gap, gives (T).

*Proof of the corollary.*  Under unimodality the steepest-ascent map
`A_lambda` has a unique fixed point and `pibar_lambda` strictly increases along
it, so it is a valid transition function and its orbits terminate within
`2^{k_lambda} - 1` steps.  For the memoryless ensemble it generates, the demand
across `(J, A_lambda(J))` is at most `pibar[Lambda^A(J)]`, the denominator is
`min{pibar(J), pibar(A_lambda(J))} = pibar(J)`, and every precedent of `J`
carries mass at most `pibar(J)` because every step is uphill, so

    varrho_lambda <= max_J pibar[Lambda^A(J)]/pibar(J) <= 2^{k_lambda}.  (F3)

QED

Two sanity checks on (T-U).  At `k_lambda = 0` there are no block edges and it
reads `Gap >= 1/{32 p s_0 (2s*+1)}`, which is exactly `eq:gapM0` — the new
bound degrades continuously into the old one.  And measured against the two
estimates the proof makes: the funnel edges of the finished ensemble have
congestion 58 to 67 against the bound `8p(2s*+1) = 784`, the block edges have
congestion 59 to 99 against the bound `8 p s* (2s*+1) = 2352`, and the
precedent ratio in (F3) is 1.97 to 3.14 against the bound `2^k <= 8`.  At
`s* = 5` the funnel congestion is 60 to 69 against `8p(2s*+1) = 1232`, and the
precedent ratio is 4.02 to 4.52 against `2^k <= 32`.

### Where (T-U) beats `eq:gapM`

(T-U) beats `eq:gapM` iff

    max(1, s* 2^{k}) (4 s_0 + 2^{k+1})  <  4 s_0 max(1,Theta_lambda)
                                            (1 + k Theta_lambda)^{k}.

For `k >= 1` and `2^{k+1} <= 4 s_0` the left side is `~ 8 s_0 s* 2^k`, so the
condition is essentially `Theta_lambda >~ 2 s* 2^{k}`, i.e.

    lambda psibar_j  >=  kappa log p + log(2 s*) + k log 2 + O(1)
                         for the worst j in W_lambda.

Note what this says: a covariate is in `W_lambda` precisely because
`lambda psiunderbar_j < kappa log p + log(2 s*)`, so the condition is that
`psibar_j` exceeds `psiunderbar_j` by more than `k log 2 + O(1)` — that is,
(T-U) wins exactly when the undecided log-gains are *context dependent*, which
is the case `eq:Thetabound` describes with `R = psibar_j / psiunderbar_j > 1`.
In terms of `R`, `eq:Thetabound` gives `(1 + k Theta)^k <= p^{kappa (R-1) s*}
(s*)^{O(s*)}`, and it is a *lower* bound of the same shape that makes the
current theorem lose; (T) wins as soon as

    kappa (R - 1) s* log p  >  2 s* log 2 + log s*,
    i.e.  R  >  1 + (2 log 2)/(kappa log p) + o(1).

For `p = 14`, `kappa = 2` that is `R > 1.27`, and the measured `R` is 3.08
already in the orthogonal design.  Numerically, (T-U) versus `eq:gapM`:

| design | `eq:gapM` | (T-U) | (T-U)/`eq:gapM` | exact cong. of that ensemble | true gap |
|---|---|---|---|---|---|
| rho=0.00, max k = 2 | 8.99e-23 | 4.43e-6 | 4.9e16 | 8.62e-4 | 1.79e-2 |
| rho=0.50, max k = 3 | 4.79e-29 | 1.66e-6 | 3.5e22 | 8.69e-4 | 1.66e-2 |
| rho=0.85, max k = 3 | 7.07e-127 | 1.66e-6 | 2.3e120 | 8.66e-4 | 1.62e-2 |
| rho=0.95, max k = 3 | 1.84e-143 | 1.66e-6 | 9.0e136 | 1.31e-3 | 1.80e-2 |

((T-U) is evaluated at each temperature and minimised, as the other columns
are.  The fifth column is `1/(rho ell)` for the *actual* ensemble (T-U)
analyses, so (T-U) itself is loose by a factor 195 to 786 — all of it in the two `2^k`
estimates and in the constant `2(2s*+1)`, none of it in the construction.)

The one regime where (T-U) is worse is `k_lambda >= 1` with
`Theta_lambda = O(1)` — genuinely neutral covariates with context-free
log-gains.  There `eq:gapM` costs `(1+k)^k <= (s*+1)^{s*}` and (T-U) costs
`s* 4^{s*}`; the two are comparable, and one can of course take the maximum of
the two bounds.

---

## 3. Assumption (U): needed, testable, and not free

**(U) holds throughout the primary experiment.**  `pibar_lambda` has exactly
one local maximum for flips+exchanges at all 40 temperatures and all four
designs.  With flips only it has up to **3** local maxima at rho = 0.85 and
0.95.  So the exchange move is what makes (U) true, which is the same message
as `unimodal_test.txt` one level up.

**(U) cannot be dropped.**  Take the cancelling pair: `X_0, X_1` correlated at
`r`, `beta* = (b, -b, 0, ..., 0)`.  Neither column explains anything alone,
the pair explains everything, so `pibar` has a genuine bottleneck between
`J = {}` and `J = {0,1}`.  Exact gaps (`p=14, s_0=4`, 40 temperatures):

| r | b | min Gap | min Gap x 4p | # local maxima of `pi^(lambda)` on `M(s_0)` |
|---|---|---|---|---|
| 0.90 | 1 | 1.78e-2 | 0.995 | 1 |
| 0.90 | 3 | **9.77e-6** | **0.0005** | 2 |
| 0.99 | 1 | 1.79e-2 | 1.001 | 1 |
| 0.99 | 3 | **3.96e-6** | **0.0002** | 2 |
| 0.999 | 1 | 1.79e-2 | 1.000 | 1 |
| 0.999 | 3 | 1.79e-2 | 1.000 | 1 |

The gap really falls by a factor 2000–4500 exactly when the landscape becomes
bimodal.  No canonical-path or flow argument can avoid an assumption of this
kind, because the conclusion is false without it.

**(U) fails at `s* = 5` for two of the four designs, and the gap degrades with
it.**  Second experiment: `p = 14, s_0 = 6, s* = 5`,
`beta* = (1.20, 1.00, 0.85, 0.70, 0.45)`, `|M(s_0)| = 6476`, 20 temperatures
(`flow_numerics_big.txt`).

Congestion `varrho` (max over the temperatures at which the route exists):

| rho | (U) holds at | max #modes of `pibar` (flips+exch / flips only) | min Gap | `Gap(Pbar)` | H1 | H2 | H4 | H5 | H3 ascent | H6 LP |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 15/18 | **2** / 2 | 3.97e-3 | 4.04e-3 | 520 | 520 | 143 | 143 | n/a at 3 | 2.27 |
| 0.50 | 18/18 | 1 / 3 | 1.45e-2 | 1.48e-2 | 6.25e8 | 1.05e4 | 7.64e7 | 3.93e7 | **1.27** | 1.00 |
| 0.85 | 18/19 | **2** / 4 | 1.40e-2 | 1.40e-2 | 1.62e35 | 2.89e5 | 1.90e34 | 1.74e34 | n/a at 1 | 1.00 |
| 0.95 | 19/19 | 1 / 3 | 1.80e-2 | 1.80e-2 | 4.46e42 | 2.86e5 | 2.04e42 | 7.13e41 | **1.63** | 1.00 |

("(U) holds at m/n" counts the temperatures with `k_lambda > 0`.  The LP itself
failed to solve at 2, 0, 1 and 4 of those temperatures respectively; those are
excluded from its column and flagged in `flow_numerics_big.txt`.)

Gap lower bounds actually delivered (min over the temperatures at which the
route exists), against `eq:gapM` and the truth:

| rho | `eq:gapM` | two-level H2 | two-level H3 | two-level H6 | truth |
|---|---|---|---|---|---|
| 0.00 | 6.95e-60 | 1.90e-6 | 3.27e-4 | 1.17e-4 | 3.97e-3 |
| 0.50 | 3.37e-76 | 9.98e-8 | **3.14e-4** | 1.62e-4 | 1.45e-2 |
| 0.85 | 5.83e-215 | 3.63e-9 | 2.47e-4 | 2.09e-4 | 1.40e-2 |
| 0.95 | 3.73e-267 | 3.67e-9 | **3.45e-4** | 1.44e-4 | 1.80e-2 |

Four things to read off.

* Where (U) holds, the ascent route is still excellent at `k_lambda = 5`:
  `varrho = 1.27` and `1.63`, against `1.05e4` and `2.86e5` for the exchange
  pairing H2 and `6.2e8` and `4.5e42` for the increasing-order route H1.
* Where (U) fails, the gap really does degrade: at rho = 0 the true gap dips
  to `3.97e-3 = (1/4p)/4.5`, and the projected chain on `2^{W_lambda}` carries
  the whole dip (`Gap(Pbar) = 4.04e-3`).  So (U) tracks something real.
* Where (U) fails, **fractionality is the only thing that still works**: the
  optimal flow gives `varrho = 2.27` at rho = 0 and `1.00` at rho = 0.85,
  against `520` and `2.9e5` for the best single route — factors of 229 and
  `2.9e5`.  This is the one place in the study where a multicommodity flow
  does something no single path can.  I have no analytic bound on the LP
  optimum, so it does not yield a theorem; it says a theorem *could* exist if
  one could bound the conductance of `pibar_lambda` on `2^{W_lambda}`
  directly.
* Every *fixed-order* route degrades sharply with `k_lambda`: from `s* = 3` to
  `s* = 5` at rho = 0.85, H1 goes `1.07e27 -> 1.62e35` and H2 goes
  `1.00 -> 2.89e5`.  Only H3 and H6 stay `O(1)`.  Fractional order (H4, H5)
  tracks the fixed-order routes to within a factor 8.5.

**A sufficient condition for (U) in terms of the `psi`'s.**  Suppose
`W_lambda` can be ordered `j_1, ..., j_k` so that

    mu(v_J)  =  h(|J|) . prod_{j_i in J} w_i ,   w_1 >= ... >= w_k ,

with `h` log-concave (a *size-biased product measure*).  Then exchanges force a
local maximum to consist of the `|J|` heaviest coordinates, and log-concavity
of `h(m) w_m` in `m` forces the size to have a single local maximum; hence (U).
The family contains both extremes of the problem: `h(m) = c^m` is the product
case (`psi_j` context free — the near-neutral picture the draft has in mind),
and `h` concentrated on one size is the collinear-group case
(`psi_j(U) ~ 0` once a partner of `j` is in `U`).  I did **not** verify
numerically that the designs above are close to this family, and I am not
claiming they are; the observation is only that the family covers the two
regimes and gives a checkable sufficient condition.

---

## 4. Decomposition, and why the flow route is better

The alternative is the partition `{B_J}` and the decomposition theorem of
Jerrum, Son, Tetali and Vigoda (2004, Thm 1) / Madras and Randall (2002),

    Gap(P) >= (1/2) Gap(Pbar) min_J Gap(P|_{B_J}),

with `Pbar` the projection onto `2^{W_lambda}` and `P|_{B_J}` the restriction
with rejection.  I verified the inequality numerically at every temperature and
every design of both experiments (the ratio of the right side to the truth is
at most 0.025 over `{k_lambda > 0}`), so I am not misquoting the statement.  Both ingredients are
available:

* `min_J Gap(P|_{B_J})`: inside a block there are no undecided coordinates,
  `F_lambda` is a valid transition function with definite drops, and the
  argument behind `eq:gapM0` applies verbatim, giving
  `>= 1/{32 p s_0 (2s*+1)}`.  Measured: `1.786e-2 ~ 1/(4p)` in every case.
* `Gap(Pbar)`: the same `2^{k_lambda}`-state problem as before; Sinclair with
  the ascent ensemble gives `>= 2^{-(2 k_lambda + 1)}` under (U).  Measured at
  `s* = 3`: `1.6e-2` to `3.4e-2`, again `~ 1/(4p)`.  At `s* = 5` it drops to
  `4.0e-3` at rho = 0, which is where the true gap drops too.

**But this route is structurally worse by a factor `p`.**  `Gap(Pbar)` is
`Theta(1/p)` because the projected chain inherits the `1/(4p)` proposal rate
of the flips, and `min_J Gap(P|_{B_J})` is `Theta(1/p)` for the same reason, so
the product is `Theta(1/p^2)`, however sharply each factor is bounded.
Measured: `(1/2) Gap(Pbar) min_J Gap(P|_{B_J})` is between `1.4e-4` and
`3.0e-4`, against a truth of `1.6e-2` — a fixed factor of about `p^1` below
what the flow route delivers (`8.6e-4` to `1.6e-3`).  Since the whole point of
`sec:gap` is a bound of order `1/(p s_0^2)`, paying `p` twice is not
acceptable.

The plain **comparison** technique (Diaconis–Saloff-Coste) is worse still.
Comparing `P_lambda` with a chain targeting a product measure on
`2^{W_lambda}` costs `max_J mu(J)/mu'(J)`, and that ratio is exactly the
quantity `exp(k_lambda Delta^C_lambda)` that `prop:rhobound` of v2 already
charged; `window_numerics.txt` reports `Delta^C_lambda` up to 85.4 with
`k_lambda = 3`, so `k_lambda Delta^C_lambda` is up to 256 and the comparison
constant is `e^{256}`.  Comparison converts the problem into itself.

**Recommendation: the flow (really, the two-level path) route.**  It pays `p`
once, it reduces to `eq:gapM0` when `k_lambda = 0`, and its only new
ingredient is (U), which is checkable, is delivered by the exchange move, and
is known to be necessary.

---

## 5. Obstructions and open points

1. **(U) is an assumption on the target, not a consequence of
   `ass:A`–`ass:D`.**  I did not find a derivation from the paper's
   assumptions, and the cancelling-pair design shows none can exist in the
   form "always true".  What `ass:B` (sparse eigenvalue `nu`) and
   `ass:betamin` would have to buy is a *uniform* version of
   `lemma:forward`\ref{it:fs_a} — a gain floor for *every* missing influential
   covariate, not only the greedy one — and that is false at `rho = 0.95`,
   where `psiunderbar_j < 0` for all three.
2. **The `2^{k_lambda}` in (F3) is crude.**  The measured precedent ratio is
   at most 3.14 against a bound of 8.  A better bound needs a branching count
   for `A_lambda`, i.e. a bound on how many neighbours of `J` can ascend into
   `J`; that is `k + k(k-1)/2` in the worst case and gives nothing better than
   `2^k` after summing over distances.  Getting `poly(k)` here would need a
   definite drop along the ascent, which is exactly what an undecided
   coordinate does not provide.
3. **The ascent orbit length is bounded only by `2^k - 1`** (measured: 3 at
   `s* = 3`, 4 at `s* = 5`, against bounds 7 and 31).
   A route that both ascends and decreases `|J Delta J*|` would give `2k`, but
   (U) does not guarantee an ascending neighbour *towards* the mode.
4. **The LP optimum is not analytically bounded.**  `varrho_LP = 1` in every
   `s* = 3` case, and `1.00` to `2.27` at `s* = 5`, including where (U) fails.
   A theorem covering the non-unimodal case would have to bound the min-cut of
   `pibar_lambda` on `2^{W_lambda}` from below, i.e. assume a conductance
   bound on the block measure — a weaker and less interpretable hypothesis
   than (U), but strictly more general.  Worth stating as
   "`Gap >= 1/{8 p (2s*+1) max(1, s* varrho_lambda) (4 s_0 + 2^{k+1})}` for
   `varrho_lambda` the minimum congestion of any flow on `2^{W_lambda}`",
   with (U) as one sufficient condition for `varrho_lambda <= 2^{k_lambda}`.
5. **`M_lambda` never uses the exchange move, and that is the root cause.**
   `lemma:Mvalid`\ref{it:mv_i} advertises this as an economy (it saves one
   factor `s_0` in `eq:gapM`).  It costs `1e30`.  The exchange move is the
   only way for the chain to move between two collinear covariates without
   passing through a state that carries no mass, which is exactly what
   `exchange_test.txt` measures at the level of the full chain (gap ratio 97
   at rho = 0.95) and what H1 versus H2 measures at the level of the blocks
   (ratio `8.3e33` at rho = 0.95).
6. **The flow programme is numerically delicate, and I nearly published a
   wrong number because of it.**  The block masses span forty orders of
   magnitude (`pibar` down to `1.4e-38` at `s* = 5`), so with the raw
   formulation HiGHS's absolute feasibility tolerance lets it put a flow of
   `8e-8` on an edge of capacity `1.4e-38` and report an optimum of `0.61`
   for a flow whose true congestion is `5.7e30`.  `flow_numerics.py` now
   scales the columns by the capacities and the conservation row at `J` by
   `1/pibar(J)`, so all coefficients and right-hand sides lie in `[-1,1]`, and
   it reports the congestion of the flow the solver *returns* rather than the
   solver's objective value.  A self-test against the one-commodity-per-pair
   formulation, at log-spreads 2, 20 and 60, is printed at the top of
   `flow_numerics.txt`.  Even so the solver fails outright at a few of the
   `k_lambda = 5` temperatures; those are reported as UNAVAILABLE rather than
   silently dropped.  None of the H1-H5 numbers use the solver: they are
   direct floating-point evaluations of an explicit flow.
7. **Nothing here reaches the truth.**  The measured gap is `~1/(4p)`,
   flat in `lambda`, with no visible dependence on `s_0`, `s*` or
   `k_lambda`.  The best bound in these notes is `~1/(p s_0 (s*)^2 4^{s*})`.
   The `4^{s*}` and the `s_0` are artefacts of the path method, as
   `rem:gapM` already says; removing them needs tensorisation, not a better
   ensemble.

---

## 6. What I would change in the draft

Minimal, and it does not disturb `sec:main` or `sec:gap`:

* reinstate `lemma:funnel` and the block construction from v2 (they were
  correct; only segment (ii) was wrong);
* replace segment (ii) by the ascent `A_lambda` over flips **and** exchanges;
* add (U) as a hypothesis, with the cancelling-pair design as the remark
  showing it cannot be dropped, and the size-biased-product family as the
  sufficient condition;
* state (T), with (T-U) as its corollary;
* keep `rem:third` (the window is real) but retire the claim that the
  penalty is "the honest cost of a memoryless ensemble" — the ensemble here is
  memoryless too, per block, and pays `4^{s*}` instead of
  `p^{1 + kappa (R-1) s*}`.  The honest cost is of *contracting the residual
  classes to a single state* and of *refusing the exchange move*.
