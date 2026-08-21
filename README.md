# Birch–Swinnerton-Dyer: computational tooling

Pure-Python tools for computing both sides of the Birch–Swinnerton-Dyer
conjecture for elliptic curves over **Q**, built from scratch with no
dependencies — no Sage, no PARI, no sympy, no numpy.

**This does not prove BSD.** BSD is open. What is here computes the algebraic
rank by descent, computes the analytic rank from the L-function, and checks
that they agree on curves where both are reachable.

---

## What works, and how it was validated

| Tool | Does | Validation |
|---|---|---|
| `rank.py` / `descent.py` | Descent via 2-isogeny. Needs a rational point of order 2. | 10/10 exact ranks on congruent-number curves with independently known answers |
| `mw.py` | Rank *lower* bounds for any curve: point search + Néron–Tate height pairing | Regulators match published values to 5 s.f. on 37a1, 389a1, 5077a1 |
| `descent2.py` | **General** 2-descent (Birch–Swinnerton-Dyer quartic method). No 2-torsion needed. | dim Sel₂ = 1, 2, 3 exactly on 37a1, 389a1, 5077a1 |
| `complete2.py` | **Complete** 2-descent for curves with full rational 2-torsion. **No search box, no reduction theory** — the candidate set is exact. | 56/60 exact ranks; every Selmer count a power of 2; zero bug indicators |
| `bsd.py` | Runs the two engines together and cross-checks them | 60-curve sweep, **zero contradictions**, 55/60 exact ranks |
| `lfun.py` | The L-function from point counts + modular continuation | a_p matches newform 37.2.a.a; analytic ranks 1, 2, 3 recovered |
| `sha.py`, `bm.py`, `bm3.py` | Ш made explicit; the Brauer–Manin obstruction computed as a number | Σ_v inv_v(A) = 1/2 on five independent adelic points |

### Headline results

- **56/60 exact ranks** on `y² = x³ − n²x` for n = 1…60, by complete 2-descent
  (`complete2.py`) — rigorous, with no search box anywhere in the computation.
- **Zero contradictions** between independently written engines — across the
  60-curve sweep the point-search lower bound never once exceeded the descent
  upper bound, and every Selmer count came out a power of 2 as group structure
  demands.
- **One confirmed Ш obstruction, at n = 17**: complete 2-descent gives
  dim Sel₂ = 4 exactly; E(Q)[2] contributes 2; the L-function gives
  L(E,1) = 2.54 ≠ 0, so the analytic rank is 0, and Coates–Wiles (CM curves)
  forces rank = 0. Hence **dim Ш[2] = 2**, i.e. #Ш[2] = 4 — a perfect square,
  as Cassels–Tate requires.
- **17 and 42 are not congruent numbers** — no right triangle with rational
  sides has either as its area. For n = 42 complete 2-descent settles it alone,
  with no analytic input.
- **Analytic rank = algebraic rank** verified independently on 37a1 (1),
  389a1 (2), 5077a1 (3) — by two programs sharing no code and no concepts.

Telling a real obstruction from a tired search is most of the work. Of ten
apparent Ш gaps in the first sweep, five were merely a weak point search and
closed under a larger height bound. The remaining gaps at n = 37, 47, 53 are
also weak lower bounds, not Ш: those are primes ≡ 5, 7 mod 8, hence congruent
numbers of rank 1.

### A correction

An earlier version of this README claimed **two** Ш obstructions, at n = 17 and
n = 42. That was wrong. The n = 42 gap came from the weaker descent via
2-isogeny, whose excess measures Ш[φ] on the *isogenous* curve E′ = y²=x³+4n²x,
not Ш(E)[2]. Complete 2-descent shows Ш(E)[2] = 0 for n = 42. The n = 17 result
survives and is now rigorous rather than heuristic.

---

## What does not work

Stated plainly so nobody builds on sand.

1. **`descent2.py`'s upper bound is heuristic, not proof** — but only for
   curves whose 2-division cubic is irreducible. For curves with full rational
   2-torsion use `complete2.py` instead, which is exact and box-free. The
   quartic enumeration uses a search box with no proven reduction bound.
   It self-audits — Sel₂ is a group, so a class count that is not a power of 2
   proves the box was too small, and the tool escalates and refuses to report
   when it cannot stabilise. But a *missing subgroup* would still leave a power
   of 2, so completeness is not certified. Correct fix: Cremona's reduction
   theory for binary quartics. **Attempted and failed** — see `unfinished/enum2.py`.

2. **`lfun.py` must be handed the conductor.** Computing N needs Tate's
   algorithm. **Attempted and failed** — see `unfinished/`. `tate.py` does not
   terminate. `cond.py` scores 7/12 and its minimisation is wrong at p = 2, 3
   (Kraus's criterion has extra congruence conditions there); it can return
   N = 1, which is impossible.

3. Local solvability tests can return inconclusive at their depth cap. When
   that happens `complete2.py` counts the class as being **inside** Sel₂. That
   is the safe direction: it can only weaken the rank upper bound, never make
   it too small. (Dropping such classes instead — an earlier bug — produced
   `rank ≤ −1` on n = 59, which the power-of-2 self-audit caught.) Currently
   inconclusive on n = 47, 53, 59 of the sweep.
4. Canonical heights are numerical, with a heuristic threshold (1e−3 on the
   normalised Gram determinant) for independence.
5. `unfinished/enum2.py` is fast but **incomplete** — it misses quartics. Do not
   use it for upper bounds.
6. `onetors.py` handles imaginary quadratic K only. Real quadratic K needs
   fundamental units (the norm form is indefinite, so infinitely many elements
   share a norm); not implemented.

---

## Usage

```bash
python rank.py 0 0 0 -289 0      # 2-isogeny descent, Weierstrass a1 a2 a3 a4 a6
python descent2.py               # general 2-descent on the benchmark curves
python sweep2.py                 # complete 2-descent, 60-curve sweep (exact)
python mw.py                     # rank lower bounds via height pairings
python bsd.py                    # 60-curve combined sweep
python sha.py                    # two curves that violate the Hasse principle
python bm.py && python bm3.py    # the Brauer–Manin obstruction, computed
```

Python 3.8+. No dependencies.

---

## Why the two sides never meet

Everything in `descent.py`, `descent2.py` and `mw.py` computes the algebraic
rank — rational points, torsors, Selmer groups, local solvability.
Everything in `lfun.py` computes the analytic rank — point counts mod p fed
through a modular continuation. The two share no code and no concepts, and
they agree every time they are both asked.

For analytic rank ≤ 1 there is a bridge: Gross–Zagier computes L'(E,1) as the
canonical height of a Heegner point, so a simple zero hands you a generator,
and Kolyvagin's Euler system bounds the rest. There is no known analogue at
order ≥ 2. 5077a1 is the standard illustration: its L-function vanishes to
order exactly 3 — verified here — and no method extracts those three
generators from that vanishing. The points in this repository were found by
searching for them.

That gap is the Millennium Problem.

---

## Session 2 additions

### `complete2.py` — complete 2-descent (working)

For curves with full rational 2-torsion the étale algebra splits as Q×Q×Q,
the S-unit group collapses to squarefree integers supported on S, and the
candidate set is finite and exact. **No search box, no reduction theory.**
Sweep of `y² = x³ − n²x`, n = 1…60: **57/60 exact ranks**, every Selmer count
a power of 2, zero bug indicators.

Three bugs were found and fixed by the self-audits:
- p-adic branches where one square condition had already failed were being
  refined instead of pruned, manufacturing false "inconclusive" results.
- Inconclusive classes were being dropped from Sel₂ rather than counted into
  it. That shrinks the rank *upper* bound — the unsafe direction — and
  produced `rank ≤ −1` on n = 59.
- The point search required gcd(u₁,w) = 1. A projective point scales so that
  (u₁,u₂,u₃,w) are coprime *as a quadruple*; u₁ and w may share a factor.
  This cost half the image on n = 37.

### `quadfield.py`, `ks2.py` — K(S,2) for imaginary quadratic K (working)

Exact computation of K(S,2), validated against dim = |S_K| + 1 + dim Cl_S[2]
on 8 fields covering ramified, split, inert, and class number 2. Every result
a group of 2-power order.

A real bug was caught here: the first `is_square_in_K` accepted −2 as a square
in Q(i), because it only checked that a consistent (trace, norm) pair existed
for the square root and never that the root lay in the field.

### `onetors.py` — 2-descent over Q×K (working)

For curves with exactly one rational 2-torsion point, K imaginary quadratic.
The étale algebra is Q×K; the norm relation forces the Q-component, so the
descent lives in K*/K*² and Sel₂ ≤ K(S,2), computed exactly by `ks2.py` — no
search box. Each class gives a conic plus a square condition; the conic is
decided exactly by Hilbert symbols (Hasse–Minkowski), then parametrised to a
quartic.

**Scan of y² = x(x²+bx+c) over 67 curves: every Selmer count a power of 2,
zero inconclusive local tests, zero contradictions, and agreement with the
2-isogeny descent on every single curve.**

### `qpsol.py` — p-adic solubility of quartics (the fix that mattered)

Deciding z² = g(u,v) over Q_p by refining (u,v) mod pᵏ branches **p² per
level**, which is hopeless at p = 11 or 23 — it times out and reports
"inconclusive", and inconclusive classes must be counted into Sel₂, which
destroys the group structure.

Reducing to one variable — (u,v) primitive means v is a unit or p|v and u is
— branches **p per level** (8 at p=2, where the square class of a unit is only
fixed mod 8), and dividing p² of content out at each step keeps coefficients
bounded. Cross-checked against the old test: 1532 agreements at odd p and 600
at p=2, zero disagreements. On the quartics the old test could not decide it
returns an answer in under a millisecond where the old one spent 1.5 s to give
up.

This also feeds `descent2.py`, whose three benchmarks still pass.
