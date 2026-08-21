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
| `bsd.py` | Runs the two engines together and cross-checks them | 60-curve sweep, **zero contradictions**, 55/60 exact ranks |
| `lfun.py` | The L-function from point counts + modular continuation | a_p matches newform 37.2.a.a; analytic ranks 1, 2, 3 recovered |
| `sha.py`, `bm.py`, `bm3.py` | Ш made explicit; the Brauer–Manin obstruction computed as a number | Σ_v inv_v(A) = 1/2 on five independent adelic points |

### Headline results

- **55/60 exact ranks** on `y² = x³ − n²x` for n = 1…60.
- **Zero contradictions** between two independently written engines across those 60 curves — the point-search lower bound never once exceeded the descent upper bound.
- **Two confirmed Ш obstructions**, at n = 17 and n = 42: descent gives rank ∈ [0,2]; the L-function gives L(E,1) = 2.54 and 3.24, so analytic rank 0; Coates–Wiles (CM curves) then forces rank = 0, so the excess of 2 is entirely Ш. Consequently **17 and 42 are not congruent numbers**.
- **Analytic rank = algebraic rank** verified independently on 37a1 (1), 389a1 (2), 5077a1 (3) — computed by two programs sharing no code and no concepts.

Of ten apparent Ш gaps in the first sweep, five turned out to be a weak point
search and closed under a larger height bound. Distinguishing a real
obstruction from a tired search is most of the work.

---

## What does not work

Stated plainly so nobody builds on sand.

1. **`descent2.py`'s upper bound is heuristic, not proof.** The quartic
   enumeration uses a search box with no proven reduction bound behind it.
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

3. Local solvability tests can return inconclusive at their depth cap.
4. Canonical heights are numerical, with a heuristic threshold (1e−3 on the
   normalised Gram determinant) for independence.
5. `unfinished/enum2.py` is fast but **incomplete** — it misses quartics. Do not
   use it for upper bounds.

---

## Usage

```bash
python rank.py 0 0 0 -289 0      # 2-isogeny descent, Weierstrass a1 a2 a3 a4 a6
python descent2.py               # general 2-descent on the benchmark curves
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
