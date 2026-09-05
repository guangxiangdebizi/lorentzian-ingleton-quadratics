# Explicit uniform Ingleton bounds for Lorentzian quadratics

[![exact-certificates](https://github.com/guangxiangdebizi/lorentzian-ingleton-quadratics/actions/workflows/ci.yml/badge.svg)](https://github.com/guangxiangdebizi/lorentzian-ingleton-quadratics/actions/workflows/ci.yml)

**Authors:** Xingyu Chen, Yifei Sun, Xiaojing Zhu (alphabetical by family name).
All authors are affiliated with Shanghai University of Electric Power,
Shanghai, China. No equal-contribution or corresponding-author designation
is assigned.

This repository contains a self-contained research note and exact certificates
for the following result.

Let `f` be a nonzero homogeneous quadratic with nonnegative coefficients whose
Hessian has at most one positive eigenvalue. For a variable subset `S` and
`epsilon > 0`, set

```text
V_epsilon(S) = f(1_S + epsilon 1_E).
```

For arbitrary, possibly overlapping, subsets `A,B,C,D` of the variables, define

```text
R_f = V(AB)V(AC)V(AD)V(BC)V(BD)
      --------------------------------- .
      V(A)V(B)V(ABC)V(ABD)V(CD)
```

Here juxtaposition denotes union. The main theorem proves a universal bound,
independent of the number of variables, `f`, the four subsets and `epsilon`:

```text
R_f >= c_2^* >= 476656^(-5).
```

Together with an exact five-vector determinantal upper-bound family:

```text
476656^(-5) <= c_2^* <= (-107 + 51 sqrt(17))/128
4.0642e-29  <= c_2^* <= 0.806862397707...
```

The lower bound is deliberately crude and not asserted sharp. This applies
to all nonnegative Lorentzian quadratics, not only size-two PSD determinants.

## Why the proof is not just support Ingleton

The support of a Lorentzian quadratic has rank two, but coefficients and the
regularizer can degenerate at different rates. A support-only inequality does
not control such simultaneous limits. The explicit proof instead works at
the actual coefficient values:

```text
16 membership atoms
  -> at most 32 polarized clones, with exact signature preservation
  -> weighted squared distances on a sphere by row-sum normalization
  -> a bottleneck ultrametric with distance distortion at most N-1
  -> a representable rank-two coefficient template
  -> one finite-length module carrying all ten evaluations
  -> max-product Ingleton
  -> sum/max comparison: [binom(N,2) (N-1)^2]^(-5).
```

For rational coefficients, the input signature gate and the finite
coefficient comparator use exact rational arithmetic, including singular
matrices, loops and parallel classes. This is not a decision oracle for
arbitrary black-box real numbers. The paper retains the earlier
semialgebraic degeneration argument as an alternative proof of positivity;
the new finite-scale comparator replaces that argument for the explicit
bound, rather than changing the existing common-module mechanism.

## Repository map

- [`paper/main.tex`](paper/main.tex): theorem, complete proof, upper-bound
  construction, scope, and comparison with adjacent literature.
- [`paper/references.bib`](paper/references.bib): primary references.
- [`notes/novelty-audit.md`](notes/novelty-audit.md): dated search record and
  claim boundary. It records absence of a match in the searched corpus; it is
  not presented as a proof of priority.
- [`src/reproduce.py`](src/reproduce.py): exact arithmetic reconstruction of
  the finite witness, limiting formula, common-principal-minor separation, and
  finite valuated-rank-two audit, plus the two new checker reports.
- [`src/explicit_floor.py`](src/explicit_floor.py): exact polarization,
  spherical/ultrametric comparison and multiscale regression.
- [`src/hostile_audit.py`](src/hostile_audit.py): independent rational input
  enumeration, singular PSD gate, and directly measured Laurent quotients.
- [`reports/certificates.json`](reports/certificates.json): generated manifest.
- [`tests/test_reproduce.py`](tests/test_reproduce.py): regression gates.
- [`tests/test_explicit_floor.py`](tests/test_explicit_floor.py): explicit
  bound, provenance, singular support and negative-control gates.
- [`REVISION_2026-09-05.md`](REVISION_2026-09-05.md): revision scope and
  locally reproduced checks.

## Reproduce

Python 3.11 or newer is sufficient; the scripts use only the standard library.
The paper job is pinned to TeX Live 2025; the release audit also compiled the
source with Tectonic 0.17.0.

```bash
python src/reproduce.py --check reports/certificates.json
python -m unittest discover -s tests -v
```

To regenerate the report:

```bash
python src/reproduce.py --output reports/certificates.json
```

The finite enumeration is a direction and serialization regression. It is not
used as a substitute for the general proof.

## Status and claim discipline

- **Proved in the note:** the explicit universal floor `476656^(-5)`.
- **Proved in the note:** the displayed exact upper bound.
- **Proved in the note:** a two-variable exact obstruction to reducing all ten
  evaluations to principal minors of one common `4 x 4` PSD matrix through a
  common scalar and diagonal modular gauge.
- **Not proved:** the value of the optimal constant.
- **Not claimed:** a degree-three or degree-four analogue.
- **Literature status:** targeted searches through 2026-08-27 found adjacent
  results for linear ranks and principal minors, but no statement matching the
  theorem above. See the novelty audit for queries and distinctions.
- **Revision literature status:** the explicit 2026-09-05 argument was
  developed and audited offline. The old dated search does not establish
  priority for this revision, and no new literature-priority claim is made.

## License

MIT. See [`LICENSE`](LICENSE).
