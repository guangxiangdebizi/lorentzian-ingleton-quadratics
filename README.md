# Explicit uniform Ingleton bounds for Lorentzian quadratics

[![exact-certificates](https://github.com/guangxiangdebizi/lorentzian-ingleton-quadratics/actions/workflows/ci.yml/badge.svg)](https://github.com/guangxiangdebizi/lorentzian-ingleton-quadratics/actions/workflows/ci.yml)

**Authors:** Xingyu Chen, Yifei Sun, Xiaojing Zhu (alphabetical by family name).
All authors are affiliated with Shanghai University of Electric Power,
Shanghai, China. No equal-contribution or corresponding-author designation
is assigned.

**Current explicit bound (2026-09-06):** `R >= 1/16384`, or `R >= 1/16`
for disjoint queries. The proof combines a finite Schur/Fischer bound for
PSD pencils with real Clifford determinant powers. It has two internal
independent hostile audits; global priority and sharpness are not certified.

**Earlier MST attribution correction:** the qualitative floor follows from
Huang--Huh--Soskin--Wang's coefficient-to-tree comparison combined with
the shared-module transfer. This note supplies an explicit self-contained
construction; no improvement over prior quantitative tree approximation
or independent new metric mechanism is claimed. See the final section
of the paper and [the updated overlap audit](notes/novelty-audit.md).

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
R_f >= c_2^* >= 1/16384.
```

Together with an exact five-vector determinantal upper-bound family:

```text
1/16384         <= c_2^* <= (-107 + 51 sqrt(17))/128
0.00006103515625 <= c_2^* <= 0.806862397707...
```

The lower bound improves the previous local `476656^(-5)` bound by more
than `10^24`, but remains deliberately crude and is not asserted sharp. It applies
to all nonnegative Lorentzian quadratics, not only size-two PSD determinants.

## Why the proof is not just support Ingleton

The support of a Lorentzian quadratic has rank two, but coefficients and the
regularizer can degenerate at different rates. A support-only inequality does
not control such simultaneous limits. The current proof works at actual
coefficient values:

```text
Schur complements + an auxiliary covariance + conditional Fischer
  -> a four-block principal-determinant bound 4^(-dim A)
  -> numerator-only overlap duplication
  -> PSD bound 128^(-d) in every dimension
  -> a real Clifford representation det(sum z_i P_i)=(2q)^(M/2)
  -> exponent cancellation: R_q >= 128^(-2)=1/16384.
```

The Clifford matrices may have exponential size; this is an exact proof
device, not an efficient representation algorithm. The earlier rational
MST comparator and semialgebraic degeneration argument remain as supporting
alternatives. Their valid but weaker `476656^(-5)` certificate is preserved
under explicitly named legacy fields. No exact oracle for unspecified
black-box real numbers is asserted.

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
- [`src/gaussian_floor.py`](src/gaussian_floor.py): exact auxiliary-covariance
  identities and a nonunit determinant-ratio fixture.
- [`src/clifford_floor.py`](src/clifford_floor.py): signed exterior generators,
  complete low-rank symbolic determinant identities, and degenerate cases.
- [`tests/test_gaussian_clifford.py`](tests/test_gaussian_clifford.py): new
  constant, provenance, rank-zero, power-cancellation, and negative-control gates.
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
