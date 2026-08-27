# Uniform Ingleton bounds for Lorentzian quadratics

[![exact-certificates](https://github.com/guangxiangdebizi/lorentzian-ingleton-quadratics/actions/workflows/ci.yml/badge.svg)](https://github.com/guangxiangdebizi/lorentzian-ingleton-quadratics/actions/workflows/ci.yml)

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

Here juxtaposition denotes union. The main theorem proves that there is a
optimal universal constant `c_2^* > 0`, independent of the number of variables, `f`, the
four subsets, and `epsilon`, such that

```text
R_f >= c_2^*.
```

The proof does not produce the optimal constant. An exact five-vector
determinantal family proves

```text
0 < c_2^* <= (-107 + 51 sqrt(17))/128
          = 0.806862397707...
```

## Why the proof is not just support Ingleton

The support of a Lorentzian quadratic has rank two, but coefficients and the
regularizer can degenerate at different rates. A support-only inequality does
not control such simultaneous limits. The proof instead uses this chain:

```text
16 membership atoms
  -> semialgebraic curve selection
  -> a rank-two valuated matroid
  -> exact realization over a discretely valued field
  -> one finite-length module carrying all ten evaluations
  -> modular-length Ingleton
  -> no degeneration can have positive ratio valuation.
```

## Repository map

- [`paper/main.tex`](paper/main.tex): theorem, complete proof, upper-bound
  construction, scope, and comparison with adjacent literature.
- [`paper/references.bib`](paper/references.bib): primary references.
- [`notes/novelty-audit.md`](notes/novelty-audit.md): dated search record and
  claim boundary. It records absence of a match in the searched corpus; it is
  not presented as a proof of priority.
- [`src/reproduce.py`](src/reproduce.py): exact arithmetic reconstruction of
  the finite witness, limiting formula, common-principal-minor separation, and
  finite valuated-rank-two audit.
- [`reports/certificates.json`](reports/certificates.json): generated manifest.
- [`tests/test_reproduce.py`](tests/test_reproduce.py): regression gates.

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

- **Proved in the note:** existence of a universal positive constant.
- **Proved in the note:** the displayed exact upper bound.
- **Proved in the note:** a two-variable exact obstruction to reducing all ten
  evaluations to principal minors of one common `4 x 4` PSD matrix through a
  common scalar and diagonal modular gauge.
- **Not proved:** the value of the optimal constant.
- **Not claimed:** a degree-three or degree-four analogue.
- **Literature status:** targeted searches through 2026-08-27 found adjacent
  results for linear ranks and principal minors, but no statement matching the
  theorem above. See the novelty audit for queries and distinctions.
- **Priority wording:** the result is released as *apparently new after a dated
  search*, not as a proof that no unpublished or differently formulated result
  exists.

## License

MIT. See [`LICENSE`](LICENSE).
