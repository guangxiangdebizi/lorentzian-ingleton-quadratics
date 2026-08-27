# Novelty and priority audit

Audit date: **2026-08-27**

## Claim checked

The searched claim was the following precise statement, not the broad phrase
"an Ingleton inequality for Lorentzian objects":

> There is a universal positive lower bound for the multiplicative Ingleton
> ratio of the ten regularized set evaluations
> `f(1_S + epsilon 1_E)`, uniformly over all nonnegative Lorentzian
> homogeneous quadratics, all numbers of variables, arbitrary overlapping
> `A,B,C,D`, and all `epsilon > 0`.

The exact five-vector upper-bound family was searched separately.

## Search methods

The audit used targeted web, arXiv full-text/source, and Crossref searches. The
following queries were run verbatim or with harmless punctuation variants:

```text
"Lorentzian quadratic" "Ingleton"
"Lorentzian polynomials" "Ingleton inequality"
"volume polynomial" "Ingleton inequality"
"valuated matroid" "Ingleton inequality"
"quantitative Ingleton" polynomial
"regularized Ingleton ratio" polynomial
"f(1_S" Ingleton polynomial
"module length" Ingleton inequality
```

The arXiv source of the two closest recent papers was downloaded and searched
at statement level, not only by title or abstract:

- Al Ahmadieh--Rincon--Vinzant--Yu,
  [*Tropicalizing Principal Minors of Positive Definite Matrices*](https://arxiv.org/abs/2410.11220),
  version 2, 2025.
- Boege--Bouthat,
  [*Sharp Inequalities for Products of Principal Minors of Positive Definite Matrices*](https://arxiv.org/abs/2606.23632),
  2026.
- Sendov,
  [*A Bounded Determinantal Ratio for Positive Definite Matrices*](https://arxiv.org/abs/2607.27216),
  2026.

Crossref returned information-theoretic stability papers and unrelated uses of
"Lorentzian quadratic," but no record matching the claim above.

## Closest primary literature

| Source | What it establishes | Why it does not match this note |
|---|---|---|
| Branden--Huh, [*Lorentzian Polynomials*](https://arxiv.org/abs/1902.03719) (2020) | Lorentzian/M-convex theory, closure and polarization, tropicalized Lorentzian polynomials | It supplies essential tools, but does not state a ten-evaluation Ingleton lower bound. |
| Ingleton, *Representation of Matroids* (1971) | Additive Ingleton for ranks of linear subspace arrangements | It controls exact rank functions, not simultaneous coefficient and regularizer degenerations of a polynomial. |
| Speyer--Sturmfels, [*The Tropical Grassmannian*](https://arxiv.org/abs/math/0304218) (2004) | Rank-two tropical Grassmannians and tree metrics | It supplies the rank-two realization step, not the compactness theorem or ratio bound. |
| Al Ahmadieh--Rincon--Vinzant--Yu (2025) | Tropicalized principal minors and local coefficient inequalities valid for Lorentzian polynomials | The paper treats principal-minor tropicalization and lifts of local `M^natural` constraints; its statements do not contain the regularized ten-evaluation theorem. |
| Hall--Johnson, [*Bounded Ratios of Products of Principal Minors of Positive Definite Matrices*](https://arxiv.org/abs/0806.2645) (2008) | Bounded principal-minor ratios | The domain is principal-minor vectors of one positive-definite matrix, not arbitrary Lorentzian quadratics evaluated at set-dependent points. |
| Boege--Bouthat (2026) | Sharp infimum `16/27` for the Ingleton ratio of `4 x 4` positive-definite principal minors | This is the closest ratio result, but its ten inputs are principal minors of one fixed matrix. The current note allows arbitrary variable count and overlap and starts from the full nonnegative Lorentzian quadratic cone. Neither theorem is asserted to imply the other. |
| Sendov (2026) | Independent proof of the reciprocal sharp supremum `27/16` for the same `4 x 4` principal-minor ratio | It has the same principal-minor domain boundary as the preceding result and does not state a Lorentzian-polynomial evaluation theorem. |
| Matveev--Romashchenko (2025/2026) and Csirmaz (2026) | Stability of conditional Ingleton implications for entropy profiles | These concern conditional mutual information and almost-entropic functions, not Lorentzian-polynomial evaluations. |

## Reduction audit for the closest 2026 result

After relabelling the four sets, the ten-subset pattern in the 2026
principal-minor theorem is the same combinatorial Ingleton pattern used here.
This makes it the main collision risk.  To imply the present theorem, however,
one would need to encode all ten regularized polynomial evaluations as scalar
principal minors of one common `4 x 4` positive-definite matrix, up to a
modular rescaling whose factors cancel from the ratio.  No such simultaneous
encoding exists in general.  The note proves this with the exact input
`f=x0*x1`, `epsilon=1/3`, `A={0}`, `B={1}`, `C=empty`, and `D={0,1}`.
After canonical modular-gauge normalization its five free coordinates are
`(1,1,4,4,16)`.  Schur complements force the last coordinate to be at most
`13`, an exact gap of `3`.  The certificate stream is reproduced in
`reports/certificates.json`.

Several standard reductions were checked.  A degree-two determinant of a
`2 x 2` real symmetric pencil has bounded Hessian rank, whereas a nonnegative
Lorentzian quadratic such as `e_2(x_1,...,x_n)` can have arbitrarily large
Hessian rank.  Gram constructions naturally produce principal minors of a
matrix indexed by variables rather than a common scalar `4 x 4` encoding of
the ten values.  PSD-pencil and block constructions produce additive
determinant components or repeated monomials rather than the required common
principal-minor vector.  The exact separation rules out the natural common
scalar plus diagonal modular-gauge route, including singular PSD matrices. It
is not a proof that no more indirect use of a principal-minor inequality
exists.

## Search conclusion

No source in the searched corpus states the theorem of this repository or the
exact upper-bound family. The result is therefore being released as an
**apparently new, independently checkable research note**.

This wording is deliberate. Database and full-text searches cannot prove that
no one has ever obtained the result, and unpublished notes or differently
phrased statements may exist. The repository does **not** use "first" or
"never previously discovered" as a theorem. A prior-art issue with a precise
citation is sufficient reason to revise the claim and attribution.

## Claim boundary

The audit does not support any of the following stronger statements:

- the displayed upper bound is optimal;
- the same positive bound exists in degree three or four;
- every Lorentzian quadratic is a principal-minor vector;
- the finite enumeration proves the general theorem;
- absence from the searched databases proves priority.
