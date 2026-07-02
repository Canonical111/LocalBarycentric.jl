# Theory

## The second barycentric form

Given ``s`` distinct nodes ``x_1 < \dots < x_s`` and values ``f_i``, the unique
degree-``(s-1)`` interpolating polynomial can be evaluated as

```math
p(x) \;=\; \frac{\displaystyle\sum_{i=1}^{s} \frac{w_i}{x - x_i}\, f_i}
                {\displaystyle\sum_{i=1}^{s} \frac{w_i}{x - x_i}},
\qquad
w_i \;=\; \prod_{j \neq i} \frac{1}{x_i - x_j}.
```

This *second (true) barycentric form* is backward-stable and — unlike the
Newton or Vandermonde routes — evaluates in ``O(s)`` once the weights are
known. It is well behaved arbitrarily close to nodes (tested here at
``10^{-15}`` relative separation), and an exact node hit is returned directly.
The weights are computed once per window in ``O(s^2)`` and cached
(`barycentric_weights` is pure, so cached and recomputed weights are
bit-identical).

## Local windows, not global polynomials

A *global* polynomial through many equispaced nodes is destroyed by the Runge
phenomenon. This package instead slides a **local window of `order + 1`
nodes** across the grid — the same convention as Mathematica's
`Interpolation[..., InterpolationOrder -> order]` — and evaluates the window's
polynomial. For a query near the *center* of its window, the equispaced
Lebesgue constant (the factor by which input errors are amplified) stays
``O(1)``: measured **1.9** at degree 28. High order is therefore free in the
interior: interior accuracy sits at the *data* rounding floor at every tested
degree up to 120.

One weight row is precomputed per window position, so hot-loop evaluation
never recomputes weights (except in the tapered edge zone, where windows are
narrower than the cached full width and the ``O(s^2)`` recompute is cheaper
than the interior mat-vec anyway).

## The edge problem and the graded taper

Near a domain boundary the window cannot center; a full-width window *clamps*
and the query sits at the window's **edge**, where the equispaced Lebesgue
function grows like ``2^s``. Measured at degree 28, grid step ``1/101``,
Float64 data:

| query position | Lebesgue | relative error | with exact data |
|---|---|---|---|
| subinterval 1 (at the boundary) | ``1.1 \times 10^6`` | ``3.9 \times 10^{-13}`` | ``3 \times 10^{-27}`` |
| subinterval 5 | ``216`` | ``2.0 \times 10^{-15}`` | ``9 \times 10^{-31}`` |
| interior (centered) | ``1.9`` | ``2.0 \times 10^{-16}`` | ``4 \times 10^{-32}`` |

The last column shows the scheme itself is nearly exact even at the edge — the
error is amplified **input rounding**, so no amount of higher-precision
*arithmetic* recovers it. The cure is geometric: [`stencil_window`](@ref)
tapers the width to

```math
s(k) \;=\; \operatorname{clamp}\!\bigl(2k,\; s_{\min},\; \text{order}+1\bigr),
```

where ``k`` counts nodes on the thinner side of the query, so the query stays
near the window center. The floor ``s_{\min} = \texttt{edge\_order\_min} + 1``
(default degree 12, i.e. 13 nodes) balances the two error terms: amplification
``\sim 2^{s}`` falls with ``s`` while truncation ``\sim (ch)^{s}`` rises, and
at degree 12 on fine grids both sit below the Float64 data floor. Measured
effect at degree 28: worst boundary error improves from ``5 \times 10^{-11}``
to ``4 \times 10^{-15}`` (~14 000×), with interior queries bit-identical to
the untapered scheme. **More points near a boundary make edge evaluation
worse, not better** — the taper deliberately uses fewer.

## Precision semantics

Three types cooperate:

- `TN` — node storage type;
- `TW = weight_type` — the type weights are computed and stored in;
- the query type; results are `TC = promote_type(TW, typeof(x))`.

The key fact making the mixed-precision path exact: **lifting a narrower float
to a wider type is exact** (`Float64 → Double64/BigFloat` preserves the value
bit-for-bit), and lifting commutes with windowing. So a cache with `Float64`
nodes and `weight_type = Double64` has weights *identical* to those computed
from wide-typed nodes from the start, and window searches compare mixed types
exactly. Values narrower than `TC` are lifted exactly per window; values
**wider** than `TC` are demoted — choose `weight_type` at least as wide as
your data.

What extended precision buys is *arithmetic* headroom, not data repair: with
Float64 input data the error floor is the data's ``\sim 10^{-16}`` rounding
(times the small centered Lebesgue constant); with BigFloat data and weights
the measured floor at degree 28 is ``\sim 10^{-50}``.

## Order and type envelope

The unnormalized weights need exponent range:
``|w| \sim h^{-(s-1)}`` up to factorial factors. Measured on step-``1/101``
grids: Float64 is comfortable through degree 120 (``\max|w| \approx 5 \times
10^{76}``), and degree 28 works down to step ``\sim 10^{-12}``; Float32 at
degree 28 is marginal; `BigFloat`/`Double64` are effectively unconstrained.
The constructor **throws** if weights overflow, rather than propagating `Inf`.

## Complexity

| operation | cost |
|---|---|
| cache build | ``O(n \cdot s^2)`` — independent of the number of value columns |
| scalar query | ``O(s)`` after ``O(\log n)`` window search |
| `m`-column query | window + coefficients **once**, then one ``s \times m`` mat-vec |

The multi-column sharing is the structural advantage over one-interpolant-per-
column designs: at ``m = 338`` columns it was measured at ~340× less per-point
work than per-column spline objects at extended precision.
