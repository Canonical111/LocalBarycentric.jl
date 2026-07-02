# LocalBarycentric.jl

High-order **local barycentric Lagrange interpolation** of tabulated 1-D data —
precomputed sliding-window weights, any `Real` number type (Float64,
[Double64](https://github.com/JuliaMath/DoubleFloats.jl), BigFloat, …), a graded
stencil taper at the domain edges, and per-query work shared across many data
columns.

```julia
using LocalBarycentric

nodes  = 3.0 .+ (0:400) ./ 101                  # sorted, e.g. equispaced
values = @. exp(-nodes) * sin(3nodes)           # Vector, or (n × m) Matrix
cache  = LocalBarycentricCache(nodes; order = 28)

interpolate_local(cache, values, 3.505)         # ≈ exp(-3.505)sin(10.515), ~1e-16 rel
cache(values, 3.505)                            # same, callable sugar
```

## Why this package

No registered Julia package (as of mid-2026) combines all five of:

1. **High polynomial order** on tabulated data — `order` is unlimited
   (Interpolations.jl caps at cubic; Dierckx at quintic).
2. **Local windowed evaluation** — one `order+1`-node stencil per query, O(1)
   in the grid size; no global solve, no Runge blow-up on equispaced grids
   (BarycentricInterpolation.jl is global-only; BSplineKit solves a global
   collocation system).
3. **Genuine arbitrary precision** — the same code path runs Float64,
   Double64, BigFloat: a BigFloat cache reaches ~1e-50 interpolation error
   where several packages that *accept* BigFloat silently round through
   Float64 internally.
4. **Precomputed weights** — one weight row per window position, built once;
   hot-loop queries do O(order) work plus one small mat-vec.
5. **Multi-column sharing** — evaluating an `n × m` table at one point locates
   the window and computes barycentric coefficients **once** for all `m`
   columns. For wide tables this is the dominant cost factor (measured ~340×
   vs. per-column interpolant objects at m = 338).

Benchmarks against the alternatives (identical battery: 401–40 000 equispaced
nodes, degree 28, errors vs a 256-bit reference; full details in the FunBootV2
provenance below): per-point evaluation of a 338-column table costs **3–6 µs
(Float64) / ~60 µs (Double64)** here (Ryzen 32-core workstation), vs ~0.1 ms / 3.3 ms for per-column
B-spline objects (BSplineKit) and ~3.7 ms Float64 for as-shipped global
Floater–Hormann (BaryRational) at n = 15 000. At BigFloat the gap widens to
hours-vs-seconds per 10⁴-point workload.

## The edge taper

Near a domain boundary a full-width window cannot center on the query; the
query sits at the window's edge, where the equispaced degree-`n` Lebesgue
function is ~2ⁿ, amplifying *data* rounding (≈10⁶× at degree 28 — i.e. Float64
inputs lose 5–6 digits precisely at the boundary). `stencil_window` therefore
grades the width down (`clamp(2k, edge_order_min + 1, order + 1)` nodes, `k` =
nodes on the thinner side), keeping the query near the window center. Measured
on degree 28, step 1/101: worst-case boundary error improves from ~5e-11 to
~4e-15 with Float64 data. Interior queries are untouched. The floor
`edge_order_min = 12` keeps truncation below the Float64 data floor for grid
steps ≲ 0.05; on much coarser grids consider raising it.

Queries that hit a node exactly return the tabulated row unchanged. Queries
outside the node hull extrapolate with the boundary window — avoid relying on
that.

**Order/type envelope**: the unnormalized weights need exponent range; the
constructor throws if they overflow. Measured on step-1/101 equispaced grids:
Float64 is comfortable through degree 120 (max |w| ≈ 5e76) and degree 28 works
down to step ~1e-12; Float32 at degree 28 is marginal (~4 orders of headroom);
BigFloat/Double64 are effectively unconstrained. Interior accuracy stays at the
data floor at every tested degree — centered windows do not suffer Runge
oscillation.

## Conventions

- `order` = polynomial degree; a window uses `order + 1` nodes, matching
  Mathematica's `InterpolationOrder -> order` node count. (For even orders
  Mathematica aligns its window one node to the *right* of this package's
  choice — values agree to interpolation accuracy on smooth data but are not
  bit-interchangeable.)
- `weight_type` decouples weight precision from node storage:
  `LocalBarycentricCache(f64_nodes; weight_type = Double64)` computes weights
  from exactly-lifted nodes, giving a true extended-precision path over
  Float64-stored grids.
- Result type: `promote_type(weight_type, typeof(x))`; narrower `values` are
  lifted exactly per window.
- Second (true) barycentric form; numerically stable arbitrarily close to
  nodes (tested at 1e-15 relative separation). Unnormalized weights: safe for
  equispaced grids down to step ~1e-12 at degree 28 (Float64 range).

## Provenance

Extracted from the conformal-bootstrap package **FunBootV2** (Zechuan Zheng),
where it interpolates 2 GB tables of crossing-equation functionals inside a
custom simplex solver at Float64/Double64 precision, and where every numeric
claim above was measured. The extraction is bit-for-bit faithful to the
original implementation (verified by a golden-run regression producing
identical results down to the solver's full pivot trajectory).

## License

MIT.
