# LocalBarycentric.jl

High-order **local barycentric Lagrange interpolation** of tabulated 1-D data:
precomputed sliding-window weights, any `Real` number type (`Float64`,
`Double64`, `BigFloat`, even `Rational`), a graded stencil taper at the domain
edges, and per-query work shared across many data columns.

## Installation

Until registration:

```julia
using Pkg
Pkg.add(url = "https://github.com/Canonical111/LocalBarycentric.jl")
```

## Quickstart

```@example quick
using LocalBarycentric

nodes  = 3.0 .+ (0:400) ./ 101            # sorted (here: equispaced)
values = @. exp(-nodes) * sin(3nodes)     # Vector, or an (n × m) Matrix
cache  = LocalBarycentricCache(nodes; order = 28)

interpolate_local(cache, values, 3.505)
```

The cache is built once per node grid; evaluation is `O(order)` per query
regardless of grid size, plus one small matrix–vector product when `values`
is a matrix — the window and coefficients are shared across all its columns.

## When to use it

- Your data is **given** on a fixed (typically equispaced) grid — you cannot
  resample at Chebyshev points.
- You need **degree ≳ 10** accuracy (piecewise-cubic packages top out around
  `h⁴`).
- You evaluate **many times** (hot loops) and/or **many columns** sharing one
  grid.
- You need the *same* code path at `Float64`, `Double64`, or `BigFloat`, with
  a genuine extended-precision error floor.

If instead you can sample a function adaptively, use a Chebyshev package; if
you need smoothing of noisy data, use a spline package. See
[Theory](theory.md) for how this scheme relates to the alternatives.

## API

```@docs
LocalBarycentric
LocalBarycentricCache
interpolate_local
barycentric_weights
stencil_window
```
