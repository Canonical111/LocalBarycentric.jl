# Examples

## Basic interpolation

```@example basic
using LocalBarycentric

h = 1 / 101
nodes  = collect(3.0 .+ h .* (0:400))
f(x)   = exp(-x) * sin(3x)
values = f.(nodes)

cache = LocalBarycentricCache(nodes; order = 28)
x = 3.505
interpolate_local(cache, values, x), f(x)
```

## Many columns, one query

An `n × m` value matrix shares the window search and barycentric coefficients
across all `m` columns — build one cache per grid, not per column:

```@example basic
values2 = hcat(values, cos.(nodes), values .^ 2)   # 401 × 3
interpolate_local(cache, values2, x)
```

## Extended precision: Double64 weights over Float64 data

`weight_type` decouples the arithmetic from the node storage. This is the
pattern for a high-precision evaluation path over an existing Float64 table:

```@example basic
using DoubleFloats

cache_hp = LocalBarycentricCache(nodes; order = 28, weight_type = Double64)
xd = Double64(nodes[200]) + Double64(1) / 303      # query below Float64 resolution
r  = interpolate_local(cache_hp, values, xd)
typeof(r), r
```

The result responds to query changes far below `eps(Float64)` — the window
search compares `Float64` nodes against the `Double64` query exactly, and only
the stencil window of `values` is lifted (exactly) per call.

## Fully BigFloat

With wide-typed data *and* weights the error floor follows the working
precision (measured ``\sim 10^{-50}`` at degree 28, 256 bits):

```@example basic
setprecision(BigFloat, 256)
bnodes  = BigFloat.(nodes)
bvalues = f.(bnodes)
bcache  = LocalBarycentricCache(bnodes; order = 28)
bx = (bnodes[200] + bnodes[201]) / 2
abs(interpolate_local(bcache, bvalues, bx) - f(bx))
```

## Edge behavior

Queries near the domain boundary automatically use narrower centered windows
(see [Theory](theory.md)); nothing is required from the caller:

```@example basic
using LocalBarycentric: stencil_window
[stencil_window(nodes, nodes[1] + k * h + h / 2, 29; edge_min = 13) for k in (0, 3, 7, 20)]
```

## Exact node hits and the overflow guard

```@example basic
interpolate_local(cache, values, nodes[57]) === values[57]
```

```@example guard
using LocalBarycentric
tiny = Float32.(1.0 .+ (0:40) .* 1.0f-6)   # degree 28 weights overflow Float32
try
    LocalBarycentricCache(tiny; order = 28)
catch e
    e
end
```

## Practical guidance

- **Reuse the cache.** Building is `O(n · order²)`; do it once per grid.
- **Batch columns** into a matrix rather than looping `interpolate_local`
  per column.
- **Choose `order`** by your accuracy target and grid spacing; on fine grids
  the interior is data-floor-limited already at moderate order, and the
  default edge floor (degree 12) is usually right. Raise `edge_order_min`
  only on coarse grids where degree-12 truncation is visible.
- **Avoid extrapolation** — outside the node hull the boundary window's
  polynomial is evaluated, and its accuracy degrades rapidly.
