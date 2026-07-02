"""
    LocalBarycentric

High-order **local barycentric Lagrange interpolation** of tabulated 1-D data,
with precomputed sliding-window weights, arbitrary-precision number types, a
graded stencil taper at the domain edges, and work-sharing across many data
columns per query point.

The scheme matches Mathematica's `Interpolation[..., InterpolationOrder -> n]`
node-count convention: polynomial degree `n` uses a local window of `n + 1`
nodes centered on the query (Mathematica's window alignment differs by one node
for even orders; see the README).

Designed for the workload where no existing Julia package fits: degrees ≳ 10 on
equispaced tabulated grids, evaluated ~10⁴–10⁶ times in a hot loop, in Float64,
Double64, BigFloat, or any `Real` subtype — with the same code path.

```julia
nodes  = 3.0 .+ (0:400) ./ 101
values = @. exp(-nodes) * sin(3nodes)          # or an (n × m) matrix of columns
cache  = LocalBarycentricCache(nodes; order = 28)
interpolate_local(cache, values, 3.505)
```
"""
module LocalBarycentric

using LinearAlgebra

export LocalBarycentricCache, barycentric_weights, stencil_window, interpolate_local

"""
    barycentric_weights(xs) -> Vector

Unnormalized barycentric weights `wᵢ = 1 / ∏_{j≠i} (xᵢ - xⱼ)` for the node
vector `xs` (distinct values). O(n²); element type follows `eltype(xs)`.
Pure and deterministic: repeated calls on the same input are bit-identical.
"""
function barycentric_weights(xs::AbstractVector{T}) where T <: Real
    n = length(xs)
    weights = Vector{T}(undef, n)
    for i in 1:n
        prod = one(T)
        xi = xs[i]
        for j in 1:n
            j == i && continue
            prod *= (xi - xs[j])
        end
        weights[i] = inv(prod)
    end
    return weights
end

"""
    stencil_window(nodes, x, stencil; edge_min = stencil) -> (start, width)

Select the local window for a query `x` over sorted `nodes`. Interior queries
get the full centered `stencil`. Near a domain edge a full-width window would
clamp and evaluate `x` at the window edge, where the degree-(stencil-1)
equispaced Lebesgue function is ~2^stencil and amplifies data rounding (~10⁶ at
stencil 29). The width is instead tapered to `clamp(2k, min(edge_min, stencil),
stencil)`, `k` being the number of nodes on the thinner side of `x`, keeping
the query near the window center.

`nodes` and `x` may have different `Real` types; comparisons follow Julia's
exact mixed-type promotion (e.g. `Float64` nodes vs `Double64`/`BigFloat`
query compare exactly).
"""
function stencil_window(xnodes::AbstractVector{<:Real}, x::Real, stencil::Int;
                        edge_min::Int = stencil)
    n = length(xnodes)
    j = searchsortedfirst(xnodes, x)     # x lies between nodes j-1 and j
    k = min(j - 1, n - j + 1)
    s = clamp(2k, min(edge_min, stencil), stencil)
    lo = max(1, j - cld(s, 2))
    hi = min(n, lo + s - 1)
    return max(1, hi - s + 1), s
end

"""
    LocalBarycentricCache(nodes; order = 10, weight_type = eltype(nodes),
                          edge_order_min = 12)

Precompute barycentric weights for every sliding window of `order + 1` nodes
over the sorted vector `nodes` (Mathematica `InterpolationOrder` convention:
`order` is the polynomial degree). One weight row is stored per window start,
so hot-loop evaluation never recomputes weights except in the tapered edge
zone.

`weight_type` sets the arithmetic type of the stored weights independently of
the node storage type: e.g. `Float64` nodes with `weight_type = Double64`
gives an extended-precision evaluation path whose weights are computed from
exactly-lifted nodes (lifting `Float64 → Double64/BigFloat` is exact, so this
equals computing in the wide type from the start).

`edge_order_min` floors the tapered edge windows at `edge_order_min + 1` nodes
(default degree 12 — at that size the edge Lebesgue constant is ~10², small
enough that data rounding stays near the working-precision floor while
truncation remains negligible on fine grids).
"""
struct LocalBarycentricCache{TN <: Real, TW <: Real}
    nodes::Vector{TN}
    stencil::Int
    edge_min::Int
    weights::Matrix{TW}
end

function LocalBarycentricCache(nodes::AbstractVector{TN};
                               order::Int = 10,
                               weight_type::Type{TW} = TN,
                               edge_order_min::Int = 12) where {TN <: Real, TW <: Real}
    issorted(nodes) || throw(ArgumentError("nodes must be sorted"))
    allunique(nodes) || throw(ArgumentError("nodes must be distinct"))
    stencil = clamp(order + 1, 2, length(nodes))
    edge_min = min(edge_order_min + 1, stencil)
    # Weights from the wide-type-lifted full node vector; lifting is elementwise,
    # so windowing before or after the lift gives identical bits.
    wnodes = TW === TN ? nodes : TW.(nodes)
    starts = max(length(nodes) - stencil + 1, 1)
    weights = Matrix{TW}(undef, starts, stencil)
    for start in 1:starts
        weights[start, :] .= barycentric_weights(@view wnodes[start:(start + stencil - 1)])
    end
    all(isfinite, weights) || throw(ArgumentError(
        "barycentric weights overflow $(TW)'s exponent range at order $(stencil - 1) " *
        "on this grid; reduce the order, coarsen the grid, or use a wider weight_type"))
    return LocalBarycentricCache{TN, TW}(collect(nodes), stencil, edge_min, weights)
end

# Shared core: window, coefficients, denominator. Returns (rng, coeff, den) with
# coeff/den in TC = promote_type(TW, typeof(x)). Sequential accumulation of den
# is deliberate (deterministic, matches the weight-cache convention).
@inline function _window_coeffs(cache::LocalBarycentricCache{TN, TW}, x::Real) where {TN, TW}
    TC = promote_type(TW, typeof(x))
    nodes = cache.nodes
    start, s = stencil_window(nodes, x, cache.stencil; edge_min = cache.edge_min)
    rng = start:(start + s - 1)
    xs = TC === TN ? (@view nodes[rng]) : TC.(@view nodes[rng])
    weights = s == cache.stencil ? (@view cache.weights[start, :]) :
                                   barycentric_weights(xs)
    coeff = Vector{TC}(undef, s)
    den = zero(TC)
    for i in 1:s
        coeff[i] = weights[i] / (x - xs[i])
        den += coeff[i]
    end
    return rng, coeff, den
end

@inline function _exact_node_index(nodes::AbstractVector, x::Real)
    j = searchsortedfirst(nodes, x)
    return (j <= length(nodes) && nodes[j] == x) ? j : 0
end

"""
    interpolate_local(cache, values, x)

Evaluate the local barycentric interpolant of `values` at `x`.

`values` is either a length-`n` vector (one data column; returns a scalar) or
an `n × m` matrix (`m` columns sharing the same nodes; returns a length-`m`
vector). The window location, weights, and barycentric coefficients are
computed **once per query** and shared across all `m` columns — for many-column
tables this is the dominant cost advantage over per-column interpolant objects.

The result type is `promote_type(weight_type, typeof(x))`; `values` of a
narrower type are lifted exactly per window, and values of a **wider** type are
demoted to it — choose `weight_type` (or the query type) at least as wide as
your data to avoid silent precision loss. Queries that hit a node exactly
return the tabulated row unchanged (no interpolation error, no division).
Queries outside the node hull evaluate the boundary window's polynomial
(extrapolation — accuracy degrades rapidly; avoid by construction).
"""
function interpolate_local(cache::LocalBarycentricCache{TN, TW},
                           values::AbstractMatrix, x::Real) where {TN, TW}
    size(values, 1) == length(cache.nodes) ||
        throw(DimensionMismatch("values has $(size(values, 1)) rows for $(length(cache.nodes)) nodes"))
    TC = promote_type(TW, typeof(x))
    j = _exact_node_index(cache.nodes, x)
    j != 0 && return TC.(@view values[j, :])

    rng, coeff, den = _window_coeffs(cache, x)
    block = TC === eltype(values) ? (@view values[rng, :]) : TC.(@view values[rng, :])
    return transpose(block) * coeff / den
end

function interpolate_local(cache::LocalBarycentricCache{TN, TW},
                           values::AbstractVector, x::Real) where {TN, TW}
    length(values) == length(cache.nodes) ||
        throw(DimensionMismatch("values has $(length(values)) entries for $(length(cache.nodes)) nodes"))
    TC = promote_type(TW, typeof(x))
    j = _exact_node_index(cache.nodes, x)
    j != 0 && return TC(values[j])

    rng, coeff, den = _window_coeffs(cache, x)
    num = zero(TC)
    @inbounds for (i, r) in enumerate(rng)
        num += TC(values[r]) * coeff[i]
    end
    return num / den
end

(cache::LocalBarycentricCache)(values, x::Real) = interpolate_local(cache, values, x)

end # module
