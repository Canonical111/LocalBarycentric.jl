using Test
using LocalBarycentric
using DoubleFloats
using Random
using LinearAlgebra

setprecision(BigFloat, 256)

# Reference battery grid: equispaced, step 1/101 (typical tabulated-data spacing).
const H = 1 / 101
const N = 401
const NODES = [3.0 + k * H for k in 0:(N - 1)]
fexact(x) = exp(-x) * sin(3x) + exp(-x * log(oftype(x, 4.0)))
relerr(v, x) = abs(Float64((BigFloat(v) - fexact(BigFloat(x))) / fexact(BigFloat(x))))

@testset "LocalBarycentric" begin

@testset "weights: correctness and purity" begin
    xs = [0.0, 1.0, 3.0]
    w = barycentric_weights(xs)
    # w_i = 1/prod(x_i - x_j): [1/((0-1)(0-3)), 1/((1-0)(1-3)), 1/((3-0)(3-1))]
    @test w ≈ [1 / 3, -1 / 2, 1 / 6]
    # purity: bit-identical on repeated calls
    @test all(barycentric_weights(NODES[1:29]) .=== barycentric_weights(NODES[1:29]))
end

@testset "window selection invariants (sweep)" begin
    for n in (2, 3, 5, 13, 29, 30, 100), stencil in (2, min(5, n), min(13, n), min(29, n))
        nodes = collect(range(0.0, 1.0, length = n))
        h = n > 1 ? nodes[2] - nodes[1] : 1.0
        for x in vcat(nodes .+ 1e-3 * h, nodes .- 1e-3 * h,
                      [nodes[1] - h, nodes[end] + h],
                      [(nodes[i] + nodes[i + 1]) / 2 for i in 1:(n - 1)])
            start, s = stencil_window(nodes, x, stencil; edge_min = min(13, stencil))
            @test 1 <= start
            @test start + s - 1 <= n
            @test 2 <= s <= stencil || stencil < 2
        end
    end
    # interior windows are full-width and centered
    start, s = stencil_window(NODES, NODES[200] + 0.4H, 29; edge_min = 13)
    @test s == 29
    @test start == 201 - 15   # 15 below, 14 above searchsortedfirst index
    # left edge tapers, right edge mirrors
    _, sl = stencil_window(NODES, NODES[1] + 0.5H, 29; edge_min = 13)
    _, sr = stencil_window(NODES, NODES[end] - 0.5H, 29; edge_min = 13)
    @test sl == 13 && sr == 13
end

@testset "degree-d polynomial exactness" begin
    for order in (12, 28)
        cache = LocalBarycentricCache(NODES; order = order)
        p(x) = ((x - 3.7) / 2)^order + 3 * ((x - 3.7) / 2)^(order ÷ 2) + 1
        vals = p.(NODES)
        for i in (100, 200, 300)
            x = (NODES[i] + NODES[i + 1]) / 2
            @test abs(interpolate_local(cache, vals, x) - p(x)) <= 1e-12 * abs(p(x))
        end
    end
end

@testset "accuracy floors vs BigFloat reference" begin
    vals = fexact.(NODES)                       # Float64 data
    cache = LocalBarycentricCache(NODES; order = 28)
    # interior: Float64 data floor
    xs_mid = [(NODES[i] + NODES[i + 1]) / 2 for i in 150:250]
    @test maximum(relerr(interpolate_local(cache, vals, x), x) for x in xs_mid) < 5e-14
    # tapered edge: stays near the data floor (was ~1e-10 with full clamped windows)
    xs_edge = [(NODES[i] + NODES[i + 1]) / 2 for i in 1:14]
    @test maximum(relerr(interpolate_local(cache, vals, x), x) for x in xs_edge) < 1e-13

    # genuine extended precision: BigFloat nodes+values+weights → far below 1e-16
    bnodes = BigFloat.(NODES)
    bvals = fexact.(bnodes)
    bcache = LocalBarycentricCache(bnodes; order = 28)
    x = (bnodes[200] + bnodes[201]) / 2
    @test abs((interpolate_local(bcache, bvals, x) - fexact(x)) / fexact(x)) < BigFloat(10)^-45

    # Double64-weight path on Float64 data: limited by the data, not arithmetic
    dcache = LocalBarycentricCache(NODES; order = 28, weight_type = Double64)
    xd = Double64(NODES[200]) + Double64(0.5) * H
    r = interpolate_local(dcache, vals, xd)
    @test r isa Double64
    @test relerr(r, Float64(xd)) < 5e-14
end

@testset "exact node hits" begin
    vals = fexact.(NODES)
    cache = LocalBarycentricCache(NODES; order = 28)
    @test interpolate_local(cache, vals, NODES[57]) === vals[57]
    m = hcat(vals, 2 .* vals)
    @test interpolate_local(cache, m, NODES[57]) == [vals[57], 2vals[57]]
    # Double64 query exactly on a Float64 node (exact lift ⇒ exact hit)
    dcache = LocalBarycentricCache(NODES; order = 28, weight_type = Double64)
    @test Float64(interpolate_local(dcache, vals, Double64(NODES[57]))) === vals[57]
end

@testset "matrix path ≡ vector path (per column, same window/coeffs)" begin
    Random.seed!(1)
    m = randn(N, 7)
    cache = LocalBarycentricCache(NODES; order = 28)
    for x in (NODES[3] + 0.3H, NODES[200] + 0.7H, NODES[end - 1] + 0.5H)
        vm = interpolate_local(cache, m, x)
        for c in 1:7
            @test isapprox(vm[c], interpolate_local(cache, m[:, c], x); rtol = 1e-13, atol = 1e-300)
        end
    end
end

@testset "small orders / small grids degenerate correctly" begin
    # order ≤ edge_order_min ⇒ taper is a no-op (edge_min == stencil)
    c12 = LocalBarycentricCache(NODES; order = 12)
    @test c12.edge_min == c12.stencil == 13
    # linear interpolation on two nodes
    c1 = LocalBarycentricCache([0.0, 1.0]; order = 1)
    @test interpolate_local(c1, [0.0, 2.0], 0.25) ≈ 0.5
    # stencil clamps to the node count
    c = LocalBarycentricCache([0.0, 0.5, 1.0]; order = 28)
    @test c.stencil == 3
end

@testset "mixed query types promote exactly" begin
    vals = fexact.(NODES)
    cache = LocalBarycentricCache(NODES; order = 28)
    x64 = NODES[200] + 0.37H
    rF = interpolate_local(cache, vals, x64)
    rD = interpolate_local(cache, vals, Double64(x64))       # F64 weights, D64 coeffs
    @test rD isa Double64
    @test Float64(rD) ≈ rF rtol = 1e-14
end

@testset "argument validation" begin
    @test_throws ArgumentError LocalBarycentricCache([1.0, 0.5, 2.0]; order = 2)   # unsorted
    @test_throws ArgumentError LocalBarycentricCache([1.0, 1.0, 2.0]; order = 2)   # duplicate
    # weight overflow guard: degree 28 on a ~1e-6-step Float32 grid overflows its exponent range
    tiny = Float32.(1.0 .+ (0:40) .* 1.0f-6)
    @test_throws ArgumentError LocalBarycentricCache(tiny; order = 28)
    cache = LocalBarycentricCache(NODES; order = 12)
    @test_throws DimensionMismatch interpolate_local(cache, zeros(N - 1), 3.5)
    @test_throws DimensionMismatch interpolate_local(cache, zeros(N - 1, 3), 3.5)
end

@testset "callable sugar" begin
    vals = fexact.(NODES)
    cache = LocalBarycentricCache(NODES; order = 28)
    x = NODES[100] + 0.6H
    @test cache(vals, x) === interpolate_local(cache, vals, x)
end

end
