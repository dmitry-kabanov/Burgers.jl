using Test
using Burgers

@testset "Sanity checks" begin
    N = 3
    x, dx, u0, times, dt = compute_init(A = 1.0, theta = pi / 2, T = 0.4, N = N)
    @test length(x) == N
    @test length(u0) == N

    u = solve(x, dx, u0, times, dt)
    @test length(u) == N
end

@testset "Test observed order of accuracy" begin
    N_best = 4096
    x, dx, u0, times, dt = compute_init(A = 1.0, theta = pi / 2, T = 0.4, N = N_best)
    u_best = solve(x, dx, u0, times, dt)

    resolutions = [32, 64, 128, 256, 512]
    errors_Linf = Vector{Float64}(undef, length(resolutions))
    for i = 1:length(resolutions)
        N = resolutions[i]
        x, dx, u0, times, dt = compute_init(A = 1.0, theta = pi / 2, T = 0.1, N = N)
        u = solve(x, dx, u0, times, dt)
        step = trunc(Int, (N_best - 1) / (N - 1))
        println("Step = ", step)
        errors_Linf[i] = maximum(abs.(u - u_best[1:step:end]))
    end

    obs_orders = -log.(errors_Linf[1:end-1] / errors_Linf[2:end]) / log(2)

    for order in obs_orders
        @test 1.8 <= order
        @test order <= 2.2
    end
end
