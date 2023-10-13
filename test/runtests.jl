using Test
using Burgers

@testset "Sanity checks" begin
    N = 3
    x, dx, u0, times, dt = compute_init(A=1.0, theta=pi/2, T=0.4, N=N)
    @test length(x) == N
    @test length(u0) == N
    
    u = solve(x, dx, u0, times, dt)
    @test length(u) == N
end