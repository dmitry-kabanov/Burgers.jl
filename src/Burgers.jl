module Burgers

using LaTeXStrings
using Plots

export compute_init, solve, run, plot_init_and_final_solutions, runAndPlot

function compute_init(; A=1.0, theta=np.pi / 2, T=0.4, N=101)
x = range(0, 1.0, length=N)
dx = 1.0 / N

u0 = A * sin.(2 * pi * x .+ theta)

cfl = 0.5

dt = dx * cfl
times = 0:dt:T

return x, dx, u0, times, dt
end

function solve(x, dx, u0, times, dt)
    u = u0
    n = length(u0)
    u_new = similar(u)
    for i = 2:length(times)
        f = 0.5 * u.^2
        local_ss = maximum(abs.(u))

        f_hat = 0.5 * (f[1:n-1] + f[2:n]) - 0.5 * local_ss * (u[2:n] - u[1:n-1])
        f_hat_plus = f_hat[2:n-1]
        f_hat_minus = f_hat[1:n-2]
        u_new[2:n-1] = u[2:n-1] - dt/dx * (f_hat_plus - f_hat_minus)
        
        local_ss_bdy = max(u[1], u[n])
        f_rb = 0.5 * (f[1] + f[n]) - 0.5 * local_ss_bdy * (u[1] - u[n])
        f_lb = f_rb
        
        u_new[1] = u[1] - dt / dx * (f_hat_minus[1] - f_lb)
        u_new[n] = u[n] - dt / dx * (f_rb - f_hat_plus[n-2])
        
        u = copy(u_new)
    end
    
    return u
end

function run(; A=1.0, theta=pi / 2.0, T=0.4, N=501)
x, dx, u0, times, dt= compute_init(A=A, theta=theta, T=T, N=N)
u = solve(x, dx, u0, times, dt)
return x, u
end

function plot_init_and_final_solutions(x, u0, u, T)
    plot(x, u0, label=L"t = 0", linestyle=:dash, linewidth=3)
    plot!(x, u, label=L"t = %$T", linewidth=3)
end

function runAndPlot(A=1.0, theta=pi / 2.0, T=0.4, N=501)
x, dx, u0, times, dt= compute_init(A=A, theta=theta, T=T, N=N)
u = solve(x, dx, u0, times, dt)
plot_init_and_final_solutions(x, u0, u, T)
end

end # module Burgers
