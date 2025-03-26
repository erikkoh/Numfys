using Plots
using LinearAlgebra
using SparseArrays
using Printf

function euler_explicit(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/nx
    r = (λ^2/(τ)) * (dt/(dx^2))
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    for t in 2:nt
        solution[t, 1] = solution[t-1,1]*(1-2*r-dt/(2τ)) +2*r*solution[t-1,2]
        for x in 2:nx-1
            solution[t,x] = solution[t, x] = solution[t-1,x]*(1-2*r-dt/(2τ)) +r*solution[t-1,x+1] +r*solution[t-1,x-1]
        end
        solution[t,end] = solution[t-1,end]*(1-2*r-dt/(2τ)) +2*r*solution[t-1,end-1]
    end
    return solution
end 


#Each step is calculated with the A*x_n+1 = x_n which can than be used to solve for the next step with: x_n+1 = A^-1*x_n
function euler_implicit(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/(nx)
    r = (λ^2/(τ)) * (dt/(dx^2))
    α = 1 + 2*r + dt/(2τ)
    β = -r
    # @printf("α = %g, β = %g\n", α, β)

    diagonals = [β*ones(nx-1),  α*ones(nx),  β*ones(nx-1)]
    A = spdiagm(-1 => diagonals[1], 0 => diagonals[2], 1 => diagonals[3])
    A[1,2] = -2*r
    A[end,end-1] = -2*r
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    for t in 2:nt
        b_int = solution[t-1, :]
        interior_solution = A \ b_int
        solution[t,:] = interior_solution
    end
    return solution
end



#Each step is calculated with the A_1*x_n+1 = A_2*x_n which can than be used to solve for the next step with: x_n+1 = A_1^-1*A_2*x_n
function crank_nicolson(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/(nx)
    # println("dt = $dt and dx = $dx ")
    r = (λ^2/(τ)) * (dt/(2*dx^2)) 
    α = 1 + 2*r + dt/(2τ)
    β = -r
    γ = 1 - 2*r - dt/(2τ)
    δ = r
    
    diagonals_next = [fill(β, nx-1), fill(α,nx),  fill(β,nx-1)]
    digonals_previous = [δ*ones(nx-1), γ*ones(nx), δ*ones(nx-1)]
    A_next = spdiagm(-1 => diagonals_next[1], 0 => diagonals_next[2], 1 => diagonals_next[3])
    A_next[1, 2] = -2*r 
    A_next[end,end-1] = -2*r
    A_previous = spdiagm( -1 => digonals_previous[1], 0 => digonals_previous[2], 1 => digonals_previous[3] )
    A_previous[1, 2] = 2*r
    A_previous[end,end-1] = 2*r
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    for t in 2:nt
        b_int = A_previous*solution[t-1,:]
        interior_solution = A_next \ b_int
        solution[t,:] = interior_solution
    end
    return solution
end

function analytical_solution(V0, V_0_tilde, nx, nt, a, b, time_stop, x0=0.5, λ=1.0, τ=1.0)
    dt = time_stop/nt
    x_list = range(a, b, length=nx)
    t_list = range(dt, time_stop, length=nt)
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    analytical = (x, t) -> V_0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
                           exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
    for t_step in 2:nt
        t_current = t_list[t_step]
        solution[t_step, :] = analytical.(x_list, t_current)
    end
    
    return solution
end

a = 0.0
b = 1.0
nx = 100
nt = 2_000
time_stop = 1.0
dt = time_stop/nt
λ = 1.0
τ = 1.0
x0 = 0.5

# Function to analyze the change in error with decreasing dx
function analyze_error_with_dx(a, b, x0, λ, τ, nx_list, nt, t_end)
    errors_crank = zeros(Float64, length(nx_list))
    errors_implicit = zeros(Float64, length(nx_list))
    errors_explicit = zeros(Float64, length(nx_list))
    dx_values = zeros(Float64, length(nx_list))
    for (i, nx) in enumerate(nx_list)
        println("Started on nx = $nx")
        dx = (b - a) / nx
        dx_values[i] = dx
        
        # Define initial condition
        x = range(a, b, length=nx)
        V0_tilde = 1.0
        analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
                               exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
        V0 = analytical.(x, t_end / nt)
        
        # Compute numerical solution using Crank-Nicolson
        sol_crank = crank_nicolson(V0, nx, nt, t_end, a, b, λ, τ)
        sol_implicit = euler_implicit(V0, nx, nt, t_end, a, b, λ, τ)
        sol_explicit = euler_explicit(V0, nx, nt, t_end, a, b, λ, τ)
        
        # Compute analytical solution directly using the lambda function
        sol_analytical = analytical.(x, t_end)
        
        # Compute error as the L2 norm of the difference at the final time step
        errors_crank[i] = norm(sol_crank[end, :] - sol_analytical, 2) / sqrt(nx)
        errors_implicit[i] = norm(sol_implicit[end, :] - sol_analytical, 2) / sqrt(nx)
        errors_explicit[i] = norm(sol_explicit[end, :] - sol_analytical, 2) / sqrt(nx)
    end
    
    # Plot error vs dx
    plot(dx_values, errors_crank, lw=2, marker=:o, xlabel="Δx", ylabel="Error", title="Error vs Δx", xscale=:log10, yscale=:log10, label="Crank-Nicolson")
    plot!(dx_values, errors_implicit, lw=2, marker=:o, label="Implicit Euler")
    plot!(dx_values, errors_explicit, lw=2, marker=:o, label="Explicit Euler")
    
    # Add trend line for Δx^2
    trend_line = dx_values .^ 1.5
    plot!(dx_values, trend_line , lw=2, ls=:dash, label="Δx² trend")
    savefig("error_vs_dx.png")
    println("Done!")
end

# Function to analyze the change in error with decreasing dt
function analyze_error_with_dt(a, b, x0, λ, τ, nx, nt_list, t_end)
    errors_crank = zeros(Float64, length(nt_list))
    errors_implicit = zeros(Float64, length(nt_list))
    errors_explicit = zeros(Float64, length(nt_list))
    dt_values = zeros(Float64, length(nt_list))
    for (i, nt) in enumerate(nt_list)
        println("Started on nt = $nt")
        dt = t_end / nt
        dt_values[i] = dt
        
        # Define initial condition
        x = range(a, b, length=nx)
        V0_tilde = 1.0
        analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
                               exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
        V0 = analytical.(x, dt)
        
        # Compute numerical solution using Crank-Nicolson
        sol_crank = crank_nicolson(V0, nx, nt, t_end, a, b, λ, τ)
        sol_implicit = euler_implicit(V0, nx, nt, t_end, a, b, λ, τ)
        sol_explicit = euler_explicit(V0, nx, nt, t_end, a, b, λ, τ)
        
        # Compute analytical solution directly using the lambda function
        sol_analytical = analytical.(x, t_end)
        
        # Compute error as the L2 norm of the difference at the final time step
        errors_crank[i] = norm(sol_crank[end, :] - sol_analytical, 2) / sqrt(nx)
        errors_implicit[i] = norm(sol_implicit[end, :] - sol_analytical, 2) / sqrt(nx)
        errors_explicit[i] = norm(sol_explicit[end, :] - sol_analytical, 2) / sqrt(nx)
    end
    
    # Plot error vs dt
    plot(dt_values, errors_crank, lw=2, marker=:o, xlabel="Δt", ylabel="Error", title="Error vs Δt", xscale=:log10, yscale=:log10, label="Crank-Nicolson")
    plot!(dt_values, errors_implicit, lw=2, ls=:dash, marker=:o, label="Implicit Euler")
    plot!(dt_values, errors_explicit, lw=2, marker=:o, label="Explicit Euler")
    
    # Add trend lines for Δt and Δt²
    trend_line_dt = dt_values 
    trend_line_dt2 = dt_values .^ 2
    plot!(dt_values, trend_line_dt, lw=2, ls=:dash, label="Δt trend")
    plot!(dt_values, trend_line_dt2, lw=2, ls=:dash, label="Δt² trend")
    savefig("error_vs_dt.png")
    println("Done!")
end




nx_list = [200, 300, 400, 500, 600, 700, 800, 1_000]
nt_list = [2_000, 20_000, 40_000, 80_000]#To ensure stability

analyze_error_with_dx(a, b, x0, λ, τ, nx_list, nt, time_stop)
analyze_error_with_dt(a, b, x0, λ, τ, nx, nt_list ,time_stop)

# t = range(0,time_stop, length = nt)
# x = range(a,b, length=nx)
# # V0 = exp.(- (x .- x0).^2 ./ (4*λ^2/τ)*dt .- dt/τ) 
# # V0 = V0 ./ (σ*sqrt(4π))
# # V0_tilde =  sqrt(4 * π * (λ^2/τ))
# V0_tilde = 1.0
# analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
#                            exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
# V0 = analytical.(x, dt)


# sol_explicit = euler_explicit(V0, nx, nt, time_stop, a, b)
# sol_implicit = euler_implicit(V0, nx, nt, time_stop, a, b)
# sol_crank = crank_nicolson(V0, nx, nt, time_stop, a, b)
# sol_analytical = analytical_solution(V0, V0_tilde, nx, nt, a, b, time_stop)

# d1 = plot(x, sol_explicit[end,:], lw=2, ls=:dot, label="explicit euler", xlabel="x", ylabel="V", title="Cable Equation: Final Profiles (Neumann BCs)")
# plot!(x, sol_implicit[end,:], lw=2, ls =:dashdot, label="Euler implicit")
# plot!(x,sol_crank[end,:], lw=2, label = "Crank-Nicolson")
# plot!(x, analytical.(x,time_stop), lw=2, ls=:dash, label="Analytical")
# savefig(d1, "cable_eq_final_profiles.png")

# d2 = plot(t,sol_explicit[:,50], ls=:dot,  xlabel = "time", ylabel = "V" ,title="Potential for middle point over time", label= "Euler explicit")
# plot!(t,sol_implicit[:,50], ls =:dashdot,label="Euler implicit")
# plot!(t,sol_crank[:,50], label = "Crank-Nicolson")
# plot!(t,sol_analytical[:,50], ls=:dash, label = "Analytical")
# savefig(d2, "cable_eq_first_point.png")

# # heatmap(x,t, sol_explicit)
# # heatmap(x,t, sol_implicit)
# heatmap(x, t, sol_crank, xlabel="x", ylabel="t", zlabel="V", title="Surface plot of V(x, t)")
# heatmap(x,t, sol_analytical,)
