using Plots
using LinearAlgebra
using SparseArrays
using Printf
using ProgressBars
using Statistics

function euler_explicit(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/nx
    r = (λ^2/(τ)) * (dt/(dx^2))
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    for t in 2:nt
        solution[t, 1] = solution[t-1,1]*(1-2*r-dt/(τ)) +2*r*solution[t-1,2]
        for x in 2:nx-1
            solution[t,x] = solution[t-1,x]*(1-2*r-dt/(τ)) +r*solution[t-1,x+1] +r*solution[t-1,x-1]
        end
        solution[t,end] = solution[t-1,end]*(1-2*r-dt/(τ)) +2*r*solution[t-1,end-1]
    end
    return solution
end 


#Each step is calculated with the A*x_n+1 = x_n which can than be used to solve for the next step with: x_n+1 = A^-1*x_n
function euler_implicit(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/(nx)
    r = (λ^2/(τ)) * (dt/(dx^2))
    α = 1 + 2*r + dt/(τ)
    β = -r

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
    r = (λ^2/(τ)) * (dt/(2*dx^2)) 
    α = 1 + 2*r + dt/(2τ)
    β = -r
    γ = 1 - 2*r - dt/(2τ)
    δ = r
    
    diagonals_next = [fill(β, nx-1), fill(α,nx),  fill(β,nx-1)]
    digonals_previous = [fill(δ,nx-1), fill(γ,nx), fill(δ,nx-1)]
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
    t_list = range(1.0, time_stop+1.0, length=nt)
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


# Function to analyze the change in error with decreasing dx
function analyze_error_with_dx(a, b, x0, λ, τ, nx_list, nt, t_end)
    errors_crank = zeros(Float64, length(nx_list))
    errors_implicit = zeros(Float64, length(nx_list))
    errors_explicit = zeros(Float64, length(nx_list))
    dx_values = zeros(Float64, length(nx_list))
    V0_tilde = 1.0
    analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
    exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
    # best_numerical_solution = golden_numerical_approximation_in_pos(analytical.(x,0.001), nt, t_end, a, b)
    for (i, nxs) in enumerate(nx_list)
        println("Started on nx = $nxs")
        dx = (b - a) / nxs
        dx_values[i] = dx
        
        # Define initial condition
        x = range(a, b, length=nxs)
        V0 = analytical.(x, 1.0)
        
        # Compute numerical solution using Crank-Nicolson
        sol_crank = crank_nicolson(V0, nxs, nt, t_end, a, b, λ, τ)
        sol_implicit = euler_implicit(V0, nxs, nt, t_end, a, b, λ, τ)
        sol_explicit = euler_explicit(V0, nxs, nt, t_end, a, b, λ, τ)
        
        # approx_solution = golden_numerical_approximation(V0, nxs, t_end, a, b)
        sol_analytical = analytical.(x, time_stop+1.0)
        
        # Compute error  at the final time step
        errors_crank[i] = sqrt(mean((sol_crank[end,:] - sol_analytical).^2))
        errors_explicit[i] = sqrt(mean((sol_explicit[end,:] - sol_analytical).^2))
        errors_implicit[i] = sqrt(mean((sol_implicit[end,:] - sol_analytical).^2))
    end
    # Plot error vs dx
    plot(dx_values, errors_crank, lw=2, marker=:o, xlabel="Log(Δx)", ylabel="Log(Error)", title="Error vs Δx after $t_end time, and $nt time steps", xscale=:log10, yscale=:log10, label="Crank-Nicolson")
    plot!(dx_values, errors_implicit, lw=2, marker=:o, label="Implicit Euler")
    plot!(dx_values, errors_explicit, lw=2, marker=:o, label="Explicit Euler")
    
    # Add trend line for Δx^2
    trend_line = dx_values .^ 2
    plot!(dx_values, trend_line*errors_crank[1]/trend_line[1] , lw=2, ls=:dash, label="Δx² trend")
    savefig("error_vs_dx.png")
    println("Done!")
end

#The Crank Nicolson doesn't seem to converge with Δt^2 but rather with Δt at lower time values and converges to Δt^2 at higher time values (t_end = 100).
# Function to analyze the change in error with decreasing dt
function analyze_error_with_dt(a, b, x0, λ, τ, nx, nt_list, t_end, late::Bool=false)
    # t_end = t_end + 1.0
    errors_crank = zeros(Float64, length(nt_list))
    errors_implicit = zeros(Float64, length(nt_list))
    errors_explicit = zeros(Float64, length(nt_list))
    dt_values = zeros(Float64, length(nt_list))
    V0_tilde = 1.0
    analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
                            exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
    x = range(a, b, length=nx)
    
    V0 = analytical.(x, 1.0)

    best_numerical_solution = analytical.(x, t_end+1.0)
    for (i, nt) in enumerate(nt_list)
        dt = t_end / nt
        dt_values[i] = dt
        println("Started on nt = $nt and Δt = $dt" )
        if !late    
            sol_implicit = euler_implicit(V0, nx, nt, t_end, a, b, λ, τ)
            sol_explicit = euler_explicit(V0, nx, nt, t_end, a, b, λ, τ)
            errors_implicit[i] = sqrt(mean((sol_implicit[end,:] - best_numerical_solution).^2))
            errors_explicit[i] = sqrt(mean((sol_explicit[end,:] - best_numerical_solution).^2))
        end

        sol_crank = crank_nicolson(V0, nx, nt, t_end, a, b, λ, τ)
        errors_crank[i] = sqrt(mean((sol_crank[end,:] - best_numerical_solution).^2))
    end
    println(log10.(errors_crank))

    plot(log10.(dt_values), log10.(errors_crank), lw=2, marker=:o, xlabel="Log(Δt)", ylabel="Log(Error)", title="Error vs Δt after $t_end and $nx position steps ", label="Crank-Nicolson")
    if !late
        plot!(log10.(dt_values), log10.(errors_implicit), lw=2, ls=:dash, marker=:o, label="Implicit Euler")
        plot!(log10.(dt_values), log10.(errors_explicit), lw=2, marker=:o, label="Explicit Euler")
    end
    # Add trend lines for Δt and Δt²
    trend_line_dt = dt_values 
    trend_line_dt2 = dt_values .^ 2
    if !late
        plot!(log10.(dt_values), log10.(trend_line_dt*errors_crank[1]/trend_line_dt[1]), lw=2, ls=:dash, label="Δt trend")
    end
    plot!(log10.(dt_values), log10.(trend_line_dt2*errors_crank[1]/trend_line_dt2[1]), lw=2, ls=:dash, label="Δt² trend")
    if !late
    savefig("error_vs_dt.png")
    else
        savefig("error_vs_dt_late.png")
    end
    println("Done!")
end


function analyse_error_with_time(a,b,x0,nx, nt, t_end, λ=1.0, τ=1.0)
    println("Analysing error with time using Δt = $(t_end/nt) and Δx = $((b-a)/nx)") 
    V0_tilde = 1.0
    dt = t_end/nt
    analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
    exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
    x = range(a, b, length=nx)
    V0 = analytical.(x, 1e-6)
    sol_crank = crank_nicolson(V0, nx, nt, t_end, a, b, λ, τ)
    sol_implicit = euler_implicit(V0, nx, nt, t_end, a, b, λ, τ)
    sol_explicit = euler_explicit(V0, nx, nt, t_end, a, b, λ, τ)
    
    # Compute analytical solution directly using the lambda function

    errors_crank = zeros(Float64, nt)
    errors_implicit = zeros(Float64, nt)
    errors_explicit = zeros(Float64, nt)

    for i in 2:nt
        errors_crank[i] = sqrt(mean((sol_crank[i,:] - analytical.(x,i*dt)).^2))
        errors_implicit[i] = sqrt(mean((sol_implicit[i,:] - analytical.(x,i*dt)).^2))
        errors_explicit[i] = sqrt(mean((sol_explicit[i,:] - analytical.(x,i*dt)).^2))
    end
    t = range(0.0, t_end, length=nt)
    println("Started plotting")
    plot(t , log10.(errors_crank), lw=2 , label="Crank-Nicolson", title ="Error vs Time , with $nt time and position steps", xlabel="Time", ylabel="Log10(Error)", legend=:topright)
    plot!(t, log10.(errors_implicit), lw=2 , label="Implicit Euler")
    plot!(t, log10.(errors_explicit), lw=2 , label="Explicit Euler")
    savefig("error_vs_time.png")
    println("Done!")
end


nx_list = [5, 10 ,20, 50, 100, 200]
dt_list = [0.02,0.05,0.2,0.5, 1, 2, 5, 10]

a = 0.0
b = 100.0
nx = 200
nt = 200
time_stop = 1.0
time_stop_2 = 125.0
λ = 1.0
τ = 1.0
x0 = (b-a)/2
dt_n_list =Int.(round.(time_stop ./ dt_list[1:end-3]))# the higher dt values reulsts in vectors of length 0
dt_n_list_2 =Int.(round.(time_stop_2 ./ dt_list))

t = range(0,time_stop, length = nt)
x = range(a,b, length=nx)
V0_tilde = 1.0
analytical = (x, t) -> V0_tilde / sqrt(4 * π * (λ^2/τ) * t) *
exp( - ((x - x0)^2 / (4 * (λ^2/τ) * t)) - t/τ)
V0 = analytical.(x, 1.0)


analyse_error_with_time(a,b,x0,nx, nt, 1.0, λ, τ)
analyze_error_with_dx(a, b, x0, λ, τ, nx_list, nt, time_stop)
analyze_error_with_dt(a, b, x0, λ, τ, nx, dt_n_list ,time_stop)
analyze_error_with_dt(a, b, x0, λ, τ, nx, dt_n_list_2 ,time_stop_2, true)



sol_explicit = euler_explicit(V0, nx, nt, time_stop, a, b)
sol_implicit = euler_implicit(V0, nx, nt, time_stop, a, b)
sol_crank = crank_nicolson(V0, nx, nt, time_stop, a, b)
sol_analytical = analytical_solution(V0, V0_tilde, nx, nt, a, b, time_stop)

d1 = plot(x, sol_crank[end,:], lw=2, ls=:dot, label="Crank-Nicolson", xlabel="x", ylabel="V", title="Cable Equation: Final Profiles after time $time_stop")
plot!(x, sol_implicit[end,:], lw=2, ls =:dashdot, label="Euler implicit")
plot!(x,sol_explicit[end,:], lw=2, label = "Euler explicitn")
plot!(x, analytical.(x,time_stop+1.0), lw=2, ls=:dash, label="Analytical")
savefig(d1, "cable_eq_final_profiles.png")

d2 = plot(t,sol_crank[:,50], ls=:dot,  xlabel = "time", ylabel = "V" ,title="Potential for middle point over time", label= "Crank-Nicolson")
plot!(t,sol_implicit[:,50], ls =:dashdot,label="Euler implicit")
# plot!(t,sol_explicit[:,50], label = "Euler explicit")
plot!(t,sol_analytical[:,50], ls=:dash, label = "Analytical")
savefig(d2, "cable_eq_first_point.png")

