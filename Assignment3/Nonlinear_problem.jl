using Plots
using LinearAlgebra
using SparseArrays

function g_Na(V,γ =0.5, V_star = -40.0)
    result =  100 ./(1 .+ exp.(γ.*(V_star .- V))) .+ 1/5
    # println("g_Na = ", result)
    return result
end

function reaction_term(V, g_K = 5.0, nu_Na = 56, nu_K = -76)
    return (((g_Na(V))*(V- nu_Na) + g_K*(V - nu_K))) 
end

function g_K_reaction(V, V_star = -40.0, γ = 0.5)
    return 100 ./(1 .+ exp.((V_star .- V)*γ)) + 1/5
end

function updating_h(h, V, gamma = 1.0, V_mem = -70.0)
    V_m = -V_mem - V
    α_h = gamma*(0.07 * exp(-(V_m) / 20))
    β_h = gamma*(1.0 / (1 + exp((30-V_m) / 10)))
    return (α_h * (1 - h) - β_h * h)
end

function g_Na_2(V, γ = 0.5, V_star = -40.0)	
    V_m = -V_mem - V
    conduct =  (100.0/(1 + exp(γ*(V_star - V_m))) + 1/5)
    return conduct 
end


function reaction_term_h_dependent(V, h, g_K = 5.0,nu_Na = +56, nu_K = -76)
    # println(g_Na_2(V, x)/g_K) 
    # h here is an relaxation parameter
    return g_Na_2(V)*h*(V - nu_Na) + g_K * (V -  nu_K)
end


#todo: fixe point iteration and implement crank nickelson for relaxation perameter
function crank_nicolson(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 0.18 , τ = 2.0)	
    dt = t_end/nt
    dx = (b-a)/(nx)
    r = (λ^2/(τ)) * (dt/(2*dx^2)) 
    α = 1 + 2*r 
    β = -r
    γ_2 = 1 - 2*r 
    δ = r
    h_array = fill(0.0, nx)
    calculate_NA_channles = Int(0.10/dx)
    h_array[calculate_NA_channles:end] .= 1.0
    diagonals_next = [β*ones(nx-1),  α*ones(nx),  β*ones(nx-1)]
    digonals_previous = [δ*ones(nx-1), γ_2*ones(nx), δ*ones(nx-1)]
    A_next = spdiagm(-1 => diagonals_next[1], 0 => diagonals_next[2], 1 => diagonals_next[3])
    A_next[1,1] = 1+2*r
    A_next[1,2] = -2*r
    A_next[end,end] = 1+2*r
    A_next[end,end-1] = -2*r
    A_previous = spdiagm( -1 => digonals_previous[1], 0 => digonals_previous[2], 1 => digonals_previous[3] )
    A_previous[1,1] = 1-2*r
    A_previous[1,2] = 2*r
    A_previous[end,end-1] = 2*r
    A_previous[end,end] = 1-2*r
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    for t in 2:nt
        h_array = h_array .+ dt*updating_h.(h_array, solution[t-1, :])
        reaction = reaction_term_x_dependent.(solution[t-1, :], h_array)
        b_int = A_previous*solution[t-1, :] .- dt/τ.* reaction
        interior_solution = A_next \ b_int
        solution[t,:] = interior_solution
    end
    return solution
end

function crank_nicolson_fixed_point(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 0.18 , τ = 2.0)
    max_iter = 100
    dt = t_end/nt
    dx = (b-a)/(nx)
    tol = dx*dt*10^(-6)
    r = (λ^2/(τ)) * (dt/(2*dx^2)) 
    α = 1 + 2*r 
    β = -r
    γ_2 = 1 - 2*r 
    δ = r
    N = nx-2
    diagonals_next = [β*ones(nx-1),  α*ones(nx),  β*ones(nx-1)]
    digonals_previous = [δ*ones(nx-1), γ_2*ones(nx), δ*ones(nx-1)]
    A_next = spdiagm(-1 => diagonals_next[1], 0 => diagonals_next[2], 1 => diagonals_next[3])
    A_next[1,1] = 1+2*r
    A_next[1,2] = -2*r
    A_next[end,end] = 1+2*r
    A_next[end,end-1] = -2*r
    A_previous = spdiagm( -1 => digonals_previous[1], 0 => digonals_previous[2], 1 => digonals_previous[3] )
    A_previous[1,1] = 1-2*r
    A_previous[1,2] = 2*r
    A_previous[end,end-1] = 2*r
    A_previous[end,end] = 1-2*r
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    h_array = fill(0.0, nx)
    calculate_NA_channles = Int(0.10/dx)
    h_array[calculate_NA_channles:end] .= 1.0
    for t in 2:nt
        h_array = h_array .+ dt*updating_h.(h_array, solution[t-1, :])
        previus_reaction = reaction_term_h_dependent.(solution[t-1, :], h_array)
        b_int = A_previous*solution[t-1, :] .- dt/τ.*previus_reaction# the average of the nonlinear term is simply the nonlinear for the previous value 
        interior_solution = A_next \ b_int
        v_temp = interior_solution 
        #Fixed point iteration
        for j in 1:max_iter
            b_int = A_previous*solution[t-1, :] .- dt/(2*τ).*(previus_reaction + reaction_term_h_dependent.(v_temp, h_array))
            interior_solution = A_next \ b_int
            v_next  = interior_solution
            if norm(v_next-v_temp,Inf) < tol
                # println("Iteration $j, and at $t norm = ", norm(v_next - v_temp))
                v_temp = v_next
                break
            end
            v_temp = v_next
        end
        solution[t,:] = v_temp
    end
    return solution
end


# Define parameters for action potential modeling
a = 0.0
b = 1.0
λ = 0.18
V_mem = -70
V_appl = -30
nx = 1_00
nt = 20_000
dx = (b-a)/nx
time_stop = 10.0
x0 = 0.0
t = range(0, time_stop, length=nt)
x = range(a, b, length=nx)
critical_index = Int(0.25/dx)
V0 = (V_appl - V_mem) .* exp.(-((x .- x0).^2) ./ (2*λ^2)) .+ V_mem


# Solve using Crank-Nicolson fixed point method
solution_crank = crank_nicolson_fixed_point(V0, nx, nt, time_stop, a, b, λ)

# Plot results
d3 = plot(t, solution_crank[:, critical_index], xlabel="time", ylabel="V", title="Potential for first point over time", label="Crank Nicolson")
savefig(d3, "./Assignment3/plots/First_point_over_time")

d1 = heatmap(x, t, solution_crank, xlabel="x", ylabel="time", title="Heatmap of the Crank-Nicolson Fixed Point Simulation", colorbar_title="V")
savefig(d1, "./Assignment3/plots/Heatplot_test")

d2 = surface(x, t, solution_crank)
savefig(d2, "./Assignment3/plots/Surface_plot")


# a = 0.0
# b = 1.0
# λ = 0.18
# V_mem = -70
# V_appl = -50
# nx = 1_000
# nt = 2_000
# time_stop = 4.0
# x0 = 0.5
# t = range(0,time_stop, length = nt)
# x = range(a,b, length=nx)
# V0 = (V_appl - V_mem).* exp.( -((x .- x0).^2)./(2*λ^2)) .+ V_mem


# solution_crank = crank_nicolson_fixed_point(V0, nx, nt,time_stop, a, b, λ )
# d3 = plot(t,solution_crank[:,1], xlabel = "time", ylabel = "V" ,title="Potential for first point over time", label= "Crank Nicolson")
# savefig(d3, "./Assignment3/plots/First_point_over_time")

# d1 = heatmap(x, t, solution_crank, xlabel="x", ylabel="time", title="Heatmap of the Crank-Nicolson Fixed Point Simulation", colorbar_title="V")
# savefig(d1, "./Assignment3/plots/Heatplot_test")
# d2 = surface(x,t,solution_crank)
# savefig(d2, "./Assignment3/plots/Surface_plot")
