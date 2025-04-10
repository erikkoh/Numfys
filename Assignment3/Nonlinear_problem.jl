using Plots
using LinearAlgebra
using SparseArrays


if !isdir("plots")
    mkdir("plots")
end



function g_K_reaction(V, V_star = -40.0, γ = 0.5)
    return 100 ./(1 .+ exp.((V_star .- V)*γ)) + 1/5
end

function updating_h(h, V, gamma = 1.0, V_mem = -70.0)
    V_m = -V_mem - V
    α_h = gamma*(0.07 * exp(-(V_m) / 20))
    β_h = gamma*(1.0 / (1 + exp((40-V_m) / 10)))
    return (α_h * (1 - h) - β_h * h)
end

function g_Na(V, γ = 0.5, V_star = -40.0)	
    conduct =  100.0/(1 + exp(γ*(V_star - V))) + 1/5
    return conduct 
end


function reaction_term_h_dependent(V, h, g_K = 5.0,nu_Na = +56, nu_K = -76)
    # h is ment to either disable or enable the sodium channel
    return g_Na(V)/g_K*h*(V - nu_Na) + (V -  nu_K)
end


#todo: fixed point iteration and implement crank nickelson for relaxation perameter
function crank_nicolson(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, x0::Float64, λ = 0.18 , τ = 2.0)	
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
        reaction = reaction_term_h_dependent.(solution[t-1, :], h_array)
        b_int = A_previous*solution[t-1, :] .- dt/τ.* reaction
        interior_solution = A_next \ b_int
        solution[t,:] = interior_solution
    end
    return solution
end

function crank_nicolson_fixed_point(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, x_0::Float64, λ = 0.18 , τ = 2.0)
    max_iter = 100
    dt = t_end/nt
    dx = (b-a)/(nx)
    tol = dx*dt*10^(-6)
    r = (λ^2/(τ)) * (dt/(2*dx^2)) 
    α = 1 + 2*r 
    β = -r
    γ_2 = 1 - 2*r 
    δ = r
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
    calculate_NA_channles_plus = Int(round((x_0+0.25)/dx))
    calculate_NA_channles_neg = Int(round((x_0+0.25)/dx))
    h_array[calculate_NA_channles_plus:end] .= 1.0
    h_array[1:calculate_NA_channles_neg] .= 1.0
    for t in 2:nt
        previus_reaction = reaction_term_h_dependent.(solution[t-1, :], h_array)
        b_int = A_previous*solution[t-1, :] .- dt/τ.*previus_reaction# the average of the nonlinear term is simply the nonlinear for the previous value 
        interior_solution = A_next \ b_int
        v_temp = interior_solution 
        previous_vec = A_previous*solution[t-1, :] 
        # Fixed point iteration
        for _j in 1:max_iter
            b_int = previous_vec .- dt/(2*τ).*(previus_reaction + reaction_term_h_dependent.(v_temp, h_array))
            interior_solution = A_next \ b_int
            v_next  = interior_solution
            if norm(v_next-v_temp,2) < tol
                v_temp = v_next 
                break
            end
            v_temp = v_next
        end
        solution[t,:] = v_temp
    end
    return solution
end

function find_lowest_inpulse(V_mem::Int , nx::Int, t_end::Float64, a::Float64, b::Float64, x_0::Float64, λ = 0.18, τ = 2.0)
    for V_appl in -70.0:0.1:0.0
        V0 = (V_appl - V_mem) .* exp.(-((x .- x_0).^2) ./ (2*λ^2)) .+ V_mem
        solution = crank_nicolson_fixed_point(V0, nx, Int(t_end/dx), t_end, a, b, x_0, λ, τ)
        if maximum(solution[:, critical_index]) > 0.0
            println("Lowest V_appl to induce action potential: ", V_appl)
            return V_appl
        end
    end
    println("No action potential induced within the tested range.")
    return nothing
end




# Define parameters for action potential modeling
a = 0.0
b = 10.0
λ = 0.18
V_mem = -70
nx = 1_00
nt = 2_000
dx = (b-a)/nx
time_stop = 10.0
x0 = (b-a)/2
t = range(0, time_stop, length=nt)
x = range(a, b, length=nx)
critical_index = Int(round(0.25/dx))
V_appl = find_lowest_inpulse(V_mem, nx, time_stop, a, b, x0, λ)
V0 = (V_appl - V_mem) .* exp.(-((x .- x0).^2) ./ (2*λ^2)) .+ V_mem



# Solve using Crank-Nicolson fixed point method
solution_crank = crank_nicolson_fixed_point(V0, nx, nt, time_stop, a, b, x0, λ)


# Plot results
d3 = plot(t, solution_crank[:, critical_index], xlabel="time", ylabel="V", title="Potential for first point over time", label="Crank Nicolson")
savefig(d3, "./Assignment3/plots/First_point_over_time")

d1 = heatmap(x, t, solution_crank, xlabel="x", ylabel="time", title="Heatmap of the Crank-Nicolson Fixed Point Simulation", colorbar_title="V")
savefig(d1, "./Assignment3/plots/Heatplot_test")


