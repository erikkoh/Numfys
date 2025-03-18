using Plots
using LinearAlgebra
using SparseArrays
using Printf

function central_difference(V, pos, h)
    return (V[pos+1] - 2*V[pos] + V[pos-1])/(h^2)
end


#Applying the Neumann boundary condition by ussing a central difference where dV/dt = 0 
function apply_neumann!(V)
    V[1] = V[2]
    V[end] = V[end-1]
end


function euler_explicit(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/nx
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    apply_neumann!(view(solution,1, :))
    for t in 2:nt
        for x in 2:nx-1
            solution[t,x] = solution[t-1,x] + dt/τ*(λ^2*central_difference(solution[t-1,:], x, dx) - solution[t-1,x])
        end
        apply_neumann!(view(solution,t,:))
    end
    return solution
end 


#Each step is calculated with the A*x_n+1 = x_n which can than be used to solve for the next step with: x_n+1 = A^-1*x_n
function euler_implicit(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    pos_array = range(a,b, length=nx)
    V0 = exp.(-100 .* (pos_array .- 0.5).^2)
    dt = t_end/nt
    dx = (b-a)/(nx)
    α = 1 + dt/τ + λ^2*2*dt/(dx^2*τ)
    β = -λ^2/τ*dt/dx^2
    N = nx - 2
    @printf("α = %g, β = %g\n", α, β)

    diagonals = [β*ones(N-1),  α*ones(N),  β*ones(N-1)]
    A = spdiagm(-1 => diagonals[1], 0 => diagonals[2], 1 => diagonals[3])
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    apply_neumann!(view(solution,1, :))
    solution_array = Vector{Float64}(undef, nx)
    for t in 2:nt
        b_int = solution[t-1, 2:end-1]
        interior_solution = A \ b_int
        solution_array[2:end-1] = interior_solution
        apply_neumann!(solution_array)
        solution[t,:] = solution_array
    end
    return solution
end

#Each step is calculated with the A_1*x_n+1 = A_2*x_n which can than be used to solve for the next step with: x_n+1 = A_1^-1*A_2*x_n
function crank_nicolson(V0, nx::Int, nt::Int,t_end::Float64, a::Float64, b::Float64, λ = 1.0 , τ = 1.0 )
    dt = t_end/nt
    dx = (b-a)/(nx)
    α = 1/dt + 1/(2*τ) + λ^2/(τ*dx^2)
    β = -1/(2*τ)*λ^2/dx^2
    γ = 1/dt - 1/(2*τ) - λ^2/dx^2
    δ = -β
    N = nx-2
    diagonals_next = [β*ones(N-1),  α*ones(N),  β*ones(N-1)]
    digonals_previous = [δ*ones(N-1), γ*ones(N), δ*ones(N-1)]
    A_next = spdiagm(-1 => diagonals_next[1], 0 => diagonals_next[2], 1 => diagonals_next[3])
    A_previous = spdiagm( -1 => digonals_previous[1], 0 => digonals_previous[2], 1 => digonals_previous[3] )
    solution = zeros(Float64, nt, nx)
    solution[1, :] = V0
    apply_neumann!(view(solution,1, :))
    solution_array = Vector{Float64}(undef, nx)
    for t in 2:nt
        b_int = A_previous*solution[t-1, 2:end-1]
        interior_solution = A_next \ b_int
        solution_array[2:end-1] = interior_solution
        apply_neumann!(solution_array)
        solution[t,:] = solution_array
    end
    return solution
end

function analytical_solution(x, t, λ = 1.0 , τ = 1.0)
    V_0_tilde = 1
    return 1/sqrt(4*π*(λ^2/τ)*t)*exp(-(x-x_0)^2/(4*λ^2/τ*t)-t/τ)
end


a = 0.0
b = 1.0
nx = 100
nt = 20_000
time_stop = 1.0
x = range(a,b, length=nx)
t = range(0,time_stop, length = nt)
V0 = exp.(-10 .* (x .- 0.5).^2)
sol_explicit = euler_explicit(V0, nx, nt, time_stop, a, b)
sol_implicit = euler_implicit(V0, nx, nt, time_stop, a, b)
sol_crank = crank_nicolson(V0, nx, nt, time_stop, a, b)

plot(x, sol_explicit[1,:], lw=2, ls=:dash, label="Initial", xlabel="x", ylabel="V", title="Cable Equation: Final Profiles (Neumann BCs)")
plot!(x, sol_explicit[end,:], lw=2, label="Euler explicit")
plot!(x, sol_implicit[end,:], lw=2, label="Euler implicit")

plot(t,sol_explicit[:,1], xlabel = "time", ylabel = "V" ,title="Potential for first point over time", label= "Euler explicit")
plot!(t,sol_implicit[:,1], label="Euler implicit")
plot!(t,sol_crank[:,1], label = "Crank-Nicolson")

surface(x,t, sol_explicit)
surface(x,t, sol_implicit)
surface(x, t, sol_crank, xlabel="x", ylabel="t", zlabel="V",
        title="Surface plot of V(x, t)", camera = (80, 30))

