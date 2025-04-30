using Plots
using Random
using ProgressBars
using Statistics
using StatsBase
using GLM                
using DataFrames         

Random.seed!(1234242)  # Set a random seed for reproducibility
mkpath("./plots")  # Create a directory for plots if it doesn't exist

"Generates a random system of spins with dimensions x_dim and y_dim."
function generate_system(x_dim, y_dim)
    system = zeros(x_dim, y_dim)
    possible_spins = [-1, 1]
    for i in eachindex(system)
        system[i] = rand(possible_spins)
    end
    return system
end



"Calculates the energy difference when flipping a spin at (i, j) in the system."
function calc_dE(system, i, j, H)
    x_dim, y_dim = size(system)  
    # Calculate the energy difference when flipping the spin at (i, j)
    i_up = mod1(i-1, x_dim)
    i_down = mod1(i+1, x_dim)
    j_left = mod1(j-1, y_dim)
    j_right = mod1(j+1, y_dim)
    ΔE = 2 * system[i, j] * (system[i_up, j] + system[i_down, j] + system[i, j_left] + system[i, j_right]) + 2* H *system[i,j]
    return ΔE
end

"Calculates the total energy of the system, by calculating the energy of each spin and its neighbors."
function get_energy(system, H=0.0)
    # Calculate the total energy of the system
    x_dim, y_dim = size(system)  
    energy = 0.0
    for i in 1:x_dim
        for j in 1:y_dim
            energy += -system[i, j] * (system[mod1(i-1, x_dim), j] + system[mod1(i+1, x_dim), j] + system[i, mod1(j-1, y_dim)] + system[i, mod1(j+1, y_dim)])
        end
    end
    return energy / 2.0 - H*sum(system)  # Each pair is counted twice
end


"Performs the Mc sweep of the system, flipping spins with a probability determined by the energy difference and temperature."
function sweep(system, num_sweeps, T, possible_indexes, H = 0.0)
    # Perform a sweep of the system
    magnetization = zeros(num_sweeps)  # Initialize magnetization with the first value
    magnetization[1] = sum(system)  # Initial magnetization
    energy = zeros(num_sweeps)  # Initialize energy with the first value
    energy[1] = get_energy(system)  # Initial energy

    #built in shuffle function in julia to ensure all spins are considered
    for m in 2:num_sweeps+1
        possible_indexes = shuffle!(possible_indexes)
        total_deltaE = 0.0
        for s in eachindex(possible_indexes)
            i, j = possible_indexes[s]
            ΔE = calc_dE(system, i, j, H)
            r = rand()  # Generate a random number between 0 and 1
            if r < exp(-ΔE / T)'
                total_deltaE += ΔE  # Accumulate the energy difference
                system[i, j] *= -1  # Flip the spin
            end
        end
        energy[mod1(m,num_sweeps)] = energy[m-1] + total_deltaE  # Update energy after each sweep    
        magnetization[mod1(m,num_sweeps)] = sum(system)  # Update magnetization after each sweep
    end
    return system, magnetization, energy
end

"Simple function to visualize the ising model"
function plot_system(system, ttl)
    # Plot the system using a heatmap
    heatmap(system, c=:grays, title = ttl, xlabel="x", ylabel="y", aspect_ratio=1)
end


function ploting_magnetization()
    temps = [2.1, 2.3, 2.5]
    x_dim, y_dim = 40, 40
    num_sweeps = 5_000
    initial_system = generate_system(x_dim, y_dim)
    p1 = plot_system(initial_system, "")
    savefig(p1, "./plots/Initial_System.png")
    possible_indeces = [(i, j) for i in 1:x_dim for j in 1:y_dim]
    system_1, magnetization_1, _ = sweep(initial_system, num_sweeps, temps[1], possible_indeces)
    p2 = plot(magnetization_1, xlabel="Sweep", ylabel="Magnetization", label="T=$(temps[1])")
    p3 = plot_system(system_1,"T=$(temps[1])")
    savefig(p3, "./plots/Final_System$(temps[1]).png")
    for ts in temps[2:end]
        system_1, magnetization_1, _ = sweep(initial_system, num_sweeps, ts, possible_indeces)
        p3 = plot_system(system_1, "T=$ts")
        savefig(p3, "./plots/Final_System$ts.png")
        plot!(p2, magnetization_1,legend=:topright, label="T=$ts")
    end
    savefig(p2,"./plots/Magnetization_vs_Sweep.png")
end


"Runns sweep function for different temperatures and calculates the energy, magnetization, heat capacity, and critical temperature"
function tempreture_trend_for_energy_and_mag(;L::Int=40, num_sweeps::Int=5_000, temps = range(1.0, 3.0, length= 100),possible_indexes = [])
    x_dim, y_dim = L, L
    initial_system = generate_system(x_dim, y_dim)
    energy = zeros(length(temps))
    magnetization = zeros(length(temps))
    variance_dE = zeros(length(temps))
    heat_capacity = zeros(length(temps))
    cut_off = 1000  # Cutoff for the first 1000 sweeps
    if possible_indexes == []  # If no possible indexes are provided, generate them
        possible_indexes = [(i, j) for i in 1:x_dim for j in 1:y_dim]
    end
    for (i, ts) in ProgressBar(enumerate(temps))
        system_1, magnetization_1, energies = sweep(initial_system, num_sweeps, ts, possible_indexes)
        energy[i] = energies[end]  # Final energy after all sweeps
        magnetization[i] = abs(magnetization_1[end])
        variance_dE[i] = mean(energies[cut_off:end].^2) - (mean(energies[cut_off:end]))^2  # Calculate variance of energy
        heat_capacity[i] = variance_dE[i] / (ts^2)  # Calculate heat capacity
    end
    argmax_heat_capacity = argmax(heat_capacity[2:end])  # Find the index of the maximum heat capacity
    #find index of argmax_heat_capacity in temps
    critical_temp = temps[argmax_heat_capacity]  # Find the temperature corresponding to the maximum heat capacity
    return temps, energy, magnetization, heat_capacity, critical_temp
end

function ploting_heat_mag_and_energy()
    temps, energy, magnetization, heat_capacity, critical_temp = tempreture_trend_for_energy_and_mag()
    println("Critical Temperature: ", critical_temp)
    scatter(temps, energy, xlabel="Temperature", ylabel="Energy", label="Energy")
    scatter!(legend=:topright)
    savefig("./plots/Energy_vs_Temperature.png")
    scatter(temps, magnetization, xlabel="Temperature", ylabel="Magnetization", label="Magnetization")
    scatter!(legend=:topright)
    savefig("./plots/Magnetization_vs_Temperature.png")
    scatter(temps, heat_capacity, xlabel="Temperature", ylabel="Heat Capacity", label="Heat Capacity")
    scatter!(legend=:topleft)
    savefig("./plots/Heat_Capacity_vs_Temperature.png")
    temps, _, magnetization,_,_ = tempreture_trend_for_energy_and_mag(temps = range(0.0, 1.0, length=100))
    scatter(temps, magnetization, xlabel="Temperature", ylabel="Magnetization", label="Magnetization")
    savefig("./plots/Appendicmagnetisation")
end

"Runs the sweep function for different L values and calculates the critical temperature."
function find_L_dependence_of_alpha(Ls::Vector{Int} = [10, 20, 40, 80], num_sweeps::Int=5_000)
    critical_temp_L = zeros(length(Ls))
    heat_capacity_L = Vector{Vector{Float64}}(undef, length(Ls))
    for (i, L) in ProgressBar(enumerate(Ls))
        _,_,_, heat_capacity, critical_temp = tempreture_trend_for_energy_and_mag(L=L, num_sweeps=num_sweeps)
        critical_temp_L[i] = critical_temp
        heat_capacity_L[i] = heat_capacity
    end
    return critical_temp_L, heat_capacity_L
end

function plot_L_dependence_of_alpha(num_sweeps::Int=10_000)
    temps = range(1.0, 3.0, length=100)
    Ls = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
    critical_temp_L, heat_capacity_L = find_L_dependence_of_alpha(Ls, num_sweeps)
    p1 = plot(temps[10:end], heat_capacity_L[1][10:end], xlabel="Temperature", ylabel="Heat Capacity", label="L=$(Ls[1])")
    p2 = plot(Ls, critical_temp_L, label="T_c", xlabel="L", ylabel="Temprature")
    for i in 2:length(Ls)
        plot!(p1, temps[10:end], heat_capacity_L[i][10:end], label="L=$(Ls[i])")
    end
    plot!(p1, legend=:topright)
    plot!(p2, legend=:topright)
    scatter!(p3, legend=:topright)
    savefig(p1, "./plots/Heat_Capacity_vs_Temperature_L.png")
    savefig(p2, "./plots/critical_t_L.png")
    savefig(p3, "./plots/Scaling_of_max_C")
end

"Modified sweep function to include constant points."
function ising_model_with_constant_points(x_dim, y_dim, num_sweeps, p = 0.01)
    num_constant_points = Int(round(p * x_dim * y_dim))  # Number of constant points
    possible_indeces = [(i, j) for i in 1:x_dim for j in 1:y_dim]
    constant_points_indexes = sample(1:length(possible_indeces), num_constant_points, replace=false)  # Randomly select constant points
    deleteat!(possible_indeces, sort(constant_points_indexes))  # Remove constant points from possible indexes

    (temps, energy, magnetization, heat_capacity, critical_temp) = tempreture_trend_for_energy_and_mag(L=40, num_sweeps=num_sweeps,possible_indexes = possible_indeces)
    return temps, magnetization, energy, heat_capacity
end

function plot_different_constant_points()
    x_dim, y_dim = 40, 40
    num_sweeps = 5_000
    p_values = [0.01, 0.05, 0.1]  # Different probabilities for constant points
    (temps, energy, magnetization, heat_capacity, critical_temp) = tempreture_trend_for_energy_and_mag(L=x_dim, num_sweeps=num_sweeps)

    p1 = scatter(temps, magnetization, xlabel="Temprature", ylabel="Magnetization", label="p=0.0")
    p2 = scatter(temps, energy, xlabel="Temprature", ylabel="Energy", label="p=0.0")
    p3 = scatter(heat_capacity, xlabel="Temprature", ylabel="Heat Capacity", label="p=0.0")
    for p in p_values
        temps, magnetization, energy, heat_capacity = ising_model_with_constant_points(x_dim, y_dim, num_sweeps, p)
        scatter!(p1, temps, magnetization, label="p=$p")
        scatter!(p2, temps, energy, label="p=$p")
        scatter!(p3, heat_capacity, label="p=$p")
        
    end
    scatter!(p1, legend=:topright)
    scatter!(p2, legend=:topright)
    scatter!(p3, legend=:topright)

    savefig(p1, "./plots/Magnetization_vs_temprature_with_p.png")
    savefig(p2, "./plots/Energy_vs_temprature_with_p.png")
    savefig(p3, "./plots/Heat_Capacity_vs_temprature_with_p.png")
end


"simulated_anehiling_with_constant_points function to run the sweep function with constant points."
function simulated_anehiling_with_constant_points(system, possible_indeces, constant_points_indexes, iterations, num_sweeps, p = 0.01)
    temps = range(0.5, 3.0, length= 1000)
    temps = reverse(temps)  # Reverse the order of temperatures
    energies = zeros()
    magnetization = zeros(length(temps))  # Initialize magnetization with the first value
    for T in temps
        system, magnetization, energies = sweep(system, num_sweeps, T, possible_indeces)
    end
    system[constant_points_indexes] .*= 2
    heatmap(system, c=:grays, title = "$p impureties and $iterations attempt", aspect_ratio=1)
    mkpath("./plots/SA$(p)")
    savefig("./plots/SA$(p)/SA_with_constant_points_$p _$iterations.png")
    return magnetization[end], energies[end]
end



# Todo make the initial condition the same for all iterations
function simulated_anehiling_with_constant_points_all(iterations::Int=3)
    p_values = [0.01, 0.05, 0.1]  # Different probabilities for constant points
    x_dim, y_dim = 40, 40
    for p in p_values
        num_sweeps = 1
        magnetization = zeros(iterations)
        energies = zeros(iterations)
        num_constant_points = Int(round(p * x_dim * y_dim))  # Number of constant points
        for i in 1:iterations
            possible_indeces = [(i, j) for i in 1:x_dim for j in 1:y_dim]
            constant_points_indexes = sample(1:length(possible_indeces), num_constant_points, replace=false)  # Randomly select constant points
            deleteat!(possible_indeces, sort(constant_points_indexes))  # Remove constant points from possible indexes
            system_copy = generate_system(x_dim, y_dim)  # Create a copy of the initial system for each iteration
            magnetization[i], energies[i] = simulated_anehiling_with_constant_points(system_copy, possible_indeces, constant_points_indexes, i,num_sweeps, p)
            num_sweeps *= 10
        end
        println("p = $p: Final Magnetization: ", magnetization, " Final Energy: ", energies)
    end
end



#This marks the begining of the second part of the main assignment


"Ploting fuction"
function plot_with_h_field(T, H)
    x_dim, y_dim = 40, 40
    num_sweeps = 5_000
    initial_system = generate_system(x_dim, y_dim)
    possible_indeces = [(i, j) for i in 1:x_dim for j in 1:y_dim]
    system_1, magnetization_1, _ = sweep(initial_system, num_sweeps, T, possible_indeces, H)
    plot(magnetization_1, title="Magnetization vs. Sweep at ", xlabel="Sweep", ylabel="Magnetization", label="T=$T and H=$H")
    plot!(legend=:bottomleft, title="Magnetization vs. Sweep at T=$T and H=$H")
    savefig("./plots/Magnetization_vs_Sweep_$(T)_$H.png")
    heatmap(system_1, c=:grays, title="Final System after $num_sweeps sweeps for T=$T and H=$H", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig("./plots/H_field_effect_at_$(T)_$H.png")
end


function calc_suseptibility(H, temps = range(1.0, 4.0, length=200),num_sweeps::Int=5_000)
    x_dim, y_dim = 40, 40
    initial_system = generate_system(x_dim, y_dim)
    cut_off = 1000  # Cutoff for the first 1000 sweeps
    susceptibility = zeros(length(temps))
    magnetization_final = zeros(length(temps))  # Initialize magnetization with the first value
    for (T, i) in ProgressBar(zip(temps, 1:length(temps)))
        possible_indeces = [(i, j) for i in 1:x_dim for j in 1:y_dim]
        _, magnetization, _ = sweep(initial_system, num_sweeps, T, possible_indeces, H)
        susceptibility[i] = (mean(magnetization[cut_off:end].^2) - (mean(magnetization[cut_off:end]))^2)/T  # Calculate variance of energy
        magnetization_final[i] = abs(magnetization[end])
    end
    return temps, magnetization_final, susceptibility
end


function plot_susceptibility_for_different_H()
    Hs = [0.01,0.02,0.03, 0.05]  # Different probabilities for constant points
    temps, magnetization, susceptibility = calc_suseptibility(Hs[1])
    p1 = scatter(temps, susceptibility, xlabel="Temperature", ylabel="Susceptibility", label="H=$(Hs[1])")
    p2 = scatter(temps, magnetization, xlabel="Temperature", ylabel="Magnetization", label="H=$(Hs[1])")
    for H in Hs[2:end-1]
        temps, magnetization, susceptibility = calc_suseptibility(H)
        scatter!(p1, temps, susceptibility, label="H=$H)")
        scatter!(p2, temps, magnetization, label="H=$H")
    end
    temps, magnetization, susceptibility = calc_suseptibility(Hs[end])
    scatter!(p1, temps, susceptibility, label="H=$(Hs[end])")
    scatter!(p2, temps, magnetization, label="H=$(Hs[end])")
    p3 = scatter(temps, magnetization, xlabel="Temprature", ylabel = "Magnetization", label="M", title="H=$(Hs[end])")
    scatter!(p3, temps, magnetization*Hs[end], label="X*H")


    scatter!(p1, legend=:topright)
    scatter!(p2, legend=:topright)
    scatter!(p3, legend=:topright)
    savefig(p1, "./plots/Susceptibility_vs_Temperature.png")
    savefig(p2, "./plots/Magnetization_vs_Temperature.png")
    savefig(p3, "./plots/Linearmagnetic.png")
end

function scaling_of_susceptibility(Hs = [0.01, 0.02 , 0.03, 0.05])
    critical_temp = 2.3
    gamma = 7/4
    beta = 1/8
    delta = 15
    temps = range(2.0,2.6, length=100)
    mkpath("./plots/Susceptibility")
    for (i,H) in enumerate(Hs)
        temps, _, susceptibility = calc_suseptibility(H, temps)
        t_h = abs.((temps .- critical_temp)) ./ critical_temp
        susceptibility_scaled = susceptibility .* t_h .^ gamma
        h_scaled = H .* t_h .^ (-beta * delta)
        
        # Split data into before and after critical_temp
        before_critical = findall(t -> t < critical_temp, temps)
        after_critical = findall(t -> t > critical_temp, temps)
        
        scatter((-1).*(h_scaled[before_critical]), (susceptibility_scaled[before_critical]),
            label="T < Critical Temp",
            xlabel="H scaled", ylabel="Susceptibility scaled", color=:blue)

        scatter!((h_scaled[after_critical]), (susceptibility_scaled[after_critical]),
            label="T > Critical Temp",
            xlabel="H scaled", ylabel="Susceptibility scaled", color=:red, title="$H")
            scatter!(legend=:topright)
        savefig("./plots/Susceptibility/scaling_of_susceptibility_after$i.png")

        # Log-log plot
        scatter(log.(h_scaled[before_critical]), log.(susceptibility_scaled[before_critical]),
            label="T < Critical Temp",
            xlabel="log(H scaled)", ylabel="log(Susceptibility scaled)", color=:blue)

        scatter!(log.(h_scaled[after_critical]), log.(susceptibility_scaled[after_critical]),
            label="T > Critical Temp",
            xlabel="log(H scaled)", ylabel="log(Susceptibility scaled)", color=:red, title="$H")
        savefig("./plots/Susceptibility/loglog_scaling_of_susceptibility_$i.png")
    end
end
