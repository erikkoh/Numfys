using JSON
using ProgressBars
using Plots
using Combinatorics
using LinearAlgebra
using Base.Threads
using Distributions
include("JSON_functions.jl")

using .JSONFunctions: find_folder, write_to_JSON


macro Name(args)
    string(args)
end


function calculating_binomial_coefficients_striling(number_of_bonds::Int64)
    stirlings_approximation = x -> x * log(x) - x + 0.5 * log(2 * π * x)
        n_log = stirlings_approximation(number_of_bonds)
        coefficient = []
        println("Generating binomial coefficients:")
        for i  in ProgressBar(0:number_of_bonds-1)
            n_i_log = stirlings_approximation(number_of_bonds-i)
            i_log = stirlings_approximation(i)
            push!(coefficient,((n_log - n_i_log - i_log)))
        end
        coefficient[end] = 1.0
        coefficient[1] = 1.0
        println("Done!")
        return coefficient
end

function calculating_binomial_coefficients(number_of_bonds::Int64)
    M = number_of_bonds
    coefficients = Vector{Float64}(undef, M)
    coefficients[1] = 0.0
    for n in 1:(M-1)
        coefficients[n+1] = coefficients[n] + log(M - n + 1) - log(n)
    end
    return coefficients
end


    
function convolution_for_prop_using_dist(prop)
    M = length(prop)
    q_list = 0.001:0.0001:1.0-0.001  
    result = Vector{Float64}(undef, length(q_list))
    for (i, q) in ProgressBar(enumerate(q_list))
        binomial_dist = Binomial(M, q)
        pmf_values = pdf(binomial_dist, 0:(M-1))  # Get PMF values for 0 to M
        conv_result = sum(pmf_values .* prop)
        result[i] =  conv_result
    end
    return result
end


function convolution_for_property_log(prop, binom_coef)
    M = length(prop)
    q_list = 0.001:0.0001:1.0-0.001
    property = Vector{Float64}(undef, length(q_list))
    
    log_q = log.(q_list)
    log_1_q = log.(1 .- q_list)
    
    ns = 0:(M-1)
    
    for (i, q) in ProgressBar(enumerate(q_list))
        binom = (binom_coef) .+ (ns .* log_q[i]) .+ ((M .-ns) .* log_1_q[i])
        binom = exp.(binom)  
        property[i] = sum(binom .* prop) 
    end
    return property
end




function convelute_data(N::Int,grid_stype::String)
    JSON_info = JSON.parsefile(find_folder("JSON_files") * "$(grid_stype)_Grid_simulation_$N"*".json")
    p_inf = Float64.(JSON_info["p_inf"])
    s = Float64.(JSON_info["s"])
    suseptibility = JSON_info["susept"]
    p_list = JSON_info["p"]
    binomial_coefficients = calculating_binomial_coefficients(length(p_inf))
    println("Started convolution for p_inf")
    p_inf_convolution = convolution_for_property_log(p_inf, binomial_coefficients)
    # p_inf_convolution = convolution_for_prop_using_dist(p_inf)
    println("Started convolution for s")
    s_convolution = convolution_for_property_log(s, binomial_coefficients)
    println("Started convolution for susceptiblitly")
    susept_convolution = convolution_for_property_log(suseptibility, binomial_coefficients)
    # s_convolution = convolution_for_prop_using_dist(s)
    return Dict("p_inf"=> p_inf_convolution, "s"=> s_convolution, "suscept" => susept_convolution)
end



function plotting_all_convoluted_values()
    
    Ns = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    Ns = [100, 200, 300, 400]
    p_inf_dic = Dict()
    s_dic = Dict()
    q_list = 0.001:0.0001:1.0-0.001
    # q_list = 0.001:0.001:1.0-0.01  

    for n in Ns 
        JSON_info = JSON.parsefile(find_folder("JSON_files") * "convolution_square_$(n).json")
        p_inf_dic[n] = Float64.(JSON_info["p_inf"])
        s_dic[n] = Float64.(JSON_info["s"])
    end
    
    plot(title= "P_Inf", xlabel="q", ylabel="P_inf")

    for k in keys(p_inf_dic)
        p_list = p_inf_dic[k]
        plot!(q_list, p_list, label="N = $k")
    end
    savefig(find_folder("Convolution_plots")*"all_convoluted")

end
