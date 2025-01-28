using JSON
using ProgressBars
using Plots
using Combinatorics
include("JSON_functions.jl")

using .JSONFunctions: find_folder, write_to_JSON

# JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_Grid_simulation_100.json")
# JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_Grid_simulation_200.json")
# JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_Grid_simulation_300.json")



macro Name(args)
    string(args)
end



function calculating_binomial_coefficients_striling(number_of_bonds::Int64)
stirlings_approximation = x -> x * log10(x) - x + 0.5 * log10(2 * π * x)
    n_log = stirlings_approximation(number_of_bonds)
    coefficient = []
    println("Generating binomial coefficients:")
    for i  in ProgressBar(0:number_of_bonds-1)
        n_i_log = stirlings_approximation(number_of_bonds-i)
        i_log = stirlings_approximation(i)
        push!(coefficient,big(10^(n_log - n_i_log - i_log)))
    end
    coefficient[end] = 1.0
    println("Done!")
    return coefficient
end

function binomial_coefficient(n::Int, k::Int)
    if k > n
        return BigInt(0)
    end
    if k == 0 || k == n
        return BigInt(1)
    end
    k = min(k, n - k)  # Take advantage of symmetry
    coeff = BigInt(1)
    for i in 1:k
        coeff *= (n - i + 1)
        coeff //= i
    end
    return coeff
end

function calculating_binomial_coefficients(n::Int)
    coefficients = Vector{BigInt}(undef, n + 1)
    println("Generating binomial coefficients:")
    for k in ProgressBar(0:n)
        coefficients[k + 1] = binomial_coefficient(n, k)
    end
    println("Done!")
    return coefficients
end


function calculating_binomial_coefficients_log_2(number_of_bonds::Int64)
    println("Generating binomial coefficients:")
    binomial_coefficients = [binomial(BigInt(number_of_bonds),BigInt(i)) for i in ProgressBar(0:(number_of_bonds-1))]
    println("Done!")
return binomial_coefficients
end
    

function convolution_for_property(prop,binom_coef)
    M = length(prop)-1
    q_list = range(0, stop=1.0, step=0.001)
    # q_list = p_list
    property = [sum([ binom_coef[n+1]*q^n*(1-q)^(M-n)*prop[n+1] for n in 0:M]) for q in ProgressBar(q_list)]
    return property, q_list
end

function ploting_convolution_property(property)
    property_plot, q_plot = convolution_for_property(property, binomial_coefficient)
    plot(q_plot, property_plot, label = "$property", xlabel = "q", ylabel = "p_inf", title = "$property as a function of q")
    savefig(find_folder("Plots") * "convolution_"*"$(property)"*".png")
end

function save_data(prop, name::String )
    property, q_list = convolution_for_property(prop)
    value_dic = Dict("$name" => property, "q" => q_list)
    write_to_JSON(value_dic, "convolution_"*"$name"*"_"*"$(length(property))")
end


function convelute_data(N::Int)
    JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_Grid_simulation_$N"*".json")
    p_inf = JSON_info["p_inf"]
    s = JSON_info["s"]
    suseptibility = JSON_info["susept"]
    p_list = JSON_info["p"]
    println("Started convolution for p_inf")
    binomial_coefficients = calculating_binomial_coefficients_log_2(length(p_inf))
    println("Started convolution for p_inf")
    p_inf_convolution = convolution_for_property(p_inf, binomial_coefficients)
    println("Started convolution for s")
    s_convolution = convolution_for_property(s, binomial_coefficients)
    return Dict("p_inf"=> p_inf_convolution, "s"=> s_convolution)
end