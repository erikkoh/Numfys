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
            push!(coefficient,BigFloat((n_log - n_i_log - i_log)))
        end
        coefficient[end] = 1.0
        coefficient[1] = 1.0
        coefficient = [BigFloat(exp(i)) for i in coefficient]
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
    
function convolution_for_prop_using_dist(prop,step_size)
    M = length(prop)-1
    setprecision(128)
    q_list = (range(0.0, 1.0-step_size, step=step_size))
    result = []
    for q in ProgressBar(q_list)
        binomial_dist = Binomial(M, q)
        pmf_values = pdf(binomial_dist, 0:M)  # Get PMF values for 0 to M
        conv_result = sum(pmf_values .* prop[1:length(pmf_values)])
        push!(result, conv_result)
    end
    return result
end

function convolution_for_property(prop,binom_coef,step_size)
    M = length(prop)-1
    q_list = BigFloat.(range(0.0, stop=1.0-step_size, step=step_size))
    property = []
    for q in ProgressBar(q_list)
        push!(property,sum([BigFloat(binom_coef[n+1]*q^n*(1-q)^(M-n)*prop[n+1]) for n in 0:M]))
    end
    return property
end

function convolution_for_property_log(prop, binom_coef)
    M = length(prop)-1
    q_list = BigFloat.(range(0.001, stop=1.0-0.001, step=0.001))
    property = [sum([(exp((log(binom_coef[n+1])+n*log(q)+(M-n)*log(1-q)+log(prop[n+1])))) for n in 0:M]) for q in ProgressBar(q_list)]
    return property
end

function convolution_for_property_optimized(prop,binom_coef,step_size)
    M = length(prop) - 1
    q_list = range(0.0, stop=1.0-step_size, step=step_size)
    
    binom_coef = BigFloat.(binom_coef)
    prop = BigFloat.(prop)
    
    property = Vector{BigFloat}(undef, length(q_list))

    Threads.@threads for i in ProgressBar(eachindex(q_list))
        q = q_list[i]
        q_pow = BigFloat(q) .^ (0:M)  # Precompute q^n for n = 0 to M
        one_minus_q_pow = reverse((1 - BigFloat(q)) .^ (0:M))  # Precompute (1-q)^(M-n) for n = 0 to M
        
        # Compute the property for this q
        property[i] = sum(binom_coef .* q_pow .* one_minus_q_pow .* prop)
    end

    return property
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
    binomial_coefficients = calculating_binomial_coefficients_striling(length(p_inf))
    println("Started convolution for p_inf")
    p_inf_convolution = convolution_for_prop_using_dist(p_inf,0.001)
    println("Started convolution for s")
    s_convolution = convolution_for_prop_using_dist(s,0.001)
    return Dict("p_inf"=> p_inf_convolution, "s"=> s_convolution)
end


function compare_time()
    JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_Grid_simulation_300"*".json")
    p_inf = JSON_info["p_inf"]
    s = JSON_info["s"]
    binomial_coefficient = calculating_binomial_coefficients_striling(length(p_inf))
    @time begin
        convolution_for_property_optimized(p_inf, binomial_coefficient)
    end
    @time begin
        convolution_for_property(p_inf, binomial_coefficient)
    end
end


function plotting_all_convoluted_values()
    
    Ns = [100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 900, 1000]
    p_inf_dic = Dict()
    s_dic = Dict()
    q_list = [i for i in range(0.0, 1.0-0.001, step=0.001)]

    for n in Ns 
        JSON_info = JSON.parsefile(find_folder("JSON_files") * "convolution_square_$n.json")
        p_inf_dic[n] = BigFloat.(JSON_info["p_inf"])
        s_dic[n] = BigFloat.(JSON_info["s"])
    end
    
    plot(title= "P_Inf", xlabel="q", ylabel="P_inf")

    for k in keys(p_inf_dic)
        p_list = p_inf_dic[k]
        plot!(q_list[300:600], p_list[300:600], label="N = $k")
    end
    savefig("all_convoluted")

end