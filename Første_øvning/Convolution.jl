using JSON
using ProgressBars
using Plots
include("JSON_functions.jl")

using .JSONFunctions: find_folder

JSON_info = JSON.parsefile(find_folder("JSON_files") * "avarage_values_20000.json")


p_inf = JSON_info["p_inf"]
s = JSON_info["s"]
suseptibility = JSON_info["susept"]
p_list = JSON_info["p"]

println(length(p_inf))

stirlings_approximation = x -> x * log(x) - x + 0.5 * log(2 * π * x)
function calculating_binomial_coefficients_log_1(number_of_bonds::Int64)
    n_log = stirlings_approximation(number_of_bonds)
    coefficient = []
    println("Generating binomial coefficients:")
    for i  in ProgressBar(0:number_of_bonds-1)
        n_i_log = stirlings_approximation(number_of_bonds-i)
        i_log = stirlings_approximation(i)
        push!(coefficient,exp(n_log - n_i_log - i_log))
    end
    coefficient[end] = 1.0
    println("Done!")
    return coefficient
end


function calculating_binomial_coefficients_log_2(number_of_bonds::Int64)
    log_coefficient = [0.0]
    for i in 1:number_of_bonds
        next_log_coeff = log_coefficient[end] + log(number_of_bonds - i + 1) - log(i)
        push!(log_coefficient, next_log_coeff)
    end
    coefficient = exp.(log_coefficient)
    problem = findfirst(x -> x == Inf, coefficient)
    println(problem)
    println(stirlings_approximation(problem))
    return coefficient
end

    
binomial_coefficients = calculating_binomial_coefficients_log_1(length(p_inf))
println(binomial_coefficients)

function convolution_for_property(prop)
    M = length(prop)-1
    q_list = range(0, stop=1.0, step=0.001)
    property = [sum([ binomial_coefficients[n+1]*q^n*(1-q)^(M-n)*prop[n+1] for n in 0:M]) for q in ProgressBar(q_list)]
    return property, q_list
end
p_inf_plot, q_plot = convolution_for_property(p_inf)


plot(q_plot, p_inf_plot, label = "p_inf", xlabel = "q", ylabel = "p_inf", title = "p_inf as a function of q")

