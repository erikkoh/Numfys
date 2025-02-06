using JSON
using Plots
using GLM
using DataFrames
using ProgressBars
using Statistics

include("JSON_functions.jl")

import .JSONFunctions: find_folder

Ns = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 900, 1000]
# Ns = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
# Ns = [100, 200, 250, 350, 450, 500, 600, 700, 900, 1000]
# Ns = [150, 250, 350, 450, 550]
p_inf_dic = Dict()
s_dic = Dict()

for n in Ns 
    JSON_info = JSON.parsefile(find_folder("JSON_files") * "convolution_square_$n.json")
    p_inf_dic[n] = BigFloat.(JSON_info["p_inf"])
    s_dic[n] = BigFloat.(JSON_info["s"])
end


epsilon_list = log.(Ns)

function extract_p_inf(q)
    p_inf_list = [(log((p[q]))) for p in values(p_inf_dic)]
    return p_inf_list
end 


function find_q_and_prop()
    q_list = [i for i in range(0.0,1.0-0.001,step=0.001)]
    slope = 0.0
    highest_cor = 0.0
    right_q = nothing
    for q in 1:1000

        y_data = Float64.([log(p[q]) for p in values(p_inf_dic)])
        x_data = log.(Ns)
        
        data = DataFrame(y= y_data, x = x_data)
        
        model = lm(@formula(y ~ x), data)
        
        coef_model = coef(model)[2] #extrating the coeficent for x as this would the be exponensial for epsilon
        
        current_cor = r2(model)

        if current_cor > highest_cor
            highest_cor = current_cor
            slope = coef_model
            right_q = q
        end
        
    end
    println("Hieghest cor: ", highest_cor)
    println("Right q value is: ", q_list[right_q])
    return slope, q_list[right_q]
end

function test_for_know_q(func)
    q = 500
    q_list = [i for i in range(0.0,1.0-0.001,step=0.001)]
    println("Test q: ", q_list[q])

    data = DataFrame(y= Float64.(func(q)), x = Float64.(epsilon_list))
        
    model = lm(@formula(y ~ x), data)
    
    coef_model = coef(model)[2]
    return coef_model
end

function find_gamma_nu(s_max_list)
    data = DataFrame(y = Float64.(log.(s_max_list)), x = Float64.(epsilon_list))
    model = lm(@formula(y ~ x), data)
    
    coef_model = coef(model)[2]
    return coef_model
end

function extract_max_s()
    s_list = [maximum(s_dic[n]) for n in Ns]
    # debug = findfirst(x-> x==s_list[1], s_dic[Ns[1]])
    # println(debug, s_dic[100][debug])
    return s_list
end

function find_nu()
    q_list = [i for i in range(0.0, 1.0-0.001, step=0.001)]
    s_list = [maximum(s_dic[n]) for n in Ns]
    max_q_index = []
    for i in eachindex(s_list)
        push!(max_q_index, findfirst(x-> x==s_list[i], s_dic[Ns[i]]))
    end
    println(max_q_index)


    q_max = find_q_and_prop()[2]
    q_q_max_list = [log(abs(q_list[q] - q_max)) for q in max_q_index if q !== nothing]
    
    data = DataFrame(y = Float64.(q_q_max_list), x = Float64.(epsilon_list))
    model = lm(@formula(y ~ x), data)

    coef_model = coef(model)[2]
    return coef_model
end

max_s = extract_max_s()

beta_nu = find_q_and_prop()[1]
test_beta_nu = test_for_know_q(extract_p_inf)

s_list = extract_max_s()

gamma_nu = find_gamma_nu(s_list)
inverse_nu = find_nu()

nu = (-1)/inverse_nu
println("and nu is:", nu)
println("Gamma is: ", gamma_nu*nu)
println("Beta is: ", beta_nu*(-nu))


# p_inf_list = extract_p_inf(right_q)


q_list = [i for i in range(0, 1.0-0.001, step= 0.001)]

# p_list = [BigFloat(log(i)) for i in p_inf_dic[100]]

# plot(epsilon_list, p_inf_list, label="p_inf", title="Overlayed Plots", xlabel="q", ylabel="Values", legend=:topright)

plot(q_list, p_inf_dic[700], label="epsilon_0_list")