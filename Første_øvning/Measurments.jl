using JSON
using Plots
using GLM
using DataFrames
using ProgressBars
using Statistics
include("JSON_functions.jl")

import .JSONFunctions: find_folder

Ns = [100,200,300,400,500,700,900,1000]
# Ns = [100, 200, 300, 400, 500, 700, 900, ]

p_inf_dic = Dict()
s_dic = Dict()
q_list = range(0.0,1.0-0.01,step=0.01)

for n in Ns 
    p_inf_dic[n] =JSON_info = JSON.parsefile(find_folder("JSON_files") * "convolution_square_$n.json")["p_inf"]
    p
end


epsilon_list = [log(i) for i in Ns]
println(length(p_inf_dic[100]))

function extract_p_inf(q)
    p_inf_list = [(log(BigFloat(p[q]))) for p in values(p_inf_dic)]
    return p_inf_list
end 

function find_q()
    q_list = [i for i in range(0.0,1.0-0.01,step=0.01)]
    highest_coef = 0
    right_q = nothing
    for q in ProgressBar(1:100)

        data = DataFrame(y= Float64.(extract_p_inf(q)), x = Float64.(epsilon_list))

    # Fit a linear model
        model = lm(@formula(y ~ x), data)

        coef_model = coef(model)[2]

    # Check if this is the lowest std
        if coef_model > highest_coef
            highest_coef = coef_model
            right_q = q
        end

    end
    println(right_q)
    println(highest_coef)
    println(q_list[right_q])
    return highest_coef, right_q
end

right_q = find_q()[2]

p_inf_list = extract_p_inf(right_q)

# p_list = [BigFloat(log(i)) for i in p_inf_dic[100]]

plot(epsilon_list, p_inf_list, label="p_inf", title="Overlayed Plots", xlabel="q", ylabel="Values", legend=:topright)

# plot!(q_list_0, epsilon_0_list, label="epsilon_0_list")