using JSON
using Plots
using GLM
using Statistics
include("JSON_functions.jl")

import .JSONFunctions: find_folder

N_0 = 100
N_1 = 10000

raw_data_0 = JSON.parsefile(find_folder("JSON_files")*"convolution_p_inf_201.json")
raw_data_1 = JSON.parsefile(find_folder("JSON_files")*"convolution_p_inf_1001.json")

p_inf_0 = raw_data_0["p_inf"]
q_list_0 = raw_data_0["q"]
p_inf_1 = raw_data_1["p_inf"]
q_list_1 = raw_data_1["q"]

p_inf = [log(i) for i in p_inf_1]
q_list =[log(i) for i in q_list_1]

epsilon_0 = sqrt(N_0)
epsilon_0_list = [log(epsilon_0) for i in 1:length(p_inf_1)]

epsilon_1 = sqrt(N_1)
#Find what q the growth of p_inf to q is similar to beta/v*log(epsilon) wtf
# test_mesurment = findfirst(x, y-> x/y==)


ratio = p_inf./(log(epsilon_0))

plot(q_list_1, p_inf_1, label="p_inf", title="Overlayed Plots", xlabel="q", ylabel="Values", legend=:topright)
plot!(q_list_0, p_inf_0)
# plot!(q_list_0, epsilon_0_list, label="epsilon_0_list")