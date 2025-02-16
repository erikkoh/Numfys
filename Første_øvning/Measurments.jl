using JSON
using Plots
using GLM
using DataFrames
using ProgressBars
using Statistics

include("JSON_functions.jl")

import .JSONFunctions: find_folder, write_to_JSON


function find_critical_exponents(grid_type::String)

    Ns = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    # Ns = [100, 200, 300, 400]
    p_inf_dic = Dict()
    s_dic = Dict()

    for n in Ns 
        JSON_info = JSON.parsefile(find_folder("JSON_files") * "convolution_$(grid_type)_$n.json")
        p_inf_dic[n] = round.((JSON_info["p_inf"]),digits=2)
        s_dic[n] = (JSON_info["s"])
    end


    epsilon_list = log.(Ns)

    function extract_p_inf(q)
        p_inf_list = [(log((p[q]))) for p in values(p_inf_dic)]
        return p_inf_list
    end 


    function find_q_and_prop()
        q_list = 0.001:0.0001:1.0-0.001
        # q_list = 0.001:0.001:1.0-0.01  
        slope = 0.0
        highest_cor = 0.0
        right_q = NaN
        r2_list = []
        for (i,q) in enumerate(q_list)

            y_data = Float64.([log(p[i]) for p in values(p_inf_dic)])
            x_data = log.(Ns)
            
            data = DataFrame(y= y_data, x = x_data)
            
            model = lm(@formula(y ~ x), data)
            
            coef_model = coef(model)[2] #extrating the coeficent for x as this would the be exponensial for epsilon
            
            y_pred = predict(model)
        
            SS_res = sum((y_data .- y_pred).^2)
            SS_tot = sum((y_data .- mean(y_data)).^2)
            
            current_cor = 1 - (SS_res / SS_tot)

            current_cor = r2(model)
            # current_cor = if current_cor<0 0.0 else current_cor end
            # current_cor = cor(data.y, data.x)

            push!(r2_list, current_cor)
        
            if abs(current_cor) > abs(highest_cor)
                highest_cor = current_cor
                slope = coef_model
                right_q = q
            end
            
        end
        plot!(q_list, r2_list)
        savefig("q_r2_plot")
        println("Hieghest cor: ", highest_cor)
        println("Right q value is: ", right_q)
        return slope, right_q
    end

    function test_for_know_q(func)
        q = 5000-9
        # q = 500
        q_list = 0.001:0.0001:1.0-0.001
        # q_list = 0.001:0.001:1.0-0.01  
        println("Test q: ", q_list[q])
        y_values = Float64.(func(q))
        x_values = Float64.(log.(Ns)) 
        data = DataFrame(y= y_values, x = x_values)
            
        model = lm(@formula(y ~ x), data)
        r2_score = r2(model)
        
        println("Test corelation: ",r2_score)
        coef_model = coef(model)[2]
        scatter(data.x, data.y)
        plot!(data.x, predict(model))
        savefig("linarregress2")
        return coef_model
    end

    function find_gamma_nu(s_max_list)
        data = DataFrame(y = Float64.(log.(s_max_list)), x = Float64.(epsilon_list))
        model = lm(@formula(y ~ x), data)
        
        scatter(data.x, data.y)
        plot!(data.x, predict(model), xlabel="log(N)/2", ylabel="log(Smax)", title= "Smax")
        savefig("linarregress_s_max")
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
        q_list = 0.001:0.0001:1.0-0.001
        # q_list = 0.001:0.001:1.0-0.01  
        s_list = [maximum(s_dic[n]) for n in Ns]
        max_q_index = []
        for i in eachindex(s_list)
            push!(max_q_index, findfirst(x-> x==s_list[i], s_dic[Ns[i]]))
        end


        # q_max = find_q_and_prop()[2]
        q_max=0.5
        q_q_max_list = [log(abs(q_list[q] - q_max)) for q in max_q_index if q !== nothing]
        println(q_q_max_list)
        
        data = DataFrame(y = Float64.(q_q_max_list), x = Float64.(epsilon_list))
        model = lm(@formula(y ~ x), data)

        scatter(data.x, data.y)
        plot!(data.x, predict(model))
        savefig("linarregress_q_s_max")
        coef_model = coef(model)[2]
        return coef_model
    end

    max_s = extract_max_s()

    beta_nu, critical_prob = find_q_and_prop()
    test_beta_nu = test_for_know_q(extract_p_inf)

    s_list = extract_max_s()

    gamma_nu = find_gamma_nu(s_list)
    inverse_nu = find_nu()

    nu = (-1)/inverse_nu
    gamma = gamma_nu*nu
    beta = beta_nu*(-nu)
    println("and nu is:", nu)
    println("Gamma is: ", gamma)
    println("Beta is: ", beta)
    println("Test beta is: " ,test_beta_nu*(-nu))

    results_dic = Dict("P_crit" => critical_prob, "Nu" => nu, "Gamma" => gamma, "Beta" => beta)
    write_to_JSON(results_dic, "Results_$(grid_type)_grid")

end

find_critical_exponents("Square")


JSON_info = JSON.parsefile(find_folder("JSON_files") * "Results_Square_grid.json")
println(JSON_info)