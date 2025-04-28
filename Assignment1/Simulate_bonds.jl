import JSON
using Random
using ProgressBars
using Base.Threads
include("JSON_functions.jl")

using  .JSONFunctions: find_folder, write_to_JSON

# rng = Random.MersenneTwister(123456789)


function swaping_bonds(bonds_list)

    for i in (1:(length(bonds_list)-1))
        candidate_number = rand(i+1:(length(bonds_list)))
        temp = bonds_list[i]
        bonds_list[i] = bonds_list[candidate_number]
        bonds_list[candidate_number] = temp
    end
    return bonds_list
end


function has_duplicates(list)
    seen = Set()
    for item in list
        if item in seen
            return true
        end
        push!(seen, item)
    end
    return false
end





function find_root_node(j, sites)
    if sites[j] < 0
        return j
    else
        sites[j] = find_root_node(sites[j],sites)
        return sites[j]
    end
end

function simulate_bonds(n::Int)
    JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_grid_bonds_$n.json")
    bonds = JSON_info["Bonds"]
    num_bonds =JSON_info["Number of bonds"]
    num_nodes = bonds[end][1]
    bonds = swaping_bonds(bonds)
    N = num_nodes
    sites = [-1 for i in 1:(N)]
    largest_cluster = Dict("index" => 1, "size" => 1)
    number_of_activated_bonds = 0
    p_0 = [0.0]
    p_inf = [1/N]
    p_inf_2 = [(1/N)^2]
    susept = [N*(p_inf_2[1]- p_inf[1]^2)^(1/2)]
    avarage_s = N
    s = [0.0]
    s_step = (avarage_s-(N*p_inf[end])^2)/(N*(1-p_inf[end]))
    steps = 0
    println("Start loop")
    for i in (1:num_bonds)
        node1, node2 = bonds[i]
        root1 = find_root_node(node1,sites)
        root2 = find_root_node(node2,sites)
        number_of_activated_bonds += 1
        if root1 != root2
            avarage_s = avarage_s -sites[root1]^2 - sites[root2]^2
            if sites[root1] < sites[root2]
                sites[root1] += sites[root2]
                sites[root2] = root1
                new_root = find_root_node(root1,sites)
            else
                sites[root2] += sites[root1]
                sites[root1] = root2
                new_root = find_root_node(root2,sites)
            end
            if sites[new_root] < largest_cluster["size"]
                largest_cluster = Dict("index" => new_root, "size" => sites[new_root])
            else
                s_step = (avarage_s - (N * p_inf[end])^2) / (N * (1 - p_inf[end]))
                s_step = isfinite(s_step) ? s_step : 0.0  # Ensure s_step is finite
            end
            avarage_s += sites[new_root]^2
        end
        push!(s,s_step)
        push!(p_0, number_of_activated_bonds/num_bonds)
        push!(p_inf,-largest_cluster["size"]/N)
        push!(p_inf_2,(largest_cluster["size"]/N)^2)
        push!(susept,N*(sum(p_inf_2)/length(p_inf_2) - (sum(p_inf)/length(p_inf))^2)^(1/2))
        steps += 1
    end
    println("Number of steps: ", steps)
    println("P_0: ", p_0[end])
    println("Number of bonds: ", num_bonds)
    println("Number of activated bonds: ", number_of_activated_bonds)
    println("Susceptiblitly: ", susept[end])
    println(length(s) == length(p_0)==length(p_inf)==length(p_inf_2)==length(susept))
    println("Started looking for largest cluster:")
    largest_cluster_list = [i for i in ProgressBar(eachindex(sites)) if find_root_node(i,sites) == largest_cluster["index"]]
    println("Done!")
    return largest_cluster, sites, largest_cluster_list, p_0, p_inf, p_inf_2, susept, s
end

function simulate_bonds_optimized(bonds, num_bonds, n::Int,iteration::Int, gird_type::String)
    N = Float64(n)
    sites = (fill(-1,Int64(N)))
    largest_cluster = Int64[1,-1]
    number_of_activated_bonds = Int64(0)
    p_0 = Vector{Float64}(undef, num_bonds)
    p_0[1] = 0
    p_inf = Vector{Float64}(undef, num_bonds)
    p_inf[1] = 1/N
    p_inf_2 = Vector{Float64}(undef, num_bonds)
    p_inf_2[1] = (1/N)^2
    susept = Vector{Float64}(undef, num_bonds)
    susept[1] = N*(p_inf_2[1]- p_inf[1]^2)^(1/2)
    avarage_s = Float64(N)
    s = Vector{Float64}(undef,num_bonds)
    s[1] = 0.0
    s_step = Float64(avarage_s-(N*p_inf[1])^2)/(N*(1-p_inf[1]))
    steps = Int64(0)
    sum_p_inf = Float64(0.0)
    sum_p_inf_2 = Float64(0.0)
    for i in (1:num_bonds)
        node1, node2 = bonds[i]
        root1 = find_root_node(node1,sites)
        root2 = find_root_node(node2,sites)
        number_of_activated_bonds += 1
        if root1 != root2
            avarage_s = avarage_s -sites[root1]^2 - sites[root2]^2
            if sites[root1] < sites[root2]
                sites[root1] += sites[root2]
                sites[root2] = root1
                new_root = root1
            else
                sites[root2] += sites[root1]
                sites[root1] = root2
                new_root = root2
            end
            if sites[new_root] < largest_cluster[2]
                largest_cluster[1] = new_root
                largest_cluster[2] = sites[new_root]
            else
                s_step = (avarage_s - (N * p_inf[i-1])^2) / (N * (1 - p_inf[i-1]))
                s_step = isfinite(s_step) ? s_step : 0.0
            end
            avarage_s += sites[new_root]^2
        end
        s[i] = s_step
        p_0[i] = number_of_activated_bonds/num_bonds
        p_inf[i] = abs(largest_cluster[2]/N)
        p_inf_2[i] = (largest_cluster[2]/N)^2
        sum_p_inf += p_inf[i]
        sum_p_inf_2 += p_inf_2[i]
        mean_p_inf = sum_p_inf / i
        mean_p_inf_2 = sum_p_inf_2 / i
        susept[i] = N * sqrt(abs(mean_p_inf_2 - mean_p_inf^2)) #The value of p_inf_2 is always larger og equal to p_inf ^2 but due to rounding errors abs is used to avoid complex roots
        steps += 1
    end
    # largest_cluster_list = [i for i in (eachindex(sites)) if find_root_node(i,sites) == largest_cluster[1]]
        largest_cluster_list = [1]
    return largest_cluster, sites, largest_cluster_list, p_0, p_inf, p_inf_2, susept, s, bonds
end

function testing_speed()
    @time begin
        simulate_bonds_optimized(300)
    end
    @time begin
        simulate_bonds(300)
    end
end


function p_images(p_list)
    p_dic = Dict{Float64,Any}()
    for p in p_list
        i = Int(floor(p*num_bonds))
        largest_cluster, sites, largest_cluster_list, p_0 = simulate_bonds(i)
        p_dic[p] = Dict("largest_cluster" => largest_cluster, "sites" => sites, "largest_cluster_list" => largest_cluster_list, "p_0" => p_0[end])
    end
    return p_dic
end


function avarage_values(iterations::Int64, n::Int64, grid_type::String)
    rng = Random.MersenneTwister(123456789)
    JSON_info = JSON.parsefile(find_folder("JSON_files") * "$(grid_type)_grid_bonds_$n.json")
    bonds = JSON_info["Bonds"]
    num_bonds =JSON_info["Number of bonds"]
    num_nodes = Int64(bonds[end][1])
    p_inf_list = zeros(num_bonds)
    p_inf_list = zeros(num_bonds)
    susept_list = zeros(num_bonds)
    s_list = zeros(num_bonds)
    p_list = zeros(num_bonds)
    p_inf_2_list = zeros(num_bonds)
    lock = Threads.SpinLock()
    pbar = ProgressBar(total = iterations)
    println("Starting simulation with $iterations iterations")
    Threads.@threads for i in (1:iterations)
        bonds = swaping_bonds(bonds)
        result = simulate_bonds_optimized(bonds, num_bonds, num_nodes, i, grid_type)
        p = result[4]
        p_inf = result[5]
        p_inf_2 = result[6]
        susept = result[7]
        s = result[8]
        Threads.lock(lock)
        update(pbar)
        p_inf_list = p_inf_list + p_inf/iterations
        susept_list = susept_list + susept/iterations
        s_list = s_list + s/iterations
        p_list = p_list + p/iterations
        p_inf_2_list = p_inf_2_list + p_inf_2/iterations
        Threads.unlock(lock)
        # push!(p_list,p)
        # push!(p_inf_list,p_inf)
        # push!(p_inf_2_list,p_inf_2)
        # push!(susept_list,susept)
        # push!(s_list,s)
    end
    println("Done!")
    println("started avraging values")
    # p_inf_avarage = [sum([p[i] for p in p_inf_list ])/length(p_inf_list) for i in eachindex(p_inf_list[1])]
    # susept_avarage = [sum([p[i] for p in susept_list ])/length(susept_list) for i in eachindex(susept_list[1])]
    # s_avarage = [sum([p[i] for p in s_list ])/length(s_list) for i in eachindex(s_list[1])]
    # p_avarage = [sum([p[i] for p in p_list])/length(p_list) for i in eachindex(p_list[1])]
    # p_inf_2_avarage = [sum([p[i] for p in p_inf_2_list ])/length(p_inf_2_list) for i in eachindex(p_inf_2_list[1])]
    values_dic = Dict("p_inf" => p_inf_list,"p_inf_2" => p_inf_2_list, "susept" => susept_list, "s" => s_list, "p" => p_list)
    println("Done!")
    return values_dic
end
