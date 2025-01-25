import JSON
import Random
using ProgressBars
include("JSON_functions.jl")

using  .JSONFunctions: find_folder

Random.seed!(144)

JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_grid.json")

bonds = JSON_info["Bonds"]

num_bonds =JSON_info["Number of bonds"]

function swaping_bonds(bonds_list)
    for i in 1:(eachindex(bonds_list)-1)
        candidate_number = rand(i+1:num_bonds)
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

function simulate_bonds(p)
    max_iter = 2_000_000
    pgb = ProgressBar(total=max_iter)
    N = bonds[end][1]
    sites = [-1 for i in 1:(N+1)]
    largest_cluster = Dict("index" => 1, "size" => 1)
    number_of_activated_bonds = 0
    p_0 = [0.0]
    p_inf = [1/N]
    p_inf_2 = [(1/N)^2]
    susept = [N*(p_inf_2[1]- p_inf[1]^2)^(1/2)]
    avarage_s = N
    s = [0.0]
    steps = 0
    println("Start loop")
    while p_0[end] <= p && steps < max_iter && largest_cluster["size"] < N
        # lowest_value = minimum(sites) this is way to ineficent for 1000^2 this would take a stupid amount of time
        
        selected_bond = rand(1:num_bonds)
        node1, node2 = bonds[selected_bond]
        root1 = find_root_node(node1,sites)
        root2 = find_root_node(node2,sites)
        if root1 != root2
            avarage_s = avarage_s -sites[root1]^2 - sites[root2]^2
            update(pgb)
            number_of_activated_bonds += 1
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
            end
            avarage_s += sites[new_root]^2
            s_step = (avarage_s-(N*p_inf[end])^2)/(N*(1-p_inf[end])) 
            push!(s,s_step)
            push!(p_0, number_of_activated_bonds/num_bonds)
            push!(p_inf,-largest_cluster["size"]/N)
            push!(p_inf_2,(largest_cluster["size"]/N)^2)
            push!(susept,N*(sum(p_inf_2)/length(p_inf_2) - (sum(p_inf)/length(p_inf))^2)^(1/2))
        else
        end
        steps += 1
    end
    println("Number of steps: ", steps)
    println("P_0: ", p_0[end])
    println("Number of bonds: ", num_bonds)
    println("Number of activated bonds: ", number_of_activated_bonds)
    println("Susceptiblitly: ", susept[end])
    println("Started looking for largest cluster:")
    largest_cluster_list = [i for i in ProgressBar(eachindex(sites)) if find_root_node(i,sites) == largest_cluster["index"]]
    println("Done!")
    return largest_cluster, sites, largest_cluster_list, p_0, p_inf, p_inf_2, susept, s
end

function p_images(p_list)
    p_dic = Dict{Float64,Any}()
    for i in p_list
        p = i
        largest_cluster, sites, largest_cluster_list, p_0 = simulate_bonds(p)
        p_dic[p] = Dict("largest_cluster" => largest_cluster, "sites" => sites, "largest_cluster_list" => largest_cluster_list, "p_0" => p_0[end])
    end
    return p_dic
end