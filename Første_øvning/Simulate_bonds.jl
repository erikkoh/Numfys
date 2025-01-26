import JSON
import Random
using ProgressBars
include("JSON_functions.jl")

using  .JSONFunctions: find_folder, write_to_JSON

Random.seed!(144)

JSON_info = JSON.parsefile(find_folder("JSON_files") * "Square_grid.json")

bonds = JSON_info["Bonds"]

num_bonds =JSON_info["Number of bonds"]
num_nodes = bonds[end][1]

function swaping_bonds(bonds_list)
    println("Swaping bonds")
    for i in ProgressBar(1:(length(bonds_list)-1))
        candidate_number = rand(i+1:num_bonds)
        temp = bonds_list[i]
        bonds_list[i] = bonds_list[candidate_number]
        bonds_list[candidate_number] = temp
    end
    println("Done!")
    return bonds_list
end

bonds = swaping_bonds(bonds)

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

function simulate_bonds(n=num_bonds, bonds_list=bonds)
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
    for i in ProgressBar(1:n)
        node1, node2 = bonds_list[i]
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
            end
            avarage_s += sites[new_root]^2
            s_step = (avarage_s-(N*p_inf[end])^2)/(N*(1-p_inf[end])) 
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

function p_images(p_list)
    p_dic = Dict{Float64,Any}()
    for p in p_list
        i = Int(floor(p*num_bonds))
        largest_cluster, sites, largest_cluster_list, p_0 = simulate_bonds(i, bonds )
        p_dic[p] = Dict("largest_cluster" => largest_cluster, "sites" => sites, "largest_cluster_list" => largest_cluster_list, "p_0" => p_0[end])
    end
    return p_dic
end


function avarage_values(iterations::Int64, bond_list)
    p_inf_list = []
    susept_list = []
    s_list = []
    p_list = []
    for i in ProgressBar(1:iterations)
        bond_list = swaping_bonds(bond_list)
        result = simulate_bonds(num_bonds, bond_list)
        p = result[4]
        p_inf = result[5]
        susept = result[7]
        s = result[8]
        println("length of p_inf: ", length(p_inf))
        push!(p_list,p)
        push!(p_inf_list,p_inf)
        push!(susept_list,susept)
        push!(s_list,s)
    end

    p_inf_avarage = [sum([p[i] for p in p_inf_list ])/length(p_inf_list) for i in eachindex(p_inf_list[1])]
    susept_avarage = [sum([p[i] for p in susept_list ])/length(susept_list) for i in eachindex(susept_list[1])]
    s_avarage = [sum([p[i] for p in s_list ])/length(s_list) for i in eachindex(s_list[1])]
    p_avarage = [sum([p[i] for p in p_list])/length(p_list) for i in eachindex(p_list[1])]
    values_dic = Dict("p_inf" => p_inf_avarage, "susept" => susept_avarage, "s" => s_avarage, "p" => p_avarage)
    write_to_JSON(values_dic,"avarage_values_$num_bonds")
end
