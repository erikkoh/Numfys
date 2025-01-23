import JSON
import Random
include("JSON_functions.jl")

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
    end
end

function simulate_bonds(p)
    sites = [-1 for i in 1:(bonds[end][1])]
    clusters_comparison = []
    largest_cluster = [1,-1]
    number_of_activated_bonds = 0
    p_0 = 0
    steps = 0
    while p_0 <= p && steps < 100_000
        largest_cluster = Dict("index" => findfirst(x ->x == minimum(sites),sites),"size" => minimum(sites))
        selected_bond = rand(1:num_bonds)
        node1, node2 = bonds[selected_bond]
        root1 = find_root_node(node1,sites)
        root2 = find_root_node(node2,sites)
        if root1 != root2
            number_of_activated_bonds += 1
            if sites[root1] < sites[root2]
                sites[root1] += sites[root2]
                sites[root2] = root1
                new_root = find_root_node(root1,sites)
                push!(clusters_comparison, sites[new_root]/largest_cluster["size"])
            else
                sites[root2] += sites[root1]
                sites[root1] = root2
                new_root = find_root_node(root2,sites)
                push!(clusters_comparison, sites[new_root]/largest_cluster["size"])
            end
        end
        p_0 = number_of_activated_bonds/num_bonds
        steps += 1
        # Apperently this is a known issue and is not ment to be solved(yet)
        # for j in 1:length(sites)
        #     find_root_node(j)
        # end
    end
    println("P_0: ", p_0)
    println("Number of bonds: ", num_bonds)
    println("Number of activated bonds: ", number_of_activated_bonds)
    largest_cluster_list = [i for i in eachindex(sites) if find_root_node(i,sites) == largest_cluster["index"]]
    return clusters_comparison, largest_cluster, sites, largest_cluster_list, p_0
end

function p_images(p_list)
    p_dic = Dict{Float64,Any}()
    for i in p_list
        p = i
        cluster_comp, largest_cluster, sites, largest_cluster_list, p_0 = simulate_bonds(p)
        p_dic[p] = Dict("cluster_comp" => cluster_comp, "largest_cluster" => largest_cluster, "sites" => sites, "largest_cluster_list" => largest_cluster_list)
    end
    return p_dic
end
