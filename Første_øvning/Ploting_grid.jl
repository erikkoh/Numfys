using Plots
include("Simulate_bonds.jl")
include("JSON_functions.jl")

using .JSONFunctions: find_folder

gr()

p_list = [0.25, 0.4, 0.45, 0.5]
function plot_images_points(p_values)
    p_dic = p_images(p_values)
    L = Int((length(p_dic[p_values[end]]["sites"]))^(1/2))
    x_cord = [mod1(i, L) for i in 1:L^2]
    y_cord = [div(i-1, L) + 1 for i in 1:L^2]

    for i in p_values
        largest_cluster_list = p_dic[i]["largest_cluster_list"]
        largest_cluster_y = [div(i-1, L) + 1 for i in (largest_cluster_list)]
        largest_cluster_x = [mod1(i, L) for i in largest_cluster_list]
        points = scatter(x_cord, y_cord, label = "Grid", xlabel = "X-coordinate", ylabel = "Y-coordinate", title = "p= $(p_dic[i]["p_0"])", legend = :topleft)
        scatter!(points, largest_cluster_x, largest_cluster_y, color = :red, label = "Largest Cluster")
        savefig(find_folder("Plots")*"Largest_cluster_p_$(i).png")
    end
end


function plot_images(p_values)
    p_dic = p_images(p_values)
    L = Int((length(p_dic[p_values[end]]["sites"]))^(1/2))
    
    for p in p_values
        sites = p_dic[p]["sites"]
        largest_cluster_list = p_dic[p]["largest_cluster_list"]
        
        grid = fill(0, L, L)
        
        for i in largest_cluster_list
            x = mod1(i, L)
            y = div(i-1, L) + 1
            grid[y, x] = 1
        end
        
        # Plot the grid using a heatmap
        heatmap(grid, c=:grays, xlabel="X-coordinate", ylabel="Y-coordinate", title="p= $(p_dic[p]["p_0"])", color=:reds, label="Largest Cluster")
        savefig(find_folder("Plots") * "Largest_cluster_p_$(p).png")
    end
end

plot_images(p_list)
