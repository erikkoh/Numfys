include("generate_koch_kurve.jl")
using Plots
using HDF5
using ProgressBars




function does_intersect(p1::Vector{Float64}, p2::Vector{Float64}, point::Vector{Float64})
    if p1[2] > p2[2]
        p1, p2 = p2, p1
    end
    if p1[2] <= point[2] < p2[2]
        return (point[2] + 0.001 - p1[1]) * (p2[1] - p1[1]) / (p2[2] - p1[2]) + p1[1] > point[1]
    end 
    return false
end

function in_boundary(point::Vector{Float64}, boundry::Vector{Vector{Float64}})
    count = 0
    list_length = length(boundry)
    if point in boundry
        return false
    else
        for i in eachindex(boundry)
            if does_intersect(boundry[i], boundry[((i)%(list_length-1))+1], point)
                count += 1
            end
        end
    end
    return count % 2 == 1
end



function in_boundary_2(point::Vector{Float64}, boundary::Vector{Vector{Float64}})
    count = 0
    list_length = length(boundary)

    function angle(a, b)
        return acosd(clamp(a ⋅ b / (norm(a) * norm(b)), -1, 1))
    end
    function cross_product_sign(a, b)
        return a[1] * b[2] - a[2] * b[1]
    end
    for i in 1:list_length
        prev_point = boundary[i]
        next_point = boundary[(i % list_length) + 1]
        vec_prev = prev_point - point
        vec_next = next_point - point
        angle_val = angle(vec_next, vec_prev)
        cross_sign = cross_product_sign(vec_prev, vec_next)
        count += cross_sign < 0 ? -angle_val : angle_val
    end
     return abs(count/360 - 1) <  0.1 
end

function generate_grid(l::Int, checking_func = in_boundary)
    grid_size = max(4^l, 4^3)
    boundary = generate_square(l, grid_size)
    x_boundary = [point[1] for point in boundary]
    y_boundary = [point[2] for point in boundary]
    x_min, x_max = Int(floor(minimum(x_boundary))), Int(ceil(maximum(x_boundary)))   
    y_min, y_max = Int(floor(minimum(y_boundary))), Int(ceil(maximum(y_boundary))) 
    grid_x = collect(x_min:x_max)
    grid_y = collect(y_min:y_max)
    grid = [[Float64(x), Float64(y)] for y in grid_y for x in grid_x]
    grid_z = [checking_func(Float64.(point), boundary) for point in grid]
    grid = [[point[1], point[2], in_boundry] for (point, in_boundry) in zip(grid, grid_z)]
    return hcat(grid...), boundary
end



function comparing_speed()
    for l in 1:4
        time1 = @elapsed begin
            generate_grid(l,in_boundary)
        end
        time2 = @elapsed begin
            generate_grid(l,in_boundary_2)
        end
        println("Ray casting took $time1 and Winding number took $time2")
    end
end

function plotting_grid()
    grid, _ = generate_grid(3, in_boundary_2)
    grid = convert(Array{Float64, 2}, grid)
    z = grid[3,:]
    grid_size = Int(sqrt(length(z)))
    z_matrix = reshape(z, grid_size, grid_size)
    p1 = heatmap(z_matrix, xlabel="X-axis", ylabel="Y-axis", title="Grid")
    display(p1)
end

function checking_if_boundary_is_included()
    grid, boundary = generate_grid(1)
    grid = convert(Array{Float64, 2}, grid)
    println(grid[:,9])
    println(boundary[1])
    similar_points = [grid[i, j] for i in 1:size(grid, 1), j in 1:size(grid, 2) if grid[i, j] in boundary]
    return similar_points, length(similar_points) == length(boundary)
end


function write_grid_to_file(l::Int)
    grid, boundary = generate_grid(l)
    
    grid = convert(Array{Float64, 2}, grid)
    num_rows = length(boundary)
    num_cols = length(boundary[1])    
    flat_vector = vcat(boundary...)
    boundary_matrix = reshape(flat_vector, num_cols, num_rows)

    h5open("./Assignment2/HDF5_files/grid_boundary$l.h5", "w") do file
        write(file, "grid", grid)
        write(file, "boundary", boundary_matrix)
    end
end

comparing_speed()

for i in 1:4
    println("started writing for l=$i")    
    write_grid_to_file(i)
end