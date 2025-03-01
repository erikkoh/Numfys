include("generate_koch_kurve.jl")
import .KochCurve: generate_square



function does_intersect(p1::Vector{Float64}, p2::Vector{Float64}, point::Vector{Float64})
    if p1[2] > p2[2]
        p1, p2 = p2, p1
    end
    if p1[2] <= point[2] < p2[2]
        x_intersect = (point[2] + 0.001 - p1[1]) * (p2[1] - p1[1]) / (p2[2] - p1[2]) + p1[1]
        return x_intersect > point[1]
    end 
    return false
end

function in_boundry(point::Vector{Float64}, boundry::Vector{Vector{Float64}})
    count = 0
    for i in eachindex(boundry)
        if does_intersect(boundry[i], boundry[mod((i+1), length(boundry))+1], point)
            count += 1
        end
    end
    return count % 2 == 1
end


function generate_grid(l)
    grid_size = 4^l
    if l< 3
        grid_size = 4^3
    end
    boundary = generate_square(l, grid_size)
    x_boundry = [point[1] for point in boundary]
    y_boundry = [point[2] for point in boundary]
    x_min, x_max = minimum(x_boundry), maximum(x_boundry)
    y_min, y_max = minimum(y_boundry), maximum(y_boundry)
    println((x_min, x_max, y_min, y_max))
    grid_x = range(x_min, x_max, length=grid_size)
    grid_y = range(y_min, y_max, length=grid_size)
    grid = [[x, y] for x in grid_x for y in grid_y]
    grid_z = [in_boundry(point, boundary) for point in grid]
    grid_x_coords = [point[1] for point in grid]
    grid_y_coords = [point[2] for point in grid]
    return hcat(grid_x_coords, grid_y_coords, grid_z)
end

function checking_speed()
    @time begin
        grid = generate_grid(4)
        println(size(grid))
        println(grid[end, :])
    end
end