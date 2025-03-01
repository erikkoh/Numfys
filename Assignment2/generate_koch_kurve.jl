using LinearAlgebra
using Plots



function rotate_vec(vector, n)
    rotation_matrices = [
        [1 0; 0 1],    # 0 degrees
        [0 -1; 1 0],   # 90 degrees
        [-1 0; 0 -1],  # 180 degrees
        [0 1; -1 0]    # 270 degrees
    ]
    rotated_vector = rotation_matrices[n+1] * vector
    return rotated_vector
end

function segment_point(p1, p2)
    v1 = p2 .- p1
    relative_length = norm(v1) / 4
    v1 = v1 / norm(v1)
    angle = atan(v1[2], v1[1])
    angle_degrees = rad2deg(angle) % 360
    n = mod(round(Int, angle_degrees / 90), 4)
    vectors_list = [
        [relative_length, 0.0], 
        [0.0, relative_length], 
        [relative_length, 0.0], 
        [0.0, -relative_length], 
        [0.0, -relative_length], 
        [relative_length, 0.0], 
        [0.0, relative_length]
    ]
    point = p1
    points_list = [p1]
    for vectors in vectors_list
        rotated_vec = rotate_vec(vectors, n)
        point = point .+ rotated_vec
        push!(points_list, point)
    end
    return points_list
end

function new_line(line)
    new_line = []
    for i in 1:(length(line)-1)
        segment_points = segment_point(line[i], line[i+1])
        append!(new_line, segment_points)
    end
    push!(new_line, line[end])
    return new_line
end

function fract(l, line)
    for _ in 1:l
        line = new_line(line)
    end
    return line
end


function plotting_line(l, size)
    line = [[0.0, 0.0], [size, 0.0]]
    p1 = plot()
    for i in 1:l
        line = fract(i, line)
        x_cords = [p[1] for p in line]
        y_cords = [p[2] for p in line]
        plot!(x_cords, y_cords, label="")
    end
    display(p1)
end


function plotting_square(l, size)
    p1 = plot()
    for i in 0:l
        total_points = generate_square_without_shift(i, size)
        x_cords = [p[1] for p in total_points]
        y_cords = [p[2] for p in total_points]
        plot!(x_cords, y_cords, label="")
    end
    display(p1)
end

function generate_square_without_shift(l, size)
    lines = [
        [[0.0, 0.0], [size, 0.0]],  
        [[size, 0.0], [size, size]],
        [[size, size], [0.0, size]],
        [[0.0, size], [0.0, 0.0]]
    ]
    total_points = []
    for line in lines
        line = fract(l, line)
        append!(total_points, line)
    end
    return total_points
end

function generate_square(l, size)
    lines = [
        [[0.0, 0.0], [size, 0.0]],  
        [[size, 0.0], [size, size]],
        [[size, size], [0.0, size]],
        [[0.0, size], [0.0, 0.0]]
    ]
    total_points = []
    for line in lines
        line_points = fract(l, line)
        append!(total_points, line_points)
    end
    min_x = minimum([p[1] for p in total_points])
    min_y = minimum([p[2] for p in total_points])
    adjusted_points = [[p[1] - min_x, p[2] - min_y] for p in total_points]
    return adjusted_points
end

plotting_square(6, 4^3)
