include("checking_points.jl")
using SparseArrays
using Arpack





function sparse_descrete_lapclacian(vector_length::Int, row_length::Int, grid::Matrix{Float64})
laplacian_matrix = spzeros(vector_length, vector_length)
for i in 1:size(grid, 1)
    if grid[i, 3] == 1
        next_x = i + 1
        next_y = i + row_length
        previous_x = i - 1
        previous_y = i - row_length

        if next_x <= vector_length && grid[mod(i + 1, vector_length), 3] == 1
            laplacian_matrix[i, next_x] = -1
        end
        if next_y <= vector_length &&  grid[mod(i + row_length, vector_length), 3] == 1
            laplacian_matrix[i, next_y] = -1
        end
        if previous_x >= 1 && grid[previous_x, 3] == 1
            laplacian_matrix[i, previous_x] = -1
        end
        if previous_y >= 1  && grid[previous_y, 3] == 1
            laplacian_matrix[i, previous_y] = -1
        end

        laplacian_matrix[i, i] = 4
    end
end

return laplacian_matrix
end


function find_eigenvalues_eigenvectors(laplacian_matrix)
    λ, ϕ = eigs(laplacian_matrix, nev = 10, which=:SM);
    return λ, ϕ
end


function plotting_eigenvectors(laplacian_matrix, boundry)
    x_collum = [point[1] for point in boundry]
    y_collum = [point[2] for point in boundry]
    λ, ϕ = find_eigenvalues_eigenvectors(laplacian_matrix)
    println(λ)
    p1 = plot(x_collum, y_collum)
    max_dim = Int(sqrt(length(ϕ[:,1])))
    p1 = plot!(heatmap(reshape(Real.(ϕ[:,1]), (max_dim, max_dim))))
    display(p1)
end
grid, boundry = generate_grid(4)
laplacian_matrix = sparse_descrete_lapclacian(length(grid[:,1]), Int(maximum(grid[:,1])), grid)
plotting_eigenvectors(laplacian_matrix, boundry)


