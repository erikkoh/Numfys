import JSON

include("JSON_functions.jl")


function Grid_Bonds(L::Int)
    bonds = []
    for i in 1:L^2
        row_pos = (i-1)%L + 1
        collum_pos = div(i-1,L) + 1
        row_neighbour = ((row_pos) % L + 1) + (collum_pos-1)*L
        collum_neighbour = (row_pos) + (collum_pos % L  )*L
        neigbours = [row_neighbour, collum_neighbour ]
        push!(bonds,neigbours)
    end
    return bonds
end

square_grid = Grid_Bonds(10)


write_to_JSON(square_grid, "Square_grid")


function Triangular_bonds(L::Int)
    bonds = []
    for i in 1:L*(L+1)/2
        row_num = floor(Int, (sqrt(8 * (i - 1) + 1) - 1) / 2) + 1
        row_start = (row_num * (row_num - 1)) ÷ 2 + 1
        row_length = row_num
        collum_pos = i - row_start + 1

        right_neighbour = if collum_pos < row_length i + 1 else 0 end
        down_left_neighbour = if row_num < L i + row_length else 0 end
        down_right_neighbour = if row_num < L i + row_length + 1 else 0 end

        neighbor = [right_neighbour, down_left_neighbour, down_right_neighbour]
        
        push!(bonds, neighbor)
    end
    return bonds
end

triangular_grid = Triangular_bonds(10)

write_to_JSON(triangular_grid, "Triangular_grid")