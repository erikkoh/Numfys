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


function Triangular_bonds(L::Int)
    bonds = []
    row_start_func = x -> (x*(x-1)) ÷ 2 + 1
    artimic_sum = x -> (x*(x+1)) / 2
    for i in 1:Int(artimic_sum(L))
        row_num = floor(Int, (sqrt(8 * (i - 1) + 1) - 1) / 2) + 1
        row_start = Int(row_start_func(row_num))
        row_length = row_num
        row_pos = i - row_start + 1
        
        
        right_neighbour = if row_pos < row_length i + 1 else row_start end
        down_left_neighbour = if row_num < L i + row_length else Int(artimic_sum(row_pos)) end
        down_right_neighbour = if row_num < L i + row_length + 1 else Int(row_start_func(L-row_pos + 1)) end
        
        neighbor = [right_neighbour, down_left_neighbour, down_right_neighbour]
        
        push!(bonds, neighbor)
    end
    return bonds
end


function Honey_comb_grid(L::Int)

end

square_grid = Grid_Bonds(10)
triangular_grid = Triangular_bonds(3)

write_to_JSON(square_grid, "Square_grid", length(square_grid)*2)
write_to_JSON(triangular_grid, "Triangular_grid", length(triangular_grid)*3)
