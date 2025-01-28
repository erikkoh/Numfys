include("JSON_functions.jl")
using .JSONFunctions: write_to_JSON 
#Grid:
# x---x---x
# |   |   |
# x---x---x
# |   |   |
# x---x---x

function Grid_Bonds(L::Int)
    bonds = []
    for i in 1:L^2
        row_pos = (i-1)%L + 1
        collum_pos = div(i-1,L) + 1
        row_neighbour = [i,( (row_pos) % L + 1) + (collum_pos-1)*L]
        collum_neighbour = [i, (row_pos) + (collum_pos % L  )*L]
        push!(bonds,row_neighbour)
        push!(bonds,collum_neighbour)
    end
    return Dict("Number of bonds" => length(bonds), "Bonds" => bonds)
end


#Grid:
#     x
#    / \
#   x---x
#  / \ / \
# x---x---x

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
    return Dict("Number of bonds" => length(bonds), "Bonds" => bonds)
end






#Grid:
#    x---x
#   /     \
#  x       x
#   \     /
#    x---x
function Honey_comb_grid(L::Int)
    bonds = []
    function double_neighbour(i::Int, L::Int)
        row_pos = (i-1)%L + 1
        collum_pos = div(i-1,L) + 1
        right_neighbour = [i,((row_pos) % L + 1) + (collum_pos-1)*L]
        down_neighbour = [i, if collum_pos < L i+L else row_pos end]
        neigbours = [right_neighbour,down_neighbour]
        return neigbours
    end
    
    function single_neighbour(i::Int, L::Int)
        row_pos = (i-1)%L + 1
        collum_pos = div(i-1,L) + 1
        down_neighbour = [i, if collum_pos < L i+L else row_pos end]
        return down_neighbour
    end
    for i in 1:L^2
        row_pos = (i-1)%L + 1
        collum_pos = div(i-1,L) + 1

        if isodd(collum_pos)
            if isodd(i)
                neigbours = single_neighbour(i,L)
            else
                neigbours = double_neighbour(i,L)
            end
        else
            if isodd(i)
                neigbours = double_neighbour(i,L)
            else
                neigbours = single_neighbour(i,L)
            end
        end

        push!(bonds,neigbours)

    end
    return Dict("Number of bonds" => length(bonds), "Bonds" => bonds)
end


