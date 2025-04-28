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

# grid = Grid_Bonds(3)["Bonds"]
# println(grid)
# seen = Set()
# number = Dict()
# for i in 1:length(grid)
#     for j in 1:2
#         if grid[i][j] in seen
#             number[grid[i][j]] += 1
#         else
#             push!(seen, grid[i][j])
#             number[grid[i][j]] = 1
#         end
#     end
# end
# println(number)
# println(seen)

#Grid:
# x---x---x 
# | \ | \ |
# x---x---x
# | \ | \ |
# x---x---x

function Triangular_bonds(L::Int)
    bonds = []
    for i in 1:L^2
        row_pos = (i-1)%L + 1
        collum_pos = div(i-1,L) + 1
        row_neighbour = [i,( (row_pos) % L + 1) + (collum_pos-1)*L]
        down_neighbour = [i, (row_pos) + (collum_pos % L  )*L]
        if collum_pos != L
            down_right_neighbour = [i, if row_pos != L down_neighbour[2] + 1 else L-collum_pos + 1 end]
        else
            down_right_neighbour = [i, 1 + L*(L-row_pos) ]
        end
        push!(bonds,row_neighbour)
        push!(bonds,down_neighbour)
        push!(bonds, down_right_neighbour)
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
