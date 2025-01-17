import JSON


function Grid_Bonds(L)
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


function write_to_JSON(bond_list, name)
    json_data = JSON.json(bond_list)
    
    file_path = "Første_øvning/JSON_files/" * name * ".json"

    open(file_path, "w") do file
        write(file, json_data)
    end
    print("Written to file path $file_path")
end

write_to_JSON(square_grid, "Square_grid")