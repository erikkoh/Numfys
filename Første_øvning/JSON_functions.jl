

function write_to_JSON(bond_list,name, number_of_bonds)
    json_data = JSON.json(bond_list)
    
    file_path = "Første_øvning/JSON_files/" * name * ".json"

    open(file_path, "w") do file
        JSON.print(file, number_of_bonds)
        write(file, "\n")
        write(file, json_data)
    end
    println("Written to file path $file_path")
end