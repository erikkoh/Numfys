

function write_to_JSON(bond_list, name)
    json_data = JSON.json(bond_list)
    
    file_path = "Første_øvning/JSON_files/" * name * ".json"

    open(file_path, "w") do file
        write(file, json_data)
    end
    print("Written to file path $file_path")
end