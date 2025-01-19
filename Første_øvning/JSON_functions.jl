using JSON

function write_to_JSON(bond_list,name, number_of_bonds)
    data_dic = Dict("Number of bonds" => number_of_bonds, "Bonds" => bond_list)
    json_data = JSON.json(data_dic)
    
    file_path = "Første_øvning/JSON_files/" * name * ".json"

    open(file_path, "w") do file
        write(file, json_data)
    end
    println("Written to file path $file_path")
end