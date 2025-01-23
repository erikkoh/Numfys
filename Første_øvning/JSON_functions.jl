using JSON

module JSON_functions

export write_to_JSON, find_folder

function write_to_JSON(bond_list,name, number_of_bonds)
    data_dic = Dict("Number of bonds" => number_of_bonds, "Bonds" => bond_list)
    json_data = JSON.json(data_dic)
    
    file_path = find_folder("JSON_files")* name * ".json"
    println(file_path)

    open(file_path, "w") do file
        write(file, json_data)
    end
    println("Written to file path $file_path")
end

function find_folder(folder_name::String)
    folder_path = ""
    for (root, dirs, files) in walkdir(@__DIR__)
        if occursin(folder_name, root)
            return joinpath(root, folder_name)
        end
    end
    error("Folder not found")
end
end