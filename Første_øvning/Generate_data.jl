using Plots

include("Generationg_bonds.jl")
include("Simulate_bonds.jl")


root_num_nodes_list =[100,200,300,400,500,700,900,1000]

for i in root_num_nodes_list
    println("Starting with $i nodes")
    JSON_info = Grid_Bonds(i)
    write_to_JSON(JSON_info, "Square_grid_bonds_$i")
    simulation = simulate_bonds(i^2)
    values_dic = Dict("p_inf" => simulation[5],"p_inf_2" => simulation[6], "susept" => simulation[7], "s" => simulation[8], "p" => simulation[4])
    write_to_JSON(values_dic, "Square_Grid_simulation_$i")
    convolution_dic = convelute_data(i)
    write_to_JSON(convolution_dic, "convolution_$i")
    println("Finished with $i nodes")
end


