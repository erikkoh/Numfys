using Plots

include("Generationg_bonds.jl")
include("Simulate_bonds.jl")
include("Convolution.jl")

@time begin
root_num_nodes_list =[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
root_num_nodes_list = [100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 900, 1000]

root_num_nodes_list =[650, 750, 850, 950]



for i in root_num_nodes_list
    println("Starting with $i nodes")
    JSON_info = Grid_Bonds(i)
    write_to_JSON(JSON_info, "Square_grid_bonds_$i")
    simulation = avarage_values(10,i)
    write_to_JSON(simulation, "Square_Grid_simulation_$i")
    convolution_dic = convelute_data(i)
    write_to_JSON(convolution_dic, "convolution_square_$i")
    println("Finished with $i nodes")
end

end


#7762.171387 seconds (68.98 G allocations: 3.346 TiB, 8.03% gc time, 0.03% compilation time)
#12779.706205 seconds (21.06 G allocations: 967.585 GiB, 3.52% gc time, 0.02% compilation time)