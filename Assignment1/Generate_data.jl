using Plots

include("Generationg_bonds.jl")
include("Simulate_bonds.jl")
include("Convolution.jl")



function Generate_data(grid_type::String)
    @time begin
        Ns =[900, 1000]
        function_dic =Dict("Square" =>(i::Int64)-> Grid_Bonds(i), "Triangular" => (i::Int64)-> Triangular_bonds(i), "Honey_comb" => (i::Int64)-> Triangular_bonds(i) )
        for i in Ns
            println("Starting with $i nodes")
            JSON_info = function_dic[grid_type](i)
            write_to_JSON(JSON_info, "$(grid_type)_grid_bonds_$i")
            simulation = avarage_values(1000,i, grid_type)
            write_to_JSON(simulation, "$(grid_type)_Grid_simulation_$i")
            convolution_dic = convelute_data(i, grid_type)
            write_to_JSON(convolution_dic, "convolution_$(grid_type)_$(i)")
            println("Finished with $i nodes")
        end
    end
end

Generate_data("Square")

#7762.171387 seconds (68.98 G allocations: 3.346 TiB, 8.03% gc time, 0.03% compilation time)
#12779.706205 seconds (21.06 G allocations: 967.585 GiB, 3.52% gc time, 0.02% compilation time)
# 770.194204 seconds (8.05 G allocations: 191.440 GiB, 25.72% gc time, 0.27% compilation time) Mayor improvment -> push! is horrible