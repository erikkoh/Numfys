include("functions.jl")


function main()
    seed = 1234242
    "Problem 3.1"
    "3.1.1"
    Random.seed!(seed)
    ploting_magnetization()

    "3.1.2 and 3.1.3"
    Random.seed!(seed)
    ploting_heat_mag_and_energy()
    Random.seed!(seed)
    plot_L_dependence_of_alpha()
    
    "3.1.4"
    Random.seed!(seed)
    plot_different_constant_points()
    "Used for the appendix A"
    Random.seed!(seed)
    plot_different_constant_points(range(1.5,3.0, length=50))

    "3.1.5"
    Random.seed!(seed)
    simulated_anehiling_with_constant_points_all()

    "Problem 3.2"
    "3.2.1"
    Random.seed!(seed)
    plot_with_h_field(2.3, 0.01)
    
    "3.2.2"
    Random.seed!(seed*13) #The different seed was mainly due to a very extreme artifact in the data which made the rest of the plot unusable
    plot_susceptibility_for_different_H()
    
    "3.2.3"
    Random.seed!(seed)
    scaling_of_susceptibility()
    Random.seed!(seed)
    data_collaps_of_susptibility()
end
main()