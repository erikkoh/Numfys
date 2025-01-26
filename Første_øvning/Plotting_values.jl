using Plots

include("Simulate_bonds.jl")


simulation = simulate_bonds()

p_0 = simulation[4]
p_inf = simulation[5]
s = simulation[8]
println("Started plotting")
p1 = plot(p_0, s, label = "<s>", xlabel = "p_0", ylabel = "s", title = "<s> as a function of p_0")
p2 = plot(p_0, p_inf, label = "p_inf", xlabel = "p_0", ylabel = "p_inf", title = "p_inf as a function of p_0")
plot(p1, p2, layout = (2,1))
savefig(find_folder("Plots") * "values_plot.png")
println("Finished plotting")