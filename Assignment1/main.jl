include("Generate_data.jl")
include("Measurments.jl")

function main()
    grids = ["Square", "Triangular", "Honey_comb"]
    for grid in girds
        Generate_data(grid)
        find_critical_exponents(grid)    
    end
end

if __name__ == "main.jl"
    main()
end