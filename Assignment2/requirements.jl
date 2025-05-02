# requirements.jl
# This script installs the required packages for the project.

using Pkg

# List of required packages
required_packages = [
    "Plots",
    "HDF5",
    "ProgressBars",
    "SparseArrays",
    "Arpack"
    ]

# Install missing packages
for pkg in required_packages
    if !haskey(Pkg.dependencies(), pkg)!
        println("Installing $pkg...")
        Pkg.add(pkg)
    else
        println("$pkg is already installed.")
    end
end

println("All required packages are installed.")