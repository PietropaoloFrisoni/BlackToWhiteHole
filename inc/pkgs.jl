using SL2Cfoam
using JLD2
using Distributed
using DataFrames
using Dates
using CSV
using HalfIntegers
using LoopVectorization
using LinearAlgebra

# useful code script to install all required pkgs in a single loop

#=
vec = ["Distributions", "HalfIntegers", "LoopVectorization", "JLD2", "LinearAlgebra", "Distributed", "CSV", "DataFrames", "Dates"]
import Pkg
for p in vec
Pkg.add("$p")
end
=#
