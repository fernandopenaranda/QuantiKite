module QuantiKite

module ChebyshevExpansions
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) QuantiKite
    using LinearAlgebra, SparseArrays, Quantica, CairoMakie
    export 
    include("structures.jl")   
end

end
