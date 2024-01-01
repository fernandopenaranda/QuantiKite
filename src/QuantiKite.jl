module QuantiKite

module ChebyshevExpansions
    # Use README as the docstring of the module:
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) QuantiKite
    using LinearAlgebra, SparseArrays, Quantica, CairoMakie
    export h5gen
    export Configuration, Dos, Ldos, Conductivity_optical, Conductivity_dc, 
        Conductivity_optical_non_linear, Singleshot_conductivity_dc, Arpes
    include("structures.jl")
    INCLUDE("")
end

end
