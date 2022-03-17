__precompile__()

module FiniteElementHermite

using QuadGK

include("hermite.jl")
export FE_rep, FE, D, I
export νe, D_νe, I_νe, νo, D_νo, I_νo
export inner_product

end
