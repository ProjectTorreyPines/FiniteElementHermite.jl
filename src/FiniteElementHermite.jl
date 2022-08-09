__precompile__()

module FiniteElementHermite

import QuadGK: quadgk
import FastGaussQuadrature: gausslegendre
using BandedMatrices
using StaticArrays

include("hermite.jl")
export FE_rep, FE, D, I
export νe, D_νe, I_νe, νo, D_νo, I_νo
export inner_product

include("matrices.jl")
export mass_matrix

end
