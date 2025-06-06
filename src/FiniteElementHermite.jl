module FiniteElementHermite

import QuadGK: quadgk
import FastGaussQuadrature: gausslegendre
using BandedMatrices
using StaticArrays

include("hermite.jl")
export FE_rep, FE, D, I
export compute_bases, compute_D_bases, compute_both_bases, evaluate, evaluate_inbounds
export νe, D_νe, I_νe, νo, D_νo, I_νo

include("integration.jl")
export inner_product, dual_inner_product, get_quadrature

include("matrices.jl")
export mass_matrix

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
