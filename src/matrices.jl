"""
    mass_matrix(N::Integer, ρ::AbstractVector{<:Real})

NxN matrix of inner products of all Hermite cubic finite elements on grid ρ
Returns a BandedMatrix as only neighboring finite elements overlap
N.B. - N should be determined from ρ and will be removed later
"""
function mass_matrix(N::Integer, ρ::AbstractVector{<:Real})
    Ms = zeros(eltype(ρ), 2N, 7)
    for j in 1:N
        # First odd with all nearest neighbors
        Ms[2j-1, 1] = 0.0
        if j > 1
            Ms[2j-1, 2] = inner_product(νo, j, νo, j - 1, ρ)
            Ms[2j-1, 3] = inner_product(νo, j, νe, j - 1, ρ)
        end
        Ms[2j-1, 4] = inner_product(νo, j, νo, j, ρ)
        Ms[2j-1, 5] = inner_product(νo, j, νe, j, ρ)
        if j < N
            Ms[2j-1, 6] = inner_product(νo, j, νo, j + 1, ρ)
            Ms[2j-1, 7] = inner_product(νo, j, νe, j + 1, ρ)
        end

        # Then even with all nearest neighbors
        if j > 1
            Ms[2j, 1] = inner_product(νe, j, νo, j - 1, ρ)
            Ms[2j, 2] = inner_product(νe, j, νe, j - 1, ρ)
        end
        Ms[2j, 3] = inner_product(νe, j, νo, j, ρ)
        Ms[2j, 4] = inner_product(νe, j, νe, j, ρ)
        if j < N
            Ms[2j, 5] = inner_product(νe, j, νo, j + 1, ρ)
            Ms[2j, 6] = inner_product(νe, j, νe, j + 1, ρ)
        end
        Ms[2j, 7] = 0.0
    end
    return BandedMatrix(-3 => Ms[4:end, 1], -2 => Ms[3:end, 2], -1 => Ms[2:end, 3], 0 => Ms[:, 4],
        1 => Ms[1:end-1, 5], 2 => Ms[1:end-2, 6], 3 => Ms[1:end-3, 7])
end