#==============================================================
Define Hermite cubic finite elements
===============================================================#

"""
    quad_deriv(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Real}

Fit a parabola going through three points (x1, y1), (x2, y2), and (x3, y3)
  and return derivative at x2
"""
function quad_deriv(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Real}
    hl = x2 - x1
    hu = x3 - x2
    return ((y3 - y2) * hl^2 + (y2 - y1) * hu^2) / (hu * hl * (hu + hl))
end

"""
    fit_derivative(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Returns dy/dx at every point in x, based on local quadratic fit
"""
function fit_derivative(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    dy_dx = similar(x)
    dy_dx[1] = quad_deriv(x[3], x[1], x[2], y[3], y[1], y[2])
    dy_dx[end] = quad_deriv(x[end-1], x[end], x[end-2], y[end-1], y[end], y[end-2])
    for k in 2:length(x)-1
        dy_dx[k] = quad_deriv(x[k-1], x[k], x[k+1], y[k-1], y[k], y[k+1])
    end
    return dy_dx
end

function which_region(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    ρk = ρ[k]
    if x <= ρk
        if k > 1
            ρkl = ρ[k-1]
            region = x > ρkl ? :in_low : :out_low
            return region, ρkl, ρk
        elseif x < ρk
            return :out_low, -Inf, ρk
        end
    end
    if x >= ρk
        if k < length(ρ)
            ρku = ρ[k+1]
            region = x < ρku ? :in_up : :out_up
            return region, ρk, ρku
        elseif x != ρk
            return :out_up, ρk, Inf
        end
    end
    return :error, -Inf, Inf
end

# Even element

function νel(x::Real, ρkl::Real, ρk::Real)
    t = (x - ρk) / (ρk - ρkl)
    return t^2 * (-2.0 * t - 3.0) + 1.0
end

function νeu(x::Real, ρk::Real, ρku::Real)
    t = (x - ρk) / (ρku - ρk)
    return t^2 * (2.0 * t - 3.0) + 1.0
end

function νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return νel(x, ρ1, ρ2)
    region === :in_up  && return νeu(x, ρ1, ρ2)
    return 0.0
end

function D_νel(x::Real, ρkl::Real, ρk::Real)
    ihl = 1.0 / (ρk - ρkl)
    t = (x - ρk) * ihl
    return -6.0 * ihl * t * (t + 1.0)
end

function D_νeu(x::Real, ρk::Real, ρku::Real)
    ihu = 1.0 / (ρku - ρk)
    t = (x - ρk) * ihu
    return 6.0 * ihu * t * (t - 1.0)
end

function D_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return D_νel(x, ρ1, ρ2)
    region === :in_up  && return D_νeu(x, ρ1, ρ2)
    return 0.0
end

function I_νel(x::Real, ρkl::Real, ρk::Real)
    hl = ρk - ρkl
    t = (x - ρk) / hl
    return hl * (-0.5 * t^4 - t^3 + t + 0.5)
end

function I_νeu(x::Real, ρk::Real, ρku::Real)
    hu = ρku - ρk
    t = (x - ρk) / hu
    return hu * t * (t^2 * (0.5 * t - 1.0) + 1.0)
end

function I_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    ρk = ρ[k]

    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :out_low && return 0.0
    region === :in_low  && return I_νel(x, ρ1, ρ2)

    Il = k > 1 ? 0.5 * (ρk - ρ[k-1]) : 0.0
    region === :in_up  && return Il + I_νeu(x, ρ1, ρ2)
    Iu = k < length(ρ) ? 0.5 * (ρ[k+1] - ρk) : 0.0
    region === :out_up && return Il + Iu
end

# Odd element

function νol(x::Real, ρkl::Real, ρk::Real)
    hl = ρk - ρkl
    t = (x - ρk) / hl
    return hl * t * (t + 1.0)^2
end

function νou(x::Real, ρk::Real, ρku::Real)
    hu = ρku - ρk
    t = (x - ρk) / hu
    return hu * t * (t - 1.0)^2
end

function νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return νol(x, ρ1, ρ2)
    region === :in_up  && return νou(x, ρ1, ρ2)
    return 0.0
end

function D_νol(x::Real, ρkl::Real, ρk::Real)
    t = (x - ρk) / (ρk - ρkl)
    return t * (3.0 * t + 4.0) + 1.0
end

function D_νou(x::Real, ρk::Real, ρku::Real)
    t = (x - ρk) / (ρku - ρk)
    return t * (3.0 * t - 4.0) + 1.0
end

function D_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return D_νol(x, ρ1, ρ2)
    region === :in_up  && return D_νou(x, ρ1, ρ2)
    return 0.0
end

const one_twelf  = 1.0 / 12.0
const two_thirds = 2.0 / 3.0

function I_νol(x::Real, ρkl::Real, ρk::Real)
    hl = ρk - ρkl
    t = (x - ρk) / hl
    return hl^2 * (0.25 * t^4 + two_thirds * t^3 + 0.5 * t^2 - one_twelf)
end

function I_νou(x::Real, ρk::Real, ρku::Real)
    hu = ρku - ρk
    t = (x - ρk) / hu
    return (hu * t)^2 * (t * (0.25 * t - two_thirds) + 0.5)
end

function I_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    ρk = ρ[k]

    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :out_low && return 0.0
    region === :in_low  && return I_νol(x, ρ1, ρ2)

    hl = k > 1 ? ρk - ρ[k-1] : 0.0
    Il = -one_twelf * hl^2
    region === :in_up  && return Il + I_νou(x, ρ1, ρ2)

    hu = k < length(ρ) ? ρ[k+1] - ρk : 0.0
    Iu = one_twelf * hu^2
    return Il + Iu
end


# FE Structure

function hermite_coeffs(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    dy_dx = fit_derivative(x, y)
    C = Vector{eltype(y)}(undef, 2 * length(x))
    @inbounds for i in eachindex(x)
        ti = 2i
        C[ti-1] = dy_dx[i]
        C[ti] = y[i]
    end
    return C
end

struct FE_rep{S <: AbstractVector{<:Real}, T <: AbstractVector{<:Real}}
    x::S
    coeffs::T
    function FE_rep{S, T}(x::S, coeffs::T) where {S <: AbstractVector{<:Real}, T<:AbstractVector{<:Real}}
        return length(coeffs) == 2length(x) ? new{S, T}(x, coeffs) : throw(DimensionMismatch)
    end
end
function FE_rep(x::S, coeffs::T) where {S <: AbstractVector{<:Real}, T<:AbstractVector{<:Real}}
    return FE_rep{S, T}(x, coeffs)
end

FE(x, y) = FE_rep(x, hermite_coeffs(x, y))

# create FE_rep for data (x, y), but mapped to grid
function FE(grid::AbstractVector{<:Real}, xy::Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}})
    f = FE(xy...)
    C = Vector{eltype(grid)}(undef, 2 * length(grid))
    @inbounds for (i, g) in enumerate(grid)
        ti = 2i
        C[ti-1] = D(f, g)
        C[ti] = f(g)
    end
    return FE_rep(grid, C)
end


@inline function compute_bases(X::AbstractVector{<:Real}, x::Real)
    k = searchsortedlast(X, x)
    k == length(X) && (k -= 1)

    @inbounds ρk = X[k]
    @inbounds ρku = X[k+1]
    nu_ou = νou(x, ρk, ρku)
    nu_eu = νeu(x, ρk, ρku)
    nu_ol = νol(x, ρk, ρku)
    nu_el = νel(x, ρk, ρku)
    return k, nu_ou, nu_eu, nu_ol, nu_el
end

@inline function evaluate(Y::FE_rep, k::Integer, nu_ou::Real, nu_eu::Real, nu_ol::Real, nu_el::Real)
    tk = 2k
    y  = Y.coeffs[tk-1] * nu_ou
    y += Y.coeffs[tk  ] * nu_eu
    y += Y.coeffs[tk+1] * nu_ol
    y += Y.coeffs[tk+2] * nu_el
    return y
end

@inline function evaluate_inbounds(Y::FE_rep, k::Integer, nu_ou::T, nu_eu::T, nu_ol::T, nu_el::T) where {T<:Real}
    tk = 2k
    @inbounds y  = Y.coeffs[tk-1] * nu_ou + Y.coeffs[tk  ] * nu_eu
    @inbounds y += Y.coeffs[tk+1] * nu_ol + Y.coeffs[tk+2] * nu_el
    return y
end

function (Y::FE_rep)(x::Real)
    k, nu_ou, nu_eu, nu_ol, nu_el = compute_bases(Y.x, x)
    y = evaluate_inbounds(Y, k, nu_ou, nu_eu, nu_ol, nu_el)
    return y
end

@inline function compute_D_bases(X::AbstractVector{<:Real}, x::Real)
    k = searchsortedlast(X, x)
    k == length(X) && (k -= 1)

    @inbounds ρk = X[k]
    @inbounds ρku = X[k+1]
    D_nu_ou = D_νou(x, ρk, ρku)
    D_nu_eu = D_νeu(x, ρk, ρku)
    D_nu_ol = D_νol(x, ρk, ρku)
    D_nu_el = D_νel(x, ρk, ρku)
    return k, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el
end

@inline function compute_both_bases(X::AbstractVector{<:Real}, x::Real)
    k = searchsortedlast(X, x)
    k == length(X) && (k -= 1)

    @inbounds ρk = X[k]
    @inbounds ρku = X[k+1]
    nu_ou = νou(x, ρk, ρku)
    nu_eu = νeu(x, ρk, ρku)
    nu_ol = νol(x, ρk, ρku)
    nu_el = νel(x, ρk, ρku)
    D_nu_ou = D_νou(x, ρk, ρku)
    D_nu_eu = D_νeu(x, ρk, ρku)
    D_nu_ol = D_νol(x, ρk, ρku)
    D_nu_el = D_νel(x, ρk, ρku)
    return k, nu_ou, nu_eu, nu_ol, nu_el, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el
end

function D(Y::FE_rep, x::Real)
    k, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el = compute_D_bases(Y.x, x)
    dy_dx = evaluate_inbounds(Y, k, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el)
    return dy_dx
end

function I(Y::FE_rep, x::Real)
    K = min(searchsortedfirst(Y.x, x), length(Y.x))
    yint = 0.0
    @inbounds for k in 1:K
        tk = 2k
        yint += Y.coeffs[tk-1] * I_νo(x, k, Y.x)
        yint += Y.coeffs[tk] * I_νe(x, k, Y.x)
    end
    return yint
end