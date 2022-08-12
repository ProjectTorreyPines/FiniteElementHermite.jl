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
    return -2.0 * t^3 - 3.0 * t^2 + 1.0
end

function νeu(x::Real, ρk::Real, ρku::Real)
    t = (x - ρk) / (ρku - ρk)
    return 2.0 * t^3 - 3.0 * t^2 + 1.0
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
    return hu * t * (0.5 * t^3 - t^2 + 1.0)
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
    return 3.0 * t^2 + 4.0 * t + 1.0
end

function D_νou(x::Real, ρk::Real, ρku::Real)
    t = (x - ρk) / (ρku - ρk)
    return 3.0 * t^2 - 4.0 * t + 1.0
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
    return (hu * t)^2 * (0.25 * t^2 - two_thirds * t + 0.5)
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
    N = length(x)
    C = Vector{typeof(y[1])}(undef, 2N)
    @inbounds for i in 1:N
        C[2i-1] = dy_dx[i]
        C[2i] = y[i]
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

function (Y::FE_rep)(x::Real)
    k = searchsortedlast(Y.x, x)
    @inbounds y  = Y.coeffs[2k-1] * νo(x, k, Y.x)
    @inbounds y += Y.coeffs[2k]   * νe(x, k, Y.x)
    @inbounds if x != Y.x[k]
        @inbounds y += Y.coeffs[2k+1] * νo(x, k+1, Y.x)
        @inbounds y += Y.coeffs[2k+2] * νe(x, k+1, Y.x)
    end
    return y
end
function D(Y::FE_rep, x::Real)
    k = searchsortedlast(Y.x, x)
    @inbounds dy_dx  =  Y.coeffs[2k-1] * D_νo(x, k, Y.x)
    @inbounds dy_dx += Y.coeffs[2k]    * D_νe(x, k, Y.x)
    @inbounds if x != Y.x[k]
        @inbounds dy_dx += Y.coeffs[2k+1] * D_νo(x, k+1, Y.x)
        @inbounds dy_dx += Y.coeffs[2k+2] * D_νe(x, k+1, Y.x)
    end
    return dy_dx
end
function I(Y::FE_rep, x::Real)
    K = min(searchsortedfirst(Y.x, x), length(Y.x))
    yint = 0.0
    @inbounds for k in 1:K
        yint += Y.coeffs[2k-1] * I_νo(x, k, Y.x)
        yint += Y.coeffs[2k] * I_νe(x, k, Y.x)
    end
    return yint
end

#==============================================================
Define inner products
===============================================================#

function gl_preallocate(N_gl::Integer)
    gξ = zeros(N_gl, N_gl)
    gw = zeros(N_gl, N_gl)
    for i in 1:N_gl
        gξ[1:i, i], gw[1:i, i] = gausslegendre(i)
    end
    return SMatrix{N_gl, N_gl}(gξ),  SMatrix{N_gl, N_gl}(gw)
end

const N_gl = 50
const gξ_pa, gw_pa = gl_preallocate(N_gl)

function integrate(f, lims::SVector; tol::Real=eps(typeof(1.0)))
    return quadgk(f, lims..., atol=tol, rtol=sqrt(tol), maxevals=100)[1]
end

integrate(f, lims::SVector, order::Nothing; tol::Real=eps(typeof(1.0))) = integrate(f, lims; tol)

function integrate(f, lims::SVector{2,<:Real}, order::Integer)
    @assert order <= 50
    I = 0.0
    dxdξ = 0.5*(lims[2] - lims[1])
    xavg = 0.5*(lims[2] + lims[1])
    for k in 1:order
        I += f(dxdξ * gξ_pa[k, order] + xavg) * gw_pa[k, order] * dxdξ
    end
    return I
end

function integrate(f, lims::SVector{3,<:Real}, order::Integer)
    return integrate(f, SVector(lims[1], lims[2]), order) + integrate(f, SVector(lims[2], lims[3]), order)
end

function limits(k1::Integer, k2::Integer, ρ::AbstractVector{<:Real})
    k1 != k2        && return SVector(min(ρ[k1], ρ[k2]), max(ρ[k1], ρ[k2]))
    k1 == 1         && return SVector(ρ[1], ρ[2])
    k1 == length(ρ) && return SVector(ρ[end-1], ρ[end])
    return SVector(ρ[k1-1], ρ[k1], ρ[k1+1])
end

function limits(k::Integer, ρ::AbstractVector{<:Real})
    k == 1         && return SVector(ρ[1], ρ[2])
    k == length(ρ) && return SVector(ρ[end-1], ρ[end])
    return SVector(ρ[k-1], ρ[k], ρ[k+1])
end

function inner_product(nu1, k1::Integer, nu2, k2::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing, Integer}=nothing)
    abs(k1 - k2) > 1 && return 0.0
    integrand(x) = nu1(x, k1, ρ) * nu2(x, k2, ρ)
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function inner_product(f, nu1, k1::Integer, nu2, k2::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing, Integer}=nothing)
    abs(k1 - k2) > 1 && return 0.0
    integrand(x) = f(x) * nu1(x, k1, ρ) * nu2(x, k2, ρ)
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function inner_product(nu1, k1::Integer, f, fnu2, g, gnu2, k2::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing, Integer}=nothing)
    abs(k1 - k2) > 1 && return 0.0
    integrand(x) = nu1(x, k1, ρ) * (f(x) * fnu2(x, k2, ρ) + g(x) * gnu2(x, k2, ρ))
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function inner_product(f::Function, nu::Function, k::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing, Integer}=nothing)
    integrand(x) = f(x) * nu(x, k, ρ)
    return integrate(integrand, limits(k, ρ), order)
end