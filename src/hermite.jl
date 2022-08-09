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

function νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    if k > 1 && x >= ρ[k-1] && x <= ρ[k]
        # we're in the lower half
        t = (x - ρ[k]) / (ρ[k] - ρ[k-1])
        return -2.0 * t^3 - 3.0 * t^2 + 1.0
    elseif k < length(ρ) && x >= ρ[k] && x <= ρ[k+1]
        # we're in the upper half
        t = (x - ρ[k]) / (ρ[k+1] - ρ[k])
        return 2.0 * t^3 - 3.0 * t^2 + 1.0
    end
    return 0.0
end

function D_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    if k > 1 && x >= ρ[k-1] && x <= ρ[k]
        # we're in the lower half
        ihl = 1.0 / (ρ[k] - ρ[k-1])
        t = (x - ρ[k]) * ihl
        return -6.0 * ihl * t * (t + 1.0)
    elseif k < length(ρ) && x >= ρ[k] && x <= ρ[k+1]
        # we're in the upper half
        ihu = 1.0 / (ρ[k+1] - ρ[k])
        t = (x - ρ[k]) * ihu
        return 6.0 * ihu * t * (t - 1.0)
    end
    return 0.0
end

function I_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    k == 1 ? hl = 0.0 : hl = ρ[k] - ρ[k-1]
    k == length(ρ) ? hu = 0.0 : hu = ρ[k+1] - ρ[k]

    if x <= ρ[k]
        Iu = 0.0
        if k > 1 && x >= ρ[k-1]
            # we're in the lower half
            t = (x - ρ[k]) / hl
            Il = hl * (-0.5 * t^4 - t^3 + t + 0.5)
        else
            Il = 0.0
        end
    elseif x >= ρ[k]
        Il = 0.5 * hl
        if k < length(ρ) && x <= ρ[k+1]
            # we're in the upper half
            t = (x - ρ[k]) / hu
            Iu = hu * t * (0.5 * t^3 - t^2 + 1.0)
        else
            Iu = 0.5 * hu
        end
    end
    return Il + Iu
end

function νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    if k > 1 && x >= ρ[k-1] && x <= ρ[k]
        # we're in the lower half
        hl = ρ[k] - ρ[k-1]
        t = (x - ρ[k]) / hl
        return hl * t * (t + 1.0)^2
    elseif k < length(ρ) && x >= ρ[k] && x <= ρ[k+1]
        # we're in the upper half
        hu = ρ[k+1] - ρ[k]
        t = (x - ρ[k]) / hu
        return hu * t * (t - 1.0)^2
    end
    return 0.0
end

function D_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    if k > 1 && x >= ρ[k-1] && x <= ρ[k]
        # we're in the lower half
        t = (x - ρ[k]) / (ρ[k] - ρ[k-1])
        return 3.0 * t^2 + 4.0 * t + 1.0
    elseif k < length(ρ) && x >= ρ[k] && x <= ρ[k+1]
        # we're in the upper half
        t = (x - ρ[k]) / (ρ[k+1] - ρ[k])
        return 3.0 * t^2 - 4.0 * t + 1.0
    end
    return 0.0
end

function I_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    k == 1 ? hl = 0.0 : hl = ρ[k] - ρ[k-1]
    if x <= ρ[k]
        Iu = 0.0
        if k > 1 && x >= ρ[k-1]
            # we're in the lower half
            t = (x - ρ[k]) / hl
            Il = hl^2 * (0.25 * t^4 + (2.0 / 3.0) * t^3 + 0.5 * t^2 - 1.0 / 12.0)
        else
            Il = 0.0
        end
    elseif x >= ρ[k]
        Il = -hl^2 / 12.0
        k == length(ρ) ? hu = 0.0 : hu = ρ[k+1] - ρ[k]
        if k < length(ρ) && x <= ρ[k+1]
            # we're in the upper half
            t = (x - ρ[k]) / hu
            Iu = (hu * t)^2 * (0.25 * t^2 - (2.0 / 3.0) * t + 0.5)
        else
            Iu = hu^2 / 12.0
        end
    end
    return Il + Iu
end

function hermite_coeffs(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    dy_dx = fit_derivative(x, y)
    N = length(x)
    C = Vector{typeof(y[1])}(undef, 2N)
    for i in 1:N
        C[2i-1] = dy_dx[i]
        C[2i] = y[i]
    end
    return C
end

struct FE_rep{S <: AbstractVector{<:Real}, T <: AbstractVector{<:Real}}
    x::S
    coeffs::T
end
FE(x, y) = FE_rep(x, hermite_coeffs(x, y))

function (Y::FE_rep)(x::Real)
    y = 0.0
    for k in 1:length(Y.x)
        y += Y.coeffs[2k-1] * νo(x, k, Y.x)
        y += Y.coeffs[2k] * νe(x, k, Y.x)
    end
    return y
end
function D(Y::FE_rep, x::Real)
    dy_dx = 0.0
    for k in 1:length(Y.x)
        dy_dx += Y.coeffs[2k-1] * D_νo(x, k, Y.x)
        dy_dx += Y.coeffs[2k] * D_νe(x, k, Y.x)
    end
    return dy_dx
end
function I(Y::FE_rep, x::Real)
    yint = 0.0
    for k in 1:length(Y.x)
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