#==============================================================
Define Hermite cubic finite elements
===============================================================#

"""
    quad_deriv(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Real}

Fit a parabola going through three points `(x1, y1)`, `(x2, y2)`, and `(x3, y3)`
and return derivative at `x2`
"""
function quad_deriv(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Real}
    hl = x2 - x1
    hu = x3 - x2
    return ((y3 - y2) * hl^2 + (y2 - y1) * hu^2) / (hu * hl * (hu + hl))
end

"""
    fit_derivative(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Returns dy/dx at every point in `x`, based on local quadratic fit
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

"""
    νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})

Returns the value of the even basis element centered at `ρ[k]` at `x`
"""
function νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return νel(x, ρ1, ρ2)
    region === :in_up && return νeu(x, ρ1, ρ2)
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

"""
    D_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})

Returns the derivative of the even basis element centered at `ρ[k]` at `x`
"""
function D_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return D_νel(x, ρ1, ρ2)
    region === :in_up && return D_νeu(x, ρ1, ρ2)
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

"""
    I_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})

Returns the integral of the even basis element centered at `ρ[k]` from 0 to `x`
"""
function I_νe(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    ρk = ρ[k]

    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :out_low && return 0.0
    region === :in_low && return I_νel(x, ρ1, ρ2)

    Il = k > 1 ? 0.5 * (ρk - ρ[k-1]) : 0.0
    region === :in_up && return Il + I_νeu(x, ρ1, ρ2)
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

"""
    νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})

Returns the value of the odd basis element centered at `ρ[k]` at `x`
"""
function νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return νol(x, ρ1, ρ2)
    region === :in_up && return νou(x, ρ1, ρ2)
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

"""
    D_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})

Returns the derivative of the odd basis element centered at `ρ[k]` at `x`
"""
function D_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :in_low && return D_νol(x, ρ1, ρ2)
    region === :in_up && return D_νou(x, ρ1, ρ2)
    return 0.0
end

const one_twelf = 1.0 / 12.0
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

"""
    I_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})

Returns the integral of the odd basis element centered at `ρ[k]` from 0 to `x`
"""
function I_νo(x::Real, k::Integer, ρ::AbstractVector{<:Real})
    ρk = ρ[k]

    region, ρ1, ρ2 = which_region(x, k, ρ)
    region === :out_low && return 0.0
    region === :in_low && return I_νol(x, ρ1, ρ2)

    hl = k > 1 ? ρk - ρ[k-1] : 0.0
    Il = -one_twelf * hl^2
    region === :in_up && return Il + I_νou(x, ρ1, ρ2)

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

struct FE_rep{S<:AbstractVector{<:Real},T<:AbstractVector{<:Real}}
    x::S
    coeffs::T
    function FE_rep{S,T}(x::S, coeffs::T) where {S<:AbstractVector{<:Real},T<:AbstractVector{<:Real}}
        return length(coeffs) == 2length(x) ? new{S,T}(x, coeffs) : throw(DimensionMismatch)
    end
end

"""
    FE_rep(x::S, coeffs::T) where {S<:AbstractVector{<:Real},T<:AbstractVector{<:Real}}

Create Hermite-cubic finite-element representation on grid `x` with cofficients `coeffs`
`coeffs[2k]` is the value of the function at `x[k]`
`coeffs[2k-1]` is the derivative of the function at `x[k]`

An `FE_rep` is callable, giving the value of the finite-element representation, for example:
```
Y = FE_rep(x, coeffs)
Y(a) # gives value of the finite-element representation at x=a
```
"""
function FE_rep(x::S, coeffs::T) where {S<:AbstractVector{<:Real},T<:AbstractVector{<:Real}}
    return FE_rep{S,T}(x, coeffs)
end

"""
    FE(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Returns finite element representation of data `(x, y)`, fitting local quadratics to determine derivative
"""
function FE(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    return FE_rep(x, hermite_coeffs(x, y))
end

"""
    FE(grid::AbstractVector{<:Real}, xy::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}})

Returns finite element representation of data `(x, y)`, mapped to `grid`
"""
function FE(grid::AbstractVector{<:Real}, xy::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}})
    f = FE(xy...)
    C = Vector{eltype(grid)}(undef, 2 * length(grid))
    @inbounds for (i, g) in enumerate(grid)
        ti = 2i
        C[ti-1] = D(f, g)
        C[ti] = f(g)
    end
    return FE_rep(grid, C)
end

"""
    compute_bases(X::AbstractVector{<:Real}, x::Real)

For grid `X`, find the value of the four basis functions at `x`
This can be used with `evaluate` or `evaluate_inbounds` to compute the value of multiple
  `FE_rep`s efficiently if they share the same grid

Returns `(k, nu_ou, nu_eu, nu_ol, nu_el)`
where `nu_ou` is odd  basis element centered at `X[k]`
      `nu_eu` is even basis element centered at `X[k]`
      `nu_ol` is odd  basis element centered at `X[k+1]`
      `nu_el` is even basis element centered at `X[k+1]`
"""
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

"""
    evaluate(Y::FE_rep, k::Integer, nu_ou::Real, nu_eu::Real, nu_ol::Real, nu_el::Real)

Evaluate `FE_rep` `Y` assuming basis function values/derivatives have been computed by
    `compute_bases`, `compute_D_bases`, or `compute_both_bases`

    `nu_ou` is value or derivative of odd  basis element centered at `X[k]`
    `nu_eu` is value or derivative of even basis element centered at `X[k]`
    `nu_ol` is value or derivative of odd  basis element centered at `X[k+1]`
    `nu_el` is value or derivative of even basis element centered at `X[k+1]`
"""
@inline function evaluate(Y::FE_rep, k::Integer, nu_ou::Real, nu_eu::Real, nu_ol::Real, nu_el::Real)
    tk = 2k
    y = Y.coeffs[tk-1] * nu_ou
    y += Y.coeffs[tk] * nu_eu
    y += Y.coeffs[tk+1] * nu_ol
    y += Y.coeffs[tk+2] * nu_el
    return y
end

"""
    evaluate_inbounds(Y::FE_rep, k::Integer, nu_ou::T, nu_eu::T, nu_ol::T, nu_el::T) where {T<:Real}

Evaluate `FE_rep` `Y` without bounds checking, assuming basis function values/derivatives
     have been computed by `compute_bases`, `compute_D_bases`, or `compute_both_bases`

    `nu_ou` is value or derivative of odd  basis element centered at `X[k]`
    `nu_eu` is value or derivative of even basis element centered at `X[k]`
    `nu_ol` is value or derivative of odd  basis element centered at `X[k+1]`
    `nu_el` is value or derivative of even basis element centered at `X[k+1]`
"""
@inline function evaluate_inbounds(Y::FE_rep, k::Integer, nu_ou::T, nu_eu::T, nu_ol::T, nu_el::T) where {T<:Real}
    tk = 2k
    @inbounds y = Y.coeffs[tk-1] * nu_ou + Y.coeffs[tk] * nu_eu
    @inbounds y += Y.coeffs[tk+1] * nu_ol + Y.coeffs[tk+2] * nu_el
    return y
end

"""
    (Y::FE_rep)(x::Real)

Functor for `FE_rep`, giving the value of the finite-element representation at `x`
"""
function (Y::FE_rep)(x::Real)
    k, nu_ou, nu_eu, nu_ol, nu_el = compute_bases(Y.x, x)
    y = evaluate_inbounds(Y, k, nu_ou, nu_eu, nu_ol, nu_el)
    return y
end

"""
    compute_D_bases(X::AbstractVector{<:Real}, x::Real)

For grid `X`, find the derivative of the four basis functions at `x`
This can be used with `evaluate` or `evaluate_inbounds` to compute the derivative of multiple
  `FE_rep`s efficiently if they share the same grid

Returns `(k, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el)`
where `D_nu_ou` is derivative of odd  basis element centered at `X[k]`
      `D_nu_eu` is derivative of even basis element centered at `X[k]`
      `D_nu_ol` is derivative of odd  basis element centered at `X[k+1]`
      `D_nu_el` is derivative of even basis element centered at `X[k+1]`
"""
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

"""
    compute_both_bases(X::AbstractVector{<:Real}, x::Real)

For grid `X`, find the value and derivative of the four basis functions at `x`
This can be used with `evaluate` or `evaluate_inbounds` to compute the value and derivative of multiple
  `FE_rep` efficiently if they share the same grid

Returns `(k, nu_ou, nu_eu, nu_ol, nu_el, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el)`
where `nu_ou` is odd  basis element centered at `X[k]`
      `nu_eu` is even basis element centered at `X[k]`
      `nu_ol` is odd  basis element centered at `X[k+1]`
      `nu_el` is even basis element centered at `X[k+1]`
      `D_nu_ou` is derivative of odd  basis element centered at `X[k]`
      `D_nu_eu` is derivative of even basis element centered at `X[k]`
      `D_nu_ol` is derivative of odd  basis element centered at `X[k+1]`
      `D_nu_el` is derivative of even basis element centered at `X[k+1]`
"""
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

"""
    D(Y::FE_rep, x::Real)

Return derivative of FE_rep `Y` at location `x`
"""
function D(Y::FE_rep, x::Real)
    k, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el = compute_D_bases(Y.x, x)
    dy_dx = evaluate_inbounds(Y, k, D_nu_ou, D_nu_eu, D_nu_ol, D_nu_el)
    return dy_dx
end

"""
    I(Y::FE_rep, x::Real)

Return integral of FE_rep `Y` from 0 to `x`
"""
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