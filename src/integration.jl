#==============================================================
Define inner products
===============================================================#

function gl_preallocate(N_gl::Integer)
    gξ = zeros(N_gl, N_gl)
    gw = zeros(N_gl, N_gl)
    for i in 1:N_gl
        gξ[1:i, i], gw[1:i, i] = gausslegendre(i)
    end
    return SMatrix{N_gl,N_gl}(gξ), SMatrix{N_gl,N_gl}(gw)
end

const N_gl = 50
const gξ_pa, gw_pa = gl_preallocate(N_gl)

function get_quadrature(x::AbstractVector, order::Integer)
    N = length(x)
    xq = zeros(eltype(x), order * (N - 1))
    wq = similar(xq)
    for j in eachindex(x)
        j == N && break
        dxdξ = 0.5 * (x[j+1] - x[j])
        xavg = 0.5 * (x[j+1] + x[j])
        for k in 1:order
            xq[k+order*(j-1)] = dxdξ * gξ_pa[k, order] + xavg
            wq[k+order*(j-1)] = gw_pa[k, order] * dxdξ
        end
    end
    return xq, wq
end

function integrate(f, lims::SVector; tol::Real=eps(typeof(1.0)))
    return quadgk(f, lims...; atol=tol, rtol=sqrt(tol), maxevals=100)[1]
end

integrate(f, lims::SVector, order::Nothing; tol::Real=eps(typeof(1.0))) = integrate(f, lims; tol)

function integrate(f, lims::SVector{2,<:Real}, order::Integer)
    @assert order <= N_gl
    I = 0.0
    dxdξ = 0.5 * (lims[2] - lims[1])
    xavg = 0.5 * (lims[2] + lims[1])
    for k in 1:order
        I += f(dxdξ * gξ_pa[k, order] + xavg) * gw_pa[k, order] * dxdξ
    end
    return I
end

function integrate(f, lims::SVector{3,<:Real}, order::Integer)
    return integrate(f, SVector(lims[1], lims[2]), order) + integrate(f, SVector(lims[2], lims[3]), order)
end

function dual_integrate(fs, lims::SVector{2,<:Real}, order::Integer)
    @assert order <= N_gl
    I1 = 0.0
    I2 = 0.0
    dxdξ = 0.5 * (lims[2] - lims[1])
    xavg = 0.5 * (lims[2] + lims[1])
    for k in 1:order
        v1, v2 = fs(dxdξ * gξ_pa[k, order] + xavg)
        w = gw_pa[k, order] * dxdξ
        I1 += v1 * w
        I2 += v2 * w
    end
    return I1, I2
end

function dual_integrate(fs, lims::SVector{3,<:Real}, order::Integer)
    hl1 = 0.5 * lims[1]
    hl2 = 0.5 * lims[2]
    hl3 = 0.5 * lims[3]
    dxdξl = hl2 - hl1
    xavgl = hl2 + hl1
    dxdξu = hl3 - hl2
    xavgu = hl3 + hl2
    return dual_integrate(fs, dxdξl, xavgl, dxdξu, xavgu, order)
end

function dual_integrate(fs, dxdξl::Real, xavgl::Real, dxdξu::Real, xavgu::Real, order::Integer)
    @assert order <= N_gl
    I1 = 0.0
    I2 = 0.0
    for k in 1:order
        @inbounds gξ = gξ_pa[k, order]
        @inbounds gw = gw_pa[k, order]

        v1, v2 = fs(muladd(dxdξl, gξ, xavgl))
        w = gw * dxdξl
        I1 += v1 * w
        I2 += v2 * w

        v1, v2 = fs(muladd(dxdξu, gξ, xavgu))
        w = gw * dxdξu
        I1 += v1 * w
        I2 += v2 * w
    end
    return I1, I2
end

function limits(k1::Integer, k2::Integer, ρ::AbstractVector{<:Real})
    k1 != k2 && return SVector(min(ρ[k1], ρ[k2]), max(ρ[k1], ρ[k2]))
    k1 == 1 && return SVector(ρ[1], ρ[2])
    k1 == length(ρ) && return SVector(ρ[end-1], ρ[end])
    return SVector(ρ[k1-1], ρ[k1], ρ[k1+1])
end

function limits(k::Integer, ρ::AbstractVector{<:Real})
    k == 1 && return SVector(ρ[1], ρ[2])
    k == length(ρ) && return SVector(ρ[end-1], ρ[end])
    return SVector(ρ[k-1], ρ[k], ρ[k+1])
end

function inner_product(nu1, k1::Integer, nu2, k2::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing,Integer}=nothing)
    abs(k1 - k2) > 1 && return 0.0
    integrand = x -> nu1(x, k1, ρ) * nu2(x, k2, ρ)
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function inner_product(f, nu1, k1::Integer, nu2, k2::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing,Integer}=nothing)
    abs(k1 - k2) > 1 && return 0.0
    integrand = x -> f(x) * nu1(x, k1, ρ) * nu2(x, k2, ρ)
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function inner_product(nu1, k1::Integer, f, fnu2, g, gnu2, k2::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing,Integer}=nothing)
    abs(k1 - k2) > 1 && return 0.0
    integrand = x -> nu1(x, k1, ρ) * (f(x) * fnu2(x, k2, ρ) + g(x) * gnu2(x, k2, ρ))
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function inner_product(f::Function, nu::Function, k::Integer, ρ::AbstractVector{<:Real}, order::Union{Nothing,Integer}=nothing)
    integrand = x -> f(x) * nu(x, k, ρ)
    return integrate(integrand, limits(k, ρ), order)
end

function inner_product(integrand::Function, k1, k2, ρ::AbstractVector{<:Real}, order::Union{Nothing,Integer}=nothing)
    return integrate(integrand, limits(k1, k2, ρ), order)
end

function dual_inner_product(integrands::Function, k1, k2, ρ::AbstractVector{<:Real}, order::Union{Nothing,Integer}=nothing)
    return dual_integrate(integrands, limits(k1, k2, ρ), order)
end