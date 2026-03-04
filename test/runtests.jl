using FiniteElementHermite
using Test



@testset "even_elements" begin

    ρ = [0.0, 0.33, 0.5, 1.0]

    # Test elements
    for k in eachindex(ρ)
        ρk = ρ[k]
        hl = 0.0

        if k > 1
            ρkl = ρ[k-1]
            hl = ρk - ρkl

            @test νe(ρkl - 0.1, k, ρ) == 0.0
            @test D_νe(ρkl - 0.1, k, ρ) == 0.0
            @test I_νe(ρkl - 0.1, k, ρ) == 0.0

            @test νe(ρkl, k, ρ) == 0.0
            @test D_νe(ρkl, k, ρ) == 0.0
            @test I_νe(ρkl, k, ρ) == 0.0

            @test νe(0.5 * (ρkl + ρk), k, ρ) ≈ 0.5
            @test D_νe(0.5 * (ρkl + ρk), k, ρ) ≈ 1.5 / hl
            @test I_νe(0.5 * (ρkl + ρk), k, ρ) ≈ hl * 3.0 / 32.0
        else
            @test νe(ρk - 0.1, k, ρ) == 0.0
            @test D_νe(ρk - 0.1, k, ρ) == 0.0
            @test I_νe(ρk - 0.1, k, ρ) == 0.0
        end

        @test νe(ρk, k, ρ) == 1.0
        @test D_νe(ρk, k, ρ) == 0.0

        Il = 0.5 * hl

        @test I_νe(ρk, k, ρ) == Il

        if k < length(ρ)
            ρku = ρ[k+1]
            hu = ρku - ρk

            @test νe(0.5 * (ρk + ρku), k, ρ) ≈ 0.5
            @test D_νe(0.5 * (ρk + ρku), k, ρ) ≈ -1.5 / hu
            @test I_νe(0.5 * (ρk + ρku), k, ρ) ≈ Il + hu * 13.0 / 32.0

            Iu = 0.5 * hu

            @test νe(ρku, k, ρ) == 0.0
            @test D_νe(ρku, k, ρ) == 0.0
            @test I_νe(ρku, k, ρ) ≈ Il + Iu

            @test νe(ρku + 0.1, k, ρ) == 0.0
            @test D_νe(ρku + 0.1, k, ρ) == 0.0
            @test I_νe(ρku + 0.1, k, ρ) ≈ Il + Iu
        else
            @test νe(ρk + 0.1, k, ρ) == 0.0
            @test D_νe(ρk + 0.1, k, ρ) == 0.0
            @test I_νe(ρk + 0.1, k, ρ) ≈ Il
        end
    end
end

@testset "odd_elements" begin

    ρ = [0.0, 0.33, 0.5, 1.0]

    for k in eachindex(ρ)
        ρk = ρ[k]
        hl = 0.0

        if k > 1
            ρkl = ρ[k-1]
            hl = ρk - ρkl

            @test νo(ρkl - 0.1, k, ρ) == 0.0
            @test D_νo(ρkl - 0.1, k, ρ) == 0.0
            @test I_νo(ρkl - 0.1, k, ρ) == 0.0

            @test νo(ρkl, k, ρ) == 0.0
            @test D_νo(ρkl, k, ρ) == 0.0
            @test I_νo(ρkl, k, ρ) == 0.0

            @test νo(0.5 * (ρkl + ρk), k, ρ) ≈ -0.125 * hl
            @test D_νo(0.5 * (ρkl + ρk), k, ρ) ≈ -0.25
            @test I_νo(0.5 * (ρkl + ρk), k, ρ) ≈ -(5.0 / 192.0) * hl^2
        else
            @test νo(ρk - 0.1, k, ρ) == 0.0
            @test D_νo(ρk - 0.1, k, ρ) == 0.0
            @test I_νo(ρk - 0.1, k, ρ) == 0.0
        end

        @test νo(ρk, k, ρ) == 0.0
        @test D_νo(ρk, k, ρ) == 1.0

        Il = -hl^2 / 12.0

        @test I_νo(ρk, k, ρ) ≈ Il

        if k < length(ρ)
            ρku = ρ[k+1]
            hu = ρku - ρk

            @test νo(0.5 * (ρk + ρku), k, ρ) ≈ 0.125 * hu
            @test D_νo(0.5 * (ρk + ρku), k, ρ) ≈ -0.25
            @test I_νo(0.5 * (ρk + ρku), k, ρ) ≈ Il + (11.0 / 192.0) * hu^2

            Iu = hu^2 / 12.0

            @test νo(ρku, k, ρ) == 0.0
            @test D_νo(ρku, k, ρ) == 0.0
            @test I_νo(ρku, k, ρ) ≈ Il + Iu

            @test νo(ρku + 0.1, k, ρ) == 0.0
            @test D_νo(ρku + 0.1, k, ρ) == 0.0
            @test I_νo(ρku + 0.1, k, ρ) ≈ Il + Iu
        else
            @test νo(ρk + 0.1, k, ρ) == 0.0
            @test D_νo(ρk + 0.1, k, ρ) == 0.0
            @test I_νo(ρk + 0.1, k, ρ) ≈ Il
        end
    end
end

@testset "FE_representation" begin
    If(x) = x^4 + x^3 + x^2 + x
    f(x) = 4x^3 + 3x^2 + 2x + 1
    Df(x) = 12x^2 + 6 * x + 2

    N = 100

    grid(n) = (3.0 * range(0.0, 1.0, n) .^ 2 .- 1.0)

    ρ = grid(N)

    Ifρ = If.(ρ)
    fρ = f.(ρ)
    Dfρ = Df.(ρ)

    C = zeros(2N)
    C[1:2:end] .= Dfρ
    C[2:2:end] .= fρ
    F = FE_rep(ρ, C)

    x = grid(3N)
    Ifx = If.(x)
    fx = f.(x)
    Dfx = Df.(x)

    # These should be exact for a cubic
    @test all(F.(x) .≈ fx)
    @test all(D.(Ref(F), x) .≈ Dfx)
    @test all(I.(Ref(F), x) .≈ Ifx)

    # These are approximate
    F = FE(ρ, fρ)
    @test all(F.(ρ) .≈ fρ)
    @test all(isapprox.(D.(Ref(F), ρ), Dfρ, rtol=1.0 / N))
    @test all(isapprox.(I.(Ref(F), ρ), Ifρ, rtol=1.0 / N))

end

@testset "second_derivative_elements" begin

    ρ = [0.0, 0.33, 0.5, 1.0]

    for k in eachindex(ρ)
        ρk = ρ[k]
        hl = 0.0

        if k > 1
            ρkl = ρ[k-1]
            hl = ρk - ρkl

            # DD_νe should be 0 outside the element
            @test DD_νe(ρkl - 0.1, k, ρ) == 0.0

            # At boundary ρkl (t=-1), DD_νe = -(12*(-1) + 6)/hl^2 = 6/hl^2
            @test DD_νe(ρkl, k, ρ) ≈ 6.0 / hl^2

            # At midpoint of lower element (t=-0.5), DD_νe = -(12*(-0.5) + 6)/hl^2 = 0
            @test DD_νe(0.5 * (ρkl + ρk), k, ρ) ≈ 0.0 atol=1e-12
        else
            @test DD_νe(ρk - 0.1, k, ρ) == 0.0
        end

        if k < length(ρ)
            ρku = ρ[k+1]
            hu = ρku - ρk

            # At midpoint of upper element (t=0.5), DD_νe = (12*0.5 - 6)/hu^2 = 0
            @test DD_νe(0.5 * (ρk + ρku), k, ρ) ≈ 0.0 atol=1e-12

            # At boundary ρku (t=1), DD_νe = (12*1 - 6)/hu^2 = 6/hu^2
            @test DD_νe(ρku, k, ρ) ≈ 6.0 / hu^2

            # DD_νe should be 0 outside the element
            @test DD_νe(ρku + 0.1, k, ρ) == 0.0
        else
            @test DD_νe(ρk + 0.1, k, ρ) == 0.0
        end
    end

    # Test DD_νo
    for k in eachindex(ρ)
        ρk = ρ[k]
        hl = 0.0

        if k > 1
            ρkl = ρ[k-1]
            hl = ρk - ρkl

            # DD_νo should be 0 outside the element
            @test DD_νo(ρkl - 0.1, k, ρ) == 0.0

            # At boundary ρkl (t=-1), DD_νo = (6*(-1) + 4)/hl = -2/hl
            @test DD_νo(ρkl, k, ρ) ≈ -2.0 / hl

            # At midpoint of lower element (t=-0.5), DD_νo = (6*(-0.5) + 4)/hl = 1/hl
            @test DD_νo(0.5 * (ρkl + ρk), k, ρ) ≈ 1.0 / hl
        else
            @test DD_νo(ρk - 0.1, k, ρ) == 0.0
        end

        if k < length(ρ)
            ρku = ρ[k+1]
            hu = ρku - ρk

            # At midpoint of upper element (t=0.5), DD_νo = (6*0.5 - 4)/hu = -1/hu
            @test DD_νo(0.5 * (ρk + ρku), k, ρ) ≈ -1.0 / hu

            # At boundary ρku (t=1), DD_νo = (6*1 - 4)/hu = 2/hu
            @test DD_νo(ρku, k, ρ) ≈ 2.0 / hu

            # DD_νo should be 0 outside the element
            @test DD_νo(ρku + 0.1, k, ρ) == 0.0
        else
            @test DD_νo(ρk + 0.1, k, ρ) == 0.0
        end
    end
end

@testset "DD_FE_representation" begin
    # Use a cubic polynomial: f(x) = 4x³ + 3x² + 2x + 1
    # First derivative: Df(x) = 12x² + 6x + 2
    # Second derivative: DDf(x) = 24x + 6
    f(x) = 4x^3 + 3x^2 + 2x + 1
    Df(x) = 12x^2 + 6 * x + 2
    DDf(x) = 24x + 6

    N = 100
    grid(n) = (3.0 * range(0.0, 1.0, n) .^ 2 .- 1.0)

    ρ = grid(N)
    fρ = f.(ρ)
    Dfρ = Df.(ρ)

    C = zeros(2N)
    C[1:2:end] .= Dfρ
    C[2:2:end] .= fρ
    F = FE_rep(ρ, C)

    # Test at interior points (not exactly on grid points to avoid discontinuities)
    x_test = grid(3N)
    # Filter to interior points only (avoid element boundaries where DD may be discontinuous)
    x_interior = filter(xi -> !any(isapprox.(xi, ρ, atol=1e-10)), x_test)

    DDfx = DDf.(x_interior)

    # DD should be exact for a cubic polynomial
    @test all(DD.(Ref(F), x_interior) .≈ DDfx)
end

@testset "extrapolation" begin
    # Use a quadratic polynomial: f(x) = ax² + bx + c
    # Extrapolation should be exact for a quadratic
    a, b, c = 2.0, -3.0, 5.0
    f(x) = a * x^2 + b * x + c
    Df(x) = 2a * x + b
    DDf(x) = 2a

    N = 50
    ρ = range(0.0, 1.0, N)
    fρ = f.(ρ)
    Dfρ = Df.(ρ)

    C = zeros(2N)
    C[1:2:end] .= Dfρ
    C[2:2:end] .= fρ
    F = FE_rep(collect(ρ), C)

    # Test low boundary extrapolation
    x_low = [-0.5, -0.2, -0.1, -0.01]
    for x in x_low
        @test extrapolate(F, x) ≈ f(x) atol=1e-10
    end

    # Test high boundary extrapolation
    x_high = [1.01, 1.1, 1.5, 2.0]
    for x in x_high
        @test extrapolate(F, x) ≈ f(x) atol=1e-10
    end

    # Test continuity at boundaries
    # Value should match at boundary
    @test extrapolate(F, ρ[1]) ≈ F(ρ[1]) atol=1e-10
    @test extrapolate(F, ρ[end]) ≈ F(ρ[end]) atol=1e-10

    # Test multiple FE_reps with shared grid (efficient extrapolation)
    g(x) = -x^2 + 4x - 1
    Dg(x) = -2x + 4
    gρ = g.(ρ)
    Dgρ = Dg.(ρ)

    C2 = zeros(2N)
    C2[1:2:end] .= Dgρ
    C2[2:2:end] .= gρ
    G = FE_rep(collect(ρ), C2)

    # Pre-compute extrapolation bases
    x_test = -0.3
    side, k, Δx, DD_bases... = compute_extrapolation_bases(F.x, x_test)

    # Extrapolate both FE_reps using shared bases
    @test extrapolate(F, side, k, Δx, DD_bases...) ≈ f(x_test) atol=1e-10
    @test extrapolate(G, side, k, Δx, DD_bases...) ≈ g(x_test) atol=1e-10

    # Fix #1: in-bounds x should throw
    @test_throws ArgumentError compute_extrapolation_bases(collect(ρ), 0.5)
    @test_throws ArgumentError extrapolate(F, 0.5)

    # Fix #2: invalid order should throw
    @test_throws ArgumentError extrapolate(F, -0.1; order=0)
    @test_throws ArgumentError extrapolate(F, -0.1; order=3)
end
