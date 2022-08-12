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
            
            @test νe(0.5*(ρkl + ρk), k, ρ) ≈ 0.5
            @test D_νe(0.5*(ρkl + ρk), k, ρ) ≈ 1.5/hl
            @test I_νe(0.5*(ρkl + ρk), k, ρ) ≈ hl * 3.0 / 32.0
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
            
            @test νe(0.5*(ρk + ρku), k, ρ) ≈ 0.5
            @test D_νe(0.5*(ρk + ρku), k, ρ) ≈ -1.5/hu
            @test I_νe(0.5*(ρk + ρku), k, ρ) ≈ Il + hu * 13.0 / 32.0
            
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
            
            @test νo(0.5*(ρkl + ρk), k, ρ) ≈ -0.125 * hl
            @test D_νo(0.5*(ρkl + ρk), k, ρ) ≈ -0.25
            @test I_νo(0.5*(ρkl + ρk), k, ρ) ≈ -(5.0/192.0) * hl^2 
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
            
            @test νo(0.5*(ρk + ρku), k, ρ) ≈ 0.125 * hu
            @test D_νo(0.5*(ρk + ρku), k, ρ) ≈ -0.25
            @test I_νo(0.5*(ρk + ρku), k, ρ) ≈ Il + (11.0 / 192.0) * hu^2
            
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
    f(x)  = 4x^3 + 3x^2 + 2x + 1
    Df(x) = 12x^2 + 6*x + 2

    N = 100

    grid(n) = (3.0 * range(0.0, 1.0, n).^ 2 .- 1.0)
    
    ρ = grid(N)

    Ifρ = If.(ρ)
    fρ  = f.(ρ)
    Dfρ = Df.(ρ)

    C = zeros(2N)
    C[1:2:end] .= Dfρ
    C[2:2:end] .= fρ
    F = FE_rep(ρ, C)

    x = grid(3N)
    Ifx = If.(x)
    fx  = f.(x)
    Dfx = Df.(x)

    # These should be exact for a cubic
    @test all(F.(x) .≈ fx)
    @test all(D.(Ref(F), x) .≈ Dfx)
    @test all(I.(Ref(F), x) .≈ Ifx)

    # These are approximate 
    F = FE(ρ, fρ)
    @test all(F.(ρ) .≈ fρ)
    @test all(isapprox.(D.(Ref(F), ρ), Dfρ, rtol = 1.0/N))
    @test all(isapprox.(I.(Ref(F), ρ), Ifρ, rtol = 1.0/N))
    
end