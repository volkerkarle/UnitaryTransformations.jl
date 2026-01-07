@testset "Magnus Expansion" begin
    using QuantumAlgebra
    using UnitaryTransformations: magnus_expansion, FourierHamiltonian, check_hermiticity
    using Symbolics

    @testset "FourierHamiltonian construction" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        modes =
            Dict(0 => Δ / 2 * σz(), 1 => Ω / 4 * (σp() + σm()), -1 => Ω / 4 * (σp() + σm()))

        H = FourierHamiltonian(modes, ω)

        @test H[0] == Δ / 2 * σz()
        @test H[1] == Ω / 4 * (σp() + σm())
        @test H[2] == zero(QuExpr)  # Non-existent mode
        @test haskey(H, 0)
        @test haskey(H, 1)
        @test !haskey(H, 2)

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Hermiticity check" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # Hermitian case: H₋ₙ = Hₙ† for σx = σ⁺ + σ⁻ (self-adjoint)
        modes_hermitian =
            Dict(0 => Δ / 2 * σz(), 1 => Ω / 4 * (σp() + σm()), -1 => Ω / 4 * (σp() + σm()))
        @test check_hermiticity(modes_hermitian) == true

        # Non-Hermitian case: H₁ ≠ H₋₁†
        modes_non_hermitian = Dict(
            0 => Δ / 2 * σz(),
            1 => Ω * σp(),
            -1 => Ω * σp(),  # Should be σm() for Hermiticity
        )
        @test_throws ArgumentError check_hermiticity(modes_non_hermitian)
        @test check_hermiticity(modes_non_hermitian; warn_only = true) == false

        # Correct non-Hermitian fixed
        modes_fixed = Dict(
            0 => Δ / 2 * σz(),
            1 => Ω * σp(),
            -1 => Ω * σm(),  # Now H₋₁ = H₁†
        )
        @test check_hermiticity(modes_fixed) == true

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Magnus order 1 - time average" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # H(t) = Δ/2 σz + Ω cos(ωt) σx
        modes =
            Dict(0 => Δ / 2 * σz(), 1 => Ω / 4 * (σp() + σm()), -1 => Ω / 4 * (σp() + σm()))

        result = magnus_expansion(modes, ω; order = 1)

        # Order 1 is just the time average = H₀
        @test normal_form(result.H_eff) == normal_form(Δ / 2 * σz())
        @test normal_form(result.Ω1) == normal_form(Δ / 2 * σz())
        @test result.Ω2 == zero(QuExpr)

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Magnus order 2 - driven qubit" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # H(t) = Δ/2 σz + Ω cos(ωt) σx
        # cos(ωt) = (e^{iωt} + e^{-iωt})/2
        # So H₁ = H₋₁ = Ω/2 × σx/2 = Ω/4 (σ⁺ + σ⁻)

        modes =
            Dict(0 => Δ / 2 * σz(), 1 => Ω / 4 * (σp() + σm()), -1 => Ω / 4 * (σp() + σm()))

        result = magnus_expansion(modes, ω; order = 2)

        # Second order contribution comes from [H₁, H₋₁]
        # [Ω/4(σ⁺+σ⁻), Ω/4(σ⁺+σ⁻)] = (Ω/4)² × ([σ⁺,σ⁻] + [σ⁻,σ⁺])
        #                            = (Ω²/16) × (σz - σz) = 0
        # So for this symmetric case, Ω₂ = 0!

        @test result.Ω2 == zero(QuExpr)

        # H_eff should just be H₀
        @test normal_form(result.H_eff) == normal_form(Δ / 2 * σz())

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Magnus order 2 - asymmetric drive" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # Circularly polarized drive: H₁ = Ω σ⁺, H₋₁ = Ω σ⁻
        # This has non-zero [H₁, H₋₁] = Ω² [σ⁺, σ⁻] = Ω² σz

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())

        result = magnus_expansion(modes, ω; order = 2)

        # [H₁, H₋₁] = (Ω/2)² [σ⁺, σ⁻] = (Ω²/4) σz
        # Ω₂ coefficient is -1/(nω) = -1/ω for n=1
        # So Ω₂ = -1/ω × (Ω²/4) σz = -Ω²/(4ω) σz

        # H_eff = Δ/2 σz - Ω²/(4ω) σz = (Δ/2 - Ω²/(4ω)) σz
        @test result.H_eff isa QuExpr
        @test !isempty(result.Ω2.terms)

        # In σpm basis, σz = 2σ⁺σ⁻ - 1
        # The result should contain only identity and σ⁺σ⁻ terms (equivalent to σz)
        for (term, coeff) in result.H_eff.terms
            ops = term.bares.v
            # Should be either identity (length 0) or σ⁺σ⁻ (length 2)
            if length(ops) == 2
                # Check it's σ⁺σ⁻
                @test ops[1].t == QuantumAlgebra.TLSCreate_
                @test ops[2].t == QuantumAlgebra.TLSDestroy_
            else
                @test length(ops) == 0  # identity term
            end
        end

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Magnus with bosonic drive" begin
        @variables ω_c Ω ω

        # Driven cavity: H(t) = ω_c a†a + Ω cos(ωt) (a + a†)
        # H₀ = ω_c a†a
        # H₁ = H₋₁ = Ω/2 (a + a†)

        modes = Dict(
            0 => ω_c * a'() * a(),
            1 => Ω / 2 * (a() + a'()),
            -1 => Ω / 2 * (a() + a'()),
        )

        result = magnus_expansion(modes, ω; order = 2)

        # [H₁, H₋₁] = (Ω/2)² [a + a†, a + a†] = 0 (commutes with itself)
        @test result.Ω2 == zero(QuExpr)

        # H_eff = H₀
        @test normal_form(result.H_eff) == normal_form(ω_c * a'() * a())
    end

    @testset "Magnus order 3" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # For third order, we need modes where n + m + k = 0
        # With modes 0, 1, -1: possible combinations are (1, -1, 0), (0, 1, -1), etc.

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())

        result = magnus_expansion(modes, ω; order = 3)

        # Third order should produce some contribution
        # Just check it runs without error and produces valid QuExpr
        @test result.H_eff isa QuExpr
        @test result.Ω3 isa QuExpr

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Magnus expansion - multiple harmonics" begin
        @variables Δ Ω₁ Ω₂ ω

        QuantumAlgebra.use_σpm(true)

        # Drive with fundamental and second harmonic
        # H(t) = Δ/2 σz + Ω₁ cos(ωt) σx + Ω₂ cos(2ωt) σx

        σx_expr = σp() + σm()

        modes = Dict(
            0 => Δ / 2 * σz(),
            1 => Ω₁ / 4 * σx_expr,
            -1 => Ω₁ / 4 * σx_expr,
            2 => Ω₂ / 4 * σx_expr,
            -2 => Ω₂ / 4 * σx_expr,
        )

        result = magnus_expansion(modes, ω; order = 2)

        # Should run without error
        @test result.H_eff isa QuExpr

        # Second order has contributions from [H₁, H₋₁] and [H₂, H₋₂]
        # But σx commutes with itself, so Ω₂ should still be zero
        @test result.Ω2 == zero(QuExpr)

        QuantumAlgebra.use_σpm(false)
    end

    @testset "FourierHamiltonian method" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())

        H = FourierHamiltonian(modes, ω)

        # Both methods should give same result
        result1 = magnus_expansion(modes, ω; order = 2)
        result2 = magnus_expansion(H; order = 2)

        @test normal_form(result1.H_eff) == normal_form(result2.H_eff)
        @test normal_form(result1.Ω2) == normal_form(result2.Ω2)

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Error handling" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # Test that order=0 is rejected
        modes = Dict(0 => Δ / 2 * σz())
        @test_throws ArgumentError magnus_expansion(modes, ω; order = 0)

        # Test that a static Hamiltonian (only mode 0) works correctly
        result = magnus_expansion(modes, ω; order = 3)
        @test result.Ω1 == modes[0]
        @test result.Ω2 == zero(QuExpr)
        @test result.Ω3 == zero(QuExpr)

        # Test that non-Hermitian modes throw an error
        non_hermitian = Dict(0 => Δ / 2 * σz(), 1 => Ω * σp())  # Missing mode -1
        @test_throws ArgumentError magnus_expansion(non_hermitian, ω; order = 2)

        QuantumAlgebra.use_σpm(false)
    end
end
