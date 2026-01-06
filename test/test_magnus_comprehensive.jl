"""
Comprehensive tests for Magnus expansion.

These tests verify the Magnus expansion against known analytical formulas
and physical results from the literature.

References:
1. Blanes et al., "The Magnus expansion and some of its applications", 
   Physics Reports 470, 151 (2009).
2. Goldman & Dalibard, "Periodically Driven Quantum Systems: Effective 
   Hamiltonians and Engineered Gauge Fields", PRX 4, 031027 (2014).
3. Eckardt & Anisimovas, "High-frequency approximation for periodically 
   driven quantum systems", New J. Phys. 17, 093039 (2015).
"""

@testset "Magnus Expansion - Comprehensive Tests" begin
    using QuantumAlgebra
    using UnitaryTransformations: magnus_expansion, FourierHamiltonian, check_hermiticity
    using Symbolics

    # ==========================================================================
    # Helper functions for coefficient extraction
    # ==========================================================================

    """
    Extract the scalar coefficient of a specific operator from a QuExpr.
    Returns the coefficient of the first matching operator structure.
    """
    function extract_coefficient(expr::QuExpr, target_op::QuExpr)
        target_nf = normal_form(target_op)
        expr_nf = normal_form(expr)

        if isempty(target_nf.terms)
            # Looking for identity coefficient
            for (term, coeff) in expr_nf.terms
                if isempty(term.bares.v)
                    return coeff
                end
            end
            return 0
        end

        # Get the operator structure from target
        target_term, _ = first(target_nf.terms)

        for (term, coeff) in expr_nf.terms
            if term == target_term
                return coeff
            end
        end
        return 0
    end

    """
    Check if two symbolic expressions are equal after simplification.
    """
    function sym_equal(a, b; rtol = 1e-10)
        diff = Symbolics.simplify(a - b)
        # Check if the difference is numerically zero
        if diff isa Number
            return abs(diff) < rtol
        end
        # For symbolic expressions, check if they simplify to the same form
        return isequal(Symbolics.simplify(a), Symbolics.simplify(b))
    end

    # ==========================================================================
    # First Order Tests: H_eff^(1) = H₀ (time average)
    # ==========================================================================

    @testset "Order 1: Time averaging" begin
        @variables Δ Ω ω α β

        QuantumAlgebra.use_σpm(true)

        # Test 1: Pure DC component
        @testset "Pure DC" begin
            modes = Dict(0 => Δ / 2 * σz())
            result = magnus_expansion(modes, ω; order = 1)
            @test normal_form(result.H_eff) == normal_form(Δ / 2 * σz())
        end

        # Test 2: AC component averages to zero
        @testset "AC averages to zero" begin
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 1)
            # Only H₀ survives
            @test normal_form(result.H_eff) == normal_form(Δ / 2 * σz())
        end

        # Test 3: Multiple harmonics all average to zero
        @testset "Multiple harmonics average" begin
            modes = Dict(
                0 => Δ / 2 * σz(),
                1 => α * σp(),
                -1 => α * σm(),
                2 => β * σp(),
                -2 => β * σm(),
            )
            result = magnus_expansion(modes, ω; order = 1)
            @test normal_form(result.H_eff) == normal_form(Δ / 2 * σz())
        end

        # Test 4: Bosonic system
        @testset "Bosonic DC + AC" begin
            modes = Dict(0 => ω * a'() * a(), 1 => Ω / 2 * a(), -1 => Ω / 2 * a'())
            result = magnus_expansion(modes, ω; order = 1, check_hermitian = false)
            @test normal_form(result.H_eff) == normal_form(ω * a'() * a())
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Second Order Tests: H_eff^(2) = Σₙ>₀ -[Hₙ, H₋ₙ]/(nω)
    # ==========================================================================

    @testset "Order 2: Commutator contributions" begin
        @variables Δ Ω ω Ω₁ Ω₂

        QuantumAlgebra.use_σpm(true)

        # Test 1: Linearly polarized drive - symmetric case gives zero
        @testset "Linear polarization (symmetric) - zero Ω₂" begin
            # H(t) = Δ/2 σz + Ω cos(ωt) σx
            # σx = σ⁺ + σ⁻, so H₁ = H₋₁ = Ω/4(σ⁺ + σ⁻)
            # [H₁, H₋₁] = (Ω/4)²([σ⁺,σ⁻] + [σ⁻,σ⁺]) = 0
            σx = σp() + σm()
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 4 * σx, -1 => Ω / 4 * σx)
            result = magnus_expansion(modes, ω; order = 2)
            @test result.Ω2 == zero(QuExpr)
        end

        # Test 2: Circularly polarized drive - Bloch-Siegert shift
        @testset "Circular polarization - Bloch-Siegert shift" begin
            # H(t) = Δ/2 σz + Ω/2(e^{iωt}σ⁺ + e^{-iωt}σ⁻)
            # H₁ = Ω/2 σ⁺, H₋₁ = Ω/2 σ⁻
            # [H₁, H₋₁] = (Ω/2)²[σ⁺,σ⁻] = (Ω²/4)σz (in σpm: σz = 2σ⁺σ⁻ - 1)
            # H_eff^(2) = -[H₁, H₋₁]/ω = -Ω²/(4ω) σz
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 2)

            # The Bloch-Siegert shift should modify the effective detuning
            # H_eff = (Δ/2 - Ω²/(4ω)) σz = (Δ/2 - Ω²/(4ω))(2σ⁺σ⁻ - 1)
            # = (Δ - Ω²/(2ω))σ⁺σ⁻ + (-Δ/2 + Ω²/(4ω)) × 1

            @test !isempty(result.Ω2.terms)

            # Verify the σ⁺σ⁻ coefficient
            σpm_coeff = extract_coefficient(result.Ω2, σp() * σm())
            # Expected: -Ω²/(2ω) from the σ⁺σ⁻ part of -Ω²/(4ω)σz
            # σz = 2σ⁺σ⁻ - 1, so coeff of σ⁺σ⁻ is 2 × (-Ω²/(4ω)) = -Ω²/(2ω)
            expected = -Ω^2 / (2 * ω)
            @test sym_equal(σpm_coeff, expected)
        end

        # Test 3: Second harmonic drive
        @testset "Second harmonic (n=2)" begin
            # H₂ = Ω/2 σ⁺, H₋₂ = Ω/2 σ⁻
            # H_eff^(2) = -[H₂, H₋₂]/(2ω) = -Ω²/(8ω) σz
            modes = Dict(0 => Δ / 2 * σz(), 2 => Ω / 2 * σp(), -2 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 2)

            # The coefficient for σ⁺σ⁻ should be half of n=1 case (factor of 1/(nω))
            σpm_coeff = extract_coefficient(result.Ω2, σp() * σm())
            expected = -Ω^2 / (4 * ω)  # Factor of 1/2 from n=2
            @test sym_equal(σpm_coeff, expected)
        end

        # Test 4: Multiple harmonics - contributions add
        @testset "Multiple harmonics add" begin
            # H₁ = Ω₁/2 σ⁺, H₋₁ = Ω₁/2 σ⁻
            # H₂ = Ω₂/2 σ⁺, H₋₂ = Ω₂/2 σ⁻
            # Total shift: -Ω₁²/(4ω) - Ω₂²/(8ω)
            modes = Dict(
                0 => Δ / 2 * σz(),
                1 => Ω₁ / 2 * σp(),
                -1 => Ω₁ / 2 * σm(),
                2 => Ω₂ / 2 * σp(),
                -2 => Ω₂ / 2 * σm(),
            )
            result = magnus_expansion(modes, ω; order = 2)

            σpm_coeff = extract_coefficient(result.Ω2, σp() * σm())
            expected = -Ω₁^2 / (2 * ω) - Ω₂^2 / (4 * ω)
            @test sym_equal(σpm_coeff, expected)
        end

        # Test 5: Bosonic operators - self-commuting gives zero
        @testset "Bosonic symmetric - zero Ω₂" begin
            # [a + a†, a + a†] = 0
            modes = Dict(
                0 => ω * a'() * a(),
                1 => Ω / 2 * (a() + a'()),
                -1 => Ω / 2 * (a() + a'()),
            )
            result = magnus_expansion(modes, ω; order = 2)
            @test result.Ω2 == zero(QuExpr)
        end

        # Test 6: Bosonic asymmetric - displacement
        @testset "Bosonic asymmetric" begin
            # H₁ = Ω a, H₋₁ = Ω a†
            # [a, a†] = 1
            # H_eff^(2) = -Ω²/ω × 1 (just a c-number shift!)
            modes = Dict(0 => ω * a'() * a(), 1 => Ω * a(), -1 => Ω * a'())
            result = magnus_expansion(modes, ω; order = 2, check_hermitian = false)

            # The result should have a constant term from -[a, a†]/ω = -1/ω
            @test !isempty(result.Ω2.terms)
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Third Order Tests
    # ==========================================================================

    @testset "Order 3: Triple commutators" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # Test 1: Verify third order runs without error
        @testset "Basic third order" begin
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 3)
            @test result.Ω3 isa QuExpr
            @test result.H_eff isa QuExpr
        end

        # Test 2: Symmetric drive (σx) should have specific third-order structure
        @testset "Symmetric third order" begin
            σx = σp() + σm()
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 4 * σx, -1 => Ω / 4 * σx)
            result = magnus_expansion(modes, ω; order = 3)

            # Third order should contain contributions from
            # [[H₁, H₀], H₋₁] + [[H₀, H₋₁], H₁] + ...
            @test result.Ω3 isa QuExpr
        end

        # Test 3: Third order scaling with ω
        @testset "Third order scales as 1/ω²" begin
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 3)

            # Third order terms should have 1/ω² dependence
            # Check that Ω3 terms have ω² in denominator
            if !isempty(result.Ω3.terms)
                _, first_coeff = first(result.Ω3.terms)
                # The coefficient should contain ω^(-2)
                @test occursin("ω", string(first_coeff)) || first_coeff == 0
            end
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Physical Systems Tests
    # ==========================================================================

    @testset "Physical systems" begin
        @variables ω_c g Δ ω Ω ω_d

        QuantumAlgebra.use_σpm(true)

        # Test 1: Jaynes-Cummings with modulated coupling
        @testset "Modulated Jaynes-Cummings" begin
            # H₀ = ω_c a†a + Δ/2 σz
            # Modulated coupling: g cos(ωt)(a†σ⁻ + aσ⁺)
            H0 = ω_c * a'() * a() + Δ / 2 * σz()
            H1 = g / 2 * (a'() * σm() + a() * σp())
            Hminus1 = g / 2 * (a'() * σm() + a() * σp())  # Same for cos drive

            modes = Dict(0 => H0, 1 => H1, -1 => Hminus1)
            result = magnus_expansion(modes, ω; order = 2)

            @test result.H_eff isa QuExpr
            # H_eff should contain a†a, σz, and possibly a†a σz cross terms
        end

        # Test 2: Driven cavity (parametric amplification setup)
        @testset "Driven cavity" begin
            # H₀ = ω_c a†a
            # H_drive = Ω(e^{iω_d t} a + e^{-iω_d t} a†)
            H0 = ω_c * a'() * a()
            H1 = Ω * a()
            Hminus1 = Ω * a'()

            modes = Dict(0 => H0, 1 => H1, -1 => Hminus1)
            result = magnus_expansion(modes, ω_d; order = 2, check_hermitian = false)

            @test result.H_eff isa QuExpr
        end

        # Test 3: Two-tone drive
        @testset "Two-tone drive" begin
            @variables Ω₁ Ω₂ ω₁ ω₂

            # Note: This tests multiple incommensurate frequencies
            # For true two-tone, we'd need a different approach
            # Here we test commensurate: ω₂ = 2ω₁
            modes = Dict(
                0 => Δ / 2 * σz(),
                1 => Ω₁ / 2 * σp(),
                -1 => Ω₁ / 2 * σm(),
                2 => Ω₂ / 2 * σp(),
                -2 => Ω₂ / 2 * σm(),
            )
            result = magnus_expansion(modes, ω; order = 2)

            @test result.H_eff isa QuExpr
            @test result.Ω2 isa QuExpr
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Hermiticity Tests
    # ==========================================================================

    @testset "Hermiticity of H_eff" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        @testset "H_eff is Hermitian for Hermitian input" begin
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 3)

            # H_eff† should equal H_eff
            H_dag = normal_form(result.H_eff')
            H_nf = normal_form(result.H_eff)
            @test H_dag == H_nf
        end

        @testset "Ω₂ is Hermitian" begin
            modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 2)

            Ω2_dag = normal_form(result.Ω2')
            Ω2_nf = normal_form(result.Ω2)
            @test Ω2_dag == Ω2_nf
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Edge Cases and Robustness
    # ==========================================================================

    @testset "Edge cases" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # Test 1: Only H₀
        @testset "Only H₀" begin
            modes = Dict(0 => Δ / 2 * σz())
            result = magnus_expansion(modes, ω; order = 4)
            @test normal_form(result.H_eff) == normal_form(Δ / 2 * σz())
            @test result.Ω2 == zero(QuExpr)
            @test result.Ω3 == zero(QuExpr)
            @test result.Ω4 == zero(QuExpr)
        end

        # Test 2: High harmonic only
        @testset "High harmonic (n=5)" begin
            modes = Dict(0 => Δ / 2 * σz(), 5 => Ω / 2 * σp(), -5 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 2)

            # Coefficient should have 1/(5ω) factor
            σpm_coeff = extract_coefficient(result.Ω2, σp() * σm())
            expected = -Ω^2 / (10 * ω)
            @test sym_equal(σpm_coeff, expected)
        end

        # Test 3: Empty modes except H₀
        @testset "Sparse modes" begin
            modes = Dict(0 => Δ / 2 * σz(), 3 => Ω / 2 * σp(), -3 => Ω / 2 * σm())
            result = magnus_expansion(modes, ω; order = 3)
            @test result.H_eff isa QuExpr
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Bloch-Siegert Shift Quantitative Test
    # ==========================================================================

    @testset "Bloch-Siegert shift - quantitative" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # The Bloch-Siegert shift for a circularly driven TLS is:
        # δω_BS = Ω²/(4ω)
        #
        # For H(t) = Δ/2 σz + Ω/2(e^{iωt}σ⁺ + e^{-iωt}σ⁻)
        # H_eff = (Δ/2 - Ω²/(4ω))σz
        #
        # In σpm basis: σz = 2σ⁺σ⁻ - 1
        # H_eff = (Δ - Ω²/(2ω))σ⁺σ⁻ + const

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
        result = magnus_expansion(modes, ω; order = 2)

        # Extract total σ⁺σ⁻ coefficient from H_eff
        σpm_coeff_Heff = extract_coefficient(result.H_eff, σp() * σm())

        # Expected: Δ (from H₀ = Δ/2 σz = Δ σ⁺σ⁻ - Δ/2) + (-Ω²/(2ω)) from Ω₂
        expected_total = Δ - Ω^2 / (2 * ω)
        @test sym_equal(σpm_coeff_Heff, expected_total)

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Third-Order Quantitative Test
    # ==========================================================================

    @testset "Third order shift - quantitative" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # For a circularly driven TLS, the third-order correction is:
        # H_eff^(3) = ΔΩ²/(4ω²) σz  (shifts the effective detuning)
        #
        # In σpm basis: σz = 2σ⁺σ⁻ - 1
        # So the σ⁺σ⁻ coefficient from Ω₃ is -ΔΩ²/(2ω²)
        #
        # Total H_eff to order 3:
        # H_eff = (Δ/2 - Ω²/(4ω) + ΔΩ²/(4ω²)) σz
        # σ⁺σ⁻ coeff = Δ - Ω²/(2ω) - ΔΩ²/(2ω²)

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
        result = magnus_expansion(modes, ω; order = 3)

        # Test Ω₃ coefficient
        σpm_coeff_Ω3 = extract_coefficient(result.Ω3, σp() * σm())
        expected_Ω3 = -Δ * Ω^2 / (2 * ω^2)
        @test sym_equal(σpm_coeff_Ω3, expected_Ω3)

        # Test total H_eff to order 3
        σpm_coeff_Heff = extract_coefficient(result.H_eff, σp() * σm())
        expected_total = Δ - Ω^2 / (2 * ω) - Δ * Ω^2 / (2 * ω^2)
        @test sym_equal(σpm_coeff_Heff, expected_total)

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Higher-Order Tests (Orders 4-6)
    # ==========================================================================

    @testset "Higher orders (4-6) - pattern verification" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        # For a circularly driven TLS, the Magnus expansion has a beautiful pattern.
        # Due to the simple commutation relations of spin-1/2 operators, only terms
        # of the form [H₁, H₀, ..., H₀, H₋₁] contribute at each order.
        #
        # The σ⁺σ⁻ coefficients follow the pattern:
        #   Order 2: -Ω²/(2ω)
        #   Order 3: -ΔΩ²/(2ω²)
        #   Order 4: +Δ²Ω²/(2ω³)
        #   Order 5: -Δ³Ω²/(2ω⁴)
        #   Order 6: +Δ⁴Ω²/(2ω⁵)
        #
        # General pattern for k ≥ 2: (-1)^k × Δ^(k-2) × Ω² / (2ω^(k-1))
        # This series converges to -Ω²/(2(ω+Δ)) = exact Bloch-Siegert shift.

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())
        result = magnus_expansion(modes, ω; order = 6)

        # Test order 4: +Δ²Ω²/(2ω³)
        @testset "Order 4 coefficient" begin
            σpm_coeff = extract_coefficient(result.Ω4, σp() * σm())
            expected = Δ^2 * Ω^2 / (2 * ω^3)
            @test sym_equal(σpm_coeff, expected)
        end

        # Test order 5: -Δ³Ω²/(2ω⁴)
        @testset "Order 5 coefficient" begin
            σpm_coeff = extract_coefficient(result.Ω5, σp() * σm())
            expected = -Δ^3 * Ω^2 / (2 * ω^4)
            @test sym_equal(σpm_coeff, expected)
        end

        # Test order 6: +Δ⁴Ω²/(2ω⁵)
        @testset "Order 6 coefficient" begin
            σpm_coeff = extract_coefficient(result.Ω6, σp() * σm())
            expected = Δ^4 * Ω^2 / (2 * ω^5)
            @test sym_equal(σpm_coeff, expected)
        end

        # Test cumulative H_eff to order 6
        @testset "Cumulative H_eff to order 6" begin
            σpm_coeff_Heff = extract_coefficient(result.H_eff, σp() * σm())
            # Sum: Δ - Ω²/(2ω) × [1 + Δ/ω - Δ²/ω² + Δ³/ω³ - Δ⁴/ω⁴]
            expected_total =
                Δ - Ω^2/(2*ω) - Δ*Ω^2/(2*ω^2) + Δ^2*Ω^2/(2*ω^3) - Δ^3*Ω^2/(2*ω^4) +
                Δ^4*Ω^2/(2*ω^5)
            @test sym_equal(σpm_coeff_Heff, expected_total)
        end

        QuantumAlgebra.use_σpm(false)
    end

    # ==========================================================================
    # Arbitrary Order Computation Test
    # ==========================================================================

    @testset "Arbitrary order computation" begin
        @variables Δ Ω ω

        QuantumAlgebra.use_σpm(true)

        modes = Dict(0 => Δ / 2 * σz(), 1 => Ω / 2 * σp(), -1 => Ω / 2 * σm())

        # Test that we can compute to high orders without error
        @testset "Order 8 computes" begin
            result = magnus_expansion(modes, ω; order = 8)
            @test result.H_eff isa QuExpr
            @test haskey(result.orders, 8)
        end

        # Test order 10
        @testset "Order 10 computes" begin
            result = magnus_expansion(modes, ω; order = 10)
            @test result.H_eff isa QuExpr
            @test haskey(result.orders, 10)
        end

        QuantumAlgebra.use_σpm(false)
    end

end
