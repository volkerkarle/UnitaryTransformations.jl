@testset "Commutator Series" begin
    using QuantumAlgebra
    using UnitaryTransformations: nested_commutator, commutator_series, bch_transform

    @testset "Nested commutators" begin
        # [a, a†] = 1
        @test normal_form(nested_commutator(a(), a'(), 1)) == one(QuExpr)

        # n=0 returns the original operator
        @test nested_commutator(a(), a'()*a(), 0) == a'()*a()

        # [a, a†a] = a (number operator commutator)
        @test normal_form(nested_commutator(a(), a'()*a(), 1)) == a()

        # [a†, a†a] = -a†
        @test normal_form(nested_commutator(a'(), a'()*a(), 1)) == -a'()

        # Double nested: [a, [a, a†a]] = [a, a] = 0
        @test normal_form(nested_commutator(a(), a'()*a(), 2)) == zero(QuExpr)
    end

    @testset "BCH expansion" begin
        # For S = ε * a and H = a†a, compute e^S H e^{-S}
        # e^{εa} (a†a) e^{-εa} = a†a + ε[a, a†a] + (ε²/2)[a,[a,a†a]] + ...
        #                      = a†a + εa + 0 + ... = a†a + εa

        @variables ε
        S = ε * a()
        H = a'() * a()

        # Order 1: H + [S,H] = a†a + ε*a
        result = commutator_series(S, H, 1)
        expected = a'()*a() + ε*a()
        @test normal_form(result) == normal_form(expected)

        # Order 2: Same result since [a,[a,a†a]] = 0
        result = commutator_series(S, H, 2)
        @test normal_form(result) == normal_form(expected)
    end

    @testset "bch_transform alias" begin
        @variables ε
        S = ε * a()
        H = a'() * a()

        @test bch_transform(S, H; order = 2) == commutator_series(S, H, 2)
    end

    @testset "Spin commutators" begin
        # Use σ± basis for cleaner commutators
        QuantumAlgebra.use_σpm(true)

        # [σz, σ+] = 2σ+  (σ+ is eigenoperator of ad_σz)
        @test normal_form(nested_commutator(σz(), σp(), 1)) == 2*σp()

        # [σz, σ-] = -2σ-
        @test normal_form(nested_commutator(σz(), σm(), 1)) == -2*σm()

        # [σ+, σ-] = σz
        @test normal_form(nested_commutator(σp(), σm(), 1)) == σz()

        QuantumAlgebra.use_σpm(false)
    end
end
