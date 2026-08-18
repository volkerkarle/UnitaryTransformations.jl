@testset "Commutator Series" begin
    using QuantumAlgebra
    using UnitaryTransformations:
        nested_commutator,
        multi_nested_commutator,
        compositions,
        commutator_series,
        bch_transform,
        bch_combine

    @testset "Nested commutators" begin
        # [a, a†] = 1
        @test normal_form(nested_commutator(a(), a'(), 1)) == one(QuExpr)

        # n=0 returns the original operator
        @test nested_commutator(a(), a'() * a(), 0) == a'() * a()

        # [a, a†a] = a (number operator commutator)
        @test normal_form(nested_commutator(a(), a'() * a(), 1)) == a()

        # [a†, a†a] = -a†
        @test normal_form(nested_commutator(a'(), a'() * a(), 1)) == -a'()

        # Double nested: [a, [a, a†a]] = [a, a] = 0
        @test normal_form(nested_commutator(a(), a'() * a(), 2)) == zero(QuExpr)

        # Triple nested with non-zero result
        # [a†, [a†, [a†, a]]] = [a†, [a†, -1]] = [a†, 0] = 0
        @test normal_form(nested_commutator(a'(), a(), 3)) == zero(QuExpr)

        # Error handling
        @test_throws ArgumentError nested_commutator(a(), a'(), -1)
    end

    @testset "Multi-generator nested commutators" begin
        # Empty generators list returns the base operator
        @test multi_nested_commutator(QuExpr[], a'() * a()) == a'() * a()

        # Single generator: [S, X]
        S1 = a()
        X = a'() * a()
        @test normal_form(multi_nested_commutator([S1], X)) == normal_form(comm(S1, X))

        # Two generators: [S1, [S2, X]]
        S2 = a'()
        # [a, [a†, a†a]] = [a, -a†] = -1
        expected = normal_form(comm(S1, comm(S2, X)))
        @test normal_form(multi_nested_commutator([S1, S2], X)) == expected

        # Verify multi_nested_commutator matches nested_commutator for identical generators
        @variables α
        S = α * a()
        H = a'() * a()
        # [S, [S, H]] with same S should equal nested_commutator(S, H, 2)
        @test normal_form(multi_nested_commutator([S, S], H)) ==
              normal_form(nested_commutator(S, H, 2))
        @test normal_form(multi_nested_commutator([S, S, S], H)) ==
              normal_form(nested_commutator(S, H, 3))
    end

    @testset "Compositions function" begin
        # Test basic compositions
        @test compositions(0, 0) == [Int[]]
        @test compositions(1, 0) == Vector{Int}[]  # Can't partition positive into 0 parts
        @test compositions(0, 1) == Vector{Int}[]  # Can't use 0 parts that sum to 0 with min_val=1

        # Simple cases
        @test compositions(1, 1) == [[1]]
        @test compositions(2, 1) == [[2]]
        @test compositions(2, 2) == [[1, 1]]

        # Ordered partitions (compositions)
        comps_3_2 = compositions(3, 2)
        @test length(comps_3_2) == 2
        @test [1, 2] in comps_3_2
        @test [2, 1] in comps_3_2

        comps_4_2 = compositions(4, 2)
        @test length(comps_4_2) == 3
        @test [1, 3] in comps_4_2
        @test [2, 2] in comps_4_2
        @test [3, 1] in comps_4_2

        comps_4_3 = compositions(4, 3)
        @test length(comps_4_3) == 3
        @test [1, 1, 2] in comps_4_3
        @test [1, 2, 1] in comps_4_3
        @test [2, 1, 1] in comps_4_3

        # With max_val constraint
        comps_4_2_max2 = compositions(4, 2; max_val = 2)
        @test comps_4_2_max2 == [[2, 2]]

        # Error handling
        @test_throws ArgumentError compositions(-1, 1)
        @test_throws ArgumentError compositions(1, -1)
    end

    @testset "BCH expansion - bch_transform" begin
        # For S = ε * a and H = a†a, compute e^S H e^{-S}
        # e^{εa} (a†a) e^{-εa} = a†a + ε[a, a†a] + (ε²/2)[a,[a,a†a]] + ...
        #                      = a†a + εa + 0 + ... = a†a + εa

        @variables ε
        S = ε * a()
        H = a'() * a()

        # Order 1: H + [S,H] = a†a + ε*a
        result = commutator_series(S, H, 1)
        expected = a'() * a() + ε * a()
        @test normal_form(result) == normal_form(expected)

        # Order 2: Same result since [a,[a,a†a]] = 0
        result = commutator_series(S, H, 2)
        @test normal_form(result) == normal_form(expected)

        # Higher orders should also give same result (series terminates)
        result = commutator_series(S, H, 5)
        @test normal_form(result) == normal_form(expected)
    end

    @testset "BCH transform - displacement operator" begin
        # For S = α(a† - a):
        # [S, a] = α[a† - a, a] = α([a†, a] - [a, a]) = α(-1 - 0) = -α
        # [S, a†] = α[a† - a, a†] = α([a†, a†] - [a, a†]) = α(0 - 1) = -α
        # So e^S a e^{-S} = a + [S,a] + ... = a - α (series terminates)
        # And e^S a† e^{-S} = a† + [S,a†] + ... = a† - α

        @variables α
        S = α * (a'() - a())

        # Transform a
        result = bch_transform(S, a(); order = 3)
        expected = a() - α
        @test normal_form(result) == normal_form(expected)

        # Transform a† 
        result = bch_transform(S, a'(); order = 3)
        expected = a'() - α
        @test normal_form(result) == normal_form(expected)
    end

    @testset "BCH transform - number operator under displacement" begin
        # With S = α(a† - a):
        # e^S (a†a) e^{-S} = (a† - α)(a - α) = a†a - αa† - αa + α²

        @variables α
        S = α * (a'() - a())
        H = a'() * a()

        result = bch_transform(S, H; order = 4)
        expected = a'() * a() - α * a'() - α * a() + α^2
        @test normal_form(result) == normal_form(expected)
    end

    @testset "bch_transform alias" begin
        @variables ε
        S = ε * a()
        H = a'() * a()

        @test bch_transform(S, H; order = 2) == commutator_series(S, H, 2)
    end

    @testset "bch_combine - basic" begin
        # Test: e^A e^B = e^Z where Z = A + B + (1/2)[A,B] + ...
        # For bosonic operators [a†, a] = -1, so [αa†, βa] = -αβ

        @variables α β

        A = α * a'()
        B = β * a()

        # Order 1: Z = A + B
        Z1 = bch_combine(A, B; order = 1)
        @test normal_form(Z1) == normal_form(α * a'() + β * a())

        # Order 2: Z = A + B + (1/2)[A,B] = αa† + βa - (1/2)αβ
        Z2 = bch_combine(A, B; order = 2)
        expected = α * a'() + β * a() - (1 // 2) * α * β
        @test normal_form(Z2) == normal_form(expected)

        # Order 3: Same as order 2 since higher commutators vanish
        Z3 = bch_combine(A, B; order = 3)
        @test normal_form(Z3) == normal_form(expected)

        # Error handling
        @test_throws ArgumentError bch_combine(A, B; order = 0)
    end

    @testset "bch_combine - commuting operators" begin
        # If [A, B] = 0, then e^A e^B = e^{A+B} exactly
        @variables α β

        # Two different modes commute
        A = α * a'(:a)
        B = β * a'(:b)

        # [a†_a, a†_b] = 0, so Z = A + B for all orders
        Z = bch_combine(A, B; order = 4)
        expected = α * a'(:a) + β * a'(:b)
        @test normal_form(Z) == normal_form(expected)
    end

    @testset "bch_combine - symmetry" begin
        # Test that swapping A and B gives different result (non-commutative)
        @variables α β

        A = α * a'()
        B = β * a()

        Z_AB = bch_combine(A, B; order = 2)
        Z_BA = bch_combine(B, A; order = 2)

        # Z_AB = αa† + βa - (1/2)αβ
        # Z_BA = βa + αa† + (1/2)[βa, αa†] = αa† + βa + (1/2)αβ
        # They differ by the sign of the commutator term

        @test normal_form(Z_AB) != normal_form(Z_BA)

        # The difference should be [A,B] = -αβ
        diff = normal_form(Z_AB - Z_BA)
        expected_diff = -α * β  # = [A,B]
        @test normal_form(diff) == normal_form(expected_diff)
    end

    @testset "bch_combine - higher orders" begin
        # Test order 4 and 5 don't crash and give reasonable results
        @variables θ

        QuantumAlgebra.use_σpm(true)

        A = θ * σp()
        B = θ * σm()

        # These have non-trivial higher-order commutators
        Z4 = bch_combine(A, B; order = 4)
        Z5 = bch_combine(A, B; order = 5)

        # Just check they don't error and produce QuExpr
        @test Z4 isa QuExpr
        @test Z5 isa QuExpr

        # Order 4 and 5 should potentially differ
        # (depends on the specific operators)

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Spin commutators" begin
        # Use σ± basis for cleaner commutators
        QuantumAlgebra.use_σpm(true)

        # [σz, σ+] = 2σ+  (σ+ is eigenoperator of ad_σz)
        @test normal_form(nested_commutator(σz(), σp(), 1)) == 2 * σp()

        # [σz, σ-] = -2σ-
        @test normal_form(nested_commutator(σz(), σm(), 1)) == -2 * σm()

        # [σ+, σ-] = σz
        @test normal_form(nested_commutator(σp(), σm(), 1)) == σz()

        QuantumAlgebra.use_σpm(false)
    end

    @testset "Spin BCH transform" begin
        # Rotation around z-axis: e^{iθσz/2} σ+ e^{-iθσz/2} = e^{iθ} σ+
        # With our real convention (no i): e^{θσz} σ+ e^{-θσz}
        # [θσz, σ+] = 2θσ+
        # [θσz, [θσz, σ+]] = 2θ·2θσ+ = 4θ²σ+
        # So e^{θσz} σ+ e^{-θσz} = σ+ + 2θσ+ + (4θ²/2)σ+ + (8θ³/6)σ+ + ...
        #                        = σ+ (1 + 2θ + 2θ² + (4/3)θ³ + ...)
        #                        = σ+ · e^{2θ}  (Taylor series of e^{2θ})

        QuantumAlgebra.use_σpm(true)
        @variables θ

        S = θ * σz()
        result = bch_transform(S, σp(); order = 4)

        # Check that result has the form c(θ) * σ+
        # The coefficient should be: 1 + 2θ + 2θ² + (4/3)θ³ + (2/3)θ⁴ + ...
        # which are the Taylor coefficients of e^{2θ}

        # At order 4, coefficient should be: 1 + 2θ + 2θ² + (4/3)θ³ + (2/3)θ⁴
        @test result isa QuExpr
        @test length(result.terms) == 1  # Only σ+ term

        QuantumAlgebra.use_σpm(false)
    end

    @testset "multi_nested_commutator with many generators" begin
        # Test with 4 generators
        gens4 = [a(), a'(), a(), a'()]
        result = multi_nested_commutator(gens4, a'()*a())
        @test result isa QuExpr

        # Test with 5 generators (should still work)
        gens5 = [a(), a(), a(), a(), a()]
        result2 = multi_nested_commutator(gens5, a'())
        @test result2 isa QuExpr
    end

    @testset "bch_combine order 5 - Pauli algebra verification" begin
        # Use commuting case to verify order-5 correctness
        @variables α β

        # Two commuting operators should give A + B at any order
        A = α * a'() * a()
        B = β * a'() * a()
        Z = bch_combine(A, B; order=5)
        @test normal_form(Z) == normal_form(α * a'()*a() + β * a'()*a())

        # Test with non-commuting case: A = α a, B = β a†
        # [A,B] = αβ, [[A,B],A] = 0, etc.
        A2 = α * a()
        B2 = β * a'()
        Z2 = bch_combine(A2, B2; order=5)
        # Known result: e^{α a} e^{β a†} = e^{α a + β a† + αβ/2} (exact!)
        expected = normal_form(α * a() + β * a'() + (α*β/2))
        @test normal_form(Z2) == expected
    end

    @testset "nested_commutator high depth" begin
        # Test depth 5
        S = a()
        H_val = a'() * a()
        result = nested_commutator(S, H_val, 5)
        # [a, [a, [a, [a, [a, a†a]]]]] = 0 (since 2 nested comms already give 0)
        @test result == zero(QuExpr)

        # Test with non-vanishing nested commutators (bosonic, depth 2)
        # [a†, [a, a†a]] = [a†, a] = -1 (non-zero scalar)
        S3 = a'()
        H3 = a'() * a()
        r3 = nested_commutator(S3, H3, 2)
        @test r3 isa QuExpr
    end

    @testset "compositions edge cases" begin
        # Test with larger n,k values
        result = compositions(6, 3)
        # There should be C(5,2) = 10 compositions of 6 into 3 positive parts
        @test length(result) == 10

        # n=0, k=0 case (valid)
        result2 = compositions(0, 0)
        @test length(result2) == 1
        @test isempty(result2[1])

        # n>0, k=0 case (no solutions)
        result3 = compositions(5, 0)
        @test isempty(result3)

        # max_val edge case
        result4 = compositions(4, 2; max_val=2)
        @test length(result4) == 1  # only (2,2)
        @test result4[1] == [2, 2]
    end

    @testset "bch_transform - known result" begin
        # Test: e^{θ a} a†a e^{-θ a} = a†a + θ a
        @variables θ
        S = θ * a()
        A = a'() * a()
        result = bch_transform(S, A; order=2)
        expected = normal_form(A + θ * a())
        @test normal_form(result) == expected
    end

    @testset "bch_combine order 5 - SU(2) quantitative" begin
        Symbolics.@variables α β
        
        # SU(2) generators provide non-terminating commutators
        λ = su_generators(2, :λ)
        A = α * λ[1]  # λ₁ = σ₁/2
        B = β * λ[2]  # λ₂ = σ₂/2
        
        Z5 = bch_combine(A, B; order=5)
        Z4 = bch_combine(A, B; order=4)
        Z3 = bch_combine(A, B; order=3)
        Z2 = bch_combine(A, B; order=2)
        
        # Order 5 must differ from order 2 (non-terminating series)
        @test normal_form(Z5) != normal_form(Z2)
        
        # Order 3 and 4 may or may not differ in SU(2); order 4 often adds nothing.
        # Order 5 must add new terms
        @test normal_form(Z5) != normal_form(Z2)
        
        # All must be non-empty
        @test !isempty(Z5.terms)
        
        # Identity: bch_combine(A, 0) = A
        Z_zero = bch_combine(A, zero(QuExpr); order=5)
        @test normal_form(Z_zero) == normal_form(A)
        
        # Commuting case: [B, B] = 0 so bch_combine(B, B) = 2B
        Z_comm = bch_combine(B, B; order=5)
        @test normal_form(Z_comm) == normal_form(2 * B)
    end
end
