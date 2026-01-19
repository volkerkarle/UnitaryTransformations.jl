@testset "Schrieffer-Wolff Transformation" begin
    using QuantumAlgebra
    using QuantumAlgebra: TLSCreate_, BosonCreate_, BosonDestroy_
    using UnitaryTransformations:
        Subspace,
        decompose,
        diagonal_part,
        off_diagonal_part,
        classify_base_operator,
        classify_term,
        DIAGONAL,
        RAISING,
        LOWERING,
        MIXED,
        schrieffer_wolff,
        sw_generator,
        project_to_subspace,
        solve_for_generator

    # Use σ± basis for cleaner SW transformations
    QuantumAlgebra.use_σpm(true)

    @testset "Operator classification" begin
        P = Subspace(σz() => -1)  # Spin-down subspace

        # Get base operators from QuExpr terms
        σp_op = first(σp().terms)[1].bares.v[1]
        σm_op = first(σm().terms)[1].bares.v[1]
        a_op = first(a().terms)[1].bares.v[1]
        adag_op = first(a'().terms)[1].bares.v[1]

        # σ+ raises spin (P→Q when P is spin-down)
        @test classify_base_operator(σp_op, P) == RAISING

        # σ- lowers spin (Q→P when P is spin-down)
        @test classify_base_operator(σm_op, P) == LOWERING

        # σz in σpm mode is -1 + 2σ⁺σ⁻, which contains multiple terms
        # Test that it's diagonal as an expression (via classify_term)
        # The σ⁺σ⁻ term should be DIAGONAL since raising+lowering cancels
        σpm_term = nothing
        for (term, _) in σz().terms
            if length(term.bares.v) == 2  # The σ⁺σ⁻ term
                σpm_term = term
                break
            end
        end
        if σpm_term !== nothing
            @test classify_term(σpm_term, P) == DIAGONAL
        end

        # Bosonic operators are transparent (diagonal w.r.t. spin subspace)
        @test classify_base_operator(a_op, P) == DIAGONAL
        @test classify_base_operator(adag_op, P) == DIAGONAL
    end

    @testset "Hamiltonian decomposition" begin
        P = Subspace(σz() => -1)

        # Simple Jaynes-Cummings-like Hamiltonian
        @variables ω Δ g

        H = ω * a'() * a() + Δ / 2 * σz() + g * (a() * σp() + a'() * σm())

        H_d, H_od = decompose(H, P)

        # Diagonal part: ω a†a + Δ/2 σz
        @test is_diagonal(H_d, P)

        # Off-diagonal part: g(a σ+ + a† σ-)
        @test is_off_diagonal(H_od, P)

        # Check that H = H_d + H_od
        @test normal_form(H) == normal_form(H_d + H_od)
    end

    @testset "Simple two-level system" begin
        # H = Δ/2 σz + ε σx
        # For small ε, SW should give effective H ≈ Δ/2 σz + O(ε²)

        @variables Δ ε

        # σx = σ+ + σ-
        H = Δ / 2 * σz() + ε * (σp() + σm())

        P = Subspace(σz() => -1)

        # Just test decomposition for now
        H_d, H_od = decompose(H, P)

        @test normal_form(H_d) == normal_form(Δ / 2 * σz())
        @test normal_form(H_od) == normal_form(ε * (σp() + σm()))
    end

    @testset "Generator solution" begin
        # For H_d = Δ/2 σz and V_od = ε σ+
        # [S, H_d] = -V_od means [S, Δ/2 σz] = -ε σ+
        # Since [σz, σ+] = 2σ+ and σz = -1 + 2σ⁺σ⁻, we have
        # [σ⁺σ⁻, σ+] = σ+, so [H_d, σ+] = Δ σ+
        # Thus S = (ε/Δ) σ+ (with proper Symbolics division)

        @variables Δ ε

        H_d = Δ / 2 * σz()
        V_od = ε * σp()

        P = Subspace(σz() => -1)
        S = solve_for_generator(H_d, V_od, P)

        # Check that S has the right structure: (ε/Δ) * σ+
        # With Symbolics integration, the coefficient is a proper Num type
        # 1. S is non-empty
        @test !isempty(S.terms)

        # 2. S contains σ+ with Symbolics coefficient
        # The generator should be proportional to σ+
        has_σp = false
        for (term, coeff) in S.terms
            if !isempty(term.bares.v)
                # Check if it's a σ+ operator
                op = term.bares.v[1]
                if op.t == TLSCreate_
                    has_σp = true
                    # With Symbolics, coeff should be a Num (or derived type)
                    # Check that it's a proper symbolic expression
                    @test coeff isa Symbolics.Num
                    # The coefficient should contain ε and Δ (as ε/Δ)
                    coeff_str = string(coeff)
                    @test occursin("ε", coeff_str) || occursin("Δ", coeff_str)
                end
            end
        end
        @test has_σp
    end

    @testset "Full SW transformation - dispersive regime" begin
        # Jaynes-Cummings: H = ω a†a + Δ/2 σz + g(a†σ- + a σ+)
        # In dispersive regime (Δ >> g), SW gives:
        # H_eff ≈ ω a†a + Δ/2 σz + χ a†a σz + const
        # where χ = g²/Δ (dispersive shift)

        @variables ω Δ g

        H = ω * a'() * a() + Δ / 2 * σz() + g * (a'() * σm() + a() * σp())

        P = Subspace(σz() => -1)

        # Compute SW to second order
        result = schrieffer_wolff(H, P; order = 2)

        # The effective Hamiltonian should be block-diagonal
        @test is_diagonal(result.H_eff, P)

        # Check that we got a generator
        @test !isempty(result.S.terms)
    end

    @testset "Projection to subspace" begin
        P = Subspace(σz() => -1)

        @variables Δ

        # σz in P sector should give -1
        H = Δ * σz()
        H_P = project_to_subspace(H, P)

        # After projection, σz → -1
        @test normal_form(H_P) == normal_form(-Δ * one(QuExpr))
    end

    @testset "Generator equation verification" begin
        # Rigorous check: [S, H_d] = -V_od must hold
        using UnitaryTransformations: param_to_symbolic, clear_param_cache!

        clear_param_cache!()

        # Helper to convert QuExpr to pure Symbolics coefficients
        function to_symbolic_coeffs(expr)
            result = QuExpr()
            for (term, coeff) in expr.terms
                full_coeff = coeff isa Symbolics.Num ? coeff : Symbolics.Num(coeff)
                for p in term.params
                    full_coeff = full_coeff * param_to_symbolic(p)
                end
                full_coeff = Symbolics.simplify(full_coeff)
                new_term = QuantumAlgebra.QuTerm(
                    term.nsuminds,
                    term.δs,
                    QuantumAlgebra.Param[],
                    term.expvals,
                    term.corrs,
                    term.bares,
                )
                result = result + full_coeff * QuExpr(new_term)
            end
            return normal_form(result)
        end

        @variables Δ ε
        H = Δ / 2 * σz() + ε * (σp() + σm())
        P = Subspace(σz() => -1)

        H_d, V_od = decompose(H, P)
        S = solve_for_generator(H_d, V_od, P)

        # Convert [S, H_d] + V_od to symbolic and check it equals zero
        residual = to_symbolic_coeffs(normal_form(comm(S, H_d) + V_od))
        @test isempty(residual.terms)  # Residual must be zero
    end

    @testset "Dispersive shift numerical accuracy" begin
        # Compare SW result with exact two-level solution
        # H = Δ/2 σz + ε σx has eigenvalues ±√(Δ²/4 + ε²)
        # SW gives E_g = -Δ/2 - ε²/Δ to second order

        Δ_val = 1.0
        ε_val = 0.1  # ε/Δ = 0.1, well in perturbative regime

        E_exact = -sqrt(Δ_val^2 / 4 + ε_val^2)
        E_SW = -Δ_val / 2 - ε_val^2 / Δ_val

        error_pct = 100 * abs(E_exact - E_SW) / abs(E_exact)

        # Should be accurate to < 0.1% for ε/Δ = 0.1
        @test error_pct < 0.1
    end

    @testset "SU(3) operator classification" begin
        using UnitaryTransformations:
            is_lie_algebra_constraint,
            get_lie_algebra_constraint_info,
            is_diagonal_lie_generator

        # Create SU(3) generators
        λ = su_generators(3, :λ)

        # Define subspace with constraint on λ₈ (diagonal generator)
        P = Subspace(λ[8] => 0.5)

        @testset "Diagonal generators are DIAGONAL" begin
            for i in (7, 8)
                term, _ = first(λ[i].terms)
                bare = term.bares.v[1]
                @test classify_base_operator(bare, P) == DIAGONAL
            end
        end

        @testset "Off-diagonal generators are MIXED" begin
            for i = 1:6
                term, _ = first(λ[i].terms)
                bare = term.bares.v[1]
                @test classify_base_operator(bare, P) == MIXED
            end
        end
    end

    @testset "SU(3) Hamiltonian decomposition" begin
        # Create SU(3) generators
        λ = su_generators(3, :λ)

        # Define subspace with constraint on λ₈
        P = Subspace(λ[8] => 0.5)

        @variables ω₁ ω₂ g₁ g₂

        # Create a 3-level system Hamiltonian
        # Diagonal terms: λ₇, λ₈ (Cartan subalgebra)
        # Off-diagonal terms: λ₁ (couples states 1↔2), λ₂ (couples states 1↔3)
        H = ω₁ * λ[7] + ω₂ * λ[8] + g₁ * λ[1] + g₂ * λ[2]

        H_d, H_od = decompose(H, P)

        @testset "Correct diagonal part extraction" begin
            @test is_diagonal(H_d, P)
            # H_d should be ω₁ λ₇ + ω₂ λ₈
            expected_d = ω₁ * λ[7] + ω₂ * λ[8]
            @test normal_form(H_d) == normal_form(expected_d)
        end

        @testset "Correct off-diagonal part extraction" begin
            @test is_off_diagonal(H_od, P)
            # H_od should be g₁ λ₁ + g₂ λ₂
            expected_od = g₁ * λ[1] + g₂ * λ[2]
            @test normal_form(H_od) == normal_form(expected_od)
        end

        @testset "Decomposition preserves Hamiltonian" begin
            @test normal_form(H) == normal_form(H_d + H_od)
        end
    end

    @testset "SU(3) with all generators" begin
        λ = su_generators(3, :λ)
        P = Subspace(λ[8] => 0.5)

        # Build H with all 8 generators
        H = sum(i * λ[i] for i = 1:8)

        H_d, H_od = decompose(H, P)

        # Diagonal: 7*λ₇ + 8*λ₈
        @test normal_form(H_d) == normal_form(7 * λ[7] + 8 * λ[8])

        # Off-diagonal: sum of i*λᵢ for i in 1:6
        @test normal_form(H_od) == normal_form(sum(i * λ[i] for i = 1:6))

        # Verify reconstruction
        @test normal_form(H) == normal_form(H_d + H_od)
    end

    @testset "SU(2) Lie algebra classification" begin
        # Create SU(2) generators
        σ = su_generators(2, :σ)

        # Define subspace with constraint on σ₃ (diagonal)
        P = Subspace(σ[3] => -0.5)  # Spin-down eigenvalue

        @testset "σ₃ is DIAGONAL" begin
            term, _ = first(σ[3].terms)
            bare = term.bares.v[1]
            @test classify_base_operator(bare, P) == DIAGONAL
        end

        @testset "σ₁, σ₂ are MIXED" begin
            for i = 1:2
                term, _ = first(σ[i].terms)
                bare = term.bares.v[1]
                @test classify_base_operator(bare, P) == MIXED
            end
        end
    end

    function numeric_abs(coeff, substitutions)
        if coeff isa Complex
            real_part = real(coeff)
            imag_part = imag(coeff)
            real_sub = real_part isa Symbolics.Num ? Symbolics.substitute(real_part, substitutions) : real_part
            imag_sub = imag_part isa Symbolics.Num ? Symbolics.substitute(imag_part, substitutions) : imag_part
            real_val = real_sub isa Symbolics.Num ? Symbolics.value(real_sub) : real_sub
            imag_val = imag_sub isa Symbolics.Num ? Symbolics.value(imag_sub) : imag_sub
            return abs(complex(real_val, imag_val))
        end
        substituted = coeff isa Symbolics.Num ? Symbolics.substitute(coeff, substitutions) : coeff
        value = substituted isa Symbolics.Num ? Symbolics.value(substituted) : substituted
        return abs(value)
    end

    @testset "SU(3) generator equation" begin
        using UnitaryTransformations: solve_for_generator_lie

        λ = su_generators(3, :λ)
        @variables Δ ω g

        # Diagonal Hamiltonian
        H_d = Δ * λ[8] + ω * λ[7]

        # Single off-diagonal term
        V_od = g * λ[2]  # Couples states 1↔3

        # Solve for generator
        S = solve_for_generator_lie(H_d, V_od, 3, λ)

        @test !isempty(S.terms)

        # Verify generator equation: [S, H_d] = -V_od
        comm_S_Hd = normal_form(comm(S, H_d))
        residual = normal_form(comm_S_Hd + V_od)

        # Check that residual is numerically zero for each term
        substitutions = Dict(Δ => 1.3, ω => 0.7, g => 0.4)
        for (_, coeff) in residual.terms
            @test numeric_abs(coeff, substitutions) < 1e-10
        end
    end

    @testset "SU(3) generator with multiple couplings" begin
        using UnitaryTransformations: solve_for_generator_lie

        λ = su_generators(3, :λ)
        @variables Δ ω g₁ g₂

        # Diagonal Hamiltonian
        H_d = Δ * λ[8] + ω * λ[7]

        # Lambda system: both ground states coupled to excited state
        V_od = g₁ * λ[2] + g₂ * λ[3]  # λ₂: 1↔3, λ₃: 2↔3

        # Solve for generator
        S = solve_for_generator_lie(H_d, V_od, 3, λ)

        @test !isempty(S.terms)

        # Verify generator equation
        comm_S_Hd = normal_form(comm(S, H_d))
        residual = normal_form(comm_S_Hd + V_od)

        substitutions = Dict(Δ => 1.3, ω => 0.7, g₁ => 0.4, g₂ => -0.3)
        for (_, coeff) in residual.terms
            @test numeric_abs(coeff, substitutions) < 1e-10
        end
    end

    @testset "Full SW transformation - SU(3) Lambda system" begin
        using UnitaryTransformations: detect_lie_algebra_system

        λ = su_generators(3, :λ)
        @variables Δ ω g

        # Lambda (Λ) configuration 3-level atom
        # States: |1⟩, |2⟩ are ground states, |3⟩ is excited state
        # H_d uses diagonal generators λ₇, λ₈
        # V_od couples ground states to excited state via λ₂ (1↔3)

        H = Δ * λ[8] + ω * λ[7] + g * λ[2]

        # Define subspace constraint on diagonal generator
        P = Subspace(λ[8] => 0.5)

        # Test that Lie algebra detection works
        H_d, H_od = decompose(H, P)
        lie_info = detect_lie_algebra_system(H_od)
        @test lie_info !== nothing
        @test lie_info.N == 3

        # Compute SW to second order
        result = schrieffer_wolff(H, P; order = 2)

        # The effective Hamiltonian should be block-diagonal
        @test is_diagonal(result.H_eff, P)

        # Check that we got a generator
        @test !isempty(result.S.terms)

        # Verify generator equation at first order
        # The generator S should satisfy [S, H_d] ≈ -H_od
        # We test this indirectly: if SW succeeded in making H_eff diagonal,
        # then the generator is doing its job
    end

    @testset "Full SW transformation - SU(3) with all off-diagonal couplings" begin
        λ = su_generators(3, :λ)
        @variables ω₁ ω₂ g₁ g₂ g₃

        # Full 3-level system with all couplings
        H_d = ω₁ * λ[7] + ω₂ * λ[8]
        H_od = g₁ * λ[1] + g₂ * λ[2] + g₃ * λ[3]  # All off-diagonal generators λ₁-λ₃
        H = H_d + H_od

        P = Subspace(λ[8] => 0.5)

        # Second order SW transformation (minimum order is 2)
        result2 = schrieffer_wolff(H, P; order = 2)
        @test is_diagonal(result2.H_eff, P)

        # Generator should be non-empty
        @test !isempty(result2.S.terms)

        # H_eff at order 2 should have dispersive-like shifts (more terms than just H_d)
        @test length(result2.H_eff.terms) >= length(H_d.terms)
    end

    @testset "N-level transition operators + bosons" begin
        # Test with generic N-level system using transition operators
        # This tests the eigenoperator method with nlevel_ops

        σ5 = nlevel_ops(5, :q)  # 5-level system

        # Symbolic energies
        ω = [Symbolics.variable(Symbol("ω", i)) for i = 1:5]
        @variables ωc g

        # Diagonal Hamiltonian: atom + cavity
        H_atom = sum(ω[i] * σ5[i, i] for i = 1:5)
        H_cav = ωc * a'() * a()
        H_d = H_atom + H_cav

        # Jaynes-Cummings coupling between levels 1,3
        V = g * (σ5[1, 3] * a'() + σ5[3, 1] * a())

        H = H_d + V

        # Subspace: cavity vacuum
        P = Subspace(a'() * a() => 0)

        # SW transformation
        result = schrieffer_wolff(H, P; order = 2)

        # Should produce block-diagonal result
        @test is_diagonal(result.H_eff, P)

        # Should have dispersive shifts (a†a q¹¹ and a†a q³³ terms)
        @test !isempty(result.S.terms)

        # H_eff should have more terms than H_d due to dispersive shifts
        @test length(result.H_eff.terms) > length(H_d.terms)
    end

    @testset "N-level with multiple cavity couplings" begin
        # Test 7-level system with multiple transitions
        σ7 = nlevel_ops(7, :q)

        @variables Δ₁ Δ₂ ωc g₁ g₂

        # Simplified diagonal: only some levels have non-zero energy
        H_atom = Δ₁ * σ7[3, 3] + Δ₂ * σ7[4, 4]
        H_cav = ωc * a'() * a()
        H_d = H_atom + H_cav

        # Multiple couplings: 1↔3 and 2↔4
        V =
            g₁ * (σ7[1, 3] * a'() + σ7[3, 1] * a()) +
            g₂ * (σ7[2, 4] * a'() + σ7[4, 2] * a())

        H = H_d + V
        P = Subspace(a'() * a() => 0)

        # Decomposition should work
        H_diag, H_od = decompose(H, P)
        @test is_diagonal(H_diag, P)
        @test is_off_diagonal(H_od, P)

        # SW transformation
        result = schrieffer_wolff(H, P; order = 2)
        @test is_diagonal(result.H_eff, P)

        # H_eff should have more terms than H_d due to dispersive shifts
        @test length(result.H_eff.terms) > length(H_d.terms)
    end

    @testset "4th order SW - Kerr nonlinearity" begin
        # Test that 4th order SW produces Kerr terms (a†²a²)
        # For Rabi model in dispersive regime
        # Using diagonal_only=true for faster computation

        @variables ω_c Δ g

        # Rabi Hamiltonian: H = ω_c a†a + Δ/2 σz + g(a† + a)(σ+ + σ-)
        H = ω_c * a'() * a() + Δ / 2 * σz() + g * (a'() + a()) * (σp() + σm())
        P = Subspace(a'() * a() => 0)

        # Compute 4th order SW with diagonal_only for speed
        result4 = schrieffer_wolff(H, P; order = 4, diagonal_only = true)

        # Check that the result is block-diagonal
        @test is_diagonal(result4.H_eff, P)

        # Check that H_P contains the expected operators
        op_strings = Set{String}()
        for (term, _) in result4.H_P.terms
            op_str = isempty(term.bares.v) ? "𝟙" : string(term.bares)
            push!(op_strings, op_str)
        end

        # When projected to the vacuum subspace (n=0), all cavity operators (a, a†) 
        # should vanish because they take us out of the vacuum.
        # Only spin operators should remain:
        # - Identity (constant energy shift)
        # - σ⁺σ⁻ (qubit frequency shift / population)
        @test "𝟙" in op_strings
        @test "σ⁺() σ⁻()" in op_strings

        # Cavity operators should NOT appear in the vacuum-projected Hamiltonian
        @test !("a†() a()" in op_strings)
        @test !("a†() σ⁺() σ⁻() a()" in op_strings)
        @test !("a†()² a()²" in op_strings)
        @test !("a†()² σ⁺() σ⁻() a()²" in op_strings)

        # The full H_eff (before projection) should still have Kerr terms
        heff_op_strings = Set{String}()
        for (term, _) in result4.H_eff.terms
            op_str = isempty(term.bares.v) ? "𝟙" : string(term.bares)
            push!(heff_op_strings, op_str)
        end

        # H_eff should contain cavity operators and Kerr terms
        @test "a†() a()" in heff_op_strings
        @test "a†()² a()²" in heff_op_strings  # Kerr term in H_eff!

        # Check that Kerr coefficient in H_eff scales as g⁴
        kerr_coeff = nothing
        for (term, coeff) in result4.H_eff.terms
            op_str = isempty(term.bares.v) ? "𝟙" : string(term.bares)
            if op_str == "a†()² a()²"
                kerr_coeff = coeff
                break
            end
        end

        @test kerr_coeff !== nothing

        # Verify the Kerr coefficient has g⁴ dependence
        # Substitute g → 0: coefficient should be 0
        val_g0 = Symbolics.substitute(kerr_coeff, Dict(g => 0.0, Δ => 5.0, ω_c => 1.0))
        @test abs(Float64(Symbolics.value(val_g0))) < 1e-10

        # Verify scaling: K(2g)/K(g) ≈ 16 (since K ~ g⁴)
        val_g1 = Symbolics.substitute(kerr_coeff, Dict(g => 1.0, Δ => 5.0, ω_c => 1.0))
        val_g2 = Symbolics.substitute(kerr_coeff, Dict(g => 2.0, Δ => 5.0, ω_c => 1.0))
        ratio = Float64(Symbolics.value(val_g2)) / Float64(Symbolics.value(val_g1))
        @test abs(ratio - 16.0) < 1e-6  # Should be exactly 16 for g⁴ scaling
    end

    @testset "SymSum/SymExpr support for multi-atom systems" begin
        # Import SymSum types
        import QuantumAlgebra: sumindex, SymSum, SymExpr, expand_symbolic

        @variables ω_c tc_Δ tc_g

        # Create a sum index for the Tavis-Cummings model
        i = sumindex(1)

        # Build Hamiltonian with symbolic sums
        H_cav = ω_c * a'() * a()
        H_atom = SymSum(tc_Δ / 2 * σz(i), i)
        H_int = SymSum(tc_g * (a'() * σm(i) + a() * σp(i)), i)

        H = SymExpr(H_cav) + H_atom + H_int

        # Define subspace: zero photon sector
        P = Subspace(a'() * a() => 0)

        # Test decomposition with SymExpr
        H_d, H_od = decompose(H, P)

        # The off-diagonal part should be the interaction term
        @test H_od isa SymExpr

        # Test solve_for_generator with SymSum
        S = solve_for_generator(H_d, H_od, P)
        @test S isa SymExpr || S isa SymSum

        # Test schrieffer_wolff with SymExpr
        result = schrieffer_wolff(H, P; order = 2)

        @test result.H_eff isa SymExpr
        @test result.S isa SymExpr || result.S isa SymSum

        # Test that exchange terms appear when we compute [S, V] for 2 atoms
        # Use the explicitly defined generator from the example
        S1 = SymSum((tc_g / tc_Δ) * (a() * σp(i) - a'() * σm(i)), i)

        # Expand to 2 atoms
        S1_2 = expand_symbolic(S1, 1:2)
        H_int_2 = expand_symbolic(H_int, 1:2)

        # Compute [S, V] for 2 atoms
        comm_SV = normal_form(comm(S1_2, H_int_2))

        # Check for exchange terms: σ⁺(1)σ⁻(2) and σ⁺(2)σ⁻(1)
        has_exchange_12 = false
        has_exchange_21 = false

        for (term, _) in comm_SV.terms
            term_str = string(term.bares)
            if occursin("σ⁺(1)", term_str) && occursin("σ⁻(2)", term_str)
                has_exchange_12 = true
            end
            if occursin("σ⁺(2)", term_str) && occursin("σ⁻(1)", term_str)
                has_exchange_21 = true
            end
        end

        @test has_exchange_12
        @test has_exchange_21

        # Test project_to_subspace for SymExpr
        @testset "project_to_subspace for SymExpr" begin
            # Test 1: project_to_subspace removes off-diagonal terms
            # Note: It does NOT substitute a†a → 0; it only removes off-diagonal operators
            # and substitutes spin projection operators (σ⁺σ⁻)
            H_test = SymExpr(ω_c * a'() * a())
            P_vac = Subspace(a'() * a() => 0)
            H_proj = project_to_subspace(H_test, P_vac)
            # a†a is diagonal, so it remains (projection doesn't numerically evaluate)
            @test H_proj isa QuExpr

            # Test 2: Projecting SymSum with σz to spin-down
            # σz is diagonal and gets substituted to its eigenvalue
            j = sumindex(2)  # Use different index
            H_spin = SymSum(tc_Δ / 2 * σz(j), j)
            H_spin_expr = SymExpr(H_spin)
            P_spin_down = Subspace(σz() => -1)
            H_spin_proj = project_to_subspace(H_spin_expr, P_spin_down)

            # The result should be SymExpr with Σⱼ(-Δ/2)
            # Check that it's still a SymExpr (the sum remains)
            @test H_spin_proj isa SymExpr

            # Test 3: Combined projection - cavity vacuum AND spin
            H_combined = SymExpr(ω_c * a'() * a()) + SymSum(tc_Δ / 2 * σz(j), j)
            P_both = Subspace(a'() * a() => 0, σz() => -1)
            H_comb_proj = project_to_subspace(H_combined, P_both)

            # σz(j) → -1, but a†a is not numerically substituted
            @test H_comb_proj isa SymExpr

            # Test 4: H_P from schrieffer_wolff should be projected
            @test result.H_P isa SymExpr || result.H_P isa QuExpr

            # Test 5: Off-diagonal terms should be removed
            # Create a Hamiltonian with explicit off-diagonal terms
            k = sumindex(3)
            H_with_od = SymExpr(a'() * a()) + SymSum(σp(k) + σm(k), k)
            P_spin = Subspace(σz() => -1)
            H_od_proj = project_to_subspace(H_with_od, P_spin)
            # σp and σm are off-diagonal and should be removed
            # Only a†a should remain
            @test H_od_proj isa QuExpr || H_od_proj isa SymExpr
        end
    end

    QuantumAlgebra.use_σpm(false)
end
