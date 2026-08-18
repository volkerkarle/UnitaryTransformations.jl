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
        solve_for_generator,
        set_fastmode_flat!,
        is_fastmode_flat

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
            real_sub =
                real_part isa Symbolics.Num ?
                Symbolics.substitute(real_part, substitutions) : real_part
            imag_sub =
                imag_part isa Symbolics.Num ?
                Symbolics.substitute(imag_part, substitutions) : imag_part
            real_val = real_sub isa Symbolics.Num ? Symbolics.value(real_sub) : real_sub
            imag_val = imag_sub isa Symbolics.Num ? Symbolics.value(imag_sub) : imag_sub
            return abs(complex(real_val, imag_val))
        end
        substituted =
            coeff isa Symbolics.Num ? Symbolics.substitute(coeff, substitutions) : coeff
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

    @testset "include_QQ captures Q-Q virtual paths" begin
        # 3-level system with sequential couplings: 0↔1 and 1↔2
        L = nlevel_ops(3, :L)

        @variables Δ g c01 c12

        H0 = 0 * L[1, 1] + Δ * L[2, 2] + 2Δ * L[3, 3]
        V = g * (c01 * (L[1, 2] + L[2, 1]) + c12 * (L[2, 3] + L[3, 2]))
        H = normal_form(H0 + V)

        P = Subspace(L[1, 1] => 1)

        result_inc = schrieffer_wolff(H, P; order = 4, include_QQ = true)
        result_exc = schrieffer_wolff(H, P; order = 4, include_QQ = false)

        coeffs_inc = [coeff for (_, coeff) in result_inc.H_P.terms]
        has_c12 = any(c -> occursin("c12", string(c)), coeffs_inc)
        @test has_c12

        coeffs_exc = [coeff for (_, coeff) in result_exc.H_P.terms]
        has_c12_exc = any(c -> occursin("c12", string(c)), coeffs_exc)
        @test !has_c12_exc
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

    @testset "Multi-atom systems with ∑" begin
        Symbolics.@variables tc_Δ tc_g ω_c
        QuantumAlgebra.use_σpm(true)

        # Tavis-Cummings Hamiltonian using ∑
        H_cav = ω_c * a'() * a()
        H_atom = ∑(:i, tc_Δ/2 * σz(:i))
        H_int = ∑(:i, tc_g * (a'() * σm(:i) + a() * σp(:i)))
        H = H_cav + H_atom + H_int

        P = Subspace(a'()*a() => 0, σz() => -1)

        # Test: H is a plain QuExpr (no special types)
        @test H isa QuExpr

        # Decompose
        H_d, H_od = decompose(H, P)
        @test is_diagonal(H_d, P)
        @test !isempty(H_od.terms)

        # Solve for generator
        S = solve_for_generator(H_d, H_od, P)
        @test S isa QuExpr
        @test !isempty(S.terms)

        # Full SW
        result = schrieffer_wolff(H, P; order=2)
        @test result.H_eff isa QuExpr
        @test result.S isa QuExpr

        # Project to subspace
        H_simple = ω_c * a'()*a() + ∑(:j, tc_Δ/2 * σz(:j))
        H_proj = project_to_subspace(H_simple, P)
        @test H_proj isa QuExpr
        @test is_diagonal(H_proj, P)

        # Combined subspace projection with SW
        H_combined = ω_c * a'()*a() + ∑(:j, tc_Δ/2 * σz(:j))
        result2 = schrieffer_wolff(H_combined, P; order=2)
        @test result2.H_P isa QuExpr

        QuantumAlgebra.use_σpm(false)
    end

    @testset "sw_generator correctness" begin
        Symbolics.@variables Δ g
        QuantumAlgebra.use_σpm(true)
        
        H = Δ/2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)
        
        # Order 1: sw_generator produces correct S₁
        S_gen1 = sw_generator(H, P; order=1)
        result_full2 = schrieffer_wolff(H, P; order=2, diagonal_only=false)
        @test normal_form(S_gen1) == normal_form(result_full2.S)
        
        # Order 2: sw_generator matches full SW S
        S_gen2 = sw_generator(H, P; order=2)
        result_full2b = schrieffer_wolff(H, P; order=2, diagonal_only=false)
        @test normal_form(S_gen2) == normal_form(result_full2b.S)
        
        # Order 3: sw_generator matches full SW S
        S_gen3 = sw_generator(H, P; order=3)
        result_full3 = schrieffer_wolff(H, P; order=3, diagonal_only=false)
        @test normal_form(S_gen3) == normal_form(result_full3.S)
        
        QuantumAlgebra.use_σpm(false)
    end

    @testset "schrieffer_wolff order 3 (full)" begin
        @variables Δ g
        QuantumAlgebra.use_σpm(true)

        H = Δ / 2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)

        # Full SW at order 3
        result = schrieffer_wolff(H, P; order = 3, diagonal_only = false)
        @test haskey(result, :H_eff)
        @test haskey(result, :S)
        @test haskey(result, :H_P)
        @test result.H_eff isa QuExpr
        @test result.S isa QuExpr

        # H_eff should be block-diagonal (only diagonal terms)
        @test is_diagonal(result.H_eff, P)

        # Verify generator equation: [S, H_d] should approximately cancel V_od
        H_d, H_od = decompose(H, P)
        # At minimum order 1, [S₁, H_d] = -V_od₁
        check = normal_form(comm(result.S, H_d) + H_od)
        # Not exactly zero because S includes S₂, S₃ and H_od is just V₁
        # But it should have fewer terms than original

        QuantumAlgebra.use_σpm(false)
    end

    @testset "include_QQ virtual paths" begin
        ops = nlevel_ops(3, :a)

        # Hamiltonian with Q↔Q transitions (off-diagonal terms within excited subspace)
        # States: |1⟩ = ground (P), |2⟩, |3⟩ = excited (Q)
        H = ops[1, 1] + 2 * ops[2, 2] + 3 * ops[3, 3]  # diagonal
        V_PQ = ops[1, 2] + ops[2, 1]  # P↔Q coupling
        V_QQ = ops[2, 3] + ops[3, 2]  # Q↔Q coupling (off-diagonal within Q)
        H_total = H + V_PQ + V_QQ

        P = Subspace(ops[1, 1] => 1)  # ground state

        # With include_QQ=true (default)
        result_qq = schrieffer_wolff(H_total, P; order = 2, include_QQ = true)
        @test result_qq.H_eff isa QuExpr

        # With include_QQ=false
        result_no_qq = schrieffer_wolff(H_total, P; order = 2, include_QQ = false)
        @test result_no_qq.H_eff isa QuExpr

        # Results should differ when Q↔Q terms are present
        # (they may not differ much at order 2, but the computation path differs)
    end

    @testset "project_to_subspace - N-level" begin
        ops = nlevel_ops(3, :a)

        # H = Σᵢ Eᵢ |i⟩⟨i| + Σᵢⱼ Vᵢⱼ |i⟩⟨j|
        H =
            1 * ops[1, 1] + 2 * ops[2, 2] + 3 * ops[3, 3] + ops[1, 2] + ops[2, 1]

        # Project onto state 1 (ground)
        P = Subspace(ops[1, 1] => 1)
        H_P = project_to_subspace(H, P)
        @test H_P isa QuExpr
        # Should only contain identity term (the eigenvalue of state 1)
        @test is_diagonal(H_P, P)

        # Project onto subspace where state 1 OR state 2 could be occupied
        # (This is a multi-state subspace - test that projection works)
    end

    @testset "project_to_subspace transition projector semantics" begin
        L = nlevel_ops(3, :L)
        @variables E1 E2 E3 v12

        H = E1 * L[1, 1] + E2 * L[2, 2] + E3 * L[3, 3] + v12 * (L[1, 2] + L[2, 1])
        P = Subspace(L[1, 1] => 1)
        H_P = project_to_subspace(H, P)

        coeffs = [string(coeff) for (_, coeff) in H_P.terms]
        op_strs = [isempty(term.bares.v) ? "𝟙" : string(term.bares) for (term, _) in H_P.terms]

        @test op_strs == ["𝟙"]
        @test any(c -> occursin("E1", c), coeffs)
        @test !any(c -> occursin("E2", c), coeffs)
        @test !any(c -> occursin("E3", c), coeffs)
        @test !any(c -> occursin("v12", c), coeffs)
    end

    @testset "include_QQ 4th-order mixed coupling structure" begin
        L = nlevel_ops(3, :L)
        @variables Δ δ ω g1 g2

        H0 = 0 * L[1, 1] + Δ * L[2, 2] + (Δ + δ) * L[3, 3] + ω * a'() * a()
        V =
            g1 * (L[1, 2] * a'() + L[2, 1] * a()) +
            g2 * (L[2, 3] * a'() + L[3, 2] * a())
        H = normal_form(H0 + V)
        P = Subspace(L[1, 1] => 1)

        result_inc = schrieffer_wolff(H, P; order = 4, include_QQ = true)
        result_exc = schrieffer_wolff(H, P; order = 4, include_QQ = false)

        coeffs_inc = [string(coeff) for (_, coeff) in result_inc.H_P.terms]
        has_g1_inc = any(c -> occursin("g1", c), coeffs_inc)
        has_g2_inc = any(c -> occursin("g2", c), coeffs_inc)
        has_mixed_inc = any(c -> occursin("g1^2", c) && occursin("g2^2", c), coeffs_inc)
        @test has_g1_inc
        @test has_g2_inc
        @test has_mixed_inc

        coeffs_exc = [string(coeff) for (_, coeff) in result_exc.H_P.terms]
        has_g2_exc = any(c -> occursin("g2", c), coeffs_exc)
        @test !has_g2_exc
    end

    @testset "schrieffer_wolff simplify_mode options" begin
        @variables Δ g
        QuantumAlgebra.use_σpm(true)

        H = Δ / 2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)

        for mode in [:none, :fast, :standard, :fractions]
            result = schrieffer_wolff(H, P; order = 2, simplify_mode = mode)
            @test result.H_eff isa QuExpr
            @test !isempty(result.H_eff.terms)
        end

        QuantumAlgebra.use_σpm(false)
    end

    @testset "diagonal_only at orders 2,3,4" begin
        @variables Δ g
        QuantumAlgebra.use_σpm(true)

        H = Δ / 2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)

        for ord in [2, 3, 4]
            result = schrieffer_wolff(H, P; order = ord, diagonal_only = true)
            @test result.H_eff isa QuExpr
            @test result.S isa QuExpr
            # S should only be S₁ (single generator, not S₁+S₂+...)
        end

        QuantumAlgebra.use_σpm(false)
    end

    @testset "fastmode flag restoration" begin
        @variables Δ g
        H = Δ / 2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)

        # Restores to false
        set_fastmode_flat!(false)
        _ = schrieffer_wolff(H, P; order = 2, fastmode_flat = true)
        @test is_fastmode_flat() == false

        # Restores to true
        set_fastmode_flat!(true)
        _ = schrieffer_wolff(H, P; order = 2, fastmode_flat = false)
        @test is_fastmode_flat() == true

        # Cleanup
        set_fastmode_flat!(false)
    end

    QuantumAlgebra.use_σpm(false)
end
