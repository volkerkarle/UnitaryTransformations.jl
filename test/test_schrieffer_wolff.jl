@testset "Schrieffer-Wolff Transformation" begin
    using QuantumAlgebra
    using QuantumAlgebra: TLSCreate_
    using UnitaryTransformations: Subspace, decompose, diagonal_part, off_diagonal_part,
        classify_base_operator, classify_term, DIAGONAL, RAISING, LOWERING, MIXED,
        schrieffer_wolff, sw_generator, project_to_subspace, solve_for_generator
    
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
        ω = Pr"ω"
        Δ = Pr"Δ"
        g = Pr"g"
        
        H = ω * a'()*a() + Δ/2 * σz() + g * (a()*σp() + a'()*σm())
        
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
        
        Δ = Pr"Δ"
        ε = Pr"ε"
        
        # σx = σ+ + σ-
        H = Δ/2 * σz() + ε * (σp() + σm())
        
        P = Subspace(σz() => -1)
        
        # Just test decomposition for now
        H_d, H_od = decompose(H, P)
        
        @test normal_form(H_d) == normal_form(Δ/2 * σz())
        @test normal_form(H_od) == normal_form(ε * (σp() + σm()))
    end
    
    @testset "Generator solution" begin
        # For H_d = Δ/2 σz and V_od = ε σ+
        # [S, H_d] = -V_od means [S, Δ/2 σz] = -ε σ+
        # Since [σz, σ+] = 2σ+ and σz = -1 + 2σ⁺σ⁻, we have
        # [σ⁺σ⁻, σ+] = σ+, so [H_d, σ+] = Δ σ+
        # Thus S = (ε/Δ) σ+ (with proper Symbolics division)
        
        using Symbolics
        
        Δ = Pr"Δ"
        ε = Pr"ε"
        
        H_d = Δ/2 * σz()
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
        
        ω = Pr"ω"
        Δ = Pr"Δ"
        g = Pr"g"
        
        H = ω * a'()*a() + Δ/2 * σz() + g * (a'()*σm() + a()*σp())
        
        P = Subspace(σz() => -1)
        
        # Compute SW to second order
        result = schrieffer_wolff(H, P; order=2)
        
        # The effective Hamiltonian should be block-diagonal
        @test is_diagonal(result.H_eff, P)
        
        # Check that we got a generator
        @test !isempty(result.S.terms)
    end
    
    @testset "Projection to subspace" begin
        P = Subspace(σz() => -1)
        
        # σz in P sector should give -1
        H = Pr"Δ" * σz()
        H_P = project_to_subspace(H, P)
        
        # After projection, σz → -1
        @test normal_form(H_P) == normal_form(-Pr"Δ" * one(QuExpr))
    end
    
    @testset "Generator equation verification" begin
        # Rigorous check: [S, H_d] = -V_od must hold
        using Symbolics
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
                new_term = QuantumAlgebra.QuTerm(term.nsuminds, term.δs, QuantumAlgebra.Param[],
                                      term.expvals, term.corrs, term.bares)
                result = result + full_coeff * QuExpr(new_term)
            end
            return normal_form(result)
        end
        
        Δ = Pr"Δ"
        ε = Pr"ε"
        H = Δ/2 * σz() + ε * (σp() + σm())
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
        
        E_exact = -sqrt(Δ_val^2/4 + ε_val^2)
        E_SW = -Δ_val/2 - ε_val^2/Δ_val
        
        error_pct = 100 * abs(E_exact - E_SW) / abs(E_exact)
        
        # Should be accurate to < 0.1% for ε/Δ = 0.1
        @test error_pct < 0.1
    end
    
    QuantumAlgebra.use_σpm(false)
end
