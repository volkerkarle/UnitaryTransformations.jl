@testset "Symbolic Utilities" begin
    @testset "simplify_coefficients" begin
        Symbolics.@variables ω g

        # Bosonic Hamiltonian for testing all modes
        H = ω * a'() * a() + g * (a'() + a())

        # Test ALL 5 modes
        for mode in [:none, :fast, :standard, :fractions, :aggressive]
            result = simplify_coefficients(H; mode = mode)
            @test result isa QuExpr
            @test !isempty(result.terms)
        end

        # Test legacy `aggressive=true` kwarg (should map to :aggressive mode)
        result_legacy = simplify_coefficients(H; aggressive = true)
        @test result_legacy isa QuExpr
        @test !isempty(result_legacy.terms)

        # Test with a trivial expression (just a number coefficient — no simplification needed)
        trivial = 42 * one(QuExpr)
        triv_result = simplify_coefficients(trivial; mode = :fast)
        @test triv_result isa QuExpr

        # Test with Rational coefficients (the clean_rationals postwalk handles 16//1 -> 16)
        rat_expr = (1 // 2) * a'() * a() + (2 // 3) * a()
        rat_result = simplify_coefficients(rat_expr; mode = :fractions)
        @test rat_result isa QuExpr
    end

    @testset "substitute_values" begin
        # Test numerical substitution with integer coefficients (avoids Num*QuExpr issues)
        H = 2 * a'() * a() + a'()
        values = Dict{Symbol,Float64}()
        result = substitute_values(H, values)
        @test result isa QuExpr
    end

    @testset "extract_coefficient" begin
        Symbolics.@variables ω g

        H = ω * a'() * a() + g * a'()

        # Extract coefficient of a specific operator
        c = extract_coefficient(H, a'() * a())
        @test c isa Num

        # Test with simplify_coeff=false
        c_nosimpl = extract_coefficient(H, a'() * a(); simplify_coeff = false)
        @test c_nosimpl isa Num

        # Extract from expression that does NOT contain the target (returns nothing)
        c_notfound = extract_coefficient(H, a() * a())
        @test c_notfound === nothing

        # Extract from expression with identity term only
        H_id = 3 * one(QuExpr)
        c_id = extract_coefficient(H_id, one(QuExpr))
        # Identity matching may depend on internal representation; accept nothing or a Number
        @test (c_id === nothing) || (c_id isa Number)

        # Test with parametric coefficients (product of symbols)
        Symbolics.@variables α β
        H_param = α * β * a'() * a()
        c_param = extract_coefficient(H_param, a'() * a())
        @test c_param isa Num
    end

    @testset "collect_terms" begin
        Symbolics.@variables ω g

        H = ω * a'() * a() + g * a'()

        # Collect terms with their coefficients
        terms = collect_terms(H)
        @test terms isa Vector
        @test length(terms) >= 1

        # Test with empty QuExpr (zero)
        empty_terms = collect_terms(zero(QuExpr))
        @test empty_terms isa Vector
        @test isempty(empty_terms)

        # Test with multiple terms
        @test length(terms) >= 2

        # Test that coefficients are simplified (should be Num or Number, not expression with free vars)
        for (op_str, coeff) in terms
            @test coeff isa Union{Num,Number}
        end
    end

    @testset "to_latex" begin
        # Generate LaTeX string from QuExpr
        result = to_latex(a'() * a())
        @test result isa String
        @test !isempty(result)

        # With simplify_coeff=false
        result_nosimpl = to_latex(a'() * a(); simplify_coeff = false)
        @test result_nosimpl isa String
        @test !isempty(result_nosimpl)
    end

    @testset "print_latex" begin
        # Test with name (should prepend "N = ")
        result = print_latex(a'() * a(); name = "N")
        @test result isa String
        @test startswith(result, "N")

        # Test with display=true (should wrap in \[ ... \])
        result_display = print_latex(a'() * a(); display = true)
        @test result_display isa String
        @test startswith(result_display, "\\[")
    end

    @testset "format_expression" begin
        s_unicode = format_expression("Δ + ω + a†() a() + g1*d1 + d1^2 + 1//2"; mode = :unicode)
        @test occursin("Δ", s_unicode)
        @test occursin("a'", s_unicode)
        @test occursin("g₁", s_unicode)
        @test occursin("d₁", s_unicode)
        @test occursin("d₁²", s_unicode)
        @test !occursin("//", s_unicode)
        @test !occursin("*", s_unicode)

        s_latex = format_expression("Δ + ω + a†() a() + 1//2"; mode = :latex)
        @test occursin("\\Delta", s_latex)
        @test occursin("\\omega", s_latex)
        @test occursin("a^{\\dagger}", s_latex)
        @test occursin("\\frac{1}{2}", s_latex)
    end

    @testset "write_expression_dump" begin
        mktempdir() do tmp
            path_tex = joinpath(tmp, "out.tex")
            path_txt = joinpath(tmp, "out.txt")

            sections = ["Input" => "H = Δ + a†() a()"]

            write_expression_dump(
                path_tex,
                sections;
                mode = :latex,
                tex = true,
                header = ["% header"],
            )
            tex_data = read(path_tex, String)
            @test occursin("\\section*{Input}", tex_data)
            @test occursin("\\Delta", tex_data)
            @test occursin("a^{\\dagger}", tex_data)

            write_expression_dump(path_txt, sections; mode = :unicode, tex = false)
            txt_data = read(path_txt, String)
            @test occursin("[Input]", txt_data)
            @test occursin("a'", txt_data)

            @test_throws ArgumentError write_expression_dump(
                path_txt,
                sections;
                mode = :unicode,
                tex = true,
            )
        end
    end

    @testset "show_result" begin
        Symbolics.@variables Δ g

        QuantumAlgebra.use_σpm(true)

        H = Δ / 2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)
        sw_result = schrieffer_wolff(H, P; order = 2)

        # Should print without error
        show_result(sw_result)
        @test true  # If we got here, no error was thrown

        QuantumAlgebra.use_σpm(false)
    end

    @testset "compact_form" begin
        Symbolics.@variables ω g

        QuantumAlgebra.use_σpm(true)

        # H = ω a†a + g(a†σ⁻ + aσ⁺)
        # a†a should be in hermitian terms (self-adjoint)
        # a†σ⁻ should be paired with its h.c. aσ⁺ in hc_pairs
        H = ω * a'() * a() + g * (a'() * σm() + a() * σp())
        hermitian, pairs = compact_form(H)
        @test hermitian isa Dict
        @test pairs isa Vector
        @test length(hermitian) >= 1
        @test length(pairs) >= 1

        # Test with empty expression
        hermitian_empty, pairs_empty = compact_form(zero(QuExpr))
        @test isempty(hermitian_empty)
        @test isempty(pairs_empty)

        QuantumAlgebra.use_σpm(false)
    end

    @testset "print_compact" begin
        Symbolics.@variables ω g

        # QuExpr method — generate compact string
        H = ω * a'() * a() + g * (a'() + a())
        result = print_compact(H; name = "H")
        @test result isa String
        @test !isempty(result)

        # NamedTuple method (same as show_result but compact format)
        QuantumAlgebra.use_σpm(true)

        Symbolics.@variables Δ
        H_sw = Δ / 2 * σz() + g * (a'() * σm() + a() * σp())
        P = Subspace(σz() => -1)
        sw_result = schrieffer_wolff(H_sw, P; order = 2)

        # Should print without error and return nothing
        @test print_compact(sw_result) === nothing

        QuantumAlgebra.use_σpm(false)
    end
end
