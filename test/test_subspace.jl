@testset "Subspace" begin
    using QuantumAlgebra
    using UnitaryTransformations:
        Subspace,
        OperatorConstraint,
        is_spin_constraint,
        is_number_constraint,
        is_lie_algebra_constraint,
        get_lie_algebra_constraint_info,
        is_diagonal_lie_generator,
        is_off_diagonal_lie_generator,
        get_lie_generator_state_info,
        get_constraints

    @testset "Subspace construction" begin
        # Spin-down subspace
        P = Subspace(σz() => -1)
        @test length(get_constraints(P)) == 1
        @test get_constraints(P)[1].eigenvalue == -1

        # Vacuum subspace
        P = Subspace(a'()*a() => 0)
        @test length(get_constraints(P)) == 1
        @test get_constraints(P)[1].eigenvalue == 0

        # Product state: spin down AND vacuum
        P = Subspace(σz() => -1, a'()*a() => 0)
        @test length(get_constraints(P)) == 2
    end

    @testset "Constraint classification" begin
        # Spin constraint
        c_spin = OperatorConstraint(σz(), -1)
        @test is_spin_constraint(c_spin)
        @test !is_number_constraint(c_spin)

        # Number constraint (bosonic)
        c_num = OperatorConstraint(a'()*a(), 0)
        @test !is_spin_constraint(c_num)
        @test is_number_constraint(c_num)
    end

    @testset "Indexed subspaces" begin
        # Indexed spin
        P = Subspace(σz(:i) => -1)
        @test length(get_constraints(P)) == 1

        # Multiple indexed modes
        P = Subspace(σz(:i) => -1, a'(:k)*a(:k) => 0)
        @test length(get_constraints(P)) == 2
    end

    @testset "SU(3) Lie algebra constraints" begin
        # Create SU(3) generators
        λ = su_generators(3, :λ)

        # Test diagonal generator identification
        @testset "Diagonal generator identification" begin
            for i = 1:8
                term, _ = first(λ[i].terms)
                bare = term.bares.v[1]
                if i in (7, 8)
                    @test is_diagonal_lie_generator(bare)
                    @test !is_off_diagonal_lie_generator(bare)
                else
                    @test !is_diagonal_lie_generator(bare)
                    @test is_off_diagonal_lie_generator(bare)
                end
            end
        end

        # Test Lie algebra constraint
        @testset "Lie algebra constraint recognition" begin
            # Constraint on diagonal generator λ₈
            c_lie = OperatorConstraint(λ[8], 0.5)
            @test is_lie_algebra_constraint(c_lie)
            @test !is_spin_constraint(c_lie)
            @test !is_number_constraint(c_lie)

            # Off-diagonal generators should NOT be valid constraints
            c_offdiag = OperatorConstraint(λ[1], 0.0)
            @test !is_lie_algebra_constraint(c_offdiag)

            # Constraint on λ₇ (also diagonal)
            c_lie7 = OperatorConstraint(λ[7], -0.5)
            @test is_lie_algebra_constraint(c_lie7)
        end

        # Test constraint info extraction
        @testset "Lie algebra constraint info" begin
            c_lie = OperatorConstraint(λ[8], 0.5)
            info = get_lie_algebra_constraint_info(c_lie)
            @test info !== nothing
            @test Symbol(info.name) == :λ  # name is QuOpName, not Symbol
            @test info.gen_idx == 8
            @test info.eigenvalue == 0.5
        end

        # Test state coupling info for off-diagonal generators
        @testset "Off-diagonal generator state coupling" begin
            for i = 1:8
                term, _ = first(λ[i].terms)
                bare = term.bares.v[1]
                state_info = get_lie_generator_state_info(bare)

                if i == 1 || i == 4
                    @test state_info == (1, 2)
                elseif i == 2 || i == 5
                    @test state_info == (1, 3)
                elseif i == 3 || i == 6
                    @test state_info == (2, 3)
                else
                    @test state_info === nothing  # Diagonal generators
                end
            end
        end
    end

    @testset "SU(2) Lie algebra constraints" begin
        # Create SU(2) generators
        σ = su_generators(2, :σ)

        # Test diagonal generator identification
        @testset "SU(2) diagonal generator identification" begin
            for i = 1:3
                term, _ = first(σ[i].terms)
                bare = term.bares.v[1]
                if i == 3
                    @test is_diagonal_lie_generator(bare)
                else
                    @test !is_diagonal_lie_generator(bare)
                end
            end
        end

        # Test SU(2) constraint
        @testset "SU(2) constraint recognition" begin
            c_su2 = OperatorConstraint(σ[3], -0.5)  # Spin down eigenvalue
            @test is_lie_algebra_constraint(c_su2)
        end

        # Test state coupling for off-diagonal generators
        @testset "SU(2) off-diagonal generator state coupling" begin
            for i = 1:3
                term, _ = first(σ[i].terms)
                bare = term.bares.v[1]
                state_info = get_lie_generator_state_info(bare)

                if i in (1, 2)
                    @test state_info == (1, 2)
                else
                    @test state_info === nothing
                end
            end
        end
    end
end
