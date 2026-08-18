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
        P = Subspace(a'() * a() => 0)
        @test length(get_constraints(P)) == 1
        @test get_constraints(P)[1].eigenvalue == 0

        # Product state: spin down AND vacuum
        P = Subspace(σz() => -1, a'() * a() => 0)
        @test length(get_constraints(P)) == 2
    end

    @testset "Constraint classification" begin
        # Spin constraint
        c_spin = OperatorConstraint(σz(), -1)
        @test is_spin_constraint(c_spin)
        @test !is_number_constraint(c_spin)

        # Number constraint (bosonic)
        c_num = OperatorConstraint(a'() * a(), 0)
        @test !is_spin_constraint(c_num)
        @test is_number_constraint(c_num)
    end

    @testset "Indexed subspaces" begin
        # Indexed spin
        P = Subspace(σz(:i) => -1)
        @test length(get_constraints(P)) == 1

        # Multiple indexed modes
        P = Subspace(σz(:i) => -1, a'(:k) * a(:k) => 0)
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

    @testset "get_spin_constraint_info" begin
        using UnitaryTransformations: get_spin_constraint_info, OperatorConstraint, is_spin_constraint
        QuantumAlgebra.use_σpm(true)
        c1 = OperatorConstraint(σz(), -1)
        info1 = get_spin_constraint_info(c1)
        @test info1 !== nothing
        name1, inds1, is_down1 = info1
        @test is_down1 == true
        c2 = OperatorConstraint(σz(), 1)
        info2 = get_spin_constraint_info(c2)
        @test info2 !== nothing
        _, _, is_down2 = info2
        @test is_down2 == false
        QuantumAlgebra.use_σpm(false)
    end

    @testset "is_spin_constraint edge cases" begin
        # Non-spin operators
        c = UnitaryTransformations.OperatorConstraint(a'() * a(), 0)
        @test !UnitaryTransformations.is_spin_constraint(c)
        
        # Boson creation operator is not a spin constraint
        c2 = UnitaryTransformations.OperatorConstraint(a'(), 1)
        @test !UnitaryTransformations.is_spin_constraint(c2)
    end

    @testset "Fermion number constraints" begin
        using UnitaryTransformations: OperatorConstraint, is_number_constraint
        QuantumAlgebra.@fermion_ops e
        c = OperatorConstraint(e'() * e(), 1)
        @test is_number_constraint(c)
        c2 = OperatorConstraint(e'() * e(), 0)
        @test is_number_constraint(c2)
        c3 = OperatorConstraint(e'(), 1)
        @test !is_number_constraint(c3)
    end

    @testset "Transition operator constraints" begin
        ops = nlevel_ops(3, :a)
        
        # |1⟩⟨1| projector → eigenvalue 1
        c = UnitaryTransformations.OperatorConstraint(ops[1, 1], 1)
        @test UnitaryTransformations.is_transition_constraint(c)
        
        info = UnitaryTransformations.get_transition_constraint_info(c)
        @test info !== nothing
        @test info.state == 1
        
        # Non-projector transition operator (|1⟩⟨2|) is not a constraint
        c2 = UnitaryTransformations.OperatorConstraint(ops[1, 2], 1)
        @test !UnitaryTransformations.is_transition_constraint(c2)
        
        # Bosonic operator is not a transition constraint
        c3 = UnitaryTransformations.OperatorConstraint(a'() * a(), 0)
        @test !UnitaryTransformations.is_transition_constraint(c3)
        @test UnitaryTransformations.get_transition_constraint_info(c3) === nothing
    end

    @testset "Lie algebra constraint info" begin
        gens = su_generators(2, :λ)
        
        # λ₃ (diagonal, σz/2) constraint
        c = UnitaryTransformations.OperatorConstraint(gens[3], -1//2)
        @test UnitaryTransformations.is_lie_algebra_constraint(c)
        
        info = UnitaryTransformations.get_lie_algebra_constraint_info(c)
        @test info !== nothing
        
        # Non-diagonal Lie generator is not a constraint
        c2 = UnitaryTransformations.OperatorConstraint(gens[1], 1//2)
        @test !UnitaryTransformations.is_lie_algebra_constraint(c2)
    end

    @testset "Subspace edge cases" begin
        # Empty subspace
        P_empty = UnitaryTransformations.Subspace()
        @test length(P_empty.constraints) == 0
        
        # Multiple independent constraints
        P_multi = UnitaryTransformations.Subspace(σz() => -1, a'() * a() => 0)
        @test length(P_multi.constraints) == 2
        
        # Indexed operators
        P_idx = UnitaryTransformations.Subspace(σz(:i) => -1, σz(:j) => 1)
        @test length(P_idx.constraints) == 2
    end
end  # @testset "Subspace"
