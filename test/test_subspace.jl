@testset "Subspace" begin
    using QuantumAlgebra
    using UnitaryTransformations:
        Subspace,
        OperatorConstraint,
        is_spin_constraint,
        is_number_constraint,
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
end
