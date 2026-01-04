using UnitaryTransformations
using QuantumAlgebra
using Test

@testset "UnitaryTransformations.jl" begin
    @testset "Module loads correctly" begin
        @test isdefined(UnitaryTransformations, :comm)
        @test isdefined(UnitaryTransformations, :normal_form)
        @test isdefined(UnitaryTransformations, :Subspace)
        @test isdefined(UnitaryTransformations, :schrieffer_wolff)
        @test isdefined(UnitaryTransformations, :commutator_series)
    end

    @testset "QuantumAlgebra integration" begin
        # Basic commutator test: [a, a†] = 1
        @test normal_form(comm(a(), a'())) == one(QuExpr)

        # Bosonic operators
        @test normal_form(a() * a'()) == a'() * a() + 1
    end

    # Include specific test files
    include("test_subspace.jl")
    include("test_commutator_series.jl")
    include("test_schrieffer_wolff.jl")
end
