"""
    UnitaryTransformations

A Julia package for performing unitary transformations on quantum Hamiltonians,
including Schrieffer-Wolff transformations and Magnus expansions.

Built on top of QuantumAlgebra.jl for symbolic quantum operator algebra.
Uses Symbolics.jl for proper symbolic manipulation of energy denominators.
"""
module UnitaryTransformations

using QuantumAlgebra
using QuantumAlgebra:
    comm,
    normal_form,
    QuExpr,
    QuTerm,
    BaseOperator,
    BaseOpProduct,
    Param,
    TLSx_,
    TLSy_,
    TLSz_,
    TLSCreate_,
    TLSDestroy_,
    BosonCreate_,
    BosonDestroy_,
    FermionCreate_,
    FermionDestroy_,
    QuOpName,
    # Symbolic sum types for multi-atom/multi-site systems
    SymSum,
    SymProd,
    SymExpr,
    AbstractSymbolicAggregate,
    expand_symbolic,
    sumindex

using Symbolics: Num, unwrap

# Re-export commonly used QuantumAlgebra functions
export comm, normal_form

# Include submodules in dependency order
include("subspace.jl")
include("decompose.jl")
include("commutator_series.jl")
include("inverse_liouvillian.jl")
include("symbolic_utils.jl")  # Must come before schrieffer_wolff.jl (imports simplify_coefficients)
include("schrieffer_wolff.jl")
include("magnus.jl")

end # module UnitaryTransformations
