"""
    UnitaryTransformations

A Julia package for performing unitary transformations on quantum Hamiltonians,
including Schrieffer-Wolff and Lang-Firsov transformations.

Built on top of QuantumAlgebra.jl for symbolic quantum operator algebra.
Uses Symbolics.jl for proper symbolic manipulation of energy denominators.
"""
module UnitaryTransformations

using QuantumAlgebra
using QuantumAlgebra: comm, normal_form, QuExpr, QuTerm, BaseOperator, BaseOpProduct,
    Param, TLSx_, TLSy_, TLSz_, TLSCreate_, TLSDestroy_,
    BosonCreate_, BosonDestroy_, FermionCreate_, FermionDestroy_,
    QuOpName

# Re-export commonly used QuantumAlgebra functions
export comm, normal_form

# Include submodules in dependency order
include("subspace.jl")
include("decompose.jl")
include("commutator_series.jl")
include("inverse_liouvillian.jl")
include("schrieffer_wolff.jl")
include("symbolic_utils.jl")

end # module UnitaryTransformations
