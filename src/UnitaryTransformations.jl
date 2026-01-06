"""
    UnitaryTransformations

A Julia package for performing unitary transformations on quantum Hamiltonians,
including Schrieffer-Wolff transformations and Magnus expansions.

Built on top of QuantumAlgebra.jl for symbolic quantum operator algebra.
Uses Symbolics.jl for proper symbolic manipulation of energy denominators.
"""
module UnitaryTransformations

# Disable precompilation due to iszero method override
__precompile__(false)

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
    QuOpName

using Symbolics: Num, unwrap

# Override iszero for Symbolics.Num to avoid expensive polynomial computations
# that occur when QuantumAlgebra's normal_form checks for zero coefficients.
# This is safe because we conservatively return false for symbolic expressions,
# meaning we may keep some zero terms but won't incorrectly discard non-zero terms.
function Base.iszero(x::Num)
    u = unwrap(x)
    # Only return true for literal numeric zero
    if u isa Number
        return iszero(u)
    end
    # For symbolic expressions, conservatively say not zero
    # The final simplify_coefficients will handle actual zeros
    return false
end

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
