# TODO: Symbolic Sum Support for Exchange Interactions

## STATUS: COMPLETED ✓

The SymSum functionality has been implemented and is working correctly.

## Summary of Work Done

### QuantumAlgebra (upstream) - COMPLETED ✓
- [x] Implemented `SymSum` type with bound variable semantics
- [x] Implemented `replace_index(expr, old_idx, new_idx)` utility
- [x] Implemented commutator rules for `SymSum` with same-site and cross-site contributions
- [x] Handle double sums `Σᵢ Σⱼ` and constrained sums `Σᵢ≠ⱼ`
- [x] Add arithmetic operations (`+`, `*`, scalar multiplication)
- [x] Add pretty printing (fixed display bug with nested sums showing `Σ#₁≠#₁`)
- [x] File: `/home/vsilv/.nextcloud/src/QuantumAlgebra.jl/src/symbolic_sums.jl`

### UnitaryTransformations.jl - COMPLETED ✓
- [x] Extended `decompose` to work with `SymSum/SymExpr`
- [x] Extended `solve_for_generator` to work with `SymSum/SymExpr`
- [x] Added `schrieffer_wolff(H::SymExpr, P; ...)` method
- [x] Updated Tavis-Cummings example to use new API
- [x] Added test case verifying exchange interaction appears
- [x] All 259 tests passing

## What Was Fixed

### Problem
The Schrieffer-Wolff transformation did not correctly generate **exchange interactions** (e.g., `J_ij ∝ gᵢgⱼ σ⁺(i)σ⁻(j)`) when using symbolic sum indices.

Using `sumindex(1)` created a dummy index `#₁`, but when computing `[Σᵢ Sᵢ, Σⱼ Vⱼ]`, it only produced the `i=j` contribution because both operands shared the same index.

### Solution
Implemented `SymSum` type in QuantumAlgebra that correctly handles bound variable semantics:

```julia
i = sumindex(1)
V = SymSum(g * (a'()*σm(i) + a()*σp(i)), i)

# Commutator [Σᵢ Aᵢ, Σⱼ Bⱼ] now correctly produces:
# - Same-site: Σᵢ [Aᵢ, Bᵢ]
# - Cross-site: Σᵢ Σⱼ≠ᵢ [Aᵢ, Bⱼ]
```

The exchange terms `σ⁺(1)σ⁻(2) + σ⁺(2)σ⁻(1)` now appear correctly in the Tavis-Cummings effective Hamiltonian.

## Usage Example

```julia
using QuantumAlgebra, UnitaryTransformations, Symbolics
import QuantumAlgebra: sumindex, SymSum, SymExpr

@variables ω_c Δ g

i = sumindex(1)

# Build Tavis-Cummings Hamiltonian with symbolic sums
H = SymExpr(ω_c * a'()*a()) + 
    SymSum(Δ/2 * σz(i), i) + 
    SymSum(g * (a'()*σm(i) + a()*σp(i)), i)

# Run Schrieffer-Wolff transformation
P = Subspace(a'()*a() => 0)  # Zero photon sector
result = schrieffer_wolff(H, P; order=2)

# result.H_eff contains the exchange interaction terms!
```

## Remaining Limitations

1. **Order > 2**: The `schrieffer_wolff(H::SymExpr, ...)` method currently only fully supports order=2. Higher orders may miss some cross-site contributions.

2. **H_d with SymSum terms**: When the diagonal Hamiltonian contains symbolic sums (e.g., `Σᵢ (Δ/2)σz(i)`), only the scalar (QuExpr) part is used for computing energy denominators. A warning is issued.

3. **Project to subspace**: `project_to_subspace` is not yet implemented for SymExpr. The full H_eff is returned without projection.

## References

- Tavis-Cummings example: `examples/tavis_cummings.jl`
- QuantumAlgebra SymSum: `/home/vsilv/.nextcloud/src/QuantumAlgebra.jl/src/symbolic_sums.jl`
- Test case: `test/test_schrieffer_wolff.jl` (last testset)
