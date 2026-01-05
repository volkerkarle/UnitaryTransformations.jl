# Session Summary: 2026-01-04

## Documentation Overhaul Completed

### New Documentation Structure

Restructured the documentation to present UnitaryTransformations.jl as a general unitary transformations toolkit, with Schrieffer-Wolff as the currently implemented method and others (Lang-Firsov, Bogoliubov, Holstein-Primakoff) planned.

**Files Updated:**

| File | Changes |
|------|---------|
| `docs/make.jl` | Added `theory.md` to page structure |
| `docs/src/index.md` | Landing page (updated previously) |
| `docs/src/theory.md` | NEW - Mathematical background (BCH, SW derivation, methods) |
| `docs/src/tutorial.md` | Rewritten with N-level and SU(N) examples, links to theory |
| `docs/src/examples.md` | Added N-level atom + cavity, multi-level atoms, SU(3) Lambda system |
| `docs/src/api.md` | Added new functions, better organization |
| `README.md` | Complete rewrite with tables, expanded planned transformations |

### Documentation Highlights

- **Tutorial** now covers:
  - Step-by-step SW transformation
  - N-level atoms with `nlevel_ops()`
  - SU(N) systems with automatic detection
  
- **Examples** now includes:
  - Jaynes-Cummings dispersive regime
  - Two-level system with transverse field
  - Rabi model / Bloch-Siegert shift
  - **NEW**: 5-level atom + cavity
  - **NEW**: 7-level atom with multiple couplings
  - **NEW**: Three-level Lambda system (SU(3))

- **API Reference** documents:
  - `solve_for_generator_eigenoperator`, `solve_for_generator_lie`
  - `detect_lie_algebra_system`, `compute_energy_eigenvalues`
  - Parameter conversion utilities

---

## CI Fixes

### Problem
CI was failing because:
1. QuantumAlgebra v1.3.1 (registry) lacks `nlevel_ops`, `su_generators`, SU(N) support
2. `docs/Project.toml` had hardcoded absolute path to local UnitaryTransformations

### Solution
Used Julia's `[sources]` section in Project.toml to specify dependencies from non-registry locations:

**`Project.toml`:**
```toml
[sources]
QuantumAlgebra = {url = "https://github.com/volkerkarle/QuantumAlgebra.jl"}
```

**`docs/Project.toml`:**
```toml
[sources]
QuantumAlgebra = {url = "https://github.com/volkerkarle/QuantumAlgebra.jl"}
UnitaryTransformations = {path = ".."}
```

**`.github/workflows/CI.yml`:**
- Simplified - no longer needs manual `Pkg.add(url=...)` steps
- Tests run on Julia 1.10, 1.12, and nightly

---

## Code Formatting

Ran `JuliaFormatter.format(".")` to fix format check CI. 

**Files formatted:**
- `src/commutator_series.jl`
- `src/decompose.jl`
- `src/inverse_liouvillian.jl`
- `src/schrieffer_wolff.jl`
- `src/subspace.jl`
- `test/test_commutator_series.jl`
- `test/test_schrieffer_wolff.jl`
- `test/test_subspace.jl`

---

## Test Status

All **147 tests passing** after all changes.

---

## Files Changed This Session

```
.github/workflows/CI.yml          # Simplified, added Julia 1.10
Project.toml                      # Added [sources], lowered Julia compat to 1.10
README.md                         # Complete rewrite
docs/Project.toml                 # Fixed paths in [sources]
docs/make.jl                      # Added theory.md
docs/src/api.md                   # Updated with new functions
docs/src/examples.md              # Added N-level examples
docs/src/theory.md                # NEW file
docs/src/tutorial.md              # Rewritten
src/commutator_series.jl          # Formatted
src/decompose.jl                  # Formatted
src/inverse_liouvillian.jl        # Formatted
src/schrieffer_wolff.jl           # Formatted
src/subspace.jl                   # Formatted
test/test_commutator_series.jl    # Formatted
test/test_schrieffer_wolff.jl     # Formatted
test/test_subspace.jl             # Formatted
```

---

## Next Steps

1. **Push changes** and verify CI passes
2. **Future development:**
   - Lang-Firsov transformation (polaron)
   - Bogoliubov transformation (quadratic bosonic)
   - Holstein-Primakoff transformation (spin-to-boson mapping)
