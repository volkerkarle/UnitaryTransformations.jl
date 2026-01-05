# [Tutorial](@id tutorial)

This tutorial walks through the Schrieffer-Wolff transformation step by step. For the mathematical background, see the [Theory](@ref theory) page.

## Overview

The Schrieffer-Wolff (SW) transformation finds an effective Hamiltonian that acts within a chosen subspace by perturbatively eliminating couplings to other subspaces. This is useful when:

- You have a system with well-separated energy scales
- You want to derive an effective low-energy theory
- You need analytical expressions for perturbative corrections

## Step 1: Set Up the System

Let's work with a concrete example: a two-level system (qubit) coupled to a harmonic oscillator (the Jaynes-Cummings model).

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Use σ± basis - important for SW to work correctly with spins
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ g  # Δ = qubit-oscillator detuning, g = coupling strength

# Jaynes-Cummings Hamiltonian (in rotating frame)
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())
```

The Hamiltonian describes:
- A qubit with splitting ``\Delta`` (the ``\sigma_z`` term)
- Coupling to an oscillator mode (the ``a^\dagger \sigma^-`` and ``a \sigma^+`` terms)

## Step 2: Define the Subspace

We need to specify which states belong to the low-energy subspace ``P``. For this example, we choose the qubit ground state:

```julia
# P = states where σz = -1 (qubit in ground state |g⟩)
P = Subspace(σz() => -1)
```

The `Subspace` type specifies expectation values of operators in the subspace. Here, we say that in subspace ``P``, the operator ``\sigma_z`` has eigenvalue ``-1``.

### Multiple Constraints

For more complex systems, you can specify multiple constraints:

```julia
# Subspace with qubit ground state AND zero photons
P = Subspace(σz() => -1, a'()*a() => 0)
```

## Step 3: Understand the Decomposition

The SW transformation requires splitting ``H`` into diagonal and off-diagonal parts with respect to ``P``:

```julia
H_d, V_od = decompose(H, P)

println("Diagonal:     ", H_d)
println("Off-diagonal: ", V_od)
```

Output:
```
Diagonal:     -0.5Δ + Δ σ⁺σ⁻
Off-diagonal: g a†σ⁻ + g a σ⁺
```

- **Diagonal** (``H_d``): Operators that don't change the subspace (like ``\sigma^+\sigma^-``, ``a^\dagger a``)
- **Off-diagonal** (``V_{od}``): Operators that connect ``P`` and ``Q`` subspaces (like ``\sigma^+``, ``\sigma^-``)

## Step 4: Perform the Transformation

Now apply the Schrieffer-Wolff transformation:

```julia
result = schrieffer_wolff(H, P; order=2)
```

This returns a named tuple with:
- `result.S` - The generator of the unitary transformation ``e^S``
- `result.H_eff` - The block-diagonal effective Hamiltonian
- `result.H_P` - The effective Hamiltonian projected onto subspace ``P``

### Higher Orders

You can go to higher orders for more accuracy:

```julia
result_4th = schrieffer_wolff(H, P; order=4)
```

Note: Higher orders produce more complex expressions and take longer to compute.

## Step 5: Analyze the Results

### The Generator

```julia
println("Generator S = ", result.S)
# S = (g/Δ) a†σ⁻ + (-g/Δ) a σ⁺
```

The generator ``S`` is anti-Hermitian (``S^\dagger = -S``) and satisfies:

```math
[S, H_d] = -V_{od}
```

### The Effective Hamiltonian

```julia
# Collect and display all terms with simplified coefficients
terms = collect_terms(result.H_eff)
for (op, coeff) in terms
    println("  ", coeff, "  ", op)
end
```

Output:
```
  -0.5Δ        𝟙
  -(g²)/Δ      a†a
  Δ + (g²)/Δ   σ⁺σ⁻
  ...
```

### The Projected Hamiltonian

For many applications, we only care about the subspace ``P``:

```julia
println("H_P = ", result.H_P)
# H_P = -0.5Δ + (-(g²)/Δ) a†a
```

This is the **dispersive Hamiltonian**: the cavity frequency is shifted by ``-g^2/\Delta`` when the qubit is in the ground state!

## Step 6: Extract Physical Parameters

Use utility functions to extract specific coefficients:

```julia
# Get the dispersive shift (coefficient of a†a)
χ = extract_coefficient(result.H_P, a'()*a())
println("Dispersive shift χ = ", χ)
# Output: -(g²)/Δ
```

## Step 7: Numerical Evaluation

Substitute numerical values to get concrete numbers:

```julia
H_numeric = substitute_values(result.H_P, Dict(:g => 0.1, :Δ => 1.0))
println("H_P with g=0.1, Δ=1.0: ", H_numeric)
```

## Step 8: LaTeX Output

For publications and documentation, you can output results in LaTeX:

```julia
# Convert a single expression to LaTeX
println(to_latex(result.H_P))
# Output: - \frac{1}{2} \Delta + \frac{-g^{2}}{\Delta} {a}^{\dagger} {a}

# Pretty-print with a name
print_latex(result.H_P; name="H_P")
# Output: H_P = - \frac{1}{2} \Delta + \frac{-g^{2}}{\Delta} {a}^{\dagger} {a}

# Show all components of the result
show_result(result)
```

The `show_result` function prints the generator ``S``, effective Hamiltonian ``H_{\text{eff}}``, and projected Hamiltonian ``H_P`` in LaTeX format.

## Complete Example

Here's the full code:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Setup
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters and Hamiltonian
@variables Δ g
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())

# Define subspace and transform
P = Subspace(σz() => -1)
result = schrieffer_wolff(H, P; order=2)

# Analyze results
println("Effective Hamiltonian in ground state subspace:")
for (op, coeff) in collect_terms(result.H_P)
    println("  ", coeff, "  ", op)
end

# Extract dispersive shift
χ = extract_coefficient(result.H_P, a'()*a())
println("\nDispersive shift: χ = ", χ)
println("Expected: -g²/Δ")
```

---

## N-Level Atoms

The package supports N-level atomic systems using QuantumAlgebra's `nlevel_ops`:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Create 5-level atom operators: σ[i,j] = |i⟩⟨j|
σ5 = nlevel_ops(5, :q)

# Define level energies and coupling
ω = [Symbolics.variable(Symbol("ω", i)) for i in 1:5]
@variables ωc g

# Hamiltonian: 5-level atom + cavity, coupling levels 1↔3
H = sum(ω[i] * σ5[i,i] for i in 1:5) + 
    ωc * a'()*a() + 
    g * (σ5[1,3] * a'() + σ5[3,1] * a())

# Zero-photon subspace
P = Subspace(a'()*a() => 0)

result = schrieffer_wolff(H, P; order=2)
println("Effective Hamiltonian:")
println(result.H_eff)
```

The result contains dispersive shifts and AC Stark corrections for all levels.

---

## SU(N) Systems

For systems described by SU(N) Lie algebras, the package automatically detects and uses the matrix-element method:

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# SU(3) generators (Gell-Mann matrices)
λ = su_generators(3, :λ)

@variables Δ ω g

# Three-level Lambda system
H = Δ * λ[8] + ω * λ[7] + g * λ[2]

# Subspace defined by λ₈ eigenvalue
P = Subspace(λ[8] => 0.5)

# Automatically uses matrix-element method for SU(3)
result = schrieffer_wolff(H, P; order=2)
```

---

## Key Points

1. **Always use `QuantumAlgebra.use_σpm(true)`** for SW transformations with spins
2. **Define subspace carefully** - this determines what "diagonal" means
3. **Use `collect_terms`** to see simplified coefficients
4. **The physics is in the coefficients** - extract them with `extract_coefficient`
5. **N-level and SU(N) systems** are automatically handled with appropriate methods

## Next Steps

- See [Examples](@ref examples) for complete physics applications
- See [Theory](@ref theory) for mathematical details
- See [API Reference](@ref api) for function documentation
