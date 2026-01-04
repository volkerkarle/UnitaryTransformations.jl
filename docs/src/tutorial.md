# [Tutorial](@id tutorial)

This tutorial walks through the Schrieffer-Wolff transformation step by step.

## The Problem

Consider a quantum system with Hamiltonian:

```math
H = H_0 + V
```

where ``H_0`` has well-separated energy eigenspaces and ``V`` is a perturbation that couples them. We want to find an effective Hamiltonian ``H_{\text{eff}}`` that:

1. Acts only within a chosen low-energy subspace ``P``
2. Captures the effects of ``V`` to a given order in perturbation theory

## Step 1: Set Up the System

Let's work with a concrete example: a two-level system (qubit) coupled to a harmonic oscillator.

```julia
using UnitaryTransformations
using QuantumAlgebra
using Symbolics

# Use σ± basis - this is important for SW to work correctly!
QuantumAlgebra.use_σpm(true)

# Define symbolic parameters
@variables Δ g  # Δ = qubit-oscillator detuning, g = coupling strength

# Jaynes-Cummings Hamiltonian (in rotating frame)
H = Δ/2 * σz() + g * (a'()*σm() + a()*σp())
```

The Hamiltonian describes:
- A qubit with splitting ``Δ`` (the ``σ_z`` term)
- Coupling to an oscillator mode (the ``a^\dagger σ^-`` and ``a σ^+`` terms)

## Step 2: Define the Subspace

We need to specify which states belong to the low-energy subspace ``P``. For this example, we choose the qubit ground state:

```julia
# P = states where σz = -1 (qubit in ground state |g⟩)
P = Subspace(σz() => -1)
```

The `Subspace` type specifies expectation values of operators in the subspace. Here, we say that in subspace ``P``, the operator ``σ_z`` has eigenvalue ``-1``.

## Step 3: Decompose the Hamiltonian

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

- **Diagonal** (``H_d``): Operators that don't change the subspace (like ``σ^+σ^-``, ``a^\dagger a``)
- **Off-diagonal** (``V_{od}``): Operators that connect ``P`` and ``Q`` subspaces (like ``σ^+``, ``σ^-``)

## Step 4: Perform the Transformation

Now we apply the Schrieffer-Wolff transformation:

```julia
result = schrieffer_wolff(H, P; order=2)
```

This returns a named tuple with:
- `result.S` - The generator of the unitary transformation
- `result.H_eff` - The block-diagonal effective Hamiltonian
- `result.H_P` - The effective Hamiltonian projected onto subspace ``P``

## Step 5: Analyze the Results

### The Generator

```julia
println("Generator S = ", result.S)
# S = (g/Δ) a†σ⁻ + (-g/Δ) a σ⁺
```

The generator ``S`` is anti-Hermitian (``S^\dagger = -S``) and satisfies the fundamental equation:

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

Use the utility functions to extract specific coefficients:

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
println("Expected: -g²/Δ ✓")
```

## Key Points

1. **Always use `QuantumAlgebra.use_σpm(true)`** for SW transformations with spins
2. **Define subspace carefully** - this determines what "diagonal" means
3. **Use `collect_terms`** to see simplified coefficients
4. **The physics is in the coefficients** - extract them with `extract_coefficient`
