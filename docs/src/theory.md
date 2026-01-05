# [Theory](@id theory)

This section provides the mathematical background for the unitary transformations implemented in this package.

## Unitary Transformations in Quantum Mechanics

A unitary transformation maps a Hamiltonian ``H`` to a new Hamiltonian ``\tilde{H}`` via:

```math
\tilde{H} = U H U^\dagger
```

where ``U`` is a unitary operator (``U^\dagger U = \mathbb{1}``). The physics is unchanged—the eigenvalues are preserved—but the representation may be simpler or more convenient.

When ``U = e^S`` for some anti-Hermitian generator ``S`` (i.e., ``S^\dagger = -S``), the transformed Hamiltonian can be expanded using the **Baker-Campbell-Hausdorff (BCH) formula**:

```math
e^S H e^{-S} = H + [S, H] + \frac{1}{2!}[S, [S, H]] + \frac{1}{3!}[S, [S, [S, H]]] + \cdots
```

This infinite series becomes tractable when ``S`` is small (perturbative) or when the series terminates after a finite number of terms.

---

## [Schrieffer-Wolff Transformation](@id sw_transformation)

The **Schrieffer-Wolff (SW) transformation** is a perturbative method for block-diagonalizing Hamiltonians with well-separated energy scales. It was introduced by Schrieffer and Wolff in 1966 to derive the Kondo exchange interaction from the Anderson impurity model.

### The Problem

Consider a Hamiltonian of the form:

```math
H = H_0 + V
```

where:
- ``H_0`` is the **unperturbed Hamiltonian** with known eigenstates grouped into low-energy (``P``) and high-energy (``Q``) sectors
- ``V`` is a **perturbation** that couples the ``P`` and ``Q`` sectors

We want to find an effective Hamiltonian ``H_{\text{eff}}`` that:
1. Acts only within the low-energy sector ``P``
2. Captures the effects of ``V`` to a given order in perturbation theory

### The Transformation

We seek a unitary ``U = e^S`` such that the transformed Hamiltonian:

```math
H_{\text{eff}} = e^S H e^{-S}
```

is **block-diagonal** with respect to the ``P`` and ``Q`` subspaces. This means ``H_{\text{eff}}`` has no matrix elements connecting ``P`` and ``Q``.

### Determining the Generator

Decompose the Hamiltonian and generator into block-diagonal and off-block-diagonal parts:

```math
H = H_d + V_{od}, \quad S = S_{od}
```

where:
- ``H_d = P H P + Q H Q`` (block-diagonal)
- ``V_{od} = P H Q + Q H P`` (off-block-diagonal)
- ``S_{od}`` is purely off-block-diagonal (anti-Hermitian)

The **generator equation** at first order is:

```math
[S, H_d] = -V_{od}
```

This is the fundamental equation that determines ``S``. It states that the commutator of ``S`` with the diagonal Hamiltonian must cancel the off-diagonal perturbation.

### Solving the Generator Equation

For operators ``O`` that are **eigenoperators** of the adjoint action of ``H_d``—meaning ``[H_d, O] = \varepsilon \cdot O`` for some energy ``\varepsilon``—the solution is:

```math
S = \sum_\alpha \frac{V_\alpha}{\varepsilon_\alpha}
```

where ``V_\alpha`` are the components of ``V_{od}`` and ``\varepsilon_\alpha`` are the corresponding energy denominators.

**Example**: For a two-level system with ``H_d = \frac{\Delta}{2}\sigma_z`` and ``V_{od} = g\,\sigma^+``:
- ``[\sigma_z, \sigma^+] = 2\sigma^+``, so ``[H_d, \sigma^+] = \Delta\,\sigma^+``
- Therefore ``S = \frac{g}{\Delta}\sigma^+`` (plus Hermitian conjugate for the ``\sigma^-`` term)

### The Effective Hamiltonian

Using the BCH expansion:

```math
H_{\text{eff}} = H_d + \frac{1}{2}[S, V_{od}] + O(V^3)
```

The key second-order contribution is:

```math
H^{(2)} = \frac{1}{2}[S, V_{od}]
```

This generates effective interactions within the low-energy sector that arise from virtual transitions to high-energy states.

### Physical Interpretation

The SW transformation captures the physics of **virtual processes**:

1. The system starts in the low-energy sector ``P``
2. The perturbation ``V`` virtually excites it to high-energy sector ``Q``
3. The system returns to ``P`` via another application of ``V``

This virtual excitation costs energy ``\Delta E`` and contributes to the effective Hamiltonian as ``\sim V^2/\Delta E``.

### Order-by-Order Expansion

At higher orders, the SW transformation proceeds iteratively:

| Order | Contribution | Energy dependence |
|-------|--------------|-------------------|
| 0 | ``H_d`` | Original diagonal |
| 2 | ``\frac{1}{2}[S_1, V]`` | ``\sim g^2/\Delta`` |
| 3 | ``\frac{1}{2}[S_1, [S_1, H_d]] + [S_2, V]`` | ``\sim g^3/\Delta^2`` |
| 4 | Higher nested commutators | ``\sim g^4/\Delta^3`` |

Each order adds terms suppressed by additional powers of ``g/\Delta``.

---

## Eigenoperator Method

This package implements two methods for solving the generator equation. The **eigenoperator method** works for operators that satisfy:

```math
[H_d, O] = \varepsilon \cdot O
```

Such operators are called **eigenoperators** of the Liouvillian ``\mathcal{L}_{H_d}(\cdot) = [H_d, \cdot]``.

### Examples of Eigenoperators

| System | Operator | Eigenvalue |
|--------|----------|------------|
| TLS: ``H_d = \frac{\Delta}{2}\sigma_z`` | ``\sigma^+`` | ``+\Delta`` |
| | ``\sigma^-`` | ``-\Delta`` |
| Cavity: ``H_d = \omega a^\dagger a`` | ``a^\dagger`` | ``+\omega`` |
| | ``a`` | ``-\omega`` |
| N-level: ``H_d = \sum_i E_i \vert i\rangle\langle i\vert`` | ``\vert i\rangle\langle j\vert`` | ``E_i - E_j`` |

For these operators, the generator is simply:

```math
S = \frac{O}{\varepsilon}
```

### Composite Operators

Products of eigenoperators are also eigenoperators with additive eigenvalues:

```math
[H_d, O_1 O_2] = (\varepsilon_1 + \varepsilon_2) O_1 O_2
```

**Example**: For the Jaynes-Cummings interaction ``a^\dagger \sigma^-``:
- ``[H_d, a^\dagger \sigma^-] = (\omega_c - \Delta) a^\dagger \sigma^-``
- Energy denominator: ``\omega_c - \Delta`` (the detuning)

---

## Matrix-Element Method for Lie Algebras

For SU(N) systems expressed in the Gell-Mann basis, the generators are **not eigenoperators** of the diagonal Hamiltonian. For example:

```math
[\lambda_8, \lambda_2] \neq c \cdot \lambda_2
```

Instead, the commutator produces a linear combination of generators.

### Cartan-Weyl Basis

The solution is to work in the **Cartan-Weyl basis**, where the off-diagonal generators are replaced by transition operators:

| Gell-Mann | Cartan-Weyl |
|-----------|-------------|
| ``\lambda_1, \lambda_4`` | ``E_{12} = \vert 1\rangle\langle 2\vert``, ``E_{21} = \vert 2\rangle\langle 1\vert`` |
| ``\lambda_2, \lambda_5`` | ``E_{13} = \vert 1\rangle\langle 3\vert``, ``E_{31} = \vert 3\rangle\langle 1\vert`` |
| ``\lambda_3, \lambda_6`` | ``E_{23} = \vert 2\rangle\langle 3\vert``, ``E_{32} = \vert 3\rangle\langle 2\vert`` |

The transition operators **are** eigenoperators:

```math
[H_d, \vert i\rangle\langle j\vert] = (E_i - E_j)\vert i\rangle\langle j\vert
```

### Algorithm

1. **Compute energy eigenvalues** ``E_i`` from the diagonal Hamiltonian
2. **Convert** ``V_{od}`` from Gell-Mann to Cartan-Weyl basis
3. **Apply inverse Liouvillian**: ``S_{ij} = V_{ij} / (E_i - E_j)``
4. **Convert** ``S`` back to Gell-Mann basis

This is implemented in `solve_for_generator_lie()`.

---

## General Method: Liouvillian Linear System

When operators are **not** eigenoperators of the adjoint action—i.e., when ``[H_d, O]`` produces a **linear combination** of operators rather than a scalar multiple of ``O``—a more general approach is needed.

### The Problem

We want to solve the generator equation:

```math
[S, H_d] = -V_{\text{od}}
```

Expand both ``S`` and ``V_{\text{od}}`` in a basis of operators ``\{O_1, O_2, \ldots, O_n\}``:

```math
S = \sum_j s_j O_j, \quad V_{\text{od}} = \sum_i v_i O_i
```

The commutator ``[O_j, H_d]`` may not be proportional to ``O_j``, but is a linear combination:

```math
[O_j, H_d] = \sum_i L_{ij} O_i
```

where ``L`` is the **Liouvillian matrix** (the matrix representation of the adjoint action ``\mathcal{L}_{H_d}``).

### The Linear System

Substituting into the generator equation:

```math
[S, H_d] = \sum_j s_j [O_j, H_d] = \sum_j s_j \sum_i L_{ij} O_i = \sum_i \left(\sum_j L_{ij} s_j\right) O_i
```

For this to equal ``-V_{\text{od}} = -\sum_i v_i O_i``, we need:

```math
\sum_j L_{ij} s_j = -v_i \quad \Rightarrow \quad L \cdot \mathbf{s} = -\mathbf{v}
```

This is a **linear system** that can be solved symbolically.

### Algorithm

The implementation in `solve_for_generator_general()`:

1. **Extract basis**: Collect all unique operator structures from ``V_{\text{od}}``
2. **Compute Liouvillian matrix**: For each basis operator ``O_j``, compute ``[O_j, H_d]`` and extract coefficients
3. **Solve linear system**: ``\mathbf{s} = -L^{-1}\mathbf{v}`` using Cramer's rule (small systems) or Gaussian elimination
4. **Reconstruct generator**: ``S = \sum_j s_j O_j``

### When to Use

The general method is automatically used when:
- The eigenoperator method fails to find a scalar eigenvalue
- Operators mix under commutation with ``H_d``
- You explicitly request it via `method=:general`

### Example

For coupled oscillators where ``[H_d, a_1]`` produces both ``a_1`` and ``a_2`` terms, the eigenoperator method fails but the general method handles it correctly by solving the coupled linear system.

---

## Energy Denominators

The energy denominators in SW transformations have important physical meaning:

### Resonance Condition

When an energy denominator approaches zero (``\varepsilon \to 0``), the perturbation theory breaks down. This indicates a **resonance** where the two sectors are no longer well-separated.

### Example: Dispersive Regime

In circuit QED, the Jaynes-Cummings Hamiltonian has energy denominator ``\Delta = \omega_q - \omega_c``. The dispersive approximation is valid when:

```math
|g| \ll |\Delta|
```

The effective Hamiltonian contains the dispersive shift:

```math
\chi = -\frac{g^2}{\Delta}
```

This diverges as ``\Delta \to 0`` (resonance), signaling the breakdown of the perturbative treatment.

---

## Comparison with Other Methods

### vs. Löwdin Partitioning

Löwdin partitioning (also called quasi-degenerate perturbation theory) achieves the same goal but works directly with the Hamiltonian matrix rather than through a unitary transformation. SW provides the explicit generator ``S``, which can be useful for understanding the transformation and computing other observables.

### vs. Adiabatic Elimination

Adiabatic elimination assumes fast variables equilibrate instantly. SW is more systematic and provides higher-order corrections, but requires a perturbative expansion.

### vs. Numerical Diagonalization

Numerical methods give exact eigenvalues but not analytical expressions. SW produces symbolic results like ``g^2/\Delta`` that reveal the parameter dependence and scaling.

---

## References

1. J. R. Schrieffer and P. A. Wolff, "Relation between the Anderson and Kondo Hamiltonians," *Phys. Rev.* **149**, 491 (1966).

2. S. Bravyi, D. P. DiVincenzo, and D. Loss, "Schrieffer-Wolff transformation for quantum many-body systems," *Ann. Phys.* **326**, 2793 (2011).

3. C. Cohen-Tannoudji, J. Dupont-Roc, and G. Grynberg, *Atom-Photon Interactions* (Wiley, 1998), Chapter 3.

4. M. Wagner, *Unitary Transformations in Solid State Physics* (North-Holland, 1986).
