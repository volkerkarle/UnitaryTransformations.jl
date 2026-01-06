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

5. W. Magnus, "On the exponential solution of differential equations for a linear operator," *Comm. Pure Appl. Math.* **7**, 649 (1954).

6. S. Blanes et al., "The Magnus expansion and some of its applications," *Physics Reports* **470**, 151 (2009).

---

## [Magnus Expansion](@id magnus_expansion)

The **Magnus expansion** is a technique for solving time-dependent Schrödinger equations and computing effective Hamiltonians for periodically driven (Floquet) systems.

### The Problem

Consider a time-dependent Schrödinger equation:

```math
i\frac{\partial U}{\partial t} = H(t) U(t)
```

where ``H(t)`` is periodic with period ``T = 2\pi/\omega``. We want to find an effective time-independent Hamiltonian ``H_{\text{eff}}`` such that after one period:

```math
U(T) = e^{-i H_{\text{eff}} T}
```

### Fourier Representation

A periodic Hamiltonian can be written in Fourier form:

```math
H(t) = \sum_n H_n e^{in\omega t}
```

where the Hermiticity condition requires ``H_{-n} = H_n^\dagger``.

### The Magnus Series

The solution to the Schrödinger equation can be written as:

```math
U(t) = e^{\Omega(t)}
```

where ``\Omega(t)`` is the Magnus series:

```math
\Omega(t) = \Omega_1(t) + \Omega_2(t) + \Omega_3(t) + \cdots
```

The effective Hamiltonian is:

```math
H_{\text{eff}} = i\Omega(T)/T = \sum_{k=1}^\infty \Omega_k
```

### Orders of the Expansion

**First order (k=1):**
```math
\Omega_1 = H_0
```

The leading term is simply the time-averaged Hamiltonian.

**Second order (k=2):**
```math
\Omega_2 = \sum_{n>0} \frac{-[H_n, H_{-n}]}{n\omega}
```

This captures effects like the **Bloch-Siegert shift** from counter-rotating drive terms.

**Higher orders (k≥3):**

For order ``k``, the Magnus term involves nested commutators with ``k`` Fourier components:

```math
\Omega_k = \sum_{\substack{n_1,\ldots,n_k \\ \sum_j n_j = 0}} C(n_1,\ldots,n_k) \, [[\cdots[[H_{n_1}, H_{n_2}], H_{n_3}],\ldots], H_{n_k}]
```

The coefficients are:

```math
C(n_1,\ldots,n_k) = \frac{1}{\omega^{k-1} \prod_{j=1}^{k-1} s_j}
```

where ``s_j = n_1 + \cdots + n_j`` are partial sums.

### Resonance Condition

The sum ``\sum_j n_j = 0`` is the **resonance condition**. It ensures that only certain combinations of Fourier modes contribute to the effective Hamiltonian.

### Reducible Terms

A term is **reducible** if any intermediate partial sum ``s_j = 0``. Such terms are excluded because they factorize into products of lower-order terms.

### Example: Circularly Driven Qubit

For a qubit driven by a circular field:

```math
H(t) = \frac{\Delta}{2}\sigma_z + \frac{\Omega}{2}(e^{i\omega t}\sigma^+ + e^{-i\omega t}\sigma^-)
```

The Fourier modes are:
- ``H_0 = \frac{\Delta}{2}\sigma_z``
- ``H_1 = \frac{\Omega}{2}\sigma^+``
- ``H_{-1} = \frac{\Omega}{2}\sigma^-``

The second-order correction is:

```math
\Omega_2 = -\frac{[H_1, H_{-1}]}{\omega} = -\frac{\Omega^2}{4\omega}[\sigma^+, \sigma^-] = -\frac{\Omega^2}{4\omega}\sigma_z
```

This is the **Bloch-Siegert shift** — a correction to the qubit frequency proportional to ``\Omega^2/\omega``.

### Convergence

The Magnus expansion converges when:

```math
\int_0^T \|H(t)\| \, dt < \pi
```

For high-frequency driving (``\omega \gg \|H\|``), convergence is typically rapid.

### Physical Applications

| Application | Key Effect |
|-------------|------------|
| NMR | Rotating-wave corrections |
| Trapped ions | Micromotion effects |
| Circuit QED | AC Stark shifts |
| Floquet engineering | Artificial gauge fields |

### Comparison with Floquet Theory

The Magnus expansion provides a **high-frequency expansion** of Floquet theory. While exact Floquet diagonalization gives the full quasienergy spectrum, the Magnus expansion produces analytical expressions valid in the ``\omega \gg \|V\|`` regime.
