# Changelog

All notable changes to UnitaryTransformations.jl are documented in this file.

## Unreleased

### Added

- Added an explicitly graded Schrieffer-Wolff interface for Hamiltonians with
  perturbations at different physical orders:
  `schrieffer_wolff(H0, Dict(1 => V, 2 => W), P; order=4)`.
- Added physical-order BCH collection, including direct perturbations,
  perturbation-based nested commutators, and the corresponding unperturbed
  Hamiltonian contributions.

### Fixed

- Fixed incorrect fourth-order results from `fastmode_flat=true`. Deep nested
  commutators are now normal-ordered term by term without passing their large
  symbolic coefficients through `normal_form`.
- Resolve both transition contractions and bosonic commutators before
  decomposition, generator solving, and projection. This restores contributions
  such as the positive two-level fourth-order Kerr term while retaining the
  low-memory flat commutator path.
- Ensure projected flat-mode Hamiltonians are returned in canonical operator
  form instead of requiring a later external `normal_form` call.

### Tests

- Added fast-mode versus stepwise fourth-order Kerr regression coverage.
- Added equivalence coverage for the legacy and explicitly graded interfaces.
- Added second- and fourth-order dipole self-energy tests, including invariance
  of quartic photon channels under an order-two matter self-energy.
- Removed stale test-runner includes for test files that are no longer present.
