"""
    Decomposition of operators into diagonal and off-diagonal parts.

Given a Subspace P, decompose any operator H into:
- H_d: block-diagonal part (P·H·P + Q·H·Q)
- H_od: off-block-diagonal part (P·H·Q + Q·H·P)
"""

export decompose, diagonal_part, off_diagonal_part, is_diagonal, is_off_diagonal
export classify_operator, OperatorClass

using QuantumAlgebra
using QuantumAlgebra:
    QuExpr,
    QuTerm,
    BaseOperator,
    BaseOpProduct,
    BaseOpType,
    TLSx_,
    TLSy_,
    TLSz_,
    TLSCreate_,
    TLSDestroy_,
    BosonCreate_,
    BosonDestroy_,
    FermionCreate_,
    FermionDestroy_,
    normal_form,
    QuOpName

import ..UnitaryTransformations:
    is_spin_constraint, is_number_constraint, get_spin_constraint_info

"""
    OperatorClass

Classification of an operator with respect to a subspace:
- `DIAGONAL`: Operator preserves the subspace (stays in P or stays in Q)
- `RAISING`: Operator takes P → Q (increases quantum number)
- `LOWERING`: Operator takes Q → P (decreases quantum number)
- `MIXED`: Operator has both diagonal and off-diagonal components
"""
@enum OperatorClass DIAGONAL RAISING LOWERING MIXED

"""
    classify_base_operator(op::BaseOperator, P::Subspace)

Classify a single base operator with respect to the subspace P.
Returns (class, affected_constraint_index) where class is the OperatorClass
and affected_constraint_index indicates which constraint is affected (0 if none/transparent).
"""
function classify_base_operator(op::BaseOperator, P::Subspace)
    t = op.t

    for (i, constraint) in enumerate(P.constraints)
        # Check if this is a spin constraint
        spin_info = get_spin_constraint_info(constraint)
        if spin_info !== nothing
            spin_name, spin_inds, is_spin_down = spin_info

            # Check if operator has same name and indices as the constraint
            if op.name == spin_name && op.inds == spin_inds
                # σ+ raises spin: takes ↓ to ↑ (P → Q if P is spin down)
                if t == TLSCreate_
                    return is_spin_down ? RAISING : LOWERING
                    # σ- lowers spin: takes ↑ to ↓ (Q → P if P is spin down)
                elseif t == TLSDestroy_
                    return is_spin_down ? LOWERING : RAISING
                    # σx and σy are off-diagonal in σz basis (they flip spin)
                elseif t == TLSx_ || t == TLSy_
                    return MIXED  # Contains both raising and lowering
                # σz is diagonal
                elseif t == TLSz_
                    return DIAGONAL
                end
            end
        end

        # For number constraints (a†a eigenvalue)
        if is_number_constraint(constraint)
            constraint_ops = first(constraint.operator.terms)[1].bares.v
            constraint_name = constraint_ops[1].name
            constraint_inds = constraint_ops[1].inds

            # Check if operator has same name and indices
            if op.name == constraint_name && op.inds == constraint_inds
                # a† creates a boson: increases n
                if t == BosonCreate_
                    return constraint.eigenvalue == 0 ? RAISING : MIXED
                    # a annihilates a boson: decreases n (but can't go below 0)
                elseif t == BosonDestroy_
                    return constraint.eigenvalue == 0 ? LOWERING : MIXED
                    # f† creates a fermion
                elseif t == FermionCreate_
                    return constraint.eigenvalue == 0 ? RAISING : MIXED
                    # f annihilates a fermion
                elseif t == FermionDestroy_
                    return constraint.eigenvalue == 0 ? LOWERING : MIXED
                end
            end
        end
    end

    # Operator doesn't affect any constraint - it's "transparent" (diagonal)
    return DIAGONAL
end

"""
    classify_term(term::QuTerm, P::Subspace)

Classify a QuTerm based on how its bare operators affect the subspace.
A term is diagonal if its net effect on the subspace is zero (stays in P or Q).
A term is off-diagonal if it has a net raising or lowering effect.
"""
function classify_term(term::QuTerm, P::Subspace)
    ops = term.bares.v

    # Track net raising/lowering count PER constraint (per degree of freedom)
    # A term is diagonal if the net effect on ALL constraints is zero

    has_mixed = false
    constraint_effects = Dict{Int,Int}()  # constraint_index => net_raising (+1 for raising, -1 for lowering)

    for op in ops
        class = classify_base_operator(op, P)

        if class == MIXED
            has_mixed = true
        elseif class == RAISING
            # Find which constraint this affects
            idx = find_affected_constraint(op, P)
            constraint_effects[idx] = get(constraint_effects, idx, 0) + 1
        elseif class == LOWERING
            idx = find_affected_constraint(op, P)
            constraint_effects[idx] = get(constraint_effects, idx, 0) - 1
        end
        # DIAGONAL operators don't affect any constraint
    end

    # If any operator is truly MIXED (like σx), the term is mixed
    if has_mixed
        return MIXED
    end

    # Check net effects on each constraint
    # If any constraint has non-zero net effect, the term is off-diagonal
    for (idx, net_effect) in constraint_effects
        if net_effect > 0
            return RAISING
        elseif net_effect < 0
            return LOWERING
        end
    end

    # All net effects are zero -> diagonal
    return DIAGONAL
end

"""
    find_affected_constraint(op::BaseOperator, P::Subspace)

Find which constraint index the operator affects.
Returns 0 if no constraint is affected (transparent).
"""
function find_affected_constraint(op::BaseOperator, P::Subspace)
    t = op.t

    for (i, constraint) in enumerate(P.constraints)
        # Check spin constraints
        spin_info = get_spin_constraint_info(constraint)
        if spin_info !== nothing
            spin_name, spin_inds, _ = spin_info
            if op.name == spin_name && op.inds == spin_inds
                if t in (TLSCreate_, TLSDestroy_, TLSx_, TLSy_, TLSz_)
                    return i
                end
            end
        end

        # Check number constraints
        if is_number_constraint(constraint)
            constraint_ops = first(constraint.operator.terms)[1].bares.v
            constraint_name = constraint_ops[1].name
            constraint_inds = constraint_ops[1].inds

            if op.name == constraint_name && op.inds == constraint_inds
                if t in (BosonCreate_, BosonDestroy_, FermionCreate_, FermionDestroy_)
                    return i
                end
            end
        end
    end

    return 0  # No constraint affected
end

"""
    is_diagonal(expr::QuExpr, P::Subspace)

Check if an expression is purely block-diagonal with respect to subspace P.
"""
function is_diagonal(expr::QuExpr, P::Subspace)
    for (term, _) in expr.terms
        class = classify_term(term, P)
        if class != DIAGONAL
            return false
        end
    end
    return true
end

"""
    is_off_diagonal(expr::QuExpr, P::Subspace)

Check if an expression is purely off-block-diagonal with respect to subspace P.
"""
function is_off_diagonal(expr::QuExpr, P::Subspace)
    for (term, _) in expr.terms
        class = classify_term(term, P)
        if class == DIAGONAL
            return false
        end
    end
    return true
end

"""
    diagonal_part(H::QuExpr, P::Subspace)

Extract the block-diagonal part of H with respect to subspace P.
This includes terms that preserve P (P·H·P) and terms that preserve Q (Q·H·Q).

Note: Terms classified as MIXED (like σx when not using σpm mode) are skipped.
For proper handling, ensure Hamiltonians are expressed in the σ± basis.
"""
function diagonal_part(H::QuExpr, P::Subspace)
    result = QuExpr()
    for (term, coeff) in H.terms
        class = classify_term(term, P)
        if class == DIAGONAL
            result = result + coeff * QuExpr(term)
        end
        # MIXED terms are skipped - they need to be expanded at a higher level
    end
    return normal_form(result)
end

"""
    off_diagonal_part(H::QuExpr, P::Subspace)

Extract the off-block-diagonal part of H with respect to subspace P.
This includes terms that couple P to Q (P·H·Q + Q·H·P).

Note: Terms classified as MIXED (like σx when not using σpm mode) are skipped.
For proper handling, ensure Hamiltonians are expressed in the σ± basis.
"""
function off_diagonal_part(H::QuExpr, P::Subspace)
    result = QuExpr()
    for (term, coeff) in H.terms
        class = classify_term(term, P)
        if class == RAISING || class == LOWERING
            result = result + coeff * QuExpr(term)
        end
        # MIXED terms are skipped
    end
    return normal_form(result)
end

"""
    decompose(H::QuExpr, P::Subspace)

Decompose H into diagonal and off-diagonal parts with respect to subspace P.

Returns (H_d, H_od) where:
- H_d is the block-diagonal part
- H_od is the off-block-diagonal part
- H = H_d + H_od
"""
function decompose(H::QuExpr, P::Subspace)
    H_d = diagonal_part(H, P)
    H_od = off_diagonal_part(H, P)
    return (H_d, H_od)
end
