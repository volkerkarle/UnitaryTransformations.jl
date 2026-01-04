"""
    Subspace definition for Schrieffer-Wolff transformation.

A Subspace defines a low-energy sector P by specifying constraints on diagonal operators.
For example, `Subspace(σz() => -1)` defines the spin-down subspace.
"""

export Subspace, OperatorConstraint, get_spin_constraint_info

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
    comm

"""
    OperatorConstraint

A single constraint defining a sector: an operator and its eigenvalue in that sector.
"""
struct OperatorConstraint
    operator::QuExpr
    eigenvalue::Number
end

"""
    Subspace(constraints...)

Define a low-energy subspace P by specifying constraints on operators.

# Examples
```julia
# Spin-down subspace
P = Subspace(σz() => -1)

# Vacuum subspace (zero bosons)
P = Subspace(a'()*a() => 0)

# Product state: spin-down AND zero bosons
P = Subspace(σz() => -1, a'()*a() => 0)

# Indexed systems
P = Subspace(σz(:i) => -1)  # All spins down
```
"""
struct Subspace
    constraints::Vector{OperatorConstraint}

    function Subspace(pairs::Pair{<:QuExpr,<:Number}...)
        constraints = [OperatorConstraint(op, val) for (op, val) in pairs]
        new(constraints)
    end
end

# Allow construction with single constraint
Subspace(op::QuExpr, val::Number) = Subspace(op => val)

"""
    get_constraints(P::Subspace)

Return the list of operator constraints defining the subspace.
"""
get_constraints(P::Subspace) = P.constraints

"""
    is_spin_constraint(c::OperatorConstraint)

Check if a constraint is on a spin operator (σz or its σ± representation).
Handles both direct σz constraints and the expanded form -1 + 2σ⁺σ⁻.
"""
function is_spin_constraint(c::OperatorConstraint)
    op = c.operator

    # Case 1: Direct σz operator (when not in σpm mode)
    if length(op.terms) == 1
        term, coeff = first(op.terms)
        if coeff == 1 && length(term.bares.v) == 1 && term.bares.v[1].t == TLSz_
            return true
        end
    end

    # Case 2: σz in σpm mode: σz = -1 + 2σ⁺σ⁻
    # Check for two terms: constant -1 and 2σ⁺σ⁻
    if length(op.terms) == 2
        has_const = false
        has_pm = false
        spin_name = nothing
        spin_inds = nothing

        for (term, coeff) in op.terms
            ops = term.bares.v
            if isempty(ops) && coeff == -1
                has_const = true
            elseif length(ops) == 2 && coeff == 2
                if ops[1].t == TLSCreate_ && ops[2].t == TLSDestroy_
                    if ops[1].name == ops[2].name && ops[1].inds == ops[2].inds
                        has_pm = true
                        spin_name = ops[1].name
                        spin_inds = ops[1].inds
                    end
                end
            end
        end

        if has_const && has_pm
            return true
        end
    end

    return false
end

"""
    get_spin_constraint_info(c::OperatorConstraint)

Extract spin name and indices from a spin constraint.
Returns (name, inds, is_spin_down) or nothing if not a spin constraint.
"""
function get_spin_constraint_info(c::OperatorConstraint)
    is_spin_constraint(c) || return nothing

    op = c.operator

    # Case 1: Direct σz
    if length(op.terms) == 1
        term, coeff = first(op.terms)
        if coeff == 1 && length(term.bares.v) == 1 && term.bares.v[1].t == TLSz_
            spin_op = term.bares.v[1]
            # eigenvalue -1 means spin down, +1 means spin up
            return (spin_op.name, spin_op.inds, c.eigenvalue == -1)
        end
    end

    # Case 2: σz = -1 + 2σ⁺σ⁻ form
    if length(op.terms) == 2
        for (term, coeff) in op.terms
            ops = term.bares.v
            if length(ops) == 2 && coeff == 2
                if ops[1].t == TLSCreate_ && ops[2].t == TLSDestroy_
                    # eigenvalue -1 for σz means σ⁺σ⁻ = 0 (spin down)
                    return (ops[1].name, ops[1].inds, c.eigenvalue == -1)
                end
            end
        end
    end

    return nothing
end

"""
    is_number_constraint(c::OperatorConstraint)

Check if a constraint is on a number operator (a'*a or f'*f).
"""
function is_number_constraint(c::OperatorConstraint)
    op = c.operator
    # Check if it's a†a form
    length(op.terms) == 1 || return false
    term, coeff = first(op.terms)
    coeff == 1 || return false
    ops = term.bares.v
    length(ops) == 2 || return false
    # Check for creation followed by annihilation
    return (ops[1].t == BosonCreate_ && ops[2].t == BosonDestroy_) ||
           (ops[1].t == FermionCreate_ && ops[2].t == FermionDestroy_)
end

"""
    get_operator_name(op::BaseOperator)

Extract the name of an operator.
"""
get_operator_name(op::BaseOperator) = op.name

"""
    get_operator_indices(op::BaseOperator)

Extract the indices of an operator.
"""
get_operator_indices(op::BaseOperator) = op.inds

"""
    operators_match_indices(op1::BaseOperator, op2::BaseOperator)

Check if two operators have the same name and indices.
"""
function operators_match_indices(op1::BaseOperator, op2::BaseOperator)
    return op1.name == op2.name && op1.inds == op2.inds
end
