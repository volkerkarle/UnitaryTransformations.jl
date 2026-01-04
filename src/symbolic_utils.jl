"""
    Symbolic utilities for working with Symbolics.jl coefficients in QuExpr.

Provides functions to simplify, substitute, and extract coefficients from
quantum expressions that use Symbolics.jl for symbolic manipulation.
"""

export simplify_coefficients, substitute_values, extract_coefficient, collect_terms

using QuantumAlgebra
using QuantumAlgebra: QuExpr, QuTerm, Param, normal_form

using Symbolics
using Symbolics: Num, simplify, simplify_fractions

"""
    simplify_coefficients(expr::QuExpr)

Simplify all Symbolics coefficients in a QuExpr.
Returns a new QuExpr with simplified coefficients.

Uses `simplify_fractions` which is much faster than full `simplify`
for the rational expressions that arise in Schrieffer-Wolff.
"""
function simplify_coefficients(expr::QuExpr)
    # Build result directly to avoid expensive iszero checks in + operator
    result_terms = Dict{QuTerm,Number}()
    for (term, coeff) in expr.terms
        if coeff isa Num
            # Use simplify_fractions - much faster for rational expressions
            simplified_coeff = simplify_fractions(coeff)
            result_terms[term] = simplified_coeff
        else
            result_terms[term] = coeff
        end
    end
    return QuExpr(result_terms)
end

"""
    substitute_values(expr::QuExpr, values::Dict{Symbol, Number})

Substitute numerical values for symbolic parameters in a QuExpr.

# Arguments
- `expr`: A QuExpr with symbolic coefficients
- `values`: Dict mapping parameter symbols to numerical values

# Example
```julia
H_P = result.H_P
values = Dict(:g => 0.1, :Δ => 1.0)
H_numeric = substitute_values(H_P, values)
```
"""
function substitute_values(expr::QuExpr, values::Dict{Symbol,T}) where {T<:Number}
    result = QuExpr()

    for (term, coeff) in expr.terms
        # Handle Symbolics coefficient
        new_coeff = coeff
        if coeff isa Num
            # Substitute in the Symbolics expression
            sub_dict = Dict(Symbolics.variable(k) => v for (k, v) in values)
            new_coeff = Symbolics.substitute(coeff, sub_dict)
            # Try to extract numerical value
            if Symbolics.isconstant(new_coeff)
                new_coeff = Symbolics.unwrap(new_coeff)
            end
        end

        # Handle remaining Params in the term
        remaining_params = Param[]
        param_factor = 1
        for p in term.params
            if haskey(values, p.name)
                param_factor *= values[p.name]
            else
                push!(remaining_params, p)
            end
        end

        # Create new term with remaining params
        new_term = QuTerm(
            term.nsuminds,
            term.δs,
            remaining_params,
            term.expvals,
            term.corrs,
            term.bares,
        )

        result = result + (new_coeff * param_factor) * QuExpr(new_term)
    end

    return normal_form(result)
end

"""
    extract_coefficient(expr::QuExpr, target_ops::QuExpr)

Extract the coefficient of a specific operator structure from a QuExpr.

# Arguments
- `expr`: The full QuExpr to search
- `target_ops`: The operator structure to match (e.g., `a'()*a()`)

# Returns
- The coefficient (Num or Number) if found, nothing otherwise
"""
function extract_coefficient(expr::QuExpr, target_ops::QuExpr)
    target_norm = normal_form(target_ops)

    # Get the bare operator structure from target
    if isempty(target_norm.terms)
        return nothing
    end
    target_bares = first(target_norm.terms)[1].bares

    for (term, coeff) in expr.terms
        if term.bares == target_bares
            # Build full coefficient including params
            full_coeff = coeff
            for p in term.params
                full_coeff = full_coeff * param_to_symbolic(p)
            end
            return full_coeff isa Num ? simplify(full_coeff) : full_coeff
        end
    end

    return nothing
end

"""
    collect_terms(expr::QuExpr)

Collect and display all terms in a QuExpr with their simplified coefficients.
Returns a vector of (operator_string, simplified_coefficient) pairs.
"""
function collect_terms(expr::QuExpr)
    results = Tuple{String,Any}[]

    for (term, coeff) in expr.terms
        # Build full coefficient
        full_coeff = coeff
        if coeff isa Num || !isempty(term.params)
            full_coeff = coeff isa Num ? coeff : Num(coeff)
            for p in term.params
                full_coeff = full_coeff * param_to_symbolic(p)
            end
            full_coeff = simplify(full_coeff)
        end

        # Get operator string
        op_str = isempty(term.bares.v) ? "𝟙" : string(term.bares)

        push!(results, (op_str, full_coeff))
    end

    return results
end
