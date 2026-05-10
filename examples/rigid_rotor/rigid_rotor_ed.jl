using LinearAlgebra
using SparseArrays

# this returns a matrix in the Hilbert space of constant m
function gen_cosθ(l_max::Int; m=0, l_min=0)
  l_all = l_min:l_max
  N = length(l_all)
  function cosθ_element(l, m, l_)
    if l_ == l + 1
      return sqrt((l + 1)^2 - m^2) / sqrt((2l + 1) * (2l + 3))
    elseif l_ == l - 1
      return sqrt(l^2 - m^2) / sqrt((2l + 1) * (2l - 1))
    else
      return 0.0
    end
  end

  V = spzeros(Float64, N, N)
  for (il1, l1) ∈ enumerate(l_all) # incoming l
    for l2 ∈ (l1 - 1, l1 + 1) # outgoing l
      if l2 < l_min || l2 > l_max
        continue
      end
      il2 = l2 - l_min + 1
      val = cosθ_element(l1, m, l2)
      if val != 0.0
        V[il1, il2] = val
      end
    end
  end
  return V
end

# this returns a matrix in the Hilbert space of constant m
function gen_cosθ²(l_max::Int; m=0, l_min=0)
  cosθ = gen_cosθ(l_max; m = m, l_min = l_min)
  return cosθ * cosθ
end


function gen_a(n_max::Int)
  return spdiagm(+1 => sqrt.(1:n_max))  # places coeffs on the first superdiagonal
end

function get_H(g::Float64, ωc::Float64; B=1.0, l_max=20, n_max=20, gauge_term=true)
  l_all = 0:l_max # angular momentum sector
  cosθ = gen_cosθ(l_max)
  cosθ² = gen_cosθ²(l_max)
  L² = spdiagm(0 => l_all .* (l_all .+ 1))
  a = gen_a(n_max)
  id_ph = spdiagm(0 => ones(n_max + 1))
  id_l = spdiagm(0 => ones(l_max + 1))
  term = g^2 / ωc * kron(id_ph, cosθ²)
  if !gauge_term
    term = 0 * term
  end
  return B * kron(id_ph, L²) + ωc * kron(a' * a, id_l) + g * kron(a + a', cosθ) + term
end
