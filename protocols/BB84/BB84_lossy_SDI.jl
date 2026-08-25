module BB84_SDI

export vN_entropy_lossy, rate_lossy

using CairoMakie, LaTeXStrings
using FastGaussQuadrature, DataStructures, LinearAlgebra
using JuMP, MosekTools, Dualization

const Ent_BFF = include("../common/ent_BFF.jl").Ent_BFF
using .Ent_BFF

# const Ent_KS = include("../common/ent_KS.jl").Ent_KS
# using .Ent_KS

function cond_prob(; x, y, a, b, q, eta)
    phi = vec(Matrix(I, 2, 2)/sqrt(2))
    rho = (1 - q)*phi*phi' + q*Matrix(I, 4, 4)/4
    Ms = [
        [[1 0; 0 0], [0 0; 0 1]],
        [0.5*[1 1; 1 1], 0.5*[1 -1; -1 1]],
    ]
    if a == 3
        return (1 - eta)*0.5
    else
        return eta*tr(kron(Ms[x][a], Ms[y][b])*rho)
    end
end

function vN_entropy_lossy(q, eta)
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)

    Ms = [
        [
            [1 0 0; 0 0 0; 0 0 0],
            [0 0 0; 0 1 0; 0 0 0],
        ],
        [
            0.5*[1 1 0; 1 1 0; 0 0 0],
            0.5*[1 -1 0; -1 1 0; 0 0 0],
        ],
    ]
    function constraints(Ns)
        # eps = 1e-3
        eps = 0
        id = Matrix(I, 3, 3)
        cs = []
        for x in 1:2, y in 1:2, a in 1:2, b in 1:2
            c = kron(Ms[x][a], Ns[y][b]) - id*cond_prob(; x, y, a, b, q, eta)
            push!(cs, +c + id*eps)
            push!(cs, -c + id*eps)
        end
        return cs
    end

    return vN_entropy(; Ms=Ms[1], dY=2, dB=2, objective=objective_A, constraints, optimizer)[1]
end

function HAgB(q, eta)
    HAB = 0
    for a in 1:2, b in 1:2
        p = cond_prob(; x=1, y=1, a, b, q, eta)
        if iszero(p)
            continue
        end
        HAB += -p*log2(p)
    end

    HB = 0
    for b in 1:2
        p = 0
        for a in 1:2
            p += cond_prob(; x=1, y=1, a, b, q, eta)
        end
        if iszero(p)
            continue
        end
        HB += -p*log2(p)
    end

    return HAB - HB
end

function rate_lossy(q, eta)
    H_AE = vN_entropy_lossy(q, eta)
    H_AB = HAgB(q, eta)
    return H_AE - H_AB
end

# const T = Float64
# # const T = ComplexF64
# 
# function get_quad(m)
#     ts, ws = gaussradau(m)
#     return 0.5*(1 .- reverse(ts)), 0.5*reverse(ws)
# end
# 
# function subs(As, Zs)
#     subs = Vector{Tuple{Monome{T}, Monome{T}}}()
# 
#     # Projective constraints
#     for A in As
#         push!(subs, (A*A, A))
#     end
# 
#     # Commutation with Eve's operators
#     for A in As, Z in Zs
#         push!(subs, (A*Z, Z*A))
#         push!(subs, (A*Z', Z'*A))
#     end
# 
#     return subs
# end
# 
# function objective(; t, As, Ns, Zs)::Matrix{Polynome{T}}
#     Ms = [As[1], (1 - As[1])]
#     N_tot = Ns[1][1] + Ns[1][2]
# 
#     obj = zeros(3, 3)
#     for a in 1:2
#         obj += kron(N_tot, Ms[a] + Ms[a]*(Zs[a] + Zs[a]' + (1 - t)*Zs[a]'*Zs[a]) + t*Zs[a]*Zs[a]')
#     end
#     return obj
# end
# 
# function constraints(; As, Ns, q, eta)
#     Ms = [[As[1], (1 - As[1])], [As[2], (1 - As[2])]]
#     id = Matrix(I, 3, 3)
#     cs = []
#     # eps = 1e-4
#     eps = 0
#     for x in 1:2, y in 1:2, a in 1:2, b in 1:3
#         c = kron(Ns[y][b], Ms[x][a]) - id*cond_prob(; x, y, a, b, q, eta)
#         push!(cs, +(c + id*eps))
#         push!(cs, -(c - id*eps))
#     end
#     return cs
# end
# 
# function vN_entropy_lossy(q, eta)
#     m = 8
#     vars = Variables{T}()
#     As = [create_var!(vars; hermitian=true) for x in 1:2]
#     Ns = [
#         [
#             [1 0 0; 0 0 0; 0 0 0],
#             [0 0 0; 0 1 0; 0 0 0],
#             [0 0 0; 0 0 0; 0 0 1],
#         ],
#         [
#             0.5*[1 1 0; 1 1 0; 0 0 0],
#             0.5*[1 -1 0; -1 1 0; 0 0 0],
#             [0 0 0; 0 0 0; 0 0 1],
#         ],
#     ]
# 
#     Zs = [create_var!(vars; hermitian=false) for i in 1:m, a in 1:2]
#     id = Matrix(I, 3, 3)
# 
#     optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
#     optimizer = dual_optimizer(optimizer)
#     model = Model(optimizer)
# 
#     substitutions = subs(As, Zs)
#     relaxation = NPARelaxation{T, 3}(;
#         model,
#         vars,
#         substitutions,
#     )
# 
#     obj = zeros(3, 3)
#     ts, ws = get_quad(m)
#     for (i, (t_i, w_i)) in enumerate(zip(ts, ws))
#         obj_i = objective(; t=t_i, As, Ns, Zs=Zs[i, :])
#         obj += w_i/(t_i*log(2))*obj_i
# 
#         # Add moment constraint
#         basis = build_basis([As; Zs[i, :]], 2, relaxation.substitutions)
#         add_moment_matrix_constraint!(relaxation, basis)
#     end
# 
#     set_objective!(relaxation, obj)
# 
#     for c in constraints(; As, Ns, q, eta)
#         add_moment_vector_constraint!(relaxation, c)
#     end
# 
#     JuMP.optimize!(relaxation.model)
#     return objective_value(relaxation.model)
# end
# 
# h_bin(p) = -p*log2(p) - (1-p)*log2(1-p)
# q = 0.2
# eta = 0.1
# @show vN_entropy_lossy(q, eta)/eta
# @show 1 - h_bin(q/2)

q = 0.1
eta = 0.01
@show (vN_entropy_lossy(q, eta) - HAgB(q, eta))/eta

end # module BB84_SDI
