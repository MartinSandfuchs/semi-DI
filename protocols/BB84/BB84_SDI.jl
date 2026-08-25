module BB84_SDI

export cond_prob, vN_entropy_BB84, HAgB, objective_A, objective_B

using LinearAlgebra
using JuMP, Dualization, MosekTools
using CairoMakie, LaTeXStrings

const Ent_BFF = include("../common/ent_BFF.jl").Ent_BFF
using .Ent_BFF

# const Ent_KS = include("../common/ent_KS.jl").Ent_KS
# using .Ent_KS

function cond_prob(; x, y, a, b, q)
    Ms = [
        [[1 0; 0 0], [0 0; 0 1]],
        [0.5*[1 1; 1 1], 0.5*[+1 -1; -1 +1]],
    ]
    Ns = copy(Ms)
    phi = [1, 0, 0, 1]/sqrt(2)
    rho = (1 - q)*phi*phi' + q*Matrix(I, 4, 4)/4

    return tr(kron(Ms[x][a], Ns[y][b]) * rho)
end

function vN_entropy_BB84(objective, q)
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)

    Ms = [
        [[1 0; 0 0], [0 0; 0 1]],
        [0.5*[1 1; 1 1], 0.5*[+1 -1; -1 +1]],
    ]
    function constraints(Ns)
        # eps = 1e-3
        eps = 0
        id = Matrix(I, 2, 2)
        cs = []
        for x in 1:2, y in 1:2, a in 1:2, b in 1:2
            c = kron(Ms[x][a], Ns[y][b]) - id*cond_prob(; x, y, a, b, q)
            push!(cs, +c + id*eps)
            push!(cs, -c + id*eps)
        end
        # for x in 1:2
        #     y = x
        #     c = kron(Ms[x][1], Ns[y][2]) + kron(Ms[x][2], Ns[y][1])
        #     push!(cs, q/2*id - c)
        # end
        return cs
    end

    return vN_entropy(; Ms=Ms[1], dB=2, dY=2, objective, constraints, optimizer)[1]
end

function h2(p)
    if p == 0 || p == 1
        return 0.0
    else
        return -p*log2(p) - (1-p)*log2(1-p)
    end
end

HAgB(q) = h2(q/2)

# q = 0.2
# @show vN_entropy_BB84(objective_A, q)
# @show vN_entropy_BB84(objective_B, q)
# @show 1 - h2(q/2)

end # module BB84_SDI

