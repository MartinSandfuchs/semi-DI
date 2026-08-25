module CHSH_QKD_SDI

using LinearAlgebra
using JuMP, Dualization, MosekTools

export vN_entropy_CHSH, objective_A, objective_B

const Ent_BFF = include("../common/ent_BFF.jl").Ent_BFF
using .Ent_BFF

# const Ent_KS = include("../common/ent_KS.jl").Ent_KS
# using .Ent_KS

function cond_prob(; x, y, a, b, q)
    psi = [1, 0, 0, 1]/sqrt(2)
    rho = (1 - q)*psi*psi' + q*Matrix(I, 4, 4)/4
    id = Matrix(I, 2, 2)
    sx = [0 1; 1 0]
    sz = [1 0; 0 -1]
    Ms = [
        [0.5*(id + sx), 0.5*(id - sx)],
        [0.5*(id + sz), 0.5*(id - sz)],
    ]
    Ns = [
        [0.5*(id + (sx + sz)/sqrt(2)), 0.5*(id - (sx + sz)/sqrt(2))],
        [0.5*(id + (sx - sz)/sqrt(2)), 0.5*(id - (sx - sz)/sqrt(2))],
    ]
    return tr(kron(Ms[x][a], Ns[y][b])*rho)
end

function winning_prob(q)
    w = 0.0
    for x in 1:2, y in 1:2, a in 1:2, b in 1:2
        if (x - 1)*(y - 1) == mod((a - 1) + (b - 1), 2)
            w += 0.25*cond_prob(; x, y, a, b, q)
        end
    end
    return w
end

function chsh_constraints(; Ms, Ns, w)
    id = Matrix(I, 2, 2)

    cs = []
    c = zeros(2, 2)
    for x in 1:2, y in 1:2, a in 1:2, b in 1:2
        if (x - 1)*(y - 1) == mod((a - 1) + (b - 1), 2)
            c += 0.25*kron(Ms[x][a], Ns[y][b])
        end
    end
    push!(cs, c - w*id)

    return cs
end

function vN_entropy_CHSH(objective, q)
    w = winning_prob(q)
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)

    Ms = [
        [0.5*[1 1; 1 1], 0.5*[1 -1; -1 1]],
        [[1 0; 0 0], [0 0; 0 1]],
    ]
    function constraints(Ns)
        return chsh_constraints(; Ms, Ns, w)
    end

    return vN_entropy(; Ms=Ms[1], dY=2, dB=2, objective, constraints, optimizer)[1]
end

end # module CHSH_QKD_SDI

