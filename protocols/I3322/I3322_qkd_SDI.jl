module I3322_QKD_SDI

export vN_entropy_I3322, objective_A, objective_B, HAgB

using LinearAlgebra
using JuMP, Dualization
using MosekTools, SCS

# const Ent_BFF = include("../common/ent_BFF.jl").Ent_BFF
# using .Ent_BFF

const Ent_KS = include("../common/ent_KS.jl").Ent_KS
using .Ent_KS

function bloch_mmt(alpha::Number)
    id = [1 0; 0 1]
    sx = [0 1; 1 0]
    sz = [1 0; 0 -1]
    A = cos(alpha)*sz + sin(alpha)*sx
    return [0.5*(id + A), 0.5*(id - A)]
end

function cond_prob(; x, y, a, b, q)
    z = [1, 0]
    o = [0, 1]
    phi = (kron(z, o) - kron(o, z))/sqrt(2)
    rho = (1 - q)*phi*phi' + q*Matrix(I, 4, 4)/4
    angles_A = [0.0, pi/3, 2pi/3]
    angles_B = [4pi/3, pi, 2pi/3]

    Ms = bloch_mmt(angles_A[x])
    Ns = bloch_mmt(angles_B[y])
    return tr(kron(Ms[a], Ns[b])*rho)
end

function vN_entropy_I3322(objective, q)
    optimizer = optimizer_with_attributes(SCS.Optimizer, "eps_abs" => 1e-5, "verbose" => false)
    # optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    # optimizer = dual_optimizer(optimizer)

    dA = 2
    Ms = [bloch_mmt(angles_A[x]) for x in 1:3]
    function constraints(Ns)
        id = Matrix(I, dA, dA)
    
        cs = []
        for x in 1:3, y in 1:3, a in 1:dA, b in 1:dA
            c = kron(Ms[x][a], Ns[y][b]) - cond_prob(; x, y, a, b, q)*id
            push!(cs, +c)
            push!(cs, -c)
        end
    
        return cs
    end

    angles_A = [0.0, pi/3, 2pi/3]
    return vN_entropy(ComplexF64; Ms=Ms[1], dY=3, dB=2, objective, constraints, optimizer)[1]
end

function HAgB(q::Number)
    HAB = 0.0
    for a in 1:dA, b in 1:dA
        p = cond_prob(; x=1, y=1, a, b, q)
        HAB += (p != 0) ? -p*log2(p) : 0
    end

    HB = 0.0
    for b in 1:dA
        p = 0.0
        for a in 1:dA
            p += cond_prob(; x=1, y=1, a, b, q)
        end
        HB += (p != 0) ? -p*log2(p) : 0
    end

    return HAB - HB
end

end # module I3322_QKD_SDI

