module MUB_QKD_SDI

export cond_prob, vN_entropy_MUB, HAgB, objective_A, objective_B

using LinearAlgebra
using JuMP, Dualization
using MosekTools, SCS

# const Ent_BFF = include("../common/ent_BFF.jl").Ent_BFF
# using .Ent_BFF

const Ent_KS = include("../common/ent_KS.jl").Ent_KS
using .Ent_KS

const dA = 3

function bases()::Vector{Vector{Vector{ComplexF64}}}
    ξ = exp(2*π*im/3)
    b = []
    z = ComplexF64[1, 0, 0]
    o = ComplexF64[0, 1, 0]
    t = ComplexF64[0, 0, 1]

    push!(b, [z, o, t])
    push!(b, [(z + o + t)/√(3), (z + ξ^2*o + ξ*t)/√(3), (z + ξ*o + ξ^2*t)/√(3)])
    push!(b, [(z + o + ξ*t)/√(3), (z + ξ^2*o + ξ^2*t)/√(3), (z + ξ*o + t)/√(3)])
    push!(b, [(z + o + ξ^2*t)/√(3), (z + ξ^2*o + t)/√(3), (z + ξ*o + ξ*t)/√(3)])

    return b
end

function cond_prob(; x, y, a, b, q)
    z = [1, 0, 0]
    o = [0, 1, 0]
    t = [0, 0, 1]
    psi = (kron(z, z) + kron(o, o) + kron(t, t))/sqrt(3)
    rho = (1 - q)*psi*psi' + q*Matrix(I, 3^2, 3^2)/3^2
    Bs = bases()
    v1 = Bs[x][a]
    v2 = Bs[y][b]
    return tr(kron(v1*v1', v2*v2') * rho) |> real
end

function vN_entropy_MUB(objective, q)
    optimizer = optimizer_with_attributes(SCS.Optimizer, "eps_abs" => 1e-5, "verbose" => false)
    # optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    # optimizer = dual_optimizer(optimizer)

    Ms = map(vec -> map(v -> v*v', vec), bases())
    function constraints(Ns)
        id = Matrix(I, dA, dA)

        cs = []
        for x in 1:4, y in 1:4, a in 1:dA, b in 1:dA
            c = kron(Ms[x][a], Ns[y][b]) - cond_prob(; x, y, a, b, q)*id
            push!(cs, +c)
            push!(cs, -c)
        end

        return cs
    end

    return vN_entropy(ComplexF64; Ms=Ms[1], dY=4, dB=3, objective, constraints, optimizer)[1]
end

function HAgB(q::Number)
    HAB = 0.0
    for a in 1:3, b in 1:3
        p = cond_prob(; x=1, y=1, a, b, q)
        if p != 0
            HAB += -p*log2(p)
        end
    end

    HB = 0.0
    for b in 1:3
        p = 0.0
        for a in 1:3
            p += cond_prob(; x=1, y=1, a, b, q)
        end
        if p != 0
            HB += -p*log2(p)
        end
    end

    return HAB - HB
end

end # module MUB_QKD_SDI
