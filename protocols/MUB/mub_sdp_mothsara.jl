module MUB_SDP

using JuMP
using Mosek
using MosekTools
using SCS
using Dualization
using LinearAlgebra
using Ket

const MUB_Common = include("mub_common.jl").MUB_Common
using .MUB_Common

# Setup Bob's measurements
function Ns(d::Int)
    B = Ket.mub(d)
    return [
        begin
            v = B[y][:, b]
            v*v'
        end
        for y in 1:(d+1), b in 1:d
    ]
end

# Compute min-entropy in dimension d=3
function hmin(d::Int; q)
    GC.gc()
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    # optimizer = optimizer_with_attributes(SCS.Optimizer)
    optimizer = dual_optimizer(optimizer)
    model = Model(optimizer)

    σ = [
        begin
            @variable(model, [1:d, 1:d] in HermitianPSDCone())
        end
        for x in 1:(d+1), a in 1:d, e in 1:d
    ]

    N = Ns(d)

    # Set the objective
    obj = 0
    for a in 1:d, b in 1:d
        obj += tr(N[1,b]*σ[1,a,b])
    end
    @objective(model, Max, real(obj))

    # Add non-signaling constraints
    for x1 in 1:(d+1), x2 in 1:(d+1), e in 1:d
        lhs = zeros(d, d)
        rhs = zeros(d, d)
        for a in 1:d
            lhs += σ[x1,a,e]
            rhs += σ[x2,a,e]
        end
        @constraint(model, lhs .== rhs)
    end

    # Add normalization constraint
    for x in 1:(d+1)
        total = 0
        for a in 1:d, e in 1:d
            total += tr(σ[x,a,e])
        end
        @constraint(model, real(total) <= 1)
    end

    # Add statistics constraints
    for x in 1:(d+1), y in 1:(d+1), a in 1:d, b in 1:d
        p = 0
        for e in 1:d
            p += tr(N[y,b]*σ[x,a,e])
        end
        @constraint(model, p == cond_prob(; d, x, y, a, b, q))
    end

    optimize!(model)
    min_ent = -log2(objective_value(model))
    return min_ent
end

function key_rate(d::Int=3; q)
    h = hmin(d; q)
    ec = HAgB(d; q)
    return h - ec
end

end # module MUB_SDP
