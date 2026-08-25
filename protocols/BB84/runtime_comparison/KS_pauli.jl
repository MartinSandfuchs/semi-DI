module KS_Pauli

using DataStructures, LinearAlgebra
using JuMP, MathOptInterface, Dualization
using MosekTools, SCS
using InteractiveUtils

const MatrixNPA = include("../../../NPA-MP/src/MatrixNPA.jl").MatrixNPA
using .MatrixNPA

const T = ComplexF64

function alpha(t, k)
    r = length(t)
    if k == 0
        return -1
    elseif k == 1
        return -((1 + t[1] / (t[2] - t[1])) * log(t[2] / t[1]) - 1)
    elseif k == r
        return -(1 - t[r-1] / (t[r] - t[r-1]) * log(t[r] / t[r-1]))
    end
    return -((1 + t[k] / (t[k+1] - t[k])) * log(t[k+1] / t[k]) - t[k-1] / (t[k] - t[k-1]) * log(t[k] / t[k-1]))
end

function beta(t, k, μ)
    r = length(t)
    if k == 0
        return μ
    elseif k == 1
        return ((1 + t[1] / (t[2] - t[1])) * log(t[2] / t[1]) - 1) * t[1]
    elseif k == r
        return (1 - t[r-1] / (t[r] - t[r-1]) * log(t[r] / t[r-1])) * t[r]
    end
    return ((1 + t[k] / (t[k+1] - t[k])) * log(t[k+1] / t[k]) - t[k-1] / (t[k] - t[k-1]) * log(t[k] / t[k-1])) * t[k]
end

function grid_function(c, ϵ, μ, λ)
    grid = [μ]
    f = μ
    i = 1
    while f < λ
        f = grid[i] + sqrt(grid[i] * ϵ / c)
        push!(grid, f)
        i += 1
    end
    grid[end] = λ  # Ensure the last element is exactly λ
    return grid
end

const mus = ["id", "X", "Y", "Z"]

function build_variables!(vars::Variables)
    S = Dict{Tuple{String,Int,Int},Monome{T}}()
    for mu in mus, b in 1:2, y in 1:2
        var = create_var!(vars; hermitian=true, name="S[$mu, $b, $y]")
        S[(mu, b, y)] = var
    end
    return S
end

function subs(; vars::Variables, S, P)
    subs = []
    r = size(P, 1)

    # Projective constraints for P
    for a in 1:2, k in 1:r
        push!(subs, (P[k, a] * P[k, a], P[k, a]))
    end

    # Commutation between P and S
    for mu in mus, a in 1:2, k in 1:r, b in 1:2, y in 1:2
        push!(subs, (P[k, a] * S[(mu, b, y)], S[(mu, b, y)] * P[k, a]))
    end

    ### Constraints on S
    for b1 in 1:2, b2 in 1:2, y in 1:2, mu in mus
        if b1 == b2
            push!(subs, (S[("id", b1, y)] * S[(mu, b2, y)], S[(mu, b1, y)]))
            push!(subs, (S[(mu, b2, y)] * S[("id", b1, y)], S[(mu, b1, y)]))
        else
            push!(subs, (S[("id", b1, y)] * S[(mu, b2, y)], vars.zero))
            push!(subs, (S[(mu, b2, y)] * S[("id", b1, y)], vars.zero))
        end
    end
    # Pauli constraints
    for mu in ["X", "Y", "Z"], b in 1:2, y in 1:2
        push!(subs, (S[(mu, b, y)] * S[(mu, b, y)], S[("id", b, y)]))
    end
    for b1 in 1:2, b2 in 1:2, y in 1:2
        if b1 == b2
            push!(subs, (S[("X", b1, y)] * S[("Y", b2, y)], 1im * S[("Z", b1, y)]))
            push!(subs, (S[("Y", b1, y)] * S[("Z", b2, y)], 1im * S[("X", b1, y)]))
            push!(subs, (S[("Z", b1, y)] * S[("X", b2, y)], 1im * S[("Y", b1, y)]))
        else
            push!(subs, (S[("X", b1, y)] * S[("Y", b2, y)], vars.zero))
            push!(subs, (S[("Y", b1, y)] * S[("Z", b2, y)], vars.zero))
            push!(subs, (S[("Z", b1, y)] * S[("X", b2, y)], vars.zero))
        end
    end

    return subs
end

function op_constraints(S)::Vector{Polynome{T}}
    cs = []
    for y in 1:2
        c = 0
        for b in 1:2
            c += S[("id", b, y)]
        end
        push!(cs, c - 1)
        push!(cs, 1 - c)
    end

    for mu in mus, y1 in 1:2, y2 in 1:2
        c = 0
        for b in 1:2
            c += S[(mu, b, y1)] - S[(mu, b, y2)]
        end
        push!(cs, c)
        push!(cs, -c)
    end

    return cs
end

function qber_constraints(; lambda, S, Q)::Vector{Polynome{T}}
    cs = []
    for x in 1:2
        y = x

        c = 0
        for a in 1:2, b in 1:2
            if a == b
                continue
            end
            for mu in mus
                c += S[(mu, b, y)] * lambda[(mu, a, x)]
            end
        end
        push!(cs, Q - c)
    end
    return cs
end

function objective(; lambda, S, P, α, β)
    obj = 0
    x0 = 1
    y0 = 1
    for a in 1:2
        for mu in mus, b in 1:2
            obj += -α * lambda[(mu, a, x0)] * S[(mu, b, y0)] * P[a]
        end
        obj += -β * P[a]
    end

    return obj
end

function compute_lambda()
    lambda = Dict{Tuple{String,Int,Int},T}()

    # x=1 => Measure in Z basis
    lambda[("id", 1, 1)] = 0.5
    lambda[("Z", 1, 1)] = 0.5
    lambda[("X", 1, 1)] = lambda[("Y", 1, 1)] = 0.0

    lambda[("id", 2, 1)] = 0.5
    lambda[("Z", 2, 1)] = -0.5
    lambda[("X", 2, 1)] = lambda[("Y", 2, 1)] = 0.0


    # x=2 => Measure in X basis
    lambda[("id", 1, 2)] = 0.5
    lambda[("X", 1, 2)] = 0.5
    lambda[("Y", 1, 2)] = lambda[("Z", 1, 2)] = 0.0

    lambda[("id", 2, 2)] = 0.5
    lambda[("X", 2, 2)] = -0.5
    lambda[("Y", 2, 2)] = lambda[("Z", 2, 2)] = 0.0

    return lambda
end

function vN_entropy(q)
    vars = Variables{T}()

    dA = 2
    μ = 1e-2
    λ = 1
    ts = grid_function(1, 0.01, μ, λ)
    alphas = [alpha(ts, k) for k in 0:length(ts)]
    betas = [beta(ts, k, μ) for k in 0:length(ts)]
    r = length(alphas)

    S = build_variables!(vars)
    S_flat = collect(values(S))
    P = [create_var!(vars; hermitian=true, name="P[$k,$a]") for k in 1:r, a in 1:2]
    lambda = compute_lambda()

    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)
    # optimizer = optimizer_with_attributes(SCS.Optimizer, "eps_abs" => 1e-5, "verbose" => false)
    model = Model(optimizer)

    relaxation = NPARelaxation{T}(;
        model,
        vars,
        subs=subs(; vars, S, P),
    )

    obj = 0
    for (k, (α, β)) in enumerate(zip(alphas, betas))
        obj += objective(; lambda, S, P=P[k, :], α, β)

        b1 = build_basis(S_flat, 1)
        b2 = build_basis(P[k, :], 1)
        basis = prod_basis(b1, b2; subs=relaxation.subs)
        add_moment_matrix_constraint!(relaxation, basis)
    end

    basis = build_basis(S_flat, 1)
    for q in op_constraints(S)
        add_localizing_matrix_constraint!(relaxation, q, basis)
    end

    for r in qber_constraints(; lambda, S, Q=q/2)
        add_moment_vector_constraint!(relaxation, r)
    end

    set_objective!(relaxation, obj)

    JuMP.optimize!(relaxation.model)
    ent = 1 / log(2) * (dA - 1 + objective_value(relaxation.model))
    return ent
end

function test()
    q = 0.2
    @show vN_entropy(q)
end

end # module KS_Pauli

