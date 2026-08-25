module KS_Mel

using CairoMakie, LaTeXStrings
using FastGaussQuadrature, DataStructures, LinearAlgebra
using JuMP, MosekTools, Dualization

const MatrixNPA = include("../../../NPA-MP/src/MatrixNPA.jl").MatrixNPA
using .MatrixNPA

const T = Float64

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

function build_variables!(vars::Variables)
    η::Array{Monome{T},4} = fill(vars.id, (2, 2, 2, 2))
    for b in 1:2, y in 1:2
        for i in 1:2, j in i:2
            η[i, j, b, y] = create_var!(vars; hermitian=(i == j), name="η_($i,$j,$b,$y)")
        end
        η[2, 1, b, y] = η[1, 2, b, y]'
    end

    return η
end

function subs(; vars, η, P)
    subs = []
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
        for b in 1:2, y in 1:2
            if j == k
                push!(subs, (η[i, j, b, y] * η[k, l, b, y], η[i, l, b, y]))
            end
        end
    end

    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
        for y in 1:2, b1 in 1:2, b2 in 1:2
            if b1 == b2
                continue
            end
            push!(subs, (η[i, j, b1, y] * η[k, l, b2, y], vars.zero))
        end
    end

    r = size(P, 1)

    # Projective constraints for P
    for a in 1:2, k in 1:r
        push!(subs, (P[k, a] * P[k, a], P[k, a]))
    end

    # Commutation between P and η
    for i in 1:2, j in 1:2
        for b in 1:2, y in 1:2, k in 1:r, a in 1:2
            push!(subs, (η[i, j, b, y] * P[k, a], P[k, a] * η[i, j, b, y]))
        end
    end

    return subs
end

function op_constraints(η)::Vector{Polynome{T}}
    cs = []
    for y in 1:2
        c = 0
        for i in 1:2, b in 1:2
            c += η[i, i, b, y]
        end
        push!(cs, c - 1)
        push!(cs, 1 - c)
    end

    for i in 1:2, j in 1:2
        c = 0
        for b in 1:2
            c += η[i, j, b, 1] - η[i, j, b, 2]
        end
        c = c + c'
        push!(cs, c)
        push!(cs, -c)
    end

    return cs
end

function objective(; Ms, P, η, α, β)
    obj = 0
    y0 = 1
    for a in 1:2
        for i in 1:2, j in 1:2, b in 1:2
            obj += -α * η[i, j, b, y0] * Ms[1][a][i, j] * P[a]
        end
        obj += -β * P[a]
    end
    return obj
end

function qber_constraints(; Ms, η, Q)::Vector{Polynome{T}}
    cs = []
    for x in 1:2
        y = x

        c = 0
        for a in 1:2, b in 1:2
            if a == b
                continue
            end
            for i in 1:2, j in 1:2
                c += η[i, j, b, y] * Ms[x][a][i, j]
            end
        end
        push!(cs, Q - c)
    end
    return cs
end

function vN_entropy(q)
    vars = Variables{T}()
    Ms = [
        [[1 0; 0 0], [0 0; 0 1]],
        [0.5 * [1 1; 1 1], 0.5 * [1 -1; -1 1]],
    ]
    η = build_variables!(vars)
    η_flat = collect(vec(η))

    dA = 2
    μ = 1e-2
    λ = 1
    ts = grid_function(1, 0.01, μ, λ)
    alphas = [alpha(ts, k) for k in 0:length(ts)]
    betas = [beta(ts, k, μ) for k in 0:length(ts)]
    r = length(alphas)

    P = [create_var!(vars; hermitian=true, name="P[$k,$a]") for k in 1:r, a in 1:2]

    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)
    model = Model(optimizer)
    relaxation = NPARelaxation{T}(;
        model,
        vars,
        subs=subs(; vars, η, P),
    )
    obj = 0
    for (k, (α, β)) in enumerate(zip(alphas, betas))
        obj += objective(; Ms, P=P[k, :], η, α, β)

        b1 = build_basis(η_flat, 1)
        b2 = build_basis(P[k, :], 1)
        basis = prod_basis(b1, b2; subs=relaxation.subs)
        add_moment_matrix_constraint!(relaxation, basis)
    end

    basis = build_basis(η_flat, 1)
    for q in op_constraints(η)
        add_localizing_matrix_constraint!(relaxation, q, basis)
    end
    for r in qber_constraints(; Ms, η, Q=q/2)
        add_moment_vector_constraint!(relaxation, r)
    end

    set_objective!(relaxation, obj)

    JuMP.optimize!(relaxation.model)
    ent = 1 / log(2) * (dA - 1 + objective_value(relaxation.model))
    return ent
end

function test()
    @show vN_entropy(0.2)
end

end # module KS_Mel

