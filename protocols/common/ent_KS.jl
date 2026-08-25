module Ent_KS

export vN_entropy_sdp, vN_entropy, objective_A, objective_B

using LinearAlgebra
using JuMP, Dualization
using MosekTools, SCS

const MatrixNPA = include("../../NPA-MP/src/MatrixNPA.jl").MatrixNPA
using .MatrixNPA

# const μ = 1e-5
const μ = 1e-2

function alpha(t, k)
    r = length(t)
    if k == 0
        return -1
    end
    if k == 1
        return -((1 + t[1]/(t[2]-t[1]))*log(t[2]/t[1]) - 1)
    end
    if k == r
        return -(1 - t[r-1]/(t[r]-t[r-1])*log(t[r]/t[r-1]))
    end
    return -((1 + t[k]/(t[k+1]-t[k]))*log(t[k+1]/t[k]) - t[k-1]/(t[k]-t[k-1])*log(t[k]/t[k-1]))
end

function beta(t, k)
    r = length(t)
    if k == 0
        return μ
    end
    if k == 1
        return ((1 + t[1]/(t[2]-t[1]))*log(t[2]/t[1]) - 1)*t[1]
    end
    if k == r
        return (1 - t[r-1]/(t[r]-t[r-1])*log(t[r]/t[r-1]))*t[r]
    end
    return ((1 + t[k]/(t[k+1]-t[k]))*log(t[k+1]/t[k]) - t[k-1]/(t[k]-t[k-1])*log(t[k]/t[k-1]))*t[k]
end

function grid_function(c, ϵ, μ, λ)
    grid = [μ]
    f = μ
    i = 1
    while f < λ
        f = grid[i] + sqrt(grid[i] * ϵ/c)
        push!(grid, f)
        i += 1
    end
    grid[end] = λ  # Ensure the last element is exactly λ
    return grid
end

function create_proj_mmt!(vars::Variables, subs, n::Int)
    Ms = [create_var!(vars; hermitian=true) for i in 1:(n-1)]
    for i in 1:(n-1), j in 1:(n-1)
        if i == j
            push!(subs, (Ms[i]*Ms[i], Ms[i]))
        else
            push!(subs, (Ms[i]*Ms[j], vars.zero))
        end
    end

    return [Ms; 1 - sum(Ms)]
end

function objective_A(; Ms, Ns, Ps, alpha, beta)
    m, n = size(Ms[1])
    obj = zeros(m, n)
    dA = length(Ms)
    M_tot = sum(Ms)

    for a in 1:dA
        obj += -alpha*kron(Ms[a], Ps[a]) - beta*kron(M_tot, Ps[a])
    end
    return obj
end

function objective_B(; Ms, Ns, Ps, alpha, beta)
    obj = 0
    dA = length(Ms)
    M_tot = sum(Ms)

    for b in 1:dA
        obj += -alpha*Ns[1][b]*Ps[b] - beta*Ps[b]
    end
    return kron(M_tot, obj)
end

function substitutions!(subs; dA, Ns, Ps)
    dY = length(Ns)
    dB = length(Ns[1])

    # Commutation between N and P
    r = size(Ps, 1)
    for y in 1:dY, i in 1:r, a in 1:dA
        for b in 1:(dB - 1)
            push!(subs, (Ns[y][b]*Ps[i,a], Ps[i,a]*Ns[y][b]))
        end
    end

    # Projective constraints for P
    for i in 1:r, a in 1:dA
        push!(subs, (Ps[i,a]*Ps[i,a], Ps[i,a]))
    end
    
    return subs
end

const SCS_optimizer = optimizer_with_attributes(SCS.Optimizer, "eps_abs" => 1e-5, "verbose" => false)

function vN_entropy(
    ::Type{T}=Float64;
    Ms, dY, dB,
    objective, constraints, offset=nothing,
    optimizer = SCS_optimizer,
) where {T}
    @assert size(Ms[1], 1) == size(Ms[1], 2)
    dQ = size(Ms[1], 1)
    dA = length(Ms)

    λ = 1
    # const ts = LinRange(μ, λ, 60)
    ts = grid_function(1, 0.01, μ, λ)
    # ts = grid_function(1, 0.1, μ, λ)
    alphas = [alpha(ts, k) for k in 0:length(ts)]
    betas = [beta(ts, k) for k in 0:length(ts)]
    r = length(alphas)

    vars = Variables{T}()
    subs = []

    Ns = [create_proj_mmt!(vars, subs, dB) for y in 1:dY]
    Ps = [create_var!(vars; hermitian=true) for i in 1:r, a in 1:dA]
    substitutions!(subs; dA, Ns, Ps)

    Ns_flat = vec([Ns[y][b] for y in 1:dY, b in 1:(dA - 1)])

    obj = zeros(dQ, dQ)
    relaxation = NPARelaxation{T, dQ}(;
        model = Model(optimizer),
        vars,
        subs,
    )

    for (k, (alpha_k, beta_k)) in enumerate(zip(alphas, betas))
        obj += objective(; Ms, Ns, Ps=Ps[k, :], alpha=alpha_k, beta=beta_k)

        b1 = build_basis(Ns_flat, 1)
        b2 = build_basis(Ps[k, :], 1)
        basis = prod_basis(b1, b2; subs=relaxation.subs)
        add_moment_matrix_constraint!(relaxation, basis)
    end

    M_tot = sum(Ms)
    delta = isnothing(offset) ? zeros(dQ, dQ) : offset(Ns)
    set_objective!(relaxation, 1/log(2)*(kron(M_tot, vars.id)*(dA - 1) + obj) + delta)

    c_refs = []
    for c in constraints(Ns)
        c_ref = add_moment_vector_constraint!(relaxation, c)
        push!(c_refs, c_ref)
    end

    JuMP.optimize!(relaxation.model)
    ent = objective_value(relaxation.model)
    grad = dual.(c_refs)
    return ent, grad
end

end # module Ent_KS
