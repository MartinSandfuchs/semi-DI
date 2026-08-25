module Ent_BFF

export vN_entropy, objective_A, objective_B

using FastGaussQuadrature, DataStructures, LinearAlgebra
using JuMP, Dualization
using SCS

const MatrixNPA = include("../../NPA-MP/src/MatrixNPA.jl").MatrixNPA
using .MatrixNPA

function get_quad(m)
    ts, ws = gaussradau(m)
    return 0.5*(1 .- reverse(ts)), 0.5*reverse(ws)
end


function substitutions!(subs; Ns, Zs)
    dX = dY = length(Ns)
    dA = dB = length(Ns[1])

    # Commutation between N and Z
    m = size(Zs, 1)
    for y in 1:dY, i in 1:m, a in 1:dA
        for b in 1:(dB - 1)
            push!(subs, (Ns[y][b]*Zs[i,a], Zs[i,a]*Ns[y][b]))
            push!(subs, (Ns[y][b]*Zs[i,a]', Zs[i,a]'*Ns[y][b]))
        end
    end

    return subs
end

function extra_monomials(Ns, Zs)
    monos = []
    dA = dB = length(Ns[1])

    # Add monomials appearing in objective function
    for Z in Zs
        for b in 1:(dB - 1)
            push!(monos, Ns[1][b]*Z'*Z)
        end
    end

    return monos
end

function objective_A(; t, Ms, Ns, Zs)
    m, n = size(Ms[1])
    obj = zeros((m, n))
    dA = length(Ms)
    M_tot = sum(Ms)

    for a in 1:dA
        obj += Ms[a] + kron(Ms[a], Zs[a] + Zs[a]' + (1 - t)*Zs[a]'*Zs[a]) + kron(M_tot, t*Zs[a]*Zs[a]')
    end
    return obj
end

function objective_B(; t, Ms, Ns, Zs)
    dA = length(Ms)
    M_tot = sum(Ms)

    obj = 0
    for a in 1:dA
        obj += Ns[1][a] + Ns[1][a]*(Zs[a] + Zs[a]' + (1 - t)*Zs[a]'*Zs[a]) + t*Zs[a]*Zs[a]'
    end
    return kron(M_tot, obj)
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

    m = 8
    vars = Variables{T}()
    subs = []

    Ns = [create_proj_mmt!(vars, subs, dB) for y in 1:dY]
    Zs = [create_var!(vars; hermitian=false) for i in 1:m, a in 1:dA]
    substitutions!(subs; Ns, Zs)

    Ns_flat = vec([Ns[y][b] for y in 1:dY, b in 1:(dA - 1)])

    relaxation = NPARelaxation{T, dQ}(;
        model = Model(optimizer),
        vars,
        subs,
    )

    obj = zeros((dQ, dQ))
    ts, ws = get_quad(m)
    for (i, (t_i, w_i)) in enumerate(zip(ts, ws))
        obj_i = objective(; t=t_i, Ms, Ns, Zs=Zs[i, :])
        obj += w_i/(t_i*log(2))*obj_i

        # Add moment constraint
        extra::Vector{Monome{T}} = extra_monomials(Ns, Zs[i, :])
        basis = build_basis([Ns_flat; Zs[i, :]], 2, relaxation.subs, extra)
        # b1 = build_basis(As, 1)
        # b2 = build_basis(Zs[i, :], 1)
        # basis = prod_basis(b1, b2, relaxation.subs)
        add_moment_matrix_constraint!(relaxation, basis)
    end

    delta = isnothing(offset) ? zeros(dQ, dQ) : offset(Ns)
    set_objective!(relaxation, obj + delta)

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

end # module Ent_BFF


